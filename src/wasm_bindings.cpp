#include <emscripten.h>
#include <emscripten/bind.h>

#include <chrono>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "Expression.h"
#include "Optimization.h"

using namespace emscripten;

struct NewtonSystem {
  std::string lhs;
  std::string rhs;
  std::string rhsShorthand;
  std::string variables;
  std::string variableDefinitions;
};

struct NewtonSystemTriplet {
  NewtonSystem full;
  NewtonSystem augmented;
  NewtonSystem normal;
};

namespace {
Optimization::Settings changePenaltyToPenaltyWithExtraVariables(
    Optimization::Settings settings) {
  if (settings.equalityHandling ==
      Optimization::EqualityHandling::PenaltyFunction) {
    settings.equalityHandling =
        Optimization::EqualityHandling::PenaltyFunctionWithExtraVariable;
  }
  return settings;
}

std::string replaceAll_(std::string str, const std::string& from,
                        const std::string& to) {
  size_t start_pos = 0;
  while ((start_pos = str.find(from, start_pos)) != std::string::npos) {
    str.replace(start_pos, from.length(), to);
    start_pos +=
        to.length();  // Handles case where 'to' is a substring of 'from'
  }
  return str;
}

std::string addLineBreaks(std::string str,
                          const std::string breakStr = " \\\\\n \\qquad ") {
  return str;
  const auto breaks = std::set<std::string>{" - ", " + "};
  int numParenthesis = 0;

  size_t numBreaks = 0;
  for (size_t i = 0; i < str.length(); ++i) {
    for (const auto& sub : breaks) {
      numBreaks += (str.substr(i, sub.length()) == sub) ? 1 : 0;
    }
  }

  size_t currentBreakIndex = 0;
  size_t lastBreak = 0;
  for (size_t i = 1; i + 2 < str.length(); ++i) {
    const auto c = str[i];
    numParenthesis += c == '(' ? 1 : 0;
    numParenthesis -= c == ')' ? 1 : 0;
    const auto isBreak = breaks.contains(str.substr(i, 3));
    currentBreakIndex += isBreak ? 1 : 0;
    if (isBreak && currentBreakIndex + 1 < numBreaks) {
      if (numParenthesis == 0 && currentBreakIndex > lastBreak + 3) {
        str.insert(i, breakStr);
        lastBreak = currentBreakIndex;
      } else if (numParenthesis == 1 && currentBreakIndex > lastBreak + 5) {
        str.insert(i, breakStr);
        lastBreak = currentBreakIndex;
      }
    }
  }
  return str;
}

NewtonSystem formatNewtonSystemStrings_(const auto& lhs, const auto& rhs,
                                        const auto& variables,
                                        const auto& variableDefinitions) {
  using namespace Expression::ExprFactory;
  const auto unity = number(1);
  const auto I = namedConstant("I");
  const auto mu = namedConstant("\\mu");
  const auto muI = product({mu, I}).simplify();
  std::string lhsStr = "";
  const auto condensed = true;
  for (const auto& row : lhs) {
    for (size_t i = 0; i < row.size(); ++i) {
      auto rowExpr = row[i];
      rowExpr = rowExpr.replaceSubexpression(unity, I);
      rowExpr = rowExpr.replaceSubexpression(mu, muI);
      lhsStr +=
          rowExpr.toString(condensed) + (i + 1 == row.size() ? "" : " & ");
    }
    lhsStr += " \\\\\n ";
  }

  std::string rhsStr = "";
  for (const auto& row : rhs) {
    auto rowStr = addLineBreaks(row.toString(condensed));
    rhsStr += rowStr + " \\\\\n ";
  }

  const auto rhsShorthand = Optimization::getShorthandRhs(variables);
  std::string rhsShorthandStr = "";
  for (const auto& row : rhsShorthand) {
    rhsShorthandStr += addLineBreaks(row.toString(condensed)) + " \\\\\n ";
  }

  std::string variablesStr = "";
  for (size_t i = 0; i < variables.size(); ++i) {
    variablesStr += "\\Delta " + variables[i].toString(condensed) +
                    (i + 1 == variables.size() ? "\n" : " \\\\\n ");
  }

  std::string definitionsStr = "";
  for (size_t i = variableDefinitions.size(); i > 0; --i) {
    definitionsStr +=
        variableDefinitions.at(i - 1).first.toString(condensed) + " &= " +
        addLineBreaks(variableDefinitions.at(i - 1).second.toString(condensed),
                      " \\\\\n & \\qquad ") +
        (i - 1 == 0 ? "\n" : " \\\\\n ");
  }

  return {lhsStr, rhsStr, rhsShorthandStr, variablesStr, definitionsStr};
}

enum class NewtonSystemType { Full, Augmented, Normal };

NewtonSystem getNewtonSystem_(const Optimization::Settings& settingsIn,
                              const NewtonSystemType type) {
  const auto settings = changePenaltyToPenaltyWithExtraVariables(settingsIn);
  const auto begin = std::chrono::high_resolution_clock::now();
  const auto variableNames = Optimization::VariableNames();
  auto [lagrangian, variables] =
      Optimization::getLagrangian(variableNames, settings);
  auto [lhs, rhs] = Optimization::getNewtonSystem(lagrangian, variables);
  const auto endNewton = std::chrono::high_resolution_clock::now();

  const auto augmentedSize = [&]() {
    const auto zero = Expression::ExprFactory::number(0);
    const auto unity = Expression::ExprFactory::number(1);
    const auto negUnity =
        Expression::ExprFactory::negate(Expression::ExprFactory::number(1));
    const auto reducibles = std::set{zero, unity, negUnity};
    size_t i = 0;
    while (i < lhs.size() && !reducibles.contains(lhs.at(0).at(i))) {
      ++i;
    }
    return i;
  }();
  const auto endAugmented = std::chrono::high_resolution_clock::now();

  std::vector<std::pair<Expression::Expr, Expression::Expr>>
      variableDefinitions;
  std::chrono::milliseconds gaussianCount =
      std::chrono::duration_cast<std::chrono::milliseconds>(endNewton -
                                                            endNewton);
  std::chrono::milliseconds deltaCount =
      std::chrono::duration_cast<std::chrono::milliseconds>(endNewton -
                                                            endNewton);
  if (type != NewtonSystemType::Full) {
    rhs = Optimization::getShorthandRhs(variables);
    while (lhs.size() > augmentedSize) {
      auto deltaVariable = Expression::ExprFactory::variable(
          "\\Delta " + variables.at(lhs.size() - 1).getName());
      const auto beginDelta = std::chrono::high_resolution_clock::now();
      auto deltaDefinition =
          Optimization::deltaDefinition(lhs, rhs, variables, lhs.size() - 1);
      const auto endDelta = std::chrono::high_resolution_clock::now();
      variableDefinitions.push_back({deltaVariable, deltaDefinition});
      Optimization::gaussianElimination(lhs, rhs, lhs.size() - 1);
      const auto endGaussian = std::chrono::high_resolution_clock::now();
      gaussianCount += std::chrono::duration_cast<std::chrono::milliseconds>(
          endGaussian - endDelta);
      deltaCount += std::chrono::duration_cast<std::chrono::milliseconds>(
          endDelta - beginDelta);

      variables.pop_back();
    }

    if (type == NewtonSystemType::Normal) {
      auto deltaVariable = Expression::ExprFactory::variable(
          "\\Delta " + variables.at(0).getName());
      auto deltaDefinition =
          Optimization::deltaDefinition(lhs, rhs, variables, 0);
      variableDefinitions.push_back({deltaVariable, deltaDefinition});
      Optimization::gaussianElimination(lhs, rhs, 0);
      variables.erase(variables.begin());
    }
  }
  const auto end = std::chrono::high_resolution_clock::now();
  std::cout << "time: newton: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(endNewton -
                                                                     begin)
                   .count()
            << ", augmented: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(
                   endAugmented - endNewton)
                   .count()
            << ", delta: " << deltaCount << ", gaussain: " << gaussianCount
            << ", full: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                     begin)
                   .count()
            << std::endl;

  return formatNewtonSystemStrings_(lhs, rhs, variables, variableDefinitions);
}

std::tuple<NewtonSystem, NewtonSystem, NewtonSystem> getNewtonSystems_(
    const Optimization::Settings& settingsIn) {
  const auto settings = changePenaltyToPenaltyWithExtraVariables(settingsIn);
  const auto begin = std::chrono::high_resolution_clock::now();
  const auto variableNames = Optimization::VariableNames();
  auto [lagrangian, variables] =
      Optimization::getLagrangian(variableNames, settings);
  auto [lhs, rhs] = Optimization::getNewtonSystem(lagrangian, variables);
  const auto endNewton = std::chrono::high_resolution_clock::now();

  const auto augmentedSize = [&]() {
    const auto zero = Expression::ExprFactory::number(0);
    const auto unity = Expression::ExprFactory::number(1);
    const auto negUnity =
        Expression::ExprFactory::negate(Expression::ExprFactory::number(1));
    const auto reducibles = std::set{zero, unity, negUnity};
    size_t i = 0;
    while (i < lhs.size() && !reducibles.contains(lhs.at(0).at(i))) {
      ++i;
    }
    return i;
  }();
  const auto endAugmented = std::chrono::high_resolution_clock::now();

  std::vector<std::pair<Expression::Expr, Expression::Expr>>
      variableDefinitions;
  std::chrono::milliseconds gaussianCount =
      std::chrono::duration_cast<std::chrono::milliseconds>(endNewton -
                                                            endNewton);
  std::chrono::milliseconds deltaCount =
      std::chrono::duration_cast<std::chrono::milliseconds>(endNewton -
                                                            endNewton);

  auto newtonSystem =
      formatNewtonSystemStrings_(lhs, rhs, variables, variableDefinitions);

  rhs = Optimization::getShorthandRhs(variables);
  while (lhs.size() > augmentedSize) {
    auto deltaVariable = Expression::ExprFactory::variable(
        "\\Delta " + variables.at(lhs.size() - 1).getName());
    const auto beginDelta = std::chrono::high_resolution_clock::now();
    auto deltaDefinition =
        Optimization::deltaDefinition(lhs, rhs, variables, lhs.size() - 1);
    const auto endDelta = std::chrono::high_resolution_clock::now();
    variableDefinitions.push_back({deltaVariable, deltaDefinition});
    Optimization::gaussianElimination(lhs, rhs, lhs.size() - 1);
    const auto endGaussian = std::chrono::high_resolution_clock::now();
    gaussianCount += std::chrono::duration_cast<std::chrono::milliseconds>(
        endGaussian - endDelta);
    deltaCount += std::chrono::duration_cast<std::chrono::milliseconds>(
        endDelta - beginDelta);

    variables.pop_back();
  }
  auto augmentedSystem =
      formatNewtonSystemStrings_(lhs, rhs, variables, variableDefinitions);

  auto deltaVariable =
      Expression::ExprFactory::variable("\\Delta " + variables.at(0).getName());
  auto deltaDefinition = Optimization::deltaDefinition(lhs, rhs, variables, 0);
  variableDefinitions.push_back({deltaVariable, deltaDefinition});
  Optimization::gaussianElimination(lhs, rhs, 0);
  variables.erase(variables.begin());

  auto normalEquations =
      formatNewtonSystemStrings_(lhs, rhs, variables, variableDefinitions);
  const auto end = std::chrono::high_resolution_clock::now();
  std::cout << "time: newton: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(endNewton -
                                                                     begin)
                   .count()
            << ", augmented: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(
                   endAugmented - endNewton)
                   .count()
            << ", delta: " << deltaCount << ", gaussian: " << gaussianCount
            << ", full: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                     begin)
                   .count()
            << std::endl;

  return {newtonSystem, augmentedSystem, normalEquations};
}

}  // namespace

std::string getLagrangian(const Optimization::Settings& settings) {
  const auto variableNames = Optimization::VariableNames();
  const auto [lagrangian, variables] =
      Optimization::getLagrangian(variableNames, settings);

  const auto condensed = true;
  auto str = "& " + lagrangian.toString(condensed);
  auto pos = 0;
  const auto addNewlines = [&](const auto& term) {
    pos = str.find(term);
    while (pos != std::string::npos) {
      str.insert(pos - 1, " \\\\\n & ");
      for (size_t i = 0; i < 3; ++i) {
        pos = pos == std::string::npos ? std::string::npos
                                       : str.find(term, pos + 1);
      }
      if (pos != std::string::npos) {
        pos = str.find(term, pos + 1);
      }
    }
  };
  addNewlines("\\lambda");
  addNewlines("- \\mu");
  return str;
}

std::string getFirstOrderOptimalityConditions(
    const Optimization::Settings& settingsIn) {
  const auto variableNames = Optimization::VariableNames();
  const auto settings = changePenaltyToPenaltyWithExtraVariables(settingsIn);
  const auto [lagrangian, variables] =
      Optimization::getLagrangian(variableNames, settings);
  auto firstOrder =
      Optimization::getFirstOrderOptimalityConditions(lagrangian, variables);
  const auto condensed = true;
  std::string str = "";
  for (const auto& c : firstOrder) {
    str += c.toString(condensed) + " &= 0 \\\\";
  }
  return str;
}

NewtonSystem getNewtonSystem(const Optimization::Settings& settings) {
  return getNewtonSystem_(settings, NewtonSystemType::Full);
}

NewtonSystem getAugmentedSystem(const Optimization::Settings& settings) {
  return getNewtonSystem_(settings, NewtonSystemType::Augmented);
}

NewtonSystem getNormalEquation(const Optimization::Settings& settings) {
  return getNewtonSystem_(settings, NewtonSystemType::Normal);
}

NewtonSystemTriplet getNewtonSystems(const Optimization::Settings& settings) {
  auto [full, augmented, normal] = getNewtonSystems_(settings);
  return {full, augmented, normal};
}

// Alternatively, use embind for more direct JS-to-C++ bindings
EMSCRIPTEN_BINDINGS(symbolic_optimization_module) {
  class_<Optimization::VariableNames>("VariableNames")
      .constructor<>()
      .property("x", &Optimization::VariableNames::x)
      .property("A_eq", &Optimization::VariableNames::A_eq)
      .property("A_ineq", &Optimization::VariableNames::A_ineq)
      .property("s_A", &Optimization::VariableNames::s_A)
      .property("s_Al", &Optimization::VariableNames::s_Al)
      .property("s_Au", &Optimization::VariableNames::s_Au)
      .property("s_xl", &Optimization::VariableNames::s_xl)
      .property("s_xu", &Optimization::VariableNames::s_xu)
      .property("l_A", &Optimization::VariableNames::l_A)
      .property("u_A", &Optimization::VariableNames::u_A)
      .property("l_x", &Optimization::VariableNames::l_x)
      .property("u_x", &Optimization::VariableNames::u_x);

  enum_<Optimization::Bounds>("Bounds")
      .value("None", Optimization::Bounds::None)
      .value("Lower", Optimization::Bounds::Lower)
      .value("Upper", Optimization::Bounds::Upper)
      .value("Both", Optimization::Bounds::Both);

  enum_<Optimization::InequalityHandling>("InequalityHandling")
      .value("Slacks", Optimization::InequalityHandling::Slacks)
      .value("SimpleSlacks", Optimization::InequalityHandling::SimpleSlacks);

  enum_<Optimization::EqualityHandling>("EqualityHandling")
      .value("None", Optimization::EqualityHandling::None)
      .value("Slacks", Optimization::EqualityHandling::Slacks)
      .value("SimpleSlacks", Optimization::EqualityHandling::SimpleSlacks)
      .value("PenaltyFunction",
             Optimization::EqualityHandling::PenaltyFunction);

  class_<Optimization::Settings>("Settings")
      .constructor<>()
      .property("inequalities", &Optimization::Settings::inequalities)
      .property("variableBounds", &Optimization::Settings::variableBounds)
      .property("equalities", &Optimization::Settings::equalities)
      .property("equalityHandling", &Optimization::Settings::equalityHandling)
      .property("inequalityHandling",
                &Optimization::Settings::inequalityHandling);

  value_object<NewtonSystem>("NewtonSystem")
      .field("lhs", &NewtonSystem::lhs)
      .field("rhs", &NewtonSystem::rhs)
      .field("rhsShorthand", &NewtonSystem::rhsShorthand)
      .field("variables", &NewtonSystem::variables)
      .field("variableDefinitions", &NewtonSystem::variableDefinitions);

  value_object<NewtonSystemTriplet>("NewtonSystemTriplet")
      .field("full", &NewtonSystemTriplet::full)
      .field("augmented", &NewtonSystemTriplet::augmented)
      .field("normal", &NewtonSystemTriplet::normal);

  // Register free functions
  function("getLagrangian", &getLagrangian);
  function("getFirstOrderOptimalityConditions",
           &getFirstOrderOptimalityConditions);
  function("getNewtonSystem", &getNewtonSystem);
  function("getAugmentedSystem", &getAugmentedSystem);
  function("getNormalEquation", &getNormalEquation);
  function("getNewtonSystems", &getNewtonSystems);
}