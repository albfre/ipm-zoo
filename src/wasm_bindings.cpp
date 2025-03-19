#include <emscripten.h>
#include <emscripten/bind.h>

#include <chrono>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "Expression.h"
#include "Optimization.h"
#include "Timer.h"

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

NewtonSystem formatNewtonSystemStrings_(const auto& lhs, const auto& rhs,
                                        const auto& variables,
                                        const auto& variableDefinitions,
                                        const auto& variableNames) {
  using namespace Expression::ExprFactory;
  const auto unity = number(1);
  const auto I = namedVector("I");
  const auto mu = namedScalar("\\mu");
  const auto muI = product({mu, I});
  const auto delta = namedScalar(variableNames.delta_eq);
  const auto delta2 = namedScalar(variableNames.delta_eq + "^2");
  const auto deltaI = product({delta, I});
  const auto deltaI2 = product({deltaI, deltaI}).simplify();
  const auto delta2I = product({delta2, I}).simplify();
  std::string lhsStr = "";
  const auto condensed = true;
  for (const auto& row : lhs) {
    for (size_t i = 0; i < row.size(); ++i) {
      auto rowExpr = row[i];
      rowExpr = rowExpr.replaceSubexpression(unity, I);
      rowExpr = rowExpr.replaceSubexpression(mu, muI);
      rowExpr = rowExpr.replaceSubexpression(delta, deltaI).simplify();
      rowExpr = rowExpr.replaceSubexpression(deltaI2, delta2I);
      lhsStr +=
          rowExpr.toString(condensed) + (i + 1 == row.size() ? "" : " & ");
    }
    lhsStr += " \\\\\n ";
  }

  std::string rhsStr = "";
  for (const auto& row : rhs) {
    auto rowStr = row.toString(condensed);
    rhsStr += rowStr + " \\\\\n ";
  }

  const auto rhsShorthand = Optimization::getShorthandRhs(variables);
  std::string rhsShorthandStr = "";
  for (const auto& row : rhsShorthand) {
    rhsShorthandStr += row.toString(condensed) + " \\\\\n ";
  }

  std::string variablesStr = "";
  for (size_t i = 0; i < variables.size(); ++i) {
    variablesStr += "\\Delta " + variables[i].toString(condensed) +
                    (i + 1 == variables.size() ? "\n" : " \\\\\n ");
  }

  std::string definitionsStr = "";
  for (size_t i = variableDefinitions.size(); i > 0; --i) {
    definitionsStr +=
        variableDefinitions.at(i - 1).first.toString(condensed) +
        " &= " + variableDefinitions.at(i - 1).second.toString(condensed) +
        (i - 1 == 0 ? "\n" : " \\\\\n ");
  }

  return {lhsStr, rhsStr, rhsShorthandStr, variablesStr, definitionsStr};
}

enum class NewtonSystemType { Full, Augmented, Normal };

std::tuple<NewtonSystem, NewtonSystem, NewtonSystem> getNewtonSystems_(
    const Optimization::Settings& settingsIn,
    const Optimization::VariableNames& variableNames) {
  Timer timer;
  timer.start("getNewtonSystems");
  const auto settings = changePenaltyToPenaltyWithExtraVariables(settingsIn);
  timer.start("getLagrangian");
  auto [lagrangian, variables] =
      Optimization::getLagrangian(variableNames, settings);
  timer.stop("getLagrangian");
  timer.start("getNewtonSystem");
  auto [lhs, rhs] = Optimization::getNewtonSystem(lagrangian, variables);
  timer.stop("getNewtonSystem");

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

  auto newtonSystem = formatNewtonSystemStrings_(
      lhs, rhs, variables, variableDefinitions, variableNames);

  rhs = Optimization::getShorthandRhs(variables);
  while (lhs.size() > augmentedSize) {
    auto deltaVariable = Expression::ExprFactory::variable(
        "\\Delta " + variables.at(lhs.size() - 1).toString());
    timer.start("deltaDefinition");
    auto deltaDefinition =
        Optimization::deltaDefinition(lhs, rhs, variables, lhs.size() - 1);
    timer.stop("deltaDefinition");
    variableDefinitions.push_back({deltaVariable, deltaDefinition});
    timer.start("gaussianElimination");
    Optimization::gaussianElimination(lhs, rhs, lhs.size() - 1);
    timer.stop("gaussianElimination");
    variables.pop_back();
  }
  auto augmentedSystem = formatNewtonSystemStrings_(
      lhs, rhs, variables, variableDefinitions, variableNames);

  auto deltaVariable = Expression::ExprFactory::variable(
      "\\Delta " + variables.at(0).toString());
  auto deltaDefinition = Optimization::deltaDefinition(lhs, rhs, variables, 0);
  variableDefinitions.push_back({deltaVariable, deltaDefinition});
  Optimization::gaussianElimination(lhs, rhs, 0);
  variables.erase(variables.begin());

  timer.start("formatStrings");
  auto normalEquations = formatNewtonSystemStrings_(
      lhs, rhs, variables, variableDefinitions, variableNames);
  timer.stop("formatStrings");
  timer.stop("getNewtonSystems");
  timer.report();

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

NewtonSystemTriplet getNewtonSystems(const Optimization::Settings& settings) {
  const auto variableNames = Optimization::VariableNames();
  auto [full, augmented, normal] = getNewtonSystems_(settings, variableNames);
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
      .value("PenaltyFunction", Optimization::EqualityHandling::PenaltyFunction)
      .value("Regularization", Optimization::EqualityHandling::Regularization);

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
  function("getNewtonSystems", &getNewtonSystems);
}