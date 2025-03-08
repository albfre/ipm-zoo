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
  std::vector<std::pair<Expression::Expr, Expression::Expr>>
      variableDefinitions;
};

namespace {
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

NewtonSystem formatNewtonSystemStrings_(const auto& lhs, const auto& rhs,
                                        const auto& variables) {
  const auto unity = Expression::ExprFactory::number(1);
  const auto I = Expression::ExprFactory::namedConstant("I");
  std::string lhsStr = "";
  const auto condensed = true;
  for (const auto& row : lhs) {
    for (size_t i = 0; i < row.size(); ++i) {
      lhsStr += row[i].replaceSubexpression(unity, I).toString(condensed) +
                (i + 1 == row.size() ? "" : " & ");
    }
    lhsStr += " \\\\ ";
  }

  std::string rhsStr = "";
  const auto breaks = std::set<std::string>{" + ", " - "};
  for (const auto& row : rhs) {
    auto rowStr = row.toString(condensed);
    int numParenthesis = 0;
    size_t numTerms = 0;
    size_t lastBreak = 0;
    for (size_t i = 0; i + 2 < rowStr.length(); ++i) {
      const auto c = rowStr[i];
      numParenthesis += c == '(' ? 1 : 0;
      numParenthesis -= c == ')' ? 1 : 0;
      numTerms += c == '+' ? 1 : 0;
      numTerms += c == '-' ? 1 : 0;
      const auto sub = rowStr.substr(i, 3);
      if (numParenthesis == 0 && breaks.contains(sub) &&
          numTerms > lastBreak + 3) {
        rowStr.insert(i, " \\\\\n \\qquad");
        lastBreak = numTerms;
      }
    }

    rhsStr += rowStr + " \\\\ ";
  }

  const auto rhsShorthand = Optimization::getShorthandRhs(variables);
  std::string rhsShorthandStr = "";
  for (const auto& row : rhsShorthand) {
    rhsShorthandStr += row.toString(condensed) + " \\\\ ";
  }

  std::string variablesStr = "";
  for (size_t i = 0; i < variables.size(); ++i) {
    variablesStr += "\\Delta " + variables[i].toString(condensed) +
                    (i + 1 == variables.size() ? "" : " \\\\ ");
  }

  return {lhsStr, rhsStr, rhsShorthandStr, variablesStr};
}

enum class NewtonSystemType { Full, Augmented, Normal };

NewtonSystem getNewtonSystem_(const Optimization::Settings& settings,
                              const NewtonSystemType type) {
  auto begin = std::chrono::high_resolution_clock::now();
  const auto variableNames = Optimization::VariableNames();
  auto [lagrangian, variables] =
      Optimization::getLagrangian(variableNames, settings);
  auto [lhs, rhs] = Optimization::getNewtonSystem(lagrangian, variables);
  auto end = std::chrono::high_resolution_clock::now();
  std::cout << "time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                     begin)
                   .count()
            << std::endl;

  auto i = lhs.size();
  if (type != NewtonSystemType::Full) {
    i = settings.equalities &&
                settings.equalityHandling ==
                    Optimization::EqualityHandling::IndefiniteFactorization
            ? 2
            : 1;
    rhs = Optimization::getShorthandRhs(variables);
  }
  if (type == NewtonSystemType::Augmented) {
    const auto zero = Expression::ExprFactory::number(0);
    const auto unity = Expression::ExprFactory::number(1);
    const auto negUnity =
        Expression::ExprFactory::negate(Expression::ExprFactory::number(1));
    const auto reducibles = std::set{zero, unity, negUnity};
    while (i < lhs.size() && !reducibles.contains(lhs.at(0).at(i))) {
      ++i;
    }
  }

  std::vector<std::pair<Expression::Expr, Expression::Expr>>
      variableDefinitions;
  while (lhs.size() > i) {
    variableDefinitions.push_back(
        Optimization::gaussianElimination(lhs, rhs, lhs.size() - 1, variables));
    variables.pop_back();
  }

  auto ns = formatNewtonSystemStrings_(lhs, rhs, variables);
  ns.variableDefinitions = std::move(variableDefinitions);
  return ns;
}

}  // namespace

std::string getLagrangian(const Optimization::Settings& settings) {
  const auto variableNames = Optimization::VariableNames();
  const auto [lagrangian, variables] =
      Optimization::getLagrangian(variableNames, settings);

  const auto condensed = true;
  auto str = "& " + lagrangian.toString(condensed);
  auto pos = 0;
  while (pos != std::string::npos) {
    pos = str.find("\\lambda", pos + 1);
    if (pos != std::string::npos && pos > 2) {
      str.insert(pos - 2, "\\\\\n & ");
    }
    for (size_t i = 0; i < 3; ++i) {
      pos = pos == std::string::npos ? std::string::npos
                                     : str.find("\\lambda", pos + 1);
    }
  }

  pos = str.find("- \\mu");
  if (pos != std::string::npos) {
    str.insert(pos, "\\\\\n & ");
  }
  return str;
}

std::string getFirstOrderOptimalityConditions(
    const Optimization::Settings& settings) {
  const auto variableNames = Optimization::VariableNames();
  const auto [lagrangian, variables] =
      Optimization::getLagrangian(variableNames, settings);
  const auto firstOrder =
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
      .value("IndefiniteFactorization",
             Optimization::EqualityHandling::IndefiniteFactorization)
      .value("Slacks", Optimization::EqualityHandling::Slacks)
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

  // Register free functions
  function("getLagrangian", &getLagrangian);
  function("getFirstOrderOptimalityConditions",
           &getFirstOrderOptimalityConditions);
  value_object<NewtonSystem>("NewtonSystem")
      .field("lhs", &NewtonSystem::lhs)
      .field("rhs", &NewtonSystem::rhs)
      .field("rhsShorthand", &NewtonSystem::rhsShorthand)
      .field("variables", &NewtonSystem::variables)
      .field("variableDefinitions", &NewtonSystem::variableDefinitions);
  function("getNewtonSystem", &getNewtonSystem);
  function("getAugmentedSystem", &getAugmentedSystem);
  function("getNormalEquation", &getNormalEquation);
}