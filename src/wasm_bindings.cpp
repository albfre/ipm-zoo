#include <emscripten.h>
#include <emscripten/bind.h>

#include <chrono>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "Expr.h"
#include "ExprFactory.h"
#include "SymbolicOptimization.h"
#include "Timer.h"

using namespace emscripten;

struct NewtonSystemStr {
  std::string lhs;
  std::string rhs;
  std::string rhsShorthand;
  std::string variables;
  std::string deltaDefinitions;
};

struct NewtonSystemTriplet {
  NewtonSystemStr full;
  NewtonSystemStr augmented;
  NewtonSystemStr normal;
};

namespace {
SymbolicOptimization::Settings changePenaltyToPenaltyWithExtraVariables(
    SymbolicOptimization::Settings settings) {
  if (settings.equalityHandling ==
      SymbolicOptimization::EqualityHandling::PenaltyFunction) {
    settings.equalityHandling = SymbolicOptimization::EqualityHandling::
        PenaltyFunctionWithExtraVariable;
  }
  return settings;
}

NewtonSystemStr formatNewtonSystemStrings_(const auto& newtonSystem,
                                           const auto& variableNames) {
  const auto& [lhs, rhs, variables, deltaDefinitions] = newtonSystem;
  using EF = Expression::ExprFactory;
  const auto unity = EF::number(1);
  const auto I = EF::namedVector("I");
  const auto mu = EF::namedScalar("\\mu");
  const auto muI = EF::product({mu, I});
  const auto delta = EF::namedScalar(variableNames.delta_eq);
  const auto delta2 = EF::namedScalar(variableNames.delta_eq + "^2");
  const auto deltaI = EF::product({delta, I});
  const auto deltaI2 = EF::product({deltaI, deltaI})->simplify();
  const auto delta2I = EF::product({delta2, I})->simplify();
  std::string lhsStr = "";
  const auto condensed = true;
  for (const auto& row : lhs) {
    for (size_t i = 0; i < row.size(); ++i) {
      auto rowExpr = row[i];
      rowExpr = rowExpr->replaceSubexpression(unity, I);
      rowExpr = rowExpr->replaceSubexpression(mu, muI);
      rowExpr = rowExpr->replaceSubexpression(delta, deltaI)->simplify();
      rowExpr = rowExpr->replaceSubexpression(deltaI2, delta2I);
      lhsStr +=
          rowExpr->toString(condensed) + (i + 1 == row.size() ? "" : " & ");
    }
    lhsStr += " \\\\\n ";
  }

  std::string rhsStr = "";
  for (const auto& row : rhs) {
    auto rowStr = row->toString(condensed);
    rhsStr += rowStr + " \\\\\n ";
  }

  const auto rhsShorthand = SymbolicOptimization::getShorthandRhs(variables);
  std::string rhsShorthandStr = "";
  for (const auto& row : rhsShorthand) {
    rhsShorthandStr += row->toString(condensed) + " \\\\\n ";
  }

  std::string variablesStr = "";
  for (size_t i = 0; i < variables.size(); ++i) {
    variablesStr += "\\Delta " + variables[i]->toString(condensed) +
                    (i + 1 == variables.size() ? "\n" : " \\\\\n ");
  }

  std::string definitionsStr = "";
  for (size_t i = deltaDefinitions.size(); i > 0; --i) {
    definitionsStr +=
        deltaDefinitions.at(i - 1).first->toString(condensed) +
        " &= " + deltaDefinitions.at(i - 1).second->toString(condensed) +
        (i - 1 == 0 ? "\n" : " \\\\\n ");
  }

  return {lhsStr, rhsStr, rhsShorthandStr, variablesStr, definitionsStr};
}

enum class NewtonSystemType { Full, Augmented, Normal };

std::tuple<NewtonSystemStr, NewtonSystemStr, NewtonSystemStr> getNewtonSystems_(
    const SymbolicOptimization::Settings& settingsIn,
    const SymbolicOptimization::VariableNames& variableNames) {
  Util::Timer timer;
  timer.start("getNewtonSystems");

  // Get correct Lagrangian
  const auto settings = changePenaltyToPenaltyWithExtraVariables(settingsIn);
  auto [lagrangian, variables] =
      SymbolicOptimization::getLagrangian(variableNames, settings);

  // Full Newton system
  auto newtonSystem =
      SymbolicOptimization::getNewtonSystem(lagrangian, variables);
  auto newtonSystemStr =
      formatNewtonSystemStrings_(newtonSystem, variableNames);

  // Use shorthand for right-hand side
  newtonSystem.rhs = SymbolicOptimization::getShorthandRhs(variables);

  // Augmented system
  auto augmentedSystem = SymbolicOptimization::getAugmentedSystem(newtonSystem);
  auto augmentedSystemStr =
      formatNewtonSystemStrings_(augmentedSystem, variableNames);

  // Normal equations
  auto normalEquations =
      SymbolicOptimization::getNormalEquations(augmentedSystem);
  auto normalEquationsStr =
      formatNewtonSystemStrings_(normalEquations, variableNames);

  timer.stop("getNewtonSystems");
  timer.report();

  return {newtonSystemStr, augmentedSystemStr, normalEquationsStr};
}
}  // namespace

std::string getLagrangian(const SymbolicOptimization::Settings& settings) {
  const auto variableNames = SymbolicOptimization::VariableNames();
  const auto [lagrangian, variables] =
      SymbolicOptimization::getLagrangian(variableNames, settings);

  const auto condensed = true;
  auto str = "& " + lagrangian->toString(condensed);
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
    const SymbolicOptimization::Settings& settingsIn) {
  const auto variableNames = SymbolicOptimization::VariableNames();
  const auto settings = changePenaltyToPenaltyWithExtraVariables(settingsIn);
  const auto [lagrangian, variables] =
      SymbolicOptimization::getLagrangian(variableNames, settings);
  auto firstOrder = SymbolicOptimization::getFirstOrderOptimalityConditions(
      lagrangian, variables);
  const auto condensed = true;
  std::string str = "";
  for (const auto& c : firstOrder) {
    str += c->toString(condensed) + " &= 0 \\\\";
  }
  return str;
}

NewtonSystemTriplet getNewtonSystems(
    const SymbolicOptimization::Settings& settings) {
  const auto variableNames = SymbolicOptimization::VariableNames();
  auto [full, augmented, normal] = getNewtonSystems_(settings, variableNames);
  return {full, augmented, normal};
}

// Alternatively, use embind for more direct JS-to-C++ bindings
EMSCRIPTEN_BINDINGS(symbolic_optimization_module) {
  class_<SymbolicOptimization::VariableNames>("VariableNames")
      .constructor<>()
      .property("x", &SymbolicOptimization::VariableNames::x)
      .property("A_eq", &SymbolicOptimization::VariableNames::A_eq)
      .property("A_ineq", &SymbolicOptimization::VariableNames::A_ineq)
      .property("s_A", &SymbolicOptimization::VariableNames::s_A)
      .property("s_Al", &SymbolicOptimization::VariableNames::s_Al)
      .property("s_Au", &SymbolicOptimization::VariableNames::s_Au)
      .property("s_xl", &SymbolicOptimization::VariableNames::s_xl)
      .property("s_xu", &SymbolicOptimization::VariableNames::s_xu)
      .property("l_A", &SymbolicOptimization::VariableNames::l_A)
      .property("u_A", &SymbolicOptimization::VariableNames::u_A)
      .property("l_x", &SymbolicOptimization::VariableNames::l_x)
      .property("u_x", &SymbolicOptimization::VariableNames::u_x);

  enum_<SymbolicOptimization::Bounds>("Bounds")
      .value("None", SymbolicOptimization::Bounds::None)
      .value("Lower", SymbolicOptimization::Bounds::Lower)
      .value("Upper", SymbolicOptimization::Bounds::Upper)
      .value("Both", SymbolicOptimization::Bounds::Both);

  enum_<SymbolicOptimization::InequalityHandling>("InequalityHandling")
      .value("Slacks", SymbolicOptimization::InequalityHandling::Slacks)
      .value("SimpleSlacks",
             SymbolicOptimization::InequalityHandling::SimpleSlacks);

  enum_<SymbolicOptimization::EqualityHandling>("EqualityHandling")
      .value("None", SymbolicOptimization::EqualityHandling::None)
      .value("Slacks", SymbolicOptimization::EqualityHandling::Slacks)
      .value("SimpleSlacks",
             SymbolicOptimization::EqualityHandling::SimpleSlacks)
      .value("PenaltyFunction",
             SymbolicOptimization::EqualityHandling::PenaltyFunction)
      .value("Regularization",
             SymbolicOptimization::EqualityHandling::Regularization);

  class_<SymbolicOptimization::Settings>("Settings")
      .constructor<>()
      .property("inequalities", &SymbolicOptimization::Settings::inequalities)
      .property("variableBounds",
                &SymbolicOptimization::Settings::variableBounds)
      .property("equalities", &SymbolicOptimization::Settings::equalities)
      .property("equalityHandling",
                &SymbolicOptimization::Settings::equalityHandling)
      .property("inequalityHandling",
                &SymbolicOptimization::Settings::inequalityHandling);

  value_object<NewtonSystemStr>("NewtonSystem")
      .field("lhs", &NewtonSystemStr::lhs)
      .field("rhs", &NewtonSystemStr::rhs)
      .field("rhsShorthand", &NewtonSystemStr::rhsShorthand)
      .field("variables", &NewtonSystemStr::variables)
      .field("deltaDefinitions", &NewtonSystemStr::deltaDefinitions);

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