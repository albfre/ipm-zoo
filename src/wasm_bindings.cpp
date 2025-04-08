#include <emscripten.h>
#include <emscripten/bind.h>

#include <chrono>
#include <iostream>
#include <memory>
#include <ranges>
#include <string>
#include <vector>

#include "Expr.h"
#include "ExprFactory.h"
#include "SymbolicOptimization.h"
#include "Utils/Timer.h"

using namespace emscripten;

struct NewtonSystemStr {
  std::string lhs;
  std::string rhs;
  std::string rhs_shorthand;
  std::string variables;
  std::string delta_definitions;
};

struct NewtonSystemTriplet {
  NewtonSystemStr full;
  NewtonSystemStr augmented;
  NewtonSystemStr normal;
};

namespace {
NewtonSystemStr format_newton_system_strings_(const auto& newton_system,
                                              const auto& variable_names) {
  const auto& [lhs, rhs, variables, delta_definitions] = newton_system;
  using EF = Expression::ExprFactory;
  const auto unity = EF::number(1);
  const auto I = EF::named_vector("I");
  const auto mu = EF::named_scalar("\\mu");
  const auto muI = EF::product({mu, I});
  const auto delta = EF::named_scalar(variable_names.delta_eq);
  const auto delta2 = EF::named_scalar(variable_names.delta_eq + "^2");
  const auto deltaI = EF::product({delta, I});
  const auto deltaI2 = EF::product({deltaI, deltaI})->simplify();
  const auto delta2I = EF::product({delta2, I})->simplify();
  std::string lhs_str = "";
  const auto condensed = true;
  for (const auto& row : lhs) {
    for (size_t i = 0; i < row.size(); ++i) {
      auto row_expr = row[i];
      row_expr = row_expr->replace_subexpression(unity, I);
      row_expr = row_expr->replace_subexpression(mu, muI);
      row_expr = row_expr->replace_subexpression(delta, deltaI)->simplify();
      row_expr = row_expr->replace_subexpression(deltaI2, delta2I);
      lhs_str +=
          row_expr->to_string(condensed) + (i + 1 == row.size() ? "" : " & ");
    }
    lhs_str += " \\\\\n ";
  }

  std::string rhs_str = "";
  for (const auto& row : rhs) {
    auto row_str = row->to_string(condensed);
    rhs_str += row_str + " \\\\\n ";
  }

  const auto rhs_shorthand =
      SymbolicOptimization::get_shorthand_rhs(newton_system).shorthand_rhs;
  std::string rhs_shorthand_str = "";
  for (const auto& row : rhs_shorthand) {
    rhs_shorthand_str += row->to_string(condensed) + " \\\\\n ";
  }

  std::string variables_str = "";
  for (size_t i = 0; i < variables.size(); ++i) {
    variables_str += "\\Delta " + variables[i]->to_string(condensed) +
                     (i + 1 == variables.size() ? "\n" : " \\\\\n ");
  }

  std::string definitions_str = "";
  for (size_t i = delta_definitions.size(); i > 0; --i) {
    definitions_str +=
        delta_definitions.at(i - 1).first->to_string(condensed) +
        " &= " + delta_definitions.at(i - 1).second->to_string(condensed) +
        (i - 1 == 0 ? "\n" : " \\\\\n ");
  }

  return {lhs_str, rhs_str, rhs_shorthand_str, variables_str, definitions_str};
}
}  // namespace

std::string get_optimization_problem(
    const SymbolicOptimization::Settings& settings,
    const SymbolicOptimization::OptimizationProblemType&
        optimization_problem_type) {
  const auto variable_names = SymbolicOptimization::VariableNames();
  const auto optimization_problem =
      SymbolicOptimization::get_optimization_problem(settings, variable_names,
                                                     optimization_problem_type);

  const auto condensed = true;
  auto str = std::string("\\text{minimize} \\quad & ");
  str += optimization_problem.objective->to_string(condensed);
  auto pos = 0;
  const auto add_newlines = [&](const auto& term) {
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
  add_newlines("\\lambda");
  add_newlines("- \\mu");

  str += " \\\\\n ";
  str += "\\text{subject to} \\quad";
  for (const auto& inequality : optimization_problem.inequalities) {
    str += " & ";
    if (inequality.lower_bound && inequality.upper_bound) {
      str += inequality.lower_bound->to_string(condensed) + " \\leq " +
             inequality.expr->to_string(condensed) + " \\leq " +
             inequality.upper_bound->to_string(condensed) + " \\\\\n";
    } else if (inequality.lower_bound) {
      str += inequality.expr->to_string(condensed) + " \\geq " +
             inequality.lower_bound->to_string(condensed) + " \\\\\n";
    } else if (inequality.upper_bound) {
      str += inequality.expr->to_string(condensed) + " \\leq " +
             inequality.upper_bound->to_string(condensed) + " \\\\\n";
    }
  }

  for (const auto& equality : optimization_problem.equalities) {
    str += " & ";
    str += equality.expr->to_string(condensed) + " = " +
           equality.rhs->to_string(condensed) + " \\\\\n";
  }

  for (const auto& bound : optimization_problem.variable_bounds) {
    str += " & ";
    if (bound.lower_bound && bound.upper_bound) {
      str += bound.lower_bound->to_string(condensed) + " \\leq " +
             bound.expr->to_string(condensed) + " \\leq " +
             bound.upper_bound->to_string(condensed) + " \\\\\n";
    } else if (bound.lower_bound) {
      str += bound.expr->to_string(condensed) + " \\geq " +
             bound.lower_bound->to_string(condensed) + " \\\\\n";
    } else if (bound.upper_bound) {
      str += bound.expr->to_string(condensed) + " \\leq " +
             bound.upper_bound->to_string(condensed) + " \\\\\n";
    }
  }

  if (!optimization_problem.nonnegative_slacks.empty()) {
    str += " & ";
    str +=
        optimization_problem.nonnegative_slacks.front()->to_string(condensed);
    for (const auto& slack :
         optimization_problem.nonnegative_slacks | std::views::drop(1)) {
      str += ", " + slack->to_string(condensed);
    }
    str += " \\geq 0 \\\\\n";
  }
  return str;
}

std::string get_lagrangian(const SymbolicOptimization::Settings& settings) {
  const auto variable_names = SymbolicOptimization::VariableNames();
  const auto problem = SymbolicOptimization::get_optimization_problem(
      settings, variable_names,
      SymbolicOptimization::OptimizationProblemType::SlackedWithBarriers);
  const auto lagrangian = SymbolicOptimization::get_lagrangian(problem);

  const auto condensed = true;
  auto str = "& " + lagrangian->to_string(condensed);
  auto pos = 0;
  const auto add_newlines = [&](const auto& term) {
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
  add_newlines("\\lambda");
  add_newlines("- \\mu");
  return str;
}

std::string get_first_order_optimality_conditions(
    const SymbolicOptimization::Settings& settings) {
  const auto variable_names = SymbolicOptimization::VariableNames();
  auto [first_order, _] =
      SymbolicOptimization::get_first_order_optimality_conditions(
          settings, variable_names);
  const auto condensed = true;
  std::string str = "";
  for (const auto& c : first_order) {
    str += c->to_string(condensed) + " &= 0 \\\\";
  }
  return str;
}

NewtonSystemTriplet get_newton_systems(
    const SymbolicOptimization::Settings& settings) {
  const auto variable_names = SymbolicOptimization::VariableNames();
  Utils::Timer timer;
  timer.start("get_newton_systems");

  // Full Newton system
  auto newton_system =
      SymbolicOptimization::get_newton_system(settings, variable_names);
  auto newton_system_str =
      format_newton_system_strings_(newton_system, variable_names);

  // Use shorthand for right-hand side
  newton_system.rhs =
      SymbolicOptimization::get_shorthand_rhs(newton_system).shorthand_rhs;

  // Augmented system
  auto augmented_system =
      SymbolicOptimization::get_augmented_system(newton_system);
  auto augmented_system_str =
      format_newton_system_strings_(augmented_system, variable_names);

  // Normal equations
  auto normal_equations =
      SymbolicOptimization::get_normal_equations(augmented_system);
  auto normal_equations_str =
      format_newton_system_strings_(normal_equations, variable_names);

  timer.stop("get_newton_systems");
  timer.report();

  return {newton_system_str, augmented_system_str, normal_equations_str};
}

// Alternatively, use embind for more direct JS-to-C++ bindings
EMSCRIPTEN_BINDINGS(symbolic_optimization_module) {
  class_<SymbolicOptimization::VariableNames>("VariableNames")
      .constructor<>()
      .property("x", &SymbolicOptimization::VariableNames::x)
      .property("A_eq", &SymbolicOptimization::VariableNames::A_eq)
      .property("A_ineq", &SymbolicOptimization::VariableNames::A_ineq)
      .property("s_A_ineq", &SymbolicOptimization::VariableNames::s_A_ineq)
      .property("s_A_ineq_l", &SymbolicOptimization::VariableNames::s_A_ineq_l)
      .property("s_A_ineq_u", &SymbolicOptimization::VariableNames::s_A_ineq_u)
      .property("s_x_l", &SymbolicOptimization::VariableNames::s_x_l)
      .property("s_x_u", &SymbolicOptimization::VariableNames::s_x_u)
      .property("l_A_ineq", &SymbolicOptimization::VariableNames::l_A_ineq)
      .property("u_A_ineq", &SymbolicOptimization::VariableNames::u_A_ineq)
      .property("l_x", &SymbolicOptimization::VariableNames::l_x)
      .property("u_x", &SymbolicOptimization::VariableNames::u_x);

  enum_<SymbolicOptimization::Bounds>("Bounds")
      .value("None", SymbolicOptimization::Bounds::None)
      .value("Lower", SymbolicOptimization::Bounds::Lower)
      .value("Upper", SymbolicOptimization::Bounds::Upper)
      .value("Both", SymbolicOptimization::Bounds::Both);

  enum_<SymbolicOptimization::InequalityHandling>("InequalityHandling")
      .value("Slacks", SymbolicOptimization::InequalityHandling::Slacks)
      .value("SlackedSlacks",
             SymbolicOptimization::InequalityHandling::SlackedSlacks)
      .value("NaiveSlacks",
             SymbolicOptimization::InequalityHandling::NaiveSlacks);

  enum_<SymbolicOptimization::EqualityHandling>("EqualityHandling")
      .value("None", SymbolicOptimization::EqualityHandling::None)
      .value("Slacks", SymbolicOptimization::EqualityHandling::Slacks)
      .value("SlackedSlacks",
             SymbolicOptimization::EqualityHandling::SlackedSlacks)
      .value("NaiveSlacks", SymbolicOptimization::EqualityHandling::NaiveSlacks)
      .value("PenaltyFunction",
             SymbolicOptimization::EqualityHandling::PenaltyFunction)
      .value("Regularization",
             SymbolicOptimization::EqualityHandling::Regularization);

  class_<SymbolicOptimization::Settings>("Settings")
      .constructor<>()
      .property("inequalities", &SymbolicOptimization::Settings::inequalities)
      .property("variableBounds",
                &SymbolicOptimization::Settings::variable_bounds)
      .property("equalities", &SymbolicOptimization::Settings::equalities)
      .property("equalityHandling",
                &SymbolicOptimization::Settings::equality_handling)
      .property("inequalityHandling",
                &SymbolicOptimization::Settings::inequality_handling);

  enum_<SymbolicOptimization::OptimizationProblemType>(
      "OptimizationProblemType")
      .value("Original",
             SymbolicOptimization::OptimizationProblemType::Original)
      .value("Slacked", SymbolicOptimization::OptimizationProblemType::Slacked)
      .value(
          "SlackedWithBarriers",
          SymbolicOptimization::OptimizationProblemType::SlackedWithBarriers);

  value_object<NewtonSystemStr>("NewtonSystem")
      .field("lhs", &NewtonSystemStr::lhs)
      .field("rhs", &NewtonSystemStr::rhs)
      .field("rhsShorthand", &NewtonSystemStr::rhs_shorthand)
      .field("variables", &NewtonSystemStr::variables)
      .field("deltaDefinitions", &NewtonSystemStr::delta_definitions);

  value_object<NewtonSystemTriplet>("NewtonSystemTriplet")
      .field("full", &NewtonSystemTriplet::full)
      .field("augmented", &NewtonSystemTriplet::augmented)
      .field("normal", &NewtonSystemTriplet::normal);

  // Register free functions
  function("getOptimizationProblem", &get_optimization_problem);
  function("getLagrangian", &get_lagrangian);
  function("getFirstOrderOptimalityConditions",
           &get_first_order_optimality_conditions);
  function("getNewtonSystems", &get_newton_systems);
}