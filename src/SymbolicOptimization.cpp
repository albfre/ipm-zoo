#include "SymbolicOptimization.h"

#include <iostream>

#include "ExprFactory.h"
#include "SymbolicOperators.h"
#include "Utils/Assert.h"
#include "Utils/Helpers.h"

namespace SymbolicOptimization {
OptimizationExpressions get_optimization_expressions(
    const VariableNames& names) {
  using EF = Expression::ExprFactory;
  OptimizationExpressions oe;
  oe.Q = EF::symmetric_matrix(names.Q);
  oe.c = EF::named_vector(names.c);
  oe.A_ineq = EF::matrix(names.A_ineq);
  oe.A_eq = EF::matrix(names.A_eq);
  oe.b_eq = EF::matrix(names.b_eq);
  oe.p_eq = EF::variable(names.p_eq);
  oe.delta_eq = EF::named_scalar(names.delta_eq);
  oe.mu = EF::named_scalar("\\mu");
  oe.e_var = EF::named_vector("e_{" + names.x + "}");
  oe.e_ineq = EF::named_vector("e_{" + names.A_ineq + "}");
  oe.e_eq = EF::named_vector("e_{" + names.A_eq + "}");
  oe.x = EF::variable(names.x);
  oe.s_A_ineq = EF::variable(names.s_A_ineq);
  oe.s_A_ineq_l = EF::variable(names.s_A_ineq_l);
  oe.s_A_ineq_u = EF::variable(names.s_A_ineq_u);
  oe.s_x_l = EF::variable(names.s_x_l);
  oe.s_x_u = EF::variable(names.s_x_u);
  oe.s_A_eq = EF::variable(names.s_A_eq);
  oe.s_A_eq_l = EF::variable(names.s_A_eq_l);
  oe.s_A_eq_u = EF::variable(names.s_A_eq_u);
  oe.lambda_A_eq = EF::variable("\\lambda_{" + names.A_eq + "}");
  oe.lambda_sAeql = EF::variable("\\lambda_{" + names.s_A_eq_l + "}");
  oe.lambda_sAequ = EF::variable("\\lambda_{" + names.s_A_eq_u + "}");
  oe.lambda_A_ineq = EF::variable("\\lambda_{" + names.A_ineq + "}");
  oe.lambda_sAineql = EF::variable("\\lambda_{" + names.s_A_ineq_l + "}");
  oe.lambda_sAinequ = EF::variable("\\lambda_{" + names.s_A_ineq_u + "}");
  oe.lambda_sxl = EF::variable("\\lambda_{" + names.s_x_l + "}");
  oe.lambda_sxu = EF::variable("\\lambda_{" + names.s_x_u + "}");
  oe.l_A_ineq = EF::named_vector(names.l_A_ineq);
  oe.u_A_ineq = EF::named_vector(names.u_A_ineq);
  oe.l_x = EF::named_vector(names.l_x);
  oe.u_x = EF::named_vector(names.u_x);
  return oe;
}

OptimizationProblem get_optimization_problem(
    const Settings& settings, const VariableNames& names,
    const OptimizationProblemType& optimization_problem_type) {
  using EF = Expression::ExprFactory;
  using namespace Expression;
  const auto o = get_optimization_expressions(names);
  const auto half = EF::number(0.5);
  const auto zero = Expression::zero;
  const auto xQx = (half * EF::transpose(o.x) * o.Q * o.x)->simplify();
  const auto cx = EF::transpose(o.c) * o.x;
  const auto exT = EF::transpose(o.e_var);
  const auto eAT = EF::transpose(o.e_ineq);
  const auto eCT = EF::transpose(o.e_eq);
  const auto Ax = o.A_ineq * o.x;
  const auto Cx = o.A_eq * o.x;
  const auto CxMinusB = Cx - o.b_eq;
  const auto X = EF::diagonal_matrix(o.x);
  const auto S_ineq = EF::diagonal_matrix(o.s_A_ineq);
  const auto S_eq = EF::diagonal_matrix(o.s_A_eq);
  const auto L_x = EF::diagonal_matrix(o.l_x);
  const auto U_x = EF::diagonal_matrix(o.u_x);
  const auto L_A_ineq = EF::diagonal_matrix(o.l_A_ineq);
  const auto U_A_ineq = EF::diagonal_matrix(o.u_A_ineq);

  const auto slacked = std::set{OptimizationProblemType::Slacked,
                                OptimizationProblemType::SlackedWithBarriers}
                           .contains(optimization_problem_type);

  OptimizationProblem problem;

  problem.objective = xQx + cx;
  problem.variables = {o.x};

  const auto add_lower_inequalities =
      std::set{Bounds::Lower, Bounds::Both}.contains(settings.inequalities);
  const auto add_upper_inequalities =
      std::set{Bounds::Upper, Bounds::Both}.contains(settings.inequalities);

  if (add_lower_inequalities || add_upper_inequalities) {
    const auto lower_expr = [&](const auto& expr) -> ExprPtr {
      return add_lower_inequalities ? expr : nullptr;
    };
    const auto upper_expr = [&](const auto& expr) -> ExprPtr {
      return add_upper_inequalities ? expr : nullptr;
    };
    if (!slacked) {
      problem.inequalities.push_back(
          {Ax, lower_expr(o.l_A_ineq), upper_expr(o.u_A_ineq),
           lower_expr(-o.lambda_sAineql), upper_expr(o.lambda_sAinequ)});
    } else {
      switch (settings.inequality_handling) {
        case InequalityHandling::Slacks:
          problem.equalities.push_back(
              {Ax - o.s_A_ineq, zero, o.lambda_A_ineq});
          problem.variable_bounds.push_back(
              {o.s_A_ineq, lower_expr(o.l_A_ineq), upper_expr(o.u_A_ineq),
               lower_expr(o.lambda_sAineql), upper_expr(o.lambda_sAinequ)});
          problem.variables2.push_back(o.lambda_A_ineq);
          problem.variables3.push_back(o.s_A_ineq);
          break;
        case InequalityHandling::SlackedSlacks:
          problem.equalities.push_back(
              {Ax - o.s_A_ineq, zero, o.lambda_A_ineq});
          problem.variables2.push_back(o.lambda_A_ineq);
          problem.variables3.push_back(o.s_A_ineq);
          if (add_lower_inequalities) {
            problem.equalities.push_back(
                {o.s_A_ineq - o.s_A_ineq_l, o.l_A_ineq, -o.lambda_sAineql});
            problem.variables4.push_back(o.lambda_sAineql);
            problem.nonnegative_slacks.push_back(o.s_A_ineq_l);
          }
          if (add_upper_inequalities) {
            problem.equalities.push_back(
                {o.s_A_ineq + o.s_A_ineq_u, o.u_A_ineq, o.lambda_sAinequ});
            problem.variables4.push_back(o.lambda_sAinequ);
            problem.nonnegative_slacks.push_back(o.s_A_ineq_u);
          }
          break;
        case InequalityHandling::NaiveSlacks:
          if (add_lower_inequalities) {
            problem.equalities.push_back(
                {Ax - o.s_A_ineq_l, o.l_A_ineq, -o.lambda_sAineql});
            problem.variables2.push_back(o.lambda_sAineql);
            problem.nonnegative_slacks.push_back(o.s_A_ineq_l);
          }
          if (add_upper_inequalities) {
            problem.equalities.push_back(
                {Ax + o.s_A_ineq_u, o.u_A_ineq, o.lambda_sAinequ});
            problem.variables2.push_back(o.lambda_sAinequ);
            problem.nonnegative_slacks.push_back(o.s_A_ineq_u);
          }
          break;
        default:
          ASSERT(false);
      }
    }
  }

  if (settings.equalities) {
    if (!slacked || settings.equality_handling == EqualityHandling::None) {
      problem.equalities.push_back({Cx, o.b_eq, o.lambda_A_eq});
      problem.variables2.push_back(o.lambda_A_eq);
    } else {
      switch (settings.equality_handling) {
        case EqualityHandling::Slacks:
          problem.equalities.push_back({Cx - o.s_A_eq, zero, o.lambda_A_eq});
          problem.variable_bounds.push_back(
              {o.s_A_eq, o.b_eq, o.b_eq, -o.lambda_sAeql, o.lambda_sAequ});
          problem.variables3.push_back(o.s_A_eq);
          problem.variables2.push_back(o.lambda_sAeql);
          problem.variables2.push_back(o.lambda_sAequ);
          break;
        case EqualityHandling::SlackedSlacks:
          problem.equalities.push_back({Cx - o.s_A_eq, zero, o.lambda_A_eq});
          problem.equalities.push_back(
              {o.s_A_eq - o.s_A_eq_l, o.b_eq, -o.lambda_sAeql});
          problem.equalities.push_back(
              {o.s_A_eq + o.s_A_eq_u, o.b_eq, o.lambda_sAequ});
          problem.variables2.push_back(o.lambda_A_eq);
          problem.variables3.push_back(o.s_A_eq);
          problem.variables4.push_back(o.lambda_sAeql);
          problem.variables4.push_back(o.lambda_sAequ);
          problem.nonnegative_slacks.push_back(o.s_A_eq_l);
          problem.nonnegative_slacks.push_back(o.s_A_eq_u);
          break;
        case EqualityHandling::NaiveSlacks:
          problem.equalities.push_back(
              {Cx - o.s_A_eq_l, o.b_eq, -o.lambda_sAeql});
          problem.equalities.push_back(
              {Cx + o.s_A_eq_u, o.b_eq, o.lambda_sAequ});
          problem.variables2.push_back(o.lambda_sAeql);
          problem.variables2.push_back(o.lambda_sAequ);
          problem.nonnegative_slacks.push_back(o.s_A_eq_l);
          problem.nonnegative_slacks.push_back(o.s_A_eq_u);
          break;
        case EqualityHandling::PenaltyFunction: {
          const auto muTerm = half * EF::invert(o.mu);
          problem.objective =
              problem.objective + muTerm * EF::transpose(CxMinusB) * CxMinusB;
          break;
        }
        case EqualityHandling::PenaltyFunctionWithExtraDual:
          problem.equalities.push_back(
              {CxMinusB - half * o.mu * o.lambda_A_eq, zero, o.lambda_A_eq});
          problem.variables2.push_back(o.lambda_A_eq);
          break;
        case EqualityHandling::Regularization: {
          problem.objective =
              problem.objective +
              (half * EF::transpose(o.p_eq) * o.p_eq)->simplify();
          problem.equalities.push_back(
              {CxMinusB + o.delta_eq * o.p_eq, zero, o.lambda_A_eq});
          problem.variables2.push_back(o.lambda_A_eq);
          problem.variables3.push_back(o.p_eq);
          break;
        }
        default:
          ASSERT(false);
      }
    }
  }

  const auto add_lower_bounds =
      std::set{Bounds::Lower, Bounds::Both}.contains(settings.variable_bounds);
  const auto add_upper_bounds =
      std::set{Bounds::Upper, Bounds::Both}.contains(settings.variable_bounds);

  if (add_lower_bounds || add_upper_bounds) {
    const auto lower_expr = [&](const auto& expr) -> ExprPtr {
      return add_lower_bounds ? expr : nullptr;
    };
    const auto upper_expr = [&](const auto& expr) -> ExprPtr {
      return add_upper_bounds ? expr : nullptr;
    };
    if (!slacked ||
        settings.inequality_handling == InequalityHandling::Slacks) {
      problem.variable_bounds.push_back(
          {o.x, lower_expr(o.l_x), upper_expr(o.u_x), lower_expr(o.lambda_sxl),
           upper_expr(o.lambda_sxu)});
    } else {
      switch (settings.inequality_handling) {
        case InequalityHandling::SlackedSlacks:
        case InequalityHandling::NaiveSlacks:
          if (add_lower_bounds) {
            problem.equalities.push_back({o.x - o.s_x_l, o.l_x, -o.lambda_sxl});
            problem.variables4.push_back(o.lambda_sxl);
            problem.nonnegative_slacks.push_back(o.s_x_l);
          }
          if (add_upper_bounds) {
            problem.equalities.push_back({o.x + o.s_x_u, o.u_x, o.lambda_sxu});
            problem.variables4.push_back(o.lambda_sxu);
            problem.nonnegative_slacks.push_back(o.s_x_u);
          }
          break;
        default:
          ASSERT(false);
      }
    }
  }

  if (optimization_problem_type ==
      OptimizationProblemType::SlackedWithBarriers) {
    // Convert all inequalities to log barriers
    ASSERT(problem.inequalities.size() <= 1);
    const auto e_ineq_set = std::set{o.s_A_ineq, o.s_A_ineq_l, o.s_A_ineq_u};
    const auto e_eq_set = std::set{o.s_A_eq, o.s_A_eq_l, o.s_A_eq_u};
    const auto e_var_set = std::set{o.x, o.s_x_l, o.s_x_u};
    const auto get_eT = [&](const auto& expr) -> Expression::ExprPtr {
      return e_var_set.contains(expr)    ? exT
             : e_ineq_set.contains(expr) ? eAT
             : e_eq_set.contains(expr)   ? eCT
                                         : nullptr;
    };
    for (const auto& bounds : {problem.inequalities, problem.variable_bounds}) {
      for (const auto& bound : bounds) {
        const auto& [expr, lower, upper, lower_dual, upper_dual] = bound;
        const auto eT = get_eT(expr);
        if (lower != nullptr) {
          problem.objective = problem.objective -
                              (o.mu * eT * EF::log(expr - lower))->simplify();
        }

        if (upper != nullptr) {
          problem.objective = problem.objective -
                              (o.mu * eT * EF::log(upper - expr))->simplify();
        }
      }
    }
    for (const auto& slack : problem.nonnegative_slacks) {
      const auto eT = get_eT(slack);
      problem.objective =
          problem.objective - (o.mu * eT * EF::log(slack))->simplify();
    }
    problem.inequalities.clear();
    problem.variable_bounds.clear();
  }
  return problem;
}

Expression::ExprPtr get_lagrangian(const OptimizationProblem& problem) {
  using EF = Expression::ExprFactory;
  auto terms = std::vector{problem.objective};

  std::vector<Expression::ExprPtr> lambda_definitions;
  for (const auto& bounds : {problem.inequalities, problem.variable_bounds}) {
    for (const auto& bound : bounds) {
      ASSERT(bound.lower_dual_variable != nullptr ||
             bound.upper_dual_variable != nullptr);
      if (bound.lower_bound) {
        ASSERT(bound.lower_dual_variable != nullptr);
        terms.push_back((-EF::transpose(bound.lower_dual_variable) *
                         (bound.expr - bound.lower_bound))
                            ->simplify());
      }
      if (bound.upper_bound) {
        ASSERT(bound.upper_dual_variable != nullptr);
        terms.push_back((-EF::transpose(bound.upper_dual_variable) *
                         (bound.upper_bound - bound.expr))
                            ->simplify());
      }
    }
  }
  for (const auto& equality : problem.equalities) {
    ASSERT(equality.dual_variable != nullptr);
    terms.push_back(
        (EF::transpose(equality.dual_variable) * (equality.expr - equality.rhs))
            ->simplify());
  }
  auto lagrangian = EF::sum(terms);
  return lagrangian;
}

std::pair<std::vector<Expression::ExprPtr>, std::vector<Expression::ExprPtr>>
get_first_order_optimality_conditions(Settings settings,
                                      const VariableNames& names) {
  using EF = Expression::ExprFactory;

  if (settings.equality_handling == EqualityHandling::PenaltyFunction) {
    settings.equality_handling = EqualityHandling::PenaltyFunctionWithExtraDual;
  }
  auto optimization_problem_type = OptimizationProblemType::SlackedWithBarriers;
  if (settings.inequality_handling == InequalityHandling::Slacks) {
    optimization_problem_type = OptimizationProblemType::Slacked;
  }
  auto problem =
      get_optimization_problem(settings, names, optimization_problem_type);
  auto lagrangian = get_lagrangian(problem);
  auto variables = problem.variables;
  for (const auto& v : {problem.variables2, problem.variables3,
                        problem.variables4, problem.nonnegative_slacks}) {
    variables.insert(variables.end(), v.begin(), v.end());
  }

  std::vector<Expression::ExprPtr> lambda_definitions;
  std::vector<Expression::ExprPtr> extra_variables;
  if (optimization_problem_type !=
      OptimizationProblemType::SlackedWithBarriers) {
    const auto o = get_optimization_expressions(names);
    for (const auto& bounds : {problem.inequalities, problem.variable_bounds}) {
      for (const auto& bound : bounds) {
        ASSERT(bound.lower_dual_variable != nullptr ||
               bound.upper_dual_variable != nullptr);
        const auto e = bound.expr == o.x ? o.e_var : o.e_ineq;
        if (bound.lower_bound) {
          ASSERT(bound.lower_dual_variable != nullptr);
          lambda_definitions.push_back(
              (EF::diagonal_matrix(bound.expr) -
               EF::diagonal_matrix(bound.lower_bound)) *
                  bound.lower_dual_variable -
              o.mu * e);
          extra_variables.push_back(bound.lower_dual_variable);
        }
        if (bound.upper_bound) {
          ASSERT(bound.upper_dual_variable != nullptr);
          lambda_definitions.push_back((EF::diagonal_matrix(bound.upper_bound) -
                                        EF::diagonal_matrix(bound.expr)) *
                                           bound.upper_dual_variable -
                                       o.mu * e);
          extra_variables.push_back(bound.upper_dual_variable);
        }
      }
    }
  }

  std::vector<Expression::ExprPtr> first_order;
  first_order.reserve(variables.size());
  for (const auto& v : variables) {
    auto diff = lagrangian->differentiate(v)->simplify();
    if (const auto invV = EF::invert(EF::diagonal_matrix(v));
        diff->contains_subexpression(invV)) {
      diff = (EF::diagonal_matrix(v) * diff)->simplify();
    }
    first_order.push_back(diff);
  }

  first_order.insert(first_order.end(), lambda_definitions.begin(),
                     lambda_definitions.end());
  variables.insert(variables.end(), extra_variables.begin(),
                   extra_variables.end());
  return {first_order, variables};
}

NewtonSystem get_newton_system(const Settings& settings,
                               const VariableNames& names) {
  using EF = Expression::ExprFactory;
  auto [first_order, variables] =
      get_first_order_optimality_conditions(settings, names);
  auto lhs = std::vector<std::vector<Expression::ExprPtr>>();
  auto rhs = std::vector<Expression::ExprPtr>();
  for (auto& c : first_order) {
    lhs.emplace_back();
    auto& row = lhs.back();
    for (auto& v : variables) {
      row.push_back(c->differentiate(v)->simplify());
    }
    for (auto& v : variables) {
      for (auto& r : row) {
        if (r == v) {
          r = EF::diagonal_matrix(v);
          break;
        }
        if (r == EF::negate(v)) {
          r = -EF::diagonal_matrix(v);
          break;
        }
      }
    }
    rhs.push_back(EF::negate(c)->simplify());
  }
  return {lhs, rhs, variables, {}};
}

namespace {
const auto get_augmented_system_size(
    const std::vector<std::vector<Expression::ExprPtr>>& lhs) {
  const auto zero = Expression::zero;
  const auto unity = Expression::unity;
  const auto negUnity = -unity;
  const auto reducibles = std::set{zero, unity, negUnity};
  size_t i = 0;
  while (i < lhs.size() && !reducibles.contains(lhs.at(0).at(i))) {
    ++i;
  }
  return i;
}

}  // namespace

NewtonSystem get_augmented_system(NewtonSystem newton_system) {
  auto& [lhs, rhs, variables, delta_definitions] = newton_system;
  const auto augmented_system_size = get_augmented_system_size(lhs);
  while (lhs.size() > augmented_system_size) {
    auto delta_variable = get_delta_variable(variables.at(lhs.size() - 1));
    auto delta_def = SymbolicOptimization::delta_definition(lhs, rhs, variables,
                                                            lhs.size() - 1);
    delta_definitions.push_back({delta_variable, delta_def});
    SymbolicOptimization::gaussian_elimination(lhs, rhs, lhs.size() - 1);
    variables.pop_back();
  }
  return std::move(newton_system);
}

NewtonSystem get_normal_equations(NewtonSystem newton_system) {
  newton_system = get_augmented_system(std::move(newton_system));
  auto& [lhs, rhs, variables, delta_definitions] = newton_system;
  auto delta_variable = Expression::ExprFactory::variable(
      "\\Delta " + variables.at(0)->to_string());
  auto delta_def =
      SymbolicOptimization::delta_definition(lhs, rhs, variables, 0);
  delta_definitions.push_back({delta_variable, delta_def});
  SymbolicOptimization::gaussian_elimination(lhs, rhs, 0);
  variables.erase(variables.begin());
  return std::move(newton_system);
}

ShorthandRhs get_shorthand_rhs(const NewtonSystem& newton_system) {
  using EF = Expression::ExprFactory;
  ShorthandRhs rhs;
  ASSERT(newton_system.variables.size() == newton_system.rhs.size());
  for (size_t i = 0; i < newton_system.variables.size(); ++i) {
    const auto& var = newton_system.variables.at(i);
    const auto vec = EF::named_vector("r_{" + var->to_string() + "}");
    rhs.shorthand_rhs.push_back(EF::negate(vec));
    rhs.vector_definitions.push_back(
        {vec, EF::negate(newton_system.rhs.at(i))->simplify()});
  }
  return rhs;
}

Expression::ExprPtr get_delta_variable(const Expression::ExprPtr& expr) {
  ASSERT(Expression::is<Expression::Variable>(expr));
  return Expression::ExprFactory::variable("\\Delta " + expr->to_string());
}

Expression::ExprPtr delta_definition(
    const std::vector<std::vector<Expression::ExprPtr>>& lhs,
    const std::vector<Expression::ExprPtr>& rhs,
    const std::vector<Expression::ExprPtr>& variables,
    const size_t source_row) {
  using EF = Expression::ExprFactory;

  ASSERT(lhs.size() == rhs.size());
  ASSERT(lhs.size() <= variables.size());
  ASSERT(source_row < lhs.size());

  const auto& lhs_source_row = lhs.at(source_row);
  const auto& source_expr = lhs_source_row.at(source_row);
  auto terms = std::vector<Expression::ExprPtr>();
  terms.reserve(lhs_source_row.size());
  ASSERT(lhs_source_row.size() <= variables.size());
  for (size_t i = 0; i < lhs_source_row.size(); ++i) {
    auto delta_variable = get_delta_variable(variables.at(i));
    terms.emplace_back(
        EF::product({lhs_source_row[i], std::move(delta_variable)}));
  }

  terms.erase(terms.begin() + source_row);
  auto sum = EF::sum(std::move(terms));

  return EF::product({EF::invert(source_expr),
                      EF::sum({rhs.at(source_row), EF::negate(sum)})})
      ->simplify();
}

void gaussian_elimination(std::vector<std::vector<Expression::ExprPtr>>& lhs,
                          std::vector<Expression::ExprPtr>& rhs,
                          const size_t source_row) {
  using EF = Expression::ExprFactory;
  using Expression::zero;
  using namespace Expression;

  ASSERT(lhs.size() == rhs.size());
  ASSERT(source_row < lhs.size());
  std::set<size_t> target_rows;
  for (size_t i = 0; i < lhs.size(); ++i) {
    if (i != source_row && lhs.at(i).at(source_row) != zero) {
      target_rows.insert(i);
    }
  }
  ASSERT(!target_rows.empty());
  for (const auto& target_row : target_rows) {
    auto& lhs_target = lhs.at(target_row);
    const auto& lhs_source = lhs.at(source_row);
    const auto& target_expr = lhs_target.at(source_row);
    const auto& source_expr = lhs_source.at(source_row);
    const auto factor = (-target_expr * EF::invert(source_expr))->simplify();
    const auto weighted_add = [&factor](const auto& source,
                                        const auto& target) {
      auto factor_times_source = (factor * source)->simplify();
      return (target + factor_times_source)->simplify();
    };
    for (size_t i = 0; i < lhs_source.size(); ++i) {
      lhs_target.at(i) = weighted_add(lhs_source.at(i), lhs_target.at(i));
    }
    rhs.at(target_row) = weighted_add(rhs.at(source_row), rhs.at(target_row));
  }

  lhs.erase(lhs.begin() + source_row);
  for (auto& lhs_row : lhs) {
    lhs_row.erase(lhs_row.begin() + source_row);
  }
  rhs.erase(rhs.begin() + source_row);
}
}  // namespace SymbolicOptimization