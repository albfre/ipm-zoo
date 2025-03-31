#include "SymbolicOptimization.h"

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

std::pair<Expression::ExprPtr, std::vector<Expression::ExprPtr>> get_lagrangian(
    const Settings& settings, const VariableNames& names) {
  using EF = Expression::ExprFactory;
  using namespace Expression;
  const auto o = get_optimization_expressions(names);

  const auto xQx = EF::number(0.5) * EF::transpose(o.x) * o.Q * o.x;
  const auto cx = EF::transpose(o.c) * o.x;
  const auto exT = EF::transpose(o.e_var);
  const auto eAT = EF::transpose(o.e_ineq);
  const auto eCT = EF::transpose(o.e_eq);
  const auto Ax = o.A_ineq * o.x;
  const auto Cx = o.A_eq * o.x;
  const auto CxMinusB = Cx - o.b_eq;

  auto terms = std::vector{xQx, cx};
  auto variables = std::vector{o.x};
  const auto addLowerInequalities = settings.inequalities == Bounds::Lower ||
                                    settings.inequalities == Bounds::Both;
  const auto addUpperInequalities = settings.inequalities == Bounds::Upper ||
                                    settings.inequalities == Bounds::Both;
  const auto addLowerVariableBounds =
      settings.variable_bounds == Bounds::Lower ||
      settings.variable_bounds == Bounds::Both;
  const auto addUpperVariableBounds =
      settings.variable_bounds == Bounds::Upper ||
      settings.variable_bounds == Bounds::Both;

  const auto hasFullySlackedEqualities =
      settings.equalities &&
      settings.equality_handling == EqualityHandling::Slacks;
  const auto hasSimplySlackedEqualities =
      settings.equalities &&
      settings.equality_handling == EqualityHandling::SimpleSlacks;
  const auto hasAnySlackedEqualities =
      hasFullySlackedEqualities || hasSimplySlackedEqualities;
  const auto hasRegularizedEqualities =
      settings.equalities &&
      settings.equality_handling == EqualityHandling::Regularization;
  const auto hasFullySlackedInequalities =
      settings.inequalities != Bounds::None &&
      settings.inequality_handling == InequalityHandling::Slacks;
  const auto hasSimplySlackedInequalities =
      settings.inequalities != Bounds::None &&
      settings.inequality_handling == InequalityHandling::SimpleSlacks;

  if (hasFullySlackedInequalities) {
    variables.push_back(o.lambda_A_ineq);
    terms.push_back(EF::transpose(o.lambda_A_ineq) * (Ax - o.s_A_ineq));
    if (addLowerInequalities) {
      terms.push_back(-EF::transpose(o.lambda_sAineql) *
                      (o.s_A_ineq - o.s_A_ineq_l - o.l_A_ineq));
    }
    if (addUpperInequalities) {
      terms.push_back(EF::transpose(o.lambda_sAinequ) *
                      (o.s_A_ineq + o.s_A_ineq_u - o.u_A_ineq));
    }
  }

  if (hasSimplySlackedInequalities && addLowerInequalities) {
    variables.push_back(o.lambda_sAineql);
    terms.push_back(-EF::transpose(o.lambda_sAineql) *
                    (Ax - o.s_A_ineq_l - o.l_A_ineq));
  }
  if (hasSimplySlackedInequalities && addUpperInequalities) {
    variables.push_back(o.lambda_sAinequ);
    terms.push_back(EF::transpose(o.lambda_sAinequ) *
                    (Ax + o.s_A_ineq_u - o.u_A_ineq));
  }

  if (settings.equalities) {
    switch (settings.equality_handling) {
      case EqualityHandling::PenaltyFunction: {
        const auto muTerm = EF::number(0.5) * EF::invert(o.mu);
        terms.insert(terms.begin() + 2,
                     muTerm * EF::transpose(CxMinusB) * CxMinusB);
        break;
      }
      case EqualityHandling::PenaltyFunctionWithExtraVariable: {
        terms.push_back(EF::transpose(o.lambda_A_eq) * CxMinusB);
        const auto muTerm = EF::number(0.5) * o.mu;
        terms.push_back(-muTerm * EF::transpose(o.lambda_A_eq) * o.lambda_A_eq);
        break;
      }
      case EqualityHandling::Slacks:
        terms.push_back(EF::transpose(o.lambda_A_eq) * (Cx - o.s_A_eq));
        terms.push_back(-EF::transpose(o.lambda_sAeql) *
                        (o.s_A_eq - o.s_A_eq_l - o.b_eq));
        terms.push_back(EF::transpose(o.lambda_sAequ) *
                        (o.s_A_eq + o.s_A_eq_u - o.b_eq));
        break;
      case EqualityHandling::SimpleSlacks:
        terms.push_back(-EF::transpose(o.lambda_sAeql) *
                        (Cx - o.s_A_eq_l - o.b_eq));
        terms.push_back(EF::transpose(o.lambda_sAequ) *
                        (Cx + o.s_A_eq_u - o.b_eq));
        break;
      case EqualityHandling::Regularization: {
        const auto pterm = EF::number(0.5) * EF::transpose(o.p_eq) * o.p_eq;
        terms.insert(terms.begin() + 2, pterm);
        terms.push_back(EF::transpose(o.lambda_A_eq) *
                        (CxMinusB + (o.delta_eq * o.p_eq)));
        break;
      }
      case EqualityHandling::None: {
        terms.push_back(EF::transpose(o.lambda_A_eq) * CxMinusB);
        break;
      }
      default:
        ASSERT(false);
    }
  }
  if (settings.equalities &&
      std::set{EqualityHandling::None,
               EqualityHandling::PenaltyFunctionWithExtraVariable,
               EqualityHandling::Regularization, EqualityHandling::Slacks}
          .contains(settings.equality_handling)) {
    variables.push_back(o.lambda_A_eq);
  }

  if (hasSimplySlackedEqualities) {
    variables.push_back(o.lambda_sAeql);
    variables.push_back(o.lambda_sAequ);
  }

  if (hasFullySlackedInequalities) {
    variables.push_back(o.s_A_ineq);
  }

  if (hasFullySlackedEqualities) {
    variables.push_back(o.s_A_eq);
  }

  if (hasRegularizedEqualities) {
    variables.push_back(o.p_eq);
  }

  if (hasFullySlackedInequalities) {
    if (addLowerInequalities) {
      variables.push_back(o.lambda_sAineql);
    }
    if (addUpperInequalities) {
      variables.push_back(o.lambda_sAinequ);
    }
  }

  if (hasFullySlackedEqualities) {
    variables.push_back(o.lambda_sAeql);
    variables.push_back(o.lambda_sAequ);
  }

  if (addLowerVariableBounds) {
    variables.push_back(o.lambda_sxl);
    terms.push_back(-EF::transpose(o.lambda_sxl) * (o.x - o.s_x_l - o.l_x));
  }
  if (addUpperVariableBounds) {
    variables.push_back(o.lambda_sxu);
    terms.push_back(EF::transpose(o.lambda_sxu) * (o.x + o.s_x_u - o.u_x));
  }

  // Add log barriers
  if (addLowerInequalities) {
    variables.push_back(o.s_A_ineq_l);
    terms.push_back(-o.mu * eAT * EF::log(o.s_A_ineq_l));
  }
  if (addUpperInequalities) {
    variables.push_back(o.s_A_ineq_u);
    terms.push_back(-o.mu * eAT * EF::log(o.s_A_ineq_u));
  }
  if (hasAnySlackedEqualities) {
    variables.push_back(o.s_A_eq_l);
    terms.push_back(-o.mu * eCT * EF::log(o.s_A_eq_l));
    variables.push_back(o.s_A_eq_u);
    terms.push_back(-o.mu * eCT * EF::log(o.s_A_eq_u));
  }
  if (addLowerVariableBounds) {
    variables.push_back(o.s_x_l);
    terms.push_back(-o.mu * exT * EF::log(o.s_x_l));
  }
  if (addUpperVariableBounds) {
    variables.push_back(o.s_x_u);
    terms.push_back(-o.mu * exT * EF::log(o.s_x_u));
  }
  terms = transform(terms, [](const auto& t) { return t->simplify(); });
  auto lagrangian = EF::sum(terms);
  return {lagrangian, variables};
}

std::pair<std::vector<Expression::ExprPtr>, std::vector<Expression::ExprPtr>>
get_first_order_optimality_conditions(Settings settings,
                                      const VariableNames& names) {
  using EF = Expression::ExprFactory;

  if (settings.equality_handling == EqualityHandling::PenaltyFunction) {
    settings.equality_handling =
        EqualityHandling::PenaltyFunctionWithExtraVariable;
  }

  auto [lagrangian, variables] = get_lagrangian(settings, names);

  std::vector<Expression::ExprPtr> first_order;
  first_order.reserve(variables.size());
  for (const auto& v : variables) {
    auto diff = lagrangian->differentiate(v)->simplify();
    if (const auto invV = EF::invert(EF::diagonal_matrix(v));
        diff->contains_subexpression(invV)) {
      diff = EF::product({EF::diagonal_matrix(v), diff})->simplify();
    }
    first_order.push_back(diff);
  }
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