#include "SymbolicOptimization.h"

#include "ExprFactory.h"
#include "SymbolicOperators.h"
#include "Utils/Assert.h"
#include "Utils/Helpers.h"

namespace SymbolicOptimization {
std::pair<Expression::ExprPtr, std::vector<Expression::ExprPtr>> get_lagrangian(
    const Settings& settings, const VariableNames& names) {
  using EF = Expression::ExprFactory;
  using namespace Expression;

  const auto Q = EF::symmetric_matrix(names.Q);
  const auto c = EF::named_vector(names.c);
  const auto A_ineq = EF::matrix(names.A_ineq);
  const auto A_eq = EF::matrix(names.A_eq);
  const auto b_eq = EF::matrix(names.b_eq);
  const auto p_eq = EF::variable(names.p_eq);
  const auto delta_eq = EF::named_scalar(names.delta_eq);
  const auto mu = EF::named_scalar("\\mu");
  const auto e = EF::named_vector("e");
  const auto x = EF::variable(names.x);
  const auto s_A = EF::variable(names.s_A);
  const auto s_Al = EF::variable(names.s_Al);
  const auto s_Au = EF::variable(names.s_Au);
  const auto s_xl = EF::variable(names.s_xl);
  const auto s_xu = EF::variable(names.s_xu);
  const auto s_C = EF::variable(names.s_C);
  const auto s_Cl = EF::variable(names.s_Cl);
  const auto s_Cu = EF::variable(names.s_Cu);
  const auto lambda_C = EF::variable("\\lambda_{" + names.A_eq + "}");
  const auto lambda_sCl = EF::variable("\\lambda_{" + names.s_Cl + "}");
  const auto lambda_sCu = EF::variable("\\lambda_{" + names.s_Cu + "}");
  const auto lambda_A = EF::variable("\\lambda_{" + names.A_ineq + "}");
  const auto lambda_sAl = EF::variable("\\lambda_{" + names.s_Al + "}");
  const auto lambda_sAu = EF::variable("\\lambda_{" + names.s_Au + "}");
  const auto lambda_sxl = EF::variable("\\lambda_{" + names.s_xl + "}");
  const auto lambda_sxu = EF::variable("\\lambda_{" + names.s_xu + "}");
  const auto l_A = EF::named_vector(names.l_A);
  const auto u_A = EF::named_vector(names.u_A);
  const auto l_x = EF::named_vector(names.l_x);
  const auto u_x = EF::named_vector(names.u_x);

  const auto xQx = EF::number(0.5) * EF::transpose(x) * Q * x;
  const auto cx = EF::transpose(c) * x;
  const auto eT = EF::transpose(e);
  const auto Ax = A_ineq * x;
  const auto Cx = A_eq * x;
  const auto CxMinusB = Cx - b_eq;

  auto terms = std::vector{xQx, cx};
  auto variables = std::vector{x};
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
    variables.push_back(lambda_A);
    terms.push_back(EF::transpose(lambda_A) * (Ax - s_A));
    if (addLowerInequalities) {
      terms.push_back(-EF::transpose(lambda_sAl) * (s_A - s_Al - l_A));
    }
    if (addUpperInequalities) {
      terms.push_back(EF::transpose(lambda_sAu) * (s_A + s_Au - u_A));
    }
  }

  if (hasSimplySlackedInequalities && addLowerInequalities) {
    variables.push_back(lambda_sAl);
    terms.push_back(-EF::transpose(lambda_sAl) * (Ax - s_Al - l_A));
  }
  if (hasSimplySlackedInequalities && addUpperInequalities) {
    variables.push_back(lambda_sAu);
    terms.push_back(EF::transpose(lambda_sAu) * (Ax + s_Au - u_A));
  }

  if (settings.equalities) {
    switch (settings.equality_handling) {
      case EqualityHandling::PenaltyFunction: {
        const auto muTerm = EF::number(0.5) * EF::invert(mu);
        terms.insert(terms.begin() + 2,
                     muTerm * EF::transpose(CxMinusB) * CxMinusB);
        break;
      }
      case EqualityHandling::PenaltyFunctionWithExtraVariable: {
        terms.push_back(EF::transpose(lambda_C) * CxMinusB);
        const auto muTerm = EF::number(0.5) * mu;
        terms.push_back(-muTerm * EF::transpose(lambda_C) * lambda_C);
        break;
      }
      case EqualityHandling::Slacks:
        terms.push_back(EF::transpose(lambda_C) * (Cx - s_C));
        terms.push_back(-EF::transpose(lambda_sCl) * (s_C - s_Cl - b_eq));
        terms.push_back(EF::transpose(lambda_sCu) * (s_C + s_Cu - b_eq));
        break;
      case EqualityHandling::SimpleSlacks:
        terms.push_back(-EF::transpose(lambda_sCl) * (Cx - s_Cl - b_eq));
        terms.push_back(EF::transpose(lambda_sCu) * (Cx + s_Cu - b_eq));
        break;
      case EqualityHandling::Regularization: {
        const auto pterm = EF::number(0.5) * EF::transpose(p_eq) * p_eq;
        terms.insert(terms.begin() + 2, pterm);
        terms.push_back(EF::transpose(lambda_C) *
                        (CxMinusB + (delta_eq * p_eq)));
        break;
      }
      case EqualityHandling::None: {
        terms.push_back(EF::transpose(lambda_C) * CxMinusB);
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
    variables.push_back(lambda_C);
  }

  if (hasSimplySlackedEqualities) {
    variables.push_back(lambda_sCl);
    variables.push_back(lambda_sCu);
  }

  if (hasFullySlackedInequalities) {
    variables.push_back(s_A);
  }

  if (hasFullySlackedEqualities) {
    variables.push_back(s_C);
  }

  if (hasRegularizedEqualities) {
    variables.push_back(p_eq);
  }

  if (hasFullySlackedInequalities) {
    if (addLowerInequalities) {
      variables.push_back(lambda_sAl);
    }
    if (addUpperInequalities) {
      variables.push_back(lambda_sAu);
    }
  }

  if (hasFullySlackedEqualities) {
    variables.push_back(lambda_sCl);
    variables.push_back(lambda_sCu);
  }

  if (addLowerVariableBounds) {
    variables.push_back(lambda_sxl);
    terms.push_back(-EF::transpose(lambda_sxl) * (x - s_xl - l_x));
  }
  if (addUpperVariableBounds) {
    variables.push_back(lambda_sxu);
    terms.push_back(EF::transpose(lambda_sxu) * (x + s_xu - u_x));
  }

  // Add log barriers
  if (addLowerInequalities) {
    variables.push_back(s_Al);
    terms.push_back(-mu * eT * EF::log(s_Al));
  }
  if (addUpperInequalities) {
    variables.push_back(s_Au);
    terms.push_back(-mu * eT * EF::log(s_Au));
  }
  if (hasAnySlackedEqualities) {
    variables.push_back(s_Cl);
    terms.push_back(-mu * eT * EF::log(s_Cl));
    variables.push_back(s_Cu);
    terms.push_back(-mu * eT * EF::log(s_Cu));
  }
  if (addLowerVariableBounds) {
    variables.push_back(s_xl);
    terms.push_back(-mu * eT * EF::log(s_xl));
  }
  if (addUpperVariableBounds) {
    variables.push_back(s_xu);
    terms.push_back(-mu * eT * EF::log(s_xu));
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
    auto delta_variable = Expression::ExprFactory::variable(
        "\\Delta " + variables.at(lhs.size() - 1)->to_string());
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

std::vector<Expression::ExprPtr> get_shorthand_rhs(
    const std::vector<Expression::ExprPtr>& variables) {
  using EF = Expression::ExprFactory;

  std::vector<Expression::ExprPtr> rhs;
  for (const auto& var : variables) {
    rhs.push_back(EF::negate(EF::named_vector("r_{" + var->to_string() + "}")));
  }
  return rhs;
}

Expression::ExprPtr delta_definition(
    const std::vector<std::vector<Expression::ExprPtr>>& lhs,
    const std::vector<Expression::ExprPtr>& rhs,
    const std::vector<Expression::ExprPtr>& variables,
    const size_t source_row) {
  using EF = Expression::ExprFactory;
  using Expression::is;
  using Expression::Variable;

  ASSERT(lhs.size() == rhs.size());
  ASSERT(lhs.size() <= variables.size());
  ASSERT(source_row < lhs.size());

  const auto& lhs_source_row = lhs.at(source_row);
  const auto& source_expr = lhs_source_row.at(source_row);
  auto terms = std::vector<Expression::ExprPtr>();
  terms.reserve(lhs_source_row.size());
  ASSERT(lhs_source_row.size() <= variables.size());
  for (size_t i = 0; i < lhs_source_row.size(); ++i) {
    ASSERT(is<Variable>(variables.at(i)));
    auto delta_variable =
        EF::variable("\\Delta " + variables.at(i)->to_string());
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