#include "SymbolicOptimization.h"

#include "Assert.h"
#include "ExprFactory.h"
#include "Helpers.h"
#include "SymbolicOperators.h"

namespace SymbolicOptimization {
std::pair<Expression::ExprPtr, std::vector<Expression::ExprPtr>> getLagrangian(
    const VariableNames& names, const Settings& settings) {
  using EF = Expression::ExprFactory;
  using namespace Expression;

  const auto Q = EF::symmetricMatrix(names.Q);
  const auto c = EF::namedVector(names.c);
  const auto A_ineq = EF::matrix(names.A_ineq);
  const auto A_eq = EF::matrix(names.A_eq);
  const auto b_eq = EF::matrix(names.b_eq);
  const auto p_eq = EF::variable(names.p_eq);
  const auto delta_eq = EF::namedScalar(names.delta_eq);
  const auto mu = EF::namedScalar("\\mu");
  const auto e = EF::namedVector("e");
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
  const auto l_A = EF::namedVector(names.l_A);
  const auto u_A = EF::namedVector(names.u_A);
  const auto l_x = EF::namedVector(names.l_x);
  const auto u_x = EF::namedVector(names.u_x);

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
      settings.variableBounds == Bounds::Lower ||
      settings.variableBounds == Bounds::Both;
  const auto addUpperVariableBounds =
      settings.variableBounds == Bounds::Upper ||
      settings.variableBounds == Bounds::Both;

  const auto hasFullySlackedEqualities =
      settings.equalities &&
      settings.equalityHandling == EqualityHandling::Slacks;
  const auto hasSimplySlackedEqualities =
      settings.equalities &&
      settings.equalityHandling == EqualityHandling::SimpleSlacks;
  const auto hasAnySlackedEqualities =
      hasFullySlackedEqualities || hasSimplySlackedEqualities;
  const auto hasRegularizedEqualities =
      settings.equalities &&
      settings.equalityHandling == EqualityHandling::Regularization;
  const auto hasFullySlackedInequalities =
      settings.inequalities != Bounds::None &&
      settings.inequalityHandling == InequalityHandling::Slacks;
  const auto hasSimplySlackedInequalities =
      settings.inequalities != Bounds::None &&
      settings.inequalityHandling == InequalityHandling::SimpleSlacks;

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
    switch (settings.equalityHandling) {
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
          .contains(settings.equalityHandling)) {
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

std::vector<Expression::ExprPtr> getFirstOrderOptimalityConditions(
    const Expression::ExprPtr& lagrangian,
    const std::vector<Expression::ExprPtr>& variables) {
  using EF = Expression::ExprFactory;

  std::vector<Expression::ExprPtr> firstOrder;
  firstOrder.reserve(variables.size());
  for (const auto& v : variables) {
    auto diff = lagrangian->differentiate(v)->simplify();
    if (const auto invV = EF::invert(EF::diagonalMatrix(v));
        diff->containsSubexpression(invV)) {
      diff = EF::product({EF::diagonalMatrix(v), diff})->simplify();
    }
    firstOrder.push_back(diff);
  }
  return firstOrder;
}

NewtonSystem getNewtonSystem(
    const Expression::ExprPtr& lagrangian,
    const std::vector<Expression::ExprPtr>& variables) {
  using EF = Expression::ExprFactory;

  auto lhs = std::vector<std::vector<Expression::ExprPtr>>();
  auto rhs = std::vector<Expression::ExprPtr>();
  auto firstOrder = getFirstOrderOptimalityConditions(lagrangian, variables);
  for (auto& c : firstOrder) {
    lhs.emplace_back();
    auto& row = lhs.back();
    for (auto& v : variables) {
      row.push_back(c->differentiate(v)->simplify());
    }
    rhs.push_back(EF::negate(c)->simplify());
  }
  return {lhs, rhs, variables};
}

namespace {
const auto getAugmentedSystemSize(
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

NewtonSystem getAugmentedSystem(NewtonSystem newtonSystem) {
  auto& [lhs, rhs, variables, deltaDefinitions] = newtonSystem;
  const auto augmentedSystemSize = getAugmentedSystemSize(lhs);
  while (lhs.size() > augmentedSystemSize) {
    auto deltaVariable = Expression::ExprFactory::variable(
        "\\Delta " + variables.at(lhs.size() - 1)->toString());
    auto deltaDefinition = SymbolicOptimization::deltaDefinition(
        lhs, rhs, variables, lhs.size() - 1);
    deltaDefinitions.push_back({deltaVariable, deltaDefinition});
    SymbolicOptimization::gaussianElimination(lhs, rhs, lhs.size() - 1);
    variables.pop_back();
  }
  return std::move(newtonSystem);
}

NewtonSystem getNormalEquations(NewtonSystem newtonSystem) {
  newtonSystem = getAugmentedSystem(std::move(newtonSystem));
  auto& [lhs, rhs, variables, deltaDefinitions] = newtonSystem;
  auto deltaVariable = Expression::ExprFactory::variable(
      "\\Delta " + variables.at(0)->toString());
  auto deltaDefinition =
      SymbolicOptimization::deltaDefinition(lhs, rhs, variables, 0);
  deltaDefinitions.push_back({deltaVariable, deltaDefinition});
  SymbolicOptimization::gaussianElimination(lhs, rhs, 0);
  variables.erase(variables.begin());
  return std::move(newtonSystem);
}

std::vector<Expression::ExprPtr> getShorthandRhs(
    const std::vector<Expression::ExprPtr>& variables) {
  using EF = Expression::ExprFactory;

  std::vector<Expression::ExprPtr> rhs;
  for (const auto& var : variables) {
    rhs.push_back(EF::negate(EF::namedVector("r_{" + var->toString() + "}")));
  }
  return rhs;
}

Expression::ExprPtr deltaDefinition(
    const std::vector<std::vector<Expression::ExprPtr>>& lhs,
    const std::vector<Expression::ExprPtr>& rhs,
    const std::vector<Expression::ExprPtr>& variables, const size_t sourceRow) {
  using EF = Expression::ExprFactory;
  using Expression::is;
  using Expression::Variable;

  ASSERT(lhs.size() == rhs.size());
  ASSERT(lhs.size() <= variables.size());
  ASSERT(sourceRow < lhs.size());

  const auto& lhsSourceRow = lhs.at(sourceRow);
  const auto& sourceExpr = lhsSourceRow.at(sourceRow);
  auto terms = std::vector<Expression::ExprPtr>();
  terms.reserve(lhsSourceRow.size());
  ASSERT(lhsSourceRow.size() <= variables.size());
  for (size_t i = 0; i < lhsSourceRow.size(); ++i) {
    ASSERT(is<Variable>(variables.at(i)));
    auto deltaVariable = EF::variable("\\Delta " + variables.at(i)->toString());
    terms.emplace_back(
        EF::product({lhsSourceRow[i], std::move(deltaVariable)}));
  }

  terms.erase(terms.begin() + sourceRow);
  auto sum = EF::sum(std::move(terms));

  return EF::product({EF::invert(sourceExpr),
                      EF::sum({rhs.at(sourceRow), EF::negate(sum)})})
      ->simplify();
}

void gaussianElimination(std::vector<std::vector<Expression::ExprPtr>>& lhs,
                         std::vector<Expression::ExprPtr>& rhs,
                         const size_t sourceRow) {
  using EF = Expression::ExprFactory;
  using Expression::zero;
  using namespace Expression;

  ASSERT(lhs.size() == rhs.size());
  ASSERT(sourceRow < lhs.size());
  std::set<size_t> targetRows;
  for (size_t i = 0; i < lhs.size(); ++i) {
    if (i != sourceRow && lhs.at(i).at(sourceRow) != zero) {
      targetRows.insert(i);
    }
  }
  ASSERT(!targetRows.empty());
  for (const auto& targetRow : targetRows) {
    auto& lhsTarget = lhs.at(targetRow);
    const auto& lhsSource = lhs.at(sourceRow);
    const auto& targetExpr = lhsTarget.at(sourceRow);
    const auto& sourceExpr = lhsSource.at(sourceRow);
    const auto factor = (-targetExpr * EF::invert(sourceExpr))->simplify();
    const auto weightedAdd = [&factor](const auto& source, const auto& target) {
      auto factorTimesSource = (factor * source)->simplify();
      return (target + factorTimesSource)->simplify();
    };
    for (size_t i = 0; i < lhsSource.size(); ++i) {
      lhsTarget.at(i) = weightedAdd(lhsSource.at(i), lhsTarget.at(i));
    }
    rhs.at(targetRow) = weightedAdd(rhs.at(sourceRow), rhs.at(targetRow));
  }

  lhs.erase(lhs.begin() + sourceRow);
  for (auto& lhsRow : lhs) {
    lhsRow.erase(lhsRow.begin() + sourceRow);
  }
  rhs.erase(rhs.begin() + sourceRow);
}
}  // namespace SymbolicOptimization