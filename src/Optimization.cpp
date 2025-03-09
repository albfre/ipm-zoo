#include "Optimization.h"

#include <cassert>
#include <iostream>

namespace Optimization {
std::pair<Expression::Expr, std::vector<Expression::Expr>> getLagrangian(
    const VariableNames& names, const Settings& settings) {
  using namespace Expression::ExprFactory;
  const auto Q = symmetricMatrix(names.Q);
  const auto c = namedConstant(names.c);
  const auto A_ineq = matrix(names.A_ineq);
  const auto A_eq = matrix(names.A_eq);
  const auto b_eq = matrix(names.b_eq);
  const auto mu = namedConstant("\\mu");
  const auto e = namedConstant("e");
  const auto x = variable(names.x);
  const auto s_A = variable(names.s_A);
  const auto s_Al = variable(names.s_Al);
  const auto s_Au = variable(names.s_Au);
  const auto s_xl = variable(names.s_xl);
  const auto s_xu = variable(names.s_xu);
  const auto lambda_C = variable("\\lambda_{" + names.A_eq + "}");
  const auto lambda_A = variable("\\lambda_{" + names.A_ineq + "}");
  const auto lambda_sAl = variable("\\lambda_{" + names.s_Al + "}");
  const auto lambda_sAu = variable("\\lambda_{" + names.s_Au + "}");
  const auto lambda_sxl = variable("\\lambda_{" + names.s_xl + "}");
  const auto lambda_sxu = variable("\\lambda_{" + names.s_xu + "}");
  const auto l_A = namedConstant(names.l_A);
  const auto u_A = namedConstant(names.u_A);
  const auto l_x = namedConstant(names.l_x);
  const auto u_x = namedConstant(names.u_x);

  const auto xQx = product({number(0.5), transpose(x), Q, x});
  const auto cx = product({c, x});
  const auto Ax = product({A_ineq, x});
  const auto Cx = product({A_eq, x});
  const auto CxMinusB = sum({Cx, negate(b_eq)});

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

  if (settings.equalities) {
    if (settings.equalityHandling == EqualityHandling::PenaltyFunction) {
      const auto muTerm = product({number(0.5), invert(variable("\\mu"))});
      terms.push_back(product({muTerm, transpose(CxMinusB), CxMinusB}));
    } else {
      variables.push_back(lambda_C);
      terms.push_back(product({transpose(lambda_C), CxMinusB}));
    }
  }

  if (settings.inequalities != Bounds::None &&
      settings.inequalityHandling == InequalityHandling::Slacks) {
    variables.push_back(lambda_A);
    variables.push_back(s_A);
    terms.push_back(product({transpose(lambda_A), sum({Ax, negate(s_A)})}));
  }
  if (addLowerInequalities) {
    variables.push_back(lambda_sAl);
    if (settings.inequalityHandling == InequalityHandling::Slacks) {
      terms.push_back(negate(product(
          {transpose(lambda_sAl), sum({s_A, negate(s_Al), negate(l_A)})})));
    } else {
      terms.push_back(negate(product(
          {transpose(lambda_sAl), sum({Ax, negate(s_Al), negate(l_A)})})));
    }
  }
  if (addUpperInequalities) {
    variables.push_back(lambda_sAu);
    if (settings.inequalityHandling == InequalityHandling::Slacks) {
      terms.push_back(
          product({transpose(lambda_sAu), sum({s_A, s_Au, negate(u_A)})}));
    } else {
      terms.push_back(
          product({transpose(lambda_sAu), sum({Ax, s_Au, negate(u_A)})}));
    }
  }
  if (addLowerVariableBounds) {
    variables.push_back(lambda_sxl);
    terms.push_back(negate(
        product({transpose(lambda_sxl), sum({x, negate(s_xl), negate(l_x)})})));
  }
  if (addUpperVariableBounds) {
    variables.push_back(lambda_sxu);
    terms.push_back(
        product({transpose(lambda_sxu), sum({x, s_xu, negate(u_x)})}));
  }

  // Add log barriers
  if (addLowerInequalities) {
    variables.push_back(s_Al);
    terms.push_back(negate(product({mu, transpose(e), log({s_Al})})));
  }
  if (addUpperInequalities) {
    variables.push_back(s_Au);
    terms.push_back(negate(product({mu, transpose(e), log({s_Au})})));
  }
  if (addLowerVariableBounds) {
    variables.push_back(s_xl);
    terms.push_back(negate(product({mu, transpose(e), log({s_xl})})));
  }
  if (addUpperVariableBounds) {
    variables.push_back(s_xu);
    terms.push_back(negate(product({mu, transpose(e), log({s_xu})})));
  }
  auto lagrangian = sum(terms);
  return {lagrangian, variables};
}

std::vector<Expression::Expr> getFirstOrderOptimalityConditions(
    const Expression::Expr& lagrangian,
    const std::vector<Expression::Expr>& variables) {
  using namespace Expression;
  std::vector<Expression::Expr> firstOrder;
  firstOrder.reserve(variables.size());
  for (const auto& v : variables) {
    auto diff = lagrangian.differentiate(v).simplify();
    if (const auto invV = ExprFactory::invert(ExprFactory::diagonalMatrix(v));
        diff.containsSubexpression(invV)) {
      diff = ExprFactory::product({ExprFactory::diagonalMatrix(v), diff})
                 .simplify();
    }
    firstOrder.push_back(diff);
  }
  return firstOrder;
}

std::pair<std::vector<std::vector<Expression::Expr>>,
          std::vector<Expression::Expr>>
getNewtonSystem(const Expression::Expr& lagrangian,
                const std::vector<Expression::Expr>& variables) {
  using namespace Expression;
  auto lhs = std::vector<std::vector<Expr>>();
  auto rhs = std::vector<Expr>();
  auto firstOrder = getFirstOrderOptimalityConditions(lagrangian, variables);
  for (auto& c : firstOrder) {
    lhs.emplace_back();
    auto& row = lhs.back();
    rhs.push_back(ExprFactory::negate(c).simplify());
    for (auto& v : variables) {
      row.push_back(c.differentiate(v).simplify());
    }
  }
  return {lhs, rhs};
}

std::vector<Expression::Expr> getShorthandRhs(
    const std::vector<Expression::Expr>& variables) {
  std::vector<Expression::Expr> rhs;
  for (const auto& var : variables) {
    rhs.push_back(Expression::ExprFactory::negate(
        Expression::ExprFactory::namedConstant("r_{" + var.toString() + "}")));
  }
  return rhs;
}

Expression::Expr deltaDefinition(
    const std::vector<std::vector<Expression::Expr>>& lhs,
    const std::vector<Expression::Expr>& rhs,
    const std::vector<Expression::Expr>& variables, const size_t sourceRow) {
  using namespace Expression;
  assert(lhs.size() == rhs.size());
  assert(lhs.size() <= variables.size());
  assert(sourceRow < lhs.size());
  const auto zero = ExprFactory::number(0.0);

  const auto sourceExpr = lhs.at(sourceRow).at(sourceRow);
  auto terms = lhs.at(sourceRow);
  assert(terms.size() <= variables.size());
  for (size_t i = 0; i < terms.size(); ++i) {
    auto deltaVariable =
        ExprFactory::variable("\\Delta " + variables.at(i).getName());
    terms[i] = ExprFactory::product({terms[i], std::move(deltaVariable)});
  }

  terms.erase(terms.begin() + sourceRow);
  auto sum = ExprFactory::sum(std::move(terms));

  return ExprFactory::product(
             {ExprFactory::invert(sourceExpr),
              ExprFactory::sum({rhs.at(sourceRow), ExprFactory::negate(sum)})})
      .simplify();
}

void gaussianElimination(std::vector<std::vector<Expression::Expr>>& lhs,
                         std::vector<Expression::Expr>& rhs,
                         const size_t sourceRow) {
  using namespace Expression;
  assert(lhs.size() == rhs.size());
  assert(sourceRow < lhs.size());
  const auto zero = ExprFactory::number(0.0);
  std::set<size_t> targetRows;
  for (size_t i = 0; i < lhs.size(); ++i) {
    if (i != sourceRow && lhs.at(i).at(sourceRow) != zero) {
      targetRows.insert(i);
    }
  }
  assert(!targetRows.empty());
  for (auto& targetRow : targetRows) {
    const auto targetExpr = lhs.at(targetRow).at(sourceRow);
    const auto sourceExpr = lhs.at(sourceRow).at(sourceRow);
    const auto factor =
        ExprFactory::negate(
            ExprFactory::product({targetExpr, ExprFactory::invert(sourceExpr)}))
            .simplify();
    const auto addRowTimesFactorToRow = [&factor](const Expr& sourceTerm,
                                                  const Expr& targetTerm) {
      const auto sourceTermTimesFactor =
          ExprFactory::product({factor, sourceTerm}).simplify();
      return ExprFactory::sum({targetTerm, sourceTermTimesFactor}).simplify();
    };
    for (size_t i = 0; i < lhs.at(sourceRow).size(); ++i) {
      lhs.at(targetRow).at(i) = addRowTimesFactorToRow(lhs.at(sourceRow).at(i),
                                                       lhs.at(targetRow).at(i));
    }

    rhs.at(targetRow) =
        addRowTimesFactorToRow(rhs.at(sourceRow), rhs.at(targetRow));
  }

  lhs.erase(lhs.begin() + sourceRow);
  for (auto& lhsRow : lhs) {
    lhsRow.erase(lhsRow.begin() + sourceRow);
  }
  rhs.erase(rhs.begin() + sourceRow);
}
}  // namespace Optimization