#include "Optimization.h"

#include <iostream>

namespace Optimization {
std::pair<Expression::Expr, std::vector<Expression::Expr>> getLagrangian(
    const VariableNames& names,
    const Settings& settings) {
  using namespace Expression::ExprFactory;
  auto Q = symmetricMatrix("Q");
  auto c = namedConstant("c");
  auto A_ineq = matrix(names.A_ineq);
  auto A_eq = matrix(names.A_eq);
  auto mu = namedConstant("\\mu");
  auto e = namedConstant("e");
  auto x = variable(names.x);
  auto s_A = variable(names.s_A);
  auto s_Al = variable(names.s_Al);
  auto s_Au = variable(names.s_Au);
  auto s_xl = variable(names.s_xl);
  auto s_xu = variable(names.s_xu);
  auto lambda_A = variable("\\lambda_{" + names.A_ineq + "}");
  auto lambda_sAl = variable("\\lambda_{" + names.s_Al + "}");
  auto lambda_sAu = variable("\\lambda_{" + names.s_Au + "}");
  auto lambda_sxl = variable("\\lambda_{" + names.s_xl + "}");
  auto lambda_sxu = variable("\\lambda_{" + names.s_xu + "}");
  auto l_A = namedConstant(names.l_A);
  auto u_A = namedConstant(names.u_A);
  auto l_x = namedConstant(names.l_x);
  auto u_x = namedConstant(names.u_x);

  auto xQx = product({number(0.5), transpose(x), Q, x});
  auto cx = product({c, x});
  auto Ax = product({A_ineq, x});

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
    if (const auto invV = ExprFactory::invert(v);
        diff.containsSubexpression(invV)) {
      diff = ExprFactory::product({v, diff}).simplify();
    }
    std::cout << "fo " << diff.toString() << std::endl;
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
    rhs.push_back(ExprFactory::negate(c));
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

}  // namespace Optimization