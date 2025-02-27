#include "Lagrangian.h"

namespace Lagrangian {
std::pair<Expression::Expr, std::vector<Expression::Expr>> getLagrangian(
    VariableNames names, Settings settings) {
  using namespace Expression::ExprFactory;
  auto Q = namedConstant("Q");
  auto c = namedConstant("c");
  auto A_ineq = namedConstant(names.A_ineq);
  auto A_eq = namedConstant(names.A_eq);
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

  auto xQx = product({number(0.5), x, Q, x});
  auto cx = product({c, x});
  auto Ax = product({A_ineq, x});

  auto terms = std::vector{xQx, cx};
  auto variables = std::vector{x};
  auto suffixVariables = std::vector<Expression::Expr>();

  if (settings.inequalities != Bounds::None) {
    variables.push_back(lambda_A);
    variables.push_back(s_A);
    terms.push_back(product({lambda_A, sum({Ax, negate(s_A)})}));
  }
  if (settings.inequalities == Bounds::Lower ||
      settings.inequalities == Bounds::Both) {
    variables.push_back(lambda_sAl);
    suffixVariables.push_back(s_Al);
    terms.push_back(
        negate(product({lambda_sAl, sum({s_A, negate(s_Al), negate(l_A)})})));
    terms.push_back(negate(product({mu, e, log({s_Al})})));
  }
  if (settings.inequalities == Bounds::Upper ||
      settings.inequalities == Bounds::Both) {
    variables.push_back(lambda_sAu);
    suffixVariables.push_back(s_Au);
    terms.push_back(product({lambda_sAu, sum({s_A, s_Au, negate(u_A)})}));
    terms.push_back(negate(product({mu, e, log({s_Au})})));
  }
  if (settings.variableBounds == Bounds::Lower ||
      settings.variableBounds == Bounds::Both) {
    variables.push_back(lambda_sxl);
    suffixVariables.push_back(s_xl);
    terms.push_back(
        negate(product({lambda_sxl, sum({x, negate(s_xl), negate(l_x)})})));
    terms.push_back(negate(product({mu, e, log({s_xl})})));
  }
  if (settings.variableBounds == Bounds::Upper ||
      settings.variableBounds == Bounds::Both) {
    variables.push_back(lambda_sxu);
    suffixVariables.push_back(s_xu);
    terms.push_back(product({lambda_sxu, sum({x, s_xu, negate(u_x)})}));
    terms.push_back(negate(product({mu, e, log({s_xu})})));
  }
  variables.insert(variables.end(), suffixVariables.begin(),
                   suffixVariables.end());
  auto lagrangian = sum(terms);
  return {lagrangian, variables};
}

}  // namespace Lagrangian