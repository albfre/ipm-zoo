#include <iostream>
#include <memory>

#include "Expression.h"

void initialTest() {
  using namespace Expression;

  auto A = ExprFactory::variable("A");
  auto X = ExprFactory::variable("X");
  auto B = ExprFactory::variable("B");
  auto C = ExprFactory::variable("C");
  auto x = ExprFactory::variable("x");
  auto y = ExprFactory::variable("y");

  auto expr1 = ExprFactory::sum(
      {ExprFactory::product({A, X}),
       ExprFactory::sum(
           {ExprFactory::product({B, X}),
            ExprFactory::product({C, ExprFactory::product({X, C})})})});

  auto expr2 = ExprFactory::sum(
      {ExprFactory::product({x, y}), ExprFactory::product({y, x}),
       ExprFactory::product({A, C}),
       ExprFactory::negate(ExprFactory::product({C, A})),
       ExprFactory::number(1.3), ExprFactory::number(2.3),
       ExprFactory::product({ExprFactory::number(1.3), x}), x, x});

  auto expr3 =
      ExprFactory::negate(ExprFactory::negate(ExprFactory::sum({A, B})));

  std::cout << "Matrix Expression: " << expr1.toString() << "\n";
  std::cout << "Simplified: " << expr1.simplify().toString() << "\n";
  std::cout << "Diff: " << expr1.differentiate(X).simplify().toString() << "\n";

  std::cout << "Algebraic Expression: " << expr2.toString() << "\n";
  std::cout << "Simplified: " << expr2.simplify().toString() << "\n";
  std::cout << "Diff: "
            << expr2.simplify().differentiate(x).simplify().toString() << "\n";

  std::cout << "Algebraic Expression: " << expr3.toString() << "\n";
  std::cout << "Simplified: " << expr3.simplify().toString() << "\n";
  std::cout << "Diff A: " << expr3.differentiate(A).simplify().toString()
            << "\n";
  std::cout << "Diff x: " << expr3.differentiate(x).simplify().toString()
            << "\n";
}

Expression::Expr getLagrangian() {
  using namespace Expression::ExprFactory;
  auto x = variable("x");
  auto H = namedConstant("H");
  auto c = namedConstant("c");
  auto A = namedConstant("A");
  auto mu = namedConstant("\\mu");
  auto e = namedConstant("e");
  auto g = variable("g");
  auto t = variable("t");
  auto y = variable("y");
  auto z = variable("z");
  auto s = variable("s");
  auto lambda_A = variable("\\lambda_A");
  auto lambda_g = variable("\\lambda_g");
  auto lambda_t = variable("\\lambda_t");
  auto lambda_y = variable("\\lambda_y");
  auto lambda_z = variable("\\lambda_z");
  auto l_A = namedConstant("l_A");
  auto u_A = namedConstant("u_A");
  auto l_x = namedConstant("l_x");
  auto u_x = namedConstant("u_x");

  auto xHx = product({number(0.5), x, H, x});
  auto cx = product({c, x});
  auto Ax = product({A, x});
  auto minusMuLogG = negate(product({mu, e, log({g})}));
  auto minusMuLogT = negate(product({mu, e, log({t})}));
  auto minusMuLogY = negate(product({mu, e, log({y})}));
  auto minusMuLogZ = negate(product({mu, e, log({z})}));
  auto lambda_A_term = product({lambda_A, sum({Ax, negate(s)})});
  auto lambda_g_term =
      negate(product({lambda_g, sum({s, negate(g), negate(l_A)})}));
  auto lambda_t_term = product({lambda_t, sum({s, t, negate(u_A)})});
  auto lambda_y_term =
      negate(product({lambda_y, sum({x, negate(y), negate(l_x)})}));
  auto lambda_z_term = product({lambda_z, sum({x, z, negate(u_x)})});

  auto lagrangian = sum({xHx, cx, minusMuLogG, minusMuLogT, minusMuLogY,
                         minusMuLogZ, lambda_A_term, lambda_g_term,
                         lambda_t_term, lambda_y_term, lambda_z_term});
  return lagrangian;
}

void runLagrangianTest() {
  auto lagrangian = getLagrangian();
  std::cout << "\nLagrangian: " << lagrangian.toString() << "\n";
  // std::cout << "Simp Lagrangian: " << lagrangian.simplify().toString() <<
  // "\n";
  auto variables = lagrangian.getVariables();

  std::cout << "simplified:\n";
  for (auto& v : variables) {
    if (v.toString() == "g" || true) {
      auto diff = lagrangian.differentiate(v);
      std::cout << "diff " << v.toString() << ": " << diff.simplify().toString()
                << std::endl;
      auto k = Expression::ExprFactory::product({v, diff});
      std::cout << "diff * v (" << v.toString()
                << "): " << k.simplify().toString() << std::endl;
    }
  }
}

int main(int argc, char* argv[]) {
  // Check if we have any command line arguments
  if (argc > 1) {
    std::string arg = argv[1];
    if (arg == "-i") {
      std::cout << "Running initial test...\n";
      initialTest();
      return 0;
    }
  }

  runLagrangianTest();

  return 0;
}