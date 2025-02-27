#include <cassert>
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

Expression::Expr getLagrangian(std::vector<Expression::Expr>& variables) {
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
  variables = {x,        lambda_A, s, lambda_g, lambda_t, lambda_y,
               lambda_z, g,        t, y,        z};

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

void printLhs(const std::vector<std::vector<Expression::Expr>>& lhs) {
  for (const auto& row : lhs) {
    for (size_t i = 0; i < row.size(); ++i) {
      std::cout << row[i].toString();
      if (i < row.size() - 1) std::cout << " & ";
    }
    std::cout << " \\\\" << std::endl;
  }
}

void printRhs(const std::vector<Expression::Expr>& rhs) {
  for (const auto& c : rhs) {
    std::cout << c.toString() << " \\\\" << std::endl;
  }
}

void gaussianElimination(std::vector<std::vector<Expression::Expr>>& lhs,
                         std::vector<Expression::Expr>& rhs,
                         const size_t sourceRow, const bool print = false) {
  using namespace Expression;
  const auto zero = ExprFactory::number(0.0);
  size_t targetRow = 0;
  for (; targetRow < lhs.size(); ++targetRow) {
    if (targetRow != sourceRow && lhs.at(targetRow).at(sourceRow) != zero) {
      break;
    }
  }
  assert(targetRow < lhs.size());
  const auto targetExpr = lhs.at(targetRow).at(sourceRow);
  const auto sourceExpr = lhs.at(sourceRow).at(sourceRow);
  const auto factor =
      ExprFactory::negate(
          ExprFactory::product({ExprFactory::invert(sourceExpr), targetExpr}))
          .simplify();
  const auto addRowTimesFactorToRow = [&factor](const Expr& sourceTerm,
                                                const Expr& targetTerm) {
    const auto sourceTermTimesFactor =
        ExprFactory::product({factor, sourceTerm}).simplify();
    return ExprFactory::sum({targetTerm, sourceTermTimesFactor}).simplify();
  };
  if (print) {
    std::cout << "source: " << sourceRow << ", target: " << targetRow
              << std::endl;
    std::cout << "factor: " << factor.toString() << std::endl;
  }

  for (size_t i = 0; i < lhs.at(sourceRow).size(); ++i) {
    if (print) {
      std::cout << i << "add "
                << ExprFactory::product({factor, lhs.at(sourceRow).at(i)})
                       .simplify()
                       .toString()
                << " to row " << targetRow << " with result "
                << addRowTimesFactorToRow(lhs.at(sourceRow).at(i),
                                          lhs.at(targetRow).at(i))
                       .toString()
                << std::endl;
    }
    lhs.at(targetRow).at(i) = addRowTimesFactorToRow(lhs.at(sourceRow).at(i),
                                                     lhs.at(targetRow).at(i));
  }
  rhs.at(targetRow) =
      addRowTimesFactorToRow(rhs.at(sourceRow), rhs.at(targetRow));

  lhs.erase(lhs.begin() + sourceRow);
  for (auto& lhsRow : lhs) {
    lhsRow.erase(lhsRow.begin() + sourceRow);
  }
  rhs.erase(rhs.begin() + sourceRow);
}

void runLagrangianTest() {
  using namespace Expression;

  std::vector<Expr> variables;
  auto lagrangian = getLagrangian(variables);
  std::cout << "\nLagrangian: " << lagrangian.toString() << "\n";

  auto lhs = std::vector<std::vector<Expr>>();
  auto rhs = std::vector<Expr>();

  for (auto& v : variables) {
    const auto invV = ExprFactory::invert(v);
    lhs.emplace_back();
    auto& row = lhs.back();
    auto diff = lagrangian.differentiate(v).simplify();
    std::cout << v.toString() << ": " << diff.toString() << std::endl;

    if (diff.containsSubexpression(invV)) {
      std::cout << "contains invV" << std::endl;
      diff = ExprFactory::product({v, diff}).simplify();
      std::cout << "diff * v: " << diff.toString() << std::endl;
    }
    rhs.push_back(ExprFactory::negate(diff).simplify());
    for (auto& v2 : variables) {
      row.push_back(diff.differentiate(v2).simplify());
    }
  }

  std::cout << "Lhs matrix:" << std::endl;
  printLhs(lhs);

  std::cout << "Rhs: " << std::endl;
  printRhs(rhs);
  const auto createRhs = [](const std::vector<Expr>& variables) {
    std::vector<Expr> rhs;
    for (const auto& var : variables) {
      rhs.push_back(ExprFactory::negate(
          ExprFactory::namedConstant("r_{" + var.toString() + "}")));
    }
    return rhs;
  };

  rhs = createRhs(variables);

  std::cout << "\n\n";

  gaussianElimination(lhs, rhs, lhs.size() - 1);
  gaussianElimination(lhs, rhs, lhs.size() - 1);
  gaussianElimination(lhs, rhs, lhs.size() - 1);
  gaussianElimination(lhs, rhs, lhs.size() - 1);
  gaussianElimination(lhs, rhs, lhs.size() - 1);
  gaussianElimination(lhs, rhs, lhs.size() - 1);

  std::cout << "Lhs matrix:" << std::endl;
  printLhs(lhs);

  gaussianElimination(lhs, rhs, lhs.size() - 1);
  gaussianElimination(lhs, rhs, lhs.size() - 1);
  gaussianElimination(lhs, rhs, lhs.size() - 1);
  std::cout << "Lhs matrix:" << std::endl;
  printLhs(lhs);

  std::cout << "Rhs: " << std::endl;
  printRhs(rhs);
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