#include <cassert>
#include <iostream>
#include <memory>

#include "Expression.h"
#include "GaussianElimination.h"
#include "Optimization.h"

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

void runLagrangianTest() {
  using namespace Expression;

  const auto [lagrangian, variables] = Optimization::getLagrangian(
      Optimization::VariableNames(), Optimization::Settings());
  std::cout << "\nLagrangian: " << lagrangian.toString() << "\n";

  auto [lhs, rhs] = Optimization::getNewtonSystem(lagrangian, variables);

  std::cout << "Lhs matrix:" << std::endl;
  printLhs(lhs);

  std::cout << "Rhs: " << std::endl;
  printRhs(rhs);
  rhs = Optimization::getShorthandRhs(variables);

  std::cout << "\n\n";

  GaussianElimination::gaussianElimination(lhs, rhs, lhs.size() - 1);
  GaussianElimination::gaussianElimination(lhs, rhs, lhs.size() - 1);
  GaussianElimination::gaussianElimination(lhs, rhs, lhs.size() - 1);
  GaussianElimination::gaussianElimination(lhs, rhs, lhs.size() - 1);
  GaussianElimination::gaussianElimination(lhs, rhs, lhs.size() - 1);
  GaussianElimination::gaussianElimination(lhs, rhs, lhs.size() - 1);

  std::cout << "Lhs matrix:" << std::endl;
  printLhs(lhs);

  GaussianElimination::gaussianElimination(lhs, rhs, lhs.size() - 1);
  GaussianElimination::gaussianElimination(lhs, rhs, lhs.size() - 1);
  GaussianElimination::gaussianElimination(lhs, rhs, lhs.size() - 1);
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