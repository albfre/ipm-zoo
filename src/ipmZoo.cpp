#include <cassert>
#include <chrono>
#include <iostream>
#include <memory>

#include "Expression.h"
#include "Optimization.h"

void initialTest() {
  using namespace Expression;

  auto A = ExprFactory::variable("A");
  auto X = ExprFactory::variable("X");
  auto B = ExprFactory::variable("B");
  auto C = ExprFactory::variable("C");
  auto x = ExprFactory::variable("x");
  auto y = ExprFactory::variable("y");
  auto z = ExprFactory::variable("z");
  auto a = ExprFactory::variable("a");

  {
    std::cout << "x^T a and a^T x" << std::endl;
    auto t1 = ExprFactory::product({ExprFactory::transpose(x), a});
    auto t2 = ExprFactory::product({ExprFactory::transpose(a), x});
    std::cout << t1.differentiate(x).simplify().toString() << std::endl;
    std::cout << t2.differentiate(x).simplify().toString() << std::endl;
  }
  {
    std::cout << "a x^T Q x and a x^T Qsym x and 0.5 x^T Qsym x" << std::endl;
    auto t1 = ExprFactory::product(
        {a, ExprFactory::transpose(x), ExprFactory::variable("Q"), x});
    auto t2 = ExprFactory::product(
        {a, ExprFactory::transpose(x), ExprFactory::symmetricMatrix("Q"), x});
    auto t3 = ExprFactory::product({ExprFactory::number(0.5),
                                    ExprFactory::transpose(x),
                                    ExprFactory::symmetricMatrix("Q"), x});
    std::cout << t1.differentiate(x).simplify().toString() << std::endl;
    std::cout << t2.differentiate(x).simplify().toString() << std::endl;
    std::cout << t3.differentiate(x).simplify().toString() << std::endl;
  }

  if (false) {
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
    std::cout << "Diff: " << expr1.differentiate(X).simplify().toString()
              << "\n";

    std::cout << "Algebraic Expression: " << expr2.toString() << "\n";
    std::cout << "Simplified: " << expr2.simplify().toString() << "\n";
    std::cout << "Diff: "
              << expr2.simplify().differentiate(x).simplify().toString()
              << "\n";

    std::cout << "Algebraic Expression: " << expr3.toString() << "\n";
    std::cout << "Simplified: " << expr3.simplify().toString() << "\n";
    std::cout << "Diff A: " << expr3.differentiate(A).simplify().toString()
              << "\n";
    std::cout << "Diff x: " << expr3.differentiate(x).simplify().toString()
              << "\n";
  }
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
  auto settings = Optimization::Settings();
  settings.inequalityHandling = Optimization::InequalityHandling::SimpleSlacks;

  const auto [lagrangian, variables] =
      Optimization::getLagrangian(Optimization::VariableNames(), settings);
  std::cout << "\nLagrangian: " << lagrangian.toString() << "\n";

  // Add timing around getNewtonSystem
  auto start = std::chrono::high_resolution_clock::now();

  auto [lhs, rhs] = Optimization::getNewtonSystem(lagrangian, variables);

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> duration = end - start;
  std::cout << "Newton system computation took " << duration.count() << " ms\n";

  std::cout << "Lhs matrix:" << std::endl;
  printLhs(lhs);

  std::cout << "Rhs: " << std::endl;
  printRhs(rhs);
  rhs = Optimization::getShorthandRhs(variables);

  std::cout << "\n\n";

  Optimization::gaussianElimination(lhs, rhs, lhs.size() - 1);
  Optimization::gaussianElimination(lhs, rhs, lhs.size() - 1);
  Optimization::gaussianElimination(lhs, rhs, lhs.size() - 1);
  Optimization::gaussianElimination(lhs, rhs, lhs.size() - 1);
  if (settings.inequalityHandling == Optimization::InequalityHandling::Slacks) {
    Optimization::gaussianElimination(lhs, rhs, lhs.size() - 1);
    Optimization::gaussianElimination(lhs, rhs, lhs.size() - 1);
  }

  std::cout << "Lhs matrix:" << std::endl;
  printLhs(lhs);

  Optimization::gaussianElimination(lhs, rhs, lhs.size() - 1);
  Optimization::gaussianElimination(lhs, rhs, lhs.size() - 1);
  if (settings.inequalityHandling == Optimization::InequalityHandling::Slacks) {
    Optimization::gaussianElimination(lhs, rhs, lhs.size() - 1);
  }

  while (lhs.size() > 1) {
    Optimization::gaussianElimination(lhs, rhs, lhs.size() - 1);
  }

  std::cout << "Lhs matrix:" << std::endl;
  printLhs(lhs);

  std::cout << "Rhs: " << std::endl;
  printRhs(rhs);
}

int main(int argc, char* argv[]) {
  // Check if we have any command line arguments
  // initialTest();

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