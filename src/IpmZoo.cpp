#include <chrono>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>

#include "Expression.h"
#include "Helpers.h"
#include "SymbolicOptimization.h"

void printHeader(const std::string& title) {
  std::cout << "\n" << std::string(80, '=') << std::endl;
  std::cout << "  " << title << std::endl;
  std::cout << std::string(80, '=') << std::endl;
}

void printSubHeader(const std::string& title) {
  std::cout << "\n" << std::string(60, '-') << std::endl;
  std::cout << "  " << title << std::endl;
  std::cout << std::string(60, '-') << std::endl;
}

void printExample(
    const std::string& description, const Expression::Expr& expr,
    const Expression::Expr& var = Expression::ExprFactory::variable("")) {
  std::cout << "\n" << description << ":" << std::endl;
  std::cout << "  Expression: " << expr.toString() << std::endl;
  std::cout << "  Simplified: " << expr.simplify().toString() << std::endl;

  if (var.toString() != "") {
    std::cout << "  Derivative wrt " << var.toString() << ": "
              << expr.differentiate(var).simplify().toString() << std::endl;
  }
}

void runBasicExamples() {
  using namespace Expression;
  printHeader("Basic Expression Examples");

  auto x = ExprFactory::variable("x");
  auto y = ExprFactory::variable("y");
  auto z = ExprFactory::variable("z");
  auto a = ExprFactory::variable("a");
  auto constant = ExprFactory::namedScalar("c");
  auto A = ExprFactory::variable("A");
  auto B = ExprFactory::variable("B");
  auto C = ExprFactory::variable("C");

  // Vector-vector dot product examples
  printSubHeader("Vector-Vector Products");

  auto dotProduct1 = ExprFactory::product({ExprFactory::transpose(x), a});
  auto dotProduct2 =
      ExprFactory::product({ExprFactory::transpose(a), x, constant});

  printExample("Dot product x^T a", dotProduct1, x);
  printExample("Scaled dot product c·a^T x", dotProduct2, x);

  // Quadratic form examples
  printSubHeader("Quadratic Forms");

  auto quadForm1 = ExprFactory::product(
      {a, ExprFactory::transpose(x), ExprFactory::variable("Q"), x});
  auto quadForm2 = ExprFactory::product(
      {a, ExprFactory::transpose(x), ExprFactory::symmetricMatrix("Q"), x});
  auto quadForm3 =
      ExprFactory::product({ExprFactory::number(0.5), ExprFactory::transpose(x),
                            ExprFactory::symmetricMatrix("Q"), x});

  printExample("Quadratic form a·x^T·Q·x", quadForm1, x);
  printExample("Quadratic form with symmetric matrix a·x^T·Q_sym·x", quadForm2,
               x);
  printExample("Standard quadratic form 0.5·x^T·Q_sym·x", quadForm3, x);

  // Complex algebraic expressions
  printSubHeader("Matrix and Algebraic Expressions");

  auto expr1 = ExprFactory::sum(
      {ExprFactory::product({A, x}),
       ExprFactory::sum(
           {ExprFactory::product({B, x}),
            ExprFactory::product({C, ExprFactory::product({x, C})})})});

  auto expr2 = ExprFactory::sum(
      {ExprFactory::product({x, y}), ExprFactory::product({y, x}),
       ExprFactory::product({A, C}),
       ExprFactory::negate(ExprFactory::product({C, A})),
       ExprFactory::number(1.3), ExprFactory::number(2.3),
       ExprFactory::product({ExprFactory::number(1.3), x}), x, x});

  auto expr3 =
      ExprFactory::negate(ExprFactory::negate(ExprFactory::sum({A, B})));

  printExample("Matrix expression with nested terms", expr1, x);
  printExample("Complex algebraic expression with coefficients", expr2, x);
  printExample("Double negation example", expr3, A);
}

void printLhs(const std::vector<std::vector<Expression::Expr>>& lhs) {
  std::stringstream ss;
  size_t maxRowSize = 0;
  for (const auto& row : lhs) {
    size_t currentRowSize = 0;
    for (size_t i = 0; i < row.size(); ++i) {
      const auto rowStr = row[i].toString();
      currentRowSize += rowStr.size() + (i + 1 == row.size() ? 0 : 1);
      ss << (i == 0 ? " " : "") << rowStr << ((i + 1 < row.size()) ? " " : "");
    }
    ss << std::endl;
    maxRowSize = std::max(maxRowSize, currentRowSize);
  }
  std::cout << "┌" << std::string(maxRowSize, '-') << "┐" << std::endl;
  std::cout << ss.str();
  std::cout << "└" << std::string(maxRowSize, '-') << "┘" << std::endl;
}

void printRhs(const std::vector<Expression::Expr>& rhs) {
  std::stringstream ss;
  size_t maxSize = 0;
  for (const auto& c : rhs) {
    const auto cStr = c.toString();
    ss << cStr << std::endl;
    maxSize = std::max(maxSize, cStr.size() + 1);
  }
  std::cout << "┌" << std::string(maxSize, '-') << "┐" << std::endl;
  std::cout << ss.str() << std::endl;
  std::cout << "└" << std::string(maxSize, '-') << "┘" << std::endl;
}

void runOptimizationExample() {
  using namespace Expression;

  printHeader("Symbolic Optimization Example");

  std::cout << "Creating a standard form optimization problem and solving "
            << "symbolically using Newton's method" << std::endl;

  auto settings = SymbolicOptimization::Settings();
  settings.inequalityHandling =
      SymbolicOptimization::InequalityHandling::SimpleSlacks;

  printSubHeader("Generating Lagrangian");
  const auto [lagrangian, variables] = SymbolicOptimization::getLagrangian(
      SymbolicOptimization::VariableNames(), settings);

  std::cout << "Lagrangian function: " << lagrangian.toString() << std::endl;
  std::cout << "\nOptimization variables:" << std::endl;
  for (const auto& var : variables) {
    std::cout << "  - " << var.toString() << std::endl;
  }

  // Add timing around getNewtonSystem
  printSubHeader("Computing Newton System");
  std::cout << "Building the Newton system (may take a moment)..." << std::endl;

  auto start = std::chrono::high_resolution_clock::now();
  auto [lhs, rhs] =
      SymbolicOptimization::getNewtonSystem(lagrangian, variables);
  auto end = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double, std::milli> duration = end - start;
  std::cout << "Newton system computation took " << duration.count() << " ms\n";

  std::cout << "\nInitial Newton system:" << std::endl;
  std::cout << "Left-hand side (Hessian matrix):" << std::endl;
  printLhs(lhs);

  rhs = SymbolicOptimization::getShorthandRhs(variables);
  std::cout << "Right-hand side (gradient vector):" << std::endl;
  printRhs(rhs);

  printSubHeader("Performing Gaussian Elimination");

  std::cout
      << "Applying symbolic Gaussian elimination to solve the Newton system..."
      << std::endl;

  // Perform a few elimination steps to demonstrate the process
  for (int i = 0; i < 3; i++) {
    std::cout << "\nElimination step " << i + 1 << ":" << std::endl;
    SymbolicOptimization::gaussianElimination(lhs, rhs, lhs.size() - 1);

    if (i % 2 == 1) {  // Only show every other step to save space
      std::cout << "Matrix after elimination:" << std::endl;
      printLhs(lhs);
    }
  }

  // Finish elimination
  std::cout << "\nCompleting elimination..." << std::endl;
  while (lhs.size() > 1) {
    SymbolicOptimization::gaussianElimination(lhs, rhs, lhs.size() - 1);
  }

  std::cout << "\nFinal reduced system:" << std::endl;
  std::cout << "Reduced matrix:" << std::endl;
  printLhs(lhs);

  std::cout << "Reduced right-hand side:" << std::endl;
  printRhs(rhs);

  std::cout << "\nThe solution can now be obtained by back-substitution."
            << std::endl;
}

void printUsage() {
  std::cout << "Usage: IpmZoo [option]" << std::endl;
  std::cout << "Options:" << std::endl;
  std::cout << "  -h, --help     : Show this help message" << std::endl;
  std::cout << "  -b, --basic    : Run basic expression examples" << std::endl;
  std::cout << "  -o, --optimize : Run symbolic optimization example"
            << std::endl;
  std::cout << "  (no options)   : Run all examples" << std::endl;
}

int main(int argc, char* argv[]) {
  if (argc > 1) {
    std::string arg = argv[1];
    if (arg == "-h" || arg == "--help") {
      printUsage();
      return 0;
    } else if (arg == "-b" || arg == "--basic") {
      runBasicExamples();
      return 0;
    } else if (arg == "-o" || arg == "--optimize") {
      runOptimizationExample();
      return 0;
    } else {
      std::cout << "Unknown option: " << arg << std::endl;
      printUsage();
      return 1;
    }
  }

  // Default: run all examples
  std::cout << "Welcome to IpmZoo: Symbolic Differentiation and Optimization "
               "Framework"
            << std::endl;
  std::cout << "This program demonstrates symbolic mathematics capabilities "
               "including:"
            << std::endl;
  std::cout << "  - Expression differentiation and simplification" << std::endl;
  std::cout << "  - Symbolic construction of optimization problems"
            << std::endl;
  std::cout << "  - Solution of Newton systems through Gaussian elimination"
            << std::endl;

  runBasicExamples();
  runOptimizationExample();

  std::cout << "\nDone!" << std::endl;

  return 0;
}