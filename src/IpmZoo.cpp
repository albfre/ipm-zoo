#include <chrono>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>

#include "Assert.h"
#include "Evaluation.h"
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

void runEvaluationExample() {
  using namespace Expression;
  using namespace Evaluation;

  printHeader("Numerical Evaluation Example");

  // Create some symbolic expressions
  auto x = ExprFactory::variable("x");
  auto y = ExprFactory::variable("y");
  auto z = ExprFactory::variable("z");
  auto A = ExprFactory::matrix("A");
  auto B = ExprFactory::symmetricMatrix("B");
  auto Q = ExprFactory::symmetricMatrix("Q");
  auto c = ExprFactory::namedScalar("c");

  printSubHeader("Setting Up Evaluation Environment");
  std::cout << "Creating an environment with concrete values for variables:"
            << std::endl;

  // Create evaluation environment and populate with values
  Environment env;

  // Vector x = [1, 2, 3]
  env[x] = valVector({1.0, 2.0, 3.0});
  std::cout << "  x = [1, 2, 3]" << std::endl;

  // Vector y = [4, 5, 6]
  env[y] = valVector({4.0, 5.0, 6.0});
  std::cout << "  y = [4, 5, 6]" << std::endl;

  // Scalar c = 2.5
  env[c] = valScalar(2.5);
  std::cout << "  c = 2.5" << std::endl;

  // Matrix A = [[1, 2, 3], [3, 4, 5], [5, 6, 7], [7, 8, 9]]
  ValMatrix matrixA = {
      {1.0, 2.0, 3.0}, {3.0, 4.0, 5.0}, {5.0, 6.0, 7.0}, {7.0, 8.0, 9.0}};
  env[A] = matrixA;
  std::cout << "  A = [[1, 2, 3], [3, 4, 5], [5, 6, 7], [7, 8, 9]]"
            << std::endl;

  // Symmetric matrix Q = [[1, 2, 3], [2, 4, 5], [3, 5, 6]]
  ValMatrix matrixQ = {{1.0, 2.0, 3.0}, {2.0, 4.0, 5.0}, {3.0, 5.0, 6.0}};
  env[Q] = matrixQ;
  std::cout << "  Q = [[1, 2, 3], [2, 4, 5], [3, 5, 6]]" << std::endl;

  // Example 1: Vector dot product x^T y
  printSubHeader("Example 1: Vector Dot Product");
  auto dotProduct = ExprFactory::product({ExprFactory::transpose(x), y});
  std::cout << "Expression: " << dotProduct.toString() << std::endl;

  auto result1 = evaluate(dotProduct, env);
  std::cout << "Result: ";
  if (is<ValScalar>(result1)) {
    std::cout << std::get<ValScalar>(result1) << std::endl;
  } else {
    std::cout << "Error: Expected scalar result but was " << result1.index()
              << std::endl;
  }
  std::cout << "Expected: 1*4 + 2*5 + 3*6 = 32" << std::endl;

  // Example 2: Scaled vector c * x
  printSubHeader("Example 2: Scaled Vector");
  auto scaledVector = ExprFactory::product({c, x});
  std::cout << "Expression: " << scaledVector.toString() << std::endl;

  auto result2 = evaluate(scaledVector, env);
  std::cout << "Result: [";
  if (is<ValVector>(result2)) {
    const auto& vec = std::get<ValVector>(result2);
    for (size_t i = 0; i < vec.size(); ++i) {
      std::cout << vec[i] << (i < vec.size() - 1 ? ", " : "");
    }
  }
  std::cout << "]" << std::endl;
  std::cout << "Expected: [2.5*1, 2.5*2, 2.5*3] = [2.5, 5.0, 7.5]" << std::endl;

  // Example 3: Quadratic form x^T Q x
  printSubHeader("Example 3: Quadratic Form");
  auto quadraticForm = ExprFactory::product({ExprFactory::transpose(x), Q, x});
  std::cout << "Expression: " << quadraticForm.toString() << std::endl;

  auto result3 = evaluate(quadraticForm, env);
  std::cout << "Result: ";
  if (is<ValScalar>(result3)) {
    std::cout << std::get<ValScalar>(result3) << std::endl;
  } else {
    std::cout << "Error: Expected scalar result" << std::endl;
  }
  std::cout << "Expected: x^T Q x = 1*1*1 + 1*2*2 + 1*3*3 + 2*2*1 + 2*4*2 + "
               "2*5*3 + 3*3*1 + 3*5*2 + 3*6*3 = 157"
            << std::endl;

  // Example 4: Matrix-vector product A * x
  printSubHeader("Example 4: Matrix-Vector Product");
  auto matrixVectorProduct = ExprFactory::product({A, x});
  std::cout << "Expression: " << matrixVectorProduct.toString() << std::endl;
  ASSERT(is<ValVector>(env.at(x)));

  auto result4 = evaluate(matrixVectorProduct, env);
  std::cout << "Result: [";
  if (is<ValVector>(result4)) {
    const auto& vec = std::get<ValVector>(result4);
    for (size_t i = 0; i < vec.size(); ++i) {
      std::cout << vec[i] << (i < vec.size() - 1 ? ", " : "");
    }
  }
  std::cout << "]" << std::endl;
  std::cout << "Expected: [1*1+2*2+3*3, 3*1+4*2+5*3, 5*1+6*2+7*3, 7*1+8*2+9*3] "
               "= [14, 26, 38, 50]"
            << std::endl;

  // Example 5: Complex expression combining multiple operations
  printSubHeader("Example 5: Complex Expression");
  // 0.5 * x^T Q x + c * y^T x
  auto complexExpr = ExprFactory::sum(
      {ExprFactory::product(
           {ExprFactory::number(0.5), ExprFactory::transpose(x), Q, x}),
       ExprFactory::product({c, ExprFactory::transpose(y), x})});

  std::cout << "Expression: " << complexExpr.toString() << std::endl;
  std::cout << "Simplified: " << complexExpr.simplify().toString() << std::endl;

  auto result5 = evaluate(complexExpr, env);
  std::cout << "Result: ";
  if (is<ValScalar>(result5)) {
    std::cout << std::get<ValScalar>(result5) << std::endl;
  } else {
    std::cout << "Error: Expected scalar result" << std::endl;
  }
  std::cout << "Expected: 0.5 * 157 + 2.5 * 32 = 78.5 + 80 = 158.5"
            << std::endl;
}

void printUsage() {
  std::cout << "Usage: IpmZoo [option]" << std::endl;
  std::cout << "Options:" << std::endl;
  std::cout << "  -h, --help     : Show this help message" << std::endl;
  std::cout << "  -b, --basic    : Run basic expression examples" << std::endl;
  std::cout << "  -o, --optimize : Run symbolic optimization example"
            << std::endl;
  std::cout << "  -e, --evaluate : Run evaluation example" << std::endl;
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
    } else if (arg == "-e" || arg == "--evaluate") {
      runEvaluationExample();
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
  runEvaluationExample();

  std::cout << "\nDone!" << std::endl;

  return 0;
}