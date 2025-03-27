#include <chrono>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>

#include "Expr.h"
#include "ExprFactory.h"
#include "NumericOptimization/Evaluation.h"
#include "SymbolicOptimization.h"
#include "Utils/Assert.h"
#include "Utils/Helpers.h"

void print_header(const std::string& title) {
  std::cout << "\n" << std::string(80, '=') << std::endl;
  std::cout << "  " << title << std::endl;
  std::cout << std::string(80, '=') << std::endl;
}

void print_sub_header(const std::string& title) {
  std::cout << "\n" << std::string(60, '-') << std::endl;
  std::cout << "  " << title << std::endl;
  std::cout << std::string(60, '-') << std::endl;
}

void print_example(
    const std::string& description, const Expression::ExprPtr& expr,
    const Expression::ExprPtr& var = Expression::ExprFactory::variable("")) {
  std::cout << "\n" << description << ":" << std::endl;
  std::cout << "  Expression: " << expr->to_string() << std::endl;
  std::cout << "  Simplified: " << expr->simplify()->to_string() << std::endl;

  if (var->to_string() != "") {
    std::cout << "  Derivative wrt " << var->to_string() << ": "
              << expr->differentiate(var)->simplify()->to_string() << std::endl;
  }
}

void run_basic_examples() {
  using namespace Expression;
  print_header("Basic Expression Examples");

  auto x = ExprFactory::variable("x");
  auto y = ExprFactory::variable("y");
  auto z = ExprFactory::variable("z");
  auto a = ExprFactory::variable("a");
  auto constant = ExprFactory::named_scalar("c");
  auto A = ExprFactory::variable("A");
  auto B = ExprFactory::variable("B");
  auto C = ExprFactory::variable("C");

  // Vector-vector dot product examples
  print_sub_header("Vector-Vector Products");

  auto dot_product1 = ExprFactory::product({ExprFactory::transpose(x), a});
  auto dot_product2 =
      ExprFactory::product({ExprFactory::transpose(a), x, constant});

  print_example("Dot product x^T a", dot_product1, x);
  print_example("Scaled dot product c·a^T x", dot_product2, x);

  // Quadratic form examples
  print_sub_header("Quadratic Forms");

  auto quad_form1 = ExprFactory::product(
      {a, ExprFactory::transpose(x), ExprFactory::variable("Q"), x});
  auto quad_form2 = ExprFactory::product(
      {a, ExprFactory::transpose(x), ExprFactory::symmetric_matrix("Q"), x});
  auto quad_form3 =
      ExprFactory::product({ExprFactory::number(0.5), ExprFactory::transpose(x),
                            ExprFactory::symmetric_matrix("Q"), x});

  print_example("Quadratic form a·x^T·Q·x", quad_form1, x);
  print_example("Quadratic form with symmetric matrix a·x^T·Q_sym·x",
                quad_form2, x);
  print_example("Standard quadratic form 0.5·x^T·Q_sym·x", quad_form3, x);

  // Complex algebraic expressions
  print_sub_header("Matrix and Algebraic Expressions");

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

  print_example("Matrix expression with nested terms", expr1, x);
  print_example("Complex algebraic expression with coefficients", expr2, x);
  print_example("Double negation example", expr3, A);
}

void print_lhs(const std::vector<std::vector<Expression::ExprPtr>>& lhs) {
  std::stringstream ss;
  size_t max_row_size = 0;
  for (const auto& row : lhs) {
    size_t current_row_size = 0;
    for (size_t i = 0; i < row.size(); ++i) {
      const auto row_str = row[i]->to_string();
      current_row_size += row_str.size() + (i + 1 == row.size() ? 0 : 1);
      ss << (i == 0 ? " " : "") << row_str << ((i + 1 < row.size()) ? " " : "");
    }
    ss << std::endl;
    max_row_size = std::max(max_row_size, current_row_size);
  }
  std::cout << "┌" << std::string(max_row_size, '-') << "┐" << std::endl;
  std::cout << ss.str();
  std::cout << "└" << std::string(max_row_size, '-') << "┘" << std::endl;
}

void print_rhs(const std::vector<Expression::ExprPtr>& rhs) {
  std::stringstream ss;
  size_t max_size = 0;
  for (const auto& c : rhs) {
    const auto c_str = c->to_string();
    ss << c_str << std::endl;
    max_size = std::max(max_size, c_str.size() + 1);
  }
  std::cout << "┌" << std::string(max_size, '-') << "┐" << std::endl;
  std::cout << ss.str() << std::endl;
  std::cout << "└" << std::string(max_size, '-') << "┘" << std::endl;
}

void run_optimization_example() {
  using namespace Expression;

  print_header("Symbolic Optimization Example");

  std::cout << "Creating a standard form optimization problem and solving "
            << "symbolically using Newton's method" << std::endl;

  auto settings = SymbolicOptimization::Settings();
  settings.inequality_handling =
      SymbolicOptimization::InequalityHandling::SimpleSlacks;
  auto names = SymbolicOptimization::VariableNames();

  print_sub_header("Generating Lagrangian");
  const auto [lagrangian, variables] =
      SymbolicOptimization::get_lagrangian(settings, names);

  std::cout << "Lagrangian function: " << lagrangian->to_string() << std::endl;
  std::cout << "\nOptimization variables:" << std::endl;
  for (const auto& var : variables) {
    std::cout << "  - " << var->to_string() << std::endl;
  }

  // Add timing around get_newton_system
  print_sub_header("Computing Newton System");
  std::cout << "Building the Newton system ..." << std::endl;

  auto start = std::chrono::high_resolution_clock::now();
  auto [lhs, rhs, newton_variables, _] =
      SymbolicOptimization::get_newton_system(settings, names);
  auto end = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double, std::milli> duration = end - start;
  std::cout << "Newton system computation took " << duration.count() << " ms\n";

  std::cout << "\nInitial Newton system:" << std::endl;
  std::cout << "Left-hand side (Hessian matrix):" << std::endl;
  print_lhs(lhs);

  rhs = SymbolicOptimization::get_shorthand_rhs(variables);
  std::cout << "Right-hand side (gradient vector):" << std::endl;
  print_rhs(rhs);

  print_sub_header("Performing Gaussian Elimination");

  std::cout
      << "Applying symbolic Gaussian elimination to solve the Newton system..."
      << std::endl;

  // Perform a few elimination steps to demonstrate the process
  for (int i = 0; i < 3; i++) {
    std::cout << "\nElimination step " << i + 1 << ":" << std::endl;
    SymbolicOptimization::gaussian_elimination(lhs, rhs, lhs.size() - 1);

    if (i % 2 == 1) {  // Only show every other step to save space
      std::cout << "Matrix after elimination:" << std::endl;
      print_lhs(lhs);
    }
  }

  // Finish elimination
  std::cout << "\nCompleting elimination..." << std::endl;
  while (lhs.size() > 1) {
    SymbolicOptimization::gaussian_elimination(lhs, rhs, lhs.size() - 1);
  }

  std::cout << "\nFinal reduced system:" << std::endl;
  std::cout << "Reduced matrix:" << std::endl;
  print_lhs(lhs);

  std::cout << "Reduced right-hand side:" << std::endl;
  print_rhs(rhs);

  std::cout << "\nThe solution can now be obtained by back-substitution."
            << std::endl;
}

void run_evaluation_example() {
  using namespace Expression;
  using namespace NumericOptimization::Evaluation;

  print_header("Numerical Evaluation Example");

  // Create some symbolic expressions
  auto x = ExprFactory::variable("x");
  auto y = ExprFactory::variable("y");
  auto z = ExprFactory::variable("z");
  auto A = ExprFactory::matrix("A");
  auto B = ExprFactory::symmetric_matrix("B");
  auto Q = ExprFactory::symmetric_matrix("Q");
  auto c = ExprFactory::named_scalar("c");

  print_sub_header("Setting Up Evaluation Environment");
  std::cout << "Creating an environment with concrete values for variables:"
            << std::endl;

  // Create evaluation environment and populate with values
  Environment env;

  // Vector x = [1, 2, 3]
  env[x] = val_vector({1.0, 2.0, 3.0});
  std::cout << "  x = [1, 2, 3]" << std::endl;

  // Vector y = [4, 5, 6]
  env[y] = val_vector({4.0, 5.0, 6.0});
  std::cout << "  y = [4, 5, 6]" << std::endl;

  // Scalar c = 2.5
  env[c] = val_scalar(2.5);
  std::cout << "  c = 2.5" << std::endl;

  // Matrix A = [[1, 2, 3], [3, 4, 5], [5, 6, 7], [7, 8, 9]]
  ValMatrix matrix_a = {
      {1.0, 2.0, 3.0}, {3.0, 4.0, 5.0}, {5.0, 6.0, 7.0}, {7.0, 8.0, 9.0}};
  env[A] = matrix_a;
  std::cout << "  A = [[1, 2, 3], [3, 4, 5], [5, 6, 7], [7, 8, 9]]"
            << std::endl;

  // Symmetric matrix Q = [[1, 2, 3], [2, 4, 5], [3, 5, 6]]
  ValMatrix matrix_q = {{1.0, 2.0, 3.0}, {2.0, 4.0, 5.0}, {3.0, 5.0, 6.0}};
  env[Q] = matrix_q;
  std::cout << "  Q = [[1, 2, 3], [2, 4, 5], [3, 5, 6]]" << std::endl;

  // Example 1: Vector dot product x^T y
  print_sub_header("Example 1: Vector Dot Product");
  auto dot_product = ExprFactory::product({ExprFactory::transpose(x), y});
  std::cout << "Expression: " << dot_product->to_string() << std::endl;

  auto result1 = evaluate(dot_product, env);
  std::cout << "Result: ";
  if (is<ValScalar>(result1)) {
    std::cout << std::get<ValScalar>(result1) << std::endl;
  } else {
    std::cout << "Error: Expected scalar result but was " << result1.index()
              << std::endl;
  }
  std::cout << "Expected: 1*4 + 2*5 + 3*6 = 32" << std::endl;

  // Example 2: Scaled vector c * x
  print_sub_header("Example 2: Scaled Vector");
  auto scaled_vector = ExprFactory::product({c, x});
  std::cout << "Expression: " << scaled_vector->to_string() << std::endl;

  auto result2 = evaluate(scaled_vector, env);
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
  print_sub_header("Example 3: Quadratic Form");
  auto quadratic_form = ExprFactory::product({ExprFactory::transpose(x), Q, x});
  std::cout << "Expression: " << quadratic_form->to_string() << std::endl;

  auto result3 = evaluate(quadratic_form, env);
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
  print_sub_header("Example 4: Matrix-Vector Product");
  auto matrix_vector_product = ExprFactory::product({A, x});
  std::cout << "Expression: " << matrix_vector_product->to_string()
            << std::endl;
  ASSERT(is<ValVector>(env.at(x)));

  auto result4 = evaluate(matrix_vector_product, env);
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
  print_sub_header("Example 5: Complex Expression");
  // 0.5 * x^T Q x + c * y^T x
  auto complex_expr = ExprFactory::sum(
      {ExprFactory::product(
           {ExprFactory::number(0.5), ExprFactory::transpose(x), Q, x}),
       ExprFactory::product({c, ExprFactory::transpose(y), x})});

  std::cout << "Expression: " << complex_expr->to_string() << std::endl;
  std::cout << "Simplified: " << complex_expr->simplify()->to_string()
            << std::endl;

  auto result5 = evaluate(complex_expr, env);
  std::cout << "Result: ";
  if (is<ValScalar>(result5)) {
    std::cout << std::get<ValScalar>(result5) << std::endl;
  } else {
    std::cout << "Error: Expected scalar result" << std::endl;
  }
  std::cout << "Expected: 0.5 * 157 + 2.5 * 32 = 78.5 + 80 = 158.5"
            << std::endl;
}

void print_usage() {
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
      print_usage();
      return 0;
    } else if (arg == "-b" || arg == "--basic") {
      run_basic_examples();
      return 0;
    } else if (arg == "-o" || arg == "--optimize") {
      run_optimization_example();
      return 0;
    } else if (arg == "-e" || arg == "--evaluate") {
      run_evaluation_example();
      return 0;
    } else {
      std::cout << "Unknown option: " << arg << std::endl;
      print_usage();
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

  run_basic_examples();
  run_optimization_example();
  run_evaluation_example();

  std::cout << "\nDone!" << std::endl;

  return 0;
}