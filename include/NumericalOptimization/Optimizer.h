#pragma once
#include <map>
#include <vector>

#include "Evaluation.h"
#include "Expr.h"
#include "SymbolicOptimization.h"

namespace NumericalOptimization {
using Matrix = std::vector<std::vector<double>>;
using Vector = std::vector<double>;

class Optimizer {
 public:
  Optimizer(Evaluation::Environment& env,
            const SymbolicOptimization::OptimizationExpressions&
                optimization_expressions,
            SymbolicOptimization::NewtonSystem newton_system);

  void solve();

 private:
  void solve_indefinite_();
  void solve_quasi_definite_();

  void update_variables_(double alpha, Matrix variables, const Matrix& delta);
  void vector_plus_eq_scalar_times_vector_(std::vector<double>& x,
                                           const double s,
                                           const std::vector<double>& y);
  double get_residual_norm_(const std::vector<Expression::ExprPtr>& rhs);
  double get_mu_(const std::vector<Expression::ExprPtr>& rhs);
  double get_max_step_(Matrix variables, Matrix delta);
  std::vector<std::vector<double>> compute_search_direction_(
      const SymbolicOptimization::NewtonSystem& newton_system, const Matrix& L,
      const std::vector<double>& D);
  double dot_(const std::vector<double>& x, const std::vector<double>& y);
  Matrix get_as_matrix_(const std::vector<std::vector<Expression::ExprPtr>>& v);
  Vector get_as_vector_(const std::vector<Expression::ExprPtr>& v);

  template <typename T>
  std::vector<T> eval_vector_of_expressions_(
      const std::vector<Expression::ExprPtr>& v);

  template <typename T>
  std::vector<std::vector<T>> eval_matrix_of_expressions_(
      const std::vector<std::vector<Expression::ExprPtr>>& v);

  std::vector<double> concatenate_vectors_(
      const std::vector<std::vector<double>>& vectors);

  Matrix concatenate_matrices_(
      const std::vector<std::vector<Matrix>>& matrices);

  Evaluation::Environment& env_;
  SymbolicOptimization::OptimizationExpressions optimization_expressions_;
  SymbolicOptimization::NewtonSystem newton_system_;
  SymbolicOptimization::ShorthandRhs shorthand_rhs_;
  SymbolicOptimization::NewtonSystem augmented_system_;
  SymbolicOptimization::NewtonSystem normal_equations_;
  Expression::ExprPtr objective_;
  std::map<Expression::ExprPtr, size_t> vector_sizes_;
  std::map<Expression::ExprPtr, size_t> variable_index_;
  std::map<Expression::ExprPtr, size_t> delta_variable_index_;
};

}  // namespace NumericalOptimization