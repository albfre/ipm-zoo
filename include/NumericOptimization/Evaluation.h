#pragma once
#include <functional>
#include <map>
#include <variant>
#include <vector>

#include "Expr.h"
#include "ExprFactory.h"

namespace NumericOptimization::Evaluation {
// Types of evaluation results
using ValScalar = double;
struct ValVector : std::vector<double> {
  using std::vector<double>::vector;
};
struct ValDiagMatrix : std::vector<double> {
  using std::vector<double>::vector;
};
using ValMatrix = std::vector<ValVector>;
using EvalResult = std::variant<ValScalar, ValVector, ValDiagMatrix, ValMatrix>;

using Environment = std::map<Expression::ExprPtr, Evaluation::EvalResult>;

std::vector<std::vector<double>> evaluate_matrix(
    const Expression::ExprPtr& expr, const Environment& env);
std::vector<double> evaluate_vector(const Expression::ExprPtr& expr,
                                    const Environment& env);

EvalResult evaluate(const Expression::ExprPtr& expr, const Environment& env);

// Functions for modifying EvalResults
EvalResult unary_op(const EvalResult& x,
                    const std::function<double(double)>& lambda);

EvalResult elementwise_op(const EvalResult& x, const EvalResult& y,
                          const std::function<double(double, double)>& lambda);

// Basic operations
EvalResult negate(const EvalResult& x);
EvalResult invert(const EvalResult& x);
EvalResult add(const EvalResult& x, const EvalResult& y);
EvalResult subtract(const EvalResult& x, const EvalResult& y);
EvalResult multiply(const EvalResult& x, const EvalResult& y);
EvalResult elementwise_multiply(const EvalResult& x, const EvalResult& y);
EvalResult elementwise_divide(const EvalResult& x, const EvalResult& y);

// Constructor functions
EvalResult val_scalar(double x);
EvalResult val_vector(const std::vector<double>& v);
EvalResult val_diag_matrix(const std::vector<double>& v);
EvalResult val_matrix(const std::vector<std::vector<double>>& m);
}  // namespace NumericOptimization::Evaluation