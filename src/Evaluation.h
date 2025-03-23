#pragma once
#include <functional>
#include <map>
#include <variant>
#include <vector>

#include "Expression.h"

namespace Evaluation {
// Type definitions for evaluation results
using ValScalar = double;
struct ValVector : std::vector<double> {
  using std::vector<double>::vector;
};
struct ValDiagMatrix : std::vector<double> {
  using std::vector<double>::vector;
};
using ValMatrix = std::vector<ValVector>;
using EvalResult = std::variant<ValScalar, ValVector, ValDiagMatrix, ValMatrix>;

// Environment mapping expressions to their evaluation results
using Environment = std::map<Expression::Expr, EvalResult>;

// Main evaluation function
EvalResult evaluate(const Expression::Expr& expr, Environment& env);

EvalResult unaryOp(const EvalResult& x,
                   const std::function<double(double)>& lambda);

// Template for element-wise binary operations on evaluation results
EvalResult elementwiseOp(const EvalResult& x, const EvalResult& y,
                         const std::function<double(double, double)>& lambda);

// Vector and matrix operations
ValScalar dot(const ValVector& x, const ValVector& y);
EvalResult product(const EvalResult& x, const EvalResult& y,
                   std::vector<EvalResult>& unhandled);

// Basic operations
EvalResult negate(const EvalResult& x);
EvalResult invert(const EvalResult& x);
EvalResult add(const EvalResult& x, const EvalResult& y);
EvalResult subtract(const EvalResult& x, const EvalResult& y);
EvalResult elementwiseMultiply(const EvalResult& x, const EvalResult& y);
EvalResult elementwiseDivide(const EvalResult& x, const EvalResult& y);

// Constructor functions
EvalResult valScalar(double x);
EvalResult valVector(const std::vector<double>& v);
EvalResult valDiagMatrix(const std::vector<double>& v);
}  // namespace Evaluation