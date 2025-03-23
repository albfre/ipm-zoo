#pragma once
#include <functional>
#include <map>
#include <variant>
#include <vector>

#include "Expr.h"

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

// Functions for modifying EvalResults
EvalResult unaryOp(const EvalResult& x,
                   const std::function<double(double)>& lambda);

EvalResult elementwiseOp(const EvalResult& x, const EvalResult& y,
                         const std::function<double(double, double)>& lambda);

// Basic operations
EvalResult negate(const EvalResult& x);
EvalResult invert(const EvalResult& x);
EvalResult add(const EvalResult& x, const EvalResult& y);
EvalResult subtract(const EvalResult& x, const EvalResult& y);
EvalResult multiply(const EvalResult& x, const EvalResult& y);
EvalResult elementwiseMultiply(const EvalResult& x, const EvalResult& y);
EvalResult elementwiseDivide(const EvalResult& x, const EvalResult& y);

// Constructor functions
EvalResult valScalar(double x);
EvalResult valVector(const std::vector<double>& v);
EvalResult valDiagMatrix(const std::vector<double>& v);
}  // namespace Evaluation