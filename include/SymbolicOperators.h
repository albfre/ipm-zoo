#pragma once
#include "Expr.h"
#include "ExprFactory.h"

namespace Expression {
inline ExprPtr operator+(const ExprPtr& lhs, const ExprPtr& rhs) {
  return ExprFactory::sum({lhs, rhs});
}

inline ExprPtr operator*(const ExprPtr& lhs, const ExprPtr& rhs) {
  return ExprFactory::product({lhs, rhs});
}

inline ExprPtr operator-(const ExprPtr& lhs, const ExprPtr& rhs) {
  return ExprFactory::sum({lhs, ExprFactory::negate(rhs)});
}

inline ExprPtr operator-(const ExprPtr& expr) {
  return ExprFactory::negate(expr);
}
}  // namespace Expression