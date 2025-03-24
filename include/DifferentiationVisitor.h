#pragma once
#include "Expr.h"
#include "Helpers.h"

namespace Expression {
struct DifferentiationVisitor {
  explicit DifferentiationVisitor(const ExprPtr& var);

  ExprPtr operator()(const auto& x) const { return zero; }
  ExprPtr operator()(const Variable& x) const;
  ExprPtr operator()(const DiagonalMatrix& x) const;
  ExprPtr operator()(const Transpose& x) const;
  ExprPtr operator()(const Negate& x) const;
  ExprPtr operator()(const Invert& x) const;
  ExprPtr operator()(const Log& x) const;
  ExprPtr operator()(const Sum& x) const;
  ExprPtr operator()(const Product& x) const;

 private:
  const ExprPtr& var_;
};

}  // namespace Expression