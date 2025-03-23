#pragma once
#include "Expression.h"
#include "Helpers.h"

namespace Expression {
struct DifferentiationVisitor {
  explicit DifferentiationVisitor(const Expr& var);

  Expr operator()(const auto& x) const { return zero; }
  Expr operator()(const Variable& x) const;
  Expr operator()(const DiagonalMatrix& x) const;
  Expr operator()(const Transpose& x) const;
  Expr operator()(const Negate& x) const;
  Expr operator()(const Invert& x) const;
  Expr operator()(const Log& x) const;
  Expr operator()(const Sum& x) const;
  Expr operator()(const Product& x) const;

 private:
  const Expr& var_;
};

}  // namespace Expression