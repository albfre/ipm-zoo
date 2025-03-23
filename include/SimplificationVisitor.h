#pragma once

#include "Expression.h"
#include "Helpers.h"

namespace Expression {
struct SimplificationVisitor {
  explicit SimplificationVisitor(bool distribute);

  Expr operator()(const auto& x) const { return Expr(x); }
  Expr operator()(const DiagonalMatrix& x);
  Expr operator()(const Transpose& x) const;
  Expr operator()(const Negate& x);
  Expr operator()(const Invert& x) const;
  Expr operator()(const Log& x) const;
  Expr operator()(const Sum& x) const;
  Expr operator()(const Product& x) const;

 private:
  bool distribute_ = true;
  template <typename T>
  void eraseCanceling_(std::vector<Expr>& terms, const Expr& replacement) const;

  template <typename T>
    requires NaryType<T>
  void associativeTransformation_(std::vector<Expr>& terms) const;
};
}  // namespace Expression