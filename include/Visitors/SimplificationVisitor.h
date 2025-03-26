#pragma once
#include "Expr.h"
#include "Utils/Helpers.h"

namespace Expression {
struct SimplificationVisitor {
  explicit SimplificationVisitor(bool distribute);

  ExprPtr operator()(const auto& x) const { return ExprFactory::get_expr(x); }
  ExprPtr operator()(const DiagonalMatrix& x);
  ExprPtr operator()(const Transpose& x) const;
  ExprPtr operator()(const Negate& x);
  ExprPtr operator()(const Invert& x) const;
  ExprPtr operator()(const Log& x) const;
  ExprPtr operator()(const Sum& x) const;
  ExprPtr operator()(const Product& x) const;

 private:
  bool distribute_ = true;
  template <typename T>
  void erase_canceling_(std::vector<ExprPtr>& terms,
                        const ExprPtr& replacement) const;

  template <typename T>
    requires NaryType<T>
  void associative_transformation_(std::vector<ExprPtr>& terms) const;
};
}  // namespace Expression