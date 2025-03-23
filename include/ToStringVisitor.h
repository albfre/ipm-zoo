#pragma once
#include <ranges>
#include <sstream>

#include "Expr.h"
#include "Helpers.h"

namespace Expression {
struct ToStringVisitor {
  ToStringVisitor(bool condensed);
  std::string operator()(const Number& x) const;
  std::string operator()(const auto& x) const {
    if constexpr (NamedNullaryType<decltype(x)>) {
      return x.name;
    } else {
      static_assert(always_false_v<decltype(x)>);
    }
  }
  std::string operator()(const DiagonalMatrix& x) const;
  std::string operator()(const Transpose& x) const;
  std::string operator()(const Negate& x) const;
  std::string operator()(const Invert& x) const;
  std::string operator()(const Log& x) const;
  std::string operator()(const Sum& x) const;
  std::string operator()(const Product& x) const;

 private:
  bool condensed_ = false;
};

struct ToExpressionStringVisitor {
  std::string operator()(const Number& x) const;
  std::string operator()(const NamedScalar& x) const;
  std::string operator()(const NamedVector& x) const;
  std::string operator()(const Variable& x) const;
  std::string operator()(const Matrix& x) const;
  std::string operator()(const SymmetricMatrix& x) const;
  std::string operator()(const DiagonalMatrix& x) const;
  std::string operator()(const Transpose& x) const;
  std::string operator()(const Negate& x) const;
  std::string operator()(const Invert& x) const;
  std::string operator()(const Log& x) const;
  std::string operator()(const Sum& x) const;
  std::string operator()(const Product& x) const;

 private:
  std::string termsToString_(const std::vector<Expr>& terms) const;
};
}  // namespace Expression