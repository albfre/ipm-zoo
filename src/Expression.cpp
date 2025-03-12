#include "Expression.h"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <map>
#include <numeric>
#include <ranges>
#include <set>
#include <sstream>

#include "DifferentiationVisitor.h"
#include "Helpers.h"
#include "SimplificationVisitor.h"
#include "ToStringVisitor.h"

namespace Expression {
std::strong_ordering operator<=>(const Expr& left, const Expr& right) {
  // Order lexicographically, first by type, then by the expression string
  if (&left == &right) {
    return std::strong_ordering::equal;
  }
  if (left.getImpl().index() != right.getImpl().index()) {
    return left.getImpl().index() <=> right.getImpl().index();
  }
  return left.toExpressionString() <=> right.toExpressionString();
}

bool operator==(const Expr& left, const Expr& right) {
  return (left <=> right) == std::strong_ordering::equal;
}

Expr::Expr(const Expr& other)
    : impl_(std::visit(
          [](const auto& val) -> ExprVariant {
            using T = std::decay_t<decltype(val)>;
            if constexpr (std::is_base_of_v<UnaryExpr, T>) {
              return T{std::make_unique<Expr>(*val.child)};
            } else {
              return val;
            }
          },
          other.impl_)) {}

Expr& Expr::operator=(const Expr& other) {
  if (this != &other) {
    Expr temp(other);
    std::swap(impl_, temp.impl_);
  }
  return *this;
}

Expr Expr::differentiate(const Expr& var) const {
  if (!containsSubexpression(var)) {
    return zero;
  }
  return std::visit(DifferentiationVisitor(var), impl_);
}

Expr Expr::simplify(const bool distribute) const {
  auto expr = *this;
  auto changed = true;
  while (changed) {
    auto simplified = std::visit(SimplificationVisitor(distribute), impl_);
    changed = simplified != expr;
    std::swap(expr, simplified);
  }
  return expr;
}

bool Expr::containsSubexpression(const Expr& expr) const {
  if (*this == expr) {
    return true;
  }
  return match(impl_, [&expr](const auto& x) {
    using T = std::decay_t<decltype(x)>;
    if constexpr (std::is_base_of_v<NaryExpr, T>) {
      return std::ranges::any_of(x.terms, [&expr](const auto& t) {
        return t.containsSubexpression(expr);
      });
    } else if constexpr (std::is_base_of_v<UnaryExpr, T>) {
      return x.child->containsSubexpression(expr);
    } else {
      return false;
    }
  });
}

std::string Expr::toString(const bool condensed) const {
  return std::visit(ToStringVisitor{condensed}, impl_);
}

std::string Expr::toExpressionString() const {
  return std::visit(ToExpressionStringVisitor{}, impl_);
}

double Expr::complexity_() const {
  return match(
      impl_, [](const Number&) { return 0.5; },
      [](const NamedConstant&) { return 1.0; },
      [](const Matrix&) { return 1.0; },
      [](const SymmetricMatrix&) { return 1.0; },
      [](const Variable&) { return 1.0; },
      [](const Transpose& x) { return 0.5 + x.child->complexity_(); },
      [](const Negate& x) { return 0.5 + x.child->complexity_(); },
      [](const Invert& x) { return 0.5 + x.child->complexity_(); },
      [](const Log& x) { return 1.0 + x.child->complexity_(); },
      [](const auto& x) {
        if constexpr (std::is_base_of_v<NaryExpr, std::decay_t<decltype(x)>>) {
          return std::transform_reduce(
              x.terms.cbegin(), x.terms.cend(), 0.0, std::plus{},
              [](const auto& t) { return t.complexity_(); });
        } else {
          static_assert(always_false_v<decltype(x)>);
        }
      });
}

Expr ExprFactory::number(const double value) { return Expr(Number{value}); }
Expr ExprFactory::namedConstant(const std::string& name) {
  return Expr(NamedConstant{name});
}
Expr ExprFactory::matrix(const std::string& name) { return Expr(Matrix{name}); }
Expr ExprFactory::symmetricMatrix(const std::string& name) {
  return Expr(SymmetricMatrix{name});
}
Expr ExprFactory::variable(const std::string& name) {
  return Expr(Variable{name});
}
Expr ExprFactory::diagonalMatrix(Expr expr) {
  return Expr(ExprType::DiagonalMatrix, {std::move(expr)});
}
Expr ExprFactory::transpose(Expr expr) {
  return Expr(Transpose{std::make_unique<Expr>(std::move(expr))});
}
Expr ExprFactory::negate(Expr expr) {
  return Expr(Negate{std::make_unique<Expr>(std::move(expr))});
}
Expr ExprFactory::invert(Expr expr) {
  return Expr(Invert{std::make_unique<Expr>(std::move(expr))});
}
Expr ExprFactory::log(Expr expr) {
  return Expr(Log{std::make_unique<Expr>(std::move(expr))});
}
Expr ExprFactory::product(std::vector<Expr> terms) {
  return Expr(Product{std::move(terms)});
}
Expr ExprFactory::sum(std::vector<Expr> terms) {
  return Expr(Sum{std::move(terms)});
}
}  // namespace Expression