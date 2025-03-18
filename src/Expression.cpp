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

Expr Expr::simplifyOnce(const bool distribute) const {
  return std::visit(SimplificationVisitor(distribute), impl_);
}

Expr Expr::simplify(const bool distribute) const {
  auto expr = *this;
  auto changed = true;
  while (changed) {
    auto simplified =
        std::visit(SimplificationVisitor(distribute), expr.getImpl());
    changed = simplified != expr;
    std::swap(expr, simplified);
  }
  return expr;
}

bool Expr::containsSubexpression(const Expr& expr) const {
  if (*this == expr) {
    return true;
  }
  return match(
      *this, ([&expr](const auto& x)
                requires is_nary_v<decltype(x)> {
                  return std::ranges::any_of(x.terms, [&expr](const auto& t) {
                    return t.containsSubexpression(expr);
                  });
                }),
              [&expr](const auto& x)
                requires is_unary_v<decltype(x)>
              { return x.child->containsSubexpression(expr); },
              [](const auto&) { return false; });
}

Expr Expr::replaceSubexpression(const Expr& expr,
                                const Expr& replacement) const {
  if (*this == expr) {
    return replacement;
  }
  return match(
      *this,
      // Handle n-ary expressions
      ([&](const auto& x)
         requires is_nary_v<decltype(x)> {
           auto terms =
               transform(x.terms, [&expr, &replacement](const auto& t) {
                 return t.replaceSubexpression(expr, replacement);
               });
           return Expr(std::decay_t<decltype(x)>{std::move(terms)});
         }),

       // Handle unary expressions
       [&expr, &replacement](const auto& x)
         requires is_unary_v<decltype(x)>
       { return x.child->replaceSubexpression(expr, replacement); },

       // Handle all other expressions
       [&](const auto& x) { return *this; });
}

std::string Expr::toString(const bool condensed) const {
  return std::visit(ToStringVisitor{condensed}, impl_);
}

std::string Expr::toExpressionString() const {
  return std::visit(ToExpressionStringVisitor{}, impl_);
}

Expr Expr::getLeadingOrEndingFactor(const bool leading) const {
  return match(
      *this,
      [&](const Negate& x) {
        return x.child->getLeadingOrEndingFactor(leading);
      },
      [&](const Sum& x) {
        const auto termFactor =
            x.terms.front().getLeadingOrEndingFactor(leading);
        if (std::ranges::all_of(x.terms, [&](const auto& t) {
              return t.getLeadingOrEndingFactor(leading) == termFactor;
            })) {
          return termFactor;
        }
        return *this;
      },
      [&](const Product& x) {
        return (leading ? x.terms.front() : x.terms.back())
            .getLeadingOrEndingFactor(leading);
      },
      [&](const auto& x) { return *this; });
}

Expr Expr::factorOut(const Expr& factor, const bool leading) const {
  if (factor == *this) {
    return unity;
  }
  assert(getLeadingOrEndingFactor(leading) == factor);

  return match(
      *this,
      [&](const Negate& x) {
        return ExprFactory::negate(x.child->factorOut(factor, leading));
      },
      [&](const Sum& x) {
        auto terms = transform(x.terms, [&](const auto& t) {
          return t.factorOut(factor, leading);
        });
        return ExprFactory::sum(std::move(terms));
      },
      [&](const Product& x) {
        auto terms = x.terms;
        for (size_t i = 0; i < terms.size(); ++i) {
          if (auto& t = terms[leading ? i : terms.size() - 1 - i];
              t.getLeadingOrEndingFactor(leading) == factor) {
            t = t.factorOut(factor, leading);
            return ExprFactory::product(std::move(terms));
          }
        }
        assert(false);  // Should not be reached
        return *this;
      },
      [&](const auto& x) {
        assert(false);
        return *this;
      });
}

double Expr::complexity() const {
  return match(
      *this, [](const Number&) { return 0.5; },
      [](const auto& x)
          -> std::enable_if_t<is_named_nullary_v<decltype(x)>, double> {
        return 1.0;
      },
      [](const auto& x) -> std::enable_if_t<is_unary_v<decltype(x)>, double> {
        return 0.5 + x.child->complexity();
      },
      [](const auto& x) -> std::enable_if_t<is_nary_v<decltype(x)>, double> {
        return std::transform_reduce(
            x.terms.cbegin(), x.terms.cend(), 0.0, std::plus{},
            [](const auto& t) { return t.complexity(); });
      });
}

std::set<Expr> Expr::getVariables() const {
  std::set<Expr> variables;
  match(
      *this, [&](const Variable& x) { variables.insert(*this); },
      ([&](const auto& x)
         requires is_nary_v<decltype(x)> {
           for (const auto& t : x.terms) {
             const auto v = t.getVariables();
             variables.insert(v.begin(), v.end());
           }
         }),
       [&](const auto& x)
         requires is_unary_v<decltype(x)> {
           const auto v = x.child->getVariables();
           variables.insert(v.begin(), v.end());
         },
       [](const auto&) {});
  return variables;
}

Expr ExprFactory::number(const double value) { return Expr(Number{value}); }
Expr ExprFactory::namedConstant(const std::string& name) {
  return Expr(NamedConstant{name});
}
Expr ExprFactory::variable(const std::string& name) {
  return Expr(Variable{name});
}
Expr ExprFactory::matrix(const std::string& name) { return Expr(Matrix{name}); }
Expr ExprFactory::symmetricMatrix(const std::string& name) {
  return Expr(SymmetricMatrix{name});
}
Expr ExprFactory::diagonalMatrix(Expr expr) {
  return Expr(DiagonalMatrix{std::make_unique<Expr>(std::move(expr))});
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
Expr ExprFactory::sum(std::vector<Expr> terms) {
  return Expr(Sum{std::move(terms)});
}
Expr ExprFactory::product(std::vector<Expr> terms) {
  return Expr(Product{std::move(terms)});
}
}  // namespace Expression