#include "Expression.h"

#include <algorithm>
#include <iostream>
#include <map>
#include <numeric>
#include <ranges>
#include <set>
#include <sstream>

#include "Assert.h"
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

std::ostream& operator<<(std::ostream& os, const Expr& expr) {
  return os << expr.toString();
}

size_t ExprHash::operator()(const Expr& expr) const {
  size_t indexHash = std::hash<size_t>{}(expr.getImpl().index());
  size_t stringHash = std::hash<std::string>{}(expr.toExpressionString());

  return indexHash ^
         (stringHash + 0x9e3779b9 + (indexHash << 6) + (indexHash >> 2));
}

Expr Expr::differentiate(const Expr& var) const {
  if (!containsSubexpression(var)) {
    return zero;
  }
  return std::visit(DifferentiationVisitor(var), getImpl());
}

Expr Expr::simplify(const bool distribute) const {
  auto expr = *this;
  auto changed = true;
  while (changed) {
    auto simplified = expr.simplifyOnce(distribute);
    changed = simplified != expr;
    std::swap(expr, simplified);
  }
  return expr;
}

Expr Expr::simplifyOnce(const bool distribute) const {
  return std::visit(SimplificationVisitor(distribute), getImpl());
}

bool Expr::containsSubexpression(const Expr& expr) const {
  if (*this == expr) {
    return true;
  }
  return match(*this).with(([&expr](const auto& x)
                              requires is_nary_v<decltype(x)> {
                                return std::ranges::any_of(
                                    x.terms, [&expr](const auto& t) {
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
  return match(*this).with((
      [&](const auto& x)
        requires is_nary_v<decltype(x)> {
          using ExprType = std::decay_t<decltype(x)>;
          auto terms = transform(x.terms, [&expr, &replacement](const auto& t) {
            return t.replaceSubexpression(expr, replacement);
          });
          return Expr(ExprType{std::move(terms)});
        }),
      [&expr, &replacement](const auto& x)
        requires is_unary_v<decltype(x)> {
          auto replaced = x.child->replaceSubexpression(expr, replacement);
          using ExprType = std::decay_t<decltype(x)>;
          return Expr(ExprType{std::make_shared<Expr>(std::move(replaced))});
        },
      [&](const auto& x) { return *this; });
}

const Expr::ExprVariant& Expr::getImpl() const { return *impl_; }

std::string Expr::toString(const bool condensed) const {
  return std::visit(ToStringVisitor{condensed}, getImpl());
}

std::string Expr::toExpressionString() const {
  return std::visit(ToExpressionStringVisitor{}, getImpl());
}

Expr Expr::getLeadingOrEndingFactor(const bool leading) const {
  return match(*this).with(
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
  ASSERT(getLeadingOrEndingFactor(leading) == factor);

  return match(*this).with(
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
        ASSERT(false);  // Should not be reached
        return *this;
      },
      [&](const auto& x) {
        ASSERT(false);
        return *this;
      });
}

double Expr::complexity() const {
  return match(*this).with(
      [](const Number&) { return 0.5; },
      ([](const auto& x)
         requires is_named_nullary_v<decltype(x)> { return 1.0; }),
       ([](const auto& x)
          requires is_unary_v<decltype(x)>
        { return 0.5 + x.child->complexity(); }),
        [](const auto& x)
          requires is_nary_v<decltype(x)> {
            return std::transform_reduce(
                x.terms.cbegin(), x.terms.cend(), 0.0, std::plus{},
                [](const auto& t) { return t.complexity(); });
          });
}

std::set<Expr> Expr::getVariables() const {
  std::set<Expr> variables;
  match(*this).with([&](const Variable& x) { variables.insert(*this); },
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

Expr ExprFactory::number(const double value) { return Expr(Number(value)); }
Expr ExprFactory::namedScalar(const std::string& name) {
  return Expr(NamedScalar(name));
}
Expr ExprFactory::namedVector(const std::string& name) {
  return Expr(NamedVector(name));
}
Expr ExprFactory::variable(const std::string& name) {
  return Expr(Variable(name));
}
Expr ExprFactory::matrix(const std::string& name) { return Expr(Matrix(name)); }
Expr ExprFactory::symmetricMatrix(const std::string& name) {
  return Expr(SymmetricMatrix(name));
}
Expr ExprFactory::diagonalMatrix(Expr expr) {
  return Expr(DiagonalMatrix(std::make_shared<Expr>(std::move(expr))));
}
Expr ExprFactory::transpose(Expr expr) {
  return Expr(Transpose(std::make_shared<Expr>(std::move(expr))));
}
Expr ExprFactory::negate(Expr expr) {
  return Expr(Negate(std::make_shared<Expr>(std::move(expr))));
}
Expr ExprFactory::invert(Expr expr) {
  return Expr(Invert(std::make_shared<Expr>(std::move(expr))));
}
Expr ExprFactory::log(Expr expr) {
  return Expr(Log(std::make_shared<Expr>(std::move(expr))));
}
Expr ExprFactory::sum(std::vector<Expr> terms) {
  if (terms.empty()) {
    return zero;
  }
  if (terms.size() == 1) {
    return std::move(terms[0]);
  }
  return Expr(Sum(std::move(terms)));
}
Expr ExprFactory::product(std::vector<Expr> terms) {
  return Expr(Product(std::move(terms)));
}
}  // namespace Expression