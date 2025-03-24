#include "Expr.h"

#include <algorithm>
#include <iostream>
#include <numeric>
#include <ranges>
#include <set>

#include "Assert.h"
#include "DifferentiationVisitor.h"
#include "ExprFactory.h"
#include "Helpers.h"
#include "SimplificationVisitor.h"
#include "ToStringVisitor.h"

namespace Expression {
std::strong_ordering operator<=>(const ExprPtr& left, const ExprPtr& right) {
  if (left.get() == right.get()) {
    return std::strong_ordering::equal;
  }
  return (*left <=> *right);
}

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

bool operator==(const ExprPtr& left, const ExprPtr& right) {
  return (left <=> right) == std::strong_ordering::equal;
}

bool operator==(const Expr& left, const Expr& right) {
  return (left <=> right) == std::strong_ordering::equal;
}

std::ostream& operator<<(std::ostream& os, const ExprPtr& expr) {
  return os << expr->toString();
}

size_t ExprHash::operator()(const Expr& expr) const {
  size_t indexHash = std::hash<size_t>{}(expr.getImpl().index());
  size_t stringHash = std::hash<std::string>{}(expr.toExpressionString());

  return indexHash ^
         (stringHash + 0x9e3779b9 + (indexHash << 6) + (indexHash >> 2));
}

ExprPtr Expr::differentiate(const ExprPtr& var) const {
  if (!containsSubexpression(var)) {
    return zero;
  }
  return std::visit(DifferentiationVisitor(var), getImpl());
}

ExprPtr Expr::simplify(const bool distribute) const {
  auto expr = ExprFactory::asPtr(*this);
  auto changed = true;
  while (changed) {
    auto simplified = expr->simplifyOnce(distribute);
    changed = *simplified != *expr;
    std::swap(expr, simplified);
  }
  return expr;
}

ExprPtr Expr::simplifyOnce(const bool distribute) const {
  return std::visit(SimplificationVisitor(distribute), getImpl());
}

bool Expr::containsSubexpression(const ExprPtr& expr) const {
  if (*this == *expr) {
    return true;
  }
  return match(*this).with(([&expr](const auto& x)
                              requires UnaryType<decltype(x)>
                            { return x.child->containsSubexpression(expr); }),
                            [&expr](const auto& x)
                              requires NaryType<decltype(x)> {
                                return std::ranges::any_of(
                                    x.terms, [&expr](const auto& t) {
                                      return t->containsSubexpression(expr);
                                    });
                              },
                            [](const auto&) { return false; });
}

ExprPtr Expr::replaceSubexpression(const ExprPtr& expr,
                                   const ExprPtr& replacement) const {
  if (*this == *expr) {
    return replacement;
  }
  auto thisPtr = ExprFactory::asPtr(*this);
  return match(*this).with((
      [&](const auto& x)
        requires NaryType<decltype(x)> {
          auto terms = transform(x.terms, [&expr, &replacement](const auto& t) {
            return t->replaceSubexpression(expr, replacement);
          });
          using ExprType = std::decay_t<decltype(x)>;
          auto variant = ExprType{std::move(terms)};
          return ExprFactory::getExpr(std::move(variant));
        }),
      [&expr, &replacement](const auto& x)
        requires UnaryType<decltype(x)> {
          auto replaced = x.child->replaceSubexpression(expr, replacement);
          using ExprType = std::decay_t<decltype(x)>;
          auto variant = ExprType{std::move(replaced)};
          return ExprFactory::getExpr(std::move(variant));
        },
      [&](const auto& x) { return thisPtr; });
}

const Expr::ExprVariant& Expr::getImpl() const { return *impl_; }

std::string Expr::toString(const bool condensed) const {
  return std::visit(ToStringVisitor{condensed}, getImpl());
}

const std::string& Expr::toExpressionString() const {
  return expressionString_;
}

ExprPtr Expr::getLeadingOrEndingFactor(const bool leading) const {
  auto thisPtr = ExprFactory::asPtr(*this);
  return match(*this).with(
      [&](const Negate& x) {
        return x.child->getLeadingOrEndingFactor(leading);
      },
      [&](const Sum& x) {
        const auto termFactor =
            x.terms.front()->getLeadingOrEndingFactor(leading);
        if (std::ranges::all_of(x.terms, [&](const auto& t) {
              return t->getLeadingOrEndingFactor(leading) == termFactor;
            })) {
          return termFactor;
        }
        return thisPtr;
      },
      [&](const Product& x) {
        return (leading ? x.terms.front() : x.terms.back())
            ->getLeadingOrEndingFactor(leading);
      },
      [&](const auto& x) { return thisPtr; });
}

ExprPtr Expr::factorOut(const ExprPtr& factor, const bool leading) const {
  if (*factor == *this) {
    return unity;
  }
  ASSERT(getLeadingOrEndingFactor(leading) == factor);

  auto thisPtr = ExprFactory::asPtr(*this);
  return match(*this).with(
      [&](const Negate& x) {
        return ExprFactory::negate(x.child->factorOut(factor, leading));
      },
      [&](const Sum& x) {
        auto terms = transform(x.terms, [&](const auto& t) {
          return t->factorOut(factor, leading);
        });
        return ExprFactory::sum(std::move(terms));
      },
      [&](const Product& x) {
        auto terms = x.terms;
        for (size_t i = 0; i < terms.size(); ++i) {
          if (auto& t = terms[leading ? i : terms.size() - 1 - i];
              t->getLeadingOrEndingFactor(leading) == factor) {
            t = t->factorOut(factor, leading);
            return ExprFactory::product(std::move(terms));
          }
        }
        ASSERT(false);  // Should not be reached
        return thisPtr;
      },
      [&](const auto& x) {
        ASSERT(false);
        return thisPtr;
      });
}

double Expr::complexity() const {
  return match(*this).with(
      [](const Number&) { return 0.5; },
      ([](const auto& x)
         requires NamedNullaryType<decltype(x)> { return 1.0; }),
       ([](const auto& x)
          requires UnaryType<decltype(x)>
        { return 0.5 + x.child->complexity(); }),
        [](const auto& x)
          requires NaryType<decltype(x)> {
            return std::transform_reduce(
                x.terms.cbegin(), x.terms.cend(), 0.0, std::plus{},
                [](const auto& t) { return t->complexity(); });
          });
}

std::set<ExprPtr> Expr::getVariables() const {
  std::set<ExprPtr> variables;
  auto thisPtr = ExprFactory::asPtr(*this);
  match(*this).with([&](const Variable& x) { variables.insert(thisPtr); },
                    ([&](const auto& x)
                       requires UnaryType<decltype(x)> {
                         const auto v = x.child->getVariables();
                         variables.insert(v.begin(), v.end());
                       }),
                     [&](const auto& x)
                       requires NaryType<decltype(x)> {
                         for (const auto& t : x.terms) {
                           const auto v = t->getVariables();
                           variables.insert(v.begin(), v.end());
                         }
                       },
                     [](const auto&) {});
  return variables;
}
}  // namespace Expression