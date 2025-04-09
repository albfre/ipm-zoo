#include "Expr.h"

#include <algorithm>
#include <numeric>

#include "ExprFactory.h"
#include "Utils/Assert.h"
#include "Utils/Helpers.h"
#include "Visitors/DifferentiationVisitor.h"
#include "Visitors/SimplificationVisitor.h"
#include "Visitors/ToStringVisitor.h"

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
  if (left.get_impl().index() != right.get_impl().index()) {
    return left.get_impl().index() <=> right.get_impl().index();
  }

  return left.to_expression_string() <=> right.to_expression_string();
}

bool operator==(const ExprPtr& left, const ExprPtr& right) {
  return (left <=> right) == std::strong_ordering::equal;
}

bool operator==(const Expr& left, const Expr& right) {
  return (left <=> right) == std::strong_ordering::equal;
}

std::ostream& operator<<(std::ostream& os, const ExprPtr& expr) {
  return os << expr->to_string();
}

size_t ExprHash::operator()(const Expr& expr) const {
  size_t index_hash = std::hash<size_t>{}(expr.get_impl().index());
  size_t string_hash = std::hash<std::string>{}(expr.to_expression_string());

  return index_hash ^
         (string_hash + 0x9e3779b9 + (index_hash << 6) + (index_hash >> 2));
}

ExprPtr Expr::differentiate(const ExprPtr& var) const {
  if (!contains_subexpression(var)) {
    return zero;
  }
  return std::visit(DifferentiationVisitor(var), get_impl());
}

ExprPtr Expr::simplify(const bool distribute) const {
  auto expr = ExprFactory::as_ptr(*this);
  auto changed = true;
  while (changed) {
    auto simplified = expr->simplify_once(distribute);
    changed = *simplified != *expr;
    std::swap(expr, simplified);
  }
  return expr;
}

ExprPtr Expr::simplify_once(const bool distribute) const {
  return std::visit(SimplificationVisitor(distribute), get_impl());
}

bool Expr::contains_subexpression(const ExprPtr& expr) const {
  if (*this == *expr) {
    return true;
  }
  return match(*this).with(([&expr](const auto& x)
                              requires UnaryType<decltype(x)>
                            { return x.child->contains_subexpression(expr); }),
                            [&expr](const auto& x)
                              requires NaryType<decltype(x)> {
                                return std::ranges::any_of(
                                    x.terms, [&expr](const auto& t) {
                                      return t->contains_subexpression(expr);
                                    });
                              },
                            [](const auto&) { return false; });
}

ExprPtr Expr::replace_subexpression(const ExprPtr& expr,
                                    const ExprPtr& replacement) const {
  if (*this == *expr) {
    return replacement;
  }
  auto this_ptr = ExprFactory::as_ptr(*this);
  return match(*this).with((
      [&](const auto& x)
        requires NaryType<decltype(x)> {
          auto terms = transform(x.terms, [&expr, &replacement](const auto& t) {
            return t->replace_subexpression(expr, replacement);
          });
          using ExprType = std::decay_t<decltype(x)>;
          auto variant = ExprType{std::move(terms)};
          return ExprFactory::get_expr(std::move(variant));
        }),
      [&expr, &replacement](const auto& x)
        requires UnaryType<decltype(x)> {
          auto replaced = x.child->replace_subexpression(expr, replacement);
          using ExprType = std::decay_t<decltype(x)>;
          auto variant = ExprType{std::move(replaced)};
          return ExprFactory::get_expr(std::move(variant));
        },
      [&](const auto& x) { return this_ptr; });
}

const Expr::ExprVariant& Expr::get_impl() const { return *impl_; }

std::string Expr::to_string(const bool condensed) const {
  return std::visit(ToStringVisitor{condensed}, get_impl());
}

const std::string& Expr::to_expression_string() const {
  return expression_string_;
}

ExprPtr Expr::get_leading_or_ending_factor(const bool leading) const {
  auto this_ptr = ExprFactory::as_ptr(*this);
  return match(*this).with(
      [&](const Negate& x) {
        return x.child->get_leading_or_ending_factor(leading);
      },
      [&](const Sum& x) {
        const auto term_factor =
            x.terms.front()->get_leading_or_ending_factor(leading);
        if (std::ranges::all_of(x.terms, [&](const auto& t) {
              return t->get_leading_or_ending_factor(leading) == term_factor;
            })) {
          return term_factor;
        }
        return this_ptr;
      },
      [&](const Product& x) {
        return (leading ? x.terms.front() : x.terms.back())
            ->get_leading_or_ending_factor(leading);
      },
      [&](const auto& x) { return this_ptr; });
}

ExprPtr Expr::factor_out(const ExprPtr& factor, const bool leading) const {
  if (*factor == *this) {
    return unity;
  }
  ASSERT(get_leading_or_ending_factor(leading) == factor);

  auto this_ptr = ExprFactory::as_ptr(*this);
  return match(*this).with(
      [&](const Negate& x) {
        return ExprFactory::negate(x.child->factor_out(factor, leading));
      },
      [&](const Sum& x) {
        auto terms = transform(x.terms, [&](const auto& t) {
          return t->factor_out(factor, leading);
        });
        return ExprFactory::sum(std::move(terms));
      },
      [&](const Product& x) {
        auto terms = x.terms;
        for (size_t i = 0; i < terms.size(); ++i) {
          if (auto& t = terms[leading ? i : terms.size() - 1 - i];
              t->get_leading_or_ending_factor(leading) == factor) {
            t = t->factor_out(factor, leading);
            return ExprFactory::product(std::move(terms));
          }
        }
        ASSERT(false);  // Should not be reached
        return this_ptr;
      },
      [&](const auto& x) {
        ASSERT(false);
        return this_ptr;
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

std::set<ExprPtr> Expr::get_variables() const {
  std::set<ExprPtr> variables;
  auto this_ptr = ExprFactory::as_ptr(*this);
  match(*this).with([&](const Variable& x) { variables.insert(this_ptr); },
                    ([&](const auto& x)
                       requires UnaryType<decltype(x)> {
                         const auto v = x.child->get_variables();
                         variables.insert(v.begin(), v.end());
                       }),
                     [&](const auto& x)
                       requires NaryType<decltype(x)> {
                         for (const auto& t : x.terms) {
                           const auto v = t->get_variables();
                           variables.insert(v.begin(), v.end());
                         }
                       },
                     [](const auto&) {});
  return variables;
}
}  // namespace Expression