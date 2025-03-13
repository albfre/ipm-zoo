#pragma once
#include <algorithm>
#include <numeric>
#include <ranges>

#include "Expression.h"
#include "Helpers.h"

namespace Expression {
struct SimplificationVisitor {
  bool distribute = true;
  explicit SimplificationVisitor(const bool distribute)
      : distribute(distribute) {}

  Expr operator()(const Number& x) const { return Expr(x); }
  Expr operator()(const NamedConstant& x) const { return Expr(x); }
  Expr operator()(const Variable& x) const { return Expr(x); }
  Expr operator()(const Matrix& x) const { return Expr(x); }
  Expr operator()(const SymmetricMatrix& x) const { return Expr(x); }
  Expr operator()(const DiagonalMatrix& x) {
    auto child = x.child->simplify(distribute);
    if (child == zero || child == unity) {
      return child;
    }
    return ExprFactory::diagonalMatrix(std::move(child));
  }

  Expr operator()(const Transpose& x) const {
    auto child = x.child->simplify(distribute);
    if (child == zero || child == unity) {
      return child;  // 0^T = 0, 1^T = 1
    }
    return match(
        child, [](const Transpose& x) { return *x.child; },  // (x^T)^T = x
        [&](const SymmetricMatrix& x) { return child; },
        [&](const Invert& x) { return child; },
        [this](const Negate& x) {
          // (-x)^T  = -x^T
          return ExprFactory::negate(ExprFactory::transpose(*x.child));
        },
        [this](const Product& x) {
          // (xyz)^T = z^T y^T x^T
          auto terms = transform(x.terms, [this](const auto& t) {
            return ExprFactory::transpose(t);
          });
          std::ranges::reverse(terms);
          return ExprFactory::product(std::move(terms));
        },
        [&](const auto& x) {
          return ExprFactory::transpose(std::move(child));
        });
  }

  Expr operator()(const Negate& x) {
    auto child = x.child->simplify(distribute);
    if (child == zero) {  // -0 = 0
      return zero;
    }
    return match(
        child, [](const Number& x) { return ExprFactory::number(-x.value); },
        [](const Negate& x) { return *x.child; },  // negate(negate(x)) = x
        [&](const Product& x) {
          const auto it = std::ranges::find_if(x.terms, is<Negate>);
          if (it != x.terms.end()) {
            auto terms = x.terms;
            terms[std::distance(x.terms.begin(), it)] =
                *std::get<Negate>(it->getImpl()).child;
            return ExprFactory::product(std::move(terms));
          }
          return ExprFactory::negate(std::move(child));
        },
        [&](const Sum& x) {
          if (static_cast<size_t>(std::ranges::count_if(x.terms, is<Negate>)) >
              x.terms.size() / 2) {
            auto terms = x.terms;
            for (auto& t : terms) {
              t = match(
                  t, [](const Negate& y) { return *y.child; },
                  [&](const auto&) { return ExprFactory::negate(t); });
            }
            return ExprFactory::sum(std::move(terms));
          }
          return ExprFactory::negate(std::move(child));
        },
        [&](const auto& x) { return ExprFactory::negate(std::move(child)); });
  }

  Expr operator()(const Invert& x) const {
    auto child = x.child->simplify(distribute);
    return match(
        child,
        [](const Invert& y) { return *y.child; },  // invert(invert(x)) = x
        [](const Negate& y) {
          // invert(negate(x)) = negate(invert(x))
          return ExprFactory::negate(ExprFactory::invert(*y.child));
        },
        [](const Product& y) {
          // invert(x * y * z) = invert(z) * invert(y) * invert(x)
          auto terms = transform(
              y.terms, [](const auto& t) { return ExprFactory::invert(t); });
          std::ranges::reverse(terms);
          return ExprFactory::product(std::move(terms));
        },
        [&](const auto&) { return ExprFactory::invert(std::move(child)); });
  }

  Expr operator()(const Log& x) const {
    return ExprFactory::log(x.child->simplify(distribute));
  }

  Expr operator()(const Sum& x) const {
    // Recursive simplification
    auto terms = transform(
        x.terms, [this](const auto& t) { return t.simplify(distribute); });

    // Associative transformation ((x + y) + z = x + y + z)
    std::vector<Expr> newTerms;
    newTerms.reserve(terms.size());
    for (auto& t : terms) {
      match(
          t,
          [&](const Sum& x) {
            newTerms.insert(newTerms.end(), x.terms.begin(), x.terms.end());
          },
          [&](auto& x) {
            if constexpr (std::is_same_v<Negate, std::decay_t<decltype(x)>>) {
              auto& child = x.child;
              if constexpr (std::is_same_v<Sum,
                                           std::decay_t<decltype(child)>>) {
                newTerms.reserve(newTerms.size() + child.terms.size());
                for (auto& ct : child.terms) {
                  newTerms.emplace_back(ExprFactory::negate(std::move(ct)));
                }
                return;
              }
            }
            newTerms.emplace_back(std::move(t));
          });
    }
    terms = std::move(newTerms);

    // Numerical transformation (1 + x + 2 = 3 + x)
    if (std::ranges::any_of(terms, is<Number>)) {
      const auto value = std::accumulate(
          terms.cbegin(), terms.cend(), 0.0, [](const double s, const auto& t) {
            return s + match(
                           t, [](const Number& x) { return x.value; },
                           [](const auto&) { return 0.0; });
          });
      std::erase_if(terms, is<Number>);
      terms.push_back(ExprFactory::number(value));
    }

    // Distributive transformation (x + y + 1.3x = 2.3x + y)
    for (size_t i = 0; i < terms.size(); ++i) {
      const auto term = terms.at(i);
      const auto negTerm = ExprFactory::negate(terms.at(i));
      const auto isNumberTimesTerm = [&term](const auto& t) -> bool {
        return match(
            t,
            [&term](const Product& x) {
              return x.terms.size() == 2 && is<Number>(x.terms.front()) &&
                     x.terms.back() == term;
            },
            [](const auto&) { return false; });
      };
      const auto isTerm = [&](const auto& t) {
        return t == term || t == negTerm || isNumberTimesTerm(t);
      };
      if (std::ranges::any_of(terms, isTerm)) {
        const auto value = std::accumulate(
            terms.cbegin(), terms.cend(), 0.0,
            [&term, &negTerm, &isNumberTimesTerm](const auto s,
                                                  const auto& t) -> double {
              return s + (t == term      ? 1.0
                          : t == negTerm ? -1.0
                          : isNumberTimesTerm(t)
                              ? std::get<Number>(std::get<Product>(t.getImpl())
                                                     .terms.front()
                                                     .getImpl())
                                    .value
                              : 0.0);
            });
        std::erase_if(terms, isTerm);
        terms.push_back(
            ExprFactory::product({ExprFactory::number(value), term}));
      }
    }

    // Identity transformations (x + 0 = x)
    std::erase_if(terms, [&](const auto& t) { return t == zero; });

    if (terms.empty()) {
      return zero;
    }

    // Canceling transformation (x - x = 0)
    eraseCanceling_<Negate>(terms, zero);

    // Commutative transformation (z + y + x = x + y + z)
    std::ranges::sort(terms);

    auto simplified = ExprFactory::sum(terms);

    // Check whether extracting common factors leads to a shorter expression
    // (xy + xz + xw = x(y + z + w)
    if (distribute) {
      for (const auto leading : std::vector{true, false}) {
        std::map<Expr, size_t> numOccurrencesOfFactors;
        std::vector<Expr> factorPerTerm;
        factorPerTerm.reserve(terms.size());
        for (auto& t : terms) {
          factorPerTerm.push_back(t.getLeadingOrEndingFactor(leading));
          ++numOccurrencesOfFactors[factorPerTerm.back()];
        }
        std::vector<std::pair<Expr, size_t>> sorted(
            numOccurrencesOfFactors.begin(), numOccurrencesOfFactors.end());
        std::ranges::sort(sorted, {}, &std::pair<Expr, size_t>::second);

        while (!sorted.empty()) {
          const auto& [factor, numOccurrences] = sorted.back();
          if (numOccurrences < 2) {
            break;
          }
          std::vector<Expr> factoredTerms;
          std::vector<Expr> unfactoredTerms;
          factoredTerms.reserve(terms.size());
          unfactoredTerms.reserve(terms.size());
          for (size_t i = 0; i < terms.size(); ++i) {
            if (factorPerTerm.at(i) == factor) {
              factoredTerms.push_back(terms[i].factorOut(factor, leading));
            } else {
              unfactoredTerms.push_back(terms[i]);
            }
          }

          auto sumFactored = ExprFactory::sum(std::move(factoredTerms));
          auto factorTimesFactored =
              leading ? ExprFactory::product({factor, std::move(sumFactored)})
                      : ExprFactory::product({std::move(sumFactored), factor});
          auto factoredExpr =
              (unfactoredTerms.empty()
                   ? std::move(factorTimesFactored)
                   : ExprFactory::sum(
                         {ExprFactory::sum(std::move(unfactoredTerms)),
                          factorTimesFactored}))
                  .simplify(false);
          match(
              factoredExpr,
              [&](auto& x)
                requires is_nary_v<decltype(x)>
              {
                associativeTransformation<std::decay_t<decltype(x)>>(x.terms);
              },
              [&](const auto&) {});
          if (factoredExpr.complexity() < simplified.complexity()) {
            return factoredExpr;
          }
          sorted.pop_back();
        }
      }
    }

    if (terms.size() == 1) {
      return terms.front();
    }

    return simplified;
  }

  Expr operator()(const Product& x) const {
    // Recursive simplification
    auto terms = transform(
        x.terms, [this](const auto& t) { return t.simplify(distribute); });

    // Associative transformation (x(yz) = xyz)
    associativeTransformation_<Product>(terms);

    // Identity transformations (x * 0 = 0, x * 1 = x)
    if (std::ranges::any_of(terms, [](const auto& t) { return t == zero; })) {
      return zero;
    }
    if (std::ranges::all_of(terms, [](const auto& t) { return t == unity; })) {
      return unity;
    }
    std::erase_if(terms, [](const auto& t) { return t == unity; });

    // x * (-1) * y = -x * y
    const auto it = std::ranges::find_if(terms, is<Negate>);
    if (it != terms.end()) {
      auto newTerms = terms;  // Seems to be necessary to make a new copy
      // for this to work in WebAssembly
      newTerms[std::distance(terms.begin(), it)] =
          *std::get<Negate>(it->getImpl()).child;
      return ExprFactory::negate(ExprFactory::product(std::move(terms)));
    }

    // Canceling transformation (power transformation for the special case
    // of inverse) (x * inv(x) = 1)
    eraseCanceling_<Invert>(terms, unity);

    // Commutative transformation (2x3y = 2*3*xy)
    std::ranges::partition(terms, is<Number>);

    // Numerical transformation (2 * x * 3 = 6 * x)
    if (std::ranges::any_of(terms, is<Number>)) {
      const auto value = std::accumulate(
          terms.cbegin(), terms.cend(), 1.0, [](const double s, const auto& t) {
            return s *
                   (is<Number>(t) ? std::get<Number>(t.getImpl()).value : 1.0);
          });
      std::erase_if(terms, is<Number>);
      terms.push_back(ExprFactory::number(value));
    }

    if (terms.size() == 1) {
      return terms.front();
    }

    auto simplified = ExprFactory::product(terms);

    // Check whether distributing factor leads to a shorter expression (x(y + z
    // + w) = xy + xz + xw)
    if (distribute && terms.size() > 1) {
      for (size_t i = 0; i < terms.size(); ++i) {
        if (match(
                terms[i],
                [&](const Sum& x) {
                  const auto init =
                      std::vector<Expr>(terms.begin(), terms.begin() + i);
                  const auto rest =
                      std::vector<Expr>(terms.begin() + i + 1, terms.end());
                  auto sumTerms =
                      transform(x.terms, [&init, &rest](const auto& t) {
                        auto factors = init;
                        factors.push_back(t);
                        factors.insert(factors.end(), rest.begin(), rest.end());
                        auto term = ExprFactory::product(factors);
                        return term;
                      });
                  if (auto distributedExpr =
                          ExprFactory::sum(std::move(sumTerms)).simplify(false);
                      distributedExpr.complexity() < simplified.complexity()) {
                    simplified = std::move(distributedExpr);
                    return true;
                  }
                  return false;
                },
                [&](const auto&) { return false; })) {
          return simplified;
        }
      }
    }

    return simplified;
  }

 private:
  template <typename T>
  void eraseCanceling_(std::vector<Expr>& terms,
                       const Expr& replacement) const {
    for (size_t i = 0; i < terms.size(); ++i) {
      const auto& t1 = terms.at(i);
      for (size_t j = i + 1; j < terms.size(); ++j) {
        const auto& t2 = terms.at(j);
        if ((is<T>(t1) && *std::get<T>(t1.getImpl()).child == t2) ||
            (is<T>(t2) && *std::get<T>(t2.getImpl()).child == t1)) {
          terms.erase(terms.begin() + j);
          terms.erase(terms.begin() + i);
          terms.insert(terms.begin() + i, replacement);
          break;
        }
      }
    }
  }

  // Associative transformation lambda
  template <typename T>
    requires is_nary_v<T>
  void associativeTransformation_(std::vector<Expr>& terms) const {
    std::vector<Expr> newTerms;
    newTerms.reserve(terms.size());
    for (auto& t : terms) {
      match(
          t,
          [&](const T& x) {
            newTerms.insert(newTerms.end(), x.terms.begin(), x.terms.end());
          },
          [&](const auto&) { newTerms.emplace_back(std::move(t)); });
    }
    terms = std::move(newTerms);
  }
};
}  // namespace Expression