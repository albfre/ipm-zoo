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

  // Associative transformation lambda
  /*
  const auto associativeTransformation = [](const ExprType type, auto& terms) {
    assert((std::set{ExprType::Sum, ExprType::Product}.contains(type)));
    std::vector<Expr> newTerms;
    newTerms.reserve(terms.size());
    for (const auto& t : terms) {
      if (t.type_ == type) {
        newTerms.insert(newTerms.end(), t.terms_.begin(), t.terms_.end());
      } else {
        newTerms.push_back(t);
      }
    }
    terms = std::move(newTerms);
  };
  */

  Expr operator()(const Number& x) const { return Expr(x); }
  Expr operator()(const NamedConstant& x) const { return Expr(x); }
  Expr operator()(const Matrix& x) const { return Expr(x); }
  Expr operator()(const SymmetricMatrix& x) const { return Expr(x); }
  Expr operator()(const Variable& x) const { return Expr(x); }

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
        [&](const auto& x) { return ExprFactory::transpose(child); });
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
          return ExprFactory::negate(child);
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

  Expr operator()(const auto& x) const { return ExprFactory::variable("p"); }

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

    // Numerical transformation (1 + x + 2 = 3 + x)
    if (std::ranges::any_of(terms,
                            [](const auto& t) { return is<Number>(t); })) {
      const auto value = std::accumulate(
          terms.cbegin(), terms.cend(), 0.0, [](const double s, const auto& t) {
            return s + match(
                           t, [](const Number& x) { return x.value; },
                           [](const auto&) { return 0.0; });
          });
      std::erase_if(terms, [](const auto& t) { return is<Number>(t); });
      terms.push_back(ExprFactory::number(value));
    }

    // Commutative transformation (z + y + x = x + y + z)
    std::ranges::sort(terms);

    auto simplified = ExprFactory::sum(terms);

    // Check whether extracting common factors leads to a shorter expression
    // (xy + xz + xw = x(y + z + w)
    /*
    if (distribute) {
      for (const auto leading : std::vector{true, false}) {
        std::map<Expr, size_t> numOccurancesOfFactors;
        for (auto& t : terms) {
          ++numOccurancesOfFactors[t.getLeadingOrEndingFactor_(leading)];
        }

        if (const auto it =
                std::ranges::max_element(numOccurancesOfFactors,
                                         [](const auto& p1, const auto& p2) {
                                           return p1.second < p2.second;
                                         });
            it->second >= 2) {
          std::vector<Expr> factoredTerms;
          std::vector<Expr> unfactoredTerms;
          const auto& factor = it->first;
          for (auto& t : terms) {
            if (t.getLeadingOrEndingFactor_(leading) == factor) {
              factoredTerms.push_back(t.factorOut(factor, leading));
            } else {
              unfactoredTerms.push_back(t);
            }
          }

          const auto sumFactored = ExprFactory::sum(std::move(factoredTerms));
          const auto factorTimesFactored =
              leading ? ExprFactory::product({factor, sumFactored})
                      : ExprFactory::product({sumFactored, factor});
          auto factoredExpr =
              (unfactoredTerms.empty()
                   ? factorTimesFactored
                   : ExprFactory::sum(
                         {ExprFactory::sum(std::move(unfactoredTerms)),
                          factorTimesFactored}))
                  .simplify_();
          associativeTransformation(factoredExpr.type_, factoredExpr.terms_);
          if (factoredExpr.complexity_() < simplified.complexity_()) {
            return factoredExpr;
          }
        }
      }
    }
    */

    if (terms.size() == 1) {
      return terms.front();
    }

    return simplified;
  }
  /*
    switch (type_) {
      case ExprType::Product: {
        // Recursive simplification
        auto terms =
            transform(terms_, [](const auto& t) { return t.simplify_(); });

        // Associative transformation (x(yz) = xyz)
        associativeTransformation(ExprType::Product, terms);

        // Identity transformations (x * 0 = 0, x * 1 = x)
        if (std::ranges::any_of(terms, [](const auto& t) { return t == zero;
    })) { return zero;
        }
        if (std::ranges::all_of(terms,
                                [](const auto& t) { return t == unity; })) {
          return unity;
        }
        std::erase_if(terms, [](const auto& t) { return t == unity; });

        // x * (-1) * y = -x * y
        if (std::ranges::any_of(terms, [](const auto& t) {
              return t.type_ == ExprType::Negate;
            })) {
          auto newTerms = terms;  // Seems to be necessary to make a new copy
          // for this to work in WebAssembly
          for (size_t i = 0; i < terms.size(); ++i) {
            if (terms[i].type_ == ExprType::Negate) {
              newTerms[i] = terms[i].getSingleChild_();
              return ExprFactory::negate(
                  ExprFactory::product(std::move(newTerms)));
            }
          }
        }

        // Canceling transformation (power transformation for the special case
        // of inverse) (x * inv(x) = 1)
        eraseCanceling(ExprType::Invert, terms, unity);

        // Commutative transformation (2x3y = 6xy)
        std::partition(terms.begin(), terms.end(),
                       [](const auto& t) { return t.type_ == ExprType::Number;
    });

        // Numerical transformation (2 * x * 3 = 6 * x)
        if (std::ranges::count_if(terms, [](const auto& t) {
              return t.type_ == ExprType::Number;
            }) > 1) {
          const auto value = std::accumulate(
              terms.cbegin(), terms.cend(), 1.0,
              [](const double s, const auto& t) {
                return s * (t.type_ == ExprType::Number ? t.value_ : 1.0);
              });
          std::erase_if(
              terms, [](const auto& t) { return t.type_ == ExprType::Number;
    }); terms.push_back(ExprFactory::number(value));
        }

        const auto simplified = ExprFactory::product(terms);

        // Check whether distributing factor leads to a shorter expression
    (x(y
        // + z + w) = xy + xz + xw)
        if (distribute && terms.size() > 1) {
          for (size_t i = 0; i < terms.size(); ++i) {
            if (terms[i].type_ == ExprType::Sum) {
              auto factors = terms;
              factors.erase(factors.begin() + i);
              auto sumTerms =
                  transform(terms[i].terms_, [&factors](const auto& t) {
                    factors.push_back(t);
                    auto term = ExprFactory::product(factors);
                    factors.pop_back();
                    return term;
                  });
              auto distributedExpr =
                  ExprFactory::sum(std::move(sumTerms)).simplify(false);
              if (distributedExpr.complexity_() < simplified.complexity_()) {
                return distributedExpr;
              }
            }
          }
        }

        if (terms.size() == 1) {
          return terms.front();
        }

        return simplified;
      }
      default:
        assert(false);
        throw std::runtime_error("Unhandled enum");
    }
    */

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
};
}  // namespace Expression