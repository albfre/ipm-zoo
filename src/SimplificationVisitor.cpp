#include "SimplificationVisitor.h"

#include <algorithm>
#include <map>
#include <numeric>
#include <ranges>

#include "Assert.h"

namespace Expression {
SimplificationVisitor::SimplificationVisitor(const bool distribute)
    : distribute_(distribute) {}

Expr SimplificationVisitor::operator()(const DiagonalMatrix& x) {
  auto child = x.child->simplifyOnce(distribute_);
  if (child == zero || child == unity) {
    return child;
  }
  return ExprFactory::diagonalMatrix(std::move(child));
}

Expr SimplificationVisitor::operator()(const Transpose& x) const {
  auto child = x.child->simplifyOnce(distribute_);
  if (child == zero || child == unity) {
    return child;  // 0^T = 0, 1^T = 1
  }
  return match(child).with(
      [](const Transpose& x) { return *x.child; },  // (x^T)^T = x
      [&](const NamedScalar& x) { return child; },
      [&](const SymmetricMatrix& x) { return child; },
      [&](const DiagonalMatrix& x) { return child; },
      [&](const Invert& x) {
        ASSERT(is<DiagonalMatrix>(*x.child));
        return child;
      },
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
      [&](const auto& x) { return ExprFactory::transpose(std::move(child)); });
}

Expr SimplificationVisitor::operator()(const Negate& x) {
  auto child = x.child->simplifyOnce(distribute_);
  if (child == zero) {  // -0 = 0
    return child;
  }
  return match(child).with(
      [](const Negate& x) { return *x.child; },  // negate(negate(x)) = x
      [&](const Product& x) {
        if (const auto it = std::ranges::find_if(x.terms, is<Negate>);
            it != x.terms.end()) {
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
          auto terms = std::vector<Expr>();
          terms.reserve(x.terms.size());
          for (const auto& t : x.terms) {
            terms.emplace_back(match(t).with(
                [](const Negate& y) { return *y.child; },
                [&](const auto&) { return ExprFactory::negate(t); }));
          }
          return ExprFactory::sum(std::move(terms));
        }
        return ExprFactory::negate(std::move(child));
      },
      [&](const auto& x) { return ExprFactory::negate(std::move(child)); });
}

Expr SimplificationVisitor::operator()(const Invert& x) const {
  auto child = x.child->simplifyOnce(distribute_);
  if (child == unity) {
    return child;
  }
  return match(child).with(
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

Expr SimplificationVisitor::operator()(const Log& x) const {
  return ExprFactory::log(x.child->simplifyOnce(distribute_));
}

Expr SimplificationVisitor::operator()(const Sum& x) const {
  // Recursive simplification
  auto terms = transform(
      x.terms, [this](const auto& t) { return t.simplifyOnce(distribute_); });

  // Associative transformation ((x + y) + z = x + y + z)
  std::vector<Expr> newTerms;
  newTerms.reserve(terms.size());
  for (auto& t : terms) {
    match(t).with(
        [&](const Sum& x) {
          newTerms.insert(newTerms.end(), x.terms.begin(), x.terms.end());
        },
        [&](auto& x) {
          if constexpr (std::is_same_v<Negate, std::decay_t<decltype(x)>>) {
            if (match(*x.child).with([](const Sum& sum) { return true; },
                                     [](const auto&) { return false; })) {
              auto& sum = std::get<Sum>(x.child->getImpl());
              newTerms.reserve(newTerms.size() + sum.terms.size());
              for (auto& ct : sum.terms) {
                newTerms.emplace_back(ExprFactory::negate(ct));
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
    if (terms.at(i) != zero) {
      auto term = terms.at(i);
      const auto negTerm = ExprFactory::negate(term);
      const auto isNumberTimesTerm = [&term](const auto& t) -> bool {
        return match(t).with(
            [&term](const Product& x) {
              return x.terms.size() == 2 && is<Number>(x.terms.front()) &&
                     x.terms.back() == term;
            },
            [](const auto&) { return false; });
      };
      const auto isTerm = [&](const auto& t) {
        return t == term || t == negTerm || isNumberTimesTerm(t);
      };
      if (std::ranges::count_if(terms, isTerm) > 1) {
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
        terms.push_back(ExprFactory::product(
            {ExprFactory::number(value), std::move(term)}));
      }
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
  if (std::ranges::count_if(terms, is<Number>) > 1) {
    const auto value = std::accumulate(
        terms.cbegin(), terms.cend(), 0.0, [](const double s, const auto& t) {
          return s + match(t).with([](const Number& x) { return x.value; },
                                   [](const auto&) { return 0.0; });
        });
    std::erase_if(terms, is<Number>);
    terms.push_back(ExprFactory::number(value));
  }

  // Commutative transformation (z + y + x = x + y + z)
  std::ranges::sort(terms);

  // -x - y = -(x + y)
  if (std::ranges::all_of(terms, is<Negate>)) {
    auto newTerms = transform(terms, [](const auto& t) {
      return *std::get<Negate>(t.getImpl()).child;
    });
    return ExprFactory::negate(ExprFactory::sum(std::move(newTerms)));
  }

  if (terms.size() == 1) {
    return terms.front();
  }

  auto simplified = ExprFactory::sum(terms);

  // Check whether extracting common factors leads to a shorter expression
  // (xy + xz + xw = x(y + z + w)
  if (distribute_) {
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
        auto [factor, numOccurrences] = sorted.back();
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
            leading ? ExprFactory::product(
                          {std::move(factor), std::move(sumFactored)})
                    : ExprFactory::product(
                          {std::move(sumFactored), std::move(factor)});
        auto factoredExpr =
            (unfactoredTerms.empty()
                 ? std::move(factorTimesFactored)
                 : ExprFactory::sum(
                       {ExprFactory::sum(std::move(unfactoredTerms)),
                        std::move(factorTimesFactored)}))
                .simplify(false);
        if (factoredExpr.complexity() < simplified.complexity()) {
          return factoredExpr;
        }
        sorted.pop_back();
      }
    }
  }

  return simplified;
}

Expr SimplificationVisitor::operator()(const Product& x) const {
  // Recursive simplification
  auto terms = transform(
      x.terms, [this](const auto& t) { return t.simplifyOnce(distribute_); });

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
  if (const auto it = std::ranges::find_if(terms, is<Negate>);
      it != terms.end()) {
    auto newTerms = terms;  // Seems to be necessary to make a new copy
    // for this to work in WebAssembly
    newTerms[std::distance(terms.begin(), it)] =
        *std::get<Negate>(it->getImpl()).child;
    return ExprFactory::negate(ExprFactory::product(std::move(newTerms)));
  }

  // Canceling transformation (power transformation for the special case
  // of inverse) (x * inv(x) = 1)
  eraseCanceling_<Invert>(terms, unity);

  // Commutative transformation (2x3y = 2*3*xy)
  std::ranges::partition(terms, is<NamedScalar>);
  std::ranges::partition(terms, is<Number>);

  // Numerical transformation (2 * x * 3 = 6 * x)
  if (std::ranges::count_if(terms, is<Number>) > 1) {
    const auto value = std::accumulate(
        terms.cbegin(), terms.cend(), 1.0, [](const double s, const auto& t) {
          return s *
                 (is<Number>(t) ? std::get<Number>(t.getImpl()).value : 1.0);
        });
    std::erase_if(terms, is<Number>);
    terms.insert(terms.begin(), ExprFactory::number(value));
  }

  if (terms.size() == 1) {
    return terms.front();
  }

  auto simplified = ExprFactory::product(terms);

  // Check whether distributing factor leads to a shorter expression (x(y + z
  // + w) = xy + xz + xw)
  if (distribute_ && terms.size() > 1) {
    for (size_t i = 0; i < terms.size(); ++i) {
      if (match(terms[i]).with(
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
                    distributedExpr.complexity() <= simplified.complexity()) {
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

template <typename T>
void SimplificationVisitor::eraseCanceling_(std::vector<Expr>& terms,
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

template <typename T>
  requires NaryType<T> void SimplificationVisitor::associativeTransformation_(
    std::vector<Expr>& terms) const {
  std::vector<Expr> newTerms;
  newTerms.reserve(terms.size());
  for (auto& t : terms) {
    match(t).with(
        [&](const T& x) {
          newTerms.insert(newTerms.end(), x.terms.begin(), x.terms.end());
        },
        [&](const auto&) { newTerms.emplace_back(std::move(t)); });
  }
  terms = std::move(newTerms);
}

// Explicit template instantiations for the templates used
template void SimplificationVisitor::eraseCanceling_<Negate>(
    std::vector<Expr>& terms, const Expr& replacement) const;
template void SimplificationVisitor::eraseCanceling_<Invert>(
    std::vector<Expr>& terms, const Expr& replacement) const;
template void SimplificationVisitor::associativeTransformation_<Sum>(
    std::vector<Expr>& terms) const;
template void SimplificationVisitor::associativeTransformation_<Product>(
    std::vector<Expr>& terms) const;

}  // namespace Expression