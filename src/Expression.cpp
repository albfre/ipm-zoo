#include "Expression.h"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <map>
#include <numeric>
#include <ranges>
#include <set>
#include <sstream>

namespace Expression {
const auto unity = ExprFactory::number(1.0);
const auto zero = ExprFactory::number(0.0);

std::strong_ordering operator<=>(const Expr& left, const Expr& right) {
  // Order lexicographically, first by type, then by the expression string
  if (&left == &right) {
    return std::strong_ordering::equal;
  }
  if (left.getType() != right.getType()) {
    return left.getType() <=> right.getType();
  }
  return left.toExpressionString() <=> right.toExpressionString();
}

bool operator==(const Expr& left, const Expr& right) {
  return (left <=> right) == std::strong_ordering::equal;
}

template <typename TLambda>
std::vector<Expr> transform(const std::vector<Expr>& terms,
                            const TLambda& lambda) {
  std::vector<Expr> transformedTerms;
  transformedTerms.reserve(terms.size());
  std::ranges::transform(terms, std::back_inserter(transformedTerms),
                         [&lambda](const auto& t) { return lambda(t); });
  return transformedTerms;
}

Expr::Expr(ExprType type, std::vector<Expr> terms)
    : type_(type), terms_(std::move(terms)) {
  assert((std::set{ExprType::Transpose, ExprType::Negate, ExprType::Invert,
                   ExprType::Log, ExprType::Product, ExprType::Sum,
                   ExprType::DiagonalMatrix}
              .contains(type_)));
  if (std::set{ExprType::Transpose, ExprType::Negate, ExprType::Invert,
               ExprType::Log, ExprType::DiagonalMatrix}
          .contains(type_)) {
    assert(terms_.size() == 1);
  }
  assert(!terms_.empty());
}

Expr::Expr(ExprType type, const std::string& name) : type_(type), name_(name) {
  assert(
      (std::set{ExprType::NamedScalar, ExprType::NamedVector,
                ExprType::Variable, ExprType::SymmetricMatrix, ExprType::Matrix}
           .contains(type_)));
}

Expr::Expr(const double value) : type_(ExprType::Number), value_(value) {}

Expr Expr::differentiate(const Expr& var) const {
  if (!containsSubexpression(var)) {
    return zero;
  }
  switch (type_) {
    case ExprType::Number:
    case ExprType::NamedScalar:
    case ExprType::NamedVector:
    case ExprType::SymmetricMatrix:
    case ExprType::Matrix:
      return zero;
    case ExprType::Variable:
      return *this == var ? unity : zero;
    case ExprType::DiagonalMatrix:
      return ExprFactory::diagonalMatrix(getSingleChild_().differentiate(var));
    case ExprType::Transpose:
      return ExprFactory::transpose(getSingleChild_().differentiate(var));
    case ExprType::Negate:
      return ExprFactory::negate(getSingleChild_().differentiate(var));
    case ExprType::Invert:
      assert(false);  // Not implemented
    case ExprType::Log: {
      const auto& child = getSingleChild_();
      return ExprFactory::product(
          {ExprFactory::invert(ExprFactory::diagonalMatrix(child)),
           child.differentiate(var)});
    }
    case ExprType::Sum: {
      auto terms = transform(
          terms_, [&var](const auto& t) { return t.differentiate(var); });
      return ExprFactory::sum(std::move(terms));
    }
    case ExprType::Product: {
      std::vector<Expr> sumTerms;
      sumTerms.reserve(terms_.size());
      for (size_t i = 0; i < terms_.size(); ++i) {
        {
          auto terms = terms_;
          terms[i] = terms_[i].differentiate(var);
          if (i + 2 == terms_.size() &&
              terms_[i].type_ == ExprType::DiagonalMatrix &&
              terms_[i + 1].type_ == ExprType::Variable) {
            terms[i + 1] = ExprFactory::diagonalMatrix(terms[i + 1]);
          }

          sumTerms.push_back(
              ExprFactory::product(std::move(terms)));  // d/dx(f(x)) g(x)
        }
        if (terms_[i].type_ == ExprType::Transpose &&
            terms_[i].getSingleChild_().type_ != ExprType::Matrix &&
            i + 1 < terms_.size()) {
          // d/dx(f(x)^T g(x)) = d/dx(f(x)^T) g(x) + d/dx(g(x))^T f(x)
          auto terms = std::vector(terms_.begin(), terms_.begin() + i);
          auto rest = std::vector(terms_.begin() + i + 1, terms_.end());
          auto restTerm = rest.size() == 1
                              ? rest.front()
                              : ExprFactory::product(std::move(rest));
          terms.push_back(ExprFactory::transpose(restTerm).differentiate(
              var));                                     // d/dx(g(x))^T
          terms.push_back(terms_[i].getSingleChild_());  // f(x)
          sumTerms.push_back(
              ExprFactory::product(std::move(terms)));  // d/dx(g(x))^T f(x)
          break;
        }
      }
      return ExprFactory::sum(std::move(sumTerms));
    }
    default:
      assert(false);
      throw std::runtime_error("Unhandled enum");
  }
}

Expr Expr::simplify(const bool distribute) const {
  auto expr = *this;
  auto changed = true;
  while (changed) {
    auto simplified = expr.simplify_(distribute);
    changed = expr != simplified;
    std::swap(expr, simplified);
  }
  return expr;
}

std::string Expr::toString(const bool condensed) const {
  switch (type_) {
    case ExprType::Number: {
      std::stringstream ss;
      ss << value_;
      return ss.str();
    }
    case ExprType::NamedScalar:
    case ExprType::NamedVector:
    case ExprType::SymmetricMatrix:
    case ExprType::Matrix:
    case ExprType::Variable:
      return name_;
    case ExprType::DiagonalMatrix: {
      const auto& child = getSingleChild_();
      if (child.type_ == ExprType::Variable) {
        auto name = child.name_;
        for (size_t i = 0; i < name.size(); ++i) {
          if (std::isalpha(name[i])) {
            name[i] = std::toupper(name[i]);
            return name;
          }
        }
      }
      return "\\diag(" + child.toString(condensed) + ")";
    }
    case ExprType::Transpose: {
      const auto& child = getSingleChild_();
      if (condensed &&
          std::set{ExprType::Sum, ExprType::Product, ExprType::Invert}.contains(
              child.type_)) {
        return "(" + child.toString(condensed) + ")^T";
      }
      return child.toString(condensed) + "^T";
    }
    case ExprType::Negate: {
      const auto& child = getSingleChild_();
      if (condensed && child.type_ == ExprType::Sum) {
        return "-(" + child.toString(condensed) + ")";
      }
      return "-" + child.toString(condensed);
    }
    case ExprType::Invert: {
      const auto& child = getSingleChild_();
      if (condensed &&
          std::set{ExprType::Sum, ExprType::Product, ExprType::Transpose}
              .contains(child.type_)) {
        return "(" + child.toString(condensed) + ")^{-1}";
      }
      return child.toString(condensed) + "^{-1}";
    }
    case ExprType::Log:
      return "\\log(" + getSingleChild_().toString(condensed) + ")";
    case ExprType::Sum: {
      std::stringstream ss;
      ss << (condensed ? "" : "(");
      ss << terms_.front().toString(condensed);
      for (const auto& t : terms_ | std::views::drop(1)) {
        const auto negate = t.type_ == ExprType::Negate;
        if (negate) {
          ss << " - " << t.getSingleChild_().toString(condensed);
        } else {
          ss << " + " << t.toString(condensed);
        }
      }
      ss << (condensed ? "" : ")");
      return ss.str();
    }
    case ExprType::Product: {
      std::stringstream ss;
      ss << (condensed ? "" : "(");
      const auto& front = terms_.front();
      if (condensed && front.type_ == ExprType::Sum) {
        ss << "(" << front.toString(condensed) << ")";
      } else {
        ss << front.toString(condensed);
      }
      const auto symbol = condensed ? " " : " * ";
      for (const auto& t : terms_ | std::views::drop(1)) {
        const auto negate = t.type_ == ExprType::Negate;
        const auto sum = t.type_ == ExprType::Sum;
        if (negate || (condensed && sum)) {
          ss << symbol << "(" << t.toString(condensed) << ")";
        } else {
          ss << symbol << t.toString(condensed);
        }
      }
      ss << (condensed ? "" : ")");
      return ss.str();
    }
    default:
      assert(false);
      throw std::runtime_error("Unhandled enum");
  }
}

std::string Expr::toExpressionString() const {
  switch (type_) {
    case ExprType::Number: {
      std::stringstream ss;
      ss << "number(" << value_ << ")";
      return ss.str();
    }
    case ExprType::NamedScalar:
      return "namedScalar(" + name_ + ")";
    case ExprType::NamedVector:
      return "namedVector(" + name_ + ")";
    case ExprType::Matrix:
      return "matrix(" + name_ + ")";
    case ExprType::SymmetricMatrix:
      return "symmetricMatrix(" + name_ + ")";
    case ExprType::Variable:
      return "variable(" + name_ + ")";
    case ExprType::DiagonalMatrix:
      return "diagonalMatrix(" + getSingleChild_().toExpressionString() + ")";
    case ExprType::Transpose:
      return "transpose(" + getSingleChild_().toExpressionString() + ")";
    case ExprType::Negate:
      return "negate(" + getSingleChild_().toExpressionString() + ")";
    case ExprType::Invert:
      return "invert(" + getSingleChild_().toExpressionString() + ")";
    case ExprType::Log:
      return "log(" + getSingleChild_().toExpressionString() + ")";
    case ExprType::Sum: {
      std::stringstream ss;
      ss << "sum(";
      ss << terms_.front().toExpressionString();
      for (const auto& t : terms_ | std::views::drop(1)) {
        ss << ", " << t.toExpressionString();
      }
      ss << ")";
      return ss.str();
    }
    case ExprType::Product: {
      std::stringstream ss;
      ss << "product(";
      ss << terms_.front().toExpressionString();
      for (const auto& t : terms_ | std::views::drop(1)) {
        ss << ", " << t.toExpressionString();
      }
      ss << ")";
      return ss.str();
    }
    default:
      assert(false);
      throw std::runtime_error("Unhandled enum");
  }
}

ExprType Expr::getType() const { return type_; }

const std::string& Expr::getName() const {
  assert(
      (std::set{ExprType::NamedScalar, ExprType::NamedVector,
                ExprType::Variable, ExprType::Matrix, ExprType::SymmetricMatrix}
           .contains(type_)));
  return name_;
}

std::set<Expr> Expr::getVariables() const {
  std::set<Expr> variables;
  for (auto& t : terms_) {
    if (t.type_ == ExprType::Variable) {
      variables.insert(t);
    } else {
      const auto v = t.getVariables();
      variables.insert(v.begin(), v.end());
    }
  }
  return variables;
}

bool Expr::containsSubexpression(const Expr& expr) const {
  if (*this == expr) {
    return true;
  }
  return std::ranges::any_of(
      terms_, [&expr](const auto& t) { return t.containsSubexpression(expr); });
}

Expr Expr::replaceSubexpression(const Expr& expr,
                                const Expr& replacement) const {
  if (*this == expr) {
    return replacement;
  }
  auto terms = std::vector<Expr>();
  terms.reserve(terms_.size());
  for (const auto& t : terms_) {
    terms.push_back(t.replaceSubexpression(expr, replacement));
  }
  auto newExpr = *this;
  newExpr.terms_ = std::move(terms);
  return newExpr;
}

Expr Expr::simplify_(const bool distribute) const {
  // Canceling transformation lambda
  const auto eraseCanceling = [](const ExprType type, auto& terms,
                                 const auto& replacement) {
    for (size_t i = 0; i < terms.size(); ++i) {
      const auto& t1 = terms.at(i);
      for (size_t j = i + 1; j < terms.size(); ++j) {
        const auto& t2 = terms.at(j);
        if ((t1.type_ == type && t1.getSingleChild_() == t2) ||
            (t2.type_ == type && t2.getSingleChild_() == t1)) {
          terms.erase(terms.begin() + j);
          terms.erase(terms.begin() + i);
          terms.insert(terms.begin() + i, replacement);
          break;
        }
      }
    }
  };

  // Associative transformation lambda
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
  switch (type_) {
    case ExprType::Number:
    case ExprType::NamedScalar:
    case ExprType::NamedVector:
    case ExprType::Matrix:
    case ExprType::SymmetricMatrix:
    case ExprType::Variable:
      return *this;
    case ExprType::DiagonalMatrix: {
      auto child = getSingleChild_().simplify_(distribute);
      if (child == zero || child == unity) {
        return child;
      }
      return ExprFactory::diagonalMatrix(std::move(child));
    }
    case ExprType::Transpose: {
      auto child = getSingleChild_().simplify_(distribute);
      if (child.type_ == ExprType::Transpose) {
        // transpose(transpose(x)) = x;
        return child.getSingleChild_();
      } else if (child == zero || child == unity ||
                 child.type_ == ExprType::NamedScalar) {
        // 0^T = 0, 1^T = 1
        return child;
      } else if (child.type_ == ExprType::SymmetricMatrix ||
                 child.type_ == ExprType::DiagonalMatrix) {
        return child;
      } else if (child.type_ == ExprType::Invert) {
        if (child.getSingleChild_().type_ == ExprType::DiagonalMatrix) {
          return child;
        }
        assert(false);  // Not implemented
      } else if (child.type_ == ExprType::Negate) {
        // (-x)^T  = -x^T
        return ExprFactory::negate(
            ExprFactory::transpose(child.getSingleChild_()));
      } else if (child.type_ == ExprType::Product) {
        // (xyz)^T = z^T y^T x^T
        auto terms = transform(child.terms_, [](const auto& t) {
          return ExprFactory::transpose(t);
        });
        std::ranges::reverse(terms);
        return ExprFactory::product(std::move(terms));
      }
      return ExprFactory::transpose(std::move(child));
    }
    case ExprType::Negate: {
      auto child = getSingleChild_().simplify_(distribute);
      if (child.type_ == ExprType::Negate) {
        // negate(negate(x)) = x
        return child.getSingleChild_();
      } else if (child == zero) {
        // -0 = 0
        return zero;
      } else if (child.type_ == ExprType::Product) {
        for (size_t i = 0; i < child.terms_.size(); ++i) {
          if (child.terms_[i].type_ == ExprType::Negate) {
            auto terms = child.terms_;
            terms[i] = child.terms_[i].getSingleChild_();
            return ExprFactory::product(std::move(terms));
          }
        }
      } else if (child.type_ == ExprType::Sum) {
        if (static_cast<size_t>(
                std::ranges::count_if(child.terms_, [](const auto& t) {
                  return t.type_ == ExprType::Negate;
                })) > child.terms_.size() / 2) {
          auto terms = std::vector<Expr>();
          terms.reserve(child.terms_.size());
          for (const auto& t : child.terms_) {
            terms.emplace_back(t.type_ == ExprType::Negate
                                   ? t.getSingleChild_()
                                   : ExprFactory::negate(t));
          }
          return ExprFactory::sum(std::move(terms));
        }
      }
      return ExprFactory::negate(std::move(child));
    }
    case ExprType::Invert: {
      auto child = getSingleChild_().simplify_(distribute);
      if (child == unity) {
        return child;
      } else if (child.type_ == ExprType::Invert) {
        // invert(invert(x)) = x
        return child.getSingleChild_();
      } else if (child.type_ == ExprType::Negate) {
        // invert(negate(x)) = negate(invert(x))
        return ExprFactory::negate(
            ExprFactory::invert(child.getSingleChild_()));
      } else if (child.type_ == ExprType::Product) {
        // invert(x * y * z) = invert(z) * invert(y) * invert(x)
        auto terms = transform(
            child.terms_, [](const auto& t) { return ExprFactory::invert(t); });
        std::ranges::reverse(terms);
        return ExprFactory::product(std::move(terms));
      }
      return ExprFactory::invert(std::move(child));
    }
    case ExprType::Log: {
      return ExprFactory::log(getSingleChild_().simplify_(distribute));
    }
    case ExprType::Sum: {
      // Recursive simplification
      auto terms = transform(terms_, [distribute](const auto& t) {
        return t.simplify_(distribute);
      });

      // Associative transformation ((x + y) + z = x + y + z)
      std::vector<Expr> newTerms;
      newTerms.reserve(terms.size());
      for (auto& t : terms) {
        if (t.type_ == ExprType::Sum) {
          newTerms.insert(newTerms.end(), t.terms_.begin(), t.terms_.end());
        } else if (t.type_ == ExprType::Negate &&
                   t.getSingleChild_().type_ == ExprType::Sum) {
          auto& childTerms = t.getSingleChild_().terms_;
          newTerms.reserve(newTerms.size() + childTerms.size());
          for (auto& ct : childTerms) {
            newTerms.emplace_back(ExprFactory::negate(std::move(ct)));
          }
        } else {
          newTerms.emplace_back(std::move(t));
        }
      }
      terms = std::move(newTerms);

      // Distributive transformation (x + y + 1.3x = 2.3x + y)
      for (size_t i = 0; i < terms.size(); ++i) {
        if (terms.at(i) != zero) {
          const auto term = terms.at(i);
          const auto negTerm = ExprFactory::negate(term);
          const auto isNumberTimesTerm = [&term](const auto& t) {
            return t.type_ == ExprType::Product && t.terms_.size() == 2 &&
                   t.terms_.at(0).type_ == ExprType::Number &&
                   t.terms_.at(1) == term;
          };
          const auto isTerm = [&term, &negTerm,
                               &isNumberTimesTerm](const auto& t) {
            return t == term || t == negTerm || isNumberTimesTerm(t);
          };
          if (std::ranges::count_if(terms, isTerm) > 1) {
            const auto value = std::accumulate(
                terms.cbegin(), terms.cend(), 0.0,
                [&term, &negTerm, &isNumberTimesTerm](const double s,
                                                      const auto& t) {
                  return s + (t == term              ? 1.0
                              : t == negTerm         ? -1.0
                              : isNumberTimesTerm(t) ? t.terms_.at(0).value_
                                                     : 0.0);
                });
            std::erase_if(terms, isTerm);
            terms.push_back(
                ExprFactory::product({ExprFactory::number(value), term}));
          }
        }
      }

      // Identity transformations (x + 0 = x)
      std::erase_if(terms, [&](const auto& t) { return t == zero; });

      if (terms.empty()) {
        return zero;
      }

      // Canceling transformation (x - x = 0)
      eraseCanceling(ExprType::Negate, terms, zero);

      // Numerical transformation (1 + x + 2 = 3 + x)
      if (std::ranges::count_if(terms, [](const auto& t) {
            return t.type_ == ExprType::Number;
          }) > 1) {
        const auto value = std::accumulate(
            terms.cbegin(), terms.cend(), 0.0,
            [](const double s, const auto& t) {
              return s + (t.type_ == ExprType::Number ? t.value_ : 0.0);
            });
        std::erase_if(
            terms, [](const auto& t) { return t.type_ == ExprType::Number; });
        terms.push_back(ExprFactory::number(value));
      }

      // Commutative transformation (z + y + x = x + y + z)
      std::ranges::sort(terms);

      // -x - y = -(x + y)
      if (std::ranges::all_of(terms, [](const auto& t) {
            return t.type_ == ExprType::Negate;
          })) {
        auto newTerms = transform(terms, [](const auto& t) {
          return t.getSingleChild_();
        });  // Seems necessary for WebAssembly
        return ExprFactory::negate(ExprFactory::sum(std::move(newTerms)));
      }

      const auto simplified = ExprFactory::sum(terms);

      // Check whether extracting common factors leads to a shorter expression
      // (xy + xz + xw = x(y + z + w)
      if (distribute) {
        for (const auto leading : std::vector{true, false}) {
          std::map<Expr, size_t> numOccurrencesOfFactors;
          std::vector<Expr> factorPerTerm;
          factorPerTerm.reserve(terms.size());
          for (auto& t : terms) {
            factorPerTerm.push_back(t.getLeadingOrEndingFactor_(leading));
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
                leading
                    ? ExprFactory::product({factor, std::move(sumFactored)})
                    : ExprFactory::product({std::move(sumFactored), factor});
            auto factoredExpr =
                (unfactoredTerms.empty()
                     ? std::move(factorTimesFactored)
                     : ExprFactory::sum(
                           {ExprFactory::sum(std::move(unfactoredTerms)),
                            factorTimesFactored}))
                    .simplify(false);
            if (std::set{ExprType::Sum, ExprType::Product}.contains(
                    factoredExpr.type_)) {
              associativeTransformation(factoredExpr.type_,
                                        factoredExpr.terms_);
            }
            if (factoredExpr.complexity_() < simplified.complexity_()) {
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
    case ExprType::Product: {
      // Recursive simplification
      auto terms = transform(terms_, [distribute](const auto& t) {
        return t.simplify_(distribute);
      });

      // Associative transformation (x(yz) = xyz)
      associativeTransformation(ExprType::Product, terms);

      // Identity transformations (x * 0 = 0, x * 1 = x)
      if (std::ranges::any_of(terms, [](const auto& t) { return t == zero; })) {
        return zero;
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
      std::ranges::partition(
          terms, [](const auto& t) { return t.type_ == ExprType::Number; });

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
            terms, [](const auto& t) { return t.type_ == ExprType::Number; });
        terms.push_back(ExprFactory::number(value));
      }

      const auto simplified = ExprFactory::product(terms);

      // Check whether distributing factor leads to a shorter expression (x(y
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
            if (distributedExpr.complexity_() <= simplified.complexity_()) {
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
}

const Expr& Expr::getSingleChild_() const {
  assert(terms_.size() == 1);
  return terms_.front();
}

Expr Expr::getLeadingOrEndingFactor_(const bool leading) const {
  switch (type_) {
    case ExprType::Number:
    case ExprType::NamedScalar:
    case ExprType::NamedVector:
    case ExprType::Matrix:
    case ExprType::SymmetricMatrix:
    case ExprType::Variable:
    case ExprType::DiagonalMatrix:
    case ExprType::Transpose:
    case ExprType::Invert:
    case ExprType::Log:
      return *this;
    case ExprType::Negate: {
      return getSingleChild_().getLeadingOrEndingFactor_(leading);
    }
    case ExprType::Sum: {
      const auto termFactor = terms_.front().getLeadingOrEndingFactor_(leading);
      if (std::ranges::all_of(terms_, [&](const auto& t) {
            return t.getLeadingOrEndingFactor_(leading) == termFactor;
          })) {
        return termFactor;
      }
      return *this;
    }
    case ExprType::Product: {
      return (leading ? terms_.front() : terms_.back())
          .getLeadingOrEndingFactor_(leading);
    }
    default:
      assert(false);
      throw std::runtime_error("Unhandled enum");
  }
}

Expr Expr::factorOut(const Expr& factor, const bool leading) const {
  assert(getLeadingOrEndingFactor_(leading) == factor);
  if (factor == *this) {
    return unity;
  }

  switch (type_) {
    case ExprType::Number:
    case ExprType::NamedScalar:
    case ExprType::NamedVector:
    case ExprType::Matrix:
    case ExprType::SymmetricMatrix:
    case ExprType::Variable:
    case ExprType::DiagonalMatrix:
    case ExprType::Transpose:
    case ExprType::Invert:
    case ExprType::Log:
      assert(false);  // Should have been handled by factor == *this above
    case ExprType::Negate: {
      {
        return ExprFactory::negate(
            getSingleChild_().factorOut(factor, leading));
      }
    }
    case ExprType::Sum: {
      auto terms = terms_;
      for (auto& t : terms) {
        t = t.factorOut(factor, leading);
      }
      return ExprFactory::sum(std::move(terms));
    }
    case ExprType::Product: {
      if (terms_.size() == 1) {
        assert(terms_.front() == factor);
        return unity;
      }
      auto terms = terms_;
      for (size_t i = 0; i < terms_.size(); ++i) {
        const auto index = leading ? i : terms_.size() - 1 - i;
        const auto& t = terms_[index];
        if (t.getLeadingOrEndingFactor_(leading) == factor) {
          terms[index] = t.factorOut(factor, leading);
          return ExprFactory::product(std::move(terms));
        }
      }
      assert(false);  // Should not be reached
    }
    default:
      break;
  }
  assert(false);
  throw std::runtime_error("Unhandled enum");
}

double Expr::complexity_() const {
  switch (type_) {
    case ExprType::Number:
      return 0.5;
    case ExprType::NamedScalar:
    case ExprType::NamedVector:
    case ExprType::Matrix:
    case ExprType::SymmetricMatrix:
    case ExprType::Variable:
      return 1.0;
    case ExprType::DiagonalMatrix:
    case ExprType::Transpose:
    case ExprType::Negate:
    case ExprType::Invert:
      return 0.5 + getSingleChild_().complexity_();
    case ExprType::Log:
      return 1.0 + getSingleChild_().complexity_();
    case ExprType::Sum:
    case ExprType::Product:
      return std::accumulate(
          terms_.cbegin(), terms_.cend(), 0.0,
          [](const double s, const auto& t) { return s + t.complexity_(); });
    default:
      assert(false);
      throw std::runtime_error("Unhandled enum");
  }
}

Expr ExprFactory::number(const double value) { return Expr(value); }
Expr ExprFactory::namedScalar(const std::string& name) {
  return Expr(ExprType::NamedScalar, name);
}
Expr ExprFactory::namedVector(const std::string& name) {
  return Expr(ExprType::NamedVector, name);
}
Expr ExprFactory::matrix(const std::string& name) {
  return Expr(ExprType::Matrix, name);
}
Expr ExprFactory::symmetricMatrix(const std::string& name) {
  return Expr(ExprType::SymmetricMatrix, name);
}
Expr ExprFactory::variable(const std::string& name) {
  return Expr(ExprType::Variable, name);
}
Expr ExprFactory::diagonalMatrix(Expr expr) {
  return Expr(ExprType::DiagonalMatrix, {std::move(expr)});
}
Expr ExprFactory::transpose(Expr expr) {
  return Expr(ExprType::Transpose, {std::move(expr)});
}
Expr ExprFactory::negate(Expr expr) {
  return Expr(ExprType::Negate, {std::move(expr)});
}
Expr ExprFactory::invert(Expr expr) {
  return Expr(ExprType::Invert, {std::move(expr)});
}
Expr ExprFactory::log(Expr expr) {
  return Expr(ExprType::Log, {std::move(expr)});
}
Expr ExprFactory::product(std::vector<Expr> terms) {
  return Expr(ExprType::Product, std::move(terms));
}
Expr ExprFactory::sum(std::vector<Expr> terms) {
  return Expr(ExprType::Sum, std::move(terms));
}
}  // namespace Expression
