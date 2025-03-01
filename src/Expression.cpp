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
  // Order lexicographically, first by type, then by length of the string
  if (left.getType() != right.getType()) {
    return left.getType() <=> right.getType();
  }
  return left.toExpressionString() <=> right.toExpressionString();
}

bool operator==(const Expr& left, const Expr& right) {
  return left.toString() == right.toString();
}

std::set<Expr> intersectMultiple_(const std::vector<std::set<Expr>>& sets) {
  if (sets.empty()) return {};

  auto result = sets.front();
  for (auto& current : sets | std::views::drop(1)) {
    std::set<Expr> intersection;
    std::ranges::set_intersection(
        result, current, std::inserter(intersection, intersection.begin()));
    result = std::move(intersection);
    if (result.empty()) {
      break;
    }
  }

  return result;
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
                   ExprType::Log, ExprType::Product, ExprType::Sum}
              .contains(type_)));
  if (std::set{ExprType::Transpose, ExprType::Negate, ExprType::Invert,
               ExprType::Log}
          .contains(type_)) {
    assert(terms_.size() == 1);
  }
  assert(!terms_.empty());
}

Expr::Expr(ExprType type, const std::string& name) : type_(type), name_(name) {
  assert(
      (std::set{ExprType::NamedConstant, ExprType::Variable}.contains(type_)));
}

Expr::Expr(const double value) : type_(ExprType::Number), value_(value) {}

Expr Expr::differentiate(const Expr& var) const {
  switch (type_) {
    case ExprType::Number:
    case ExprType::NamedConstant:
      return zero;
    case ExprType::Variable:
      return *this == var ? unity : zero;
    case ExprType::Transpose:
      return getSingleChild_().differentiate(var);
    case ExprType::Negate:
      return ExprFactory::negate(getSingleChild_().differentiate(var));
    case ExprType::Invert:
      assert(false);  // Not implemented
    case ExprType::Log:
      return ExprFactory::product({ExprFactory::invert(getSingleChild_()),
                                   getSingleChild_().differentiate(var)});
    case ExprType::Sum: {
      auto terms = transform(
          terms_, [&var](const auto& t) { return t.differentiate(var); });
      return ExprFactory::sum(std::move(terms));
    }
    case ExprType::Product: {
      std::vector<Expr> sumTerms;
      sumTerms.reserve(terms_.size());
      auto terms = terms_;
      for (size_t i = 0; i < terms.size(); ++i) {
        terms[i] = terms_[i].differentiate(var);
        sumTerms.push_back(ExprFactory::product(terms));
        terms[i] = terms_[i];
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
  auto simplified = expr.simplify_(distribute);
  while (expr != simplified) {
    expr = simplified;
    simplified = expr.simplify_(distribute);
  }
  return simplified;
}

std::string Expr::toString(const bool condensed) const {
  switch (type_) {
    case ExprType::Number: {
      std::stringstream ss;
      ss << value_;
      return ss.str();
    }
    case ExprType::NamedConstant:
    case ExprType::Variable:
      return name_;
    case ExprType::Transpose:
      return getSingleChild_().toString(condensed) + "^T";
    case ExprType::Negate:
      return "-" + getSingleChild_().toString(condensed);
    case ExprType::Invert:
      return getSingleChild_().toString(condensed) + "^{-1}";
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
      ss << terms_.front().toString(condensed);
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
    case ExprType::NamedConstant:
      return "namedConstant(" + name_ + ")";
    case ExprType::Variable:
      return "variable(" + name_ + ")";
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
    case ExprType::NamedConstant:
    case ExprType::Variable:
      return *this;
    case ExprType::Transpose: {
      const auto& child = getSingleChild_();
      if (child.type_ == ExprType::Transpose) {
        // transpose(transpose(x)) = x;
        return child.getSingleChild_().simplify_();
      } else if (child == zero) {
        // 0^T = 0
        return zero;
      } else if (child.type_ == ExprType::Negate) {
        // (-x)^T  = -x^T
        return ExprFactory::negate(
            ExprFactory::transpose(child.getSingleChild_().simplify_()));
      } else if (child.type_ == ExprType::Product) {
        // (xyz)^T = z^T y^T x^T
        auto terms = transform(child.terms_, [](const auto& t) {
          return ExprFactory::transpose(t.simplify_());
        });
        std::ranges::reverse(terms);
        return ExprFactory::product(std::move(terms));
      }
      return ExprFactory::transpose(child.simplify_());
    }
    case ExprType::Negate: {
      const auto& child = getSingleChild_();
      if (child.type_ == ExprType::Negate) {
        // negate(negate(x)) = x
        return child.getSingleChild_().simplify_();
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
        const size_t count = std::ranges::count_if(
            child.terms_,
            [](const auto& t) { return t.type_ == ExprType::Negate; });
        if (count > child.terms_.size() / 2) {
          auto terms = child.terms_;
          for (auto& t : terms) {
            if (t.type_ == ExprType::Negate) {
              t = t.getSingleChild_();
            } else {
              t = ExprFactory::negate(t);
            }
          }
          return ExprFactory::sum(std::move(terms));
        }
      }
      return ExprFactory::negate(child.simplify_());
    }
    case ExprType::Invert: {
      const auto& child = getSingleChild_();
      if (child.type_ == ExprType::Invert) {
        // invert(invert(x)) = x
        return child.getSingleChild_().simplify_();
      } else if (child.type_ == ExprType::Negate) {
        // invert(negate(x)) = negate(invert(x))
        return ExprFactory::negate(
            ExprFactory::invert(child.getSingleChild_().simplify_()));
      } else if (child.type_ == ExprType::Product) {
        // invert(x * y * z) = invert(z) * invert(y) * invert(x)
        auto terms = transform(child.terms_, [](const auto& t) {
          return ExprFactory::invert(t.simplify_());
        });
        std::ranges::reverse(terms);
        return ExprFactory::product(std::move(terms));
      }
      return ExprFactory::invert(child.simplify_());
    }
    case ExprType::Log:
      return ExprFactory::log(getSingleChild_().simplify_());
    case ExprType::Sum: {
      // Recursive simplification
      auto terms =
          transform(terms_, [](const auto& t) { return t.simplify_(); });

      // Distributive transformation (x + y + 1.3x = 2.3x + y)
      for (size_t i = 0; i < terms.size(); ++i) {
        const auto term = terms.at(i);
        const auto isNumberTimesTerm = [&term](const auto& t) {
          return t.type_ == ExprType::Product && t.terms_.size() == 2 &&
                 t.terms_.at(0).type_ == ExprType::Number &&
                 t.terms_.at(1) == term;
        };
        if (std::ranges::count_if(terms,
                                  [&term, &isNumberTimesTerm](const auto& t) {
                                    return t == term || isNumberTimesTerm(t);
                                  }) > 1) {
          const auto value = std::accumulate(
              terms.cbegin(), terms.cend(), 0.0,
              [&term, &isNumberTimesTerm](const double s, const auto& t) {
                return s + (t == term              ? 1.0
                            : isNumberTimesTerm(t) ? t.terms_.at(0).value_
                                                   : 0.0);
              });
          std::erase_if(terms, [&term, &isNumberTimesTerm](const auto& t) {
            return t == term || isNumberTimesTerm(t);
          });
          terms.push_back(
              ExprFactory::product({ExprFactory::number(value), term}));
        }
      }

      // Associative transformation (z + y + x = x + y + z)
      associativeTransformation(ExprType::Sum, terms);

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

      // Commutative transformation
      std::ranges::sort(terms);

      const auto simplified = ExprFactory::sum(terms);

      // Check whether extracting common factors leads to a shorter expression
      // (xy + xz + xw = x(y + z + w)
      if (distribute) {
        std::map<Expr, size_t> numOccurancesOfFactors;
        for (auto& t : terms) {
          for (const auto& f : t.getUniqueFactors_()) {
            ++numOccurancesOfFactors[f];
          }
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
            if (t.getUniqueFactors_().contains(factor)) {
              factoredTerms.push_back(t.factorOut(factor));
            } else {
              unfactoredTerms.push_back(t);
            }
          }

          auto factoredExpr =
              unfactoredTerms.empty()
                  ? ExprFactory::product(
                        {factor, ExprFactory::sum(std::move(factoredTerms))})
                  : ExprFactory::sum(
                        {ExprFactory::sum(std::move(unfactoredTerms)),
                         ExprFactory::product(
                             {factor,
                              ExprFactory::sum(std::move(factoredTerms))})});
          associativeTransformation(factoredExpr.type_, factoredExpr.terms_);
          if (factoredExpr.complexity_() < simplified.complexity_()) {
            return factoredExpr;
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
      auto terms =
          transform(terms_, [](const auto& t) { return t.simplify_(); });

      // Associative transformation (zyx = xyz)
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
        auto newTerms = terms;
        for (size_t i = 0; i < terms.size(); ++i) {
          if (terms[i].type_ == ExprType::Negate) {
            newTerms[i] = terms[i].getSingleChild_();
          }
          return ExprFactory::negate(ExprFactory::product(std::move(newTerms)));
        }
      }

      // Canceling transformation (power transformation for the special case
      // of inverse) (x * inv(x) = 1)
      eraseCanceling(ExprType::Invert, terms, unity);

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

      // Commutative transformation
      // std::ranges::sort(terms); // Removed to facilitate matrix algebra

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
}

const Expr& Expr::getSingleChild_() const {
  assert(terms_.size() == 1);
  return terms_.front();
}

std::set<Expr> Expr::getUniqueFactors_() const {
  switch (type_) {
    case ExprType::Number:
    case ExprType::NamedConstant:
    case ExprType::Variable:
    case ExprType::Transpose:
    case ExprType::Invert:
    case ExprType::Log:
      return std::set<Expr>{*this};
    case ExprType::Negate:
      return getSingleChild_().getUniqueFactors_();
    case ExprType::Sum: {
      std::vector<std::set<Expr>> factorsPerTerm(terms_.size());
      std::ranges::transform(terms_, factorsPerTerm.begin(), [](const auto& t) {
        return t.getUniqueFactors_();
      });
      const auto factors = intersectMultiple_(factorsPerTerm);
      return factors.empty() ? std::set<Expr>{*this} : factors;
    }
    case ExprType::Product: {
      std::set<Expr> factors;
      for (const auto& t : terms_) {
        const auto childFactors = t.getUniqueFactors_();
        factors.insert(childFactors.begin(), childFactors.end());
      }
      return factors;
    }
    default:
      assert(false);
      throw std::runtime_error("Unhandled enum");
  }
}

Expr Expr::factorOut(const Expr& factor) const {
  assert(getUniqueFactors_().contains(factor));

  switch (type_) {
    case ExprType::Number:
    case ExprType::NamedConstant:
    case ExprType::Variable:
    case ExprType::Transpose:
    case ExprType::Invert:
    case ExprType::Log:
      assert(factor == *this);
      return unity;
    case ExprType::Negate: {
      return ExprFactory::negate(getSingleChild_().factorOut(factor));
    }
    case ExprType::Sum: {
      auto terms = terms_;
      for (auto& t : terms) {
        t = t.factorOut(factor);
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
        const auto& t = terms_[i];
        if (t == factor) {
          terms.erase(terms.begin() + i);
          return ExprFactory::product(std::move(terms));
        }
      }
      for (auto& t : terms) {
        const auto childFactors = t.getUniqueFactors_();
        if (std::ranges::any_of(childFactors,
                                [&factor](const auto& childFactor) {
                                  return childFactor == factor;
                                })) {
          t = t.factorOut(factor);
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
    case ExprType::NamedConstant:
    case ExprType::Variable:
      return 1.0;
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
Expr ExprFactory::namedConstant(const std::string& name) {
  return Expr(ExprType::NamedConstant, name);
}
Expr ExprFactory::variable(const std::string& name) {
  return Expr(ExprType::Variable, name);
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
