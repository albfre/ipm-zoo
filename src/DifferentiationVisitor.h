#pragma once
#include <algorithm>
#include <cassert>
#include <ranges>

#include "Expression.h"
#include "Helpers.h"

namespace Expression {
struct DifferentiationVisitor {
  const Expr& var;

  explicit DifferentiationVisitor(const Expr& variable) : var(variable) {}
  Expr operator()(const Number&) const { return zero; }
  Expr operator()(const NamedConstant&) const { return zero; }
  Expr operator()(const Matrix&) const { return zero; }
  Expr operator()(const SymmetricMatrix&) const { return zero; }
  Expr operator()(const Variable& x) const {
    return Expr(x) == var ? unity : zero;
  }
  Expr operator()(const DiagonalMatrix& x) const {
    return ExprFactory::diagonalMatrix(x.child->differentiate(var));
  }
  // Transpose: d/dx(f(x)^T) = (d/dx(f(x)))^T
  Expr operator()(const Transpose& x) const {
    return ExprFactory::transpose(x.child->differentiate(var));
  }
  // Negation: d/dx(-f(x)) = -d/dx(f(x))
  Expr operator()(const Negate& x) const {
    return ExprFactory::negate(x.child->differentiate(var));
  }
  // Inverse: Not yet implemented
  Expr operator()(const Invert& x) const {
    assert(false);  // Not implemented
    return zero;
  }
  // Logarithm: d/dx(log(f(x))) = (1/f(x)) * d/dx(f(x))
  Expr operator()(const Log& x) const {
    const auto& child = *x.child;
    return ExprFactory::product(
        {ExprFactory::invert(ExprFactory::diagonalMatrix(child)),
         child.differentiate(var)});
  }
  // Sum: d/dx(f(x) + g(x)) = d/dx(f(x)) + d/dx(g(x))
  Expr operator()(const Sum& x) const {
    auto terms = transform(
        x.terms, [this](const auto& t) { return t.differentiate(var); });
    return ExprFactory::sum(std::move(terms));
  }
  // Product: d/dx(f(x) * g(x)) = d/dx(f(x)) * g(x) + f(x) * d/dx(g(x))
  Expr operator()(const Product& x) const {
    std::vector<Expr> sumTerms;
    sumTerms.reserve(x.terms.size());

    // Handle general product rule cases
    for (size_t i = 0; i < x.terms.size(); ++i) {
      const auto& xi = x.terms[i];
      {
        auto terms = x.terms;
        terms[i] = xi.differentiate(var);  // Differentiate one term
        if (i + 2 == terms.size() && is<DiagonalMatrix>(terms[i]) &&
            is<Variable>(terms[i + 1])) {
          terms[i + 1] = ExprFactory::diagonalMatrix(terms[i + 1]);
        }
        sumTerms.push_back(ExprFactory::product(std::move(terms)));
      }

      // Handle special case for transpose
      if (x.terms.size() > i + 1 &&
          match(
              xi,
              [&](const Transpose& y) {
                return match(
                    *y.child, [](const Matrix& z) { return false; },
                    [&](const auto& z) {
                      // xi is a transpose with a non-matrix child

                      // d/dx(f(x)^T g(x)) = d/dx(f(x)^T) g(x) + d/dx(g(x))^T
                      // f(x)
                      auto terms =
                          std::vector(x.terms.begin(), x.terms.begin() + i);
                      auto restTerm =
                          i + 2 == x.terms.size()
                              ? x.terms[i + 1]
                              : ExprFactory::product(std::vector(
                                    x.terms.begin() + i + 1, x.terms.end()));
                      terms.push_back(
                          ExprFactory::transpose(std::move(restTerm))
                              .differentiate(var));  // d/dx(g(x))^T
                      terms.push_back(*y.child);     // f(x)
                      sumTerms.push_back(ExprFactory::product(
                          std::move(terms)));  // d/dx(g(x))^T f(x)
                      return true;
                    });
              },
              [](const auto&) { return false; })) {
        break;
      }
    }

    return ExprFactory::sum(std::move(sumTerms));
  }
};

}  // namespace Expression