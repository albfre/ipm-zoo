#pragma once
#include <algorithm>
#include <cassert>
#include <ranges>
#include "Expression.h"
#include "Helpers.h"

namespace Expression {
// Visitor for differentiation
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
  // Transpose: d/dx(f(x)^T) = (d/dx(f(x)))^T
  Expr operator()(const Transpose& x) const {
    return ExprFactory::transpose(x.child->differentiate(var));
  }
  // Negation: d/dx(-f(x)) = -d/dx(f(x))
  Expr operator()(const Negate& x) const {
    return ExprFactory::negate(x.child->differentiate(var));
  }
  // Inverse: Not yet implemented
  Expr operator()(const Invert&) const {
    assert(false);  // Not implemented
    return zero;
  }
  // Logarithm: d/dx(log(f(x))) = (1/f(x)) * d/dx(f(x))
  Expr operator()(const Log& x) const {
    return ExprFactory::product(
        {ExprFactory::invert(*x.child), x.child->differentiate(var)});
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
        sumTerms.push_back(ExprFactory::product(std::move(terms)));
      }

      // Handle special case for transpose
      if (x.terms.size() > 1 &&
          std::visit(
              overloaded{
                  [&](const Transpose& y) {
                    return std::visit(
                        overloaded{
                            [](const Matrix& z) { return false; },
                            [&](const auto& z) {
                              // xi is a transpose with a non-matrix child

                              // d/dx(f(x)^T g(x)) = d/dx(f(x)^T) g(x) +
                              // d/dx(g(x))^T f(x)
                              auto terms = std::vector(x.terms.begin(),
                                                       x.terms.begin() + i);
                              auto rest = std::vector(x.terms.begin() + i + 1,
                                                      x.terms.end());

                              terms.push_back(
                                  ExprFactory::transpose(
                                      ExprFactory::product(std::move(rest)))
                                      .differentiate(var));
                              terms.push_back(*y.child);

                              sumTerms.push_back(
                                  ExprFactory::product(std::move(terms)));

                              return true;
                            }},
                        y.child->getImpl());
                  },
                  [](const auto&) { return false; }},
              xi.getImpl())) {
        break;
      }
    }

    return ExprFactory::sum(std::move(sumTerms));
  }
};

}  // namespace Expression