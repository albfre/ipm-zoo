#include "DifferentiationVisitor.h"

#include <algorithm>
#include <iostream>
#include <ranges>

#include "Assert.h"

namespace Expression {

DifferentiationVisitor::DifferentiationVisitor(const ExprPtr& var) : var_(var) {
  ASSERT(is<Variable>(var));
}

ExprPtr DifferentiationVisitor::operator()(const Variable& x) const {
  return ExprFactory::variable(x.name) == var_ ? unity : zero;
}

ExprPtr DifferentiationVisitor::operator()(const DiagonalMatrix& x) const {
  return ExprFactory::diagonalMatrix(x.child->differentiate(var_));
}

// Transpose: d/dx(f(x)^T) = (d/dx(f(x)))^T
ExprPtr DifferentiationVisitor::operator()(const Transpose& x) const {
  return ExprFactory::transpose(x.child->differentiate(var_));
}

// Negation: d/dx(-f(x)) = -d/dx(f(x))
ExprPtr DifferentiationVisitor::operator()(const Negate& x) const {
  return ExprFactory::negate(x.child->differentiate(var_));
}

// Inverse: Not yet implemented
ExprPtr DifferentiationVisitor::operator()(const Invert& x) const {
  ASSERT(false);  // Not implemented
  return zero;
}

// Logarithm: d/dx(log(f(x))) = (1/f(x)) * d/dx(f(x))
ExprPtr DifferentiationVisitor::operator()(const Log& x) const {
  const auto& child = x.child;
  return ExprFactory::product(
      {ExprFactory::invert(ExprFactory::diagonalMatrix(child)),
       child->differentiate(var_)});
}

// Sum: d/dx(f(x) + g(x)) = d/dx(f(x)) + d/dx(g(x))
ExprPtr DifferentiationVisitor::operator()(const Sum& x) const {
  auto terms = transform(
      x.terms, [this](const auto& t) { return t->differentiate(var_); });
  return ExprFactory::sum(std::move(terms));
}

// Product: d/dx(f(x) * g(x)) = d/dx(f(x)) * g(x) + f(x) * d/dx(g(x))
ExprPtr DifferentiationVisitor::operator()(const Product& x) const {
  std::vector<ExprPtr> sumTerms;
  sumTerms.reserve(x.terms.size());

  // Handle general product rule cases
  for (size_t i = 0; i < x.terms.size(); ++i) {
    const auto& xi = x.terms[i];
    {
      auto terms = x.terms;
      terms[i] = xi->differentiate(var_);  // Differentiate one term
      if (i + 2 == terms.size() && is<DiagonalMatrix>(terms[i]) &&
          is<Variable>(terms[i + 1])) {
        terms[i + 1] = ExprFactory::diagonalMatrix(terms[i + 1]);
      }
      sumTerms.push_back(ExprFactory::product(std::move(terms)));
    }

    // Handle special case for transpose
    if (x.terms.size() > i + 1 &&
        match(xi).with(
            [&](const Transpose& y) {
              return match(*y.child).with(
                  [](const Matrix& z) { return false; },
                  [&](const auto& z) {
                    // xi is a transpose with a non-matrix child
                    // Its derivative is d/dx(f(x)^T g(x)) = d/dx(f(x)^T) g(x)
                    // + d/dx(g(x))^T f(x)
                    // where the first term has already been added above.
                    auto terms =
                        std::vector(x.terms.begin(), x.terms.begin() + i);
                    auto restTerm =
                        i + 2 == x.terms.size()
                            ? x.terms[i + 1]
                            : ExprFactory::product(std::vector(
                                  x.terms.begin() + i + 1, x.terms.end()));
                    terms.push_back(ExprFactory::transpose(std::move(restTerm))
                                        ->differentiate(var_));  // d/dx(g(x))^T
                    terms.push_back(y.child);                    // f(x)
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

}  // namespace Expression