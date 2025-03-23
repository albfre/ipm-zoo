#include "ExprFactory.h"

#include "Helpers.h"

namespace Expression {
Expr ExprFactory::number(const double value) { return Expr(Number(value)); }

Expr ExprFactory::namedScalar(std::string_view name) {
  return Expr(NamedScalar(name));
}

Expr ExprFactory::namedVector(std::string_view name) {
  return Expr(NamedVector(name));
}

Expr ExprFactory::variable(std::string_view name) {
  return Expr(Variable(name));
}

Expr ExprFactory::matrix(std::string_view name) { return Expr(Matrix(name)); }

Expr ExprFactory::symmetricMatrix(std::string_view name) {
  return Expr(SymmetricMatrix(name));
}

Expr ExprFactory::diagonalMatrix(Expr expr) {
  return Expr(DiagonalMatrix(std::make_shared<Expr>(std::move(expr))));
}

Expr ExprFactory::transpose(Expr expr) {
  return Expr(Transpose(std::make_shared<Expr>(std::move(expr))));
}

Expr ExprFactory::negate(Expr expr) {
  return Expr(Negate(std::make_shared<Expr>(std::move(expr))));
}

Expr ExprFactory::invert(Expr expr) {
  return Expr(Invert(std::make_shared<Expr>(std::move(expr))));
}

Expr ExprFactory::log(Expr expr) {
  return Expr(Log(std::make_shared<Expr>(std::move(expr))));
}

Expr ExprFactory::sum(std::vector<Expr> terms) {
  if (terms.empty()) {
    return zero;
  }
  if (terms.size() == 1) {
    return std::move(terms[0]);
  }
  return Expr(Sum(std::move(terms)));
}

Expr ExprFactory::product(std::vector<Expr> terms) {
  if (terms.empty()) {
    return zero;
  }
  if (terms.size() == 1) {
    return std::move(terms[0]);
  }
  return Expr(Product(std::move(terms)));
}
}  // namespace Expression