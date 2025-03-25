#include "ExprFactory.h"

#include "Helpers.h"
#include "ToStringVisitor.h"

namespace Expression {
ExprFactory::ExprFactory() : mutex_(std::make_unique<std::mutex>()) {}

ExprFactory& ExprFactory::instance_() {
  static ExprFactory factory;
  return factory;
}

ExprPtr ExprFactory::getExpr_(Expr::ExprVariant&& variant) {
  std::scoped_lock<std::mutex> lock(*mutex_);

  const auto key = std::visit(ToExpressionStringVisitor{}, variant);
  auto it = cache_.find(key);
  if (it != cache_.end()) {
    if (auto expr = it->second.lock()) {
      return expr;
    }
    cache_.erase(it);
  }
  auto expr = ExprPtr(new Expr(std::move(variant), key));
  cache_[key] = expr;

  // Clean cache
  if (++cleanupCounter_ % 1000 == 0) {
    for (auto it = cache_.begin(); it != cache_.end();) {
      if (it->second.expired()) {
        it = cache_.erase(it);
      } else {
        ++it;
      }
    }
  }
  return expr;
}

ExprPtr ExprFactory::number(const double value) {
  return instance_().getExpr_(Number(value));
}

ExprPtr ExprFactory::namedScalar(std::string_view name) {
  return instance_().getExpr_(NamedScalar(name));
}

ExprPtr ExprFactory::namedVector(std::string_view name) {
  return instance_().getExpr_(NamedVector(name));
}

ExprPtr ExprFactory::variable(std::string_view name) {
  return instance_().getExpr_(Variable(name));
}

ExprPtr ExprFactory::matrix(std::string_view name) {
  return instance_().getExpr_(Matrix(name));
}

ExprPtr ExprFactory::symmetricMatrix(std::string_view name) {
  return instance_().getExpr_(SymmetricMatrix(name));
}

ExprPtr ExprFactory::diagonalMatrix(ExprPtr expr) {
  return instance_().getExpr_(DiagonalMatrix(std::move(expr)));
}

ExprPtr ExprFactory::transpose(ExprPtr expr) {
  return instance_().getExpr_(Transpose(std::move(expr)));
}

ExprPtr ExprFactory::negate(ExprPtr expr) {
  return instance_().getExpr_(Negate(std::move(expr)));
}

ExprPtr ExprFactory::invert(ExprPtr expr) {
  return instance_().getExpr_(Invert(std::move(expr)));
}

ExprPtr ExprFactory::log(ExprPtr expr) {
  return instance_().getExpr_(Log(std::move(expr)));
}

ExprPtr ExprFactory::sum(std::vector<ExprPtr> terms) {
  if (terms.empty()) {
    return zero;
  }
  if (terms.size() == 1) {
    return std::move(terms[0]);
  }
  return instance_().getExpr_(Sum(std::move(terms)));
}

ExprPtr ExprFactory::product(std::vector<ExprPtr> terms) {
  if (terms.empty()) {
    return unity;
  }
  if (terms.size() == 1) {
    return std::move(terms[0]);
  }
  return instance_().getExpr_(Product(std::move(terms)));
}

ExprPtr ExprFactory::asPtr(const Expr& expr) {
  auto variant = expr.getImpl();
  return instance_().getExpr_(std::move(variant));
}

ExprPtr ExprFactory::getExpr(Expr::ExprVariant&& variant) {
  return instance_().getExpr_(std::move(variant));
}
}  // namespace Expression