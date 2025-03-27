#include "ExprFactory.h"

#include "Utils/Helpers.h"
#include "Visitors/ToStringVisitor.h"

namespace Expression {
ExprFactory::ExprFactory() : mutex_(std::make_unique<std::mutex>()) {}

ExprFactory& ExprFactory::instance_() {
  static ExprFactory factory;
  return factory;
}

ExprPtr ExprFactory::get_expr_(Expr::ExprVariant&& variant) {
  std::scoped_lock<std::mutex> lock(*mutex_);

  const auto key = std::visit(ToExpressionStringVisitor{}, variant);
  if (const auto it = cache_.find(key); it != cache_.end()) {
    if (auto expr = it->second.lock()) {
      return expr;
    }
    cache_.erase(it);
  }
  auto expr = ExprPtr(new Expr(std::move(variant), key));
  cache_[key] = expr;

  // Clean cache
  if (++cleanup_counter_ >= 1000) {
    std::erase_if(cache_, [](const auto& it) { return it.second.expired(); });
    cleanup_counter_ = 0;
  }
  return expr;
}

ExprPtr ExprFactory::number(const double value) {
  return instance_().get_expr_(Number(value));
}

ExprPtr ExprFactory::named_scalar(std::string_view name) {
  return instance_().get_expr_(NamedScalar(name));
}

ExprPtr ExprFactory::named_vector(std::string_view name) {
  return instance_().get_expr_(NamedVector(name));
}

ExprPtr ExprFactory::variable(std::string_view name) {
  return instance_().get_expr_(Variable(name));
}

ExprPtr ExprFactory::matrix(std::string_view name) {
  return instance_().get_expr_(Matrix(name));
}

ExprPtr ExprFactory::symmetric_matrix(std::string_view name) {
  return instance_().get_expr_(SymmetricMatrix(name));
}

ExprPtr ExprFactory::diagonal_matrix(ExprPtr expr) {
  return instance_().get_expr_(DiagonalMatrix(std::move(expr)));
}

ExprPtr ExprFactory::transpose(ExprPtr expr) {
  return instance_().get_expr_(Transpose(std::move(expr)));
}

ExprPtr ExprFactory::negate(ExprPtr expr) {
  return instance_().get_expr_(Negate(std::move(expr)));
}

ExprPtr ExprFactory::invert(ExprPtr expr) {
  return instance_().get_expr_(Invert(std::move(expr)));
}

ExprPtr ExprFactory::log(ExprPtr expr) {
  return instance_().get_expr_(Log(std::move(expr)));
}

ExprPtr ExprFactory::sum(std::vector<ExprPtr> terms) {
  if (terms.empty()) {
    return zero;
  }
  if (terms.size() == 1) {
    return std::move(terms[0]);
  }
  return instance_().get_expr_(Sum(std::move(terms)));
}

ExprPtr ExprFactory::product(std::vector<ExprPtr> terms) {
  if (terms.empty()) {
    return unity;
  }
  if (terms.size() == 1) {
    return std::move(terms[0]);
  }
  return instance_().get_expr_(Product(std::move(terms)));
}

ExprPtr ExprFactory::as_ptr(const Expr& expr) {
  auto variant = expr.get_impl();
  return instance_().get_expr_(std::move(variant));
}

ExprPtr ExprFactory::get_expr(Expr::ExprVariant&& variant) {
  return instance_().get_expr_(std::move(variant));
}
}  // namespace Expression