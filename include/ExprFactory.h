#pragma once
#include <mutex>
#include <unordered_map>

#include "Expr.h"

namespace Expression {
class ExprFactory {
 public:
  [[nodiscard]] static ExprPtr number(const double value);
  [[nodiscard]] static ExprPtr named_scalar(std::string_view name);
  [[nodiscard]] static ExprPtr named_vector(std::string_view name);
  [[nodiscard]] static ExprPtr variable(std::string_view name);
  [[nodiscard]] static ExprPtr matrix(std::string_view name);
  [[nodiscard]] static ExprPtr symmetric_matrix(std::string_view name);
  [[nodiscard]] static ExprPtr diagonal_matrix(ExprPtr expr);
  [[nodiscard]] static ExprPtr transpose(ExprPtr expr);
  [[nodiscard]] static ExprPtr negate(ExprPtr expr);
  [[nodiscard]] static ExprPtr invert(ExprPtr expr);
  [[nodiscard]] static ExprPtr log(ExprPtr expr);
  [[nodiscard]] static ExprPtr sum(std::vector<ExprPtr> terms);
  [[nodiscard]] static ExprPtr product(std::vector<ExprPtr> terms);
  [[nodiscard]] static ExprPtr as_ptr(const Expr& expr);
  [[nodiscard]] static ExprPtr get_expr(Expr::ExprVariant&& var);

 private:
  friend class ExprFactoryTest;

  ExprFactory();
  static ExprFactory& instance_();
  ExprPtr get_expr_(Expr::ExprVariant&& variant);
  std::unordered_map<std::string, std::weak_ptr<const Expr>> cache_;
  std::unique_ptr<std::mutex> mutex_;
  size_t cleanup_counter_ = 0;
};
}  // namespace Expression