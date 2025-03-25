#pragma once
#include <gtest/gtest_prod.h>

#include <mutex>
#include <unordered_map>

#include "Expr.h"

namespace Expression {
class ExprFactory {
 public:
  [[nodiscard]] static ExprPtr number(const double value);
  [[nodiscard]] static ExprPtr namedScalar(std::string_view name);
  [[nodiscard]] static ExprPtr namedVector(std::string_view name);
  [[nodiscard]] static ExprPtr variable(std::string_view name);
  [[nodiscard]] static ExprPtr matrix(std::string_view name);
  [[nodiscard]] static ExprPtr symmetricMatrix(std::string_view name);
  [[nodiscard]] static ExprPtr diagonalMatrix(ExprPtr expr);
  [[nodiscard]] static ExprPtr transpose(ExprPtr expr);
  [[nodiscard]] static ExprPtr negate(ExprPtr expr);
  [[nodiscard]] static ExprPtr invert(ExprPtr expr);
  [[nodiscard]] static ExprPtr log(ExprPtr expr);
  [[nodiscard]] static ExprPtr sum(std::vector<ExprPtr> terms);
  [[nodiscard]] static ExprPtr product(std::vector<ExprPtr> terms);
  [[nodiscard]] static ExprPtr asPtr(const Expr& expr);
  [[nodiscard]] static ExprPtr getExpr(Expr::ExprVariant&& var);

 private:
  friend class ExprFactoryTest;

  ExprFactory();
  static ExprFactory& instance_();
  ExprPtr getExpr_(Expr::ExprVariant&& variant);
  std::unordered_map<std::string, std::weak_ptr<const Expr>> cache_;
  std::unique_ptr<std::mutex> mutex_;
  size_t cleanupCounter_ = 0;
};
}  // namespace Expression