#pragma once
#include "Evaluation.h"
#include "Expr.h"
#include "NumericOptimization/LinearSolvers.h"
#include "SymbolicOptimization.h"
#include "Utils/Assert.h"
#include "Utils/Helpers.h"

namespace NumericOptimization {
using Matrix = std::vector<std::vector<double>>;
class NumericOptimization {
  NumericOptimization(Evaluation::Environment& env,
                      SymbolicOptimization::NewtonSystem newton_system)
      : env_(env),
        newton_system_(std::move(newton_system)),
        augmented_system_(
            SymbolicOptimization::get_augmented_system(newton_system_)),
        normal_equations_(
            SymbolicOptimization::get_normal_equations(augmented_system_)) {}

  void solve() {
    const auto& lhs = augmented_system_.lhs;
    const auto is_indefinite = std::ranges::any_of(
        std::views::iota(0ul, lhs.size()),
        [&lhs](const auto i) { return lhs.at(i).at(i) == Expression::zero; });
    if (is_indefinite) {
      solve_indefinite_();
    } else {
      solve_quasi_definite_();
    }
  }

 private:
  void solve_indefinite_() { ASSERT(false); }
  void solve_quasi_definite_() {
    const auto& [lhs, rhs, variables, delta_definitions] = augmented_system_;
  }

  std::vector<Evaluation::EvalResult> evalExprVector_(
      const std::vector<Expression::ExprPtr>& v) {
    std::vector<Evaluation::EvalResult> result;
    result.reserve(v.size());
    std::ranges::transform(
        v, std::back_inserter(result),
        [this](const auto& expr) { return Evaluation::evaluate(expr, env_); });
    return result;
  }

  std::vector<std::vector<Evaluation::EvalResult>> evalExprMatrix_(
      const std::vector<std::vector<Expression::ExprPtr>>& v) {
    std::vector<std::vector<Evaluation::EvalResult>> result;
    result.reserve(v.size());
    std::ranges::transform(
        v, std::back_inserter(result),
        [this](const auto& exprVec) { return evalExprVector_(exprVec); });
    return result;
  }

  Evaluation::Environment& env_;
  Expression::ExprPtr objective_;
  SymbolicOptimization::NewtonSystem newton_system_;
  SymbolicOptimization::NewtonSystem augmented_system_;
  SymbolicOptimization::NewtonSystem normal_equations_;
};

}  // namespace NumericOptimization