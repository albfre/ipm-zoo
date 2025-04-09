#pragma once
#include "Expr.h"
#include "NumericalOptimization/Evaluation.h"
#include "SymbolicOptimization.h"

namespace NumericalOptimization {
struct Data {
  std::vector<std::vector<double>> Q;
  std::vector<double> c;
  std::vector<std::vector<double>> A_ineq;
  std::vector<double> l_A_ineq;
  std::vector<double> u_A_ineq;
  std::vector<std::vector<double>> A_eq;
  std::vector<double> b_eq;
  std::vector<double> l_x;
  std::vector<double> u_x;
};

Evaluation::Environment build_environment(
    const SymbolicOptimization::VariableNames& names, const Data& data);

struct ScopedEnvironmentOverride {
  ScopedEnvironmentOverride(Evaluation::Environment& env,
                            Expression::ExprPtr var,
                            Evaluation::EvalResult temp_val);
  ScopedEnvironmentOverride(ScopedEnvironmentOverride&) = delete;
  ~ScopedEnvironmentOverride();

 private:
  Evaluation::Environment& env_;
  Expression::ExprPtr var_;
  Evaluation::EvalResult orig_val_;
};

}  // namespace NumericalOptimization