#pragma once
#include "Expr.h"
#include "NumericOptimization/Evaluation.h"
#include "SymbolicOptimization.h"

namespace NumericOptimization {
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
    const SymbolicOptimization::VariableNames& names,
    const NumericOptimization::Data& data) {
  using namespace Expression;

  const auto o = SymbolicOptimization::get_optimization_expressions(names);

  Evaluation::Environment env;

  // Set constants provided in the data struct
  env[o.Q] = Evaluation::val_matrix(data.Q);
  env[o.c] = Evaluation::val_vector(data.c);
  env[o.A_ineq] = Evaluation::val_matrix(data.A_ineq);
  env[o.l_A_ineq] = Evaluation::val_vector(data.l_A_ineq);
  env[o.u_A_ineq] = Evaluation::val_vector(data.u_A_ineq);
  env[o.A_eq] = Evaluation::val_matrix(data.A_eq);
  env[o.b_eq] = Evaluation::val_vector(data.b_eq);
  env[o.l_x] = Evaluation::val_vector(data.l_x);

  const auto vars = std::vector(data.Q.size(), 1.0);
  const auto ineqs = std::vector(data.A_ineq.size(), 1.0);
  const auto eqs = std::vector(data.A_eq.size(), 1.0);

  // Set other constants
  env[o.delta_eq] = Evaluation::val_scalar(1e-4);
  env[o.mu] = Evaluation::val_scalar(1.0);
  env[o.e_var] = Evaluation::val_vector(vars);
  env[o.e_ineq] = Evaluation::val_vector(ineqs);
  env[o.e_eq] = Evaluation::val_vector(eqs);

  // Set initial values for variables and slacks
  env[o.x] = Evaluation::val_vector(vars);
  env[o.p_eq] = Evaluation::val_vector(eqs);
  env[o.s_A_ineq] = Evaluation::val_vector(ineqs);
  env[o.s_A_ineq_l] = Evaluation::val_vector(ineqs);
  env[o.s_A_ineq_u] = Evaluation::val_vector(ineqs);
  env[o.s_x_l] = Evaluation::val_vector(vars);
  env[o.s_x_u] = Evaluation::val_vector(vars);
  env[o.s_A_eq] = Evaluation::val_vector(eqs);
  env[o.s_A_eq_l] = Evaluation::val_vector(eqs);
  env[o.s_A_eq_u] = Evaluation::val_vector(eqs);

  env[o.lambda_A_eq] = Evaluation::val_vector(eqs);
  env[o.lambda_sAeql] = Evaluation::val_vector(eqs);
  env[o.lambda_sAequ] = Evaluation::val_vector(eqs);
  env[o.lambda_A_ineq] = Evaluation::val_vector(ineqs);
  env[o.lambda_sAineql] = Evaluation::val_vector(ineqs);
  env[o.lambda_sAinequ] = Evaluation::val_vector(ineqs);
  env[o.lambda_sxl] = Evaluation::val_vector(vars);
  env[o.lambda_sxu] = Evaluation::val_vector(vars);

  return env;
}

}  // namespace NumericOptimization