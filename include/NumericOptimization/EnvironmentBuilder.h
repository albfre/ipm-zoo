#pragma once
#include "Expr.h"
#include "NumericOptimization/Evaluation.h"
#include "SymbolicOptimization.h"

namespace NumericOptimization {
struct Data {
  std::vector<std::vector<double>> A_eq;
  std::vector<double> b_eq;
  double delta_eq = 0.0;
  std::vector<std::vector<double>> A_ineq;
  std::vector<double> l_A_ineq;
  std::vector<double> u_A_ineq;
  std::vector<double> l_x;
  std::vector<double> u_x;
  std::vector<std::vector<double>> Q;
  std::vector<double> c;
};

Evaluation::Environment build_environment(
    const SymbolicOptimization::VariableNames& names,
    const NumericOptimization::Data& data) {
  using namespace Expression;
  auto x = ExprFactory::variable(names.x);
  auto A_eq = ExprFactory::matrix(names.A_eq);
  auto b_eq = ExprFactory::named_vector(names.b_eq);
  auto p_eq = ExprFactory::variable(names.p_eq);
  auto delta_eq = ExprFactory::named_scalar(names.delta_eq);
  auto A_ineq = ExprFactory::matrix(names.A_ineq);
  auto s_A_ineq = ExprFactory::variable(names.s_A_ineq);
  auto s_A_ineq_l = ExprFactory::variable(names.s_A_ineq_l);
  auto s_A_ineq_u = ExprFactory::variable(names.s_A_ineq_u);
  auto s_x_l = ExprFactory::variable(names.s_x_l);
  auto s_x_u = ExprFactory::variable(names.s_x_u);
  auto s_A_eq = ExprFactory::variable(names.s_A_eq);
  auto s_A_eq_l = ExprFactory::variable(names.s_A_eq_l);
  auto s_A_eq_u = ExprFactory::variable(names.s_A_eq_u);
  auto l_A_ineq = ExprFactory::named_vector(names.l_A_ineq);
  auto u_A_ineq = ExprFactory::named_vector(names.u_A_ineq);
  auto l_x = ExprFactory::named_vector(names.l_x);
  auto u_x = ExprFactory::named_vector(names.u_x);
  auto Q = ExprFactory::symmetric_matrix(names.Q);
  auto c = ExprFactory::named_vector(names.c);

  Evaluation::Environment env;
  env[A_eq] = Evaluation::val_matrix(data.A_eq);
  env[b_eq] = Evaluation::val_vector(data.b_eq);
  env[delta_eq] = Evaluation::val_scalar(data.delta_eq);
  env[A_ineq] = Evaluation::val_matrix(data.A_ineq);
  env[l_A_ineq] = Evaluation::val_vector(data.l_A_ineq);
  env[u_A_ineq] = Evaluation::val_vector(data.u_A_ineq);
  env[l_x] = Evaluation::val_vector(data.l_x);
  env[Q] = Evaluation::val_matrix(data.Q);
  env[c] = Evaluation::val_vector(data.c);

  const auto n_variables = data.Q.size();
  const auto n_ineq = data.A_ineq.size();
  const auto n_eq = data.A_eq.size();
  return env;
}

}  // namespace NumericOptimization