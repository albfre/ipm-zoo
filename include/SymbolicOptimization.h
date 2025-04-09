#pragma once
#include "Expr.h"

namespace SymbolicOptimization {
struct VariableNames {
  std::string x = "x";
  std::string A_eq = "C";
  std::string b_eq = "d";
  std::string p_eq = "p";
  std::string delta_eq = "\\delta";
  std::string A_ineq = "A";
  std::string s_A_ineq = "s";
  std::string s_A_ineq_l = "g";
  std::string s_A_ineq_u = "h";
  std::string s_x_l = "y";
  std::string s_x_u = "z";
  std::string s_A_eq = "t";
  std::string s_A_eq_l = "v";
  std::string s_A_eq_u = "w";
  std::string l_A_ineq = "l_A";
  std::string u_A_ineq = "u_A";
  std::string l_x = "l_x";
  std::string u_x = "u_x";
  std::string Q = "Q";
  std::string c = "c";
};

enum class Bounds {
  None,
  Lower,
  Upper,
  Both,
};

enum class InequalityHandling {
  Slacks,
  SlackedSlacks,
  NaiveSlacks,
};

enum class EqualityHandling {
  None,
  Slacks,
  SlackedSlacks,
  NaiveSlacks,
  PenaltyFunction,
  PenaltyFunctionWithExtraDual,
  Regularization,
};

enum class OptimizationProblemType {
  Original,
  Slacked,
  SlackedWithBarriers,
  ForOptimalityConditions,
};

struct Settings {
  Bounds inequalities = Bounds::Both;
  Bounds variable_bounds = Bounds::Both;
  bool equalities = false;
  EqualityHandling equality_handling = EqualityHandling::None;
  InequalityHandling inequality_handling = InequalityHandling::SlackedSlacks;
};

struct InequalityConstraint {
  Expression::ExprPtr expr;
  Expression::ExprPtr lower_bound;
  Expression::ExprPtr upper_bound;
  Expression::ExprPtr lower_dual_variable;
  Expression::ExprPtr upper_dual_variable;
};

struct EqualityConstraint {
  Expression::ExprPtr expr;
  Expression::ExprPtr rhs;
  Expression::ExprPtr dual_variable;
};

using VariableBounds = InequalityConstraint;

struct OptimizationProblem {
  Expression::ExprPtr objective;
  std::vector<InequalityConstraint> inequalities;
  std::vector<EqualityConstraint> equalities;
  std::vector<VariableBounds> variable_bounds;
  std::vector<Expression::ExprPtr> variables;
  std::vector<Expression::ExprPtr> variables2;
  std::vector<Expression::ExprPtr> variables3;
  std::vector<Expression::ExprPtr> variables4;
  std::vector<Expression::ExprPtr> nonnegative_slacks;
};

struct NewtonSystem {
  std::vector<std::vector<Expression::ExprPtr>> lhs;
  std::vector<Expression::ExprPtr> rhs;
  std::vector<Expression::ExprPtr> variables;
  std::vector<std::pair<Expression::ExprPtr, Expression::ExprPtr>>
      delta_definitions;
};

struct ShorthandRhs {
  std::vector<Expression::ExprPtr> shorthand_rhs;
  std::vector<std::pair<Expression::ExprPtr, Expression::ExprPtr>>
      vector_definitions;
};

struct OptimizationExpressions {
  Expression::ExprPtr Q;
  Expression::ExprPtr c;
  Expression::ExprPtr A_ineq;
  Expression::ExprPtr A_eq;
  Expression::ExprPtr b_eq;
  Expression::ExprPtr p_eq;
  Expression::ExprPtr delta_eq;
  Expression::ExprPtr mu;
  Expression::ExprPtr e_var;
  Expression::ExprPtr e_ineq;
  Expression::ExprPtr e_eq;
  Expression::ExprPtr x;
  Expression::ExprPtr s_A_ineq;
  Expression::ExprPtr s_A_ineq_l;
  Expression::ExprPtr s_A_ineq_u;
  Expression::ExprPtr s_x_l;
  Expression::ExprPtr s_x_u;
  Expression::ExprPtr s_A_eq;
  Expression::ExprPtr s_A_eq_l;
  Expression::ExprPtr s_A_eq_u;
  Expression::ExprPtr lambda_A_eq;
  Expression::ExprPtr lambda_sAeql;
  Expression::ExprPtr lambda_sAequ;
  Expression::ExprPtr lambda_A_ineq;
  Expression::ExprPtr lambda_sAineql;
  Expression::ExprPtr lambda_sAinequ;
  Expression::ExprPtr lambda_sxl;
  Expression::ExprPtr lambda_sxu;
  Expression::ExprPtr l_A_ineq;
  Expression::ExprPtr u_A_ineq;
  Expression::ExprPtr l_x;
  Expression::ExprPtr u_x;
};

OptimizationExpressions get_optimization_expressions(
    const VariableNames& names);

OptimizationProblem get_optimization_problem(
    const Settings& settings, const VariableNames& names,
    const OptimizationProblemType& optimization_problem_type);

Expression::ExprPtr get_lagrangian(const OptimizationProblem& problem);

std::pair<std::vector<Expression::ExprPtr>, std::vector<Expression::ExprPtr>>
get_first_order_optimality_conditions(Settings settings,
                                      const VariableNames& names);

NewtonSystem get_newton_system(const Settings& settings,
                               const VariableNames& names);
NewtonSystem get_augmented_system(NewtonSystem newton_system);
NewtonSystem get_normal_equations(NewtonSystem newton_system);

ShorthandRhs get_shorthand_rhs(const NewtonSystem& newton_system);

Expression::ExprPtr get_delta_variable(const Expression::ExprPtr& expr);

Expression::ExprPtr delta_definition(
    const std::vector<std::vector<Expression::ExprPtr>>& lhs,
    const std::vector<Expression::ExprPtr>& rhs,
    const std::vector<Expression::ExprPtr>& variables, size_t source_row);

void gaussian_elimination(std::vector<std::vector<Expression::ExprPtr>>& lhs,
                          std::vector<Expression::ExprPtr>& rhs,
                          size_t source_row);
}  // namespace SymbolicOptimization