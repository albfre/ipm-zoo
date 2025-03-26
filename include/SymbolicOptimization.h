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
  std::string s_A = "s";
  std::string s_Al = "g";
  std::string s_Au = "h";
  std::string s_xl = "y";
  std::string s_xu = "z";
  std::string s_C = "t";
  std::string s_Cl = "v";
  std::string s_Cu = "w";
  std::string l_A = "l_A";
  std::string u_A = "u_A";
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
  SimpleSlacks,
};

enum class EqualityHandling {
  None,
  Slacks,
  SimpleSlacks,
  PenaltyFunction,
  PenaltyFunctionWithExtraVariable,
  Regularization,
};

struct Settings {
  Bounds inequalities = Bounds::Both;
  Bounds variable_bounds = Bounds::Both;
  bool equalities = false;
  EqualityHandling equality_handling = EqualityHandling::None;
  InequalityHandling inequality_handling = InequalityHandling::Slacks;
};

struct NewtonSystem {
  std::vector<std::vector<Expression::ExprPtr>> lhs;
  std::vector<Expression::ExprPtr> rhs;
  std::vector<Expression::ExprPtr> variables;
  std::vector<std::pair<Expression::ExprPtr, Expression::ExprPtr>>
      delta_definitions;
};

std::pair<Expression::ExprPtr, std::vector<Expression::ExprPtr>> get_lagrangian(
    const Settings& settings, const VariableNames& names);

std::pair<std::vector<Expression::ExprPtr>, std::vector<Expression::ExprPtr>>
get_first_order_optimality_conditions(Settings settings,
                                      const VariableNames& names);

NewtonSystem get_newton_system(const Settings& settings,
                               const VariableNames& names);
NewtonSystem get_augmented_system(NewtonSystem newton_system);
NewtonSystem get_normal_equations(NewtonSystem newton_system);

std::vector<Expression::ExprPtr> get_shorthand_rhs(
    const std::vector<Expression::ExprPtr>& variables);

Expression::ExprPtr delta_definition(
    const std::vector<std::vector<Expression::ExprPtr>>& lhs,
    const std::vector<Expression::ExprPtr>& rhs,
    const std::vector<Expression::ExprPtr>& variables, size_t source_row);

void gaussian_elimination(std::vector<std::vector<Expression::ExprPtr>>& lhs,
                          std::vector<Expression::ExprPtr>& rhs,
                          size_t source_row);
}  // namespace SymbolicOptimization