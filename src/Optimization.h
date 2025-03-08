#pragma once
#include "Expression.h"

namespace Optimization {
struct VariableNames {
  std::string x = "x";
  std::string A_eq = "C";
  std::string b_eq = "b";
  std::string A_ineq = "A";
  std::string s_A = "s";
  std::string s_Al = "g";
  std::string s_Au = "t";
  std::string s_xl = "y";
  std::string s_xu = "z";
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
  IndefiniteFactorization,
  Slacks,
  PenaltyFunction,
};

struct Settings {
  Bounds inequalities = Bounds::Both;
  Bounds variableBounds = Bounds::Both;
  bool equalities = false;
  EqualityHandling equalityHandling = EqualityHandling::IndefiniteFactorization;
  InequalityHandling inequalityHandling = InequalityHandling::Slacks;
};

std::pair<Expression::Expr, std::vector<Expression::Expr>> getLagrangian(
    const VariableNames& names, const Settings& settings);

std::vector<Expression::Expr> getFirstOrderOptimalityConditions(
    const Expression::Expr& lagrangian,
    const std::vector<Expression::Expr>& variables);

std::pair<std::vector<std::vector<Expression::Expr>>,
          std::vector<Expression::Expr>>
getNewtonSystem(const Expression::Expr& lagrangian,
                const std::vector<Expression::Expr>& variables);

std::vector<Expression::Expr> getShorthandRhs(
    const std::vector<Expression::Expr>& variables);

std::pair<Expression::Expr, Expression::Expr> gaussianElimination(
    std::vector<std::vector<Expression::Expr>>& lhs,
    std::vector<Expression::Expr>& rhs, size_t sourceRow,
    const std::vector<Expression::Expr>& variables);
}  // namespace Optimization