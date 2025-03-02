#pragma once
#include "Expression.h"

namespace GaussianElimination {
void gaussianElimination(std::vector<std::vector<Expression::Expr>>& lhs,
                         std::vector<Expression::Expr>& rhs, size_t sourceRow);
}  // namespace GaussianElimination