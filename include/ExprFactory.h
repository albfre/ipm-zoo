#pragma once
#include "Expr.h"

namespace Expression::ExprFactory {
[[nodiscard]] Expr number(const double value);
[[nodiscard]] Expr namedScalar(std::string_view name);
[[nodiscard]] Expr namedVector(std::string_view name);
[[nodiscard]] Expr variable(std::string_view name);
[[nodiscard]] Expr matrix(std::string_view name);
[[nodiscard]] Expr symmetricMatrix(std::string_view name);
[[nodiscard]] Expr diagonalMatrix(Expr expr);
[[nodiscard]] Expr transpose(Expr expr);
[[nodiscard]] Expr negate(Expr expr);
[[nodiscard]] Expr invert(Expr expr);
[[nodiscard]] Expr log(Expr expr);
[[nodiscard]] Expr sum(std::vector<Expr> terms);
[[nodiscard]] Expr product(std::vector<Expr> terms);
}  // namespace Expression::ExprFactory