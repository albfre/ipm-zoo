#include "GaussianElimination.h"

#include <cassert>
#include <iostream>

namespace GaussianElimination {
void gaussianElimination(std::vector<std::vector<Expression::Expr>>& lhs,
                         std::vector<Expression::Expr>& rhs,
                         const size_t sourceRow) {
  using namespace Expression;
  const auto zero = ExprFactory::number(0.0);
  size_t targetRow = 0;
  for (; targetRow < lhs.size(); ++targetRow) {
    if (targetRow != sourceRow && lhs.at(targetRow).at(sourceRow) != zero) {
      break;
    }
  }
  assert(targetRow < lhs.size());
  const auto targetExpr = lhs.at(targetRow).at(sourceRow);
  const auto sourceExpr = lhs.at(sourceRow).at(sourceRow);
  const auto factor =
      ExprFactory::negate(
          ExprFactory::product({targetExpr, ExprFactory::invert(sourceExpr)}))
          .simplify();
  const auto addRowTimesFactorToRow = [&factor](const Expr& sourceTerm,
                                                const Expr& targetTerm) {
    const auto sourceTermTimesFactor =
        ExprFactory::product({factor, sourceTerm}).simplify();
    return ExprFactory::sum({targetTerm, sourceTermTimesFactor}).simplify();
  };

  for (size_t i = 0; i < lhs.at(sourceRow).size(); ++i) {
    lhs.at(targetRow).at(i) = addRowTimesFactorToRow(lhs.at(sourceRow).at(i),
                                                     lhs.at(targetRow).at(i));
  }

  rhs.at(targetRow) =
      addRowTimesFactorToRow(rhs.at(sourceRow), rhs.at(targetRow));

  lhs.erase(lhs.begin() + sourceRow);
  for (auto& lhsRow : lhs) {
    lhsRow.erase(lhsRow.begin() + sourceRow);
  }
  rhs.erase(rhs.begin() + sourceRow);
}
}  // namespace GaussianElimination