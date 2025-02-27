#include "GaussianElimination.h"

#include <cassert>
#include <iostream>

namespace GaussianElimination {
void gaussianElimination(std::vector<std::vector<Expression::Expr>>& lhs,
                         std::vector<Expression::Expr>& rhs,
                         const size_t sourceRow, const bool print) {
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
          ExprFactory::product({ExprFactory::invert(sourceExpr), targetExpr}))
          .simplify();
  const auto addRowTimesFactorToRow = [&factor](const Expr& sourceTerm,
                                                const Expr& targetTerm) {
    const auto sourceTermTimesFactor =
        ExprFactory::product({factor, sourceTerm}).simplify();
    return ExprFactory::sum({targetTerm, sourceTermTimesFactor}).simplify();
  };
  if (print) {
    std::cout << "source: " << sourceRow << ", target: " << targetRow
              << std::endl;
    std::cout << "factor: " << factor.toString() << std::endl;
  }

  for (size_t i = 0; i < lhs.at(sourceRow).size(); ++i) {
    if (print) {
      std::cout << i << "add "
                << ExprFactory::product({factor, lhs.at(sourceRow).at(i)})
                       .simplify()
                       .toString()
                << " to row " << targetRow << " with result "
                << addRowTimesFactorToRow(lhs.at(sourceRow).at(i),
                                          lhs.at(targetRow).at(i))
                       .toString()
                << std::endl;
    }
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