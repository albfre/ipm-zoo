#pragma once
#include <numeric>

#include "Evaluation.h"
#include "Expr.h"
#include "NumericOptimization/LinearSolvers.h"
#include "SymbolicOptimization.h"
#include "Utils/Assert.h"
#include "Utils/Helpers.h"

namespace NumericOptimization {
using Matrix = std::vector<std::vector<double>>;
using Vector = std::vector<double>;
class NumericOptimization {
 public:
  NumericOptimization(Evaluation::Environment& env,
                      SymbolicOptimization::NewtonSystem newton_system)
      : env_(env),
        newton_system_(std::move(newton_system)),
        augmented_system_(
            SymbolicOptimization::get_augmented_system(newton_system_)),
        normal_equations_(
            SymbolicOptimization::get_normal_equations(augmented_system_)) {}

  void solve() {
    const auto& lhs = augmented_system_.lhs;
    const auto is_indefinite = std::ranges::any_of(
        std::views::iota(0ul, lhs.size()),
        [&lhs](const auto i) { return lhs.at(i).at(i) == Expression::zero; });
    if (is_indefinite) {
      solve_indefinite_();
    } else {
      solve_quasi_definite_();
    }
  }

 private:
  void solve_indefinite_() { ASSERT(false); }
  void solve_quasi_definite_() {
    const auto& [lhs, rhs, variables, delta_definitions] = augmented_system_;
    auto val_lhs = eval_expr_matrix_<Matrix>(lhs);
    auto val_rhs = eval_expr_vector_<Vector>(rhs);
    auto kkt = concatenate_matrices_(val_lhs);
    auto b = concatenate_vectors_(val_rhs);

    const auto print_mat = [](std::string name, const auto& M) {
      std::cout << name << std::endl;
      for (auto& row : M) {
        for (auto& col : row) {
          std::cout << col << ", ";
        }
        std::cout << std::endl;
      }
    };
    print_mat("kkt:", kkt);

    const auto [L, D] = LinearSolvers::ldlt_decomposition(kkt);
    print_mat("L:", L);

    LinearSolvers::overwriting_solve_ldlt(L, D, b);

    for (auto& bi : b) {
      std::cout << bi << ", ";
    }
    std::cout << std::endl;
  }

  template <typename T>
  std::vector<T> eval_expr_vector_(const std::vector<Expression::ExprPtr>& v) {
    std::vector<T> result;
    result.reserve(v.size());
    std::ranges::transform(v, std::back_inserter(result),
                           [this](const auto& expr) {
                             if constexpr (std::is_same_v<T, Matrix>) {
                               return Evaluation::evaluate_matrix(expr, env_);
                             } else {
                               return Evaluation::evaluate_vector(expr, env_);
                             }
                           });
    return result;
  }

  template <typename T>
  std::vector<std::vector<T>> eval_expr_matrix_(
      const std::vector<std::vector<Expression::ExprPtr>>& v) {
    std::vector<std::vector<T>> result;
    result.reserve(v.size());
    std::ranges::transform(
        v, std::back_inserter(result),
        [this](const auto& vec) { return eval_expr_vector_<T>(vec); });
    return result;
  }

  std::vector<double> concatenate_vectors_(
      const std::vector<std::vector<double>>& vectors) {
    const auto total_size = std::accumulate(
        vectors.begin(), vectors.end(), 0,
        [](size_t sum, const auto& vec) { return sum + vec.size(); });

    std::vector<double> result;
    result.reserve(total_size);

    for (const auto& vec : vectors) {
      result.insert(result.end(), vec.begin(), vec.end());
    }

    return result;
  }

  Matrix concatenate_matrices_(
      const std::vector<std::vector<Matrix>>& matrices) {
    if (matrices.empty()) {
      return Matrix();
    }
    {
      // Assert sizes
      const auto matrix_rows = matrices.size();
      const auto matrix_cols = matrices.at(0).size();
      ASSERT(std::ranges::all_of(matrices, [&](const auto& row) {
        return row.size() == matrix_cols;
      }));
      for (size_t i = 0; i < matrices.size(); ++i) {
        for (size_t j = 0; j < matrices.at(i).size(); ++j) {
          const auto& current = matrices.at(i).at(j);
          const auto& nextCol = matrices.at(i).at((j + 1) % matrix_cols);
          // All matrices in the same row has the same height
          ASSERT(current.size() == nextCol.size());
          const auto& nextRow = matrices.at((i + 1) % matrix_rows).at(j);
          // All matrices in the same col has the same width
          ASSERT((current.empty() ? 0 : current.at(0).size()) ==
                 (nextRow.empty() ? 0 : nextRow.at(0).size()));
        }
      }
    }

    const auto rows = std::accumulate(
        matrices.begin(), matrices.end(), 0,
        [](const auto s, const auto& m) { return s + m.size(); });
    const auto cols = std::accumulate(matrices.begin(), matrices.end(), 0,
                                      [](const auto s, const auto& m) {
                                        return s + m.empty() ? 0 : m[0].size();
                                      });
    Matrix result(rows, Vector(cols, 0.0));

    size_t row_offset = 0;
    for (const auto& row_of_matrices : matrices) {
      size_t col_offset = 0;
      for (const auto& mat : row_of_matrices) {
        for (size_t i = 0; i < mat.size(); ++i) {
          // All rows of the inner matrices have the same size
          ASSERT(mat[i].size() == mat[(i + 1) % mat.size()].size());
          std::ranges::copy(mat[i],
                            result[row_offset + i].begin() + col_offset);
        }

        col_offset += mat.empty() ? 0 : mat[0].size();
      }
      row_offset += row_of_matrices.empty() ? 0 : row_of_matrices.at(0).size();
    }

    return result;
  }

  Evaluation::Environment& env_;
  Expression::ExprPtr objective_;
  SymbolicOptimization::NewtonSystem newton_system_;
  SymbolicOptimization::NewtonSystem augmented_system_;
  SymbolicOptimization::NewtonSystem normal_equations_;
};

}  // namespace NumericOptimization