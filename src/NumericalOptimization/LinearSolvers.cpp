#include "NumericalOptimization/LinearSolvers.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <tuple>
#include <vector>

#include "Utils/Assert.h"

namespace NumericalOptimization::LinearSolvers {
using Matrix = std::vector<std::vector<double>>;

std::pair<Matrix, std::vector<double>> ldlt_decomposition(const Matrix& A) {
  const auto n = A.size();
  ASSERT(
      std::ranges::all_of(A, [n](const auto& Ai) { return Ai.size() == n; }));
  Matrix L(n, std::vector<double>(n, 0.0));
  std::vector<double> D(n, 0.0);

  for (size_t i = 0; i < n; ++i) {
    auto sum_d = A[i][i];
    for (size_t j = 0; j < i; ++j) {
      sum_d -= L[i][j] * L[i][j] * D[j];
    }
    // Modify factorization to avoid zero diagonal (Vanderbei, Symmetric
    // quasi-definite matrices, 1995)
    D[i] = sum_d == 0.0 ? 1e-8 : sum_d;

    for (size_t j = i + 1; j < n; ++j) {
      auto sum = A[j][i];
      for (size_t k = 0; k < i; ++k) {
        sum -= L[j][k] * L[i][k] * D[k];
      }
      L[j][i] = sum / D[i];
    }

    L[i][i] = 1.0;
  }

  return {L, D};
}

void overwriting_solve_ldlt(const Matrix& L, const std::vector<double>& D,
                            std::vector<double>& b) {
  if (b.empty()) {
    return;
  }
  const auto n = b.size();
  ASSERT(D.size() == n);
  ASSERT(L.size() == n);
  ASSERT(
      std::ranges::all_of(L, [n](const auto& Li) { return Li.size() == n; }));

  // Forward substitution (Ly = b)
  for (size_t i = 0; i < n; ++i) {
    b[i] -= std::inner_product(L[i].begin(), L[i].begin() + i, b.begin(),
                               0.0);  // No division by L[i][i] since it is 1.0
  }

  // Diagonal scaling (Dz = y)
  for (size_t i = 0; i < n; ++i) {
    b[i] /= D[i];
  }

  // Backward substitution (L^T x = z)
  for (int i = static_cast<int>(n) - 1; i >= 0; --i) {
    auto sum = 0.0;
    for (size_t j = i + 1; j < n; ++j) {
      sum += L[j][i] * b[j];
    }
    b[i] -= sum;
  }
}

std::pair<Matrix, std::vector<int>> symmetric_indefinite_factorization(
    const Matrix& matrix) {
  auto A = matrix;
  const auto n = A.size();
  ASSERT(
      std::ranges::all_of(A, [n](const auto& Ai) { return Ai.size() == n; }));
  const auto alpha = (1.0 + std::sqrt(17.0)) / 8.0;
  std::vector<int> ipiv(n, 0);

  const auto max_in_row_or_column = [&](const size_t i_begin,
                                        const size_t i_end,
                                        const size_t constant_index,
                                        const bool check_column) {
    size_t i_max = 0;
    auto col_max = 0.0;
    for (size_t i = i_begin; i < i_end; ++i) {
      const auto v =
          std::abs(check_column ? A[i][constant_index] : A[constant_index][i]);
      if (v > col_max) {
        col_max = v;
        i_max = i;
      }
    }
    return std::pair{i_max, col_max};
  };

  int info = 0;
  size_t k = 0;
  while (k < n) {
    int k_step = 1;
    size_t kp = 0;
    const auto absakk = std::abs(A[k][k]);
    // imax is the row-index of the largest off-diagonal element in column k,
    // and colmax is its absolute value
    const auto [i_max, col_max] = max_in_row_or_column(k + 1, n, k, true);
    if (absakk == 0.0 && col_max == 0.0) {
      // Column k is zero: set info and continue
      if (info == 0) {
        info = static_cast<int>(k);
        kp = k;
      }
    } else {
      if (absakk >= alpha * col_max) {
        // No interchange, use 1-by-1 pivot block
        kp = k;
      } else {
        const auto [_, row_max1] = max_in_row_or_column(k, i_max, i_max, false);
        const auto [__, row_max2] =
            max_in_row_or_column(i_max + 1, n, i_max, true);
        const auto row_max = std::max(row_max1, row_max2);
        if (absakk * row_max >= alpha * col_max * col_max) {
          kp = k;  // No interchange, use 1-by-1 pivot block
        } else if (std::abs(A[i_max][i_max]) >= alpha * row_max) {
          kp = i_max;  // Interchange rows and columns k and imax, use 1-by-1
                       // pivot block
        } else {
          // Interchange rows and columns k+1 and imax, use 2-by-2 pivot block
          kp = i_max;
          k_step = 2;
        }
      }

      const auto kk = k + k_step - 1;
      if (kp != kk) {
        // Interchange rows and columns kk and kp in the trailing submatrix
        // A(k:n,k:n)
        for (size_t i = kp + 1; i < n; ++i) {
          std::swap(A[i][kp], A[i][kk]);
        }
        for (size_t j = kk + 1; j < kp; ++j) {
          std::swap(A[kp][j], A[j][kk]);
        }
        std::swap(A[kp][kp], A[kk][kk]);
        if (k_step == 2) {
          std::swap(A[kk][k], A[kp][k]);
        }
      }

      // Update the trailing submatrix
      if (k_step == 1) {
        // 1-by-1 pivot block D(k): column k now holds W(k) = L(k)*D(k) where
        // L(k) is the k-th column of L Perform a rank-1 update of
        // A(k+1:n,k+1:n) as A := A - L(k)*D(k)*L(k)**T = A -
        // W(k)*(1/D(k))*W(k)**T
        const auto inv_pivot = 1.0 / A[k][k];
        for (size_t j = k + 1; j < n; ++j) {
          const auto scaled_factor = inv_pivot * A[j][k];
          for (size_t i = j; i < n; ++i) {
            A[i][j] -= scaled_factor * A[i][k];
          }
          A[j][k] *= inv_pivot;
        }
      } else {
        // 2-by-2 pivot block D(k): columns k and k+1 now hold ( W(k) W(k+1) )
        // = ( L(k) L(k+1) )*D(k) where L(k) and L(k+1) are the k-th and
        // (k+1)-th columns of L
        if (k < n - 1) {
          // Perform a rank-2 update of A(k+2:n,k+2:n) as
          // A := A - ( L(k) L(k+1) )*D(k)*( L(k) L(k+1) )**T = A - ( W(k)
          // W(k+1) )*inv(D(k))*( W(k) W(k+1) )**T where L(k) and L(k+1) are
          // the k-th and (k+1)-th columns of L
          auto d21 = A[k + 1][k];
          const auto d11 = A[k + 1][k + 1] / d21;
          const auto d22 = A[k][k] / d21;
          const auto t = 1.0 / (d11 * d22 - 1.0);
          d21 = t / d21;

          for (size_t j = k + 2; j < n; ++j) {
            const auto wk = d21 * (d11 * A[j][k] - A[j][k + 1]);
            const auto wkp1 = d21 * (d22 * A[j][k + 1] - A[j][k]);
            for (size_t i = j; i < n; ++i) {
              A[i][j] -= (A[i][k] * wk + A[i][k + 1] * wkp1);
            }
            A[j][k] = wk;
            A[j][k + 1] = wkp1;
          }
        }
      }
    }

    if (k_step == 1) {
      ipiv[k] = static_cast<int>(kp);
    } else {
      ipiv[k] = -static_cast<int>(kp);
      ipiv[k + 1] = -static_cast<int>(kp);
    }

    k += k_step;
  }

  return {std::move(A), std::move(ipiv)};
}

void overwriting_solve_indefinite(const Matrix& L, const std::vector<int>& ipiv,
                                  std::vector<double>& b) {
  if (b.empty()) {
    return;
  }
  const auto n = b.size();
  ASSERT(ipiv.size() == n);
  ASSERT(L.size() == n);
  ASSERT(
      std::ranges::all_of(L, [n](const auto& Li) { return Li.size() == n; }));

  auto apply_row_transformation = [&](size_t i_start, size_t j_index) {
    auto multiplier = -b[j_index];
    for (size_t i = i_start; i < n; ++i) {
      b[i] += L[i][j_index] * multiplier;
    }
  };

  size_t k = 0;
  while (k < n) {
    if (ipiv[k] >= 0) {
      // 1 x 1 diagonal block, interchange rows k and ipiv(k).
      const auto kp = static_cast<size_t>(ipiv[k]);
      if (kp != k) {
        std::swap(b[k], b[kp]);
      }
      // Multiply by inv(L(k)), where L(k) is the transformation stored in
      // column k of L.
      apply_row_transformation(k + 1, k);
      b[k] /= L[k][k];
      k += 1;
    } else {
      // 2 x 2 diagonal block, interchange rows k+1 and -ipiv(k).
      const auto kp = static_cast<size_t>(-ipiv[k]);
      if (kp != k + 1) {
        std::swap(b[k + 1], b[kp]);
      }

      // Multiply by inv(L(k)), where L(k) is the transformation stored in
      // columns k and k+1 of L.
      if (k < n - 1) {
        apply_row_transformation(k + 2, k);
        apply_row_transformation(k + 2, k + 1);
      }

      // Multiply by the inverse of the diagonal block.
      const auto akm1k = L[k + 1][k];
      const auto akm1 = L[k][k] / akm1k;
      const auto ak = L[k + 1][k + 1] / akm1k;
      const auto denom = akm1 * ak - 1.0;
      const auto bkm1 = b[k] / akm1k;
      const auto bk = b[k + 1] / akm1k;
      b[k] = (ak * bkm1 - bk) / denom;
      b[k + 1] = (akm1 * bk - bkm1) / denom;
      k += 2;
    }
  }

  // Next solve L**T *X = B, overwriting B with X.
  // k is the main loop index, decreasing from n-1 to 0 in steps of 1 or 2,
  // depending on the size of the diagonal blocks.
  auto compute_transpose_product = [&](size_t i_start, size_t j_index) {
    auto sum = 0.0;
    for (size_t i = i_start; i < n; ++i) {
      sum += L[i][j_index] * b[i];
    }
    b[j_index] -= sum;
  };

  k = n - 1;
  while (k < n) {  // Using underflow for unsigned comparison
    if (ipiv[k] >= 0) {
      // 1 x 1 diagonal block, multiply by inv(L**T(k)), where L(k) is the
      // transformation stored in column k of L.
      if (k < n - 1) {
        compute_transpose_product(k + 1, k);
      }

      // Interchange rows k and ipiv(k).
      const auto kp = static_cast<size_t>(ipiv[k]);
      if (kp != k) {
        std::swap(b[k], b[kp]);
      }

      if (k == 0) {
        break;
      }
      --k;
    } else {
      // 2 x 2 diagonal block, multiply by inv(L**T(k-1)), where L(k-1) is the
      // transformation stored in columns k-1 and k of L.
      if (k < n - 1) {
        compute_transpose_product(k + 1, k);
        compute_transpose_product(k + 1, k - 1);
      }

      // Interchange rows k and -ipiv(k).
      const auto kp = static_cast<size_t>(-ipiv[k]);
      if (kp != k) {
        std::swap(b[k], b[kp]);
      }

      if (k <= 1) {
        break;
      }
      k -= 2;
    }
  }
}

}  // namespace NumericalOptimization::LinearSolvers