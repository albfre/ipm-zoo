#pragma once
#include <vector>

namespace NumericOptimization::LinearSolvers {
using Matrix = std::vector<std::vector<double>>;

/**
 * Performs LDL^T decomposition on a symmetric matrix A.
 * Returns L and D where A = LDL^T.
 */
std::pair<Matrix, std::vector<double>> ldlt_decomposition(const Matrix& A);

/**
 * Solves the system Ax = b using LDL^T decomposition.
 */
void overwriting_solve_ldlt(const Matrix& L, const std::vector<double>& D,
                            std::vector<double>& b);

/**
 * Performs symmetric indefinite factorization (Bunch-Kaufman).
 * Returns the factorized matrix and pivot information.
 */
std::pair<Matrix, std::vector<int>> symmetric_indefinite_factorization(
    const Matrix& matrix);

/**
 * Solves the system Ax = b using Bunch-Kaufman factorization.
 */
void overwriting_solve_bunch_kaufman(const Matrix& L,
                                     const std::vector<int>& ipiv,
                                     std::vector<double>& b);
}  // namespace NumericOptimization::LinearSolvers