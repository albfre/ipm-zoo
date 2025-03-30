#include "NumericalOptimization/Optimizer.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>

#include "ExprFactory.h"
#include "NumericalOptimization/EnvironmentBuilder.h"
#include "NumericalOptimization/LinearSolvers.h"
#include "Utils/Assert.h"
#include "Utils/Helpers.h"

namespace NumericalOptimization {

Optimizer::Optimizer(Evaluation::Environment& env,
                     const SymbolicOptimization::OptimizationExpressions&
                         optimization_expressions,
                     SymbolicOptimization::NewtonSystem newton_system)
    : env_(env),
      optimization_expressions_(optimization_expressions),
      newton_system_(std::move(newton_system)),
      augmented_system_(
          SymbolicOptimization::get_augmented_system(newton_system_)),
      normal_equations_(
          SymbolicOptimization::get_normal_equations(augmented_system_)) {
  using namespace Expression;
  const auto Q = optimization_expressions.Q;
  const auto x = optimization_expressions.x;
  const auto cT = ExprFactory::transpose(optimization_expressions.c);
  objective_ =
      ExprFactory::sum({ExprFactory::product({ExprFactory::number(0.5),
                                              ExprFactory::transpose(x), Q, x}),
                        ExprFactory::product({cT, x})});
  const auto var_vec = eval_expr_vector_<Vector>(newton_system_.variables);
  for (size_t i = 0; i < var_vec.size(); ++i) {
    const auto& variable = newton_system_.variables.at(i);
    variable_index_[variable] = i;
    vector_sizes_[variable] = var_vec.at(i).size();
    const auto delta_variable =
        SymbolicOptimization::get_delta_variable(variable);
    delta_variable_index_[delta_variable] = i;
  }
}

void Optimizer::solve() {
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

void Optimizer::solve_indefinite_() { ASSERT(false); }

void Optimizer::solve_quasi_definite_() {
  const auto print_mat = [](std::string name, const auto& M) {
    std::cout << name << std::endl;
    for (auto& row : M) {
      if (!row.empty()) {
        for (auto& col : row) {
          std::cout << col << ", ";
        }
        std::cout << std::endl;
      }
    }
    std::cout << std::endl;
  };
  const auto print_var = [&](std::string name,
                             const std::vector<Expression::ExprPtr>& M) {
    std::cout << name << std::endl;
    for (auto& row : M) {
      auto row_val = Evaluation::evaluate_vector(row, env_);
      if (!row_val.empty()) {
        std::cout << row->to_string() << ": ";
        for (auto& col : row_val) {
          std::cout << col << ", ";
        }
        std::cout << std::endl;
      }
    }
    std::cout << std::endl;
  };
  const auto& [lhs, rhs, variables, delta_definitions] = augmented_system_;
  const auto& full_rhs = newton_system_.rhs;

  const auto sigma = 1.0;
  const auto tolerance = 1e-8;
  const size_t max_iter = 30;
  size_t iter = 0;
  for (; iter < max_iter; ++iter) {
    const auto f = Evaluation::evaluate_scalar(objective_, env_);
    const auto residual_norm = get_residual_norm_(full_rhs);
    const auto mu = get_mu_(full_rhs);
    // env_.at(optimization_expressions_.mu) = Evaluation::val_scalar(mu);
    std::cout << "iter: " << iter << ", f: " << f << ", res: " << residual_norm
              << ", gap: " << mu << std::endl;
    if (residual_norm < tolerance && mu < tolerance) {
      break;
    }

    const auto kkt = get_as_matrix_(lhs);
    print_mat("kkt:", kkt);
    const auto [L, D] = LinearSolvers::symmetric_indefinite_factorization(kkt);

    auto variables = eval_expr_vector_<Vector>(newton_system_.variables);
    const auto delta =
        compute_search_direction_(augmented_system_, L, D, sigma);
    print_var("var", newton_system_.variables);
    print_mat("delta", delta);
    const auto alpha = get_max_step_(variables, delta);
    std::cout << " alpha: " << alpha << std::endl;
    const auto fraction_to_boundary = 0.995;
    update_variables_(fraction_to_boundary * alpha, variables, delta);
  }
}

void Optimizer::update_variables_(double alpha, Matrix& variables,
                                  const Matrix& delta) {
  ASSERT(variables.size() == delta.size());
  ASSERT(variables.size() == newton_system_.variables.size());
  for (size_t i = 0; i < variables.size(); ++i) {
    vector_plus_eq_scalar_times_vector_(variables[i], alpha, delta[i]);
    env_.at(newton_system_.variables.at(i)) =
        Evaluation::val_vector(variables[i]);
  }
}

void Optimizer::vector_plus_eq_scalar_times_vector_(
    std::vector<double>& x, const double s, const std::vector<double>& y) {
  ASSERT(x.size() == y.size());
  std::ranges::transform(x, y, x.begin(),
                         [s](double xi, double yi) { return xi + s * yi; });
}

double Optimizer::get_residual_norm_(
    const std::vector<Expression::ExprPtr>& rhs) {
  const auto& mu = optimization_expressions_.mu;
  ScopedEnvironmentOverride temp(env_, mu, Evaluation::val_scalar(0.0));
  for (auto& r : rhs) {
    std::cout << r->to_string() << std::endl;
    auto rv = Evaluation::evaluate_vector(r, env_);
    for (auto& rr : rv) {
      std::cout << rr << ", ";
    }
    std::cout << std::endl;
  }

  const auto rhs_vec = get_as_vector_(rhs);
  const auto residual = std::sqrt(dot_(rhs_vec, rhs_vec));
  return residual;
}

double Optimizer::get_mu_(const std::vector<Expression::ExprPtr>& rhs) {
  const auto& mu = optimization_expressions_.mu;
  const auto& e_var = optimization_expressions_.e_var;
  const auto& e_ineq = optimization_expressions_.e_ineq;
  const auto& e_eq = optimization_expressions_.e_eq;
  auto filtered = rhs | std::views::filter([&](const auto& expr) {
                    return (expr->contains_subexpression(e_var) ||
                            expr->contains_subexpression(e_ineq) ||
                            expr->contains_subexpression(e_eq)) &&
                           expr->contains_subexpression(mu);
                  });
  const auto mu_terms =
      std::vector<Expression::ExprPtr>{filtered.begin(), filtered.end()};
  ScopedEnvironmentOverride temp(env_, mu, Evaluation::val_scalar(0.0));
  const auto mu_vec = get_as_vector_(mu_terms);
  const auto sum =
      std::accumulate(mu_vec.begin(), mu_vec.end(), 0.0, std::plus<>());
  return mu_vec.empty() ? 0.0 : -sum / mu_vec.size();
}

double Optimizer::get_max_step_(Matrix variables, Matrix delta) {
  ASSERT(variables.size() == delta.size());
  double max_step = 1.0;
  const auto& oe = optimization_expressions_;
  static const auto non_negative = std::set{
      oe.s_A_ineq_l,     oe.s_A_ineq_u,     oe.s_x_l,        oe.s_x_u,
      oe.s_A_eq_l,       oe.s_A_eq_u,       oe.lambda_sAeql, oe.lambda_sAequ,
      oe.lambda_sAineql, oe.lambda_sAinequ, oe.lambda_sxl,   oe.lambda_sxu,
  };
  ASSERT(variables.size() == newton_system_.variables.size());
  for (size_t i = 0; i < variables.size(); ++i) {
    if (!non_negative.contains(newton_system_.variables.at(i))) {
      continue;
    }
    ASSERT(variables.at(i).size() == delta.at(i).size());
    for (size_t j = 0; j < variables.at(i).size(); ++j) {
      const auto dij = delta.at(i).at(j);
      const auto vij = variables.at(i).at(j);
      if (dij < 0.0) {
        max_step = std::min(max_step, -vij / dij);
      }
    }
  }
  return max_step;
}

std::vector<std::vector<double>> Optimizer::compute_search_direction_(
    const SymbolicOptimization::NewtonSystem& newton_system, const Matrix& L,
    const std::vector<int>& D, const double sigma) {
  const auto& [lhs, rhs, variables, delta_definitions] = newton_system;
  const auto& mu = optimization_expressions_.mu;
  ScopedEnvironmentOverride temp(env_, mu,
                                 Evaluation::scale(env_.at(mu), sigma));
  auto b = get_as_vector_(rhs);
  LinearSolvers::overwriting_solve_bunch_kaufman(L, D, b);
  auto result = std::vector<std::vector<double>>(variable_index_.size());
  size_t offset = 0;
  for (const auto& var : newton_system.variables) {
    const auto index = variable_index_.at(var);
    const auto next_offset = offset + vector_sizes_.at(var);
    ASSERT(b.size() >= next_offset);
    result.at(index) = std::vector(b.begin() + offset, b.begin() + next_offset);
    const auto delta_variable = SymbolicOptimization::get_delta_variable(var);
    env_[delta_variable] = Evaluation::val_vector(result.at(index));
    offset = next_offset;
  }
  for (const auto& [delta_variable, definition] :
       std::views::reverse(newton_system.delta_definitions)) {
    const auto index = delta_variable_index_.at(delta_variable);
    result.at(index) = Evaluation::evaluate_vector(definition, env_);
    env_[delta_variable] = Evaluation::val_vector(result.at(index));
  }
  return result;
}

double Optimizer::dot_(const std::vector<double>& x,
                       const std::vector<double>& y) {
  return std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
}

Matrix Optimizer::get_as_matrix_(
    const std::vector<std::vector<Expression::ExprPtr>>& v) {
  const auto val = eval_expr_matrix_<Matrix>(v);
  return concatenate_matrices_(val);
}

Vector Optimizer::get_as_vector_(const std::vector<Expression::ExprPtr>& v) {
  const auto val = eval_expr_vector_<Vector>(v);
  return concatenate_vectors_(val);
}

template <typename T>
std::vector<T> Optimizer::eval_expr_vector_(
    const std::vector<Expression::ExprPtr>& v) {
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
std::vector<std::vector<T>> Optimizer::eval_expr_matrix_(
    const std::vector<std::vector<Expression::ExprPtr>>& v) {
  std::vector<std::vector<T>> result;
  result.reserve(v.size());
  std::ranges::transform(
      v, std::back_inserter(result),
      [this](const auto& vec) { return eval_expr_vector_<T>(vec); });
  return result;
}

std::vector<double> Optimizer::concatenate_vectors_(
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

Matrix Optimizer::concatenate_matrices_(
    const std::vector<std::vector<Matrix>>& matrices) {
  if (matrices.empty()) {
    return Matrix();
  }
  {
    // Assert sizes
    const auto matrix_rows = matrices.size();
    const auto matrix_cols = matrices.at(0).size();
    ASSERT(std::ranges::all_of(
        matrices, [&](const auto& row) { return row.size() == matrix_cols; }));
    for (size_t i = 0; i < matrices.size(); ++i) {
      for (size_t j = 0; j < matrices.at(i).size(); ++j) {
        const auto& current = matrices.at(i).at(j);
        const auto& next_col = matrices.at(i).at((j + 1) % matrix_cols);
        // All matrices in the same row has the same height
        /*
                const auto& next_row = matrices.at((i + 1) % matrix_rows).at(j);
                // All matrices in the same col has the same width
                ASSERT((current.empty() ? 0 : current.at(0).size()) ==
                       (next_row.empty() ? 0 : next_row.at(0).size()));
                       */
      }
    }
  }

  // Compute total number of rows and columns
  std::vector<size_t> max_height_per_row(matrices.size());
  std::vector<size_t> max_width_per_col(matrices.at(0).size());
  for (size_t i = 0; i < matrices.size(); ++i) {
    for (size_t j = 0; j < matrices.at(0).size(); ++j) {
      const auto& m = matrices.at(i).at(j);
      const auto height = m.size();
      const auto width = m.empty() ? 0 : m.at(0).size();
      max_height_per_row.at(i) = std::max(max_height_per_row.at(i), height);
      max_width_per_col.at(j) = std::max(max_width_per_col.at(j), width);
    }
  }
  const auto total_rows = std::accumulate(
      max_height_per_row.begin(), max_height_per_row.end(), 0, std::plus<>());
  const auto total_cols = std::accumulate(
      max_width_per_col.begin(), max_width_per_col.end(), 0, std::plus<>());

  Matrix result(total_rows, Vector(total_cols, 0.0));

  size_t row_offset = 0;
  for (size_t i = 0; i < matrices.size(); ++i) {
    size_t col_offset = 0;
    for (size_t j = 0; j < matrices.at(0).size(); ++j) {
      const auto& mat = matrices.at(i).at(j);
      for (size_t k = 0; k < mat.size(); ++k) {
        std::ranges::copy(mat[k], result[row_offset + k].begin() + col_offset);
      }
      col_offset += max_width_per_col.at(j);
    }
    row_offset += max_height_per_row.at(i);
  }

  return result;
}

}  // namespace NumericalOptimization
