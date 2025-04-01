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
      newton_system_(newton_system),
      shorthand_rhs_(SymbolicOptimization::get_shorthand_rhs(newton_system_)) {
  using namespace Expression;
  newton_system.rhs = shorthand_rhs_.shorthand_rhs;
  augmented_system_ = SymbolicOptimization::get_augmented_system(newton_system);
  normal_equations_ =
      SymbolicOptimization::get_normal_equations(augmented_system_);

  const auto Q = optimization_expressions.Q;
  const auto x = optimization_expressions.x;
  const auto cT = ExprFactory::transpose(optimization_expressions.c);
  objective_ =
      ExprFactory::sum({ExprFactory::product({ExprFactory::number(0.5),
                                              ExprFactory::transpose(x), Q, x}),
                        ExprFactory::product({cT, x})});
  const auto var_vec =
      eval_vector_of_expressions_<Vector>(newton_system_.variables);
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
  const auto print_mat = [&](std::string name, const auto& M,
                             bool print_names = true) {
    std::cout << name << std::endl;
    size_t i = 0;
    for (auto& row : M) {
      if (!row.empty()) {
        if (print_names) {
          std::cout << newton_system_.variables.at(i)->to_string() << ": ";
        }
        for (auto& col : row) {
          std::cout << col << ", ";
        }
        std::cout << std::endl;
      }
      ++i;
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
  const auto& full_rhs = newton_system_.rhs;
  const auto& augmented_lhs = augmented_system_.lhs;

  std::vector<Expression::ExprPtr> delta_aff_variables;
  std::ranges::transform(
      newton_system_.variables, std::back_inserter(delta_aff_variables),
      [](const auto& v) {
        const auto dv = SymbolicOptimization::get_delta_variable(v);
        ASSERT(Expression::is<Expression::Variable>(dv));
        const auto& var = std::get<Expression::Variable>(dv->get_impl());
        return Expression::ExprFactory::variable(var.name + "_affine");
      });

  const auto tolerance = 1e-8;
  const size_t max_iter = 20;
  size_t iter = 0;
  for (; iter < max_iter; ++iter) {
    const auto f = Evaluation::evaluate_scalar(objective_, env_);
    const auto residual_norm = get_residual_norm_(full_rhs);
    const auto mu_val = get_mu_(full_rhs);
    std::cout << "iter: " << iter << ", f: " << f << ", res: " << residual_norm
              << ", gap: " << mu_val << std::endl;
    if (residual_norm < tolerance && mu_val < tolerance) {
      break;
    }

    const auto kkt = get_as_matrix_(augmented_lhs);
    // print_mat("kkt:", kkt, false);
    const auto [L, D] = LinearSolvers::ldlt_decomposition(kkt);
    // print_mat("L:", L, false);
    /*
    std::cout << "D: ";
    for (auto& di : D) {
      std::cout << di << ", ";
    }
    std::cout << std::endl;
    */

    // Compute affine scaling step
    env_.at(optimization_expressions_.mu) = Evaluation::val_scalar(0.0);
    for (const auto& [vec, def] : shorthand_rhs_.vector_definitions) {
      env_[vec] = Evaluation::evaluate(def, env_);
    }

    /*
    auto b = get_as_vector_(augmented_system_.rhs);
    std::cout << "b: ";
    for (auto& bi : b) {
      std::cout << bi << ", ";
    }
    std::cout << std::endl;
    */

    const auto variable_values =
        eval_vector_of_expressions_<Vector>(newton_system_.variables);
    const auto delta_aff_values =
        compute_search_direction_(augmented_system_, L, D);
    print_var("var", newton_system_.variables);
    print_mat("delta affine", delta_aff_values);
    {
      // Temporarily update the variables to compute mu_affine
      const auto alpha_affine =
          get_max_step_(variable_values, delta_aff_values);
      std::cout << " alpha affine: " << alpha_affine << std::endl;
      std::vector<std::unique_ptr<ScopedEnvironmentOverride>> temps;
      temps.reserve(newton_system_.variables.size());
      for (const auto& var : newton_system_.variables) {
        auto val = Evaluation::evaluate(var, env_);
        temps.push_back(
            std::make_unique<ScopedEnvironmentOverride>(env_, var, val));
      }
      update_variables_(alpha_affine, variable_values, delta_aff_values);
      const auto mu_affine = get_mu_(full_rhs);
      const auto sigma = mu_val > 0.0 ? std::pow(mu_affine / mu_val, 3) : 0.0;
      env_.at(optimization_expressions_.mu) =
          Evaluation::val_scalar(mu_val * sigma);
      std::cout << "sigma " << sigma << ", new mu: " << mu_val * sigma
                << std::endl;
    }
    {
      // Update complementarity rhs to include the affine scaling direction
      const auto& mu = optimization_expressions_.mu;
      const auto& e_var = optimization_expressions_.e_var;
      const auto& e_ineq = optimization_expressions_.e_ineq;
      const auto& e_eq = optimization_expressions_.e_eq;
      std::cout << "residuals2" << std::endl;
      for (const auto& [shorthand, expr] : shorthand_rhs_.vector_definitions) {
        auto shorthand_val = Evaluation::evaluate(expr, env_);
        if ((expr->contains_subexpression(e_var) ||
             expr->contains_subexpression(e_ineq) ||
             expr->contains_subexpression(e_eq)) &&
            expr->contains_subexpression(mu)) {
          auto delta_expr = expr->replace_subexpression(mu, Expression::zero);
          const auto expr_variables = expr->get_variables();
          for (auto& var : expr_variables) {
            const auto index = variable_index_.at(var);
            const auto& delta_aff_variable = delta_aff_variables.at(index);
            const auto& delta_aff_value = delta_aff_values.at(index);
            env_[delta_aff_variable] = Evaluation::val_vector(delta_aff_value);
            delta_expr =
                delta_expr->replace_subexpression(var, delta_aff_variable);
          }

          const auto delta_expr_val = Evaluation::evaluate(delta_expr, env_);
          shorthand_val = Evaluation::add(shorthand_val, delta_expr_val);
        }
        env_[shorthand] = shorthand_val;
        std::cout << shorthand->to_string() << ": ";
        const auto vv = std::get<Evaluation::ValVector>(shorthand_val);
        for (const auto& v : vv) {
          std::cout << v << ", ";
        }
        std::cout << std::endl;
      }

      // Compute aggregated centering-corrector step
      const auto delta = compute_search_direction_(augmented_system_, L, D);
      print_mat("delta ", delta);
      const auto alpha = get_max_step_(variable_values, delta);
      std::cout << " alpha: " << alpha << std::endl;

      constexpr auto fraction_to_boundary = 0.995;
      update_variables_(fraction_to_boundary * alpha, variable_values, delta);
    }
  }
}

void Optimizer::update_variables_(double alpha, Matrix variables,
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
  std::cout << "residuals" << std::endl;
  size_t i = 0;
  for (auto& r : rhs) {
    std::cout << newton_system_.variables.at(i++)->to_string() << ": ";
    std::cout << "(" << r->to_string() << "):  ";
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
  const auto sum = std::accumulate(
      mu_vec.begin(), mu_vec.end(), 0.0,
      [](const auto s, const auto& v) { return s + std::abs(v); });
  return mu_vec.empty() ? 0.0 : sum / mu_vec.size();
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
    const std::vector<double>& D) {
  const auto& [lhs, rhs, variables, delta_definitions] = newton_system;
  auto b = get_as_vector_(rhs);

  LinearSolvers::overwriting_solve_ldlt(L, D, b);
  std::cout << "sol: ";
  for (auto& bi : b) {
    std::cout << bi << ", ";
  }
  std::cout << std::endl;

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
  for (auto& [rh, _] : shorthand_rhs_.vector_definitions) {
    std::cout << rh->to_string() << ": ";
    auto v = Evaluation::evaluate_vector(rh, env_);
    for (auto& vv : v) {
      std::cout << vv << ", ";
    }
    std::cout << std::endl;
  }

  for (const auto& [delta_variable, delta_definition] :
       std::views::reverse(newton_system.delta_definitions)) {
    const auto index = delta_variable_index_.at(delta_variable);
    std::cout << "computing " << delta_variable->to_string() << " as "
              << delta_definition->to_string() << std::endl;
    result.at(index) = Evaluation::evaluate_vector(delta_definition, env_);
    for (auto& r : result.at(index)) {
      std::cout << r << ", ";
    }
    std::cout << std::endl;

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
  const auto val = eval_matrix_of_expressions_<Matrix>(v);
  return concatenate_matrices_(val);
}

Vector Optimizer::get_as_vector_(const std::vector<Expression::ExprPtr>& v) {
  const auto val = eval_vector_of_expressions_<Vector>(v);
  return concatenate_vectors_(val);
}

template <typename T>
std::vector<T> Optimizer::eval_vector_of_expressions_(
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
std::vector<std::vector<T>> Optimizer::eval_matrix_of_expressions_(
    const std::vector<std::vector<Expression::ExprPtr>>& v) {
  std::vector<std::vector<T>> result;
  result.reserve(v.size());
  std::ranges::transform(
      v, std::back_inserter(result),
      [this](const auto& vec) { return eval_vector_of_expressions_<T>(vec); });
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
                const auto& next_row = matrices.at((i + 1) %
           matrix_rows).at(j);
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