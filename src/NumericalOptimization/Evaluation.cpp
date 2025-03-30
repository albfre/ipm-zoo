#include "NumericalOptimization/Evaluation.h"

#include <algorithm>
#include <iostream>
#include <numeric>

#include "Utils/Assert.h"
#include "Utils/Helpers.h"

namespace NumericalOptimization::Evaluation {

namespace {
template <typename T>
concept VectorOrDiagonal = std::is_same_v<std::decay_t<T>, ValVector> ||
                           std::is_same_v<std::decay_t<T>, ValDiagMatrix>;

ValScalar dot_(const ValVector& x, const ValVector& y) {
  ASSERT(x.size() == y.size());
  return std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
}

EvalResult multiply_(const EvalResult& x, const EvalResult& y,
                     std::vector<EvalResult>& unhandled) {
  return Expression::match(x, y).with(
      [](const ValScalar& a, const auto& b) {
        return unary_op(b, [&](const auto& bi) { return a * bi; });
      },
      [](const ValVector& a, const ValVector& b) -> EvalResult {
        return dot_(a, b);
      },
      [](const auto& a, const auto& b)
        requires VectorOrDiagonal<decltype(a)> && VectorOrDiagonal<decltype(b)>
      { return elementwise_multiply(a, b); },
      [](const ValMatrix& a, const ValVector& b) -> EvalResult {
        auto r = ValVector(a.size());
        for (size_t i = 0; i < a.size(); ++i) {
          r[i] = dot_(a[i], b);
        }
        return r;
      },
      [&](const ValVector& a, const ValMatrix& b) -> EvalResult {
        unhandled.push_back(a);
        return b;
      },
      [&](const auto&, const auto&) {
        ASSERT(false);
        return x;
      });
}
}  // namespace

std::vector<std::vector<double>> evaluate_matrix(
    const Expression::ExprPtr& expr, const Environment& env) {
  auto val = evaluate(expr, env);
  return Expression::match(val).with(
      [](const auto&) {
        ASSERT(false);
        return std::vector<std::vector<double>>();
      },
      [](const ValDiagMatrix& diag) {
        auto result = std::vector<std::vector<double>>(
            diag.size(), std::vector<double>(diag.size()));
        for (size_t i = 0; i < diag.size(); ++i) {
          result[i][i] = diag[i];
        }
        return result;
      },
      [](const ValMatrix& mat) {
        auto result = std::vector<std::vector<double>>();
        result.reserve(mat.size());
        for (const auto& row : mat) {
          result.emplace_back(row.begin(), row.end());
        }
        return result;
      });
}

std::vector<double> evaluate_vector(const Expression::ExprPtr& expr,
                                    const Environment& env) {
  auto val = evaluate(expr, env);
  return Expression::match(val).with(
      [](const auto&) {
        ASSERT(false);
        return std::vector<double>();
      },
      [](const ValVector& x) { return std::vector<double>(x); },
      [](const ValDiagMatrix& x) { return std::vector<double>(x); });
}

double evaluate_scalar(const Expression::ExprPtr& expr,
                       const Environment& env) {
  auto val = evaluate(expr, env);
  return Expression::match(val).with(
      [](const auto&) {
        ASSERT(false);
        return 0.0;
      },
      [](const ValScalar& x) { return x; });
}

EvalResult evaluate(const Expression::ExprPtr& expr, const Environment& env) {
  using namespace Expression;

  // Return cached result if available
  if (env.contains(expr)) {
    return env.at(expr);
  }

  // Evaluate the expression based on its type
  auto res = match(expr).with(
      [&](const auto&) {
        std::cout << expr->to_string() << std::endl;
        ASSERT(env.contains(expr));
        return env.at(expr);
      },
      [&](const Number& x) { return val_scalar(x.value); },
      [&](const DiagonalMatrix& x) {
        const auto vec = evaluate(x.child, env);
        return match(vec).with(
            [](const ValVector& v) { return val_diag_matrix(v); },
            [](const auto&) {
              ASSERT(false);
              return EvalResult();
            });
      },
      [&](const Transpose& x) {
        const auto r = evaluate(x.child, env);
        return match(r).with([](const auto& y) -> EvalResult { return y; },
                             [](const ValMatrix& y) -> EvalResult {
                               const auto m = y.size();
                               const auto n = y.empty() ? 0 : y.at(0).size();
                               ValMatrix xT(n, ValVector(m));
                               for (size_t i = 0; i < m; ++i) {
                                 for (size_t j = 0; j < n; ++j) {
                                   xT[j][i] = y[i][j];
                                 }
                               }
                               return xT;
                             });
      },
      [&](const Invert& x) { return invert(evaluate(x.child, env)); },
      [&](const Log& x) {
        ASSERT(false);
        return evaluate(x.child, env);
      },
      [&env](const Negate& x) { return negate(evaluate(x.child, env)); },
      [&](const Sum& x) {
        auto res = evaluate(x.terms.at(0), env);
        for (auto& term : x.terms | std::views::drop(1)) {
          res = add(res, evaluate(term, env));
        }
        return res;
      },
      [&](const Product& x) {
        auto res = evaluate(x.terms.at(0), env);
        // Unhandled is used to avoid computing, e.g., x^T H
        // when computing x^T H x. Instead, x^T is saved in the "unhandled"
        // vector and H x is computed first.
        std::vector<EvalResult> unhandled;
        for (const auto& term : x.terms | std::views::drop(1)) {
          const auto next = evaluate(term, env);
          res = multiply_(res, next, unhandled);
        }
        std::ranges::reverse(unhandled);

        // Now, multiply the result with the unhandled terms (e.g., x^T
        // is multiplied with H x).
        std::vector<EvalResult> unhandled2;
        for (auto& init : unhandled) {
          res = multiply_(init, res, unhandled2);
        }
        ASSERT(unhandled2.empty());
        return res;
      });

  return res;
}

EvalResult unary_op(const EvalResult& x,
                    const std::function<double(double)>& lambda) {
  return Expression::match(x).with(
      [&](const ValScalar& y) -> EvalResult { return lambda(y); },
      [&](const auto& y) -> EvalResult
        requires VectorOrDiagonal<decltype(y)> {
          std::decay_t<decltype(y)> r;
          r.reserve(y.size());
          std::ranges::transform(
              y, std::back_inserter(r),
              [&lambda](const auto& yi) { return lambda(yi); });
          return r;
        },
      [&](const ValMatrix& y) -> EvalResult {
        auto m = y;
        for (auto& m_row : m) {
          for (auto& elem : m_row) {
            elem = lambda(elem);
          }
        }
        return m;
      });
}

EvalResult elementwise_op(const EvalResult& x, const EvalResult& y,
                          const std::function<double(double, double)>& lambda) {
  constexpr auto elementwise_vec_op = [](auto& v, const auto& a, const auto& b,
                                         const auto& lambda) {
    ASSERT(a.size() == b.size());
    for (size_t i = 0; i < v.size(); ++i) {
      v[i] = lambda(a[i], b[i]);
    }
    return v;
  };
  constexpr auto elementwise_diag_mat_op = [](const auto& m, const auto& d,
                                              const auto& lambda) {
    ASSERT(m.size() == d.size());
    auto new_mat = m;
    for (size_t i = 0; i < m.size(); ++i) {
      new_mat[i][i] = lambda(m[i][i], d[i]);
    }
    return new_mat;
  };

  return Expression::match(x, y).with(
      [&lambda](const ValScalar& a, const ValScalar& b) -> EvalResult {
        return lambda(a, b);
      },
      [&](const ValDiagMatrix& a, const ValDiagMatrix& b) -> EvalResult {
        auto v = ValDiagMatrix(b.size());
        return elementwise_vec_op(v, a, b, lambda);
      },
      [&](const auto& a, const auto& b) -> EvalResult
        requires VectorOrDiagonal<decltype(a)> && VectorOrDiagonal<decltype(b)>
      {
        auto v = ValVector(b.size());
        return elementwise_vec_op(v, a, b, lambda);
      },
      [&](const ValMatrix& a, const ValDiagMatrix& b) -> EvalResult {
        return elementwise_diag_mat_op(a, b, lambda);
      },
      [&](const ValDiagMatrix& a, const ValMatrix& b) -> EvalResult {
        return elementwise_diag_mat_op(b, a, lambda);
      },
      [&lambda](const ValMatrix& a, const ValMatrix& b) -> EvalResult {
        ASSERT(a.size() == b.size());
        auto v = ValMatrix(a.size(), ValVector(a.empty() ? 0 : a.at(0).size()));
        for (size_t i = 0; i < a.size(); ++i) {
          ASSERT(a[i].size() == b[i].size());
          for (size_t j = 0; j < a[i].size(); ++j) {
            v[i][j] = lambda(a[i][j], b[i][j]);
          }
        }
        return v;
      },
      [&](const auto&, const auto&) {
        ASSERT(false);
        return x;
      });
}

EvalResult negate(const EvalResult& x) {
  return unary_op(x, [](const auto& xi) { return -xi; });
}

EvalResult invert(const EvalResult& x) {
  return unary_op(x, [](const auto& xi) { return 1.0 / xi; });
}

EvalResult add(const EvalResult& x, const EvalResult& y) {
  return elementwise_op(x, y, [](const auto a, const auto b) { return a + b; });
}

EvalResult subtract(const EvalResult& x, const EvalResult& y) {
  return elementwise_op(x, y, [](const auto a, const auto b) { return a - b; });
}

EvalResult multiply(const EvalResult& x, const EvalResult& y) {
  std::vector<EvalResult> unhandled;
  auto res = multiply_(x, y, unhandled);
  ASSERT(unhandled.empty());
  return res;
}

EvalResult elementwise_multiply(const EvalResult& x, const EvalResult& y) {
  return elementwise_op(x, y, [](const auto a, const auto b) { return a * b; });
}

EvalResult elementwise_divide(const EvalResult& x, const EvalResult& y) {
  return elementwise_op(x, y, [](const auto a, const auto b) { return a / b; });
}

EvalResult val_scalar(double x) { return ValScalar(x); }

EvalResult val_vector(const std::vector<double>& v) {
  return ValVector(v.begin(), v.end());
}

EvalResult val_diag_matrix(const std::vector<double>& v) {
  return ValDiagMatrix(v.begin(), v.end());
}

EvalResult val_matrix(const std::vector<std::vector<double>>& m) {
  ValMatrix result;
  result.reserve(m.size());

  for (const auto& row : m) {
    result.emplace_back(row.begin(), row.end());
  }

  return result;
}

}  // namespace NumericalOptimization::Evaluation