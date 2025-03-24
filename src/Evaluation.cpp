#include "Evaluation.h"

#include <algorithm>
#include <numeric>
#include <ranges>

#include "Assert.h"
#include "Helpers.h"

namespace Evaluation {

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
        return unaryOp(b, [&](const auto& bi) { return a * bi; });
      },
      [](const ValVector& a, const ValVector& b) -> EvalResult {
        return dot_(a, b);
      },
      [](const auto& a, const auto& b)
        requires VectorOrDiagonal<decltype(a)> && VectorOrDiagonal<decltype(b)>
      { return elementwiseMultiply(a, b); },
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

EvalResult evaluate(const Expression::ExprPtr& expr, Environment& env) {
  using namespace Expression;

  // Return cached result if available
  if (env.contains(expr)) {
    return env.at(expr);
  }

  // Evaluate the expression based on its type
  auto res = match(expr).with(
      [&](const auto&) { return env.at(expr); },
      [&](const Number& x) { return valScalar(x.value); },
      [&](const DiagonalMatrix& x) {
        const auto vec = evaluate(x.child, env);
        return match(vec).with(
            [](const ValVector& v) { return valDiagMatrix(v); },
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

  // Cache and return the result
  env[expr] = res;
  return res;
}

EvalResult unaryOp(const EvalResult& x,
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
        for (auto& mRow : m) {
          for (auto& elem : mRow) {
            elem = lambda(elem);
          }
        }
        return m;
      });
}

EvalResult elementwiseOp(const EvalResult& x, const EvalResult& y,
                         const std::function<double(double, double)>& lambda) {
  constexpr auto elementwiseVecOp = [](auto& v, const auto& a, const auto& b,
                                       const auto& lambda) {
    ASSERT(a.size() == b.size());
    for (size_t i = 0; i < v.size(); ++i) {
      v[i] = lambda(a[i], b[i]);
    }
    return v;
  };
  constexpr auto elementwiseDiagMatOp = [](const auto& m, const auto& d,
                                           const auto& lambda) {
    ASSERT(m.size() == d.size());
    auto newMat = m;
    for (size_t i = 0; i < m.size(); ++i) {
      newMat[i][i] = lambda(m[i][i], d[i]);
    }
    return newMat;
  };

  return Expression::match(x, y).with(
      [&lambda](const ValScalar& a, const ValScalar& b) -> EvalResult {
        return lambda(a, b);
      },
      [&](const ValDiagMatrix& a, const ValDiagMatrix& b) -> EvalResult {
        auto v = ValDiagMatrix(b.size());
        return elementwiseVecOp(v, a, b, lambda);
      },
      [&](const auto& a, const auto& b) -> EvalResult
        requires VectorOrDiagonal<decltype(a)> && VectorOrDiagonal<decltype(b)>
      {
        auto v = ValVector(b.size());
        return elementwiseVecOp(v, a, b, lambda);
      },
      [&](const ValMatrix& a, const ValDiagMatrix& b) -> EvalResult {
        return elementwiseDiagMatOp(a, b, lambda);
      },
      [&](const ValDiagMatrix& a, const ValMatrix& b) -> EvalResult {
        return elementwiseDiagMatOp(b, a, lambda);
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
  return unaryOp(x, [](const auto& xi) { return -xi; });
}

EvalResult invert(const EvalResult& x) {
  return unaryOp(x, [](const auto& xi) { return 1.0 / xi; });
}

EvalResult add(const EvalResult& x, const EvalResult& y) {
  return elementwiseOp(x, y, [](const auto a, const auto b) { return a + b; });
}

EvalResult subtract(const EvalResult& x, const EvalResult& y) {
  return elementwiseOp(x, y, [](const auto a, const auto b) { return a - b; });
}

EvalResult multiply(const EvalResult& x, const EvalResult& y) {
  std::vector<EvalResult> unhandled;
  auto res = multiply_(x, y, unhandled);
  ASSERT(unhandled.empty());
  return res;
}

EvalResult elementwiseMultiply(const EvalResult& x, const EvalResult& y) {
  return elementwiseOp(x, y, [](const auto a, const auto b) { return a * b; });
}

EvalResult elementwiseDivide(const EvalResult& x, const EvalResult& y) {
  return elementwiseOp(x, y, [](const auto a, const auto b) { return a / b; });
}

EvalResult valScalar(double x) { return ValScalar(x); }

EvalResult valVector(const std::vector<double>& v) {
  return ValVector(v.begin(), v.end());
}

EvalResult valDiagMatrix(const std::vector<double>& v) {
  return ValDiagMatrix(v.begin(), v.end());
}

}  // namespace Evaluation