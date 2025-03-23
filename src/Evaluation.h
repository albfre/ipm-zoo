#pragma once
#include <algorithm>
#include <map>
#include <numeric>
#include <ranges>

#include "Assert.h"
#include "Expression.h"
#include "Helpers.h"

namespace Evaluation {
using ValScalar = double;
struct ValVector : std::vector<double> {
  using std::vector<double>::vector;
};
struct ValDiagMatrix : std::vector<double> {
  using std::vector<double>::vector;
};
using ValMatrix = std::vector<ValVector>;
using EvalResult = std::variant<ValScalar, ValVector, ValDiagMatrix, ValMatrix>;
template <typename T>
constexpr bool is_vector_or_diag_v =
    std::is_same_v<std::decay_t<T>, ValVector> ||
    std::is_same_v<std::decay_t<T>, ValDiagMatrix>;

using Environment = std::map<Expression::Expr, EvalResult>;

// Forward declarations
ValScalar dot(const ValVector& x, const ValVector& y);
EvalResult product(const EvalResult& x, const EvalResult& y,
                   std::vector<EvalResult>& unhandled);
EvalResult matrixVectorProduct(const ValMatrix& m, const ValVector& v);
EvalResult negate(const EvalResult& x);
EvalResult invert(const EvalResult& x);
EvalResult add(const EvalResult& x, const EvalResult& y);
EvalResult subtract(const EvalResult& x, const EvalResult& y);
EvalResult elementwiseMultiply(const EvalResult& x, const EvalResult& y);
EvalResult elementwiseDivide(const EvalResult& x, const EvalResult& y);
EvalResult valScalar(double x);
EvalResult valVector(const std::vector<double>& v);
EvalResult valDiagMatrix(const std::vector<double>& v);

template <typename TLambda>
EvalResult unaryOp(const EvalResult& x, const TLambda& lambda) {
  return Expression::match(x).with(
      [&](const ValScalar& y) -> EvalResult { return lambda(y); },
      [&](const auto& y) -> EvalResult {
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

EvalResult evaluate(const Expression::Expr& expr, Environment& env) {
  using namespace Expression;
  if (env.contains(expr)) {
    return env.at(expr);
  }
  auto res = match(expr).with(
      [&](const auto&) { return env.at(expr); },
      [&](const Number& x) { return valScalar(x.value); },
      [&](const DiagonalMatrix& x) {
        const auto vec = evaluate(*x.child, env);
        return match(vec).with(
            [](const ValVector& v) { return valDiagMatrix(v); },
            [](const auto&) {
              ASSERT(false);
              return EvalResult();
            });
      },
      [&](const Transpose& x) {
        const auto r = evaluate(*x.child, env);
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
      [&](const Invert& x) { return invert(evaluate(*x.child, env)); },
      [&](const Log& x) {
        ASSERT(false);
        return evaluate(*x.child, env);
      },
      [&env](const Negate& x) { return negate(evaluate(*x.child, env)); },
      [&](const Sum& x) {
        auto res = evaluate(x.terms.at(0), env);
        for (auto& term : x.terms | std::views::drop(1)) {
          res = add(res, evaluate(term, env));
        }
        return res;
      },
      [&](const Product& x) {
        auto res = evaluate(x.terms.at(0), env);
        std::vector<EvalResult> unhandled;
        for (const auto& term : x.terms | std::views::drop(1)) {
          const auto next = evaluate(term, env);
          res = product(res, next, unhandled);
        }
        std::ranges::reverse(unhandled);
        std::vector<EvalResult> unhandled2;
        for (auto& init : unhandled) {
          res = product(init, res, unhandled2);
        }
        ASSERT(unhandled2.empty());
        return res;
      });
  env[expr] = res;
  return res;
}

EvalResult product(const EvalResult& x, const EvalResult& y,
                   std::vector<EvalResult>& unhandled) {
  return Expression::match(x, y).with(
      [](const ValScalar& a, const auto& b) {
        return unaryOp(b, [&](const auto& bi) { return a * bi; });
      },
      [](const ValVector& a, const ValVector& b) -> EvalResult {
        return dot(a, b);
      },
      [](const auto& a, const auto& b)
        requires(is_vector_or_diag_v<decltype(a)> &&
                 is_vector_or_diag_v<decltype(b)>)
      { return elementwiseMultiply(a, b); },
      [](const ValMatrix& a, const ValVector& b) {
        return matrixVectorProduct(a, b);
      },
      [&](const ValVector& a, const ValMatrix& b) -> EvalResult {
        unhandled.push_back(a);
        return b;
      },
      [](const auto& a, const auto& b) -> EvalResult {
        ASSERT(false);
        return a;
      });
}

EvalResult matrixVectorProduct(const ValMatrix& m, const ValVector& v) {
  auto r = ValVector(m.size());
  for (size_t i = 0; i < m.size(); ++i) {
    r[i] = dot(m[i], v);
  }
  return r;
}

ValScalar dot(const ValVector& x, const ValVector& y) {
  ASSERT(x.size() == y.size());
  return std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
}

EvalResult negate(const EvalResult& x) {
  return unaryOp(x, [](const auto& xi) { return -xi; });
}

EvalResult invert(const EvalResult& x) {
  return unaryOp(x, [](const auto& xi) { return 1.0 / xi; });
}

template <typename TLambda>
EvalResult elementwiseOp(const EvalResult& a, const EvalResult& b,
                         const TLambda& lambda) {
  const auto elementwiseVecOp = [&lambda](auto& v, const auto& av,
                                          const auto& bv) {
    ASSERT(av.size() == bv.size());
    for (size_t i = 0; i < v.size(); ++i) {
      v[i] = lambda(av[i], bv[i]);
    }
    return v;
  };
  const auto elementwiseDiagMatOp = [&lambda](const auto& m, const auto& d) {
    ASSERT(m.size() == d.size());
    auto newMat = m;
    for (size_t i = 0; i < m.size(); ++i) {
      newMat[i][i] = lambda(m[i][i], d[i]);
    }
    return newMat;
  };

  return Expression::match(a, b).with(
      [&](const ValScalar& av, const ValScalar& bv) -> EvalResult {
        return lambda(av, bv);
      },
      [&](const ValDiagMatrix& av, const ValDiagMatrix& bv) -> EvalResult {
        auto v = ValDiagMatrix(bv.size());
        return elementwiseVecOp(v, av, bv);
      },
      [&](const ValVector& av, const ValDiagMatrix& bv) -> EvalResult {
        auto v = ValVector(bv.size());
        return elementwiseVecOp(v, av, bv);
      },
      [&](const ValDiagMatrix& av, const ValVector& bv) -> EvalResult {
        auto v = ValVector(bv.size());
        return elementwiseVecOp(v, av, bv);
      },
      [&](const ValMatrix& av, const ValDiagMatrix& bv) -> EvalResult {
        return elementwiseDiagMatOp(av, bv);
      },
      [&](const ValDiagMatrix& av, const ValMatrix& bv) -> EvalResult {
        return elementwiseDiagMatOp(bv, av);
      },
      [&](const ValMatrix& av, const ValMatrix& bv) -> EvalResult {
        ASSERT(av.size() == bv.size());
        auto v =
            ValMatrix(av.size(), ValVector(av.empty() ? 0 : av.at(0).size()));
        for (size_t i = 0; i < av.size(); ++i) {
          ASSERT(av[i].size() == bv[i].size());
          for (size_t j = 0; j < av[i].size(); ++j) {
            v[i][j] = lambda(av[i][j], bv[i][j]);
          }
        }
        return v;
      },
      [&](const auto&, const auto&) {
        ASSERT(false);
        return b;
      });
}

EvalResult add(const EvalResult& x, const EvalResult& y) {
  return elementwiseOp(x, y, [](const auto a, const auto b) { return a + b; });
}

EvalResult subtract(const EvalResult& x, const EvalResult& y) {
  return elementwiseOp(x, y, [](const auto a, const auto b) { return a - b; });
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