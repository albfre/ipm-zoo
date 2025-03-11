#pragma once
#include <ranges>
#include <sstream>
#include "Expression.h"
#include "Helpers.h"

namespace Expression {
struct ToStringVisitor {
  bool condensed = false;
  std::string operator()(const Number& x) const {
    std::stringstream ss;
    ss << x.value;
    return ss.str();
  }

  std::string operator()(const auto& x) const {
    if constexpr (std::is_base_of_v<NamedNullaryExpr,
                                    std::decay_t<decltype(x)>>) {
      return x.name;
    } else {
      static_assert(always_false_v<decltype(x)>);
    }
  }

  std::string operator()(const Transpose& x) const {
    const auto& c = *x.child;
    if (condensed && (is<Sum>(c) || is<Product>(c) || is<Invert>(c))) {
      return "(" + c.toString(condensed) + ")^T";
    }
    return c.toString(condensed) + "^T";
  }

  std::string operator()(const Negate& x) const {
    const auto& c = *x.child;
    if (condensed && is<Sum>(c)) {
      return "-(" + c.toString(condensed) + ")";
    }
    return "-" + c.toString(condensed);
  }

  std::string operator()(const Invert& x) const {
    const auto& c = *x.child;
    if (condensed && (is<Sum>(c) || is<Product>(c) || is<Transpose>(c))) {
      return "(" + c.toString(condensed) + ")^{-1}";
    }
    return c.toString(condensed) + "^{-1}";
  }

  std::string operator()(const Log& x) const {
    return "\\log(" + x.child->toString(condensed) + ")";
  }

  std::string operator()(const Sum& x) const {
    std::stringstream ss;
    ss << x.terms.front().toString(condensed);
    for (const auto& t : x.terms | std::views::drop(1)) {
      ss << (is<Negate>(t) ? " - " : " + ") << t.toString(condensed);
    }
    return condensed ? "(" + ss.str() + ")" : ss.str();
  }

  std::string operator()(const Product& x) const {
    std::stringstream ss;
    ss << x.terms.front().toString(condensed);
    const auto symbol = condensed ? " " : " * ";
    for (const auto& t : x.terms | std::views::drop(1)) {
      if (is<Negate>(t) || (condensed && is<Sum>(t))) {
        ss << symbol << "(" << t.toString(condensed) << ")";
      } else {
        ss << symbol << t.toString(condensed);
      }
    }
    return condensed ? "(" + ss.str() + ")" : ss.str();
  }
};

struct ToExpressionStringVisitor {
  std::string operator()(const Number& x) const {
    std::stringstream ss;
    ss << "number(" << x.value << ")";
    return ss.str();
  }

  std::string operator()(const NamedConstant& x) const {
    return "namedConstant(" + x.name + ")";
  }

  std::string operator()(const Matrix& x) const {
    return "matrix(" + x.name + ")";
  }

  std::string operator()(const SymmetricMatrix& x) const {
    return "symmetricMatrix(" + x.name + ")";
  }

  std::string operator()(const Variable& x) const {
    return "variable(" + x.name + ")";
  }

  std::string operator()(const Transpose& x) const {
    return "transpose(" + x.child->toExpressionString() + ")";
  }

  std::string operator()(const Negate& x) const {
    return "negate(" + x.child->toExpressionString() + ")";
  }

  std::string operator()(const Invert& x) const {
    return "invert(" + x.child->toExpressionString() + ")";
  }

  std::string operator()(const Log& x) const {
    return "log(" + x.child->toExpressionString() + ")";
  }

  std::string operator()(const Sum& x) const {
    return "sum(" + termsToString_(x.terms) + ")";
  }

  std::string operator()(const Product& x) const {
    return "product(" + termsToString_(x.terms) + ")";
  }

 private:
  std::string termsToString_(const std::vector<Expr>& terms) const {
    std::stringstream ss;
    ss << terms.front().toExpressionString();
    for (const auto& t : terms | std::views::drop(1)) {
      ss << ", " << t.toExpressionString();
    }
    return ss.str();
  }
};
}  // namespace Expression