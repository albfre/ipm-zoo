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
    if constexpr (NamedNullaryType<decltype(x)>) {
      return x.name;
    } else {
      static_assert(always_false_v<decltype(x)>);
    }
  }

  std::string operator()(const DiagonalMatrix& x) const {
    const auto& c = *x.child;
    return match(c).with(
        [](const Variable& x) {
          auto name = x.name;
          for (size_t i = 0; i < name.size(); ++i) {
            if (std::isalpha(name[i])) {
              name[i] = std::toupper(name[i]);
              return name;
            }
          }
          return "\\diag(" + name + ")";
        },
        [&](const auto&) { return "\\diag(" + c.toString(condensed) + ")"; });
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
    ss << (condensed ? "" : "(");
    ss << x.terms.front().toString(condensed);
    for (const auto& t : x.terms | std::views::drop(1)) {
      if (is<Negate>(t)) {
        ss << " - " << std::get<Negate>(t.getImpl()).child->toString(condensed);
      } else {
        ss << " + " << t.toString(condensed);
      }
    }
    ss << (condensed ? "" : ")");
    return ss.str();
  }

  std::string operator()(const Product& x) const {
    std::stringstream ss;
    ss << (condensed ? "" : "(");
    const auto& front = x.terms.front();
    if (condensed && is<Sum>(front)) {
      ss << "(" << front.toString(condensed) << ")";
    } else {
      ss << front.toString(condensed);
    }
    const auto symbol = condensed ? " " : " * ";
    for (const auto& t : x.terms | std::views::drop(1)) {
      if (is<Negate>(t) || (condensed && is<Sum>(t))) {
        ss << symbol << "(" << t.toString(condensed) << ")";
      } else {
        ss << symbol << t.toString(condensed);
      }
    }
    ss << (condensed ? "" : ")");
    return ss.str();
  }
};

struct ToExpressionStringVisitor {
  std::string operator()(const Number& x) const {
    std::stringstream ss;
    ss << "number(" << x.value << ")";
    return ss.str();
  }

  std::string operator()(const NamedScalar& x) const {
    return "namedScalar(" + x.name + ")";
  }

  std::string operator()(const NamedVector& x) const {
    return "namedVector(" + x.name + ")";
  }

  std::string operator()(const Variable& x) const {
    return "variable(" + x.name + ")";
  }

  std::string operator()(const Matrix& x) const {
    return "matrix(" + x.name + ")";
  }

  std::string operator()(const SymmetricMatrix& x) const {
    return "symmetricMatrix(" + x.name + ")";
  }

  std::string operator()(const DiagonalMatrix& x) const {
    return "diagonalMatrix(" + x.child->toExpressionString() + ")";
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