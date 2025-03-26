#include "Visitors/ToStringVisitor.h"

#include <algorithm>
#include <cctype>
#include <ranges>
#include <sstream>

#include "Utils/Helpers.h"

namespace Expression {

ToStringVisitor::ToStringVisitor(bool condensed) : condensed_(condensed) {}

std::string ToStringVisitor::operator()(const Number& x) const {
  std::stringstream ss;
  ss << x.value;
  return ss.str();
}

std::string ToStringVisitor::operator()(const DiagonalMatrix& x) const {
  const auto& c = x.child;
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
      [&](const auto&) { return "\\diag(" + c->toString(condensed_) + ")"; });
}

std::string ToStringVisitor::operator()(const Transpose& x) const {
  const auto& c = x.child;
  if (condensed_ && (is<Sum>(c) || is<Product>(c) || is<Invert>(c))) {
    return "(" + c->toString(condensed_) + ")^T";
  }
  return c->toString(condensed_) + "^T";
}

std::string ToStringVisitor::operator()(const Negate& x) const {
  const auto& c = x.child;
  if (condensed_ && is<Sum>(c)) {
    return "-(" + c->toString(condensed_) + ")";
  }
  return "-" + c->toString(condensed_);
}

std::string ToStringVisitor::operator()(const Invert& x) const {
  const auto& c = x.child;
  if (condensed_ && (is<Sum>(c) || is<Product>(c) || is<Transpose>(c))) {
    return "(" + c->toString(condensed_) + ")^{-1}";
  }
  return c->toString(condensed_) + "^{-1}";
}

std::string ToStringVisitor::operator()(const Log& x) const {
  return "\\log(" + x.child->toString(condensed_) + ")";
}

std::string ToStringVisitor::operator()(const Sum& x) const {
  std::stringstream ss;
  ss << (condensed_ ? "" : "(");
  ss << x.terms.front()->toString(condensed_);
  for (const auto& t : x.terms | std::views::drop(1)) {
    if (is<Negate>(t)) {
      ss << " - " << std::get<Negate>(t->getImpl()).child->toString(condensed_);
    } else {
      ss << " + " << t->toString(condensed_);
    }
  }
  ss << (condensed_ ? "" : ")");
  return ss.str();
}

std::string ToStringVisitor::operator()(const Product& x) const {
  std::stringstream ss;
  ss << (condensed_ ? "" : "(");
  const auto& front = x.terms.front();
  if (condensed_ && is<Sum>(front)) {
    ss << "(" << front->toString(condensed_) << ")";
  } else {
    ss << front->toString(condensed_);
  }
  const auto symbol = condensed_ ? " " : " * ";
  for (const auto& t : x.terms | std::views::drop(1)) {
    if (is<Negate>(t) || (condensed_ && is<Sum>(t))) {
      ss << symbol << "(" << t->toString(condensed_) << ")";
    } else {
      ss << symbol << t->toString(condensed_);
    }
  }
  ss << (condensed_ ? "" : ")");
  return ss.str();
}

// ToExpressionStringVisitor implementation

std::string ToExpressionStringVisitor::operator()(const Number& x) const {
  std::stringstream ss;
  ss << "number(" << x.value << ")";
  return ss.str();
}

std::string ToExpressionStringVisitor::operator()(const NamedScalar& x) const {
  return "namedScalar(" + x.name + ")";
}

std::string ToExpressionStringVisitor::operator()(const NamedVector& x) const {
  return "namedVector(" + x.name + ")";
}

std::string ToExpressionStringVisitor::operator()(const Variable& x) const {
  return "variable(" + x.name + ")";
}

std::string ToExpressionStringVisitor::operator()(const Matrix& x) const {
  return "matrix(" + x.name + ")";
}

std::string ToExpressionStringVisitor::operator()(
    const SymmetricMatrix& x) const {
  return "symmetricMatrix(" + x.name + ")";
}

std::string ToExpressionStringVisitor::operator()(
    const DiagonalMatrix& x) const {
  return "diagonalMatrix(" + x.child->toExpressionString() + ")";
}

std::string ToExpressionStringVisitor::operator()(const Transpose& x) const {
  return "transpose(" + x.child->toExpressionString() + ")";
}

std::string ToExpressionStringVisitor::operator()(const Negate& x) const {
  return "negate(" + x.child->toExpressionString() + ")";
}

std::string ToExpressionStringVisitor::operator()(const Invert& x) const {
  return "invert(" + x.child->toExpressionString() + ")";
}

std::string ToExpressionStringVisitor::operator()(const Log& x) const {
  return "log(" + x.child->toExpressionString() + ")";
}

std::string ToExpressionStringVisitor::operator()(const Sum& x) const {
  return "sum(" + termsToString_(x.terms) + ")";
}

std::string ToExpressionStringVisitor::operator()(const Product& x) const {
  return "product(" + termsToString_(x.terms) + ")";
}

std::string ToExpressionStringVisitor::termsToString_(
    const std::vector<ExprPtr>& terms) const {
  std::stringstream ss;
  ss << terms.front()->toExpressionString();
  for (const auto& t : terms | std::views::drop(1)) {
    ss << ", " << t->toExpressionString();
  }
  return ss.str();
}

}  // namespace Expression
