#include "Visitors/ToStringVisitor.h"

#include <algorithm>
#include <cctype>
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
      [](const auto& x)
        requires NamedNullaryType<decltype(x)> {
          auto name = x.name;
          for (size_t i = 0; i < name.size(); ++i) {
            if (std::isalpha(name[i])) {
              name[i] = std::toupper(name[i]);
              return name;
            }
          }
          return "\\diag(" + name + ")";
        },
      [&](const auto&) { return "\\diag(" + c->to_string(condensed_) + ")"; });
}

std::string ToStringVisitor::operator()(const Transpose& x) const {
  const auto& c = x.child;
  if (condensed_ && (is<Sum>(c) || is<Product>(c) || is<Invert>(c))) {
    return "(" + c->to_string(condensed_) + ")^T";
  }
  return c->to_string(condensed_) + "^T";
}

std::string ToStringVisitor::operator()(const Negate& x) const {
  const auto& c = x.child;
  if (condensed_ && is<Sum>(c)) {
    return "-(" + c->to_string(condensed_) + ")";
  }
  return "-" + c->to_string(condensed_);
}

std::string ToStringVisitor::operator()(const Invert& x) const {
  const auto& c = x.child;
  if (condensed_ && (is<Sum>(c) || is<Product>(c) || is<Transpose>(c))) {
    return "(" + c->to_string(condensed_) + ")^{-1}";
  }
  return c->to_string(condensed_) + "^{-1}";
}

std::string ToStringVisitor::operator()(const Log& x) const {
  return "\\log(" + x.child->to_string(condensed_) + ")";
}

std::string ToStringVisitor::operator()(const Sum& x) const {
  std::stringstream ss;
  ss << (condensed_ ? "" : "(");
  ss << x.terms.front()->to_string(condensed_);
  for (const auto& t : x.terms | std::views::drop(1)) {
    if (is<Negate>(t)) {
      ss << " - "
         << std::get<Negate>(t->get_impl()).child->to_string(condensed_);
    } else {
      ss << " + " << t->to_string(condensed_);
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
    ss << "(" << front->to_string(condensed_) << ")";
  } else {
    ss << front->to_string(condensed_);
  }
  const auto symbol = condensed_ ? " " : " * ";
  for (const auto& t : x.terms | std::views::drop(1)) {
    if (is<Negate>(t) || (condensed_ && is<Sum>(t))) {
      ss << symbol << "(" << t->to_string(condensed_) << ")";
    } else {
      ss << symbol << t->to_string(condensed_);
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
  return "named_scalar(" + x.name + ")";
}

std::string ToExpressionStringVisitor::operator()(const NamedVector& x) const {
  return "named_vector(" + x.name + ")";
}

std::string ToExpressionStringVisitor::operator()(const Variable& x) const {
  return "variable(" + x.name + ")";
}

std::string ToExpressionStringVisitor::operator()(const Matrix& x) const {
  return "matrix(" + x.name + ")";
}

std::string ToExpressionStringVisitor::operator()(
    const SymmetricMatrix& x) const {
  return "symmetric_matrix(" + x.name + ")";
}

std::string ToExpressionStringVisitor::operator()(
    const DiagonalMatrix& x) const {
  return "diagonal_matrix(" + x.child->to_expression_string() + ")";
}

std::string ToExpressionStringVisitor::operator()(const Transpose& x) const {
  return "transpose(" + x.child->to_expression_string() + ")";
}

std::string ToExpressionStringVisitor::operator()(const Negate& x) const {
  return "negate(" + x.child->to_expression_string() + ")";
}

std::string ToExpressionStringVisitor::operator()(const Invert& x) const {
  return "invert(" + x.child->to_expression_string() + ")";
}

std::string ToExpressionStringVisitor::operator()(const Log& x) const {
  return "log(" + x.child->to_expression_string() + ")";
}

std::string ToExpressionStringVisitor::operator()(const Sum& x) const {
  return "sum(" + terms_to_string_(x.terms) + ")";
}

std::string ToExpressionStringVisitor::operator()(const Product& x) const {
  return "product(" + terms_to_string_(x.terms) + ")";
}

std::string ToExpressionStringVisitor::terms_to_string_(
    const std::vector<ExprPtr>& terms) const {
  std::stringstream ss;
  ss << terms.front()->to_expression_string();
  for (const auto& t : terms | std::views::drop(1)) {
    ss << ", " << t->to_expression_string();
  }
  return ss.str();
}

}  // namespace Expression
