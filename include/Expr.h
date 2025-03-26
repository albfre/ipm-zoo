#pragma once
#include <memory>
#include <set>
#include <string>
#include <variant>
#include <vector>

namespace Expression {

struct Number;
struct NamedNullaryExpr;
struct NamedScalar;
struct NamedVector;
struct Variable;
struct Matrix;
struct SymmetricMatrix;
struct UnaryExpr;
struct DiagonalMatrix;
struct Transpose;
struct Invert;
struct Log;
struct Negate;
struct NaryExpr;
struct Sum;
struct Product;
class Expr;

using ExprPtr = std::shared_ptr<const Expr>;

class Expr {
 public:
  using ExprVariant =
      std::variant<Number, NamedScalar, NamedVector, Variable, Matrix,
                   SymmetricMatrix, DiagonalMatrix, Transpose, Invert, Log, Sum,
                   Product, Negate>;

  Expr(const Expr& other) = delete;
  Expr(Expr&& other) noexcept = delete;
  Expr& operator=(const Expr& other) = delete;

  [[nodiscard]] ExprPtr differentiate(const ExprPtr& var) const;
  [[nodiscard]] ExprPtr simplify(bool distribute = true) const;
  [[nodiscard]] ExprPtr simplify_once(bool distribute = true) const;
  [[nodiscard]] std::string to_string(bool condensed = false) const;
  [[nodiscard]] const std::string& to_expression_string() const;
  [[nodiscard]] bool contains_subexpression(const ExprPtr& expr) const;
  [[nodiscard]] ExprPtr replace_subexpression(const ExprPtr& expr,
                                              const ExprPtr& replacement) const;
  [[nodiscard]] const ExprVariant& get_impl() const;
  [[nodiscard]] ExprPtr get_leading_or_ending_factor(bool leading) const;
  [[nodiscard]] ExprPtr factor_out(const ExprPtr& factor, bool leading) const;
  [[nodiscard]] double complexity() const;
  [[nodiscard]] std::set<ExprPtr> get_variables() const;

 private:
  friend class ExprFactory;

  template <typename T>
  explicit Expr(T value, std::string_view expression_string)
      : impl_(std::make_shared<ExprVariant>(std::move(value))),
        expression_string_(std::string(expression_string)) {}
  const std::shared_ptr<const ExprVariant> impl_;
  const std::string expression_string_;
};

struct NamedNullaryExpr {
  const std::string name;
  explicit NamedNullaryExpr(std::string_view name) : name(std::string(name)) {}
  virtual ~NamedNullaryExpr() = default;
};

struct UnaryExpr {
  const ExprPtr child;
  explicit UnaryExpr(ExprPtr expr) : child(std::move(expr)) {}
  virtual ~UnaryExpr() = default;
};

struct NaryExpr {
  const std::vector<ExprPtr> terms;
  NaryExpr(std::vector<ExprPtr> terms) : terms(std::move(terms)) {}
  virtual ~NaryExpr() = default;
};

struct Number {
  const double value;
  explicit Number(double val) : value(val) {}
};
struct NamedScalar : public NamedNullaryExpr {
  using NamedNullaryExpr::NamedNullaryExpr;
};
struct NamedVector : public NamedNullaryExpr {
  using NamedNullaryExpr::NamedNullaryExpr;
};
struct Variable : public NamedNullaryExpr {
  using NamedNullaryExpr::NamedNullaryExpr;
};
struct Matrix : public NamedNullaryExpr {
  using NamedNullaryExpr::NamedNullaryExpr;
};
struct SymmetricMatrix : public NamedNullaryExpr {
  using NamedNullaryExpr::NamedNullaryExpr;
};
struct DiagonalMatrix : public UnaryExpr {
  using UnaryExpr::UnaryExpr;
};
struct Transpose : public UnaryExpr {
  using UnaryExpr::UnaryExpr;
};
struct Invert : public UnaryExpr {
  using UnaryExpr::UnaryExpr;
};
struct Log : public UnaryExpr {
  using UnaryExpr::UnaryExpr;
};
struct Negate : public UnaryExpr {
  using UnaryExpr::UnaryExpr;
};
struct Sum : public NaryExpr {
  using NaryExpr::NaryExpr;
};
struct Product : public NaryExpr {
  using NaryExpr::NaryExpr;
};

std::strong_ordering operator<=>(const ExprPtr& left, const ExprPtr& right);
std::strong_ordering operator<=>(const Expr& left, const Expr& right);
bool operator==(const ExprPtr& left, const ExprPtr& right);
bool operator==(const Expr& left, const Expr& right);
std::ostream& operator<<(std::ostream& os, const Expr& expr);
struct ExprHash {
  std::size_t operator()(const Expr& expr) const;
};
}  // namespace Expression
