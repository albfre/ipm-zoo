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

class Expr {
  using ExprVariant =
      std::variant<Number, NamedScalar, NamedVector, Variable, Matrix,
                   SymmetricMatrix, DiagonalMatrix, Transpose, Invert, Log, Sum,
                   Product, Negate>;

 public:
  template <typename T>
    requires std::is_base_of_v<UnaryExpr, std::decay_t<T>>
  explicit Expr(T value)
      : impl_(std::make_shared<ExprVariant>(std::move(value))) {}

  template <typename T>
    requires(!std::is_same_v<Expr, std::decay_t<T>> &&
             !std::is_base_of_v<UnaryExpr, std::decay_t<T>>)
  explicit Expr(T&& value)
      : impl_(std::make_shared<ExprVariant>(std::forward<T>(value))) {}

  Expr(const Expr& other) = default;
  Expr(Expr&& other) noexcept = default;
  Expr& operator=(const Expr& other) = default;

  [[nodiscard]] Expr differentiate(const Expr& var) const;
  [[nodiscard]] Expr simplify(bool distribute = true) const;
  [[nodiscard]] Expr simplifyOnce(bool distribute = true) const;
  [[nodiscard]] std::string toString(bool condensed = false) const;
  [[nodiscard]] std::string toExpressionString() const;
  [[nodiscard]] bool containsSubexpression(const Expr& expr) const;
  [[nodiscard]] Expr replaceSubexpression(const Expr& expr,
                                          const Expr& replacement) const;
  [[nodiscard]] const ExprVariant& getImpl() const;
  [[nodiscard]] Expr getLeadingOrEndingFactor(bool leading) const;
  [[nodiscard]] Expr factorOut(const Expr& factor, bool leading) const;
  [[nodiscard]] double complexity() const;
  [[nodiscard]] std::set<Expr> getVariables() const;

 private:
  std::shared_ptr<const ExprVariant> impl_;
};

struct NamedNullaryExpr {
  const std::string name;
  explicit NamedNullaryExpr(std::string_view name) : name(std::string(name)) {}
  virtual ~NamedNullaryExpr() = default;
};

struct UnaryExpr {
  const std::shared_ptr<const Expr> child;
  explicit UnaryExpr(std::shared_ptr<Expr> expr) : child(std::move(expr)) {}
  virtual ~UnaryExpr() = default;
};

struct NaryExpr {
  const std::vector<Expr> terms;
  NaryExpr(std::vector<Expr> terms) : terms(std::move(terms)) {}
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

std::strong_ordering operator<=>(const Expr& left, const Expr& right);
bool operator==(const Expr& left, const Expr& right);
std::ostream& operator<<(std::ostream& os, const Expr& expr);
struct ExprHash {
  std::size_t operator()(const Expr& expr) const;
};

namespace ExprFactory {
[[nodiscard]] Expr number(const double value);
[[nodiscard]] Expr namedScalar(std::string_view name);
[[nodiscard]] Expr namedVector(std::string_view name);
[[nodiscard]] Expr variable(std::string_view name);
[[nodiscard]] Expr matrix(std::string_view name);
[[nodiscard]] Expr symmetricMatrix(std::string_view name);
[[nodiscard]] Expr diagonalMatrix(Expr expr);
[[nodiscard]] Expr transpose(Expr expr);
[[nodiscard]] Expr negate(Expr expr);
[[nodiscard]] Expr invert(Expr expr);
[[nodiscard]] Expr log(Expr expr);
[[nodiscard]] Expr sum(std::vector<Expr> terms);
[[nodiscard]] Expr product(std::vector<Expr> terms);
}  // namespace ExprFactory
}  // namespace Expression
