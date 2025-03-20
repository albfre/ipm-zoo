#pragma once
#include <memory>
#include <set>
#include <string>
#include <variant>
#include <vector>

namespace Expression {
class Expr;

struct NamedNullaryExpr {
  std::string name;
};

struct UnaryExpr {
  std::unique_ptr<Expr> child;
};

struct NaryExpr {
  std::vector<Expr> terms;
};

struct Number {
  double value;
};
struct NamedScalar : public NamedNullaryExpr {};
struct NamedVector : public NamedNullaryExpr {};
struct Variable : public NamedNullaryExpr {};
struct Matrix : public NamedNullaryExpr {};
struct SymmetricMatrix : public NamedNullaryExpr {};
struct DiagonalMatrix : public UnaryExpr {};
struct Transpose : public UnaryExpr {};
struct Invert : public UnaryExpr {};
struct Log : public UnaryExpr {};
struct Negate : public UnaryExpr {};
struct Sum : public NaryExpr {};
struct Product : public NaryExpr {};

class Expr {
  using ExprVariant =
      std::variant<Number, NamedScalar, NamedVector, Variable, Matrix,
                   SymmetricMatrix, DiagonalMatrix, Transpose, Invert, Log, Sum,
                   Product, Negate>;

 public:
  template <typename T>
    requires std::is_base_of_v<UnaryExpr, std::decay_t<T>>
  explicit Expr(const T& value)
      : impl_(T{std::make_unique<Expr>(*value.child)}) {}

  template <typename T>
    requires(!std::is_same_v<Expr, std::decay_t<T>>)
  explicit Expr(T&& value) : impl_(std::forward<T>(value)) {}

  Expr(const Expr& other);
  Expr(Expr&& other) noexcept : impl_(std::move(other.impl_)) {}

  Expr& operator=(const Expr& other);

  Expr differentiate(const Expr& var) const;
  Expr simplify(bool distribute = true) const;
  Expr simplifyOnce(bool distribute = true) const;
  std::string toString(bool condensed = false) const;
  std::string toExpressionString() const;
  bool containsSubexpression(const Expr& expr) const;
  Expr replaceSubexpression(const Expr& expr, const Expr& replacement) const;
  const ExprVariant& getImpl() const { return impl_; }
  Expr getLeadingOrEndingFactor(bool leading) const;
  Expr factorOut(const Expr& factor, bool leading) const;
  double complexity() const;
  std::set<Expr> getVariables() const;

 private:
  ExprVariant impl_;
};

std::strong_ordering operator<=>(const Expr& left, const Expr& right);
bool operator==(const Expr& left, const Expr& right);
struct ExprHash {
  std::size_t operator()(const Expr& expr) const;
};

namespace ExprFactory {
Expr number(const double value);
Expr namedScalar(const std::string& name);
Expr namedVector(const std::string& name);
Expr variable(const std::string& name);
Expr matrix(const std::string& name);
Expr symmetricMatrix(const std::string& name);
Expr diagonalMatrix(Expr expr);
Expr transpose(Expr expr);
Expr negate(Expr expr);
Expr invert(Expr expr);
Expr log(Expr expr);
Expr sum(std::vector<Expr> terms);
Expr product(std::vector<Expr> terms);
}  // namespace ExprFactory
}  // namespace Expression
