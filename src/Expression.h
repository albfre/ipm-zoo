#pragma once
#include <memory>
#include <set>
#include <string>
#include <variant>
#include <vector>

namespace Expression {
class Expr;

struct NamedNullaryExpr {
  const std::string name;
};

struct UnaryExpr {
  const std::shared_ptr<const Expr> child;
};

struct NaryExpr {
  const std::vector<Expr> terms;
};

struct Number {
  const double value;
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

  Expr differentiate(const Expr& var) const;
  Expr simplify(bool distribute = true) const;
  Expr simplifyOnce(bool distribute = true) const;
  std::string toString(bool condensed = false) const;
  std::string toExpressionString() const;
  bool containsSubexpression(const Expr& expr) const;
  Expr replaceSubexpression(const Expr& expr, const Expr& replacement) const;
  const ExprVariant& getImpl() const { return *impl_; }
  Expr getLeadingOrEndingFactor(bool leading) const;
  Expr factorOut(const Expr& factor, bool leading) const;
  double complexity() const;
  std::set<Expr> getVariables() const;

 private:
  std::shared_ptr<ExprVariant> impl_;
};

std::strong_ordering operator<=>(const Expr& left, const Expr& right);
bool operator==(const Expr& left, const Expr& right);
std::ostream& operator<<(std::ostream& os, const Expr& expr);
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
