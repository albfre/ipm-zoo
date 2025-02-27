#pragma once
#include <memory>
#include <set>
#include <string>
#include <vector>

namespace Expression {
enum class ExprType {
  Number,
  NamedConstant,
  Variable,
  Invert,
  Log,
  Sum,
  Product,
  Negate,
};

class Expr {
 public:
  Expr(ExprType type, std::vector<Expr> terms);
  explicit Expr(ExprType type, const std::string& name);
  explicit Expr(const double value);
  Expr differentiate(const Expr& var) const;
  Expr simplify(bool distribute = true) const;
  std::string toString() const;
  ExprType getType() const;
  std::set<Expr> getVariables() const;
  bool containsSubexpression(const Expr& expr) const;

 private:
  Expr simplify_(bool distribute = true) const;
  const Expr& getSingleChild_() const;
  std::set<Expr> getUniqueFactors_() const;
  Expr factorOut(const Expr& factor) const;
  double complexity_() const;
  ExprType type_;
  std::string name_ = "";
  double value_ = 0.0;
  std::vector<Expr> terms_ = {};
};

std::strong_ordering operator<=>(const Expr& left, const Expr& right);
bool operator==(const Expr& left, const Expr& right);

namespace ExprFactory {
Expr number(const double value);
Expr namedConstant(const std::string& name);
Expr variable(const std::string& name);
Expr negate(Expr expr);
Expr invert(Expr expr);
Expr log(Expr expr);
Expr product(std::vector<Expr> terms);
Expr sum(std::vector<Expr> terms);
}  // namespace ExprFactory
}  // namespace Expression
