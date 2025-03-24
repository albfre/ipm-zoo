#include "Expr.h"

#include <gtest/gtest.h>

#include "ExprFactory.h"
#include "Helpers.h"

using namespace Expression;

class ExpressionTest : public ::testing::Test {
 protected:
  void SetUp() override {}

  ExprPtr x = ExprFactory::variable("x");
  ExprPtr y = ExprFactory::variable("y");
  ExprPtr z = ExprFactory::variable("z");
  ExprPtr a = ExprFactory::namedVector("a");
  ExprPtr b = ExprFactory::namedVector("b");
  ExprPtr c = ExprFactory::namedVector("c");
  ExprPtr zero = ExprFactory::number(0.0);
  ExprPtr one = ExprFactory::number(1.0);
  ExprPtr two = ExprFactory::number(2.0);
  ExprPtr three = ExprFactory::number(3.0);
};

// Test construction and basic operations
TEST_F(ExpressionTest, BasicConstruction) {
  // Test number
  auto num = ExprFactory::number(42.0);
  EXPECT_EQ(num->toString(), "42");
  EXPECT_TRUE(is<Number>(num));

  // Test variable
  EXPECT_EQ(x->toString(), "x");
  EXPECT_TRUE(is<Variable>(x));

  // Test named constant
  EXPECT_EQ(a->toString(), "a");
  EXPECT_TRUE(is<NamedVector>(a));
}

// Test basic arithmetic operations
TEST_F(ExpressionTest, BasicArithmetic) {
  // Test summation
  auto sum = ExprFactory::sum({x, y});
  EXPECT_EQ(sum->toString(), "(x + y)");

  // Test product
  auto product = ExprFactory::product({x, y});
  EXPECT_EQ(product->toString(), "(x * y)");

  // Test negation
  auto neg = ExprFactory::negate(x);
  EXPECT_EQ(neg->toString(), "-x");

  // Test inversion
  auto inv = ExprFactory::invert(x);
  EXPECT_EQ(inv->toString(), "x^{-1}");

  // Test logarithm
  auto log = ExprFactory::log(x);
  EXPECT_EQ(log->toString(), "\\log(x)");
}

// Test expression comparison
TEST_F(ExpressionTest, Comparison) {
  // Create identical expressions
  auto expr1 = ExprFactory::sum({x, y});
  auto expr2 = ExprFactory::sum({x, y});

  // Create different expressions
  auto expr3 = ExprFactory::sum({y, x});  // Same variables, different order
  auto expr4 = ExprFactory::sum({x, z});  // Different variable

  EXPECT_EQ(expr1, expr2);  // Should be equal
  EXPECT_NE(expr1, expr3);  // Should be different due to ordering
  EXPECT_NE(expr1, expr4);  // Should be different due to variables
}

// Test differentiation
TEST_F(ExpressionTest, DifferentiateConstant) {
  // Test constant differentiation
  EXPECT_EQ(one->differentiate(x)->toString(), "0");
}

TEST_F(ExpressionTest, DifferentiateVariable) {
  // Test variable differentiation
  EXPECT_EQ(x->differentiate(x)->toString(), "1");
  EXPECT_EQ(y->differentiate(x)->toString(), "0");
}

TEST_F(ExpressionTest, DifferentiateSum) {
  // Test sum differentiation
  auto sum = ExprFactory::sum({x, y});
  EXPECT_EQ(sum->differentiate(x)->toString(), "(1 + 0)");
}

TEST_F(ExpressionTest, DifferentiateProduct) {
  // Test product differentiation
  auto prod = ExprFactory::product({x, y});
  auto diff_prod = prod->differentiate(x);
  EXPECT_EQ(diff_prod->toString(), "((1 * y) + (x * 0))");
}

TEST_F(ExpressionTest, DifferentiateNegation) {
  // Test negate differentiation
  auto neg = ExprFactory::negate(x);
  EXPECT_EQ(neg->differentiate(x)->toString(), "-1");
}

TEST_F(ExpressionTest, DifferentiateTranspose) {
  // Test transpose differentiation
  auto t1 = ExprFactory::product({ExprFactory::transpose(x), a});
  auto t2 = ExprFactory::product({ExprFactory::transpose(a), x});
  EXPECT_EQ(t1->differentiate(x)->simplify()->toString(), "a");
  EXPECT_EQ(t2->differentiate(x)->simplify()->toString(), "a");
}

TEST_F(ExpressionTest, DifferentiateXSquared) {
  // Test more complex expressions
  auto expr = ExprFactory::product({x, x});  // x^2
  EXPECT_EQ(expr->differentiate(x)->simplify()->toString(), "(2 * x)");
}

TEST_F(ExpressionTest, DifferentiateQuadraticForm) {
  // Test quadratic form differentiation
  auto t1 = ExprFactory::product(
      {a, ExprFactory::transpose(x), ExprFactory::variable("Q"), x});
  auto t2 = ExprFactory::product(
      {a, ExprFactory::transpose(x), ExprFactory::symmetricMatrix("Q"), x});
  auto t3 =
      ExprFactory::product({ExprFactory::number(0.5), ExprFactory::transpose(x),
                            ExprFactory::symmetricMatrix("Q"), x});
  EXPECT_EQ(t1->differentiate(x)->simplify()->toString(),
            "(a * (Q + Q^T) * x)");
  EXPECT_EQ(t2->differentiate(x)->simplify()->toString(), "(2 * a * Q * x)");
  EXPECT_EQ(t3->differentiate(x)->simplify()->toString(), "(Q * x)");
}

// Test simplification
TEST_F(ExpressionTest, SimplifySumOfNumbers) {
  // Test basic numeric simplification
  auto numSum = ExprFactory::sum({one, two});
  EXPECT_EQ(numSum->simplify()->toString(), "3");
}

TEST_F(ExpressionTest, SimplifySumWithZero) {
  // Test basic numeric simplification
  auto numSum = ExprFactory::sum({x, zero});
  EXPECT_EQ(numSum->simplify()->toString(), "x");
}

TEST_F(ExpressionTest, SimplifyProductOfNumbers) {
  // Test basic product simplification
  auto numProd = ExprFactory::product({two, three});
  EXPECT_EQ(numProd->simplify()->toString(), "6");
}

TEST_F(ExpressionTest, SimplifySumAssociatively) {
  // Test associative property in sum
  auto nestedSum = ExprFactory::sum({x, ExprFactory::sum({y, z})});
  EXPECT_EQ(nestedSum->simplify()->toString(), "(x + y + z)");
}

TEST_F(ExpressionTest, SimplifySumCommutatively) {
  // Test commutative property in sum
  auto zyx = ExprFactory::sum({z, y, x});
  EXPECT_EQ(zyx->simplify()->toString(), "(x + y + z)");
}

TEST_F(ExpressionTest, SimplifySumDistributively) {
  // Test distributive property in sum
  auto sum = ExprFactory::sum({x, y, ExprFactory::product({three, x})});
  EXPECT_EQ(sum->simplify()->toString(), "(y + (4 * x))");
}

TEST_F(ExpressionTest, SimplifyProductAssociatively) {
  // Test associative property in product
  auto nestedProd = ExprFactory::product({x, ExprFactory::product({y, z})});
  EXPECT_EQ(nestedProd->simplify()->toString(), "(x * y * z)");
}

TEST_F(ExpressionTest, SimplifyDoubleNegation) {
  // Test double negation
  auto doubleNeg = ExprFactory::negate(ExprFactory::negate(x));
  EXPECT_EQ(doubleNeg->simplify()->toString(), "x");
}

TEST_F(ExpressionTest, SimplifyNegationOfProductWithNegation) {
  // Test negation of product with negation
  auto neg =
      ExprFactory::negate(ExprFactory::product({x, y, ExprFactory::negate(z)}));
  EXPECT_EQ(neg->simplify()->toString(), "(x * y * z)");
}

TEST_F(ExpressionTest, SimplifyNegationOfSumWithNegations) {
  // Test negation of sum with negations
  auto neg = ExprFactory::negate(
      ExprFactory::sum({ExprFactory::negate(x), ExprFactory::negate(y), z}));
  EXPECT_EQ(neg->simplify()->toString(), "(x + y - z)");
}

TEST_F(ExpressionTest, SimplifyProductWithNegation) {
  // Test product with negation
  auto prod = ExprFactory::product({x, y, ExprFactory::negate(z)});
  std::cout << prod->toExpressionString() << std::endl;
  std::cout << "simp: " << prod->simplify()->toExpressionString() << std::endl;
  EXPECT_EQ(prod->simplify()->toString(), "-(x * y * z)");
}

TEST_F(ExpressionTest, SimplifyExtractsFactor) {
  // Test algebraic simplification
  auto expr = ExprFactory::sum(
      {ExprFactory::product({x, y}), ExprFactory::product({x, z})});
  EXPECT_EQ(expr->simplify()->toString(), "(x * (y + z))");
}

TEST_F(ExpressionTest, SimplifyDistributesFactor) {
  // Test algebraic simplification
  auto expr =
      ExprFactory::product({x, ExprFactory::sum({ExprFactory::invert(x), z})});
  EXPECT_EQ(expr->simplify()->toString(), "(1 + (x * z))");
}

TEST_F(ExpressionTest, SimplifyRemovesCancelingTermsForSum) {
  // Test cancellation
  auto expr = ExprFactory::sum({x, ExprFactory::negate(x)});
  EXPECT_EQ(expr->simplify()->toString(), "0");
}

TEST_F(ExpressionTest, SimplifyRemovesCancelingTermsForProduct) {
  // Test cancellation
  auto expr = ExprFactory::product({x, ExprFactory::invert(x)});
  EXPECT_EQ(expr->simplify()->toString(), "1");
}

TEST_F(ExpressionTest, SimplifyInvertedProduct) {
  // Test cancellation
  auto expr =
      ExprFactory::invert(ExprFactory::product({x, ExprFactory::invert(y)}));
  EXPECT_EQ(expr->simplify()->toString(), "(y * x^{-1})");
}

TEST_F(ExpressionTest, SimplifyMultiplicationWithZero) {
  // Test multiplication with zero
  auto prodWithZero = ExprFactory::product({x, zero});
  EXPECT_EQ(prodWithZero->simplify()->toString(), "0");
}

TEST_F(ExpressionTest, SimplifyMultiplicationWithOnes) {
  // Test multiplication with ones
  auto prodWithOne = ExprFactory::product({one, one, x, one, one});
  auto prodOfOnes = ExprFactory::product({one, one, one, one});
  EXPECT_EQ(prodWithOne->simplify()->toString(), "x");
  EXPECT_EQ(prodOfOnes->simplify()->toString(), "1");
}

TEST_F(ExpressionTest, SimplifyComplex) {
  // Test more complex expression
  auto complexExpr = ExprFactory::sum(
      {ExprFactory::product({x, y}), ExprFactory::product({x, y}),
       ExprFactory::product({a, b}),
       ExprFactory::product({a, ExprFactory::invert(a)}),
       ExprFactory::negate(ExprFactory::product({a, b})), one, two,
       ExprFactory::product({three, x}), x, x});
  // x * y + x * y + a * b + a * inv(a) - a * b + 1 + 2 + 3 * x + x + x
  // = 2 * x * y + 1 + 3 + 5 * x
  // = 4 + (2 * x * y) + (5 * x))
  auto simplified = complexExpr->simplify();
  EXPECT_EQ(simplified->toString(), "(4 + (2 * x * y) + (5 * x))");
}

TEST_F(ExpressionTest, SimplifyRhsExpression) {
  // Test expression appearing in augmented system RHS
  using EF = Expression::ExprFactory;
  auto expr =
      EF::product({EF::invert(EF::diagonalMatrix(EF::variable("z"))),
                   EF::diagonalMatrix(EF::variable("\\lambda_{z}")),
                   EF::sum({EF::product({EF::invert(EF::diagonalMatrix(
                                             EF::variable("\\lambda_{z}"))),
                                         EF::namedVector("r_{z}")}),
                            EF::negate(EF::namedVector("r_{\\lambda_{z}}"))})});
  auto expr2 =
      EF::sum({EF::product({EF::invert(EF::diagonalMatrix(EF::variable("z"))),
                            EF::namedVector("r_{z}")}),
               EF::product({EF::invert(EF::diagonalMatrix(EF::variable("z"))),
                            EF::diagonalMatrix(EF::variable("\\lambda_{z}")),
                            EF::negate(EF::namedVector("r_{\\lambda_{z}}"))})});
  auto simplified = expr->simplify();
  auto simplified2 = expr2->simplify();
  std::cout << simplified->toString() << std::endl;
  EXPECT_EQ(simplified->toString(), simplified2->toString());
}

TEST_F(ExpressionTest, SimplifyComplex2) {
  // Test more complex expression
  auto complexExpr = ExprFactory::sum(
      {ExprFactory::product({a, x}),
       ExprFactory::sum(
           {ExprFactory::product({b, x}),
            ExprFactory::product({c, ExprFactory::product({x, c})})})});
  // a * x + b * x + c * (x * c)
  // = ((c * x * c) + ((a + b) * x))
  auto simplified = complexExpr->simplify();
  EXPECT_EQ(simplified->toString(), "((c * x * c) + ((a + b) * x))");
}

// Test variable extraction
TEST_F(ExpressionTest, VariableExtraction) {
  auto expr = ExprFactory::product({a, ExprFactory::sum({x, y, z}), b});

  auto vars = expr->getVariables();
  EXPECT_EQ(vars.size(), 3);
  EXPECT_TRUE(vars.find(x) != vars.end());
  EXPECT_TRUE(vars.find(y) != vars.end());
  EXPECT_TRUE(vars.find(z) != vars.end());
}

// Test the Lagrangian from ipmZoo.cpp
TEST_F(ExpressionTest, Lagrangian) {
  using EF = Expression::ExprFactory;
  auto x = EF::variable("x");
  auto H = EF::symmetricMatrix("H");
  auto c = EF::namedVector("c");
  auto A = EF::matrix("A");
  auto mu = EF::namedScalar("\\mu");
  auto e = EF::namedVector("e");
  auto g = EF::variable("g");
  auto t = EF::variable("t");
  auto y = EF::variable("y");
  auto z = EF::variable("z");
  auto s = EF::variable("s");
  auto lambda_A = EF::variable("\\lambda_A");
  auto lambda_g = EF::variable("\\lambda_g");
  auto lambda_t = EF::variable("\\lambda_t");
  auto lambda_y = EF::variable("\\lambda_y");
  auto lambda_z = EF::variable("\\lambda_z");
  auto l_A = EF::namedVector("l_A");
  auto u_A = EF::namedVector("u_A");
  auto l_x = EF::namedVector("l_x");
  auto u_x = EF::namedVector("u_x");

  auto xHx = EF::product({EF::number(0.5), x, H, x});
  auto cx = EF::product({c, x});
  auto Ax = EF::product({A, x});
  auto minusMuLogG = EF::negate(EF::product({mu, e, EF::log({g})}));
  auto minusMuLogT = EF::negate(EF::product({mu, e, EF::log({t})}));
  auto minusMuLogY = EF::negate(EF::product({mu, e, EF::log({y})}));
  auto minusMuLogZ = EF::negate(EF::product({mu, e, EF::log({z})}));
  auto lambda_A_term = EF::product({lambda_A, EF::sum({Ax, EF::negate(s)})});
  auto lambda_g_term = EF::negate(
      EF::product({lambda_g, EF::sum({s, EF::negate(g), EF::negate(l_A)})}));
  auto lambda_t_term =
      EF::product({lambda_t, EF::sum({s, t, EF::negate(u_A)})});
  auto lambda_y_term = EF::negate(
      EF::product({lambda_y, EF::sum({x, EF::negate(y), EF::negate(l_x)})}));
  auto lambda_z_term =
      EF::product({lambda_z, EF::sum({x, z, EF::negate(u_x)})});

  auto lagrangian = EF::sum({xHx, cx, minusMuLogG, minusMuLogT, minusMuLogY,
                             minusMuLogZ, lambda_A_term, lambda_g_term,
                             lambda_t_term, lambda_y_term, lambda_z_term});

  // Test that the Lagrangian can be differentiated and simplified
  auto dgL = lagrangian->differentiate(g)->simplify();
  auto dxL = lagrangian->differentiate(x)->simplify();

  // Just verify they don't crash and return something of the expected type
  EXPECT_NE(dgL->toString(), "");
  EXPECT_NE(dxL->toString(), "");

  // Verify that all variables are detected
  auto vars = lagrangian->getVariables();
  EXPECT_TRUE(vars.find(g) != vars.end());
  EXPECT_TRUE(vars.find(x) != vars.end());
  EXPECT_TRUE(vars.find(t) != vars.end());
  EXPECT_TRUE(vars.find(y) != vars.end());
  EXPECT_TRUE(vars.find(z) != vars.end());
  EXPECT_TRUE(vars.find(s) != vars.end());
  EXPECT_TRUE(vars.find(lambda_A) != vars.end());
  EXPECT_TRUE(vars.find(lambda_g) != vars.end());
  EXPECT_TRUE(vars.find(lambda_t) != vars.end());
  EXPECT_TRUE(vars.find(lambda_y) != vars.end());
  EXPECT_TRUE(vars.find(lambda_z) != vars.end());
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}