#include "Utils/Helpers.h"

#include <gtest/gtest.h>

#include "Expr.h"
#include "ExprFactory.h"

using namespace Expression;

class HelpersTest : public ::testing::Test {
 protected:
  void SetUp() override {}
};

TEST_F(HelpersTest, ConstantExprs) {
  EXPECT_EQ(unity->to_string(), "1");
  EXPECT_EQ(zero->to_string(), "0");
}

TEST_F(HelpersTest, VariantTypeCheck) {
  EXPECT_TRUE((VariantType<std::variant<int, double>>));
  EXPECT_FALSE(VariantType<int>);

  // Test the VariantType concept
  auto test_variant = [](VariantType auto var) { return true; };
  EXPECT_TRUE(test_variant(std::variant<int, double>{}));
}

TEST_F(HelpersTest, MatchFunction) {
  auto number_expr = ExprFactory::number(42.0);
  auto var_expr = ExprFactory::variable("x");

  // Test matching with Expr
  auto match_result1 = match(number_expr)
                           .with([](const Number& n) { return n.value; },
                                 [](const auto& _) { return 0.0; });
  EXPECT_EQ(match_result1, 42.0);

  auto match_result2 =
      match(var_expr).with([](const Number& n) { return "number"; },
                           [](const Variable& v) { return "variable"; },
                           [](const auto& _) { return "other"; });
  EXPECT_EQ(match_result2, "variable");
}

TEST_F(HelpersTest, VariadicMatchFunction) {
  // Test variadic match with multiple variants
  std::variant<int, double> v1 = 10;
  std::variant<std::string, bool> v2 = "hello";
  std::variant<char, float> v3 = 'A';

  const auto matcher = [&]() {
    return match(v1, v2, v3)
        .with(
            [](int i, const std::string& s, char c) {
              return std::string("int-string-char: ") + s + c +
                     std::to_string(i);
            },
            [](double d, bool b, float f) {
              return std::string("double-bool-float");
            },
            [](auto&&, auto&&, auto&&) {
              return std::string("other combination");
            });
  };

  // Match with three variants simultaneously
  auto result = matcher();
  EXPECT_EQ(result, "int-string-char: helloA10");

  // Test with different variant values
  v1 = 3.14;
  v2 = true;
  v3 = 2.5f;

  auto result2 = matcher();
  EXPECT_EQ(result2, "double-bool-float");

  // Test with other combination
  v3 = 'A';
  auto result3 = matcher();
  EXPECT_EQ(result3, "other combination");
}

TEST_F(HelpersTest, TypeChecking) {
  // Test is_named_nullary
  EXPECT_TRUE((NamedNullaryType<Variable>));
  EXPECT_FALSE((NamedNullaryType<Number>));

  // Test is_unary
  struct MockUnaryExpr : UnaryExpr {
    using UnaryExpr::UnaryExpr;
  };

  EXPECT_TRUE((UnaryType<MockUnaryExpr>));
  EXPECT_FALSE((UnaryType<Number>));

  // Test is_nary
  struct MockNaryExpr : NaryExpr {
    using NaryExpr::NaryExpr;
  };

  EXPECT_TRUE((NaryType<MockNaryExpr>));
  EXPECT_FALSE((NaryType<Number>));
}

TEST_F(HelpersTest, IsFunction) {
  auto number_expr = ExprFactory::number(42.0);
  auto var_expr = ExprFactory::variable("x");

  EXPECT_TRUE(is<Number>(number_expr));
  EXPECT_FALSE(is<Variable>(number_expr));

  EXPECT_TRUE(is<Variable>(var_expr));
  EXPECT_FALSE(is<Number>(var_expr));
}

TEST_F(HelpersTest, Transform) {
  auto terms = std::vector{ExprFactory::number(1.0), ExprFactory::number(2.0),
                           ExprFactory::number(3.0)};

  auto double_value = [](const ExprPtr& e) {
    return match(e).with(
        [](const Number& n) { return ExprFactory::number(n.value * 2.0); },
        [&](const auto&) { return e; });
  };

  auto transformed = transform(terms, double_value);

  EXPECT_EQ(transformed.size(), terms.size());
  EXPECT_EQ(match(transformed[0])
                .with([](const Number& n) { return n.value; },
                      [](const auto& _) { return 0.0; }),
            2.0);
  EXPECT_EQ(match(transformed[1])
                .with([](const Number& n) { return n.value; },
                      [](const auto& _) { return 0.0; }),
            4.0);
  EXPECT_EQ(match(transformed[2])
                .with([](const Number& n) { return n.value; },
                      [](const auto& _) { return 0.0; }),
            6.0);
}

TEST_F(HelpersTest, OverloadedVisitor) {
  std::variant<int, double, std::string> var = 3.14;

  auto result = std::visit(
      overloaded{[](int i) { return "int"; }, [](double d) { return "double"; },
                 [](const std::string& s) { return "string"; }},
      var);

  EXPECT_EQ(result, "double");
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
