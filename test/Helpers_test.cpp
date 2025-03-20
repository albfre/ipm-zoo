#include "Helpers.h"

#include <gtest/gtest.h>

#include "Expression.h"

using namespace Expression;

class HelpersTest : public ::testing::Test {
 protected:
  void SetUp() override {}
};

TEST_F(HelpersTest, ConstantExprs) {
  EXPECT_EQ(unity.toString(), "1");
  EXPECT_EQ(zero.toString(), "0");
}

TEST_F(HelpersTest, VariantTypeCheck) {
  EXPECT_TRUE((is_variant_v<std::variant<int, double>>));
  EXPECT_FALSE(is_variant_v<int>);

  // Test the VariantType concept
  auto testVariant = [](VariantType auto var) { return true; };
  EXPECT_TRUE(testVariant(std::variant<int, double>{}));
}

TEST_F(HelpersTest, MatchFunction) {
  auto numberExpr = ExprFactory::number(42.0);
  auto varExpr = ExprFactory::variable("x");

  // Test matching with Expr
  auto matchResult1 = match(
      numberExpr, [](const Number& n) { return n.value; },
      [](const auto& _) { return 0.0; });
  EXPECT_EQ(matchResult1, 42.0);

  auto matchResult2 = match(
      varExpr, [](const Number& n) { return "number"; },
      [](const Variable& v) { return "variable"; },
      [](const auto& _) { return "other"; });
  EXPECT_EQ(matchResult2, "variable");
}

TEST_F(HelpersTest, TypeChecking) {
  // Test is_named_nullary
  EXPECT_TRUE((is_named_nullary_v<Variable>));
  EXPECT_FALSE((is_named_nullary_v<Number>));

  // Test is_unary
  struct MockUnaryExpr : UnaryExpr {
    using UnaryExpr::UnaryExpr;
  };

  EXPECT_TRUE((is_unary_v<MockUnaryExpr>));
  EXPECT_FALSE((is_unary_v<Number>));

  // Test is_nary
  struct MockNaryExpr : NaryExpr {
    using NaryExpr::NaryExpr;
  };

  EXPECT_TRUE((is_nary_v<MockNaryExpr>));
  EXPECT_FALSE((is_nary_v<Number>));
}

TEST_F(HelpersTest, IsFunction) {
  auto numberExpr = ExprFactory::number(42.0);
  auto varExpr = ExprFactory::variable("x");

  EXPECT_TRUE(is<Number>(numberExpr));
  EXPECT_FALSE(is<Variable>(numberExpr));

  EXPECT_TRUE(is<Variable>(varExpr));
  EXPECT_FALSE(is<Number>(varExpr));
}

TEST_F(HelpersTest, Transform) {
  std::vector<Expr> terms = {ExprFactory::number(1.0), ExprFactory::number(2.0),
                             ExprFactory::number(3.0)};

  auto doubleValue = [](const Expr& e) {
    return match(
        e, [](const Number& n) { return ExprFactory::number(n.value * 2.0); },
        [&](const auto&) { return e; });
  };

  auto transformed = transform(terms, doubleValue);

  EXPECT_EQ(transformed.size(), terms.size());
  EXPECT_EQ(match(
                transformed[0], [](const Number& n) { return n.value; },
                [](const auto& _) { return 0.0; }),
            2.0);
  EXPECT_EQ(match(
                transformed[1], [](const Number& n) { return n.value; },
                [](const auto& _) { return 0.0; }),
            4.0);
  EXPECT_EQ(match(
                transformed[2], [](const Number& n) { return n.value; },
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
