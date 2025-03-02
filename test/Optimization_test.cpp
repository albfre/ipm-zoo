
#include "Optimization.h"

#include <gtest/gtest.h>

using namespace Optimization;

class OptimizationTest : public ::testing::Test {
 protected:
  void SetUp() override {}
};

TEST_F(OptimizationTest, GetLagrangian) {
  auto settings = Settings();
  auto names = VariableNames();
  auto [lagrangian, variables] = getLagrangian(names, settings);
  const auto str = lagrangian.toString();
  EXPECT_NE(str.find(names.s_A), std::string::npos);
  EXPECT_NE(str.find(names.s_Al), std::string::npos);
  EXPECT_NE(str.find(names.s_Au), std::string::npos);
  EXPECT_NE(str.find(names.s_xl), std::string::npos);
  EXPECT_NE(str.find(names.s_xu), std::string::npos);
}

TEST_F(OptimizationTest, GetLagrangianLower) {
  auto settings = Settings();
  settings.inequalities = Bounds::Lower;
  settings.variableBounds = Bounds::Lower;
  auto names = VariableNames();
  names.s_Au = "t123";
  auto [lagrangian, variables] = getLagrangian(names, settings);
  const auto str = lagrangian.toString();
  EXPECT_NE(str.find(names.s_A), std::string::npos);
  EXPECT_NE(str.find(names.s_Al), std::string::npos);
  EXPECT_EQ(str.find(names.s_Au), std::string::npos);
  EXPECT_NE(str.find(names.s_xl), std::string::npos);
  EXPECT_EQ(str.find(names.s_xu), std::string::npos);
}

TEST_F(OptimizationTest, GetLagrangianUpper) {
  auto settings = Settings();
  settings.inequalities = Bounds::Upper;
  settings.variableBounds = Bounds::Upper;
  auto names = VariableNames();
  names.s_Al = "g123";
  auto [lagrangian, variables] = getLagrangian(names, settings);
  const auto str = lagrangian.toString();
  EXPECT_NE(str.find(names.s_A), std::string::npos);
  EXPECT_EQ(str.find(names.s_Al), std::string::npos);
  EXPECT_NE(str.find(names.s_Au), std::string::npos);
  EXPECT_EQ(str.find(names.s_xl), std::string::npos);
  EXPECT_NE(str.find(names.s_xu), std::string::npos);
}

TEST_F(OptimizationTest, GetLagrangianNoBounds) {
  auto settings = Settings();
  settings.inequalities = Bounds::None;
  settings.variableBounds = Bounds::None;
  auto names = VariableNames();
  names.s_Al = "g123";
  names.s_Au = "t123";
  names.s_xl = "y123";
  names.s_xu = "z123";
  auto [lagrangian, variables] = getLagrangian(names, settings);
  const auto str = lagrangian.toString();
  EXPECT_EQ(str.find(names.s_A), std::string::npos);
  EXPECT_EQ(str.find(names.s_Al), std::string::npos);
  EXPECT_EQ(str.find(names.s_Au), std::string::npos);
  EXPECT_EQ(str.find(names.s_xl), std::string::npos);
  EXPECT_EQ(str.find(names.s_xu), std::string::npos);
}

TEST_F(OptimizationTest, GetFirstOrderOptimalityConditions) {
  auto settings = Settings();
  auto names = VariableNames();
  auto [lagrangian, variables] = getLagrangian(names, settings);
  auto firstOrder = getFirstOrderOptimalityConditions(lagrangian, variables);
  EXPECT_EQ(firstOrder.size(), variables.size());
}

TEST_F(OptimizationTest, GetNewtonSystem) {
  auto settings = Settings();
  auto names = VariableNames();
  auto [lagrangian, variables] = getLagrangian(names, settings);
  auto [lhs, rhs] = getNewtonSystem(lagrangian, variables);
  EXPECT_EQ(lhs.size(), rhs.size());
  EXPECT_EQ(lhs.size(), variables.size());
  for (const auto& row : lhs) {
    EXPECT_EQ(row.size(), variables.size());
  }
}

TEST_F(OptimizationTest, GetShorthandRhs) {
  auto settings = Settings();
  auto names = VariableNames();
  auto [lagrangian, variables] = getLagrangian(names, settings);
  auto rhs = getShorthandRhs(variables);
  EXPECT_EQ(variables.size(), rhs.size());
  for (size_t i = 0; i < variables.size(); ++i) {
    EXPECT_EQ("-r_{" + variables[i].toString() + "}", rhs[i].toString());
  }
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
