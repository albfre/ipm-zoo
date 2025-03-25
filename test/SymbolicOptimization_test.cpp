#include "SymbolicOptimization.h"

#include <gtest/gtest.h>

#include "ExprFactory.h"

using namespace SymbolicOptimization;

class OptimizationTest : public ::testing::Test {
 protected:
  void SetUp() override {}
};

TEST_F(OptimizationTest, GetLagrangian) {
  auto settings = Settings();
  auto names = VariableNames();
  auto [lagrangian, variables] = getLagrangian(settings, names);
  const auto str = lagrangian->toString();
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
  auto [lagrangian, variables] = getLagrangian(settings, names);
  const auto str = lagrangian->toString();
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
  auto [lagrangian, variables] = getLagrangian(settings, names);
  const auto str = lagrangian->toString();
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
  auto [lagrangian, variables] = getLagrangian(settings, names);
  const auto str = lagrangian->toString();
  EXPECT_EQ(str.find(names.s_A), std::string::npos);
  EXPECT_EQ(str.find(names.s_Al), std::string::npos);
  EXPECT_EQ(str.find(names.s_Au), std::string::npos);
  EXPECT_EQ(str.find(names.s_xl), std::string::npos);
  EXPECT_EQ(str.find(names.s_xu), std::string::npos);
}

TEST_F(OptimizationTest, GetFirstOrderOptimalityConditions) {
  auto settings = Settings();
  auto names = VariableNames();
  auto [firstOrder, variables] =
      getFirstOrderOptimalityConditions(settings, names);
  EXPECT_EQ(firstOrder.size(), variables.size());
}

TEST_F(OptimizationTest, GetNewtonSystem) {
  auto settings = Settings();
  auto names = VariableNames();
  auto [lhs, rhs, variables, _] = getNewtonSystem(settings, names);
  EXPECT_EQ(lhs.size(), rhs.size());
  EXPECT_EQ(lhs.size(), variables.size());
  for (const auto& row : lhs) {
    EXPECT_EQ(row.size(), variables.size());
  }
}

TEST_F(OptimizationTest, GetShorthandRhs) {
  auto settings = Settings();
  auto names = VariableNames();
  auto [_, variables] = getLagrangian(settings, names);
  auto rhs = getShorthandRhs(variables);
  EXPECT_EQ(variables.size(), rhs.size());
  for (size_t i = 0; i < variables.size(); ++i) {
    EXPECT_EQ("-r_{" + variables[i]->toString() + "}", rhs[i]->toString());
  }
}

TEST_F(OptimizationTest, GaussianElimination) {
  auto settings = Settings();
  auto names = VariableNames();
  auto [lhs, rhs, variables, _] = getNewtonSystem(settings, names);

  // Save original size
  size_t originalSize = lhs.size();

  // Perform one elimination step
  gaussianElimination(lhs, rhs, lhs.size() - 1);

  // Check sizes after elimination
  EXPECT_EQ(lhs.size(), originalSize - 1);
  EXPECT_EQ(rhs.size(), originalSize - 1);
  for (const auto& row : lhs) {
    EXPECT_EQ(row.size(), originalSize - 1);
  }
}

TEST_F(OptimizationTest, GetAugmentedSystem) {
  auto settings = Settings();
  settings.inequalities = Bounds::Lower;
  auto names = VariableNames();
  auto newtonSystem = getNewtonSystem(settings, names);

  // Save original size
  size_t originalSize = newtonSystem.lhs.size();
  EXPECT_GT(originalSize, 2);  // Make sure we have something to eliminate

  // Get augmented system
  auto augmentedSystem = getAugmentedSystem(std::move(newtonSystem));

  // Check that variables were eliminated and delta definitions were created
  EXPECT_LT(augmentedSystem.lhs.size(), originalSize);
  EXPECT_GT(augmentedSystem.deltaDefinitions.size(), 0);

  // Check that first delta definition references the appropriate variable
  const auto& [deltaVar, _] = augmentedSystem.deltaDefinitions.front();
  EXPECT_TRUE(deltaVar->toString().find("\\Delta") != std::string::npos);
}

TEST_F(OptimizationTest, GetNormalEquations) {
  auto settings = Settings();
  settings.inequalities = Bounds::Lower;
  auto names = VariableNames();
  auto newtonSystem = getNewtonSystem(settings, names);

  // Save original size
  size_t originalSize = newtonSystem.lhs.size();
  EXPECT_GT(originalSize, 2);  // Make sure we have something to eliminate

  // Get normal equations (should be a single equation)
  auto normalEquations = getNormalEquations(std::move(newtonSystem));

  // Check that we have only one equation left
  EXPECT_EQ(normalEquations.lhs.size(), 1);
  EXPECT_EQ(normalEquations.rhs.size(), 1);
  EXPECT_EQ(normalEquations.variables.size(), 1);

  // Check that we have delta definitions
  EXPECT_EQ(normalEquations.deltaDefinitions.size(), originalSize - 1);
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}