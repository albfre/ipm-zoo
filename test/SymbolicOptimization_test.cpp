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
  auto [lagrangian, variables] = get_lagrangian(settings, names);
  const auto str = lagrangian->to_string();
  EXPECT_NE(str.find(names.s_A), std::string::npos);
  EXPECT_NE(str.find(names.s_Al), std::string::npos);
  EXPECT_NE(str.find(names.s_Au), std::string::npos);
  EXPECT_NE(str.find(names.s_xl), std::string::npos);
  EXPECT_NE(str.find(names.s_xu), std::string::npos);
}

TEST_F(OptimizationTest, GetLagrangianLower) {
  auto settings = Settings();
  settings.inequalities = Bounds::Lower;
  settings.variable_bounds = Bounds::Lower;
  auto names = VariableNames();
  names.s_Au = "t123";
  auto [lagrangian, variables] = get_lagrangian(settings, names);
  const auto str = lagrangian->to_string();
  EXPECT_NE(str.find(names.s_A), std::string::npos);
  EXPECT_NE(str.find(names.s_Al), std::string::npos);
  EXPECT_EQ(str.find(names.s_Au), std::string::npos);
  EXPECT_NE(str.find(names.s_xl), std::string::npos);
  EXPECT_EQ(str.find(names.s_xu), std::string::npos);
}

TEST_F(OptimizationTest, GetLagrangianUpper) {
  auto settings = Settings();
  settings.inequalities = Bounds::Upper;
  settings.variable_bounds = Bounds::Upper;
  auto names = VariableNames();
  names.s_Al = "g123";
  auto [lagrangian, variables] = get_lagrangian(settings, names);
  const auto str = lagrangian->to_string();
  EXPECT_NE(str.find(names.s_A), std::string::npos);
  EXPECT_EQ(str.find(names.s_Al), std::string::npos);
  EXPECT_NE(str.find(names.s_Au), std::string::npos);
  EXPECT_EQ(str.find(names.s_xl), std::string::npos);
  EXPECT_NE(str.find(names.s_xu), std::string::npos);
}

TEST_F(OptimizationTest, GetLagrangianNoBounds) {
  auto settings = Settings();
  settings.inequalities = Bounds::None;
  settings.variable_bounds = Bounds::None;
  auto names = VariableNames();
  names.s_Al = "g123";
  names.s_Au = "t123";
  names.s_xl = "y123";
  names.s_xu = "z123";
  auto [lagrangian, variables] = get_lagrangian(settings, names);
  const auto str = lagrangian->to_string();
  EXPECT_EQ(str.find(names.s_A), std::string::npos);
  EXPECT_EQ(str.find(names.s_Al), std::string::npos);
  EXPECT_EQ(str.find(names.s_Au), std::string::npos);
  EXPECT_EQ(str.find(names.s_xl), std::string::npos);
  EXPECT_EQ(str.find(names.s_xu), std::string::npos);
}

TEST_F(OptimizationTest, GetFirstOrderOptimalityConditions) {
  auto settings = Settings();
  auto names = VariableNames();
  auto [first_order, variables] =
      get_first_order_optimality_conditions(settings, names);
  EXPECT_EQ(first_order.size(), variables.size());
}

TEST_F(OptimizationTest, GetNewtonSystem) {
  auto settings = Settings();
  auto names = VariableNames();
  auto [lhs, rhs, variables, _] = get_newton_system(settings, names);
  EXPECT_EQ(lhs.size(), rhs.size());
  EXPECT_EQ(lhs.size(), variables.size());
  for (const auto& row : lhs) {
    EXPECT_EQ(row.size(), variables.size());
  }
}

TEST_F(OptimizationTest, GetShorthandRhs) {
  auto settings = Settings();
  auto names = VariableNames();
  auto [_, variables] = get_lagrangian(settings, names);
  auto rhs = get_shorthand_rhs(variables);
  EXPECT_EQ(variables.size(), rhs.size());
  for (size_t i = 0; i < variables.size(); ++i) {
    EXPECT_EQ("-r_{" + variables[i]->to_string() + "}", rhs[i]->to_string());
  }
}

TEST_F(OptimizationTest, GaussianElimination) {
  auto settings = Settings();
  auto names = VariableNames();
  auto [lhs, rhs, variables, _] = get_newton_system(settings, names);

  // Save original size
  size_t original_size = lhs.size();

  // Perform one elimination step
  gaussian_elimination(lhs, rhs, lhs.size() - 1);

  // Check sizes after elimination
  EXPECT_EQ(lhs.size(), original_size - 1);
  EXPECT_EQ(rhs.size(), original_size - 1);
  for (const auto& row : lhs) {
    EXPECT_EQ(row.size(), original_size - 1);
  }
}

TEST_F(OptimizationTest, GetAugmentedSystem) {
  auto settings = Settings();
  settings.inequalities = Bounds::Lower;
  auto names = VariableNames();
  auto newton_system = get_newton_system(settings, names);

  // Save original size
  size_t original_size = newton_system.lhs.size();
  EXPECT_GT(original_size, 2);  // Make sure we have something to eliminate

  // Get augmented system
  auto augmented_system = get_augmented_system(std::move(newton_system));

  // Check that variables were eliminated and delta definitions were created
  EXPECT_LT(augmented_system.lhs.size(), original_size);
  EXPECT_GT(augmented_system.delta_definitions.size(), 0);

  // Check that first delta definition references the appropriate variable
  const auto& [delta_var, _] = augmented_system.delta_definitions.front();
  EXPECT_TRUE(delta_var->to_string().find("\\Delta") != std::string::npos);
}

TEST_F(OptimizationTest, GetNormalEquations) {
  auto settings = Settings();
  settings.inequalities = Bounds::Lower;
  auto names = VariableNames();
  auto newton_system = get_newton_system(settings, names);

  // Save original size
  size_t original_size = newton_system.lhs.size();
  EXPECT_GT(original_size, 2);  // Make sure we have something to eliminate

  // Get normal equations
  auto normal_equations = get_normal_equations(std::move(newton_system));

  // Check that we have only one equation left
  EXPECT_EQ(normal_equations.lhs.size(), 1);
  EXPECT_EQ(normal_equations.rhs.size(), 1);
  EXPECT_EQ(normal_equations.variables.size(), 1);

  // Check that we have delta definitions
  EXPECT_EQ(normal_equations.delta_definitions.size(), original_size - 1);
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}