#include "SymbolicOptimization.h"

#include <gtest/gtest.h>

#include "ExprFactory.h"

using namespace SymbolicOptimization;

class SymbolicOptimizationTest : public ::testing::Test {
 protected:
  void SetUp() override {}
};

TEST_F(SymbolicOptimizationTest, GetLagrangian) {
  auto settings = Settings();
  auto names = VariableNames();
  names.s_A_ineq = "s123";
  names.s_A_ineq_u = "t123";
  names.s_A_ineq_l = "g123";
  names.s_x_l = "sxl123";
  names.s_x_u = "sxu123";
  auto problem = get_optimization_problem(
      settings, names, OptimizationProblemType::SlackedWithBarriers);
  auto lagrangian = get_lagrangian(problem);
  const auto str = lagrangian->to_string();
  EXPECT_NE(str.find(names.s_A_ineq), std::string::npos);
  EXPECT_NE(str.find(names.s_A_ineq_l), std::string::npos);
  EXPECT_NE(str.find(names.s_A_ineq_u), std::string::npos);
  EXPECT_NE(str.find(names.s_x_l), std::string::npos);
  EXPECT_NE(str.find(names.s_x_u), std::string::npos);
}

TEST_F(SymbolicOptimizationTest, GetLagrangianLower) {
  auto settings = Settings();
  settings.inequalities = Bounds::Lower;
  settings.variable_bounds = Bounds::Lower;
  auto names = VariableNames();
  names.s_A_ineq_u = "t123";
  names.s_A_ineq_l = "g123";
  auto problem = get_optimization_problem(
      settings, names, OptimizationProblemType::SlackedWithBarriers);
  auto lagrangian = get_lagrangian(problem);
  const auto str = lagrangian->to_string();
  EXPECT_NE(str.find(names.s_A_ineq), std::string::npos);
  EXPECT_NE(str.find(names.s_A_ineq_l), std::string::npos);
  EXPECT_EQ(str.find(names.s_A_ineq_u), std::string::npos);
  EXPECT_NE(str.find(names.s_x_l), std::string::npos);
  EXPECT_EQ(str.find(names.s_x_u), std::string::npos);
}

TEST_F(SymbolicOptimizationTest, GetLagrangianUpper) {
  auto settings = Settings();
  settings.inequalities = Bounds::Upper;
  settings.variable_bounds = Bounds::Upper;
  auto names = VariableNames();
  names.s_A_ineq_u = "t123";
  names.s_A_ineq_l = "g123";
  auto problem = get_optimization_problem(
      settings, names, OptimizationProblemType::SlackedWithBarriers);
  auto lagrangian = get_lagrangian(problem);
  const auto str = lagrangian->to_string();
  EXPECT_NE(str.find(names.s_A_ineq), std::string::npos);
  EXPECT_EQ(str.find(names.s_A_ineq_l), std::string::npos);
  EXPECT_NE(str.find(names.s_A_ineq_u), std::string::npos);
  EXPECT_EQ(str.find(names.s_x_l), std::string::npos);
  EXPECT_NE(str.find(names.s_x_u), std::string::npos);
}

TEST_F(SymbolicOptimizationTest, GetLagrangianNoBounds) {
  auto settings = Settings();
  settings.inequalities = Bounds::None;
  settings.variable_bounds = Bounds::None;
  auto names = VariableNames();
  names.s_A_ineq_l = "g123";
  names.s_A_ineq_u = "t123";
  names.s_x_l = "y123";
  names.s_x_u = "z123";
  auto problem = get_optimization_problem(
      settings, names, OptimizationProblemType::SlackedWithBarriers);
  auto lagrangian = get_lagrangian(problem);
  const auto str = lagrangian->to_string();
  EXPECT_EQ(str.find(names.s_A_ineq), std::string::npos);
  EXPECT_EQ(str.find(names.s_A_ineq_l), std::string::npos);
  EXPECT_EQ(str.find(names.s_A_ineq_u), std::string::npos);
  EXPECT_EQ(str.find(names.s_x_l), std::string::npos);
  EXPECT_EQ(str.find(names.s_x_u), std::string::npos);
}

TEST_F(SymbolicOptimizationTest, GetFirstOrderOptimalityConditions) {
  auto settings = Settings();
  auto names = VariableNames();
  auto [first_order, variables] =
      get_first_order_optimality_conditions(settings, names);
  EXPECT_EQ(first_order.size(), variables.size());
}

TEST_F(SymbolicOptimizationTest, GetNewtonSystem) {
  auto settings = Settings();
  auto names = VariableNames();
  auto [lhs, rhs, variables, _] = get_newton_system(settings, names);
  EXPECT_EQ(lhs.size(), rhs.size());
  EXPECT_EQ(lhs.size(), variables.size());
  for (const auto& row : lhs) {
    EXPECT_EQ(row.size(), variables.size());
  }
}

TEST_F(SymbolicOptimizationTest, GetShorthandRhs) {
  const auto settings = Settings();
  const auto names = VariableNames();
  const auto newton_system = get_newton_system(settings, names);
  const auto rhs = get_shorthand_rhs(newton_system).shorthand_rhs;
  const auto& variables = newton_system.variables;
  EXPECT_EQ(variables.size(), rhs.size());
  for (size_t i = 0; i < variables.size(); ++i) {
    EXPECT_EQ("-r_{" + variables[i]->to_string() + "}", rhs[i]->to_string());
  }
}

TEST_F(SymbolicOptimizationTest, GaussianElimination) {
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

TEST_F(SymbolicOptimizationTest, GetAugmentedSystem) {
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

TEST_F(SymbolicOptimizationTest, GetNormalEquations) {
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