#include "ExprFactory.h"

#include <gtest/gtest.h>

#include "Utils/Helpers.h"

using namespace Expression;

// Test class with friend access to ExprFactory
namespace Expression {
class ExprFactoryTest : public ::testing::Test {
 protected:
  void SetUp() override {
    // Reset the ExprFactory cache for each test by recreating the instance
    resetFactory();
  }

  ExprFactory& factory() { return ExprFactory::instance_(); }

  void resetFactory() {
    auto& f = factory();
    f.cache_.clear();
    f.cleanupCounter_ = 0;
  }

  void setCleanupCounter(size_t count) {
    auto& f = factory();
    f.cleanupCounter_ = count;
  }

  bool isCached(const ExprPtr& expr) {
    const auto& f = factory();
    const auto key = expr->toExpressionString();
    auto it = f.cache_.find(key);
    if (it != f.cache_.end()) {
      return !it->second.expired();
    }
    return false;
  }

  size_t cacheSize() {
    const auto& f = factory();
    return f.cache_.size();
  }

  // Get number of valid (non-expired) entries in cache
  size_t validCacheEntries() {
    const auto& f = factory();
    size_t count = 0;
    for (const auto& [key, ptr] : f.cache_) {
      if (!ptr.expired()) {
        count++;
      }
    }
    return count;
  }
};
}  // namespace Expression

// Test basic expression creation
TEST_F(ExprFactoryTest, CreateBasicExpressions) {
  // Test creation of basic expressions
  auto num = ExprFactory::number(42.0);
  EXPECT_TRUE(is<Number>(num));
  EXPECT_EQ(std::get<Number>(num->getImpl()).value, 42.0);

  auto var = ExprFactory::variable("x");
  EXPECT_TRUE(is<Variable>(var));
  EXPECT_EQ(std::get<Variable>(var->getImpl()).name, "x");

  auto scalar = ExprFactory::namedScalar("a");
  EXPECT_TRUE(is<NamedScalar>(scalar));
  EXPECT_EQ(std::get<NamedScalar>(scalar->getImpl()).name, "a");

  auto vector = ExprFactory::namedVector("v");
  EXPECT_TRUE(is<NamedVector>(vector));
  EXPECT_EQ(std::get<NamedVector>(vector->getImpl()).name, "v");

  auto matrix = ExprFactory::matrix("M");
  EXPECT_TRUE(is<Matrix>(matrix));
  EXPECT_EQ(std::get<Matrix>(matrix->getImpl()).name, "M");

  auto symMatrix = ExprFactory::symmetricMatrix("S");
  EXPECT_TRUE(is<SymmetricMatrix>(symMatrix));
  EXPECT_EQ(std::get<SymmetricMatrix>(symMatrix->getImpl()).name, "S");
}

// Test unary operations
TEST_F(ExprFactoryTest, UnaryOperations) {
  auto x = ExprFactory::variable("x");

  auto diagX = ExprFactory::diagonalMatrix(x);
  EXPECT_TRUE(is<DiagonalMatrix>(diagX));
  EXPECT_TRUE(std::get<DiagonalMatrix>(diagX->getImpl()).child == x);

  auto transposeX = ExprFactory::transpose(x);
  EXPECT_TRUE(is<Transpose>(transposeX));
  EXPECT_TRUE(std::get<Transpose>(transposeX->getImpl()).child == x);

  auto negateX = ExprFactory::negate(x);
  EXPECT_TRUE(is<Negate>(negateX));
  EXPECT_TRUE(std::get<Negate>(negateX->getImpl()).child == x);

  auto invertX = ExprFactory::invert(x);
  EXPECT_TRUE(is<Invert>(invertX));
  EXPECT_TRUE(std::get<Invert>(invertX->getImpl()).child == x);

  auto logX = ExprFactory::log(x);
  EXPECT_TRUE(is<Log>(logX));
  EXPECT_TRUE(std::get<Log>(logX->getImpl()).child == x);
}

// Test n-ary operations
TEST_F(ExprFactoryTest, NaryOperations) {
  auto x = ExprFactory::variable("x");
  auto y = ExprFactory::variable("y");
  auto z = ExprFactory::variable("z");

  // Test sum
  auto sum = ExprFactory::sum({x, y, z});
  EXPECT_TRUE(is<Sum>(sum));
  const auto& sumTerms = std::get<Sum>(sum->getImpl()).terms;
  EXPECT_EQ(sumTerms.size(), 3);
  EXPECT_EQ(sumTerms[0], x);
  EXPECT_EQ(sumTerms[1], y);
  EXPECT_EQ(sumTerms[2], z);

  // Test product
  auto prod = ExprFactory::product({x, y, z});
  EXPECT_TRUE(is<Product>(prod));
  const auto& prodTerms = std::get<Product>(prod->getImpl()).terms;
  EXPECT_EQ(prodTerms.size(), 3);
  EXPECT_EQ(prodTerms[0], x);
  EXPECT_EQ(prodTerms[1], y);
  EXPECT_EQ(prodTerms[2], z);
}

// Test empty/single term operations
TEST_F(ExprFactoryTest, EmptyAndSingleTermOperations) {
  auto x = ExprFactory::variable("x");

  // Empty sum should return 0
  auto emptySum = ExprFactory::sum({});
  EXPECT_TRUE(is<Number>(emptySum));
  EXPECT_EQ(std::get<Number>(emptySum->getImpl()).value, 0.0);

  // Single term sum should return the term
  auto singleSum = ExprFactory::sum({x});
  EXPECT_EQ(singleSum, x);

  // Empty product should return 1
  auto emptyProd = ExprFactory::product({});
  EXPECT_TRUE(is<Number>(emptyProd));
  EXPECT_EQ(std::get<Number>(emptyProd->getImpl()).value, 1.0);

  // Single term product should return the term
  auto singleProd = ExprFactory::product({x});
  EXPECT_EQ(singleProd, x);
}

// Test cache functionality
TEST_F(ExprFactoryTest, CacheFunctionality) {
  // Initially the cache should be empty
  EXPECT_EQ(cacheSize(), 0);

  // Create some expressions
  auto x = ExprFactory::variable("x");
  auto y = ExprFactory::variable("y");

  // Check that they're cached
  EXPECT_TRUE(isCached(x));
  EXPECT_TRUE(isCached(y));
  EXPECT_EQ(validCacheEntries(), 2);

  // Create same expressions again and verify they come from cache
  auto x2 = ExprFactory::variable("x");
  EXPECT_EQ(x, x2);                   // Should be the same object
  EXPECT_EQ(x.get(), x2.get());       // Should be the same pointer
  EXPECT_EQ(validCacheEntries(), 2);  // Cache should still have 2 entries

  // Create a more complex expression
  auto sum1 = ExprFactory::sum({x, y});
  EXPECT_TRUE(isCached(sum1));
  EXPECT_EQ(validCacheEntries(), 3);

  // Create the same complex expression and verify it comes from cache
  auto sum2 = ExprFactory::sum({x, y});
  EXPECT_EQ(sum1, sum2);
  EXPECT_EQ(sum1.get(), sum2.get());
  EXPECT_EQ(validCacheEntries(), 3);

  // Create a different complex expression
  auto prod = ExprFactory::product({x, y});
  EXPECT_TRUE(isCached(prod));
  EXPECT_EQ(validCacheEntries(), 4);

  // Test cache expiration
  {
    // Create a temporary expression that will go out of scope
    auto tempVar = ExprFactory::variable("temp");
    EXPECT_TRUE(isCached(tempVar));
    EXPECT_EQ(validCacheEntries(), 5);
  }

  // Now tempVar is out of scope, but the cache entry is still there
  EXPECT_EQ(cacheSize(), 5);
  // But it should be marked as expired
  EXPECT_EQ(validCacheEntries(), 4);

  // Force cleanup by creating many expressions to trigger the cleanup counter
  setCleanupCounter(999);  // Set counter near threshold
  auto triggerCleanup = ExprFactory::variable("cleanup");

  // Cache should have been cleaned up
  EXPECT_EQ(validCacheEntries(), 5);  // 4 previous valid + 1 new
  EXPECT_EQ(cacheSize(), 5);  // At least one expired entry should be removed
}

// Test asPtr method
TEST_F(ExprFactoryTest, AsPtr) {
  auto x = ExprFactory::variable("x");
  auto y = ExprFactory::variable("y");

  // Dereference the pointer
  auto& expr = *x;

  // Convert it back to ExprPtr using asPtr
  auto ptr = ExprFactory::asPtr(expr);

  // Should be the same as the original
  EXPECT_EQ(ptr, x);
  EXPECT_EQ(ptr.get(), x.get());
}

// Test getExpr method
TEST_F(ExprFactoryTest, GetExpr) {
  // Create a variant directly
  Expr::ExprVariant variant = Variable("direct");

  // Use getExpr to create an ExprPtr from it
  auto expr = ExprFactory::getExpr(std::move(variant));

  // Verify it works and is cached
  EXPECT_TRUE(is<Variable>(expr));
  EXPECT_EQ(std::get<Variable>(expr->getImpl()).name, "direct");
  EXPECT_TRUE(isCached(expr));
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
