#pragma once
#include <stdexcept>
#include <string>

namespace Utils {

class AssertionError : public std::logic_error {
 public:
  template <typename... Args>
  explicit AssertionError(Args&&... args)
      : std::logic_error(std::forward<Args>(args)...) {}
};

// AssertionHandler handles the assertion failure and prints stack trace
void assertion_handler(bool condition, const char* condition_str,
                       const char* file, int line, const char* func,
                       const std::string& message = "");

}  // namespace Utils

// Helper macro to check if there are one or two arguments
#define ASSERT_GET_MACRO(_1, _2, NAME, ...) NAME

#define ASSERT_MESSAGE(condition)                          \
  (std::string("Assertion failed at ") + __FILE__ + ": " + \
   std::to_string(__LINE__) + " in " + __func__ + ": " + #condition)

// ASSERT macro overloads for one or two arguments
#define ASSERT1(condition)                                              \
  Utils::assertion_handler((condition), #condition, __FILE__, __LINE__, \
                           __func__)

#define ASSERT2(condition, message)                                     \
  Utils::assertion_handler((condition), #condition, __FILE__, __LINE__, \
                           __func__, message)

// Main ASSERT macro that selects the appropriate overload
#define ASSERT(...) ASSERT_GET_MACRO(__VA_ARGS__, ASSERT2, ASSERT1)(__VA_ARGS__)