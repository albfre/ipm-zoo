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
}  // namespace Utils

// Helper macro to check if there are one or two arguments
#define ASSERT_GET_MACRO(_1, _2, NAME, ...) NAME

#define ASSERT_MESSAGE(condition)                         \
  (std::string("Assertion faied at ") + __FILE__ + ": " + \
   std::to_string(__LINE__) + " in " + __func__ + ": " + #condition)

// ASSERT macro overloads for one or two arguments
#define ASSERT1(condition)                                  \
  if (!(condition)) {                                       \
    throw Utils::AssertionError(ASSERT_MESSAGE(condition)); \
  }

#define ASSERT2(condition, message)                                          \
  if (!(condition)) {                                                        \
    throw Utils::AssertionError(ASSERT_MESSAGE(condition) + "; " + message); \
  }

// Main ASSERT macro that selects the appropriate overload
#define ASSERT(...) ASSERT_GET_MACRO(__VA_ARGS__, ASSERT2, ASSERT1)(__VA_ARGS__)