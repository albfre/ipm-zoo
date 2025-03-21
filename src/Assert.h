#pragma once

#include <stdexcept>
#include <string>

// Helper macro to check if there are one or two arguments
#define ASSERT_GET_MACRO(_1, _2, NAME, ...) NAME

// ASSERT macro overloads for one or two arguments
#define ASSERT1(condition)                                                  \
  if (!(condition)) {                                                       \
    throw std::logic_error(std::string("Assertion failed at ") + __FILE__ + \
                           ":" + std::to_string(__LINE__) + " in " +        \
                           __func__ + ": " + #condition);                   \
  }

#define ASSERT2(condition, message)                                      \
  if (!(condition)) {                                                    \
    throw std::logic_error(std::string("Assertion failed: ") + message); \
  }

// Main ASSERT macro that selects the appropriate overload
#define ASSERT(...) ASSERT_GET_MACRO(__VA_ARGS__, ASSERT2, ASSERT1)(__VA_ARGS__)