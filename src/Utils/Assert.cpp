#include "Utils/Assert.h"

#include <iostream>
#include <sstream>

#include "Utils/StackTrace.h"

namespace Utils {

void assertion_handler(bool condition, const char* condition_str,
                       const char* file, int line, const char* func,
                       const std::string& message) {
  if (!condition) {
    // Create detailed error message
    std::stringstream ss;
    ss << "Assertion failed at " << file << ":" << line << " in " << func
       << ": " << condition_str;

    if (!message.empty()) {
      ss << "; " << message;
    }

    // Print the error message and stack trace
    std::cerr << "\n========== ASSERTION FAILED ==========\n";
    std::cerr << ss.str() << std::endl;
    print_stack_trace();
    std::cerr << "==========================================\n";

    // Throw an exception with the detailed message
    throw AssertionError(ss.str());
  }
}
}  // namespace Utils
