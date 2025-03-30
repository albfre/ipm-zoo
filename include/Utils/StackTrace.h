#pragma once
#include <iostream>
#include <string>
#include <vector>

namespace Utils {

struct StackFrame {
  void* address;         // Memory address
  std::string file;      // Source file
  std::string function;  // Function name
  int line;              // Line number

  StackFrame()
      : address(nullptr), file("unknown"), function("unknown"), line(0) {}
};

// Stack trace functions
void print_stack_trace(std::ostream& os = std::cerr);
std::vector<StackFrame> get_stack_trace(int skip_frames = 1);

// Sets up the global termination handler for uncaught exceptions
void initialize_termination_handler();

}  // namespace Utils
