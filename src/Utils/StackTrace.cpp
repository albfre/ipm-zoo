#include "Utils/StackTrace.h"

#include <cstdlib>
#include <iomanip>
#include <memory>
#include <sstream>

#include "Utils/Assert.h"

#ifndef __EMSCRIPTEN__
#include <cxxabi.h>    // For demangling C++ symbols
#include <dlfcn.h>     // For dladdr
#include <execinfo.h>  // For backtrace functions
#endif

namespace Utils {

std::vector<StackFrame> get_stack_trace(int skip_frames) {
#ifdef __EMSCRIPTEN__
  // WASM implementation
  std::vector<StackFrame> stack_frames;
  StackFrame frame;
  frame.address = nullptr;
  frame.function = "Stack traces not available in WebAssembly";
  frame.file = "unknown";
  frame.line = 0;
  stack_frames.push_back(frame);
  return stack_frames;
#else
  // Native implementation
  constexpr int kMaxFrames = 128;
  void* callstack[kMaxFrames];

  // Get the addresses in the stack
  const auto num_frames = backtrace(callstack, kMaxFrames);

  // Convert addresses to strings
  char** symbols = backtrace_symbols(callstack, num_frames);
  if (symbols == nullptr) {
    return {};  // Failed to get symbols
  }

  std::vector<StackFrame> stack_frames;

  // Skip the first few frames (usually the error handling functions)
  for (int i = skip_frames; i < num_frames; i++) {
    StackFrame frame;
    frame.address = callstack[i];

    std::string symbol = symbols[i];

    // Parse the raw symbol string which is typically:
    // "binary(function+offset) [address]" or "binary [address]"
    size_t pos_open = symbol.find('(');
    size_t pos_plus = symbol.find('+', pos_open);
    size_t pos_close = symbol.find(')', pos_plus);
    size_t pos_bracket = symbol.find('[');

    // Extract binary path
    auto binary_path = symbol.substr(
        0, pos_open != std::string::npos ? pos_open : pos_bracket);
    binary_path =
        binary_path.substr(0, binary_path.find_last_not_of(" \t") + 1);

    // Try to get function name
    if (pos_open != std::string::npos && pos_plus != std::string::npos &&
        pos_close != std::string::npos) {
      std::string mangled_name =
          symbol.substr(pos_open + 1, pos_plus - pos_open - 1);

      // Demangle C++ symbol
      int status;
      char* demangled =
          abi::__cxa_demangle(mangled_name.c_str(), nullptr, nullptr, &status);

      if (status == 0 && demangled != nullptr) {
        frame.function = demangled;
        free(demangled);
      } else {
        frame.function = mangled_name;
      }
    } else {
      frame.function = "??";
    }

    // Try to get file and line info using addr2line (if available)
    Dl_info info;
    if (dladdr(frame.address, &info)) {
      frame.file = info.dli_fname ? info.dli_fname : "??";
    } else {
      frame.file = binary_path;
    }

    // Save the original symbol string as well
    frame.function += " [" + symbol + "]";

    stack_frames.push_back(frame);
  }

  free(symbols);
  return stack_frames;
#endif
}

void print_stack_trace(std::ostream& os) {
  auto stack_frames =
      get_stack_trace(2);  // Skip this method and get_stack_trace

  os << "\nStack trace:" << std::endl;
  for (size_t i = 0; i < stack_frames.size(); ++i) {
    const auto& frame = stack_frames[i];
    os << "[" << i << "] ";

    // Print function name if available
    os << frame.function;

    // Print file location if available
    if (frame.file != "unknown" && frame.file != "??") {
      os << "\n    at " << frame.file;
      if (frame.line > 0) {
        os << ":" << frame.line;
      }
    }

    os << std::endl;
  }
}

void initialize_termination_handler() {
  std::set_terminate([]() {
    std::stringstream message;
    message << "\n========== UNHANDLED EXCEPTION ==========\n";
    try {
      std::rethrow_exception(std::current_exception());
    } catch (const Utils::AssertionError&) {
      // Do nothing, stack trace already printed
      std::abort();
    } catch (const std::exception& e) {
      message << "Exception type: " << typeid(e).name() << "\n";
      message << "What: " << e.what() << "\n";
    } catch (...) {
      message << "Unknown exception type\n";
    }

    // Print stack trace
    std::cerr << message.str();
    print_stack_trace();
    std::cerr << "==========================================\n";
    std::abort();
  });
}

// Automatically initialize the termination handler when this object is
// constructed
namespace {
struct AutoInitializeTerminationHandler {
  AutoInitializeTerminationHandler() { initialize_termination_handler(); }
};
static AutoInitializeTerminationHandler g_auto_init;
}  // namespace

}  // namespace Utils
