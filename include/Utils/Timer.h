#pragma once
#include <chrono>
#include <string>
#include <unordered_map>
#include <vector>

namespace Utils {
class Timer {
 public:
  Timer() = default;

  void start(const std::string& operation);
  void stop(const std::string& operation);
  void report() const;

 private:
  std::unordered_map<std::string,
                     std::chrono::high_resolution_clock::time_point>
      start_times;
  std::unordered_map<std::string, std::chrono::milliseconds> total_times;
  std::vector<std::string> operation_order;
};
}  // namespace Util