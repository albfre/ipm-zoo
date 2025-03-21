#pragma once
#include <chrono>
#include <string>
#include <unordered_map>
#include <vector>

namespace Util {

class Timer {
 public:
  Timer() = default;

  void start(const std::string& operation);
  void stop(const std::string& operation);
  void report() const;

 private:
  std::unordered_map<std::string,
                     std::chrono::high_resolution_clock::time_point>
      startTimes;
  std::unordered_map<std::string, std::chrono::milliseconds> totalTimes;
  std::vector<std::string> operationOrder;
};
}  // namespace Util