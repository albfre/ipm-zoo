#include <chrono>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

class Timer {
 private:
  std::unordered_map<std::string,
                     std::chrono::high_resolution_clock::time_point>
      startTimes;
  std::unordered_map<std::string, std::chrono::milliseconds> totalTimes;
  std::vector<std::string> operationOrder;

 public:
  void start(const std::string& operation) {
    startTimes[operation] = std::chrono::high_resolution_clock::now();
  }

  void stop(const std::string& operation) {
    auto endTime = std::chrono::high_resolution_clock::now();
    auto startTime = startTimes[operation];
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
        endTime - startTime);
    if (!totalTimes.contains(operation)) {
      operationOrder.push_back(operation);
    }
    totalTimes[operation] += duration;
  }

  void report() {
    for (const auto& operation : operationOrder) {
      std::cout << operation << ": " << totalTimes.at(operation).count()
                << " ms" << std::endl;
    }
  }
};