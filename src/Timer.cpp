#include "Timer.h"

#include <iomanip>
#include <iostream>

namespace Util {

void Timer::start(const std::string& operation) {
  startTimes[operation] = std::chrono::high_resolution_clock::now();
}

void Timer::stop(const std::string& operation) {
  auto endTime = std::chrono::high_resolution_clock::now();
  auto startTimeIter = startTimes.find(operation);

  if (startTimeIter == startTimes.end()) {
    std::cerr << "Warning: Attempting to stop timer for operation '"
              << operation << "' that was not started" << std::endl;
    return;
  }

  auto startTime = startTimeIter->second;
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
      endTime - startTime);

  if (totalTimes.find(operation) == totalTimes.end()) {
    operationOrder.push_back(operation);
    totalTimes[operation] = duration;
  } else {
    totalTimes[operation] += duration;
  }
}

void Timer::report() const {
  for (const auto& operation : operationOrder) {
    const auto& duration = totalTimes.at(operation);

    std::cout << operation << ": ";

    // Format time appropriately
    if (duration.count() < 1) {
      double microseconds = duration.count() * 1000.0;
      std::cout << std::fixed << std::setprecision(3) << microseconds << " μs";
    } else if (duration.count() < 1000) {
      std::cout << std::fixed << std::setprecision(3) << duration.count()
                << " ms";
    } else {
      double seconds = duration.count() / 1000.0;
      std::cout << std::fixed << std::setprecision(3) << seconds << " s";
    }

    std::cout << std::endl;
  }
}
}  // namespace Util
