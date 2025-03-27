#include "Utils/Timer.h"

#include <iomanip>
#include <iostream>

namespace Utils {

void Timer::start(const std::string& operation) {
  start_times[operation] = std::chrono::high_resolution_clock::now();
}

void Timer::stop(const std::string& operation) {
  auto end_time = std::chrono::high_resolution_clock::now();
  auto start_time_iter = start_times.find(operation);

  if (start_time_iter == start_times.end()) {
    std::cerr << "Warning: Attempting to stop timer for operation '"
              << operation << "' that was not started" << std::endl;
    return;
  }

  auto start_time = start_time_iter->second;
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
      end_time - start_time);

  if (total_times.find(operation) == total_times.end()) {
    operation_order.push_back(operation);
    total_times[operation] = duration;
  } else {
    total_times[operation] += duration;
  }
}

void Timer::report() const {
  for (const auto& operation : operation_order) {
    const auto& duration = total_times.at(operation);

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
}  // namespace Utils
