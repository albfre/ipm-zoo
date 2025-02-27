cmake_minimum_required(VERSION 3.14)
project(ExpressionProject VERSION 1.0)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Enable testing at the top level
enable_testing()

# Add source and test directories
add_subdirectory(src)
add_subdirectory(test)
