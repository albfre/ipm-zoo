cmake_minimum_required(VERSION 3.14)
project(SymbolicOptimizationProject VERSION 1.0)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Check if we're building for WebAssembly
if(EMSCRIPTEN)
  message(STATUS "Building for WebAssembly with Emscripten")
  
  # Set Emscripten-specific compiler flags
  set(CMAKE_EXECUTABLE_SUFFIX ".js")
  
  # Common Emscripten flags
  set(EMSCRIPTEN_COMMON_FLAGS
    "-s WASM=1"
    "-s ALLOW_MEMORY_GROWTH=1"
    "-s MODULARIZE=1"
    "-s EXPORT_ES6=0"
    "-s EXPORTED_RUNTIME_METHODS=['ccall','cwrap']"
    "-s EXPORT_NAME=SymbolicOptimizationModule"
  )
  
  # Convert the list to a string
  string(REPLACE ";" " " EMSCRIPTEN_LINK_FLAGS "${EMSCRIPTEN_COMMON_FLAGS}")
  
  # Set these flags for linking
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${EMSCRIPTEN_LINK_FLAGS}")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${EMSCRIPTEN_LINK_FLAGS}")
  
  # Disable testing for WebAssembly build
  set(BUILD_TESTING OFF)
else()
  message(STATUS "Building for native platform")
  # Enable testing for native build
  enable_testing()
endif()

# Add source directory
add_subdirectory(src)

# Only add test directory for native builds
if(NOT EMSCRIPTEN)
  add_subdirectory(test)
endif()