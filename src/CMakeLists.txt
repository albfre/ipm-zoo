# Create SymbolicOptimization library
add_library(SymbolicOptimization
    Expression.cpp
    Expression.h
    Optimization.cpp
    Optimization.h
)

# Make include directories visible to other targets
target_include_directories(SymbolicOptimization PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Create main executable for native builds
if(NOT EMSCRIPTEN)
    add_executable(ipmZoo ipmZoo.cpp)
    target_link_libraries(ipmZoo PRIVATE SymbolicOptimization)
else()
    # Create the WebAssembly bindings
    add_executable(symbolicOptimizationWasm wasm_bindings.cpp)
    target_link_libraries(symbolicOptimizationWasm PRIVATE SymbolicOptimization)
    
    # Add specific flags for this target
    target_link_options(symbolicOptimizationWasm PRIVATE
        "-s EXPORT_NAME=SymbolicOptimizationModule"
        "-s EXPORTED_FUNCTIONS=['_getLagrangian','_malloc','_free']"
        "-s EXPORTED_RUNTIME_METHODS=['ccall','cwrap'']"
        "--bind"
    )
    
    # Set output file names
    set_target_properties(symbolicOptimizationWasm PROPERTIES
        OUTPUT_NAME "symbolic-optimization-wasm"
    )
endif()