# Create SymbolicOptimization library
add_library(SymbolicOptimization
    Expression.cpp
    Expression.h
    SymbolicOptimization.cpp
    SymbolicOptimization.h
    Helpers.h
    ToStringVisitor.h
    DifferentiationVisitor.h
    DifferentiationVisitor.cpp
    SimplificationVisitor.h
    SimplificationVisitor.cpp
    Evaluation.h
    Timer.h
    Timer.cpp
)

# Make include directories visible to other targets
target_include_directories(SymbolicOptimization PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Create main executable for native builds
if(NOT EMSCRIPTEN)
    add_executable(IpmZoo IpmZoo.cpp)
    target_link_libraries(IpmZoo PRIVATE SymbolicOptimization)
else()
    # Create the WebAssembly bindings
    add_executable(symbolicOptimizationWasm wasm_bindings.cpp)
    target_link_libraries(symbolicOptimizationWasm PRIVATE SymbolicOptimization)
    
    # Add specific flags for this target
    target_link_options(symbolicOptimizationWasm PRIVATE
        "-s EXPORT_NAME=SymbolicOptimizationModule"
        "-s EXPORTED_FUNCTIONS=['_malloc','_free']"
        "-s EXPORTED_RUNTIME_METHODS=['ccall','cwrap']"
        "--bind"
    )
    
    # Set output file names
    set_target_properties(symbolicOptimizationWasm PROPERTIES
        OUTPUT_NAME "symbolic-optimization-wasm"
    )
endif()