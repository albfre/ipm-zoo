# Create SymbolicOptimization library
set(LIBRARY_SOURCES
    Expr.cpp
    ExprFactory.cpp
    SymbolicOptimization.cpp
    Visitors/DifferentiationVisitor.cpp
    Visitors/SimplificationVisitor.cpp
    Visitors/ToStringVisitor.cpp
    NumericOptimization/Evaluation.cpp
    Utils/Timer.cpp
)

# Create SymbolicOptimization library
add_library(SymbolicOptimization ${LIBRARY_SOURCES})

# Make include directories visible
target_include_directories(SymbolicOptimization 
    PUBLIC 
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
)

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

# Add install targets for better distribution support
install(TARGETS SymbolicOptimization
    EXPORT SymbolicOptimizationTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include
)

if(NOT EMSCRIPTEN)
    install(TARGETS IpmZoo
        RUNTIME DESTINATION bin
    )
endif()