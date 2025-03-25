# Find GoogleTest package
find_package(GTest QUIET)

# If GTest not found, use FetchContent to download it
if(NOT GTEST_FOUND)
    include(FetchContent)
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG release-1.12.1
    )
    FetchContent_MakeAvailable(googletest)
endif()

# Define test sources
set(TEST_SOURCES 
    Expr_test.cpp 
    ExprFactory_test.cpp
    SymbolicOptimization_test.cpp
    Helpers_test.cpp
    Evaluation_test.cpp
)

# Create test executables in a loop
foreach(TEST_SOURCE ${TEST_SOURCES})
    get_filename_component(TEST_NAME ${TEST_SOURCE} NAME_WE)
    add_executable(${TEST_NAME} ${TEST_SOURCE})

    target_link_libraries(${TEST_NAME} PRIVATE 
        SymbolicOptimization 
        GTest::gtest_main
    )

    gtest_discover_tests(${TEST_NAME})
endforeach()
