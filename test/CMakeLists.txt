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

# Add test executable
add_executable(Expression_test
    Expression_test.cpp
)

# Link test executable with Google Test and Expression library
target_link_libraries(Expression_test PRIVATE 
    Expression
    GTest::gtest_main
)

# Register tests
include(GoogleTest)
gtest_discover_tests(Expression_test)
