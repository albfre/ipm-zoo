# Create Expression library
add_library(Expression
    Expression.cpp
    Expression.h
)

# Make include directories visible to other targets
target_include_directories(Expression PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Create main executable
add_executable(ipmZoo ipmZoo.cpp)
target_link_libraries(ipmZoo PRIVATE Expression)
