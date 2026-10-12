cmake_minimum_required(VERSION 3.14)
project(NonOverlapDecomposeGrid LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)

add_library(
    grid
    STATIC
    src/vertex.cpp
    src/dim_iterator.cpp
    src/grid.cpp
)
target_link_libraries(grid PUBLIC utils)
target_include_directories(
    grid
    PUBLIC
    include
)



add_executable(
    grid_test
    tests/grid.cpp
)

target_link_libraries(grid_test PUBLIC grid)
target_link_libraries(grid_test PRIVATE Catch2::Catch2WithMain)