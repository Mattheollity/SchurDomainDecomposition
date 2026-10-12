cmake_minimum_required(VERSION 3.14)
project(PoissonApprox LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)

add_library(
    poisson
    STATIC
    src/poisson.cpp
    src/approx.cpp
)
target_link_libraries(poisson grid)

target_include_directories(
    poisson
    PUBLIC
    include
)



add_executable(
    poisson_test
    tests/poisson.cpp
)
target_link_libraries(poisson_test PUBLIC poisson)
target_link_libraries(poisson_test PRIVATE Catch2::Catch2WithMain)