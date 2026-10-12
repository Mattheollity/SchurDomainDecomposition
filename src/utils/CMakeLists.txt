cmake_minimum_required(VERSION 3.14)
project(NonOverlapDecomposeUtils LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)

add_library(
    utils
    STATIC
    src/file_io.cpp
    src/problem_conditions.cpp
)

target_include_directories(
    utils
    PUBLIC
    include
    ${EIGEN3_INCLUDE_DIR}
    ${CURSES_INCLUDE_DIR}
)

target_link_libraries(utils PUBLIC OpenMP::OpenMP_CXX)
target_link_libraries(utils PUBLIC Eigen3::Eigen)