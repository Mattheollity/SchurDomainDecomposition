cmake_minimum_required(VERSION 3.14)
project(TimeScheme LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)

find_package(eigen-expokit 0.1 REQUIRED)

add_library(
    time_scheme
    STATIC
    src/time_scheme.cpp
)
target_link_libraries(time_scheme PUBLIC poisson eigen-expokit::eexpokit)

target_include_directories(
    time_scheme
    PUBLIC
    include
)

add_executable(
    time_scheme_test
    tests/time_scheme.cpp
)
target_link_libraries(time_scheme_test PUBLIC time_scheme)
target_link_libraries(time_scheme_test PRIVATE Catch2::Catch2WithMain)

add_executable(
    schemeapprox
    src/scheme_approx.cpp
)
target_link_libraries(schemeapprox PUBLIC time_scheme)

add_executable(
    schemeconverg
    src/scheme_converg.cpp
)
target_link_libraries(schemeconverg PUBLIC time_scheme)