cmake_minimum_required(VERSION 3.14)
project(PoissonDecomposition LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)

# set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -DNDEBUG")
# set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g")

# if(MSVC)
#     set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /arch:AVX2")
# endif()

find_package(Eigen3 REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Catch2 REQUIRED)

add_subdirectory(src/grid)
add_subdirectory(src/poisson)
add_subdirectory(src/scheme)
add_subdirectory(src/utils)

add_executable(compt
    src/main.cpp
)

target_link_libraries(
    compt
    grid
    poisson
    time_scheme
    utils
)