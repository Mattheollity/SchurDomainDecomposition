#include "omp.h"
#include <catch2/catch_test_macros.hpp>

TEST_CASE("omp is active")
{
    #pragma omp parallel
    REQUIRE(omp_get_num_threads()!=1);
}