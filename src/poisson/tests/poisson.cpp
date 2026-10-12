#include "poisson.h"
#include <catch2/catch_test_macros.hpp>
#include <vector>
#include "common.h"
#include "problem_conditions.h"
#include "approx.h"
#include "file_io.h"

TEST_CASE("single subdomain poisson system")
{
    grid_conds gc(
        {
            {0,0},
            {1,1}
        },
        { 10,10},
        {{
            {0,0},
            {1,1}
        }}
    );
    approxs_conds appc = {
        approx,
        dirich_boundary_conds
    };
    grid g(gc);
    poisson_system p_system(g, appc);
    const sparse_system & p_sys = p_system.get_system();

    const spmat & diag = p_sys.system[0][0];

    int TOTAL = 8*8;
    int fssdTOTAL = 8*7;
    int main_diag = 0;
    int first_sub_diags = 0;
    int second_sub_diags = 0;
    REQUIRE(diag.cols()==TOTAL);
    REQUIRE(diag.rows()==TOTAL);
    for(int i=0; i<TOTAL; i++)
    {
        if(diag.coeff(i, i)!=0) main_diag++;

        if(i<TOTAL-8)
        if(diag.coeff(i, 8+i)!=0) first_sub_diags++;

        if(i>=8)
        if(diag.coeff(i, i-8)!=0) second_sub_diags++;
    }
    REQUIRE(main_diag==TOTAL);

    REQUIRE(first_sub_diags==fssdTOTAL);
    REQUIRE(second_sub_diags==fssdTOTAL);
    
    
    

    write_matrix(diag, "test_file_one_sbd_p_system");
}

TEST_CASE("two subdomains poisson system main part")
{
    grid_conds gc(
        {
            {0,0,0},
            {1,1,1}
        },
        { 10,10,10},
        {
            {
                {0,0,0},
                {1,0.5,1},
            },
            {
                {0,0.5,0},
                {1,1,1}
            }
        }
    );

    approxs_conds appc = {
        approx,
        dirich_boundary_conds
    };
    grid g(gc);
    poisson_system p_system(g, appc);

    write_matrix(p_system.get_system().system[0][0], "test_file_two_sbd_p_system_main_part");

    REQUIRE(p_system.get_system().system[0][0].nonZeros()!=0);
}

TEST_CASE("permutation test for two subds: 2d case")
{
    grid_conds gc(
        {
            {0,0},
            {1,1}
        },
        { 10,10},
        {
            {
                {0,0},
                {1,0.5},
            },
            {
                {0,0.5},
                {1,1}
            }
        }
    );

    approxs_conds appc = {
        approx,
        dirich_boundary_conds
    };
    grid g(gc);
    poisson_system p_system(g, appc);
    
    spmat test_perm = p_system.get_system().permutation * p_system.get_system().permutation.transpose();

    mat I(8*8, 8*8);
    I.setZero();
    for(int i=0;i<8*8;i++)
    I.coeffRef(i,i)=1;
    REQUIRE(I==test_perm.toDense());
}

TEST_CASE("subd systems correct for two: 2d case")
{
    grid_conds gc(
        {
            {0,0},
            {1,1}
        },
        { 16,16},
        {
            {
                {0,0},
                {1,0.5},
            },
            {
                {0,0.5},
                {1,1}
            }
        }
    );

    approxs_conds appc = {
        approx,
        dirich_boundary_conds
    };
    grid g(gc);
    poisson_system p_system(g, appc);



    grid_conds gc_origin(
        {
            {0,0},
            {1,1}
        },
        { 16,16},
        {
            {
                {0,0},
                {1,1}
            }
        }
    );
    grid g_origin(gc_origin);
    poisson_system p_system_origin(g_origin, appc);



    mat L_origin = p_system_origin.get_system().system[0][0].toDense();
    mat L = L_origin;

    vector<int> csum = {0};
    for(int i : p_system.get_system().sizes)
    csum.push_back(csum.back()+i);

    for(int i=0;i<p_system.get_system().amount_of_parts; i++)
    for(int j=0;j<p_system.get_system().amount_of_parts; j++)
    L.block(csum[i], csum[j], csum[i+1]-csum[i], csum[j+1]-csum[j]) = p_system.get_system().system[i][j].toDense();

    write_matrix_struct(L_origin, "L_origin");
    write_matrix_struct(L, "L");

    L_origin = p_system.get_system().permutation * L_origin * p_system.get_system().permutation.transpose();

    write_matrix_struct(L_origin, "L_origin_after");

    REQUIRE(L==L_origin);
}

TEST_CASE("permutation test for four subds: 2d case")
{
    grid_conds gc(
        {
            {0,0},
            {1,1}
        },
        { 10,10},
        {
            {
                {0,0},
                {0.5,0.5},
            },
            {
                {0.5,0},
                {1,0.5}
            },
            {
                {0,0.5},
                {0.5,1},
            },
            {
                {0.5,0.5},
                {1,1}
            }
        }
    );

    approxs_conds appc = {
        approx,
        dirich_boundary_conds
    };
    grid g(gc);
    poisson_system p_system(g, appc);
    
    spmat test_perm = p_system.get_system().permutation * p_system.get_system().permutation.transpose();

    mat I(8*8, 8*8);
    I.setZero();
    for(int i=0;i<8*8;i++)
    I.coeffRef(i,i)=1;
    REQUIRE(I==test_perm.toDense());
}

TEST_CASE("subd systems correct for four: 2d case")
{
    grid_conds gc(
        {
            {0,0},
            {1,1}
        },
        { 16,16},
        {
            {
                {0,0},
                {0.5,0.5},
            },
            {
                {0.5,0},
                {1,0.5}
            },
            {
                {0,0.5},
                {0.5,1},
            },
            {
                {0.5,0.5},
                {1,1}
            }
        }
    );

    approxs_conds appc = {
        approx,
        dirich_boundary_conds
    };
    grid g(gc);
    poisson_system p_system(g, appc);



    grid_conds gc_origin(
        {
            {0,0},
            {1,1}
        },
        { 16,16},
        {
            {
                {0,0},
                {1,1}
            }
        }
    );
    grid g_origin(gc_origin);
    poisson_system p_system_origin(g_origin, appc);



    mat L_origin = p_system_origin.get_system().system[0][0].toDense();
    mat L = L_origin;

    vector<int> csum = {0};
    for(int i : p_system.get_system().sizes)
    csum.push_back(csum.back()+i);

    for(int i=0;i<p_system.get_system().amount_of_parts; i++)
    for(int j=0;j<p_system.get_system().amount_of_parts; j++)
    L.block(csum[i], csum[j], csum[i+1]-csum[i], csum[j+1]-csum[j]) = p_system.get_system().system[i][j].toDense();

    write_matrix_struct(L_origin, "L_origin");
    write_matrix_struct(L, "L");

    L_origin = p_system.get_system().permutation * L_origin * p_system.get_system().permutation.transpose();

    write_matrix_struct(L_origin, "L_origin_after");

    REQUIRE(L==L_origin);
}