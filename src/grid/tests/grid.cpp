#include "grid.h"
#include <catch2/catch_test_macros.hpp>
#include <fstream>
#include <vector>
#include "dim_iterator.h"
#include "problem_conditions.h"
#include "file_io.h"

TEST_CASE("grid size and symmetry")
{
    cout<<"grid: size and symmetry"<<endl;

    double eps = 1e-10;
    grid g1d = grid(load_grid_conds("grid1d"));

    auto vtx1b = g1d.get_boundary_points();
    auto vtx1i = g1d.get_inner_points();

    int size = vtx1i.size();
    REQUIRE(g1d.get_vertexes().size()==10);
    REQUIRE(vtx1b.size()==2);
    REQUIRE(vtx1b[0][0]==0);
    REQUIRE(vtx1b[1][0]==1);
    REQUIRE(vtx1i.size()==8);
    for(int i=0;i<size;i++)
    REQUIRE(abs(vtx1i[i][0]-(1.0 - vtx1i[size-i-1][0]))<eps);


    grid g1d_1 = grid(load_grid_conds("grid1d_1"));

    auto vtx2b = g1d_1.get_boundary_points();
    auto vtx2i = g1d_1.get_inner_points();

    size = vtx2i.size();
    REQUIRE(g1d_1.get_h()[0]==0.1);
    REQUIRE(g1d_1.get_vertexes().size()==11);
    REQUIRE(vtx2b.size()==2);
    REQUIRE(vtx2b[0][0]==0);
    REQUIRE(vtx2b[1][0]==1);
    REQUIRE(vtx2i.size()==9);
    for(int i=0;i<size;i++)
    REQUIRE(abs(vtx2i[i][0]-(1.0 - vtx2i[size-i-1][0]))<eps);
}

TEST_CASE("grid two subdomain. division")
{
    cout<<"grid: two subdomains. division"<<endl;
    double eps = 1e-10;
    grid g1d = grid(load_grid_conds("grid1d_2sbd_sep0.65"));

    auto sbd1 = g1d.get_subdomain()[0];
    auto sbd2 = g1d.get_subdomain()[1];

    REQUIRE(sbd1.boundary.size()==1);
    REQUIRE(sbd2.boundary.size()==1);
    REQUIRE(sbd1.inner.size()==4);
    REQUIRE(sbd2.inner.size()==3);
}

TEST_CASE("grid2d two subdomain. division")
{
    cout<<"grid2d two subdomain. division"<<endl;
    double eps = 1e-10;
    grid g2d = grid(load_grid_conds("grid2d_10_20_sep0.7"));

    auto sbd1 = g2d.get_subdomain()[0];
    auto sbd2 = g2d.get_subdomain()[1];

    REQUIRE(sbd1.boundary.size()==8);
    REQUIRE(sbd2.boundary.size()==8);
    REQUIRE(sbd1.inner.size()==8*11);
    REQUIRE(sbd2.inner.size()==8*5);
}

TEST_CASE("grid2d two subdomain. interaction")
{
    cout<<"grid2d two subdomain. interaction"<<endl;
    double eps = 1e-10;
    grid g2d = grid(load_grid_conds("grid2d_10_20_sep0.7"));

    auto io1 = g2d.get_subdomain_interaction(0);
    auto io2 = g2d.get_subdomain_interaction(1);
    
    for(int i=0;i<8*5;i++)
    REQUIRE(io2[0].inputs[i].size()>1);
}

TEST_CASE("grid2d two subdomains. Permutation")
{
    cout<<"grid2d two subdomains. Permutation"<<endl;
    grid g2d = grid(load_grid_conds("grid2d_10_20_sep0.7"));

    auto pm = g2d.get_permutation_along_subd();

    REQUIRE(pm.size()==8*18);
}

// TEST_CASE("grid file: one subdomain 3d")
// {
//     grid_conds gc(
//         {
//             {0,0,0},
//             {1,1,1}
//         },
//         { 10,10,10 },
//         {{
//             {0,0,0},
//             {1,1,1}
//         }}
//     );
    
//     grid g(gc);
    
//     // auto subdomains = g.get_subdomain();
//     // REQUIRE(g.get_boundary().size()==488);
//     // REQUIRE(subdomains.size()==1);
//     // REQUIRE(subdomains[0].inner.size()==512);
//     // REQUIRE(subdomains[0].boundary.size()==0);
// }

// TEST_CASE("grid file: two subdomain 3d")
// {
//     grid_conds gc(
//         {
//             {0,0,0},
//             {1,1,1}
//         },
//         { 10,10,10 },
//         {
//             {
//                 {0,0,0},
//                 {1,0.5,1},
//             },
//             {
//                 {0,0.5,0},
//                 {1,1,1}
//             }
//         }
//     );
    
//     grid g(gc);
    
//     auto subdomains = g.get_subdomain();
//     auto vtx = g.get_vertexes();
//     REQUIRE(subdomains.size()==2);
//     REQUIRE(g.get_boundary().size()==488);
//     REQUIRE(subdomains[0].inner.size()==subdomains[1].inner.size());
//     REQUIRE(subdomains[0].boundary.size()==subdomains[1].boundary.size());
// }


// TEST_CASE("grid. One sbd, neighbors")
// {
//     grid_conds gc(
//         {
//             {0,0,0},
//             {1,1,1}
//         },
//         { 10,10,10 },
//         {{
//             {0,0,0},
//             {1,1,1}
//         }}
//     );
//     grid g(gc);

    //interaction io = g.get_subdomains_interaction(0, 0);

    // REQUIRE(io.outputs.size()<=8*8*8);
    // for(
    //     dim_iterator di({
    //         {2,2,2},
    //         {8,8,8}
    //     },
    //         {2,2,2}
    //     );
    //     !di.is_full();
    //     ++di
    // )
    // {

    // for(auto vtx : io.inputs)
    // REQUIRE
    // (
    //     io.start_input_index[di.conv({{1,1,1},{9,9,9}})+1]
    //     -
    //     io.start_input_index[
    //     di.conv({{1,1,1},{9,9,9}})
    //     ]==7
    // );
    // }
//}

// TEST_CASE("grid. 2 subd. one neighbor at boundary")
// {
//     grid_conds gc(
//         {
//             {0,0,0},
//             {1,1,1}
//         },
//         { 10,10,10 },
//         {
//             {
//                 {0,0,0},
//                 {1,0.5,1},
//             },
//             {
//                 {0,0.5,0},
//                 {1,1,1}
//             }
//         }
//     );
//     grid g(gc);

//     // interaction io = g.get_subdomains_interaction(0, 1);

//     // for(int i=0;i<io.outputs.size();i++)
//     // REQUIRE(io.start_input_index[i+1]-io.start_input_index[i]==1);
// }