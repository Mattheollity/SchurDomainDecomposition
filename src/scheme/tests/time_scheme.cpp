#include "time_scheme.h"
#include "common.h"
#include "file_io.h"
#include "poisson.h"
#include "problem_conditions.h"
#include <catch2/catch_test_macros.hpp>
#include <cstdlib>
#include <fstream>
#include <vector>
#include "approx.h"


double solution(point p, double t)
{
    return exp(-3.0*M_PI*M_PI*t)*trip_sin(p);
}

const double ppi = 3.1415926535;

TEST_CASE("1d time scheme test: theta = 0")
{
    double eps = 0.5;
    double time = 0;
    vector<double> sol_at_time;
    auto gc = load_grid_conds("grid1d");
    grid g(gc);
    auto tic = load_tic("tic_theta_0_10");
    tic.set_zero_time_data(trip_sin, g.get_inner_points());
    approxs_conds appc({approx, dirich_boundary_conds});
    poisson_system p_system(g, appc);
    tsws s_solver(p_system.get_system(), tic);
    auto points = g.get_inner_points();

    s_solver.solve();

    for(int i=0;i<tic.N;i++)
    {
        time = tic.time.first+i*tic.T;
        vect csol = s_solver.get_solution_at_time(i);

        REQUIRE(csol.rows()==points.size());
        for(int j=0;j<points.size(); j++)
        REQUIRE(std::abs
        (
            csol.coeffRef(j)
            -
            solution(points[j], time)
        )<eps);
    }
}

// TEST_CASE("single subdomain time scheme")
// {
//     grid_conds gc1(
//         {
//             {0,0,0},
//             {1,1,1}
//         },
//         {10, 10, 10},
//         {{
//             {0,0,0},
//             {1,1,1}
//         }}
//     );
//     grid_conds gc2(
//         {
//             {0,0,0},
//             {1,1,1}
//         },
//         { 20,20,20},
//         {
//             {
//                 {0,0,0},
//                 {0.3,1,1}
//             },
//             {
//                 {0.3,0,0},
//                 {0.6,1,1}
//             },
//             {
//                 {0.6,0,0},
//                 {1,1,1}
//             }
//         }
//     );
//     grid_conds gc3(
//         {
//             {0,0,0},
//             {1,1,1}
//         },
//         { 30,30,30},
//         {
//             {
//                 {0,0,0},
//                 {0.3,1,1}
//             },
//             {
//                 {0.3,0,0},
//                 {0.6,1,1}
//             },
//             {
//                 {0.6,0,0},
//                 {1,1,1}
//             }
//         }
//     );

//     grid_conds gc4(
//         {
//             {0,0,0},
//             {1,1,1}
//         },
//         { 40,40,40},
//         {
//             {
//                 {0,0,0},
//                 {0.3,1,1}
//             },
//             {
//                 {0.3,0,0},
//                 {0.6,1,1}
//             },
//             {
//                 {0.6,0,0},
//                 {1,1,1}
//             }
//         }
//     );

//     approxs_conds appc = 
//     {
//         3,
//         approx,
//         dirich_boundary_conds
//     };
//     time_integr_conds tic(
//         0,
//         0.000000001,
//         10,
//         1,
//         trip_sin
//     );

//     grid g1(gc1);
//     grid g2(gc2);
//     grid g3(gc3);
//     grid g4(gc4);

//     tic.build_y0(g1);

//     poisson_system p_system1(g1, appc);

//     tsws s_solver1(p_system1.get_system(), tic);
//     auto solution1 = s_solver1.get_solution();

//     tic.build_y0(g2);
//     poisson_system p_system2(g2, appc);
    
//     tsws s_solver2(p_system2.get_system(), tic);
//     auto solution2 = s_solver2.get_solution();

//     double L12 = L2_diffrence(solution1, solution2);
//     cout<<"L12: "<<L12<<endl;


//     tic.build_y0(g3);
//     poisson_system p_system3(g3, appc);
    
//     tsws s_solver3(p_system3.get_system(), tic);
//     auto solution3 = s_solver3.get_solution();

//     double L23 = L2_diffrence(solution2, solution3);
//     cout<<"L23: "<<L23<<endl;

//     tic.build_y0(g4);

//     poisson_system p_system4(g4, appc);

//     tsws s_solver4(p_system4.get_system(), tic);
//     auto solution4 = s_solver4.get_solution();

//     double L34 = L2_diffrence(solution3, solution4);
//     cout<<"L34: "<<L34<<endl;
// }

// TEST_CASE("two subdomains time scheme")
// {
//     grid_conds gc(
//         {
//             {0,0,0},
//             {1,1,1}
//         },
//         { 10,10,10 },
//         {{
//             {0,0,0},
//             {1,0.5,1}
//         }, {
//             {0,0.5,0},
//             {1,1,1}
//         }}
//     );
//     grid g(gc.bbox, gc.N);
//     g.create_subdomain(gc.subdomains[0]);
//     poisson_solver p_solver(g, approx, dirich_boundary_conds);
//     sparse_system p_system = p_solver.get_system();
//     vect y0_;
//     y0(g, y0_);
//     time_integr_conds tic(
//         {0,1},
//         10,
//         0.6,
//         g.get_inner().size(),
//         y0_
//     );
//     tsws s_solver(p_system, tic);

//     auto solution = s_solver.get_solution();

//     auto boundary = g.get_boundary();
//     auto inner = g.get_inner();
//     REQUIRE(solution.size()==tic.N);

//     std::ofstream file("file2.txt");
//     for(auto s : solution)
//     {
//         for(auto i : s)
//         file<<i<<endl;
//         file<<"---"<<endl;
//     }
// }