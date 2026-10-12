#ifndef SCHUR_SOLVER_
#define SCHUR_SOLVER_

#include "common.h"
#include "poisson.h"
#include "problem_conditions.h"
#include "time_scheme.h"

class schur_solver
{
    approxs_conds appc;
    grid_conds gc;
    time_integr_conds tic;

    grid * g;
    poisson_system * p_system;
    tsws * s_solver;

    public:
    schur_solver(const approxs_conds &, const grid_conds &, const time_integr_conds &);
    void reset_approxs_conds(const approxs_conds &);
    void reset_grid_conds(const grid_conds &);
    void reset_time_integr_conds(const time_integr_conds &);
    const vector<vect> solve();
    double compare_with_analytical_solution(const vector<vect> & solution);
    double compare_with_analytical_solution(double (*solution)(point, double));
    ~schur_solver();
};




#endif