#include "schur_solver.h"
#include "common.h"
#include "grid.h"
#include "poisson.h"
#include "problem_conditions.h"
#include "time_scheme.h"
#include "vertex.h"
#include <vector>

schur_solver::schur_solver(const approxs_conds & appc_, const grid_conds & gc_, const time_integr_conds & tic_) :
    appc(appc_),
    gc(gc_),
    tic(tic_),
    g(nullptr),
    p_system(nullptr),
    s_solver(nullptr)
{
    reset_grid_conds(gc);
}

void schur_solver::reset_grid_conds(const grid_conds & gc_)
{
    gc = gc_;

    if(g!=nullptr) delete g;
    g = new grid(gc);

    reset_approxs_conds(appc);
}

void schur_solver::reset_approxs_conds(const approxs_conds & appc_)
{
    appc = appc_;

    if(g==nullptr) throw simple_error("can`t create poisson system without mesh");
    if(p_system!=nullptr) delete p_system;
    p_system = new poisson_system(*g, appc);

    reset_time_integr_conds(tic);
}

void schur_solver::reset_time_integr_conds(const time_integr_conds & tic_)
{
    tic = tic_;
    if(p_system==nullptr) throw simple_error("can`t create tsws without poisson system");
    if(s_solver!=nullptr) delete s_solver;

    s_solver = new tsws(p_system->get_system(), tic);
}

const vector<vect> schur_solver::solve()
{
    if(s_solver==nullptr) throw simple_error("can`t solve withour tsws");
    return s_solver->solve();
}

double schur_solver::compare_with_analytical_solution(const vector<vect> & asol)
{
    double diffr = 0;
    vect difv(asol[0]);

    if(s_solver==nullptr) throw simple_error("can`t compare without tsws");
    const vector<vect> & sol = s_solver->get_solution();

    if(sol.size()<2) s_solver->solve();
    if(sol.size()!=asol.size()) throw simple_error("solution and analytical solution must have same sizes");

    for(int i=0; i<sol.size(); i++)
    {
        for(int j=0; j<difv.size(); j++)
        difv.coeffRef(j)=sol[i][j]-asol[i][j];

        double n_diffr = difv.norm();
        if(n_diffr>diffr) diffr=n_diffr;
    }

    return diffr;
}

double schur_solver::compare_with_analytical_solution(double (*solution)(point p, double t))
{
    if(g==nullptr) throw simple_error("can`t compare without grid");
    vector<point> points = g->get_inner_points();

    if(s_solver==nullptr) throw simple_error("can`t compare without s_solver");
    const vector<vect> & comp_solution = s_solver->get_solution();

    double max_error = 0;

    for(int i=0; i<comp_solution.size(); i++)
    for(int j=0; j<points.size(); j++)
    {
        double y = solution(points[j], tic.time.first+tic.T*i);
        double cy = comp_solution[i][j];
        double new_error = (cy-y)*(cy-y);
        if(new_error>max_error) max_error = new_error;
    }

    return max_error;
}

schur_solver::~schur_solver()
{
    if(g!=nullptr) delete g;
    if(p_system!=nullptr) delete p_system;
    if(s_solver!=nullptr) delete s_solver;
}