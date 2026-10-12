#include "problem_conditions.h"

time_integr_conds::time_integr_conds(pair<double, double> time_, double T_, double theta_, double (*y0_)(point), const vector<point> & grid) :
    time(time_),
    T(T_),
    theta(theta_),
    y0(y0_),
    zero_time_data(false),
    write_each(1)
{
    N = (time.second-time.first)/T;
    build_y0(grid);
}

time_integr_conds::time_integr_conds(pair<double, double> time_, int N_, double theta_, double (*y0_)(point), const vector<point> & grid) :
    time(time_),
    N(N_),
    theta(theta_),
    y0(y0_),
    zero_time_data(false),
    write_each(1)
{
    T = (time.second-time.first)/N;
    build_y0(grid);
}

time_integr_conds::time_integr_conds(double start_time, double T_, int N_, double theta_, double (*y0_)(point), const vector<point> & grid) :
    time({start_time, start_time+T_*N_}),
    N(N_),
    T(T_),
    theta(theta_),
    y0(y0_),
    zero_time_data(true),
    write_each(1)
{
    build_y0(grid);
}

time_integr_conds::time_integr_conds(double start_time, double T_, int N_, double theta_) :
    time({start_time, start_time+T_*N_}),
    N(N_),
    T(T_),
    theta(theta_),
    zero_time_data(false),
    write_each(1)
{}

void time_integr_conds::build_y0(const vector<point> & points)
{
    y0_vect = vect(points.size());
    
    for(int i=0;i<y0_vect.rows();i++)
    y0_vect[i]=y0(points[i]);
}

void time_integr_conds::set_zero_time_data(double (*y0_)(point), const vector<point> & grid)
{
    y0 = y0_;
    build_y0(grid);
}