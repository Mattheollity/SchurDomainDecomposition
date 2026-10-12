#ifndef PROBLEN_CONDITIONS_
#define PROBLEN_CONDITIONS_

#include "common.h"

typedef vector<double> point;
typedef vector<int> voxel_point;
class grid;

struct grid_conds
{
    pair<point, point> bbox;
    voxel_point N;
    vector<pair<point, point>> subdomains;

    grid_conds(
        pair<point, point> bbox_,
        voxel_point N_,
        vector<pair<point, point>> subdomains_
    ) : bbox(bbox_), N(N_), subdomains(subdomains_) {};
};

struct time_integr_conds
{
    double T;
    double theta;
    int N;
    pair<double, double> time;
    int inner_size;

    bool zero_time_data;
    double (*y0)(point);
    vect y0_vect;

    time_integr_conds(double start_time, double T, int N, double theta);
    time_integr_conds(pair<double, double> time_, double T, double theta, double (*y0)(point), const vector<point> &);
    time_integr_conds(pair<double, double> time_, int N, double theta, double (*y0_)(point), const vector<point> &);
    time_integr_conds(double start_time, double T, int N, double theta, double (*y0_)(point), const vector<point> &);

    void build_y0(const vector<point> &);
    void set_zero_time_data(double (*y0)(point), const vector<point> &);
    int write_each;
};

#endif