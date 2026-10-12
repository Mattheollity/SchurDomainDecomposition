#ifndef IAS_
#define IAS_

#include "common.h"
#include "poisson.h"

double trip_sin(const point p);
double real_solution(int dim, const point p, double t);
vect y0(const grid & gr, vect & v);

double approx(vector<double> h, const vertex & main, const vertex & neighbor);
double dirich_boundary_conds(point, double t);
vect get_solution_at_time_step(const grid &, const time_integr_conds & tic, int N);
vect get_time_deriv_at_time_step(const grid &, const time_integr_conds & tic, int N);
vect get_dummy_time_deriv_at_time_step(const grid &, const time_integr_conds & tic, int N);
vect get_laplace_at_time_step(const grid &, const time_integr_conds & tic, const sparse_system & p_system, int N);
double get_laplace_eqlid_at_time_step(const grid &, const time_integr_conds & tic, const sparse_system & p_system, int N);

#endif