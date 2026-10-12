#include "common.h"
#include "approx.h"
#include "poisson.h"
#include "problem_conditions.h"
#include "file_io.h"

vect get_res_at_time_step(const grid &, const time_integr_conds &, const sparse_system &, int N);
vector<vector<pair<double, double>>> check_laplace_approx(string filename);