#include "common.h"
#include "approx.h"
#include "poisson.h"
#include "problem_conditions.h"
#include <vector>
#include "eexpokit/expv.hpp"

double get_cheb_norm(const vect & v);
double get_eqlid_norm(const vect & v);
vector<double> get_eqlid_converg_at_time_step(const grid &, const time_integr_conds &, const sparse_system &, int N, string filename);
vector<vector<double>> check_space_converg(string filename);
pair<vector<string>, vector<int>> get_check_data(string filename);
vect load_vect(string filename, int N);
vect multiply_spsys_and_vect(const sparse_system & system, const vect & v, const vect & vn, const time_integr_conds & tic);
double get_eqlid_L_converg(const grid & g, const time_integr_conds & tic, const sparse_system & p_system, int N, string dir_name);
vect get_approx_solution_at_time_step(const grid & g, const time_integr_conds & tic, int N, const sparse_system & spsys);