#ifndef TIME_SCHEME_
#define TIME_SCHEME_

#include "common.h"
#include "poisson.h"
#include "vertex.h"
#include <utility>
#include <vector>
#include "problem_conditions.h"


struct one_step_conds
{
    vector<vect> y0;
    vector<vect> y;
    vector<vect> right;

    vector<int> parts_size;
    spsys left_system;
    spsys right_system;
};

struct schur_system
{
    vector<int> block_sizes;
    spmat schur_mat;
    vect y;
    vect right;
};

//time scheme with schur
class tsws
{
    time_integr_conds tic;
    sparse_system origin_system;

    schur_system ss;
    one_step_conds osc;
    
    
    vector<vect> solution;

    void expimp_matrix_manipulation(spsys & target, const spsys & origin, double coeff);

    vector<vect> split_vect(const vect & origin);

    void concate_vect(vect & target, const vector<vect> & origin);

    const schur_system & build_schur_matrix();

    const one_step_conds & set_osc();
    const schur_system & set_ss();

    const one_step_conds & solve_osc();
    const schur_system & solve_ss();

    const one_step_conds & update_osc();
    const schur_system & update_ss();

    void build_right_part_osc();

    double get_current_time() const;
    


    public:
    
    tsws(const sparse_system &, const time_integr_conds &);
    const vector<vect> & solve();

    const schur_system & get_schur_system() const {return ss;}
    const time_integr_conds & get_tic() const {return tic;}
    const one_step_conds & get_osc() const {return osc;}
    const vector<vect> & get_solution() const {return solution;}
    const vect get_solution_at_time(int step) const;
    void write_file() const;

    
};

#endif