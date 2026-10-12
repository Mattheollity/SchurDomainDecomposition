#ifndef POISSON_
#define POISSON_

#include "grid.h"
#include "common.h"
#include "problem_conditions.h"
#include <sstream>

typedef vector<vector<spmat>> spsys;

struct approxs_conds
{
    double (*approx)(const vector<double> h, const vertex &, const vertex &);
    double (*dirichlet)(point, double t);
};

struct sparse_system
{
    spsys system;
    spmat right_addition;
    spmat permutation;
    int amount_of_parts;
    vector<int> sizes;
    vector<point> boundary_points;
    double (*dirichlet_conds)(point, double t);
    double (*approx)(const vector<double> h, const vertex &, const vertex &);
    vect boundary_at_time(double t) const;
};

class poisson_system
{
    sparse_system p_system;

    void build_main_diagonal_system(const vector<interaction> &, int i);
    void build_interaction_system(const vector<interaction> &, int i, int j);
    void fill_spmat(spmat &, const interaction &);

    void write_spmat(std::ofstream &, spmat &);
    void load_spmat(std::ifstream &, spmat &);

    void build_right_addition(const grid & g);
    void add_zeros_matrixes();
    void build_permutation_matrix(const grid & g);

    std::stringstream poissonlog;

    public:

    poisson_system(const grid &, const approxs_conds &);
    poisson_system(string s) {load_from_file(s);}

    const sparse_system & build_system(const grid &);
    const sparse_system & get_system() const {return p_system;}
    const vector<int> & get_system_sizes() const;
    void write_to_file(string);
    void load_from_file(string);

    vect get_diriclet_boundary(const grid & g);

    string get_poisson_log() const {return poissonlog.str();}
    
};

#endif