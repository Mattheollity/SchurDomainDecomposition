#ifndef GRID_
#define GRID_

#include "common.h"
#include "vertex.h"
#include "dim_iterator.h"
#include "problem_conditions.h"
#include <ostream>
#include <sstream>

struct interaction
{
    vector<vertex> outputs;
    vector<list<pair<int, vertex>>> inputs;

    pair<int,int> io_size;
    vector<double> h;
};

struct subdomain
{
    unsigned i;
    pair<voxel_point, voxel_point> N_bbox;
    vector<int> boundary;
    vector<int> inner;
    set<int> neighboring_subdomains;
};

class grid
{
    unsigned dim;
    pair<point, point> bbox;
    pair<voxel_point, voxel_point> N_bbox;

    vector<double> h;

    vector<vertex> vertexes;
    vector<int> boundary;
    vector<int> inner;

    vector<subdomain> subdomains;

    vector<int> permut;

    void fill_domain();
    void find_b_i();
    void find_b_i(int i);

    void to_inner_or_to_boundary(vertex &);
    void to_inner_or_to_boundary(subdomain &, const vertex &);

    pair<voxel_point, voxel_point> truncate_subdomain_bbox(pair<point, point> &);
    void fill_subdomain(subdomain &);
    void find_neighboring_subdomains();

    void load_subdomains(std::ifstream &);
    void load_subdomain(std::ifstream &, int i);
    void load_bboxes(std::ifstream &);
    void load_vertexes(std::ifstream &);

    void write_subdomains(std::ofstream &) const;
    void write_subdomain(std::ofstream &, int i) const;
    void write_bboxes(std::ofstream &) const;
    void write_vertexes(std::ofstream &) const;

    void find_neighbors_in_part(const vertex & v, list<pair<int, vertex>> & io_input, const vector<int> & input) const;
    int find_index_in_part(int, const vector<int> &) const;

    void assert_all_vertexes_has_subdomain() const;

    void setup_subdomains();
    void setup_subdomain(int i);

    vector<int> build_permutation_along_subd();

    std::stringstream gridlog;

    public:

    grid(const grid_conds &);
    grid(string);

    void create_subdomain(int i, pair<point, point> bbox_);

    void fill_interaction_between(interaction &, const vector<int> & input, const vector<int> & output) const;
    vector<interaction> get_subdomains_interaction(int i, int j) const;
    vector<interaction> get_subdomain_interaction(int i) const;
    interaction get_inner_boundary_transpose_interaction() const;


    point to_real_coord(voxel_point) const;
    voxel_point to_voxel_coord(point) const;

    void write_to_file(string) const;
    void load_from_file(string);

    const vector<subdomain> & get_subdomain() const {return subdomains;}
    const vector<vertex> & get_vertexes() const {return vertexes;}
    const vector<int> & get_boundary() const {return boundary;}
    const vector<int> & get_inner() const {return inner;}
    const vector<double> & get_h() const {return h;}
    const pair<voxel_point, voxel_point> & get_vbbox() const {return N_bbox;}
    vector<point> get_boundary_points() const;
    vector<point> get_inner_points() const;
    vector<int> get_permutation_along_subd() const {return permut;}
    int get_dim() const {return dim;}

    string get_grid_log() const {return gridlog.str();}
};


#endif
