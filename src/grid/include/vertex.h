#ifndef VERTEX_
#define VERTEX_

#include "common.h"
#include "problem_conditions.h"

struct vertex
{
    unsigned dim;

    int i;
    int containing_subdomain;
    bool boundary;
    
    voxel_point voxel_coord;

    //neighbor index, neighbor dimension
    vector<pair<int, int>> neighbors;

    vertex(int i, voxel_point vp);
    vertex() {};

    void from_file(std::ifstream &);

    bool at_top_boundary(voxel_point vp) const;
    bool at_bottom_boundary(voxel_point vp) const;
    bool at_critical_point(pair<voxel_point, voxel_point> N_bbox) const;

    string to_string();
};

#endif