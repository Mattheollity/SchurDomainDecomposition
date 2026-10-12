#ifndef DIM_ITERATOR_
#define DIM_ITERATOR_

#include "common.h"
#include "problem_conditions.h"

struct dim_iterator
{
    unsigned dim;
    pair<voxel_point, voxel_point> bbox;
    voxel_point voxel_coord;

    voxel_point conv_sizes;

    dim_iterator(const pair<voxel_point, voxel_point> & bbox, const voxel_point & initial_position);

    dim_iterator & operator++();
    dim_iterator & operator--();
    
    bool operator<(const voxel_point & v);
    bool operator!=(const voxel_point & v);

    bool is_full();
    bool is_empty();

    int conv();
    int conv(const pair<voxel_point, voxel_point> & bbox);

    vector<pair<int, int>> get_neighbors();
    vector<pair<int, int>> get_neighbors(const pair<voxel_point, voxel_point> & bbox);
};

#endif