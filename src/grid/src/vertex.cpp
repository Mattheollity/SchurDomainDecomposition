#include "vertex.h"

vertex::vertex(int i_, voxel_point vp) :
    i(i_),
    dim(vp.size()),
    voxel_coord(vp),
    containing_subdomain(-1),
    boundary(false) {};

void vertex::from_file(std::ifstream & file)
{
    int n_size;
    file>>dim;
    file>>i;
    file>>containing_subdomain;
    file>>n_size;
    neighbors = vector<pair<int, int>> (n_size);
    voxel_coord = voxel_point(dim);
    for(int i=0;i<dim;i++) file>>voxel_coord[i];
    for(int i=0;i<n_size;i++) file>>neighbors[i].first>>neighbors[i].second;
}

bool vertex::at_top_boundary(voxel_point bbox_top_right) const
{
    for(unsigned d=0;d<dim;d++) if(voxel_coord[d]==bbox_top_right[d]-1) return true;
    return false;
}

bool vertex::at_bottom_boundary(voxel_point bbox_bottom_left) const {
    for(unsigned d=0;d<dim;d++) if(voxel_coord[d]==bbox_bottom_left[d]) return true;
    return false;
}

string vertex::to_string()
{
    std::stringstream oss;
    oss<<dim<<" "<<i<<" "<<containing_subdomain<<" "<<neighbors.size()<<"\n";
    for(auto vc : voxel_coord) oss<<vc<<"\n";
    for(auto n : neighbors) oss<<n.first<<" "<<n.second<<"\n";
    return oss.str();
}

bool vertex::at_critical_point(pair<voxel_point, voxel_point> N_bbox) const
{
    int dim_counter = 0;
    for(int d=0;d<dim;d++)
    if(voxel_coord[d]==N_bbox.second[d]-1 || voxel_coord[d]==N_bbox.first[d]) dim_counter++;

    if(dim_counter==dim) return true;
    else return false;
}