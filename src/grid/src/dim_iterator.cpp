#include "dim_iterator.h"

dim_iterator::dim_iterator(const pair<voxel_point, voxel_point> & bbox_, const voxel_point & initial_position_) : 
    bbox(bbox_),
    dim(bbox_.first.size()),
    voxel_coord(initial_position_),
    conv_sizes(dim, 1)
{
    for(unsigned d=1;d<dim;d++)
        conv_sizes[d]=conv_sizes[d-1]*(bbox.second[d-1]-bbox.first[d-1]);
}

dim_iterator & dim_iterator::operator++() {
    voxel_coord[0]++;
    for(unsigned d=0; d<dim-1;d++)
    {
        if(voxel_coord[d]==bbox.second[d])
        {
            voxel_coord[d]=bbox.first[d];
            voxel_coord[d+1]++;
        }
    }
    return *this;
}

dim_iterator & dim_iterator::operator--() {
    voxel_coord[0]--;
    for(unsigned d=0; d<dim-1;d++)
    {
        if(voxel_coord[d]==-1)
        {
            voxel_coord[d]=bbox.second[d]-1;
            voxel_coord[d+1]--;
        }
    }
    return *this;
}


bool dim_iterator::operator<(const voxel_point & v)
{
    for(unsigned d=0;d<dim;d++) if(voxel_coord[d]>=v[d]) return false;
    return true;
}

bool dim_iterator::operator!=(const voxel_point & v)
{
    for(unsigned d=0;d<dim;d++) if(voxel_coord[d]!=v[d]) return true;
    return false;
}

bool dim_iterator::is_full()
{
    return voxel_coord[dim-1]>=bbox.second[dim-1];
}

bool dim_iterator::is_empty()
{
    return voxel_coord[dim-1]<bbox.first[dim-1];
}

int dim_iterator::conv()
{
    unsigned s=0;

    for(unsigned d=0;d<dim;d++)
    s+=conv_sizes[d]*voxel_coord[d];

    return s;
}

int dim_iterator::conv(const pair<voxel_point, voxel_point> & bbox_)
{
    unsigned s=0;
    vector<unsigned> local_conv_sizes (dim,1);

    for(unsigned d=1;d<dim;d++)
    local_conv_sizes[d]=local_conv_sizes[d-1]*(bbox_.second[d-1]-bbox_.first[d-1]);

    for(unsigned d=0;d<dim;d++)
    s+=local_conv_sizes[d]*(voxel_coord[d]-bbox_.first[d]);

    return s;
}

vector<pair<int, int>> dim_iterator::get_neighbors()
{
    vector<pair<int, int>> neighbors;

    int s = conv();

    for(unsigned d=0;d<dim;d++)
    {
        if(voxel_coord[d]+1<bbox.second[d]) neighbors.push_back({s+conv_sizes[d], d});
        if(voxel_coord[d]-1>=bbox.first[d]) neighbors.push_back({s-conv_sizes[d], d});
    }

    return neighbors;
}

vector<pair<int, int>> dim_iterator::get_neighbors(const pair<voxel_point, voxel_point> & bbox_)
{
    vector<pair<int, int>> neighbors;

    int s=0;
    vector<unsigned> local_conv_sizes (dim,1);

    for(unsigned d=1;d<dim;d++)
    local_conv_sizes[d]=conv_sizes[d-1]*(bbox_.second[d]-bbox_.first[d]);

    for(unsigned d=0;d<dim;d++)
    s+=local_conv_sizes[d]*(voxel_coord[d]-bbox_.first[d]);

    for(unsigned d=0;d<dim;d++)
    {
        if(voxel_coord[d]+1<bbox_.second[d]) neighbors.push_back({s+local_conv_sizes[d], d});
        if(voxel_coord[d]-1>=bbox_.first[d]) neighbors.push_back({s-local_conv_sizes[d], d});
    }

    return neighbors;
}