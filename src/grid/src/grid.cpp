#include "grid.h"
#include <omp.h>

grid::grid(const grid_conds & gc) :
    bbox(gc.bbox),
    dim(gc.N.size()),
    N_bbox({
        vector<int>(dim,0),
        gc.N
    }),
    h(dim),
    subdomains(gc.subdomains.size()),
    gridlog(std::stringstream())
{
    double t0 = omp_get_wtime();
    cout<<"GRID CREATION"<<endl;
    for(int i=0; i<dim; i++)
    h[i] = (bbox.second[i] - bbox.first[i])/(N_bbox.second[i]-1);

    int TOTAL = 1;
    for(int i=0; i<dim; i++) TOTAL*=N_bbox.second[i];
    vertexes = vector<vertex>(TOTAL);

    double t1 = omp_get_wtime();
    fill_domain();

    double t2 = omp_get_wtime();
    cout<<"grid was filled by "<<vertexes.size()<<" vertexes with time: "<<t2-t1<<endl;

    t1=t2;
    for(int i=0; i<gc.subdomains.size(); i++)
    create_subdomain(i, gc.subdomains[i]);

    setup_subdomains();
    t2=omp_get_wtime();
    cout<<subdomains.size()<<" subdomains were created with time: "<<t2-t1<<endl;

    assert_all_vertexes_has_subdomain();

    permut = build_permutation_along_subd();

    cout<<"GRID WAS SUCCESSFULLY CREATED with time: "<<omp_get_wtime()-t0<<endl<<endl;
}

void grid::fill_domain()
{
    int i=0;
    for
    (
        dim_iterator di(N_bbox, N_bbox.first);
        !di.is_full();
        ++di
    )
    {
        vertexes[i] = {i, di.voxel_coord};

        vertex & v = vertexes[i];
        to_inner_or_to_boundary(v);
        v.neighbors = di.get_neighbors();

        i++;
    }
}

grid::grid(string s)
{
    load_from_file(s);
}

void grid::load_from_file(string s)
{
    std::ifstream file(s);
    file>>dim;

    load_subdomains(file);
    load_bboxes(file);
    load_vertexes(file);

    find_b_i();
}

void grid::load_subdomains(std::ifstream & file)
{
    int amount_subd;
    file>>amount_subd;
    subdomains = vector<subdomain> (amount_subd);
    for(int i=0; i<amount_subd; i++) load_subdomain(file, i);
}

void grid::load_subdomain(std::ifstream & file, int i)
{
    file>>subdomains[i].i;
    subdomain & s = subdomains[i];

    s.N_bbox = {vector<int> (dim), vector<int> (dim)};
    for(int d=0;d<dim;d++)
        file>>s.N_bbox.first[d];
    for(int d=0;d<dim;d++)
        file>>s.N_bbox.second[d];

    int neigh_size;
    file>>neigh_size;
    s.neighboring_subdomains = set<int> ();
    for(int j=0;j<neigh_size;j++) 
    {
        unsigned k;
        file>>k;
        s.neighboring_subdomains.insert(k);
    }
}

void grid::load_bboxes(std::ifstream & file)
{
    bbox.first = vector<double> (dim);
    bbox.second = vector<double> (dim);
    h = vector<double> (dim);

    N_bbox.first = vector<int> (dim);
    N_bbox.second = vector<int> (dim);
    for(int i=0;i<dim;i++) file>>bbox.first[i];
    for(int i=0;i<dim;i++) file>>bbox.second[i];
    for(int i=0;i<dim;i++) file>>N_bbox.first[i];
    for(int i=0;i<dim;i++) file>>N_bbox.second[i];
    for(int i=0;i<dim;i++) file>>h[i];
}

void grid::load_vertexes(std::ifstream & file)
{
    int vtx_size;
    file>>vtx_size;
    vertexes = vector<vertex> (vtx_size);
    for(int i=0;i<vtx_size;i++) vertexes[i].from_file(file);
}

void grid::find_b_i()
{
    for(int i=0;i<vertexes.size(); i++) find_b_i(i);
}

void grid::find_b_i(int i)
{
    vertex & v = vertexes[i];
    to_inner_or_to_boundary(v);
    if(v.containing_subdomain!=-1)
        to_inner_or_to_boundary(subdomains[v.containing_subdomain], v);
}

void grid::write_to_file(string s) const
{
    std::ofstream file(s);
    file<<dim<<endl;

    write_subdomains(file);
    write_bboxes(file);
    write_vertexes(file);

    file.close();
}

void grid::write_subdomains(std::ofstream & file) const
{
    file<<subdomains.size()<<endl;
    for(int i=0; i<subdomains.size(); i++)
    write_subdomain(file, i);
}

void grid::write_subdomain(std::ofstream & file, int i) const
{
    const subdomain & s = subdomains[i];
    file<<s.i<<endl;

    for(int d=0;d<dim;d++)
    file<<s.N_bbox.first[d]<<endl;

    for(int d=0;d<dim;d++)
    file<<s.N_bbox.second[d]<<endl;

    file<<s.neighboring_subdomains.size()<<endl;

    for(auto n : s.neighboring_subdomains)
    file<<n<<endl;
}

void grid::write_bboxes(std::ofstream & file) const
{
    for(auto vc : bbox.first) file<<vc<<endl;
    for(auto vc : bbox.second) file<<vc<<endl;
    for(auto vc : N_bbox.first) file<<vc<<endl;
    for(auto vc : N_bbox.second) file<<vc<<endl;
    for(auto h_ : h) file<<h_<<endl;
}

void grid::write_vertexes(std::ofstream & file) const
{
    file<<vertexes.size()<<endl;

    for(auto v : vertexes)
    file<<v.to_string();
}

void grid::create_subdomain(int i, pair<point, point> bbox_)
{
    subdomains[i] = subdomain();
    subdomains[i].i = i;
    subdomains[i].N_bbox = truncate_subdomain_bbox(bbox_);

    fill_subdomain(subdomains[i]);
}

pair<voxel_point, voxel_point> grid::truncate_subdomain_bbox(pair<point, point> & bbox_)
{
    point & a = bbox_.first;
    point & b = bbox_.second;

    voxel_point va(dim);
    voxel_point vb(dim);

    for(unsigned d=0;d<dim;d++)
    {
        if(a[d]<bbox.first[d]) a[d]=bbox.first[d];
        if(b[d]>bbox.second[d]) b[d]=bbox.second[d];

        double eps = h[d]/2;

        int A = (a[d]-bbox.first[d])/h[d];
        if( a[d] - A*h[d] < eps) va[d] = A;
        else va[d] = A+1;

        int B = (b[d]-bbox.first[d])/h[d];
        if( b[d] - B*h[d] < eps) vb[d] = B;
        else vb[d] = B+1;

        if(va[d]==N_bbox.first[d]) va[d]+=1;
        if(vb[d]==N_bbox.second[d]) vb[d]-=1;
    }

    return {va, vb};
}

void grid::fill_subdomain(subdomain & s)
{
    for
    (
        dim_iterator di(s.N_bbox, s.N_bbox.first);
        !di.is_full();
        ++di
    )
    {
        vertex & v = vertexes[di.conv(N_bbox)];

        if(v.containing_subdomain!=-1) throw simple_error("SUBDOMAINS CAN`T HAVE COMMON VERTEXES");
        v.containing_subdomain = s.i;
    } 
}

void grid::to_inner_or_to_boundary(vertex & v)
{
    if(v.at_bottom_boundary(N_bbox.first) || v.at_top_boundary(N_bbox.second))
    {
        v.boundary=true;
        boundary.push_back(v.i);
    }
    else inner.push_back(v.i);
}


void grid::to_inner_or_to_boundary(subdomain & s, const vertex & v)
{
    if (v.at_top_boundary(s.N_bbox.second) || v.at_bottom_boundary(s.N_bbox.first))
    {
        bool boundary_neight = false;
        bool another_subd = false;

        for(auto n : v.neighbors)
        if(vertexes[n.first].boundary) boundary_neight=true;
        else if (vertexes[n.first].containing_subdomain!=s.i) another_subd=true;

        if(another_subd) s.boundary.push_back(v.i);
        else if(boundary_neight) s.inner.push_back(v.i);
        else throw simple_error("UNDEFINED VERTEX BOUNDARY TYPE");
    }
    else s.inner.push_back(v.i);
}


void grid::find_neighboring_subdomains()
{
    for(int i=0;i<subdomains.size();i++)
    {
        for(auto v_i : subdomains[i].boundary)
        for(auto j : vertexes[v_i].neighbors)
        subdomains[i].neighboring_subdomains.insert(vertexes[j.first].containing_subdomain);
        
        subdomains[i].neighboring_subdomains.erase(subdomains[i].i);
    }
}

int grid::find_index_in_part(int v_index, const vector<int> & input) const
{
    int pos = 0;
    auto v = std::find
    (
        input.begin(),
        input.end(),
        v_index
    );
        
    if(v!=input.end()) pos = v-input.begin();
    else pos = -1;

    return pos;
}

point grid::to_real_coord(voxel_point vp) const
{
    point p(dim);
    for(int i=0; i<dim; i++)
    p[i] = vp[i]*h[i]+bbox.first[i];

    return p;
}

voxel_point grid::to_voxel_coord(point p) const
{
    voxel_point vp(dim);
    for(int i=0; i<dim; i++)
    vp[i] = (p[i]-bbox.first[i])/h[i];

    return vp;
}

vector<interaction> grid::get_subdomains_interaction(int i, int j) const
{
    if( !(i<subdomains.size() && j<subdomains.size()) ) throw simple_error("index of subdomain out of range");
    if(i==j) return get_subdomain_interaction(i);

    interaction io;
    io.h = h;

    const vector<int> & output = subdomains[i].boundary;
    const vector<int> & input = subdomains[j].boundary;
    fill_interaction_between(io, input, output);

    return {io};
}

vector<interaction> grid::get_subdomain_interaction(int i) const
{
    if( !(i<subdomains.size()) ) throw simple_error("index of subdomain out of range");

    const vector<int> & inner_ = subdomains[i].inner;
    const vector<int> & boundary_ = subdomains[i].boundary;
    vector<interaction> io (3);

    io[0].h = h;
    io[1].h = h;
    io[2].h = h;

    fill_interaction_between(io[0], inner_, inner_);
    fill_interaction_between(io[1],inner_,boundary_);
    fill_interaction_between(io[2],boundary_,boundary_);

    return io;
}

void grid::fill_interaction_between(interaction & io, const vector<int> & input, const vector<int> & output) const
{
    io.io_size = {input.size(),output.size()};
    io.outputs = vector<vertex> (output.size());
    io.inputs = vector<list<pair<int, vertex>>> (output.size());

    #pragma omp parallel shared(io, input, output)
    #pragma omp for
    for(int i=0; i<output.size(); i++)
    {
        io.outputs[i] = vertexes[output[i]];

        find_neighbors_in_part(io.outputs[i], io.inputs[i], input);

        if(&output==&input)
        io.inputs[i].push_back({i, vertexes[input[i]]});
    }
}

void grid::find_neighbors_in_part(const vertex & v, list<pair<int, vertex>> & io_input, const vector<int> & input) const
{
    const vector<pair<int, int>> & nbs = v.neighbors;
    for(int j=0;j<nbs.size(); j++)
    {
        int k = find_index_in_part(nbs[j].first, input);
        #pragma omp critical
        if( k != -1) io_input.push_back({k, vertexes[input[k]]});
    }
}

interaction grid::get_inner_boundary_transpose_interaction() const
{
    interaction io;
    io.h = h;
    fill_interaction_between(io, inner, boundary);
    return io;
}

vector<point> grid::get_boundary_points() const
{
    vector<point> bp(boundary.size());
    for(int i=0;i<boundary.size();i++)
    bp[i] = to_real_coord(vertexes[boundary[i]].voxel_coord);

    return bp;
}

vector<point> grid::get_inner_points() const
{
    vector<point> bp(inner.size());
    for(int i=0;i<inner.size();i++)
    bp[i] = to_real_coord(vertexes[inner[i]].voxel_coord);

    return bp;
}

vector<int> grid::build_permutation_along_subd()
{
    double t1 = omp_get_wtime();
    vector<int> permutation (inner.size());
    vector<int> paddings (subdomains.size()+1);

    paddings[0] = 0;
    for(int i=0; i<subdomains.size(); i++)
    paddings[i+1] = paddings[i] + subdomains[i].boundary.size() + subdomains[i].inner.size();

    omp_set_nested(1);
    #pragma omp parallel for shared(paddings, permutation, subdomains, inner)
    for(int i=0; i<subdomains.size(); i++)
    {
        const subdomain & s = subdomains[i];

        #pragma omp parallel for shared(paddings, permutation, subdomains, inner, i)
        for(int j=0; j<s.inner.size(); j++)
        permutation[paddings[i]+j] = std::lower_bound
        (
            inner.begin(),
            inner.end(),
            s.inner[j]
        ) - inner.begin();

        #pragma omp parallel for shared(paddings, permutation, subdomains, inner, i)
        for(int j=0; j<s.boundary.size(); j++)
        permutation[paddings[i] + s.inner.size() + j] = std::lower_bound
        (
            inner.begin(),
            inner.end(),
            s.boundary[j]
        ) - inner.begin();
    }
    omp_set_nested(0);
    double t2 = omp_get_wtime();

    cout<<"Permutation`s indices was created with time: "<<t2-t1<<endl<<endl;
    return permutation;
}

void grid::assert_all_vertexes_has_subdomain() const
{
    #pragma omp parallel for shared(vertexes, subdomains)
    for(int i=0; i<vertexes.size(); i++)
    if(vertexes[i].boundary==false && vertexes[i].containing_subdomain==-1) throw("Any vertex must has subdomain or be boundary vertex");
}

void grid::setup_subdomains()
{
    find_neighboring_subdomains();
    for(int i=0;i<vertexes.size(); i++)
    if(!vertexes[i].boundary)
    to_inner_or_to_boundary(subdomains[vertexes[i].containing_subdomain], vertexes[i]);
}