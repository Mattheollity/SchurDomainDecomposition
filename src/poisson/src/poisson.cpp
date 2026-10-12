#include "poisson.h"

poisson_system::poisson_system(const grid & g, const approxs_conds & appc) :
    poissonlog(std::stringstream())
{
    double t0, t1, t2;
    
    p_system.dirichlet_conds = appc.dirichlet;
    p_system.approx = appc.approx;

    t1 = omp_get_wtime();
    t0 = t1;
    p_system.boundary_points = g.get_boundary_points();
    t2 = omp_get_wtime();
    cout<<"Get boundary points with time: "<<t2-t1<<endl;

    t1=t2;
    build_permutation_matrix(g);
    t2=omp_get_wtime();
    cout<<"Generate permutation matrix with time: "<<t2-t1<<endl;

    t1=t2;
    build_system(g);
    t2 = omp_get_wtime();
    cout<<"System was built with sizes: "<<endl;
    for(auto i : p_system.sizes) cout<<i<<" ";
    cout<<endl<<"with time: "<<t2-t1<<endl;
    cout<<"POISSON SYSTEM WAS SUCCESSFULLY GENERATED WITH TIME: "<<omp_get_wtime()-t0<<endl<<endl;
}

void poisson_system::build_permutation_matrix(const grid & g)
{

    const vector<int> & perm_vector = g.get_permutation_along_subd();
    vector<triplet> triplets (perm_vector.size());
    p_system.permutation = spmat(perm_vector.size(), perm_vector.size());

    #pragma omp parallel for shared(perm_vector, triplets)
    for(int i=0; i<triplets.size(); i++)
    triplets[i] = triplet(i, perm_vector[i], 1);

    p_system.permutation.setFromTriplets(triplets.begin(), triplets.end());    
}

const sparse_system & poisson_system::build_system(const grid & g)
{
    p_system.system.clear();
    p_system.sizes.clear();

    int N = g.get_subdomain().size();

    p_system.amount_of_parts = 2*N;

    
    for(int i=0;i<N;i++)
    {
        p_system.system.push_back(vector<spmat> (2*N, spmat(0,0)));
        p_system.system.push_back(vector<spmat> (2*N, spmat(0,0)));
    }

    double t1 = omp_get_wtime();

    #pragma omp parallel shared(N, p_system, g)
    {
    #pragma omp for collapse(2)
    for(int i=0;i<N; i++)
    for(int j=0;j<N; j++)
    {
    double t1 = omp_get_wtime();
    if(i==j) build_main_diagonal_system(g.get_subdomain_interaction(i), i);
    else build_interaction_system(g.get_subdomains_interaction(i, j), i, j);
    double t2 = omp_get_wtime();
    cout<<"build system for "<<i<<" "<<j<<" with time: "<<t2-t1<<endl;
    }
    }

    for(int i=0; i<p_system.amount_of_parts; i++)
    p_system.sizes.push_back(p_system.system[i][i].cols());

    double t2 = omp_get_wtime();
    cout<<"All systems were built with time: "<<t2-t1<<endl<<endl;
    t1 = t2;

    add_zeros_matrixes();

    build_right_addition(g);
    t2 = omp_get_wtime();
    cout<<"Right addition was built with time: "<<t2-t1<<endl;

    return p_system;
}

void poisson_system::build_main_diagonal_system(const vector<interaction> & io, int i)
{
    if(io.size()!=3) throw simple_error("Uncorrect interaction object");

    fill_spmat(p_system.system[2*i][2*i], io[0]);
    fill_spmat(p_system.system[2*i+1][2*i], io[1]);
    fill_spmat(p_system.system[2*i+1][2*i+1], io[2]);
    p_system.system[2*i][2*i+1] = p_system.system[2*i+1][2*i].transpose();
}

void poisson_system::fill_spmat(spmat & spm, const interaction & io)
{
    spm = spmat
    (
        io.io_size.second,
        io.io_size.first
    );
    list<triplet> triplets;

    for(int i=0; i<io.outputs.size(); i++)
    {
    for(auto input : io.inputs[i])
    {
        triplets.push_back(
            triplet(
            i,
            input.first,
            p_system.approx(io.h, io.outputs[i], input.second)
            )
        );
    }}
        spm.setFromTriplets(triplets.begin(), triplets.end());
}

void poisson_system::build_interaction_system(const vector<interaction> & io, int i, int j)
{
    if(io.size()!=1) throw simple_error("uncorrect interaction object");
    fill_spmat(p_system.system[2*i+1][2*j+1], io[0]);
    p_system.system[2*j+1][2*i+1] = p_system.system[2*i+1][2*j+1].transpose();
}

void poisson_system::write_to_file(string s)
{
    std::ofstream file(s);

    file<<p_system.system.size()<<endl;
    file<<p_system.system[0].size()<<endl;

    for(auto v : p_system.system)
    for(auto sp : v)
    write_spmat(file, sp);
            
    file.close();
}

void poisson_system::write_spmat(std::ofstream & file, spmat & sp)
{
    file<<sp.rows();
    file<<sp.cols();
    file<<sp.nonZeros()<<endl;

    for (int k=0; k<sp.outerSize(); ++k)
    for (spmat::InnerIterator it(sp,k); it; ++it)
    {
        file<<it.row()<<endl;
        file<<it.col()<<endl;
        file<<it.value()<<endl;
    }
}

void poisson_system::load_from_file(string s)
{
    std::ifstream file(s);
    pair<int,int> sizes;
    file>>sizes.first;
    file>>sizes.second;

    for(int i=0; i<sizes.first; i++)
    p_system.system.push_back(vector<spmat> (sizes.second));

    for(int i=0;i<sizes.first; i++)
    for(int j=0;j<sizes.second; j++)
    load_spmat(file, p_system.system[i][j]);

    file.close();
}

void poisson_system::load_spmat(std::ifstream & file, spmat & sp)
{
    int rows,cols,nonz;
    int r,c,v;

    file>>rows>>cols>>nonz;

    sp = spmat(rows, cols);
    list<triplet> trp;

    for (int k=0;k<nonz;k++)
    {
        file>>r>>c>>v;
        trp.push_back(triplet(r,c,v));
    }

    sp.setFromTriplets(trp.begin(), trp.end());
}

void poisson_system::build_right_addition(const grid & g)
{
    double t1 = omp_get_wtime();
    interaction io = g.get_inner_boundary_transpose_interaction();

    double t2 = omp_get_wtime();
    cout<<"IO for right part object was built with time: "<<t2-t1<<endl;

    t1 = t2;

    fill_spmat(p_system.right_addition, io);
    p_system.right_addition 
    =
    p_system.permutation
    *
    p_system.right_addition.transpose();

    t2 = omp_get_wtime();
    cout<<"Right matrix was filled with time: "<<t2-t1<<endl;
}

vect sparse_system::boundary_at_time(double t) const
{
    vect bat(boundary_points.size());
    bat.setZero();
    for(int i=0; i<boundary_points.size(); i++)
    bat.coeffRef(i) = dirichlet_conds(boundary_points[i], t);

    return bat;
}

void poisson_system::add_zeros_matrixes()
{
    int N = p_system.amount_of_parts;
    #pragma omp parallel shared(N)
    {
    #pragma omp for collapse(2)
        for(int i=0; i<N; i++)
        for(int j=0; j<N; j++)
        if(p_system.system[i][j].size()==0) p_system.system[i][j] = spmat(p_system.sizes[i], p_system.sizes[j]);
    }
}