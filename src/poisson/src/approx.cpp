#include "common.h"
#include "file_io.h"
#include "poisson.h"
#include <cmath>
#include <numbers>
#include "eexpokit/expv.hpp"

double trip_sin(const point p)
{
    double c = 1.0;
    for(int d=0; d<p.size(); d++)
    c*=sin(M_PI*p[d]);

    return c;
}

double m_pi = M_PI;

vect y0(const grid & gr, vect & v)
{
    const vector<vertex> & vtx = gr.get_vertexes();
    auto inner = gr.get_inner();
    v = vect(inner.size());
    for(int i=0;i<inner.size(); i++)
    v.coeffRef(i)=trip_sin(
        gr.to_real_coord(vtx[inner[i]].voxel_coord)
    );
    return v;
}

double approx(const vector<double> h, const vertex & main, const vertex & neighbor)
{
    //Задача Дирихле, однородные гру
    //ut = vec(k)*vec(laplac)u

    if(h.size()==0) throw simple_error("h is empty!");
    int dim = h.size();
    vector<double> k (dim, 1);

    double coeff = 0;

    if(main.i==neighbor.i)
    {
        for(int d=0;d<dim;d++)
        coeff -= 2*k[d]/(h[d]*h[d]);
    }
    
    else
    {
        for(auto np : main.neighbors)
        if(np.first==neighbor.i)
        coeff=k[np.second]/(h[np.second]*h[np.second]);
    }

    return -coeff;

}

double dirich_boundary_conds(point, double t) {return 0;}

double real_solution(int dim, const point p, double t)
{
    return exp(-M_PI*M_PI*dim*t)*trip_sin(p);
}

vect get_solution_at_time_step(const grid & g, const time_integr_conds & tic, int N)
{
    auto points = g.get_inner_points();
    double time = tic.time.first+tic.T*N;

    vect sol(points.size());

    for(int i=0; i<points.size(); i++)
    sol.coeffRef(i) = real_solution(g.get_dim(), points[i], time);

    return sol;
}



vect get_time_deriv_at_time_step(const grid & g, const time_integr_conds & tic, int N)
{
    vect td = (-g.get_dim())*M_PI*M_PI*get_solution_at_time_step(g, tic, N);

    return td;
}

vect get_dummy_time_deriv_at_time_step(const grid & g, const time_integr_conds & tic, int N)
{
    vect a = get_solution_at_time_step(g, tic, N);
    vect b = get_solution_at_time_step(g, tic, N+1);
    vect td = (1.0/tic.T)*(b-a);

    return td;
}

vect get_laplace_at_time_step(const grid & g, const time_integr_conds & tic, const sparse_system & p_system, int N)
{
    auto points = g.get_inner_points();
    vect yn = get_solution_at_time_step(g, tic, N);
    vect ynn = yn;
    ynn.setZero();

    const mat P = p_system.permutation;
    yn = P*yn;

    vector<int> csum = {0};
    for(int i : p_system.sizes)
    csum.push_back(csum.back()+i);

    for(int i=0; i<p_system.amount_of_parts; i++)
    for(int j=0; j<p_system.amount_of_parts; j++)
    ynn.block(csum[i], 0, csum[i+1]-csum[i], 1) += p_system.system[i][j]*yn.block(csum[j], 0, csum[j+1]-csum[j], 1);

    return P.transpose()*ynn;
}

double get_laplace_eqlid_at_time_step(const grid & g, const time_integr_conds & tic, const sparse_system & p_system, int N)
{
    vect lpv = get_laplace_at_time_step(g, tic, p_system, N);
    double norm = 0;
    
    for(int i=0;i<lpv.rows(); i++)
    norm+=lpv.coeff(i)*lpv.coeff(i);

    norm/=lpv.rows();

    return sqrt(norm);

}