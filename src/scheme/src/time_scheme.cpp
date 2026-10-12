#include "time_scheme.h"
#include "Eigen/src/SparseCholesky/SimplicialCholesky.h"
#include "Eigen/src/SparseLU/SparseLU.h"
#include "Eigen/src/SparseQR/SparseQR.h"
#include "common.h"
#include "file_io.h"
#include "problem_conditions.h"
#include <fstream>
#include <iomanip>
#include <omp.h>
#include <ostream>
#include <string>
#include <vector>

void tsws::expimp_matrix_manipulation(spsys & target, const spsys & origin, double coeff)
{
    #pragma omp parallel for shared(target, origin, coeff)
    for(int i=0; i<origin.size(); i++)
    {
        for(int j=0; j<origin.size(); j++)
        target[i][j]*=coeff;

        int inner_size = origin[i][i].cols();
        for(int k=0; k<inner_size; k++)
        target[i][i].coeffRef(k, k)+=1;
    }
}

vector<vect> tsws::split_vect(const vect & origin)
{
    vector<int> csum = {0};
    for(int i=0;i<origin_system.sizes.size();i++)
    csum.push_back(csum.back()+origin_system.sizes[i]);

    vector<vect> target;
    for(int i=0; i<origin_system.sizes.size(); i++)
    {
        target.push_back(vect(origin_system.sizes[i]));
        target.back() = origin.block(csum[i], 0, csum[i+1]-csum[i], 1);
    }
    return target;
}

void tsws::concate_vect(vect & target, const vector<vect> & origin)
{
    vector<int> csum = {0};
    for(int i=0;i<origin_system.sizes.size();i++)
    csum.push_back(csum.back()+origin_system.sizes[i]);

    target = vect(csum.back());

    for(int i=0; i<origin_system.sizes.size(); i++)
    target.block(csum[i], 0, csum[i+1]-csum[i], 1) = origin[i];
}

const one_step_conds & tsws::set_osc()
{

    vect y0 = origin_system.permutation*tic.y0_vect;

    osc = one_step_conds();
    osc.y0 = split_vect(y0);
    osc.y = osc.y0;

    osc.left_system = origin_system.system;
    osc.right_system = origin_system.system;
    osc.parts_size = origin_system.sizes;

    expimp_matrix_manipulation(
        osc.left_system,
        origin_system.system,
        (1.0)*tic.T*tic.theta
    );
    expimp_matrix_manipulation(
        osc.right_system, 
        origin_system.system,
        (-1.0)*tic.T*(1.0-tic.theta)
    );
    build_right_part_osc();
    return osc;
}

const schur_system & tsws::set_ss()
{
    double t1 = omp_get_wtime();
    ss = schur_system();
    int size = origin_system.amount_of_parts/2;
    int schur_size = 0;

    vector<int> block_sizes(size);

    for(int i=0; i<size; i++)
    block_sizes[i] = origin_system.sizes[2*i+1];


    for(int bs : block_sizes)
    schur_size += bs;

    ss.block_sizes = block_sizes;
    ss.schur_mat = spmat(schur_size, schur_size);
    ss.right = vect::Zero(schur_size);
    ss.y = vect::Zero(schur_size);

    double t2 = omp_get_wtime();
    cout<<"set schur system parts with time: "<<t2-t1<<endl;
    t1=t2;

    build_schur_matrix();

    t2 = omp_get_wtime();
    cout<<"built schur system matrix with time: "<<t2-t1<<endl;

    return ss;
}


const one_step_conds & tsws::solve_osc()
{
    #pragma omp parallel shared(osc, ss)
    {
    #pragma omp for
    for(int i=0; i<ss.block_sizes.size(); i++)
    {
        osc.right[2*i] -= osc.left_system[2*i][2*i+1]*osc.y[2*i+1];
        Eigen::SparseLU<spmat, Eigen::COLAMDOrdering<int>> solver;
        solver.analyzePattern(osc.left_system[2*i][2*i]);
        solver.compute(osc.left_system[2*i][2*i]);
    
        vect tmpvec = solver.solve(osc.right[2*i]); 
        osc.y[2*i] = solver.solve(osc.right[2*i]);
    }
    }

    return osc;
}

const schur_system & tsws::solve_ss()
{
    double t1 = omp_get_wtime();
    cout<<"Start schur solution"<<endl;

    if(ss.schur_mat.rows()==0) return ss;
    Eigen::SimplicialLDLT<spmat> solver;
    solver.analyzePattern(ss.schur_mat);
    solver.compute(ss.schur_mat);
    ss.y = solver.solve(ss.right);

    double t2 = omp_get_wtime();

    cout<<"Schur matrix was solved with time: "<<t2-t1<<endl;
    cout<<"Start schur writting"<<endl;

    int csum = 0;
    for(int i=0;i<ss.block_sizes.size();i++)
    {
        
        osc.y[2*i+1] = ss.y.block(
            csum,
            0,
            ss.block_sizes[i],
            1
        );
        csum+=ss.block_sizes[i];
    }
    return ss;
}

const one_step_conds & tsws::update_osc()
{
    vect sol_next;
    concate_vect(sol_next, osc.y);
    solution.push_back(origin_system.permutation.transpose()*sol_next);
    int i = solution.size()-1;
    cout<<"solution at time step "<<i<< " was get"<<endl;

    if(tic.write_each != 1 && i%tic.write_each==tic.write_each-1)
    {
        std::ofstream file(to_string(i));
        for(auto d : solution[i]) file<<std::setprecision(16)<<d<<endl;
        cout<<"solution at time step "<<i<< " was wrotten"<<endl;
    }
    if(i%tic.write_each==0)
    {
        std::ofstream file(to_string(i));
        for(auto d : solution[i]) file<<std::setprecision(16)<<d<<endl;
        cout<<"solution at time step "<<i<< " was wrotten"<<endl;
    }
    
    
    osc.y0 = osc.y;
    build_right_part_osc();

    return osc;
}

const schur_system & tsws::update_ss()
{
    vector<int> csum = {0};
    for(auto i : ss.block_sizes) csum.push_back(csum.back()+i);

    #pragma omp parallel shared(ss, csum)
    {
    #pragma omp for
    for(int i=0; i<ss.block_sizes.size(); i++)
    {
        mat tmpmat;
        Eigen::SimplicialLDLT<spmat> solver;
        solver.analyzePattern(osc.left_system[2*i][2*i]);
        solver.factorize(osc.left_system[2*i][2*i]);
        tmpmat = osc.right[2*i+1]
            -
            osc.left_system[2*i+1][2*i]
            *
            solver.solve(osc.right[2*i]);
        #pragma omp critical
        {
        ss.right.block(
            csum[i],
            0,
            csum[i+1]-csum[i],
            1
        ) = tmpmat;
        }
    }
    }
    return ss;
}

const schur_system & tsws::build_schur_matrix()
{
    vector<int> csum = {0};
    for(auto i : ss.block_sizes) csum.push_back(csum.back()+i);

    double t1 = omp_get_wtime();

    int trp_size = 0;
    for(int i=0; i<ss.block_sizes.size(); i++)
    for(int j=i; j<ss.block_sizes.size(); j++)
    if(i==j) trp_size+=ss.block_sizes[i]*ss.block_sizes[i];
    else trp_size+=2*osc.left_system[2*i+1][2*j+1].nonZeros();

    list<triplet> triplets;

    #pragma omp parallel shared(ss, csum, triplets) firstprivate(t1)
    {
    #pragma omp for collapse(2)
    for(int i=0; i<ss.block_sizes.size(); i++)
    for(int j=0; j<ss.block_sizes.size(); j++)
        if(i==j && ss.block_sizes[i]!=0)
        {
            mat tmpmat;
            vect tmpvec;
            Eigen::SparseLU<spmat> solver;
            solver.analyzePattern(osc.left_system[2*i][2*i]);
            solver.factorize(osc.left_system[2*i][2*i]);

            tmpmat = osc.left_system[2*i+1][2*i+1]
                -
                osc.left_system[2*i+1][2*i]
                *
                solver.solve(osc.left_system[2*i][2*i+1]);


            tmpvec = osc.right[2*i+1]
                -
                osc.left_system[2*i+1][2*i]
                *
                solver.solve(osc.right[2*i]);

            #pragma omp critical
            {
                for(int k=0;k<tmpmat.rows();k++)
                for(int q=0;q<tmpmat.cols(); q++)
                triplets.push_back(triplet(csum[i]+k, csum[j]+q, tmpmat.coeff(k,q)));

                ss.right.block(
                    csum[i],
                    0,
                    csum[i+1]-csum[i],
                    1
                ) = tmpvec;
            }       

            double t2=omp_get_wtime();
            cout<<"build ss matrix main at "<<i<<" "<<j<<" with time: "<<t2-t1<<endl;
        }
        else
        {
            #pragma omp critical
            {
                for (int k = 0; k < osc.left_system[2*i+1][2*j+1].outerSize(); ++k)
                for (Eigen::SparseMatrix<double>::InnerIterator it(osc.left_system[2*i+1][2*j+1], k); it; ++it)
                    triplets.push_back(triplet(csum[i]+it.row(), csum[j]+it.col(), it.value()));
            }
            double t2=omp_get_wtime();
            cout<<"copy block to ss matrix at "<<i<<" "<<j<<" with time: "<<t2-t1<<endl;
        }
    }

    ss.schur_mat.setFromTriplets(triplets.begin(), triplets.end());

    return ss;
}

tsws::tsws(const sparse_system & sparse_system_, const time_integr_conds & tic_) :
    origin_system(sparse_system_),
    tic(tic_)
{
    
}

const vector<vect> & tsws::solve()
{
    solution.clear();
    solution.push_back(tic.y0_vect);

    int i = solution.size()-1;
    std::ofstream file(to_string(i));
    for(auto d : solution[i]) file<<std::setprecision(16)<<d<<endl;

    cout<<"solution at time step "<<i<<" was wrotten"<<endl<<endl;

    double t = omp_get_wtime();
    set_osc();
    double t2 = omp_get_wtime();
    cout<<"set osc with time: "<<t2-t<<endl;
    t = t2;
    t2 = omp_get_wtime();
    set_ss();
    cout<<"set osc with time: "<<t2-t<<endl<<endl;
    t = t2;
    t2 = omp_get_wtime();

    for(int i=1;i<tic.N; i++)
    {
        double time_a = omp_get_wtime();
        cout<<"time step: "<<i<<endl;
        solve_ss();
        double time_b = omp_get_wtime();
        cout<<"schur system was solved with time: "<<time_b-time_a<<endl;
        solve_osc();
        double time_c = omp_get_wtime();
        cout<<"poisson system was solved with time: "<<time_c-time_b<<endl;
        update_osc();
        double time_d = omp_get_wtime();
        cout<<"poisson system was updated with time: "<<time_d-time_c<<endl;
        update_ss();
        double time_e = omp_get_wtime();
        cout<<"schur system was updated with time: "<<time_e-time_d<<endl;
        cout<<"time step "<<i<<" was wrotten with time: "<<time_e - time_a<<endl<<endl;
    }

    return solution;
}

double tsws::get_current_time() const
{
    return tic.time.first+solution.size()*tic.T;
}

void tsws::build_right_part_osc()
{
    osc.right = osc.y0;

    for(int i=0;i<osc.parts_size.size();i++)
    osc.right[i].setZero();

    #pragma omp parallel for shared(osc)
    for(int i=0;i<osc.parts_size.size();i++)
    for(int j=0; j<osc.parts_size.size(); j++)
    osc.right[i] += osc.right_system[i][j]*osc.y0[j];

    // vect yb = (-1)*(1-tic.theta)*tic.T*origin_system.right_addition*origin_system.boundary_at_time(get_current_time()-tic.T);
    // vect yb2 = (-1)*(tic.theta)*tic.T*origin_system.right_addition*origin_system.boundary_at_time(get_current_time());

    // int csum = 0;
    // for(int i=0; i<osc.parts_size.size(); i++)
    // {
    //     osc.right[i]+=yb.block(csum, 0, osc.parts_size[i], 1);
    //     osc.right[i]+=yb2.block(csum, 0, osc.parts_size[i], 1);
    //     csum+=osc.parts_size[i];
    // }
}

const vect tsws::get_solution_at_time(int step) const
{
    std::ifstream file(to_string(step));
    vect v(solution.back().rows());
    for(int i=0;i<v.rows();i++) file>>v.coeffRef(i);
    
    return v;
}