#include "scheme_converg.h"
#include "common.h"
#include "file_io.h"
#include "approx.h"
#include "poisson.h"
#include "problem_conditions.h"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

double get_cheb_norm(const vect & v)
{
    int N = v.rows();

    double norm = 0;
    for(int i=0;i<N;i++)
    if(norm<std::abs(v.coeff(i))) norm=std::abs(v.coeff(i));

    return norm;
}

double get_eqlid_norm(const vect & v)
{
    int N = v.rows();

    double norm = 0;
    for(int i=0;i<N;i++)
    norm+=v.coeff(i)*v.coeff(i);

    return sqrt(norm/N);
}

vect load_vect(string filename, int N)
{
    if(!fs::exists(filename)) throw simple_error("CAN`T OPEN "+filename+". FILE is not exist");
    vect v(N);

    std::ifstream file(filename);
    for(int k=0;k<N;k++)
    file>>v.coeffRef(k);

    return v;
}

vect multiply_spsys_and_vect(const sparse_system & system, const vect & v, const vect & vn, const time_integr_conds & tic)
{
    vect v_ = system.permutation*v;
    vect vn_ = system.permutation*v;

    vect compt_right(v.rows());
    compt_right.setZero();

    vector<int> csum = {0};
    for(int i : system.sizes)
    csum.push_back(csum.back()+i);

    for(int a=0; a<system.amount_of_parts; a++)
    for(int b=0; b<system.amount_of_parts; b++)
    compt_right.block(csum[a], 0, csum[a+1]-csum[a], 1) 
    +=
    (1.0-tic.theta)*system.system[a][b]*v_.block(csum[b], 0, csum[b+1]-csum[b], 1)
    +
    tic.theta*system.system[a][b]*vn_.block(csum[b], 0, csum[b+1]-csum[b], 1);

    return system.permutation.transpose() * compt_right;
}

vect get_approx_solution_at_time_step(const grid & g, const time_integr_conds & tic, int N, const sparse_system & spsys)
{
    auto points = g.get_inner_points();
    double time = tic.time.first+tic.T*N;

    // double itern = 500;

    vect another = get_solution_at_time_step(g, tic, 0);
    // vect sol = another;
    // double timepow = 1.0;

    // for(int i=1;i<itern;i++)
    // {
    //     timepow*=(-time);
    //     sol+=(1.0/i)*timepow*multiply_spsys_and_vect(spsys, sol, sol, tic);
    // }
    vect sol = eexpokit::expv(-time, spsys.system[0][0], another);

    return sol;
}

vector<double> get_eqlid_converg_at_time_step(const grid & g, const time_integr_conds & tic, const sparse_system & p_system, int N, string dir_name)
{
    vect origin_sol = get_solution_at_time_step(g, tic, N);
    //vect origin_sol_next = get_solution_at_time_step(g, tic, N+1);

    vect deriv = get_time_deriv_at_time_step(g, tic, N);

    vect compt_sol = load_vect(dir_name+"/"+to_string(N), origin_sol.rows());
    vect compt_sol_next = load_vect(dir_name+"/"+to_string(N+1), origin_sol.rows());
    vect cderiv = multiply_spsys_and_vect(p_system, compt_sol, compt_sol_next, tic);

    vect sol_diffr = deriv+cderiv;
    double sol_diffr_norm = get_eqlid_norm(sol_diffr);
    double sol_norm = get_eqlid_norm(deriv);


    //double right_norm = get_eqlid_L_converg(g, tic, p_system, N, dir_name);

    return {sol_diffr_norm/sol_norm};
}

// double get_eqlid_L_converg(const grid & g, const time_integr_conds & tic, const sparse_system & p_system, int N, string dir_name)
// {
//     if(N==tic.N+1) throw simple_error("Uncorrect N value");

//     vect compt_sol = load_vect(dir_name+"/"+to_string(N), g.get_inner().size());
//     vect compt_sol_next = load_vect(dir_name+"/"+to_string(N+1), compt_sol.rows());

//     //vect compt_td = (1.0/tic.T)*(compt_sol_next-compt_sol);
//     vect compt_laplace = multiply_spsys_and_vect(p_system, compt_sol, compt_sol_next, tic);

//     vect right = compt_td-compt_laplace;
//     return get_eqlid_norm(right);
// }

vector<vector<double>> check_space_converg(string filename)
{
    auto check_data = get_check_data(filename);
    int amount_of_dirs = check_data.first.size();
    vector<string> & dir_names = check_data.first;
    vector<int> & time_steps = check_data.second;

    approxs_conds appc({approx, dirich_boundary_conds});

    vector<vector<double>> eqlid(amount_of_dirs);
    for(int i=0;i<amount_of_dirs;i++)
    {
        time_integr_conds tic = load_tic(dir_names[i]+"/tic");
        grid g(load_grid_conds(dir_names[i]+"/grid"));
        poisson_system p_system(g, appc);
        {
            int N = tic.N-tic.N%tic.write_each-tic.write_each;

            cout<<"Dir: "<<dir_names[i]<<endl;
            cout<<"time step: "<<N<<endl;

            vector<double> error = get_eqlid_converg_at_time_step(g, tic, p_system.get_system(), 0, dir_names[i]);
            for(int q=0;q<error.size(); q++)
            error[q]=log(error[q]);

            cout<<g.get_h()[0]<<endl;
            eqlid[i].push_back(log(g.get_h()[0]));
            cout<<tic.T<<endl;
            eqlid[i].push_back(log(tic.T));
            eqlid[i].insert(eqlid[i].end(), error.begin(), error.end());
            //eqlid[i][j].push_back(log(error/rnorm));
        }
    }

    return eqlid;
}

pair<vector<string>, vector<int>> get_check_data(string filename)
{
    std::ifstream file(filename);
    int amount_of_directories;
    file>>amount_of_directories;

    vector<string> dir_names(amount_of_directories);

    for(int i=0;i<amount_of_directories;i++)
    file>>dir_names[i];

    int amount_of_times;
    file>>amount_of_times;
    vector<int> time_steps(amount_of_times);
    for(int i=0;i<amount_of_times;i++)
    file>>time_steps[i];

    file.close();
    return {dir_names, time_steps};
}

int main()
{
    auto eqlid = check_space_converg("check_eqlid_converg");
    cout<<endl<<"Error by dir:"<<endl;
    for(auto dir : eqlid)
    {
        for(auto v : dir)
        cout<<std::setprecision(5)<<v<<" ";

        cout<<endl<<"----------"<<endl;
    }

    return 0;
}