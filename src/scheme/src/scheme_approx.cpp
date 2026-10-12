#include "scheme_approx.h"
#include "common.h"
#include "file_io.h"

vect get_res_at_time_step(const grid & g, const time_integr_conds & tic, const sparse_system & p_system, int N)
{
    auto points = g.get_inner_points();
    double time = tic.time.first+tic.T*N;

    vect td = get_time_deriv_at_time_step(g, tic, N);
    vect yn = get_solution_at_time_step(g, tic, N);

    vect res(td.rows());

    td = p_system.permutation*td;
    yn = p_system.permutation*yn;


    vect td_from_laplace(yn.rows());
    td_from_laplace.setZero();

    vector<int> csum {0};
    for(auto i : p_system.sizes) csum.push_back(csum.back()+i);

    for(int i=0;i<p_system.sizes.size();i++)
    for(int j=0;j<p_system.sizes.size();j++)
    {
        td_from_laplace.block(csum[i], 0, csum[i+1]-csum[i], 1) += p_system.system[i][j]*yn.block(csum[j], 0, csum[j+1]-csum[j], 1);
        
    }
    for(int i=0;i<td.rows();i++)
    {
        res.coeffRef(i)=std::abs(td_from_laplace.coeff(i)-td.coeff(i));
        if(res.coeffRef(i)>1) cout<<"Get anomal value: "<<res.coeffRef(i)<<" at index: "<<i<<endl;
    }

    return p_system.permutation.transpose()*res;
}


vector<vector<pair<double, double>>> check_laplace_approx(string filename)
{

    std::ifstream file(filename);
    int amount_of_grids;
    file>>amount_of_grids;

    vector<string> gc_names(amount_of_grids);
    vector<grid> grids;

    for(int i=0;i<amount_of_grids;i++)
    {
        file>>gc_names[i];
        grids.push_back(load_grid_conds(gc_names[i]));
    }

    time_integr_conds tic = load_tic("tic");

    int amount_of_times;
    file>>amount_of_times;
    vector<int> time_steps(amount_of_times);
    
    for(int i=0;i<amount_of_times;i++)
    file>>time_steps[i];


    vector<vector<pair<double, double>>> res(amount_of_times);

    for(int i = 0; i<amount_of_times; i++)
    {res[i] = vector<pair<double, double>> (amount_of_grids);
    for(int j = 0; j<grids.size(); j++)
    {
        approxs_conds appc({approx, dirich_boundary_conds});
        poisson_system p_system(grids[j], appc);

        double max_res = 0;
        vect res_vect = get_res_at_time_step(grids[j], tic, p_system.get_system(), time_steps[i]);

        for(int k=0;k<res_vect.rows();k++)
        {
        if(res_vect.coeff(k)>max_res) max_res=res_vect.coeff(k);
        }

        res[i][j].first=log(*grids[j].get_h().begin());
        res[i][j].second=log(max_res);
    }}

    return res;
}

int main()
{
    auto res = check_laplace_approx("check_space_approx_conds");
    for(auto time_step : res)
    {
        for(auto grid_step : time_step)
        cout<<std::setprecision(16)<<grid_step.first<<", "<<grid_step.second<<endl;

        cout<<"----------"<<endl;
    }

    return 0;
}