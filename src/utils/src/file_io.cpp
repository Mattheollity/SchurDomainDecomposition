#include "file_io.h"
#include "problem_conditions.h"

grid_conds load_grid_conds(string grid_file)
{
    if(!fs::exists(grid_file)) throw simple_error(grid_file+" file is not exist");

    std::ifstream file(grid_file);

    bool change = true;
    string block;
    
    int ibuf;
    double dbuf;

    int dim;
    pair<point, point> bbox;
    pair<voxel_point, voxel_point> v_bbox;

    int subd_size;
    vector<pair<point, point>> subdomains;

    file>>block;
    if(block=="dim")
    {   if(file>>ibuf) dim = ibuf;
        else throw simple_error("can`t set dim");
    }
    else throw simple_error("expected dim");

    file>>block;
    if(block=="bbox")
    {   
        bbox.first = point(dim);
        bbox.second = point(dim);

        for(int i=0;i<dim;i++)
        file>>bbox.first[i];

        for(int i=0;i<dim;i++)
        file>>bbox.second[i];
    }
    else throw simple_error("expected bbox");

    file>>block;
    if(block=="vbbox")
    {
        v_bbox.first = voxel_point(dim);
        v_bbox.second = voxel_point(dim);
        for(int i=0;i<dim;i++)
        v_bbox.first[i]=0;
        
        for(int i=0;i<dim;i++)
        file>>v_bbox.second[i];
    }
    else throw simple_error("expected vbbox");

    file>>block;
    if(block=="amountsubd") file>>subd_size;
    else throw simple_error("expected amountsubd");

    subdomains.resize(subd_size);

    file>>block;
    if(block=="sections")
    for(int i=0; i<subd_size; i++)
    {
        file>>block;
        if(block!="subd") throw simple_error("expected subd");
        
        pair<point, point> bbox_sbd;
        bbox_sbd.first = point(dim);
        bbox_sbd.second = point(dim);

        for(int i=0;i<dim;i++)
        file>>bbox_sbd.first[i];

        for(int i=0;i<dim;i++)
        file>>bbox_sbd.second[i];
        subdomains[i]=bbox_sbd;
    }   
    else throw simple_error("expected sections");
    
    return grid_conds(bbox, v_bbox.second, subdomains);
}

time_integr_conds load_tic(string tic_file)
{
    std::ifstream file(tic_file);
    double start_time;
    int N;
    double T;
    double theta;
    double write_each;
    file>>start_time>>N>>T>>theta>>write_each;
    
    time_integr_conds tic(start_time, T, N, theta);
    tic.write_each = write_each;
    return tic;
}

void write_matrix(const spmat & m, string s)
{
    std::ofstream file(s+".txt");
    file << std::setprecision(std::numeric_limits<double>::max_digits10);
    for(int i=0;i<m.rows();i++)
    {
        for(int j=0;j<m.cols();j++)
        if(m.coeff(i, j)!=0) file<<m.coeff(i, j)<<" ";
        else file<<" . ";
        
        file<<"\n";
    }
}

void write_matrix_struct(const spmat & m, string s)
{
    std::ofstream file(s+".txt");
    file << std::setprecision(std::numeric_limits<double>::max_digits10);
    for(int i=0;i<m.rows();i++)
    {
        for(int j=0;j<m.cols();j++)
        if(m.coeff(i, j)!=0) file<<"!"<<" ";
        else file<<" . ";
        
        file<<"\n";
    }
}

void write_matrix(mat m, string s)
{
    std::ofstream file(s+".txt");
    file << std::setprecision(std::numeric_limits<double>::max_digits10);

    for(int i=0;i<m.rows();i++)
    {for(int j=0;j<m.cols();j++) file<<m.coeff(i, j)<<" "; file<<"\n";}
}

void write_matrix_struct(mat m, string s)
{
    std::ofstream file(s+".txt");
    file << std::setprecision(std::numeric_limits<double>::max_digits10);

    for(int i=0;i<m.rows();i++)
    {
        for(int j=0;j<m.cols();j++)
        {
            if(m.coeffRef(i,j)!=0) file<<"!"<<" ";
            else file<<"."<<" ";
        }
        file<<"\n";
    }
}
void write_matrix(vect m, string s)
{
    std::ofstream file(s+".txt");
    file << std::setprecision(std::numeric_limits<double>::max_digits10);

    for(int i=0;i<m.rows();i++)
    {for(int j=0;j<m.cols();j++) file<<m.coeff(i, j)<<" "; file<<"\n";}
}