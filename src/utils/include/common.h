#ifndef COMMON_
#define COMMON_

#include <exception>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <sstream>

#include <string>
#include <vector>
#include <list>
#include <deque>
#include <set>

#include <utility>
#include <cmath>

#include <omp.h>

#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Sparse"
#include "Eigen/src/SparseLU/SparseLU.h"
#include <Eigen/Core>

typedef Eigen::SparseMatrix<double> spmat;
typedef Eigen::MatrixXd mat;
typedef Eigen::VectorXd vect;
typedef Eigen::Triplet<double> triplet;

using std::cout;
using std::cin;
using std::endl;

using std::string;
using std::to_string;
using std::vector;
using std::list;
using std::pair;
using std::deque;
using std::set;

class simple_error : public std::exception
{
    public:
    simple_error(string message)
    {
        std::cerr<<message<<endl;
        std::exit(1);
    }
};

#endif