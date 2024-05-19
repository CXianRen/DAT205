#ifndef __M_SOLVER_H__
#define __M_SOLVER_H__


#include <Eigen/Core>
#include <Eigen/Sparse>

void cudaSolve( 
    Eigen::SparseMatrix<float, Eigen::RowMajor> &A , 
    Eigen::VectorXf &b, 
    Eigen::VectorXf &x
);

#endif // __M_CUSOLVER_H__