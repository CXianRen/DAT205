#ifndef __M_SOLVER_H__
#define __M_SOLVER_H__

#include <Eigen/Core>
#include <Eigen/Sparse>

class Solver
{
public:
    Solver() = default;
    ~Solver() = default;
    virtual void setTolerance(float tol) {};

    virtual void setMaxIterations(int max_iter) {};

    virtual void getIterations(int &iter) {};

    virtual void getError(float &error) {};

    virtual void solve(
        Eigen::VectorXf &x,
        Eigen::VectorXf &b) = 0;

    virtual void solveWithGuess(
        Eigen::VectorXf &x,
        Eigen::VectorXf &b) = 0;
};

class EigenSolver : public Solver
{
public:
    EigenSolver()
    {
        // set tolerance
        ICCG.setTolerance(1e-4);
    }

    ~EigenSolver() = default;

    void setTolerance(float tol) override
    {
        ICCG.setTolerance(tol);
    }

    void setMaxIterations(int max_iter) override
    {
        ICCG.setMaxIterations(max_iter);
    }

    void getIterations(int &iter) override
    {
        iter = ICCG.iterations();
    }

    void getError(float &error) override
    {
        error = ICCG.error();
    }

    void compute(Eigen::SparseMatrix<float, Eigen::RowMajor> &A)
    {
        ICCG.compute(A);
    }

    void solve(
        Eigen::VectorXf &x,
        Eigen::VectorXf &b) override;

    void solveWithGuess(
        Eigen::VectorXf &x,
        Eigen::VectorXf &b) override;

private:
    Eigen::ConjugateGradient<
        Eigen::SparseMatrix<float>,
        Eigen::Lower | Eigen::Upper,
        // Eigen::IncompleteCholesky<float>>
        // Eigen::IncompleteLUT<float>>
        Eigen::DiagonalPreconditioner<float>>
        ICCG;
};

#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cusparse.h>

class CudaSolver : public Solver
{
public:
    CudaSolver(){};
    ~CudaSolver();

    void compute(Eigen::SparseMatrix<float, Eigen::RowMajor> &A);

    void getIterations(int &iter) override;

    void getError(float &error) override;

    void solve(
        Eigen::VectorXf &x,
        Eigen::VectorXf &b) override;

    void solveWithGuess(
        Eigen::VectorXf &x,
        Eigen::VectorXf &b) override {};

private:
    const int max_iter = 1000;
    const float tol = 1e-5f;

    int N = 0;  // size of the matrix, rows
    int nz = 0; // number of non-zero elements
    // data on the host
    int *I = NULL;     // row indices
    int *J = NULL;     // column indices
    float *val = NULL; // values

    float *x, *rhs; // solution and right hand side

    //
    float a, b, na, r0, r1;

    // data on the device
    int *d_col = NULL;
    int *d_row = NULL;
    float *d_val = NULL;
    float *d_x = NULL, *d_r = NULL;
    float *d_p = NULL, *d_Ax = NULL;

    float dot;

    int k = 0; // iteration count
    float alpha, beta, alpham1;

    // handle to the CUBLAS context
    cublasHandle_t cublasHandle = 0;
    cublasStatus_t cublasStatus;

    // handle to the CUSPARSE context
    cusparseHandle_t cusparseHandle = 0;

    // Wrap raw data into cuSPARSE generic API objects
    cusparseSpMatDescr_t matA = NULL;
    cusparseDnVecDescr_t vecx = NULL;
    cusparseDnVecDescr_t vecp = NULL;
    cusparseDnVecDescr_t vecAx = NULL;
};

// void cudaSolve(
//     Eigen::SparseMatrix<float> &A,
//     Eigen::VectorXf &b,
//     Eigen::VectorXf &x);

#endif // __M_CUSOLVER_H__