#ifndef __M_SOLVER_H__
#define __M_SOLVER_H__

#include <Eigen/Core>
#include <Eigen/Sparse>

class Solver
{
public:
    Solver() = default;
    ~Solver() = default;
    virtual void setTolerance(double tol) {};

    virtual void setMaxIterations(int max_iter) {};

    virtual void getIterations(int &iter) {};

    virtual void getError(double &error) {};

    virtual void solve(
        Eigen::VectorXd &x,
        Eigen::VectorXd &b) = 0;

    virtual void solveWithGuess(
        Eigen::VectorXd &x,
        Eigen::VectorXd &b) = 0;
};

class EigenSolver : public Solver
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    EigenSolver()
    {
        // set tolerance
        ICCG.setTolerance(1e-6);
    }

    ~EigenSolver() = default;

    void setTolerance(double tol) override
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

    void getError(double &error) override
    {
        error = ICCG.error();
    }

    void compute(Eigen::SparseMatrix<double, Eigen::RowMajor> &A)
    {
        ICCG.compute(A);
    }

    void solve(
        Eigen::VectorXd &x,
        Eigen::VectorXd &b) override;

    void solveWithGuess(
        Eigen::VectorXd &x,
        Eigen::VectorXd &b) override;

private:
    Eigen::ConjugateGradient<
        Eigen::SparseMatrix<double>,
        Eigen::Lower | Eigen::Upper,
        // Eigen::IncompleteCholesky<double>>
        // Eigen::IncompleteLUT<double>>
        Eigen::DiagonalPreconditioner<double>>
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

    void compute(Eigen::SparseMatrix<double, Eigen::RowMajor> &A);

    void getIterations(int &iter) override;

    void getError(double &error) override;

    /**
     * will first copy the data to
     * the gpu and then solve the system
     * and the data will be copied back to the host
     * store the result in x
     */
    void solve(
        Eigen::VectorXd &x,
        Eigen::VectorXd &b) override;

    void solve(
        double *x,
        double *b);

    /**
     * data from gpu
     * @todo due to current implementation
     * we need to copy the data from gpu to gpu
    */
    void solve_from_gpu(
        double *x,
        double *b);

    void solveWithGuess(
        Eigen::VectorXd &x,
        Eigen::VectorXd &b) override {};

private:
    /**
     * x and b is already on the gpu
     */
    void solve_nocp(double *x, double *b);


    const int max_iter = 1000;
    const double tol = 1e-6f;

    int N = 0;  // size of the matrix, rows
    int nz = 0; // number of non-zero elements
    // data on the host
    int *I = NULL;      // row indices
    int *J = NULL;      // column indices
    double *val = NULL; // values

    double *x, *rhs; // solution and right hand side

    //
    double a, b, na, r0, r1;

    // data on the device
    int *d_col = NULL;
    int *d_row = NULL;
    double *d_val = NULL;
    double *d_x = NULL, *d_r = NULL;
    double *d_p = NULL, *d_Ax = NULL;

    double dot;

    int k = 0; // iteration count
    double alpha, beta, alpham1;

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

#endif // __M_CUSOLVER_H__