#include "test.h"
#include "test_common.h"

#include <common/mmath.h>

#include <particle/Solver.h>

void test_eigen_solver()
{
    TEST_LOG("Test Eigen Solver");

    int N = 32;
    auto L = build_2d_laplace<float>(N, N);

    Eigen::VectorXf x = Eigen::VectorXf::Zero(N * N);
    Eigen::VectorXf b = Eigen::VectorXf::Zero(N * N);

    Eigen::VectorXf x_true = Eigen::VectorXf::Zero(N * N);
    // init x_true with sin function
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            int k = i * N + j;

            x_true(k) = std::sin(M_PI * i / (N - 1)) * std::sin(M_PI * j / (N - 1));
        }
    }
    // TEST_LOG("\t\t head(10): \n\t\t" << x_true.head(10).transpose());
    // compute b
    b = L * x_true;

    // init solver
    EigenSolver solver;
    TEST_LOG("\t\tTest compute");
    {
        TIME_START
        solver.compute(L);
        TIME_END("\t\tCompute in :")
    }

    // solve
    TEST_LOG("\t\tTest solve without guess");
    {
        TIME_START
        solver.solve(x, b);
        TIME_END("\t\tSolve in :")
    }

    int iter;
    float error;
    solver.getIterations(iter);
    solver.getError(error);

    TEST_LOG("\t\titer: " << iter << ", error: " << error);

    auto bias = 0.1 * Eigen::VectorXf::Random(N * N);

    x = x + bias;
    // show fist ten elements in x
    // TEST_LOG("x.head(10): \n\t\t" << x.head(10).transpose());
    // solve with guess

    TEST_LOG("\t\tTest solve with guess");
    {
        TIME_START
        solver.solveWithGuess(x, b);
        TIME_END("\t\tSolve in :")
    }
    solver.getIterations(iter);
    solver.getError(error);

    TEST_LOG("\t\titer: " << iter << ", error: " << error);

    // show fist ten elements in x
    // TEST_LOG("x.head(10): \n\t\t" << x.head(10).transpose());

    TEST_LOG("Eigen Solver passed");
}

void test_cuda_solver()
{
    TEST_LOG("Test Cuda Solver");

    int N = 5;
    auto L = build_2d_laplace<float>(N, N);

    Eigen::VectorXf x = Eigen::VectorXf::Zero(N * N);
    Eigen::VectorXf b = Eigen::VectorXf::Zero(N * N);

    Eigen::VectorXf x_true = Eigen::VectorXf::Zero(N * N);
    // init x_true with sin function
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            int k = i * N + j;

            x_true(k) = std::cos(M_PI * i / (N - 1));
        }
    }
    TEST_LOG("\t\t x true head(10): \n\t\t" << x_true.head(10).transpose());
    // compute b
    b = L * x_true;

    // init solver
    CudaSolver solver;
    TEST_LOG("\t\tTest compute");
    {
        TIME_START
        solver.compute(L);
        TIME_END("\t\tCompute in :")
    }

    // solve
    TEST_LOG("\t\tTest solve without guess");
    {
        TIME_START
        solver.solve(x, b);
        TIME_END("\t\tSolve in :")
    }

    int iter;
    float error;
    solver.getIterations(iter);
    solver.getError(error);

    TEST_LOG("\t\titer: " << iter << ", error: " << error);
    TEST_LOG("\t\t x head(10): \n\t\t" << x.head(10).transpose());

    TEST_LOG("Cuda Solver passed");
}

void test_solver()
{
    test_eigen_solver();
    test_cuda_solver();
}