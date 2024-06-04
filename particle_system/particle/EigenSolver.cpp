#include "Solver.h"

void EigenSolver::solve(
    Eigen::VectorXd &x,
    Eigen::VectorXd &b)
{
    /* solve sparse lenear system by ICCG */

    // check if the matrix is positive definite
    if (ICCG.info() != Eigen::Success)
    {
        printf("decomposition failed\n");
        return;
    }
    x = ICCG.solve(b);
}

void EigenSolver::solveWithGuess(
    Eigen::VectorXd &x,
    Eigen::VectorXd &b)
{
    /* solve sparse lenear system by ICCG */

    // check if the matrix is positive definite
    if (ICCG.info() != Eigen::Success)
    {
        printf("decomposition failed\n");
        return;
    }
    x = ICCG.solveWithGuess(b, x);
}