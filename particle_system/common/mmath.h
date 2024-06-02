#ifndef __MMATH_H__
#define __MMATH_H__
#include <iostream>
#include <vector>
#include <math.h>
#include <Eigen/Core>
#include <Eigen/Sparse>

#include "constants.h"
#include "Vec3.h"

#define POS(i, j, k) ((i) + Nx * (j) + Nx * Ny * (k))
#define POS_X(i, j, k) ((i) + (Nx + 1) * (j) + (Nx + 1) * (Ny) * (k))
#define POS_Y(i, j, k) ((i) + Nx * (j) + Nx * (Ny + 1) * (k))
#define POS_Z(i, j, k) ((i) + Nx * (j) + Nx * Ny * (k))

#define FOR_EACH_CELL                \
    for (int k = 0; k < Nz; ++k)     \
        for (int j = 0; j < Ny; ++j) \
            for (int i = 0; i < Nx; ++i)

#define ACCESS3D(x, y, z) ((x) + (y) * dims[0] + (z) * dims[0] * dims[1])
#define ACCESS3D_X(x, y, z) ((x) + (y) * (dims[0] + 1) + (z) * (dims[0] + 1) * dims[1])
#define ACCESS3D_Y(x, y, z) ((x) + (y) * dims[0] + (z) * dims[0] * (dims[1] + 1))
#define ACCESS3D_Z(x, y, z) ((x) + (y) * dims[0] + (z) * dims[0] * dims[1])

template <typename T>
T linearInterpolation(const T *pt, T *src, int *dims, int *maxXYZ)
{
    T pos[3];
    // clamp position
    pos[0] = std::min(
        std::max((T)0.0, pt[0]),
        (T)(maxXYZ[0]) * VOXEL_SIZE - (T)1e-6);
    pos[1] = std::min(
        std::max((T)0.0, pt[1]),
        (T)(maxXYZ[1]) * VOXEL_SIZE - (T)1e-6);
    pos[2] = std::min(
        std::max((T)0.0, pt[2]),
        (T)(maxXYZ[2]) * VOXEL_SIZE - (T)1e-6);

    int i = (int)(pos[0] / VOXEL_SIZE);
    int j = (int)(pos[1] / VOXEL_SIZE);
    int k = (int)(pos[2] / VOXEL_SIZE);

    T scale = 1.0 / VOXEL_SIZE;
    T fractx = scale * (pos[0] - i * VOXEL_SIZE);
    T fracty = scale * (pos[1] - j * VOXEL_SIZE);
    T fractz = scale * (pos[2] - k * VOXEL_SIZE);

    assert(fractx < 1.0 && fractx >= 0);
    assert(fracty < 1.0 && fracty >= 0);
    assert(fractz < 1.0 && fractz >= 0);

    // Y @ low X, low Z:
    T tmp1 = src[ACCESS3D(i, j, k)];
    T tmp2 = src[ACCESS3D(i, j + 1, k)];
    // Y @ high X, low Z:
    T tmp3 = src[ACCESS3D(i + 1, j, k)];
    T tmp4 = src[ACCESS3D(i + 1, j + 1, k)];

    // Y @ low X, high Z:
    T tmp5 = src[ACCESS3D(i, j, k + 1)];
    T tmp6 = src[ACCESS3D(i, j + 1, k + 1)];
    // Y @ high X, high Z:
    T tmp7 = src[ACCESS3D(i + 1, j, k + 1)];
    T tmp8 = src[ACCESS3D(i + 1, j + 1, k + 1)];

    // Y @ low X, low Z
    T tmp12 = ((T)(1) - fracty) * tmp1 + fracty * tmp2;
    // Y @ high X, low Z
    T tmp34 = ((T)(1) - fracty) * tmp3 + fracty * tmp4;

    // Y @ low X, high Z
    T tmp56 = ((T)(1) - fracty) * tmp5 + fracty * tmp6;
    // Y @ high X, high Z
    T tmp78 = ((T)(1) - fracty) * tmp7 + fracty * tmp8;

    // X @ low Z
    T tmp1234 = ((T)(1) - fractx) * tmp12 + fractx * tmp34;
    // X @ high Z
    T tmp5678 = ((T)(1) - fractx) * tmp56 + fractx * tmp78;

    // Z
    T tmp = ((T)(1) - fractz) * tmp1234 + fractz * tmp5678;
    return tmp;
}

template <typename T>
T linearInterpolation(const Vec3 &pt, T *src, int *dims, int *maxXYZ)
{
    return linearInterpolation<T>(pt.n, src, dims, maxXYZ);
}


#define ACC3D(x, y, z, rows, cols) ((x) + (y) * cols + (z) * cols * rows)

#define ACC2D(x, y, cols) ACC3D(x, y, 0, 0, cols)

/*
    middle difference
    du(i,j) = (u(i+1,j) - u(i-1,j)) / 2
*/

#define GET_GRADIANT_3D_X(i, j, k, data, rows, cols, dx) \
    (data[ACC3D(i + 1, j, k, rows, cols)] - data[ACC3D(i - 1, j, k, rows, cols)]) / 2 / dx

#define GET_GRADIANT_3D_Y(i, j, k, data, rows, cols, dy) \
    (data[ACC3D(i, j + 1, k, rows, cols)] - data[ACC3D(i, j - 1, k, rows, cols)]) / 2 / dy

#define GET_GRADIANT_3D_Z(i, j, k, data, rows, cols, dz) \
    (data[ACC3D(i, j, k + 1, rows, cols)] - data[ACC3D(i, j, k - 1, rows, cols)]) / 2 / dz

#define GET_GRADIANT_2D_X(i, j, data, cols, dx) \
    GET_GRADIANT_3D_X(i, j, 0, data, 0, cols, dx)

#define GET_GRADIANT_2D_Y(i, j, data, cols, dy) \
    GET_GRADIANT_3D_Y(i, j, 0, data, 0, cols, dy)

#define GET_GRADIANT_2D_Z(i, j, data, cols, dz) \
    GET_GRADIANT_3D_Z(i, j, 0, data, 0, cols, dz)

// calculate divergence \nobal /dot u

#define GET_DIVERGENCE_2D_X(i, j, u, cols, dx) \
    (u[ACC2D(i + 1, j, cols)] - u[ACC2D(i, j, cols)] / dx)

#define GET_DIVERGENCE_2D_Y(i, j, v, cols, dy) \
    (v[ACC2D(i, j + 1, cols)] - v[ACC2D(i, j, cols)] / dy)

#define GET_DIVERGENCE_2D(i, j, u, v, cols, dx, dy) \
    (GET_DIVERGENCE_2D_X(i, j, u, cols, dx) + GET_DIVERGENCE_2D_Y(i, j, v, cols, dy))

// build 2d laplace matrix
// dirichlet boundary condition
// N: number of grid points
// return: laplace matrix
/*
    a example of laplace matrix:

    2 -1  0 -1  0  0  0  0  0
   -1  3 -1  0 -1  0  0  0  0
    0 -1  2  0  0 -1  0  0  0
   -1  0  0  3 -1  0 -1  0  0
    0 -1  0 -1  4 -1  0 -1  0
    0  0 -1  0 -1  3  0  0 -1
    0  0  0 -1  0  0  2 -1  0
    0  0  0  0 -1  0 -1  3 -1
    0  0  0  0  0 -1  0 -1  2
*/
template <typename T>
Eigen::SparseMatrix<T, Eigen::RowMajor> build_2d_laplace(int Nx, int Ny)
{
    std::vector<Eigen::Triplet<T>> tripletList;
    tripletList.reserve(Nx * Ny * 5);

    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            int idx = ACC2D(i, j, Nx);
            T sum = 0.0;

            if (i > 0)
            {
                tripletList.push_back(Eigen::Triplet<T>(idx, ACC2D(i - 1, j, Nx), -1.0));
                sum += 1.0f;
            }
            if (i < Nx - 1)
            {
                tripletList.push_back(Eigen::Triplet<T>(idx, ACC2D(i + 1, j, Nx), -1.0));
                sum += 1.0f;
            }
            if (j > 0)
            {
                tripletList.push_back(Eigen::Triplet<T>(idx, ACC2D(i, j - 1, Nx), -1.0));
                sum += 1.0f;
            }
            if (j < Ny - 1)
            {
                tripletList.push_back(Eigen::Triplet<T>(idx, ACC2D(i, j + 1, Nx), -1.0));
                sum += 1.0f;
            }
            tripletList.push_back(Eigen::Triplet<T>(idx, ACC2D(i, j, Nx), sum));
        }
    }

    Eigen::SparseMatrix<T, Eigen::RowMajor> A(Nx * Ny, Nx * Ny);
    A.setFromTriplets(tripletList.begin(), tripletList.end());

    return A;
}

// build 3d laplace matrix
// dirichlet boundary condition
// N: number of grid points
// return: laplace matrix
template <typename T>
Eigen::SparseMatrix<T, Eigen::RowMajor> build_3d_laplace(int Nx, int Ny, int Nz)
{
    int M = Nx * Ny * Nz;
    std::vector<Eigen::Triplet<T>> tripletList;
    tripletList.reserve(M * 7);
    for (int k = 0; k < Nz; k++)
        for (int j = 0; j < Ny; j++)
            for (int i = 0; i < Nx; i++)
            {
                int idx = ACC3D(i, j, k, Ny, Nx);
                T sum = 0.0;

                if (i > 0)
                {
                    tripletList.push_back(
                        Eigen::Triplet<T>(
                            idx,
                            ACC3D(i - 1, j, k, Ny, Nx),
                            1.0));
                    sum += 1.0f;
                }
                if (i < Nx - 1)
                {
                    tripletList.push_back(
                        Eigen::Triplet<T>(
                            idx,
                            ACC3D(i + 1, j, k, Ny, Nx),
                            1.0));
                    sum += 1.0f;
                }
                if (j > 0)
                {
                    tripletList.push_back(
                        Eigen::Triplet<T>(
                            idx,
                            ACC3D(i, j - 1, k, Ny, Nx),
                            1.0));
                    sum += 1.0f;
                }
                if (j < Ny - 1)
                {
                    tripletList.push_back(
                        Eigen::Triplet<T>(
                            idx,
                            ACC3D(i, j + 1, k, Ny, Nx),
                            1.0));
                    sum += 1.0f;
                }
                if (k > 0)
                {
                    tripletList.push_back(
                        Eigen::Triplet<T>(
                            idx,
                            ACC3D(i, j, k - 1, Ny, Nx),
                            1.0));
                    sum += 1.0f;
                }
                if (k < Nz - 1)
                {
                    tripletList.push_back(
                        Eigen::Triplet<T>(
                            idx,
                            ACC3D(i, j, k + 1, Ny, Nx),
                            1.0));
                    sum += 1.0f;
                }
                tripletList.push_back(
                    Eigen::Triplet<T>(
                        idx,
                        ACC3D(i, j, k, Ny, Nx),
                        -sum));
            }
    auto A = Eigen::SparseMatrix<T, Eigen::RowMajor>(M, M);
    A.setFromTriplets(tripletList.begin(), tripletList.end());
    return A;
}

#endif // __MMATH_H__