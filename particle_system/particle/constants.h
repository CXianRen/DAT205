#pragma once
#include <cassert>

enum E_METHOD
{
    E_LINEAR = 0,
    E_MONOTONIC_CUBIC = 1
};

enum E_ADVECTION
{
    E_SEMI_LAGRANGE = 0,
    E_MAC_CORMACK = 1
};

enum E_EMITTER_POS
{
    E_TOP = 0,
    E_BOTTOM = 1
};

constexpr int N = 50;
constexpr int ratio[3] = {1, 2, 1}; // X, Y, Z
constexpr E_METHOD INTERPOLATION_METHOD = E_LINEAR;
constexpr E_ADVECTION ADVECTION_METHOD = E_SEMI_LAGRANGE;
constexpr E_EMITTER_POS EMITTER_POS = E_TOP;
constexpr bool OFFSCREEN_MODE = false;

// 32, 64, 32
constexpr int Nx = ratio[0] * N, Ny = ratio[1] * N, Nz = ratio[2] * N;
// the size of te sorce box
constexpr int SOURCE_SIZE_X = (int)(5);    // 8
constexpr int SOURCE_SIZE_Y = (int)(3);   // 3
constexpr int SOURCE_SIZE_Z = (int)(5);    // 8
constexpr int SOURCE_Y_MERGIN = (int)(3); // 3

constexpr int SIZE = Nx * Ny * Nz;

constexpr double DT = 0.02;
constexpr double VOXEL_SIZE = 1.0;
constexpr double INIT_DENSITY = 1.0;
constexpr double INIT_VELOCITY = 0.0;
constexpr double VORT_EPS = 1.0;
constexpr double ALPHA = 0.1;
constexpr double BETA = 5.0; 
constexpr double T_AMP = 5.0;
constexpr double T_AMBIENT = 25.0;
constexpr double EMIT_DURATION = 5.0;
constexpr double FINISH_TIME = 10.0;

// render 
constexpr float ABSORPTION = 5.0f;

/* other definitions */

#define POS(i, j, k) ((i) + Nx * (j) + Nx * Ny * (k))
#define POS_X(x, y, z) ((x) + (y) * (Nx + 1) + (z) * (Nx + 1) * Ny)
#define POS_Y(x, y, z) ((x) + (y) * Nx + (z) * Nx * (Ny + 1))
#define POS_Z(x, y, z) ((x) + (y) * Nx + (z) * Nx * Ny)


#define FOR_EACH_CELL                \
    for (int k = 0; k < Nz; ++k)     \
        for (int j = 0; j < Ny; ++j) \
            for (int i = 0; i < Nx; ++i)


#ifdef _OPENMP
#include <omp.h>
#define OPENMP_FOR _Pragma("omp parallel for")
#define OPENMP_FOR_COLLAPSE _Pragma("omp parallel for collapse(3)")
#else
#define OPENMP_FOR
#define OPENMP_FOR_COLLAPSE
#endif