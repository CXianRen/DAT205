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

constexpr int N = 32;
constexpr int ratio[3] = {1, 2, 1}; // X, Y, Z
constexpr E_METHOD INTERPOLATION_METHOD = E_MONOTONIC_CUBIC;
constexpr E_ADVECTION ADVECTION_METHOD = E_SEMI_LAGRANGE;
constexpr E_EMITTER_POS EMITTER_POS = E_TOP;
constexpr bool OFFSCREEN_MODE = false;

// 32, 64, 32
constexpr int Nx = ratio[0] * N, Ny = ratio[1] * N, Nz = ratio[2] * N;
// the size of te sorce box
constexpr int SOURCE_SIZE_X = (int)(Nx / 4);    // 8
constexpr int SOURCE_SIZE_Y = (int)(Ny / 20);   // 3
constexpr int SOURCE_SIZE_Z = (int)(Nz / 4);    // 8
constexpr int SOURCE_Y_MERGIN = (int)(Ny / 20); // 3

constexpr int SIZE = Nx * Ny * Nz;

constexpr double DT = 0.02;
constexpr double VOXEL_SIZE = 1.0;
constexpr double INIT_DENSITY = 1.0;
constexpr double INIT_VELOCITY = 80.0;
constexpr double VORT_EPS = 0.25;
constexpr double ALPHA = 9.8;
constexpr double BETA = 15.0;
constexpr double T_AMP = 5.0;
constexpr double T_AMBIENT = 50.0;
constexpr double EMIT_DURATION = 2.0;
constexpr double FINISH_TIME = 6.0;

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
