#pragma once
#include <cassert>

enum E_EMITTER_POS
{
    E_TOP = 0,
    E_BOTTOM = 1
};

constexpr int N = 64;
constexpr int ratio[3] = {1, 1, 1}; // X, Y, Z
constexpr int Nx = ratio[0] * N, Ny = ratio[1] * N, Nz = ratio[2] * N;
constexpr int SIZE = Nx * Ny * Nz;

constexpr E_EMITTER_POS EMITTER_POS = E_TOP;

// the size of te sorce box
constexpr int SOURCE_SIZE_X = (int)(10);
constexpr int SOURCE_SIZE_Y = (int)(1);
constexpr int SOURCE_SIZE_Z = (int)(10);
constexpr int SOURCE_Y_MERGIN = (int)(3);

constexpr double DT = 0.02;

constexpr double VOXEL_SIZE = 1.0;
constexpr double INIT_DENSITY = 1.0;
constexpr double INIT_VELOCITY = 10.0;
constexpr double T_AMBIENT = 25.0;

constexpr double VORT_EPS = 10.0;
constexpr double ALPHA = 0.1;
constexpr double BETA = 5.0;

constexpr double EMIT_DURATION = 5.0;
constexpr double FINISH_TIME = 10.0;

// render
constexpr float ABSORPTION = 5.0f;
