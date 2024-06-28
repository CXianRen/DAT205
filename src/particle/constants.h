#pragma once
#include <cassert>


constexpr int N = 50;
constexpr int ratio[3] = {1, 1, 1}; // X, Y, Z
constexpr int Nx = ratio[0] * N, Ny = ratio[1] * N, Nz = ratio[2] * N;
constexpr int SIZE = Nx * Ny * Nz;

// the size of te sorce box

constexpr double DT = 0.01;

constexpr double VOXEL_SIZE = 1.0;
constexpr double INIT_DENSITY = 1.0;
constexpr double INIT_VELOCITY = 0.0;
constexpr double T_AMBIENT = 25.0;

constexpr double VORT_EPS = 30.0;
constexpr double ALPHA = 0.1;
constexpr double BETA = 5.0;

constexpr double EMIT_DURATION = 5.0;
constexpr double FINISH_TIME = 10.0;