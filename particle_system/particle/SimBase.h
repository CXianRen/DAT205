#ifndef __MSIMBASE_H__
#define __MSIMBASE_H__

#include "mmath.h"
#include "constants.h"

template <typename T>
PREFIX inline void
calculateBuyancyForceBody(
    int i, int j, int k,
    int Nx, int Ny, int Nz,
    T *density, T *temperature,
    T *f_x, T *f_y, T *f_z)
{
    f_x[ACC3D(i, j, k, Ny, Nx)] = 0.0;
    f_y[ACC3D(i, j, k, Ny, Nx)] =
        -ALPHA * density[ACC3D(i, j, k, Ny, Nx)] +
        BETA * (temperature[ACC3D(i, j, k, Ny, Nx)] - T_AMBIENT);
    f_z[ACC3D(i, j, k, Ny, Nx)] = 0.0;
}

template <typename T>
PREFIX inline void
getCenter(
    int i, int j, int k,
    T *center)
{
    center[0] = 0.5 * VOXEL_SIZE + i * VOXEL_SIZE;
    center[1] = 0.5 * VOXEL_SIZE + j * VOXEL_SIZE;
    center[2] = 0.5 * VOXEL_SIZE + k * VOXEL_SIZE;
}

template <typename T>
PREFIX
    T
    getVelocityX(
        T *pos,
        T *u,
        int Nx, int Ny, int Nz)
{
    T tpos[3] = {
        pos[0],
        pos[1] - 0.5 * VOXEL_SIZE,
        pos[2] - 0.5 * VOXEL_SIZE};

    return linearInterpolation3D<T>(
        tpos,
        u,
        Nx, Ny, Nz,
        Nx, Ny - 1, Nz - 1, VOXEL_SIZE);
}

template <typename T>
PREFIX
    T
    getVelocityY(
        T *pos,
        T *v,
        int Nx, int Ny, int Nz)
{
    T tpos[3] = {
        pos[0] - 0.5 * VOXEL_SIZE,
        pos[1],
        pos[2] - 0.5 * VOXEL_SIZE};

    return linearInterpolation3D<T>(
        tpos,
        v,
        Nx, Ny, Nz,
        Nx - 1, Ny, Nz - 1, VOXEL_SIZE);
}

template <typename T>
PREFIX
    T
    getVelocityZ(
        T *pos,
        T *w,
        int Nx, int Ny, int Nz)
{
    T tpos[3] = {
        pos[0] - 0.5 * VOXEL_SIZE,
        pos[1] - 0.5 * VOXEL_SIZE,
        pos[2]};
    return linearInterpolation3D<T>(
        tpos,
        w,
        Nx, Ny, Nz,
        Nx - 1, Ny - 1, Nz, VOXEL_SIZE);
}

template <typename T>
PREFIX void getVelocity(
    T *pos,
    T *vel,
    T *u, T *v, T *w,
    int Nx, int Ny, int Nz)
{
    vel[0] = getVelocityX<T>(pos, u, Nx, Ny, Nz);
    vel[1] = getVelocityY<T>(pos, v, Nx, Ny, Nz);
    vel[2] = getVelocityZ<T>(pos, w, Nx, Ny, Nz);
}

template <typename T>
PREFIX inline void calculateAverageVelocity(
    int i, int j, int k,
    int Nx, int Ny, int Nz,
    T *u, T *v, T *w,
    T *avg_u, T *avg_v, T *avg_w)
{
    // calculate average velocity
    T t0 = 0, t1 = 0;

    t0 = u[ACC3D(i, j, k, Ny, Nx)];
    if (i < Nx - 1)
    {
        t1 = u[ACC3D(i + 1, j, k, Ny, Nx)];
    }
    avg_u[ACC3D(i, j, k, Ny, Nx)] = (t0 + t1) * 0.5;

    t0 = v[ACC3D(i, j, k, Ny, Nx)];
    t1 = 0;
    if (j < Ny - 1)
    {
        t1 = v[ACC3D(i, j + 1, k, Ny, Nx)];
    }
    avg_v[ACC3D(i, j, k, Ny, Nx)] = (t0 + t1) * 0.5;

    t0 = w[ACC3D(i, j, k, Ny, Nx)];
    t1 = 0;
    if (k < Nz - 1)
    {

        t1 = w[ACC3D(i, j, k + 1, Ny, Nx)];
    }
    avg_w[ACC3D(i, j, k, Ny, Nx)] = (t0 + t1) * 0.5;
}

template <typename T>
PREFIX inline void calculateGradient(
    int i, int j, int k,
    int Nx, int Ny, int Nz,
    T *avg_u, T *avg_v, T *avg_w,
    T *grad_u, T *grad_v, T *grad_w)
{
    // ignore boundary cells
    if (i == 0 || j == 0 || k == 0)
    {
        return;
    }
    if (i == Nx - 1 || j == Ny - 1 || k == Nz - 1)
    {
        return;
    }
    // calculate vorticity
    grad_u[ACC3D(i, j, k, Ny, Nx)] = (avg_w[ACC3D(i, j + 1, k, Ny, Nx)] - avg_w[ACC3D(i, j - 1, k, Ny, Nx)] - avg_v[ACC3D(i, j, k + 1, Ny, Nx)] + avg_v[ACC3D(i, j, k - 1, Ny, Nx)]) * 0.5 / VOXEL_SIZE;
    grad_v[ACC3D(i, j, k, Ny, Nx)] = (avg_u[ACC3D(i, j, k + 1, Ny, Nx)] - avg_u[ACC3D(i, j, k - 1, Ny, Nx)] - avg_w[ACC3D(i + 1, j, k, Ny, Nx)] + avg_w[ACC3D(i - 1, j, k, Ny, Nx)]) * 0.5 / VOXEL_SIZE;
    grad_w[ACC3D(i, j, k, Ny, Nx)] = (avg_v[ACC3D(i + 1, j, k, Ny, Nx)] - avg_v[ACC3D(i - 1, j, k, Ny, Nx)] - avg_u[ACC3D(i, j + 1, k, Ny, Nx)] + avg_u[ACC3D(i, j - 1, k, Ny, Nx)]) * 0.5 / VOXEL_SIZE;
}

template <typename T>
PREFIX inline void applyExternalForceBody(
    int i, int j, int k,
    int Nx, int Ny, int Nz,
    T *f_x, T *f_y, T *f_z,
    T *u, T *v, T *w)
{
    if (i < Nx - 1)
    {
        u[ACC3D(i + 1, j, k, Ny, Nx)] += DT * (f_x[ACC3D(i, j, k, Ny, Nx)] + f_x[ACC3D(i + 1, j, k, Ny, Nx)]) * 0.5;
    }
    if (j < Ny - 1)
    {
        v[ACC3D(i, j + 1, k, Ny, Nx)] += DT * (f_y[ACC3D(i, j, k, Ny, Nx)] + f_x[ACC3D(i, j + 1, k, Ny, Nx)]) * 0.5;
    }
    if (k < Nz - 1)
    {
        w[ACC3D(i, j, k + 1, Ny, Nx)] += DT * (f_z[ACC3D(i, j, k, Ny, Nx)] + f_x[ACC3D(i, j, k + 1, Ny, Nx)]) * 0.5;
    }
}

template <typename T>
PREFIX inline void applyPressureBody(
    int i, int j, int k,
    int Nx, int Ny, int Nz,
    T *pressure,
    T *u, T *v, T *w)
{
    // compute gradient of pressure
    // @todo optimize this (VOXEL_SIZE / DT) * DT
    if (i < Nx - 1)
    {
        u[ACC3D(i + 1, j, k, Ny, Nx)] -=
            (VOXEL_SIZE / DT) * DT * (pressure[ACC3D(i + 1, j, k, Ny, Nx)] - pressure[ACC3D(i, j, k, Ny, Nx)]) / VOXEL_SIZE;
    }
    if (j < Ny - 1)
    {
        v[ACC3D(i, j + 1, k, Ny, Nx)] -=
            (VOXEL_SIZE / DT) * DT * (pressure[ACC3D(i, j + 1, k, Ny, Nx)] - pressure[ACC3D(i, j, k, Ny, Nx)]) / VOXEL_SIZE;
    }
    if (k < Nz - 1)
    {
        w[ACC3D(i, j, k + 1, Ny, Nx)] -=
            (VOXEL_SIZE / DT)* DT * (pressure[ACC3D(i, j, k + 1, Ny, Nx)] - pressure[ACC3D(i, j, k, Ny, Nx)]) / VOXEL_SIZE;
    }
}

template <typename T>
PREFIX inline void calculateVorticityBody(
    int i, int j, int k,
    int Nx, int Ny, int Nz,
    T *omg_x, T *omg_y, T *omg_z,
    T *fx, T *fy, T *fz)
{

    // ignore boundary cells
    if (i == 0 || j == 0 || k == 0)
    {
        return;
    }
    if (i == Nx - 1 || j == Ny - 1 || k == Nz - 1)
    {
        return;
    }
    // compute gradient of vorticity
    T p, q;
    p = VEC3_NORM(omg_x[ACC3D(i + 1, j, k, Ny, Nx)], omg_y[ACC3D(i + 1, j, k, Ny, Nx)], omg_z[ACC3D(i + 1, j, k, Ny, Nx)]);
    q = VEC3_NORM(omg_x[ACC3D(i - 1, j, k, Ny, Nx)], omg_y[ACC3D(i - 1, j, k, Ny, Nx)], omg_z[ACC3D(i - 1, j, k, Ny, Nx)]);
    T grad1 = (p - q) * 0.5 / VOXEL_SIZE;

    p = VEC3_NORM(omg_x[ACC3D(i, j + 1, k, Ny, Nx)], omg_y[ACC3D(i, j + 1, k, Ny, Nx)], omg_z[ACC3D(i, j + 1, k, Ny, Nx)]);
    q = VEC3_NORM(omg_x[ACC3D(i, j - 1, k, Ny, Nx)], omg_y[ACC3D(i, j - 1, k, Ny, Nx)], omg_z[ACC3D(i, j - 1, k, Ny, Nx)]);
    T grad2 = (p - q) * 0.5 / VOXEL_SIZE;

    p = VEC3_NORM(omg_x[ACC3D(i, j, k + 1, Ny, Nx)], omg_y[ACC3D(i, j, k + 1, Ny, Nx)], omg_z[ACC3D(i, j, k + 1, Ny, Nx)]);
    q = VEC3_NORM(omg_x[ACC3D(i, j, k - 1, Ny, Nx)], omg_y[ACC3D(i, j, k - 1, Ny, Nx)], omg_z[ACC3D(i, j, k - 1, Ny, Nx)]);
    T grad3 = (p - q) * 0.5 / VOXEL_SIZE;

    T norm = VEC3_NORM(grad1, grad2, grad3);

    T ni = 0.0, nj = 0.0, nk = 0.0;

    if (norm != 0)
    {
        ni = grad1 / norm;
        nj = grad2 / norm;
        nk = grad3 / norm;
    }

    T f1, f2, f3;

    VEC3_CROSS(
        omg_x[ACC3D(i, j, k, Ny, Nx)],
        omg_y[ACC3D(i, j, k, Ny, Nx)],
        omg_z[ACC3D(i, j, k, Ny, Nx)],
        ni, nj, nk,
        f1, f2, f3);

    fx[ACC3D(i, j, k, Ny, Nx)] += VORT_EPS * VOXEL_SIZE * f1;
    fy[ACC3D(i, j, k, Ny, Nx)] += VORT_EPS * VOXEL_SIZE * f2;
    fz[ACC3D(i, j, k, Ny, Nx)] += VORT_EPS * VOXEL_SIZE * f3;
}

template <typename T>
PREFIX inline void advectVelocityBody(
    T *u, T *v, T *w,
    T *u_0, T *v_0, T *w_0,
    int i, int j, int k,
    int Nx, int Ny, int Nz)
{
    double half_dx = 0.5 * VOXEL_SIZE;

    double center[3];
    center[0] = half_dx + i * VOXEL_SIZE;
    center[1] = half_dx + j * VOXEL_SIZE;
    center[2] = half_dx + k * VOXEL_SIZE;

    double pos_u[3];
    pos_u[0] = center[0] - half_dx;
    pos_u[1] = center[1];
    pos_u[2] = center[2];

    double pos_v[3];
    pos_v[0] = center[0];
    pos_v[1] = center[1] - half_dx;
    pos_v[2] = center[2];

    double pos_w[3];
    pos_w[0] = center[0];
    pos_w[1] = center[1];
    pos_w[2] = center[2] - half_dx;

    double vel_u[3];
    getVelocity<double>(
        pos_u,
        vel_u,
        u_0,
        v_0,
        w_0,
        Nx, Ny, Nz);

    double vel_v[3];
    getVelocity<double>(
        pos_v,
        vel_v,
        u_0,
        v_0,
        w_0,
        Nx, Ny, Nz);

    double vel_w[3];
    getVelocity<double>(
        pos_w,
        vel_w,
        u_0,
        v_0,
        w_0,
        Nx, Ny, Nz);

    pos_u[0] -= DT * vel_u[0];
    pos_u[1] -= DT * vel_u[1];
    pos_u[2] -= DT * vel_u[2];

    pos_v[0] -= DT * vel_v[0];
    pos_v[1] -= DT * vel_v[1];
    pos_v[2] -= DT * vel_v[2];

    pos_w[0] -= DT * vel_w[0];
    pos_w[1] -= DT * vel_w[1];
    pos_w[2] -= DT * vel_w[2];

    u[ACC3D(i, j, k, Ny, Nx)] = getVelocityX<double>(
        pos_u,
        u_0,
        Nx, Ny, Nz);

    v[ACC3D(i, j, k, Ny, Nx)] = getVelocityY<double>(
        pos_v,
        v_0,
        Nx, Ny, Nz);

    w[ACC3D(i, j, k, Ny, Nx)] = getVelocityZ<double>(
        pos_w,
        w_0,
        Nx, Ny, Nz);
}

template <typename T>
PREFIX T getScalar(
    T *pos,
    T *src,
    int Nx, int Ny, int Nz)
{
    T tpos[3] = {
        pos[0] - 0.5 * VOXEL_SIZE,
        pos[1] - 0.5 * VOXEL_SIZE,
        pos[2] - 0.5 * VOXEL_SIZE};

    return linearInterpolation3D<T>(
        tpos,
        src,
        Nx, Ny, Nz,
        Nx - 1, Ny - 1, Nz - 1, VOXEL_SIZE);
}

template <typename T>
PREFIX inline void advectScalarBody(
    int i, int j, int k,
    int Nx, int Ny, int Nz,
    T *field, T *field_0,
    T *u_0, T *v_0, T *w_0)
{
    T pos_cell[3];
    getCenter<double>(i, j, k, pos_cell);

    T vel_cell[3];
    getVelocity<double>(
        pos_cell,
        vel_cell,
        u_0,
        v_0,
        w_0,
        Nx, Ny, Nz);

    pos_cell[0] -= vel_cell[0] * DT;
    pos_cell[1] -= vel_cell[1] * DT;
    pos_cell[2] -= vel_cell[2] * DT;

    field[ACC3D(i, j, k, Ny, Nx)] = getScalar<double>(
        pos_cell,
        field_0,
        Nx, Ny, Nz);
}

template <typename T, typename U>
PREFIX inline void applyOccupiedVoxelsBody(
    int i, int j, int k,
    int Nx, int Ny, int Nz,
    U occupied_voxels,
    T *density, T *u, T *v, T *w, T *temperature)
{
    if (occupied_voxels[ACC3D(i, j, k, Ny, Nx)])
    {
        u[ACC3D(i, j, k, Ny, Nx)] = 0.0;
        v[ACC3D(i, j, k, Ny, Nx)] = 0.0;
        w[ACC3D(i, j, k, Ny, Nx)] = 0.0;
        temperature[ACC3D(i, j, k, Ny, Nx)] = T_AMBIENT;
        density[ACC3D(i, j, k, Ny, Nx)] = 0.0;
    }
}

template <typename T>
PREFIX inline void buildRhsBody(
    int i, int j, int k,
    int Nx, int Ny, int Nz,
    T *u, T *v, T *w,
    T *b)
{
    double F[6] = {
        static_cast<double>(k > 0),
        static_cast<double>(j > 0),
        static_cast<double>(i > 0),
        static_cast<double>(i < Nx - 1),
        static_cast<double>(j < Ny - 1),
        static_cast<double>(k < Nz - 1)};

    double U[6];

    U[0] = w[ACC3D(i, j, k, Ny, Nx)];
    U[1] = v[ACC3D(i, j, k, Ny, Nx)];
    U[2] = u[ACC3D(i, j, k, Ny, Nx)];

    U[3] = 0.0;
    if (i < Nx - 1)
        U[3] = (double)(u[ACC3D(i + 1, j, k, Ny, Nx)]);

    U[4] = 0.0;
    if (j < Ny - 1)
        U[4] = (double)(v[ACC3D(i, j + 1, k, Ny, Nx)]);

    U[5] = 0.0;
    if (k < Nz - 1)
        U[5] = (double)(w[ACC3D(i, j, k + 1, Ny, Nx)]);

    for (int n = 0; n < 3; ++n)
    {
        b[ACC3D(i, j, k, Ny, Nx)] -= F[n] * U[n];
    }
    for (int n = 3; n < 6; ++n)
    {
        b[ACC3D(i, j, k, Ny, Nx)] += F[n] * U[n];
    }
}

#endif // __MSIMBASE_H__