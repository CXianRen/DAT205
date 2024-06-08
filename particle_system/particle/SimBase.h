#ifndef __MSIMBASE_H__
#define __MSIMBASE_H__

#include "mmath.h"
#include "constants.h"


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
        Nx + 1, Ny, Nz,
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
        Nx, Ny + 1, Nz,
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
        Nx, Ny, Nz + 1,
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

    u[ACC3D_X(i, j, k, Ny, Nx)] = getVelocityX<double>(
        pos_u,
        u_0,
        Nx, Ny, Nz);

    v[ACC3D_Y(i, j, k, Ny, Nx)] = getVelocityY<double>(
        pos_v,
        v_0,
        Nx, Ny, Nz);

    w[ACC3D_Z(i, j, k, Ny, Nx)] = getVelocityZ<double>(
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
#endif // __MSIMBASE_H__