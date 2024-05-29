#pragma once
#include <vector>
#include <memory>
#include <Eigen/Core>
#include <Eigen/Sparse>

#include "mmath.h"

#include "GridData.h"

#include "Solver.h"

typedef Eigen::Triplet<double> T;

class Simulator
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Simulator(double &time);
  ~Simulator();

  void update();

  const std::array<double, SIZE> &getDensity()
  {
    return density.m_data;
  }

private:
  void setEmitterVelocity();
  void addSource();

  // a cube in the center of the domain
  void setOccupiedVoxels();

  /*
   * is to calculate the external force field
   * which is the buoyancy force and the gravity force
   * and store the result in fx, fy, fz
   * fx, fy, fz are the external force field
   * fx, fy, fz are the same size as the grid
   * fx, fy, fz are the force field in x, y, z direction
   * fx, fy, fz are the force field at the center of the grid
   * f_buoyancy = alpha * density * (0, 1, 0) + beta * (T - T_ambient) * (0, 1, 0)
   */
  void calculate_external_force();
  void calculate_vorticity();
  void apply_external_force();
  void calculate_pressure();
  void apply_pressure();
  void advect_velocity();
  void advect_scalar_field();

  double &m_time;

  Vec3 getCenter(int i, int j, int k)
  {
    double half_dx = 0.5 * VOXEL_SIZE;

    double x = half_dx + i * VOXEL_SIZE;
    double y = half_dx + j * VOXEL_SIZE;
    double z = half_dx + k * VOXEL_SIZE;
    return Vec3(x, y, z);
  }

  // external force
  double fx[SIZE], fy[SIZE], fz[SIZE];
  // velocity field
  GridDataX u, u0;
  GridDataY v, v0;
  GridDataZ w, w0;
  double avg_u[SIZE], avg_v[SIZE], avg_w[SIZE];

  Vec3 getVelocity(const Vec3 &pos)
  {
    Vec3 vel;
    vel[0] = getVelocityX(pos);
    vel[1] = getVelocityY(pos);
    vel[2] = getVelocityZ(pos);
    return vel;
  }

  double getVelocityX(const Vec3 &pos)
  {
    static int dim[3] = {Nx + 1, Ny, Nz};
    static int maxIndex[3] = {Nx, Ny - 1, Nz - 1};
    static Vec3 offset = 0.5 * Vec3(0.0, VOXEL_SIZE, VOXEL_SIZE);
    return linearInterpolation<double>(
        pos - offset,
        u0.m_data.data(),
        dim,
        maxIndex);
    // return u0.interp(pos - 0.5 * Vec3(0.0, VOXEL_SIZE, VOXEL_SIZE));
  }
  double getVelocityY(const Vec3 &pos)
  {
    static int dim[3] = {Nx, Ny + 1, Nz};
    static int maxIndex[3] = {Nx - 1, Ny, Nz - 1};
    static Vec3 offset = 0.5 * Vec3(VOXEL_SIZE, 0.0, VOXEL_SIZE);
    return linearInterpolation<double>(
        pos - offset,
        v0.m_data.data(),
        dim,
        maxIndex);
    // return v0.interp(pos - 0.5 * Vec3(VOXEL_SIZE, 0.0, VOXEL_SIZE));
  }
  double getVelocityZ(const Vec3 &pos)
  {
    static int dim[3] = {Nx, Ny, Nz + 1};
    static int maxIndex[3] = {Nx - 1, Ny - 1, Nz};
    static Vec3 offset = 0.5 * Vec3(VOXEL_SIZE, VOXEL_SIZE, 0.0);
    return linearInterpolation<double>(
        pos - offset,
        w0.m_data.data(),
        dim,
        maxIndex);
    // return w0.interp(pos - 0.5 * Vec3(VOXEL_SIZE, VOXEL_SIZE, 0.0));
  }

  // vorticity field
  double omg_x[SIZE], omg_y[SIZE], omg_z[SIZE];
  double vort[SIZE];

  // pressure field
  GridDataScalar pressure;
  // double getPressure(const Vec3 &pos)
  // {
  //   return pressure.interp(
  //       pos - 0.5 * Vec3(VOXEL_SIZE, VOXEL_SIZE, VOXEL_SIZE));
  // }

  // temperature field
  GridDataScalar temperature0, temperature;
  double getTemperature(const Vec3 &pos)
  {
    return temperature0.interp(
        pos - 0.5 * Vec3(VOXEL_SIZE, VOXEL_SIZE, VOXEL_SIZE));
  }

  // density field
  GridDataScalar density, density0;
  double getDensity(const Vec3 &pos)
  {
    return density0.interp(
        pos - 0.5 * Vec3(VOXEL_SIZE, VOXEL_SIZE, VOXEL_SIZE));
  }

  // solver
  std::vector<T> tripletList;
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> ICCG;

  Eigen::SparseMatrix<double, Eigen::RowMajor> A;
  Eigen::VectorXd b;
  Eigen::VectorXd x;

  CudaSolver m_solver;
};