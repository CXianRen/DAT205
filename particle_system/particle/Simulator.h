#pragma once
#include <vector>
#include <memory>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include "MData.h"
#include "mmath.h"

#include "Solver.h"

#include "mperf.h"

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

  std::string get_performance_info()
  {
    auto time_str = get_mesaurement_info();
    // get solver info
    int iter;
    double error;
    m_solver.getIterations(iter);
    m_solver.getError(error);
    return time_str + "Solver Iterations: " +
           std::to_string(iter) + " Solver Error: " + std::to_string(error) + "\n";
  }

  void setOccupiedVoxels(std::array<bool, SIZE> &occupied_voxels)
  {
    m_occupied_voxels = occupied_voxels;
  }

  void reset(){
    density.m_data.fill(0.0);
    density0.m_data.fill(0.0);
    temperature.m_data.fill(T_AMBIENT);
    temperature0.m_data.fill(T_AMBIENT);
    u.m_data.fill(0.0);
    u0.m_data.fill(0.0);
    v.m_data.fill(0.0);
    v0.m_data.fill(0.0);
    w.m_data.fill(0.0);
    w0.m_data.fill(0.0);
    pressure.m_data.fill(0.0);
    m_time = 0.0;
  }

private:
  void setEmitterVelocity();
  void addSource();

  /**     General step in update function         **/
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

  /**   Semi-Lagarance method  **/
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
  MDataX u, u0;
  MDataY v, v0;
  MDataZ w, w0;
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
  MDataScalar pressure;

  // temperature field
  MDataScalar temperature0, temperature;
  double getTemperature(const Vec3 &pos)
  {
    static int dim[3] = {Nx, Ny, Nz};
    static int maxIndex[3] = {Nx - 1, Ny - 1, Nz - 1};
    static Vec3 offset = 0.5 * Vec3(VOXEL_SIZE, VOXEL_SIZE, VOXEL_SIZE);
    return linearInterpolation<double>(
        pos - offset,
        temperature0.m_data.data(),
        dim,
        maxIndex);
  }

  // density field
  MDataScalar density, density0;
  double getDensity(const Vec3 &pos)
  {
    static int dim[3] = {Nx, Ny, Nz};
    static int maxIndex[3] = {Nx - 1, Ny - 1, Nz - 1};
    static Vec3 offset = 0.5 * Vec3(VOXEL_SIZE, VOXEL_SIZE, VOXEL_SIZE);
    return linearInterpolation<double>(
        pos - offset,
        density0.m_data.data(),
        dim,
        maxIndex);
  }

  // solver
  EigenSolver m_e_solver;

  Eigen::SparseMatrix<double, Eigen::RowMajor> A;
  Eigen::VectorXd b;
  Eigen::VectorXd x;

  CudaSolver m_solver;

  // ocuppied voxels
  std::array<bool, SIZE> m_occupied_voxels;

  void fix_occupied_voxels();
};

std::array<double, SIZE>&
generateSphereDensity();


std::array<double, SIZE>&
generateCubeDensity();