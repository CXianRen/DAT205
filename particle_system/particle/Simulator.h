#pragma once
#include <vector>
#include <memory>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include "MData.h"
#include "mmath.h"

#include "SimBase.h"

#include "Solver.h"
#include "CudaWorker.h"

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

  void setLightPosition(
      double x, double y, double z,
      double scale_factor, double factor)
  {
    light_x = x;
    light_y = y;
    light_z = z;
    module_scale_factor = scale_factor;
    this->factor = factor;
  }

  void genTransparencyMap();

  double *getTransparency()
  {
    return transparency;
  }

  void reset()
  {
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

  // vorticity field
  double omg_x[SIZE], omg_y[SIZE], omg_z[SIZE];
  double vort[SIZE];

  // pressure field
  MDataScalar pressure;

  // temperature field
  MDataScalar temperature0, temperature;

  // density field
  MDataScalar density, density0;

  double transparency[SIZE];
  double light_x, light_y, light_z;
  double module_scale_factor;
  double factor;

  // solver
  EigenSolver m_e_solver;

  Eigen::SparseMatrix<double, Eigen::RowMajor> A;
  Eigen::VectorXd b;
  Eigen::VectorXd x;

  CudaSolver m_solver;

  MCUDA::CudaWorker CW;

  // ocuppied voxels
  std::array<bool, SIZE> m_occupied_voxels;

  void fix_occupied_voxels();
};

std::array<double, SIZE> &
generateSphereDensity();

std::array<double, SIZE> &
generateCubeDensity();