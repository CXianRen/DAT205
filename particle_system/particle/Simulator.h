#pragma once
#include <vector>
#include <memory>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include "mmath.h"
#include "SimBase.h"

#include "Solver.h"
#include "CudaSimulator.h"

#include "mperf.h"

class Simulator
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  Simulator(double &time);
  ~Simulator();

  void update();

  const double *getDensity()
  {
    return density;
  }

  const double *getTransparency()
  {
    return transparency;
  }

  const std::string getPerformanceInfo()
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

  void reset()
  {
    for (int i = 0; i < SIZE; ++i)
    {
      density[i] = 0.0;
      density0[i] = 0.0;
      temperature[i] = T_AMBIENT;
      temperature0[i] = T_AMBIENT;
      pressure[i] = 0.0;
      u[i] = 0.0;
      u0[i] = 0.0;
      v[i] = 0.0;
      v0[i] = 0.0;
      w[i] = 0.0;
      w0[i] = 0.0;
    }

    m_time = 0.0;
  }

private:
  void setEmitterVelocity();
  void addSource();

  void calculateExternalForce();
  void calculateVorticity();
  void applyExternalForce();
  void calculatePressure();
  void applyPressure();
  void advectVelocity();
  void advectScalarField();

  /**   Semi-Lagarance method  **/
  double &m_time;

  // external force
  double fx[SIZE], fy[SIZE], fz[SIZE];
  // velocity field
  double u[SIZE], u0[SIZE];
  double v[SIZE], v0[SIZE];
  double w[SIZE], w0[SIZE];
  double avg_u[SIZE], avg_v[SIZE], avg_w[SIZE];

  // vorticity field
  double omg_x[SIZE], omg_y[SIZE], omg_z[SIZE];
  double vort[SIZE];

  // pressure field
  double pressure[SIZE];

  // temperature field
  double temperature0[SIZE], temperature[SIZE];

  // density field
  double density[SIZE], density0[SIZE];

  // transparency field
  double transparency[SIZE];
  double light_x, light_y, light_z;
  double module_scale_factor;
  double factor;

  // solver
  EigenSolver m_e_solver;
  Eigen::VectorXd b;
  Eigen::VectorXd x;

  CudaSolver m_solver;
  MCUDA::CudaSimulator CW;

  // ocuppied voxels
  std::array<bool, SIZE> m_occupied_voxels;
  void fixOccupiedVoxels();

  void genTransparencyMap();
};