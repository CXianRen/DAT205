#pragma once
#include <vector>
#include <memory>
#include <Eigen/Core>
#include <Eigen/Sparse>

#include "MACGrid.h"

#include "Solver.h"

typedef Eigen::Triplet<double> T;

class Simulator
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Simulator(std::shared_ptr<MACGrid> grids, double &time);
  ~Simulator();

  void update();
  const std::array<double, SIZE> &getDensity()
  {
    return m_grids->density.m_data;
  }

private:
  void setEmitterVelocity();
  void addSource();

  // a cube in the center of the domain
  void setOccupiedVoxels();

  void resetForce();
  void calVorticity();
  void addForce();
  void calPressure();
  void applyPressureTerm();
  void advectVelocity();
  void advectScalar();

  std::shared_ptr<MACGrid> m_grids;
  double &m_time;

  // external force
  double fx[SIZE], fy[SIZE], fz[SIZE];
  // velocity field
  double avg_u[SIZE], avg_v[SIZE], avg_w[SIZE];

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

  
  // solver
  std::vector<T> tripletList;
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> ICCG;

  Eigen::SparseMatrix<double, Eigen::RowMajor> A;
  Eigen::VectorXd b;
  Eigen::VectorXd x;

  CudaSolver m_solver;
};