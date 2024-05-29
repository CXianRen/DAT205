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
    return density.m_data;
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
    return u0.interp(pos - 0.5 * Vec3(0.0, VOXEL_SIZE, VOXEL_SIZE));
  }
  double getVelocityY(const Vec3 &pos)
  {
    return v0.interp(pos - 0.5 * Vec3(VOXEL_SIZE, 0.0, VOXEL_SIZE));
  }
  double getVelocityZ(const Vec3 &pos)
  {
    return w0.interp(pos - 0.5 * Vec3(VOXEL_SIZE, VOXEL_SIZE, 0.0));
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