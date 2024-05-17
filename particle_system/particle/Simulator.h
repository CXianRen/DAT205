#pragma once
#include <vector>
#include <memory>
#include <Eigen/Core>
#include <Eigen/Sparse>

#include "constants.h"

typedef Eigen::Triplet<double> T;

class Simulator
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Simulator(double &time);
  ~Simulator();

  void update();

  std::array<double, SIZE>& getDensity() {
    return m_density;
  }

private:
  void setEmitterVelocity();
  void addSource();

  void resetForce();
  void calVorticity();
  void addForce();
  void calPressure();
  void applyPressureTerm();
  void advectVelocity();
  void advectScalar();

  // std::shared_ptr<MACGrid> m_grids;

  double &m_time;
  
  // density
  std::array<double, SIZE> m_density;
  std::array<double, SIZE> m_density0;
  // temperature
  std::array<double, SIZE> m_temperature;
  std::array<double, SIZE> m_temperature0;
  // pressure
  std::array<double, SIZE> m_pressure;
  // velocity
  std::array<double, (Nx+1) * Ny * Nz> m_u,m_u0;
  std::array<double, Nx * (Ny+1) * Nz> m_v,m_v0;
  std::array<double, Nx * Ny * (Nz+1)> m_w,m_w0;

  std::array<double, SIZE> m_avg_u, m_avg_v, m_avg_w;
  // vorticity
  std::array<double, SIZE> m_omg_x, m_omg_y, m_omg_z;
  // force
  std::array<double, SIZE> m_fx, m_fy, m_fz;


  // solver
  std::vector<T> tripletList;
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> ICCG;

  Eigen::SparseMatrix<double, Eigen::RowMajor> A;
  Eigen::VectorXd b;
  Eigen::VectorXd x;
};