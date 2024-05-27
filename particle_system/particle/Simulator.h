#pragma once
#include <vector>
#include <memory>
#include <Eigen/Core>
#include <Eigen/Sparse>

#include "MACGrid.h"

#include "Solver.h"

typedef Eigen::Triplet<float> T;

class Simulator
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Simulator(std::shared_ptr<MACGrid> grids, double &time);
  ~Simulator();

  void update();
  const std::array<double, SIZE>&  getDensity(){
    return m_grids->density.m_data;
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

  std::shared_ptr<MACGrid> m_grids;
  double &m_time;

  // solver
  std::vector<T> tripletList;
  Eigen::ConjugateGradient<Eigen::SparseMatrix<float>, Eigen::Lower | Eigen::Upper> ICCG;

  Eigen::SparseMatrix<float, Eigen::RowMajor> A;
  Eigen::VectorXf b;
  Eigen::VectorXf x;

  CudaSolver m_solver;
};