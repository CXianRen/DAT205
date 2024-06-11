#pragma once
#include <array>
#include <algorithm>
#include <cassert>
#include <cmath>
#include "constants.h"

template <int X, int Y, int Z>
class MData
{
public:
  MData() : m_data(), maxNx(X - 1), maxNy(Y - 1), maxNz(Z - 1) {}

  ~MData(){};

  double &operator()(int i, int j, int k)
  {
    assert((i >= 0 || i < X) || (j >= 0 || j < Y) || (k >= 0 || k < Z));
    return m_data[i + X * j + X * Y * k];
  }

  double *begin()
  {
    return m_data.begin();
  }

  double *end()
  {
    return m_data.end();
  }

public:
  const int maxNx;
  const int maxNy;
  const int maxNz;

  std::array<double, X * Y * Z> m_data;
};

using MDataScalar = MData<Nx, Ny, Nz>;
using MDataX = MData<Nx + 1, Ny, Nz>;
using MDataY = MData<Nx, Ny + 1, Nz>;
using MDataZ = MData<Nx, Ny, Nz + 1>;