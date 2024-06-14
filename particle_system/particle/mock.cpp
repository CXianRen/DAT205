#include "mock.h"
#include "math.h"

namespace MOCK
{
    std::array<double, SIZE> &
    generateSphereDensity()
    {
        static std::array<double, SIZE> density;
        density.fill(0.0);
        for (int k = 0; k < Nz; ++k)
            for (int j = 0; j < Ny; ++j)
                for (int i = 0; i < Nx; ++i)
                {
                    // a ball in the center, radius is 1/3 Nx
                    if (pow(i - Nx / 2, 2) + pow(j - Ny / 2, 2) + pow(k - Nz / 2, 2) < pow(Nx / 4, 2))
                    {
                        density[ACC3D(i, j, k, Ny, Nx)] = 0.5;
                    }
                }
        return density;
    }

    std::array<double, SIZE> &
    generateCubeDensity()
    {
        static std::array<double, SIZE> density;
        density.fill(0.0);
        for (int k = 0; k < Nz; ++k)
            for (int j = 0; j < Ny; ++j)
                for (int i = 0; i < Nx; ++i)
                {
                    // a cube in the center, side length is 1/3 Nx
                    if (abs(i - Nx / 2) < Nx / 3 && abs(j - Ny / 2) < 5 && abs(k - Nz / 2) < 5)
                    {
                        density[ACC3D(i, j, k, Ny, Nx)] = 0.5;
                    }
                }
        return density;
    }
}