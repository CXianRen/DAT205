#include "Vexelization.h"

#include <cmath>

#include "mmath.h"

std::array<bool, SIZE>
generate_vexelized_sphere(int radius)
{
    std::array<bool, SIZE> sphere;

    int x0 = Nx / 2 - radius;
    int x1 = Nx / 2 + radius;
    int y0 = Ny / 2 - radius;
    int y1 = Ny / 2 + radius;
    int z0 = Nz / 2 - radius;
    int z1 = Nz / 2 + radius;

    FOR_EACH_CELL
    {
        if (
            (i - Nx / 2) * (i - Nx / 2) +
                (j - Ny / 2) * (j - Ny / 2) +
                (k - Nz / 2) * (k - Nz / 2) <=
            radius * radius)
        {
            sphere[POS(i, j, k)] = true;
        }
    }
    return sphere;
}