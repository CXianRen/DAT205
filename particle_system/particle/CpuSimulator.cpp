#include "CpuSimulator.h"

#include <omp.h>

namespace CPUSIM
{

    void CpuSimulator::computeExternalForce()
    {
        #pragma omp parallel for
        FOR_EACH_CELL
        {
            computeBuyancyForceBody<double>(
                i, j, k,
                Nx_, Ny_, Nz_,
                density, temperature,
                f_x, f_y, f_z);
        }
    }

    void CpuSimulator::computeVorticity()
    {   
        #pragma omp parallel for
        FOR_EACH_CELL
        {

            computeAverageVelocity<double>(
                i, j, k,
                Nx_, Ny_, Nz_,
                u, v, w,
                avg_u, avg_v, avg_w);
        }

        #pragma omp parallel for
        FOR_EACH_CELL
        {
            computeGradient(
                i, j, k,
                Nx, Ny, Nz,
                avg_u, avg_v, avg_w,
                omg_x, omg_y, omg_z);
        }

        #pragma omp parallel for
        FOR_EACH_CELL
        {
            computeVorticityBody(
                i, j, k,
                Nx_, Ny_, Nz_,
                omg_x, omg_y, omg_z,
                f_x, f_y, f_z);
        }
    }

    void CpuSimulator::applyExternalForce()
    {
        #pragma omp parallel for
        FOR_EACH_CELL
        {
            applyExternalForceBody<double>(
                i, j, k,
                Nx_, Ny_, Nz_,
                f_x, f_y, f_z,
                u, v, w);
        }
    }

    void CpuSimulator::computePressure()
    {
        #pragma omp parallel for
        FOR_EACH_CELL
        {
            b[ACC3D(i, j, k, Ny_, Nx_)] = 0.0;
            buildRhsBody<double>(
                i, j, k,
                Nx_, Ny_, Nz_,
                u, v, w,
                b.data());
        }
        // solver_.solve(x, b);
        solver_.solveWithGuess(x, b);
        
        #pragma omp parallel for
        FOR_EACH_CELL
        {
            pressure[ACC3D(i, j, k, Ny_, Nx_)] = x(ACC3D(i, j, k, Ny_, Nx_));
        }
    }

    void CpuSimulator::applyPressure()
    {   
        #pragma omp parallel for
        FOR_EACH_CELL
        {
            applyPressureBody<double>(
                i, j, k,
                Nx_, Ny_, Nz_,
                pressure,
                u, v, w);
        }
    }

    void CpuSimulator::advectVelocityField()
    {
        std::copy(u, u + SIZE, u_0);
        std::copy(v, v + SIZE, v_0);
        std::copy(w, w + SIZE, w_0);

        #pragma omp parallel for
        FOR_EACH_CELL
        {
            advectVelocityBody<double>(
                u, v, w,
                u_0, v_0, w_0,
                i, j, k,
                Nx_, Ny_, Nz_);
        }
    }

    void CpuSimulator::advectScalarField()
    {
        //@todo should we copy it ? why not use u v w directly?
        std::copy(u, u + SIZE, u_0);
        std::copy(v, v + SIZE, v_0);
        std::copy(w, w + SIZE, w_0);

        std::copy(density, density + SIZE, density_0);
        std::copy(temperature, temperature + SIZE, temperature_0);

        #pragma omp parallel for
        FOR_EACH_CELL
        {
            advectScalarBody<double>(
                i, j, k,
                Nx_, Ny_, Nz_,
                density, density_0,
                u_0, v_0, w_0);

            advectScalarBody<double>(
                i, j, k,
                Nx, Ny, Nz,
                temperature, temperature_0,
                u_0, v_0, w_0);
        }
    }

} // namespace CPUSIM
