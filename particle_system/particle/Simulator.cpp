#include <cmath>
#include <random>
#include "Simulator.h"

#include <iostream>
#include "mmath.h"
#include "mperf.h"

Simulator::Simulator(double &time) : m_time(time), b(SIZE), x(SIZE),
                                     CW(MCUDA::CudaSimulator(SIZE, Nx, Ny, Nz))
{
    CW.init();

    static auto L = build_3d_laplace<double>(Nx, Ny, Nz);
    // m_e_solver.compute(L);
    // m_solver.compute(L);

    // initial environment temperature
    FOR_EACH_CELL
    {
        temperature[ACC3D(i, j, k, Ny, Nx)] = T_AMBIENT;
    }
}

void Simulator::update()
{
    if (m_time > FINISH_TIME)
    {
        return;
    }

    if (m_time < EMIT_DURATION)
    {
        addSource();
        setEmitterVelocity();
    }

    m_time += DT;

    clear_measurement();

    T_START("update total")

    T_START("update data to gpu")
    CW.setDensityField(density);
    CW.setPreviosDensityField(density0);
    CW.setTemperatureField(temperature);
    CW.setPreviosTemperatureField(temperature0);
    CW.setVelocityField(u, v, w);
    CW.setPreviosVelocityField(u0, v0, w0);
    T_END

    T_START("gpu computeExternalForce")
    CW.computeExternalForce();
    T_END

    T_START("gpu computeVorticity")
    CW.computeVorticity();
    T_END

    T_START("gpu applyExternalForce")
    CW.applyExternalForce();
    T_END

    T_START("gpu advectVelocity")
    CW.advectVelocityField();
    T_END

    T_START("gpu computePressure")
    CW.computePressure();
    T_END

    T_START("gpu applyPressure")
    CW.applyPressure();
    T_END

    T_START("gpu advectScalarField")
    CW.advectScalarField();
    T_END

    T_START("get data from gpu")
    CW.getDensityField(density);
    CW.getPreviosDensityField(density0);
    CW.getTemperatureField(temperature);
    CW.getPreviosTemperatureField(temperature0);
    CW.getVelocityField(u, v, w);
    CW.getPreviosVelocityField(u0, v0, w0);
    T_END

    T_START("applyOccupiedVoxels")
    applyOccupiedVoxels();
    T_END

    T_END
}

/* private */
void Simulator::addSource()
{

    std::random_device rnd;
    std::mt19937 engine(rnd());
    std::uniform_real_distribution<double> dist(800, 1000);

    switch (EMITTER_POS)
    {
    case E_TOP:
    {

        for (int k = (Nz - SOURCE_SIZE_Z) / 2; k < (Nz + SOURCE_SIZE_Z) / 2; ++k)
        {
            for (int j = SOURCE_Y_MERGIN; j < SOURCE_Y_MERGIN + SOURCE_SIZE_Y; ++j)
            {
                for (int i = (Nx - SOURCE_SIZE_X) / 2; i < (Nx + SOURCE_SIZE_X) / 2; ++i)
                // for (int i = 20; i < 20 + SOURCE_SIZE_X; ++i)
                {
                    density[ACC3D(i, j, k, Ny, Nx)] = INIT_DENSITY;
                    temperature[ACC3D(i, j, k, Ny, Nx)] = dist(engine);
                }
            }
        }
        break;
    }

    case E_BOTTOM:
    {
        // (32-8) / 2 = 12, (32+8) / 2 = 20
        for (int k = (Nz - SOURCE_SIZE_Z) / 2; k < (Nz + SOURCE_SIZE_Z) / 2; ++k)
        {
            // 64-3-3 = 58, 64-3 = 61
            for (int j = Ny - SOURCE_Y_MERGIN - SOURCE_SIZE_Y; j < Ny - SOURCE_Y_MERGIN; ++j)
            {
                // (32-8) / 2 = 12, (32+8) / 2 = 20
                for (int i = (Nx - SOURCE_SIZE_X) / 2; i < (Nx + SOURCE_SIZE_X) / 2; ++i)
                {
                    density[ACC3D(i, j, k, Ny, Nx)] = INIT_DENSITY;
                    temperature[ACC3D(i, j, k, Ny, Nx)] = dist(engine);
                }
            }
        }
        break;
    }
    }
}

void Simulator::setEmitterVelocity()
{
    switch (EMITTER_POS)
    {
    case E_TOP:
    {

        for (int k = (Nz - SOURCE_SIZE_Z) / 2; k < (Nz + SOURCE_SIZE_Z) / 2; ++k)
        {
            for (int j = SOURCE_Y_MERGIN; j < SOURCE_Y_MERGIN + SOURCE_SIZE_Y; ++j)
            {
                for (int i = (Nx - SOURCE_SIZE_X) / 2; i < (Nx + SOURCE_SIZE_X) / 2; ++i)
                {
                    // v[ACC3D(i, j, k, Ny, Nx)] = INIT_VELOCITY;
                    // v0[ACC3D(i, j, k, Ny, Nx)] = v[ACC3D(i, j, k, Ny, Nx)];
                    // random velocity
                    // v[ACC3D(i, j, k, Ny, Nx)] = INIT_VELOCITY * (rand() % 100) / 100.0;
                    // v0[ACC3D(i, j, k, Ny, Nx)] = v[ACC3D(i, j, k, Ny, Nx)];

                    // random velocity for x and z (-0.5, 0.5) * INIT_VELOCITY
                    // u(i, j, k) = (rand() % 100) / 100.0 - 0.5 * INIT_VELOCITY;
                    // u0(i, j, k) = u(i, j, k);

                    // w[ACC3D(i, j, k, Ny, Nx)] = (rand() % 100) / 100.0 - 0.5 * INIT_VELOCITY;
                    // w0[ACC3D(i, j, k, Ny, Nx)] = w[ACC3D(i, j, k, Ny, Nx)];
                }
            }
        }
        break;
    }

    case E_BOTTOM:
    {

        for (int k = (Nz - SOURCE_SIZE_Z) / 2; k < (Nz + SOURCE_SIZE_Z) / 2; ++k)
        {
            for (int j = Ny - SOURCE_Y_MERGIN - SOURCE_SIZE_Y; j < Ny - SOURCE_Y_MERGIN + 1; ++j)
            {
                for (int i = (Nx - SOURCE_SIZE_X) / 2; i < (Nx + SOURCE_SIZE_X) / 2; ++i)
                {
                    v[ACC3D(i, j, k, Ny, Nx)] = -INIT_VELOCITY;
                    v0[ACC3D(i, j, k, Ny, Nx)] = v[ACC3D(i, j, k, Ny, Nx)];
                }
            }
        }
        break;
    }
    }
}

void Simulator::computeExternalForce()
{
    FOR_EACH_CELL
    {
        computeBuyancyForceBody<double>(
            i, j, k,
            Nx, Ny, Nz,
            density, temperature,
            fx, fy, fz);
    }
}

void Simulator::computeVorticity()
{

    FOR_EACH_CELL
    {

        computeAverageVelocity<double>(
            i, j, k,
            Nx, Ny, Nz,
            u, v, w,
            avg_u, avg_v, avg_w);
    }

    FOR_EACH_CELL
    {
        computeGradient(
            i, j, k,
            Nx, Ny, Nz,
            avg_u, avg_v, avg_w,
            omg_x, omg_y, omg_z);
    }

    FOR_EACH_CELL
    {
        computeVorticityBody(
            i, j, k,
            Nx, Ny, Nz,
            omg_x, omg_y, omg_z,
            fx, fy, fz);
    }
}

void Simulator::applyExternalForce()
{
    FOR_EACH_CELL
    {
        applyExternalForceBody<double>(
            i, j, k,
            Nx, Ny, Nz,
            fx, fy, fz,
            u, v, w);
    }
}

void Simulator::computePressure()
{
    T_START("\tBuild b")
    FOR_EACH_CELL
    {
        b[ACC3D(i, j, k, Ny, Nx)] = 0.0;
        buildRhsBody<double>(
            i, j, k,
            Nx, Ny, Nz,
            u, v, w,
            b.data());
    }
    T_END

    T_START("\tSolve")
    m_solver.solve(pressure, b.data());
    T_END
}

void Simulator::applyPressure()
{

    FOR_EACH_CELL
    {
        applyPressureBody<double>(
            i, j, k,
            Nx, Ny, Nz,
            pressure,
            u, v, w);
    }
}

void Simulator::advectVelocity()
{
    T_START("\tcopy data")
    std::copy(u, u + ((Nx + 1) * Ny * Nz), u0);
    std::copy(v, v + (Nx * (Ny + 1) * Nz), v0);
    std::copy(w, w + (Nx * Ny * (Nz + 1)), w0);
    T_END

    FOR_EACH_CELL
    {
        advectVelocityBody<double>(
            u, v, w,
            u0, v0, w0,
            i, j, k,
            Nx, Ny, Nz);
    }
}

void Simulator::advectScalarField()
{
    T_START("\tcopy data")
    std::copy(u, u + ((Nx + 1) * Ny * Nz), u0);
    std::copy(v, v + (Nx * (Ny + 1) * Nz), v0);
    std::copy(w, w + (Nx * Ny * (Nz + 1)), w0);

    std::copy(density, density + SIZE, density0);
    std::copy(temperature, temperature + SIZE, temperature0);
    T_END

    FOR_EACH_CELL
    {
        advectScalarBody<double>(
            i, j, k,
            Nx, Ny, Nz,
            density, density0,
            u0, v0, w0);

        advectScalarBody<double>(
            i, j, k,
            Nx, Ny, Nz,
            temperature, temperature0,
            u0, v0, w0);
    }
}

void Simulator::applyOccupiedVoxels()
{
    FOR_EACH_CELL
    {
        applyOccupiedVoxelsBody<double>(
            i, j, k,
            Nx, Ny, Nz,
            m_occupied_voxels.data(),
            density,
            u, v, w,
            temperature);
    }
}