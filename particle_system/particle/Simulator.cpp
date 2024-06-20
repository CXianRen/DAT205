#include <cmath>
#include <random>
#include "Simulator.h"

#include <iostream>
#include "mmath.h"
#include "mperf.h"

Simulator::Simulator(double &time)
    : m_time(time),
      CudaSim(SIZE, Nx, Ny, Nz),
      CPUSim(SIZE, Nx, Ny, Nz)
{
    CPUSim.init();
    CudaSim.init();

    // initial environment temperature
    FOR_EACH_CELL
    {
        temperature[ACC3D(i, j, k, Ny, Nx)] = T_AMBIENT;
    }
}

void Simulator::update_gpu()
{
    T_START("update data to gpu")
    CudaSim.setDensityField(density);
    CudaSim.setPreviosDensityField(density0);
    CudaSim.setTemperatureField(temperature);
    CudaSim.setPreviosTemperatureField(temperature0);
    CudaSim.setVelocityField(u, v, w);
    CudaSim.setPreviosVelocityField(u0, v0, w0);
    T_END

    T_START("gpu computeExternalForce")
    CudaSim.computeExternalForce(alpha_, beta_, envTemp_);
    T_END

    T_START("gpu computeVorticity")
    CudaSim.computeVorticity(vort_eps_);
    T_END

    T_START("gpu applyExternalForce")
    CudaSim.applyExternalForce();
    T_END

    T_START("gpu advectVelocity")
    CudaSim.advectVelocityField();
    T_END

    T_START("gpu computePressure")
    CudaSim.computePressure();
    T_END

    T_START("gpu applyPressure")
    CudaSim.applyPressure();
    T_END

    T_START("gpu advectScalarField")
    CudaSim.advectScalarField();
    T_END

    T_START("get data from gpu")
    CudaSim.getDensityField(density);
    CudaSim.getPreviosDensityField(density0);
    CudaSim.getTemperatureField(temperature);
    CudaSim.getPreviosTemperatureField(temperature0);
    CudaSim.getVelocityField(u, v, w);
    CudaSim.getPreviosVelocityField(u0, v0, w0);
    T_END
}

void Simulator::update_cpu()
{
    CPUSim.setDensityField(density);
    CPUSim.setPreviosDensityField(density0);
    CPUSim.setTemperatureField(temperature);
    CPUSim.setPreviosTemperatureField(temperature0);
    CPUSim.setVelocityField(u, v, w);
    CPUSim.setPreviosVelocityField(u0, v0, w0);
    CPUSim.setforceField(fx, fy, fz);
    CPUSim.setAvgVelocityField(avg_u, avg_v, avg_w);
    CPUSim.setVorticityField(omg_x, omg_y, omg_z);
    CPUSim.setPressureField(pressure);

    T_START("cpu computeExternalForce")
    CPUSim.computeExternalForce(alpha_, beta_, envTemp_);
    T_END

    T_START("cpu computeVorticity")
    CPUSim.computeVorticity(vort_eps_);
    T_END

    T_START("cpu applyExternalForce")
    CPUSim.applyExternalForce();
    T_END

    T_START("cpu advectVelocity")
    CPUSim.advectVelocityField();
    T_END

    T_START("cpu computePressure")
    CPUSim.computePressure();
    T_END

    T_START("cpu applyPressure")
    CPUSim.applyPressure();
    T_END

    T_START("cpu advectScalarField")
    CPUSim.advectScalarField();
    T_END
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

    update_gpu();
    // update_cpu();

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
        
        // decay temperature, speed
        temperature[ACC3D(i, j, k, Ny, Nx)] *= decay_factor_;
        u[ACC3D(i, j, k, Ny, Nx)] *= decay_factor_;
        v[ACC3D(i, j, k, Ny, Nx)] *= decay_factor_;
        w[ACC3D(i, j, k, Ny, Nx)] *= decay_factor_;
        density[ACC3D(i, j, k, Ny, Nx)] *= decay_factor_;

    }
}