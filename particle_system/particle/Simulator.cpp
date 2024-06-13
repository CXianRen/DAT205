#include <cmath>
#include <random>
#include "Simulator.h"

#include <iostream>
#include "mmath.h"
#include "mperf.h"

Simulator::Simulator(double &time) : m_time(time), b(SIZE), x(SIZE),
                                     CW(MCUDA::CudaWorker(SIZE, Nx, Ny, Nz))
{
    CW.init();

    static auto L = build_3d_laplace<double>(Nx, Ny, Nz);
    m_e_solver.compute(L);
    m_solver.compute(L);

    // initial environment temperature
    FOR_EACH_CELL
    {
        // temperature[ACC3D(i, j, k, Ny, Nx)] = (j / (float)Ny) * T_AMP + T_AMBIENT;
        // temperature[ACC3D(i, j, k, Ny, Nx)] = (j / (float)Ny) * T_AMP + dist(engine) + T_AMBIENT;
        // temperature[ACC3D(i, j, k, Ny, Nx)] = dist(engine);
        temperature[ACC3D(i, j, k, Ny, Nx)] = T_AMBIENT;
    }

    addSource();
    setEmitterVelocity();
}

Simulator::~Simulator()
{
}

void Simulator::update()
{
    if (m_time > FINISH_TIME)
    {
        return;
    }

    clear_measurement();

    T_START("update total")

    T_START("calculateExternalForce")
    calculateExternalForce();
    T_END

    T_START("gpu calculateVorticity")
    CW.setforceField(fx, fy, fz);
    CW.setVelocityField(u, v, w);
    CW.calculateVorticity();
    T_END

    // T_START
    // calculateVorticity();
    // T_END("calculateVorticity")

    T_START("\tgpu applyExternalForce")
    CW.applyExternalForce();
    CW.getVelocityField(u, v, w);
    T_END

    // T_START
    // CW.getforceField(fx, fy, fz);
    // applyExternalForce();
    // T_END("applyExternalForce")

    T_START("gpu advectVelocity")
    CW.advectVelocityField();
    CW.getVelocityField(u, v, w);
    CW.getPreviosVelocityField(
        u0, v0, w0);
    T_END

    T_START("gpu calculatePressure")
    calculatePressure();
    T_END

    T_START("applyPressure")
    applyPressure();
    T_END

    T_START("gpu advectScalarField")

    T_START("\tupdate density and temperature to gpu")
    CW.setVelocityField(u, v, w);

    CW.setDensityField(density);
    CW.setPreviosDensityField(density0);

    CW.setTemperatureField(temperature);
    CW.setPreviosTemperatureField(temperature0);
    T_END

    CW.advectScalarField();

    CW.getDensityField(density);
    CW.getPreviosDensityField(density0);

    CW.getTemperatureField(temperature);
    CW.getPreviosTemperatureField(temperature0);

    T_END
    // T_START
    // advectScalarField();
    // T_END("advectScalarField")

    T_START("fixOccupiedVoxels")
    fixOccupiedVoxels();
    T_END

    T_START("gpu genTransparencyMap")
    genTransparencyMap();
    T_END

    T_END

    if (m_time < EMIT_DURATION)
    {
        addSource();
        setEmitterVelocity();
    }

    m_time += DT;
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
                    // v[ACC3D_Y(i, j, k, Ny, Nx)] = INIT_VELOCITY;
                    // v0[ACC3D_Y(i, j, k, Ny, Nx)] = v[ACC3D_Y(i, j, k, Ny, Nx)];
                    // random velocity
                    v[ACC3D_Y(i, j, k, Ny, Nx)] = INIT_VELOCITY * (rand() % 100) / 100.0;
                    v0[ACC3D_Y(i, j, k, Ny, Nx)] = v[ACC3D_Y(i, j, k, Ny, Nx)];

                    // random velocity for x and z (-0.5, 0.5) * INIT_VELOCITY
                    // u(i, j, k) = (rand() % 100) / 100.0 - 0.5 * INIT_VELOCITY;
                    // u0(i, j, k) = u(i, j, k);

                    // w[ACC3D_Y(i, j, k, Ny, Nx)] = (rand() % 100) / 100.0 - 0.5 * INIT_VELOCITY;
                    // w0[ACC3D_Y(i, j, k, Ny, Nx)] = w[ACC3D_Y(i, j, k, Ny, Nx)];
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
                    v[ACC3D_Y(i, j, k, Ny, Nx)] = -INIT_VELOCITY;
                    v0[ACC3D_Y(i, j, k, Ny, Nx)] = v[ACC3D_Y(i, j, k, Ny, Nx)];
                }
            }
        }
        break;
    }
    }
}

void Simulator::calculateExternalForce()
{
    FOR_EACH_CELL
    {
        fx[ACC3D(i, j, k, Ny, Nx)] = 0.0;
        fy[ACC3D(i, j, k, Ny, Nx)] =
            -ALPHA * density[ACC3D(i, j, k, Ny, Nx)] +
            BETA * (temperature[ACC3D(i, j, k, Ny, Nx)] - T_AMBIENT);
        fz[ACC3D(i, j, k, Ny, Nx)] = 0.0;
    }
}

void Simulator::calculateVorticity()
{

    FOR_EACH_CELL
    {

        calculateAverageVelocity<double>(
            i, j, k,
            Nx, Ny, Nz,
            u, v, w,
            avg_u, avg_v, avg_w);
    }

    FOR_EACH_CELL
    {
        calculateGradient(
            i, j, k,
            Nx, Ny, Nz,
            avg_u, avg_v, avg_w,
            omg_x, omg_y, omg_z);
    }

    FOR_EACH_CELL
    {
        calculateVorticityBody(
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

void Simulator::calculatePressure()
{
    b.setZero();
    // x.setZero();

    double coeff = 1.0;

    T_START("\tBuild b")
    FOR_EACH_CELL
    {
        double F[6] = {
            static_cast<double>(k > 0),
            static_cast<double>(j > 0),
            static_cast<double>(i > 0),
            static_cast<double>(i < Nx - 1),
            static_cast<double>(j < Ny - 1),
            static_cast<double>(k < Nz - 1)};

        static double D[6] = {-1.0, -1.0, -1.0, 1.0, 1.0, 1.0};
        static double U[6];

        U[0] = (double)(w[ACC3D_Z(i, j, k, Ny, Nx)]);
        U[1] = (double)(v[ACC3D_Y(i, j, k, Ny, Nx)]);
        U[2] = (double)(u[ACC3D_X(i, j, k, Ny, Nx)]);
        U[3] = (double)(u[ACC3D_X(i + 1, j, k, Ny, Nx)]);
        U[4] = (double)(v[ACC3D_Y(i, j + 1, k, Ny, Nx)]);
        U[5] = (double)(w[ACC3D_Z(i, j, k + 1, Ny, Nx)]);

        for (int n = 0; n < 6; ++n)
        {
            b(ACC3D(i, j, k, Ny, Nx)) += D[n] * F[n] * U[n];
        }
        b(ACC3D(i, j, k, Ny, Nx)) *= coeff;
    }

    T_END

    T_START("\tSolve")
    m_solver.solve(x, b);
    T_END

    static double err = 0.0;
    m_solver.getError(err);
    // printf("solver error: %f\n", err);
    static int it = 0;
    m_solver.getIterations(it);
    // printf("solver iterations: %d\n", it);

    // printf("#iterations:     %d \n", static_cast<int>(ICCG.iterations()));
    // printf("estimated error: %e \n", ICCG.error());
    static float t_coeff = (VOXEL_SIZE / DT);
    T_START("\tUpdate pressure")
    FOR_EACH_CELL
    {
        // pressure(i, j, k) = x( ACC3D(i, j, k, Ny, Nx)) * t_coeff;
        pressure[ACC3D(i, j, k, Ny, Nx)] = x(ACC3D(i, j, k, Ny, Nx)) * t_coeff;
    }
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
    std::copy(u, u + ((Nx+1) * Ny * Nz), u0);
    std::copy(v, v + (Nx * (Ny+1) * Nz), v0);
    std::copy(w, w + (Nx * Ny * (Nz+1)), w0);
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
    std::copy(u, u + ((Nx+1) * Ny * Nz), u0);
    std::copy(v, v + (Nx * (Ny+1) * Nz), v0);
    std::copy(w, w + (Nx * Ny * (Nz+1)), w0);

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

void Simulator::fixOccupiedVoxels()
{
    FOR_EACH_CELL
    {
        if (m_occupied_voxels[ACC3D(i, j, k, Ny, Nx)])
        {
            u[ACC3D_X(i, j, k, Ny, Nx)] = 0.0;
            v[ACC3D_Y(i, j, k, Ny, Nx)] = 0.0;
            w[ACC3D_Y(i, j, k, Ny, Nx)] = 0.0;
            temperature[ACC3D(i, j, k, Ny, Nx)] = T_AMBIENT;
            density[ACC3D(i, j, k, Ny, Nx)] = 0.0;
        }
    }
}

void Simulator::genTransparencyMap()
{
    CW.genTransparencyMap(
        light_x, light_y, light_z,
        module_scale_factor, factor);
    CW.getTransparencyMap(transparency);
}

std::array<double, SIZE> &
generateSphereDensity()
{
    static std::array<double, SIZE> density;
    density.fill(0.0);
    FOR_EACH_CELL
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
    FOR_EACH_CELL
    {
        // a cube in the center, side length is 1/3 Nx
        if (abs(i - Nx / 2) < Nx / 3 && abs(j - Ny / 2) < 5 && abs(k - Nz / 2) < 5)
        {
            density[ACC3D(i, j, k, Ny, Nx)] = 0.5;
        }
        //  density[ ACC3D(i, j, k, Ny, Nx)] = 0.5;
    }
    return density;
}