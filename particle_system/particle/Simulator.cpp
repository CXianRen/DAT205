#include <cmath>
#include <random>
#include "Simulator.h"

#include <iostream>
#include "mmath.h"

Simulator::Simulator(double &time) : m_time(time), A(SIZE, SIZE), b(SIZE), x(SIZE)
{
    // nnz size is estimated by 7*SIZE because there are 7 nnz elements in a row.(center and neighbor 6)
    tripletList.reserve(7 * SIZE);
    ICCG.setTolerance(1e-6);

    /*set temperature */
    std::random_device rnd;
    std::mt19937 engine(rnd());
    std::uniform_real_distribution<double> dist(800, 1000);

    FOR_EACH_CELL
    {
        // temperature(i, j, k) = (j / (float)Ny) * T_AMP + T_AMBIENT;
        // temperature(i, j, k) = (j / (float)Ny) * T_AMP + dist(engine) + T_AMBIENT;
        temperature(i, j, k) = dist(engine);
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

    std::cout << "m_time:" << m_time << " EEMIT_DURATION:" << EMIT_DURATION << std::endl;

    calculate_external_force();
    calculate_vorticity();
    apply_external_force();
    advect_velocity();
    calculate_pressure();
    apply_pressure();

    advect_scalar_field();
    if (m_time < EMIT_DURATION)
    {

        addSource();
        setEmitterVelocity();
    }
    setOccupiedVoxels();
    m_time += DT;
}

/* private */
void Simulator::addSource()
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
                    density(i, j, k) = INIT_DENSITY;
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
                    density(i, j, k) = INIT_DENSITY;
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
                    // v(i, j, k) = INIT_VELOCITY;
                    // v0(i, j, k) = v(i, j, k);
                    // random velocity
                    v(i, j, k) = INIT_VELOCITY * (rand() % 100) / 100.0;
                    v0(i, j, k) = v(i, j, k);

                    // random velocity for x and z (-0.5, 0.5) * INIT_VELOCITY
                    // u(i, j, k) = (rand() % 100) / 100.0 - 0.5 * INIT_VELOCITY;
                    // u0(i, j, k) = u(i, j, k);

                    // w(i, j, k) = (rand() % 100) / 100.0 - 0.5 * INIT_VELOCITY;
                    // w0(i, j, k) = w(i, j, k);
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
                    v(i, j, k) = -INIT_VELOCITY;
                    v0(i, j, k) = v(i, j, k);
                }
            }
        }
        break;
    }
    }
}

void Simulator::calculate_external_force()
{
    FOR_EACH_CELL
    {
        fx[POS(i, j, k)] = 0.0;
        fy[POS(i, j, k)] =
            -ALPHA * density(i, j, k) +
            BETA * (temperature(i, j, k) - T_AMBIENT);
        fz[POS(i, j, k)] = 0.0;
    }
}

void Simulator::calculate_vorticity()
{

    FOR_EACH_CELL
    {
        avg_u[POS(i, j, k)] = (u(i, j, k) + u(i + 1, j, k)) * 0.5;
        avg_v[POS(i, j, k)] = (v(i, j, k) + v(i, j + 1, k)) * 0.5;
        avg_w[POS(i, j, k)] = (w(i, j, k) + w(i, j, k + 1)) * 0.5;
    }

    FOR_EACH_CELL
    {
        // ignore boundary cells
        if (i == 0 || j == 0 || k == 0)
        {
            continue;
        }
        if (i == Nx - 1 || j == Ny - 1 || k == Nz - 1)
        {
            continue;
        }

        omg_x[POS(i, j, k)] = (avg_w[POS(i, j + 1, k)] - avg_w[POS(i, j - 1, k)] - avg_v[POS(i, j, k + 1)] + avg_v[POS(i, j, k - 1)]) * 0.5 / VOXEL_SIZE;
        omg_y[POS(i, j, k)] = (avg_u[POS(i, j, k + 1)] - avg_u[POS(i, j, k - 1)] - avg_w[POS(i + 1, j, k)] + avg_w[POS(i - 1, j, k)]) * 0.5 / VOXEL_SIZE;
        omg_z[POS(i, j, k)] = (avg_v[POS(i + 1, j, k)] - avg_v[POS(i - 1, j, k)] - avg_u[POS(i, j + 1, k)] + avg_u[POS(i, j - 1, k)]) * 0.5 / VOXEL_SIZE;
    }

    FOR_EACH_CELL
    {
        // ignore boundary cells
        if (i == 0 || j == 0 || k == 0)
        {
            continue;
        }
        if (i == Nx - 1 || j == Ny - 1 || k == Nz - 1)
        {
            continue;
        }
        // compute gradient of vorticity
        double p, q;
        p = Vec3(omg_x[POS(i + 1, j, k)], omg_y[POS(i + 1, j, k)], omg_z[POS(i + 1, j, k)]).norm();
        q = Vec3(omg_x[POS(i - 1, j, k)], omg_y[POS(i - 1, j, k)], omg_z[POS(i - 1, j, k)]).norm();
        double grad1 = (p - q) * 0.5 / VOXEL_SIZE;

        p = Vec3(omg_x[POS(i, j + 1, k)], omg_y[POS(i, j + 1, k)], omg_z[POS(i, j + 1, k)]).norm();
        q = Vec3(omg_x[POS(i, j - 1, k)], omg_y[POS(i, j - 1, k)], omg_z[POS(i, j - 1, k)]).norm();
        double grad2 = (p - q) * 0.5 / VOXEL_SIZE;

        p = Vec3(omg_x[POS(i, j, k + 1)], omg_y[POS(i, j, k + 1)], omg_z[POS(i, j, k + 1)]).norm();
        q = Vec3(omg_x[POS(i, j, k - 1)], omg_y[POS(i, j, k - 1)], omg_z[POS(i, j, k - 1)]).norm();
        double grad3 = (p - q) * 0.5 / VOXEL_SIZE;

        Vec3 gradVort(grad1, grad2, grad3);
        // compute N vector
        Vec3 N_ijk(0, 0, 0);
        double norm = gradVort.norm();
        if (norm != 0)
        {
            N_ijk = gradVort / gradVort.norm();
        }

        Vec3 vorticity = Vec3(omg_x[POS(i, j, k)], omg_y[POS(i, j, k)], omg_z[POS(i, j, k)]);
        Vec3 f = VORT_EPS * VOXEL_SIZE * vorticity.cross(N_ijk);
        vort[POS(i, j, k)] = f.norm();
        fx[POS(i, j, k)] += f[0];
        fy[POS(i, j, k)] += f[1];
        fz[POS(i, j, k)] += f[2];
    }
}

void Simulator::apply_external_force()
{
    FOR_EACH_CELL
    {
        if (i < Nx - 1)
        {
            u(i + 1, j, k) += DT * (fx[POS(i, j, k)] + fx[POS(i + 1, j, k)]) * 0.5;
        }
        if (j < Ny - 1)
        {
            v(i, j + 1, k) += DT * (fy[POS(i, j, k)] + fx[POS(i, j + 1, k)]) * 0.5;
        }
        if (k < Nz - 1)
        {
            w(i, j, k + 1) += DT * (fz[POS(i, j, k)] + fx[POS(i, j, k + 1)]) * 0.5;
        }
    }
}

void Simulator::calculate_pressure()
{
    tripletList.clear();
    A.setZero();
    b.setZero();
    x.setZero();

    double coeff = 1.0;

    FOR_EACH_CELL
    {
        double F[6] = {static_cast<double>(k > 0), static_cast<double>(j > 0), static_cast<double>(i > 0),
                       static_cast<double>(i < Nx - 1), static_cast<double>(j < Ny - 1), static_cast<double>(k < Nz - 1)};
        double D[6] = {-1.0, -1.0, -1.0, 1.0, 1.0, 1.0};
        double U[6];
        U[0] = (double)(w(i, j, k));
        U[1] = (double)(v(i, j, k));
        U[2] = (double)(u(i, j, k));
        U[3] = (double)(u(i + 1, j, k));
        U[4] = (double)(v(i, j + 1, k));
        U[5] = (double)(w(i, j, k + 1));
        double sum_F = 0.0;

        for (int n = 0; n < 6; ++n)
        {
            sum_F += F[n];
            b(POS(i, j, k)) += D[n] * F[n] * U[n];
        }
        b(POS(i, j, k)) *= coeff;

#pragma omp ordered
        {
            if (k > 0)
            {
                tripletList.push_back(T(POS(i, j, k), POS(i, j, k - 1), F[0]));
            }
            if (j > 0)
            {
                tripletList.push_back(T(POS(i, j, k), POS(i, j - 1, k), F[1]));
            }
            if (i > 0)
            {
                tripletList.push_back(T(POS(i, j, k), POS(i - 1, j, k), F[2]));
            }

            tripletList.push_back(T(POS(i, j, k), POS(i, j, k), -sum_F));

            if (i < Nx - 1)
            {
                tripletList.push_back(T(POS(i, j, k), POS(i + 1, j, k), F[3]));
            }
            if (j < Ny - 1)
            {
                tripletList.push_back(T(POS(i, j, k), POS(i, j + 1, k), F[4]));
            }
            if (k < Nz - 1)
            {
                tripletList.push_back(T(POS(i, j, k), POS(i, j, k + 1), F[5]));
            }
        }
    }

    {
        static bool first = true;
        if (first)
        {
            first = false;
            A.setFromTriplets(tripletList.begin(), tripletList.end());
            ICCG.compute(A);
            m_solver.compute(A);
        }
    }
    // if (ICCG.info() == Eigen::Success)
    // {
    //     printf("SUCCESS: Convergence\n");
    // }
    // else
    // {
    //     fprintf(stderr, "FAILED: No Convergence\n");
    // }
    // x = ICCG.solve(b);

    m_solver.solve(x, b);

    static double err = 0.0;
    m_solver.getError(err);
    printf("solver error: %f\n", err);
    static int it = 0;
    m_solver.getIterations(it);
    printf("solver iterations: %d\n", it);

    // printf("#iterations:     %d \n", static_cast<int>(ICCG.iterations()));
    // printf("estimated error: %e \n", ICCG.error());

    FOR_EACH_CELL
    {
        pressure(i, j, k) = x(POS(i, j, k)) * (VOXEL_SIZE / DT);
    }
}

void Simulator::apply_pressure()
{

    FOR_EACH_CELL
    {
        // compute gradient of pressure
        if (i < Nx - 1)
        {
            u(i + 1, j, k) -= DT * (pressure(i + 1, j, k) - pressure(i, j, k)) / VOXEL_SIZE;
        }
        if (j < Ny - 1)
        {
            v(i, j + 1, k) -= DT * (pressure(i, j + 1, k) - pressure(i, j, k)) / VOXEL_SIZE;
        }
        if (k < Nz - 1)
        {
            w(i, j, k + 1) -= DT * (pressure(i, j, k + 1) - pressure(i, j, k)) / VOXEL_SIZE;
        }
    }
}

void Simulator::advect_velocity()
{
    std::copy(u.begin(), u.end(), u0.begin());
    std::copy(v.begin(), v.end(), v0.begin());
    std::copy(w.begin(), w.end(), w0.begin());

    FOR_EACH_FACE_X
    {
        Vec3 pos_u = getCenter(i, j, k) - 0.5 * Vec3(VOXEL_SIZE, 0, 0);
        Vec3 vel_u = getVelocity(pos_u);
        pos_u -= DT * vel_u;
        u(i, j, k) = getVelocityX(pos_u);
    }

    FOR_EACH_FACE_Y
    {
        Vec3 pos_v = getCenter(i, j, k) - 0.5 * Vec3(0, VOXEL_SIZE, 0);
        Vec3 vel_v = getVelocity(pos_v);
        pos_v -= DT * vel_v;
        v(i, j, k) = getVelocityY(pos_v);
    }

    FOR_EACH_FACE_Z
    {
        Vec3 pos_w = getCenter(i, j, k) - 0.5 * Vec3(0, 0, VOXEL_SIZE);
        Vec3 vel_w = getVelocity(pos_w);
        pos_w -= DT * vel_w;
        w(i, j, k) = getVelocityZ(pos_w);
    }
}

void Simulator::advect_scalar_field()
{
    std::copy(density.begin(), density.end(), density0.begin());
    std::copy(temperature.begin(), temperature.end(), temperature0.begin());

    FOR_EACH_CELL
    {
        Vec3 pos_cell = getCenter(i, j, k);
        Vec3 vel_cell = getVelocity(pos_cell);
        pos_cell -= DT * vel_cell;
        density(i, j, k) = getDensity(pos_cell);
        temperature(i, j, k) = getTemperature(pos_cell);
    }
}

void Simulator::setOccupiedVoxels()
{
    // a cube in the center of the domain
    // int x0 = Nx / 2 - 4;
    // int x1 = Nx / 2 + 4;
    // int y0 = Ny / 2 - 4;
    // int y1 = Ny / 2 + 4;
    // int z0 = Nz / 2 - 4;
    // int z1 = Nz / 2 + 4;

    //
    // FOR_EACH_CELL
    // {
    //     if (i >= x0 && i <= x1 && j >= y0 && j <= y1 && k >= z0 && k <= z1)
    //     {

    //         // velocity is zero
    //         u(i, j, k) = 0.0;
    //         v(i, j, k) = 0.0;
    //         w(i, j, k) = 0.0;

    //         // density is zero
    //         density(i, j, k) = 0.0;

    //         // temperature is environment temperature
    //         temperature(i, j, k) = T_AMBIENT;
    //     }
    // }

    // a sphere in the center of the domain, R is 5
    int x0 = Nx / 2 - 8;
    int x1 = Nx / 2 + 8;
    int y0 = Ny / 2 - 8;
    int y1 = Ny / 2 + 8;
    int z0 = Nz / 2 - 8;
    int z1 = Nz / 2 + 8;

    FOR_EACH_CELL
    {
        if (
            (i - Nx / 2) * (i - Nx / 2) +
                (j - Ny / 2) * (j - Ny / 2) +
                (k - Nz / 2) * (k - Nz / 2) <=
            64)
        {
            // velocity is zero
            u(i, j, k) = 0.0;
            v(i, j, k) = 0.0;
            w(i, j, k) = 0.0;

            // density is zero
            density(i, j, k) = 0.0;

            // temperature is environment temperature
            temperature(i, j, k) = T_AMBIENT;
        }
    }
}