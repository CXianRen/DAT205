#include <cmath>
#include <random>
#include "Simulator.h"

#include <iostream>

Simulator::Simulator(std::shared_ptr<MACGrid> grids, double &time) : m_grids(grids), m_time(time), A(SIZE, SIZE), b(SIZE), x(SIZE)
{
    // nnz size is estimated by 7*SIZE because there are 7 nnz elements in a row.(center and neighbor 6)
    tripletList.reserve(7 * SIZE);
    ICCG.setTolerance(1e-6);

    /*set temperature */
    std::random_device rnd;
    std::mt19937 engine(rnd());
    std::uniform_real_distribution<double> dist(800, 1000);

    OPENMP_FOR_COLLAPSE
    FOR_EACH_CELL
    {
        // m_grids->temperature(i, j, k) = (j / (float)Ny) * T_AMP + T_AMBIENT;
        // m_grids->temperature(i, j, k) = (j / (float)Ny) * T_AMP + dist(engine) + T_AMBIENT;
        m_grids->temperature(i, j, k) = dist(engine);
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

    resetForce();
    calVorticity();
    addForce();
    advectVelocity();
    calPressure();
    applyPressureTerm();

    advectScalar();
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
        OPENMP_FOR_COLLAPSE
        for (int k = (Nz - SOURCE_SIZE_Z) / 2; k < (Nz + SOURCE_SIZE_Z) / 2; ++k)
        {
            for (int j = SOURCE_Y_MERGIN; j < SOURCE_Y_MERGIN + SOURCE_SIZE_Y; ++j)
            {
                for (int i = (Nx - SOURCE_SIZE_X) / 2; i < (Nx + SOURCE_SIZE_X) / 2; ++i)
                {
                    m_grids->density(i, j, k) = INIT_DENSITY;
                }
            }
        }
        break;
    }

    case E_BOTTOM:
    {
        OPENMP_FOR_COLLAPSE
        // (32-8) / 2 = 12, (32+8) / 2 = 20
        for (int k = (Nz - SOURCE_SIZE_Z) / 2; k < (Nz + SOURCE_SIZE_Z) / 2; ++k)
        {
            // 64-3-3 = 58, 64-3 = 61
            for (int j = Ny - SOURCE_Y_MERGIN - SOURCE_SIZE_Y; j < Ny - SOURCE_Y_MERGIN; ++j)
            {
                // (32-8) / 2 = 12, (32+8) / 2 = 20
                for (int i = (Nx - SOURCE_SIZE_X) / 2; i < (Nx + SOURCE_SIZE_X) / 2; ++i)
                {
                    m_grids->density(i, j, k) = INIT_DENSITY;
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
        OPENMP_FOR_COLLAPSE
        for (int k = (Nz - SOURCE_SIZE_Z) / 2; k < (Nz + SOURCE_SIZE_Z) / 2; ++k)
        {
            for (int j = SOURCE_Y_MERGIN; j < SOURCE_Y_MERGIN + SOURCE_SIZE_Y; ++j)
            {
                for (int i = (Nx - SOURCE_SIZE_X) / 2; i < (Nx + SOURCE_SIZE_X) / 2; ++i)
                {
                    // m_grids->v(i, j, k) = INIT_VELOCITY;
                    // m_grids->v0(i, j, k) = m_grids->v(i, j, k);
                    // random velocity
                    m_grids->v(i, j, k) = INIT_VELOCITY * (rand() % 100) / 100.0;
                    m_grids->v0(i, j, k) = m_grids->v(i, j, k);

                    // random velocity for x and z (-0.5, 0.5) * INIT_VELOCITY
                    // m_grids->u(i, j, k) = (rand() % 100) / 100.0 - 0.5 * INIT_VELOCITY;
                    // m_grids->u0(i, j, k) = m_grids->u(i, j, k);

                    // m_grids->w(i, j, k) = (rand() % 100) / 100.0 - 0.5 * INIT_VELOCITY;
                    // m_grids->w0(i, j, k) = m_grids->w(i, j, k);
                }
            }
        }
        break;
    }

    case E_BOTTOM:
    {
        OPENMP_FOR_COLLAPSE
        for (int k = (Nz - SOURCE_SIZE_Z) / 2; k < (Nz + SOURCE_SIZE_Z) / 2; ++k)
        {
            for (int j = Ny - SOURCE_Y_MERGIN - SOURCE_SIZE_Y; j < Ny - SOURCE_Y_MERGIN + 1; ++j)
            {
                for (int i = (Nx - SOURCE_SIZE_X) / 2; i < (Nx + SOURCE_SIZE_X) / 2; ++i)
                {
                    m_grids->v(i, j, k) = -INIT_VELOCITY;
                    m_grids->v0(i, j, k) = m_grids->v(i, j, k);
                }
            }
        }
        break;
    }
    }
}

void Simulator::resetForce()
{
    OPENMP_FOR_COLLAPSE
    FOR_EACH_CELL
    {
        fx[POS(i, j, k)] = 0.0;
        fy[POS(i, j, k)] =
            -ALPHA * m_grids->density(i, j, k) +
            BETA * (m_grids->temperature(i, j, k) - T_AMBIENT);
        fz[POS(i, j, k)] = 0.0;
    }
}

void Simulator::calVorticity()
{
    OPENMP_FOR_COLLAPSE
    FOR_EACH_CELL
    {
        m_grids->avg_u[POS(i, j, k)] = (m_grids->u(i, j, k) + m_grids->u(i + 1, j, k)) * 0.5;
        m_grids->avg_v[POS(i, j, k)] = (m_grids->v(i, j, k) + m_grids->v(i, j + 1, k)) * 0.5;
        m_grids->avg_w[POS(i, j, k)] = (m_grids->w(i, j, k) + m_grids->w(i, j, k + 1)) * 0.5;
    }

    OPENMP_FOR_COLLAPSE
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

        m_grids->omg_x[POS(i, j, k)] = (m_grids->avg_w[POS(i, j + 1, k)] - m_grids->avg_w[POS(i, j - 1, k)] - m_grids->avg_v[POS(i, j, k + 1)] + m_grids->avg_v[POS(i, j, k - 1)]) * 0.5 / VOXEL_SIZE;
        m_grids->omg_y[POS(i, j, k)] = (m_grids->avg_u[POS(i, j, k + 1)] - m_grids->avg_u[POS(i, j, k - 1)] - m_grids->avg_w[POS(i + 1, j, k)] + m_grids->avg_w[POS(i - 1, j, k)]) * 0.5 / VOXEL_SIZE;
        m_grids->omg_z[POS(i, j, k)] = (m_grids->avg_v[POS(i + 1, j, k)] - m_grids->avg_v[POS(i - 1, j, k)] - m_grids->avg_u[POS(i, j + 1, k)] + m_grids->avg_u[POS(i, j - 1, k)]) * 0.5 / VOXEL_SIZE;
    }

    OPENMP_FOR_COLLAPSE
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
        p = Vec3(m_grids->omg_x[POS(i + 1, j, k)], m_grids->omg_y[POS(i + 1, j, k)], m_grids->omg_z[POS(i + 1, j, k)]).norm();
        q = Vec3(m_grids->omg_x[POS(i - 1, j, k)], m_grids->omg_y[POS(i - 1, j, k)], m_grids->omg_z[POS(i - 1, j, k)]).norm();
        double grad1 = (p - q) * 0.5 / VOXEL_SIZE;

        p = Vec3(m_grids->omg_x[POS(i, j + 1, k)], m_grids->omg_y[POS(i, j + 1, k)], m_grids->omg_z[POS(i, j + 1, k)]).norm();
        q = Vec3(m_grids->omg_x[POS(i, j - 1, k)], m_grids->omg_y[POS(i, j - 1, k)], m_grids->omg_z[POS(i, j - 1, k)]).norm();
        double grad2 = (p - q) * 0.5 / VOXEL_SIZE;

        p = Vec3(m_grids->omg_x[POS(i, j, k + 1)], m_grids->omg_y[POS(i, j, k + 1)], m_grids->omg_z[POS(i, j, k + 1)]).norm();
        q = Vec3(m_grids->omg_x[POS(i, j, k - 1)], m_grids->omg_y[POS(i, j, k - 1)], m_grids->omg_z[POS(i, j, k - 1)]).norm();
        double grad3 = (p - q) * 0.5 / VOXEL_SIZE;

        Vec3 gradVort(grad1, grad2, grad3);
        // compute N vector
        Vec3 N_ijk(0, 0, 0);
        double norm = gradVort.norm();
        if (norm != 0)
        {
            N_ijk = gradVort / gradVort.norm();
        }

        Vec3 vorticity = Vec3(m_grids->omg_x[POS(i, j, k)], m_grids->omg_y[POS(i, j, k)], m_grids->omg_z[POS(i, j, k)]);
        Vec3 f = VORT_EPS * VOXEL_SIZE * vorticity.cross(N_ijk);
        m_grids->vort[POS(i, j, k)] = f.norm();
        fx[POS(i, j, k)] += f[0];
        fy[POS(i, j, k)] += f[1];
        fz[POS(i, j, k)] += f[2];
    }
}

void Simulator::addForce()
{
    OPENMP_FOR_COLLAPSE
    FOR_EACH_CELL
    {
        if (i < Nx - 1)
        {
            m_grids->u(i + 1, j, k) += DT * (fx[POS(i, j, k)] + fx[POS(i + 1, j, k)]) * 0.5;
        }
        if (j < Ny - 1)
        {
            m_grids->v(i, j + 1, k) += DT * (fy[POS(i, j, k)] + fx[POS(i, j + 1, k)]) * 0.5;
        }
        if (k < Nz - 1)
        {
            m_grids->w(i, j, k + 1) += DT * (fz[POS(i, j, k)] + fx[POS(i, j, k + 1)]) * 0.5;
        }
    }
}

void Simulator::calPressure()
{
    tripletList.clear();
    A.setZero();
    b.setZero();
    x.setZero();

    // float coeff = VOXEL_SIZE / DT;
    double coeff = 1.0;

#pragma omp parallel for collapse(3) ordered
    FOR_EACH_CELL
    {
        double F[6] = {static_cast<double>(k > 0), static_cast<double>(j > 0), static_cast<double>(i > 0),
                       static_cast<double>(i < Nx - 1), static_cast<double>(j < Ny - 1), static_cast<double>(k < Nz - 1)};
        double D[6] = {-1.0, -1.0, -1.0, 1.0, 1.0, 1.0};
        double U[6];
        U[0] = (double)(m_grids->w(i, j, k));
        U[1] = (double)(m_grids->v(i, j, k));
        U[2] = (double)(m_grids->u(i, j, k));
        U[3] = (double)(m_grids->u(i + 1, j, k));
        U[4] = (double)(m_grids->v(i, j + 1, k));
        U[5] = (double)(m_grids->w(i, j, k + 1));
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
        /* solve sparse lenear system by ICCG */
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

    // convert x to double
    // x = x.cast<double>();
    // asign x to m_grids->pressure

    // Eigen::Map<Eigen::VectorXd>(m_grids->pressure.begin(), SIZE) = x;
    for (int k = 0; k < Nz; ++k)
    {
        for (int j = 0; j < Ny; ++j)
        {
            for (int i = 0; i < Nx; ++i)
            {
                m_grids->pressure(i, j, k) = x(POS(i, j, k)) * (VOXEL_SIZE / DT);
            }
        }
    }
}

void Simulator::applyPressureTerm()
{
    OPENMP_FOR_COLLAPSE
    FOR_EACH_CELL
    {
        // compute gradient of pressure
        if (i < Nx - 1)
        {
            m_grids->u(i + 1, j, k) -= DT * (m_grids->pressure(i + 1, j, k) - m_grids->pressure(i, j, k)) / VOXEL_SIZE;
        }
        if (j < Ny - 1)
        {
            m_grids->v(i, j + 1, k) -= DT * (m_grids->pressure(i, j + 1, k) - m_grids->pressure(i, j, k)) / VOXEL_SIZE;
        }
        if (k < Nz - 1)
        {
            m_grids->w(i, j, k + 1) -= DT * (m_grids->pressure(i, j, k + 1) - m_grids->pressure(i, j, k)) / VOXEL_SIZE;
        }
    }
}

void Simulator::advectVelocity()
{
    std::copy(m_grids->u.begin(), m_grids->u.end(), m_grids->u0.begin());
    std::copy(m_grids->v.begin(), m_grids->v.end(), m_grids->v0.begin());
    std::copy(m_grids->w.begin(), m_grids->w.end(), m_grids->w0.begin());

    OPENMP_FOR_COLLAPSE
    FOR_EACH_FACE_X
    {
        Vec3 pos_u = m_grids->getCenter(i, j, k) - 0.5 * Vec3(VOXEL_SIZE, 0, 0);
        Vec3 vel_u = m_grids->getVelocity(pos_u);
        pos_u -= DT * vel_u;
        m_grids->u(i, j, k) = m_grids->getVelocityX(pos_u);
    }

    OPENMP_FOR_COLLAPSE
    FOR_EACH_FACE_Y
    {
        Vec3 pos_v = m_grids->getCenter(i, j, k) - 0.5 * Vec3(0, VOXEL_SIZE, 0);
        Vec3 vel_v = m_grids->getVelocity(pos_v);
        pos_v -= DT * vel_v;
        m_grids->v(i, j, k) = m_grids->getVelocityY(pos_v);
    }

    OPENMP_FOR_COLLAPSE
    FOR_EACH_FACE_Z
    {
        Vec3 pos_w = m_grids->getCenter(i, j, k) - 0.5 * Vec3(0, 0, VOXEL_SIZE);
        Vec3 vel_w = m_grids->getVelocity(pos_w);
        pos_w -= DT * vel_w;
        m_grids->w(i, j, k) = m_grids->getVelocityZ(pos_w);
    }
}

void Simulator::advectScalar()
{
    std::copy(m_grids->density.begin(), m_grids->density.end(), m_grids->density0.begin());
    std::copy(m_grids->temperature.begin(), m_grids->temperature.end(), m_grids->temperature0.begin());

    OPENMP_FOR_COLLAPSE
    FOR_EACH_CELL
    {
        Vec3 pos_cell = m_grids->getCenter(i, j, k);
        Vec3 vel_cell = m_grids->getVelocity(pos_cell);
        pos_cell -= DT * vel_cell;
        m_grids->density(i, j, k) = m_grids->getDensity(pos_cell);
        m_grids->temperature(i, j, k) = m_grids->getTemperature(pos_cell);
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

    // OPENMP_FOR_COLLAPSE
    // FOR_EACH_CELL
    // {
    //     if (i >= x0 && i <= x1 && j >= y0 && j <= y1 && k >= z0 && k <= z1)
    //     {

    //         // velocity is zero
    //         m_grids->u(i, j, k) = 0.0;
    //         m_grids->v(i, j, k) = 0.0;
    //         m_grids->w(i, j, k) = 0.0;

    //         // density is zero
    //         m_grids->density(i, j, k) = 0.0;

    //         // temperature is environment temperature
    //         m_grids->temperature(i, j, k) = T_AMBIENT;
    //     }
    // }

    // a sphere in the center of the domain, R is 5
    int x0 = Nx / 2 - 8;
    int x1 = Nx / 2 + 8;
    int y0 = Ny / 2 - 8;
    int y1 = Ny / 2 + 8;
    int z0 = Nz / 2 - 8;
    int z1 = Nz / 2 + 8;

    OPENMP_FOR_COLLAPSE
    FOR_EACH_CELL
    {
        if ((i - Nx / 2) * (i - Nx / 2) + (j - Ny / 2) * (j - Ny / 2) + (k - Nz / 2) * (k - Nz / 2) <= 64)
        {
            // velocity is zero
            m_grids->u(i, j, k) = 0.0;
            m_grids->v(i, j, k) = 0.0;
            m_grids->w(i, j, k) = 0.0;

            // density is zero
            m_grids->density(i, j, k) = 0.0;

            // temperature is environment temperature
            m_grids->temperature(i, j, k) = T_AMBIENT;
        }
    }
}