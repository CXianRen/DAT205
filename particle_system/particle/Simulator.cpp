#include <iostream>
#include <cmath>
#include <random>

#include <perf.h>

#include "Simulator.h"
#include "common/mmath.h"


Simulator::Simulator(double &time) : m_time(time), A(SIZE, SIZE), b(SIZE), x(SIZE)
{
    // nnz size is estimated by 7*SIZE because there are 7 nnz elements in a row.(center and neighbor 6)
    tripletList.reserve(7 * SIZE);
    ICCG.setTolerance(1e-8);

    /*set temperature */
    std::random_device rnd;
    std::mt19937 engine(rnd());
    std::uniform_real_distribution<double> dist(0, T_AMP);

    FOR_EACH_CELL
    {
        // m_grids->temperature(i, j, k) = (j / (float)Ny) * T_AMP + T_AMBIENT;
        m_temperature[POS(i, j, k)] = (j / (float)Ny) * T_AMP + dist(engine) + T_AMBIENT;
    }

    addSource();
    setEmitterVelocity();
}

Simulator::~Simulator()
{
}

void Simulator::update()
{
    labhelper::perf::Scope s("update");
    {
        labhelper::perf::Scope s("resetForce");
        resetForce();
    }
    {
        labhelper::perf::Scope s("calVorticity");
        calVorticity();
    }

    {
        labhelper::perf::Scope s("addForce");
        addForce();
    }

    {
        labhelper::perf::Scope s("calPressure");
        calPressure();
    }
    {
        labhelper::perf::Scope s("applyPressureTerm");
        applyPressureTerm();
    }
    {
        labhelper::perf::Scope s("advectVelocity");
        advectVelocity();
    }
    {
        labhelper::perf::Scope s("advectScalar");
        advectScalar();
    }
    if (m_time < EMIT_DURATION)
    {
        addSource();
        setEmitterVelocity();
    }
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
                    // m_grids->density(i, j, k) = INIT_DENSITY;
                    m_density[POS(i, j, k)] = INIT_DENSITY;
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
                    // m_grids->density(i, j, k) = INIT_DENSITY;
                    m_density[POS(i, j, k)] = INIT_DENSITY;
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
                    m_v[POS_Y(i, j, k)] = INIT_VELOCITY;
                    m_v0[POS_Y(i, j, k)] = m_v[POS_Y(i, j, k)];
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
                    m_v[POS_Y(i, j, k)] = -INIT_VELOCITY;
                    m_v0[POS_Y(i, j, k)] = m_v[POS_Y(i, j, k)];
                }
            }
        }
        break;
    }
    }
}

void Simulator::resetForce()
{
    FOR_EACH_CELL
    {
        m_fx[POS(i, j, k)] = 0.0;
        m_fy[POS(i, j, k)] = ALPHA * m_density[POS(i, j, k)];
        -BETA *(m_temperature[POS(i, j, k)] - T_AMBIENT);
        m_fz[POS(i, j, k)] = 0.0;
    }
}

void Simulator::calVorticity()
{
    {
        labhelper::perf::Scope s("calculate_3D_Field_average");
        calculate_3D_Field_average<double>(
            std::vector<double *>{m_u.data(), m_v.data(), m_w.data()},
            std::vector<double *>{m_avg_u.data(), m_avg_v.data(), m_avg_w.data()},
            std::vector<int>{Nx, Ny, Nz});
    }
    {
        labhelper::perf::Scope s2("calculate_3D_Field_vorticity");
        calculate_3D_Field_vorticity<double>(
            std::vector<double *>{m_avg_u.data(), m_avg_v.data(), m_avg_w.data()},
            std::vector<double *>{m_omg_x.data(), m_omg_y.data(), m_omg_z.data()},
            std::vector<int>{Nx, Ny, Nz}, VOXEL_SIZE);
    }
    static double vorticity_len[SIZE];
    static double vorticity_len_gradient_x[SIZE];
    static double vorticity_len_gradient_y[SIZE];
    static double vorticity_len_gradient_z[SIZE];

    {
        labhelper::perf::Scope s3("calculate vorticity length");

        for (int k = 1; k < Nz - 1; ++k)
            for (int j = 1; j < Ny - 1; ++j)
                for (int i = 1; i < Nx - 1; ++i)
                {
                    vorticity_len[POS(i, j, k)] =
                        Vec3(
                            m_omg_x[POS(i, j, k)],
                            m_omg_y[POS(i, j, k)],
                            m_omg_z[POS(i, j, k)])
                            .norm();
                }
    }
    {
        labhelper::perf::Scope s4("calculate vorticity length gradient");
        calculate_Scalar_Field_gradient<double>(
            vorticity_len,
            std::vector<double *>{vorticity_len_gradient_x, vorticity_len_gradient_y, vorticity_len_gradient_z},
            std::vector<int>{Nx, Ny, Nz}, VOXEL_SIZE);
    }
    {
        labhelper::perf::Scope s5("calculate vorticity force");
        for (int k = 1; k < Nz - 1; ++k)
            for (int j = 1; j < Ny - 1; ++j)
                for (int i = 1; i < Nx - 1; ++i)
                {
                    // compute N vector
                    Vec3 gradVort(
                        vorticity_len_gradient_x[POS(i, j, k)],
                        vorticity_len_gradient_y[POS(i, j, k)],
                        vorticity_len_gradient_z[POS(i, j, k)]);
                    // compute N vector
                    Vec3 N_ijk(0, 0, 0);
                    double norm = gradVort.norm();
                    if (norm != 0)
                    {
                        N_ijk = gradVort / gradVort.norm();
                    }

                    Vec3 vorticity = Vec3(m_omg_x[POS(i, j, k)], m_omg_y[POS(i, j, k)], m_omg_z[POS(i, j, k)]);
                    Vec3 f = VORT_EPS * VOXEL_SIZE * vorticity.cross(N_ijk);
                    // m_grids->vort[POS(i, j, k)] = f.norm();
                    m_fx[POS(i, j, k)] += f[0];
                    m_fy[POS(i, j, k)] += f[1];
                    m_fz[POS(i, j, k)] += f[2];
                }
    }
}

void Simulator::addForce()
{
    FOR_EACH_CELL
    {
        if (i < Nx - 1)
        {
            m_u[POS_X(i + 1, j, k)] += DT * (m_fx[POS(i, j, k)] + m_fx[POS(i + 1, j, k)]) * 0.5;
        }
        if (j < Ny - 1)
        {
            m_v[POS_Y(i, j + 1, k)] += DT * (m_fy[POS(i, j, k)] + m_fy[POS(i, j + 1, k)]) * 0.5;
        }
        if (k < Nz - 1)
        {
            m_w[POS_Z(i, j, k + 1)] += DT * (m_fz[POS(i, j, k)] + m_fz[POS(i, j, k + 1)]) * 0.5;
        }
    }
}

void Simulator::calPressure()
{
    tripletList.clear();
    A.setZero();
    b.setZero();
    x.setZero();

    double coeff = VOXEL_SIZE / DT;

    FOR_EACH_CELL
    {
        double F[6] = {static_cast<double>(k > 0), static_cast<double>(j > 0), static_cast<double>(i > 0),
                       static_cast<double>(i < Nx - 1), static_cast<double>(j < Ny - 1), static_cast<double>(k < Nz - 1)};
        double D[6] = {-1.0, -1.0, -1.0, 1.0, 1.0, 1.0};
        double U[6];
        U[0] = m_w[POS_Z(i, j, k)];
        U[1] = m_v[POS_Y(i, j, k)];
        U[2] = m_u[POS_X(i, j, k)];
        U[3] = m_u[POS_X(i + 1, j, k)];
        U[4] = m_v[POS_Y(i, j + 1, k)];
        U[5] = m_w[POS_Z(i, j, k + 1)];
        double sum_F = 0.0;

        for (int n = 0; n < 6; ++n)
        {
            sum_F += F[n];
            b(POS(i, j, k)) += D[n] * F[n] * U[n];
        }
        b(POS(i, j, k)) *= coeff;
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

    A.setFromTriplets(tripletList.begin(), tripletList.end());

    /* solve sparse lenear system by ICCG */
    ICCG.compute(A);
    if (ICCG.info() == Eigen::Success)
    {
        printf("SUCCESS: Convergence\n");
    }
    else
    {
        fprintf(stderr, "FAILED: No Convergence\n");
    }
    x = ICCG.solve(b);
    printf("#iterations:     %d \n", static_cast<int>(ICCG.iterations()));
    printf("estimated error: %e \n", ICCG.error());

    // asign x to m_grids->pressure
    Eigen::Map<Eigen::VectorXd>(m_pressure.begin(), SIZE) = x;
}

void Simulator::applyPressureTerm()
{
    FOR_EACH_CELL
    {
        if (i < Nx - 1)
        {
            m_u[POS_X(i + 1, j, k)] -= DT * (m_pressure[POS(i + 1, j, k)] - m_pressure[POS(i, j, k)]) / VOXEL_SIZE;
        }
        if (j < Ny - 1)
        {
            m_v[POS_Y(i, j + 1, k)] -= DT * (m_pressure[POS(i, j + 1, k)] - m_pressure[POS(i, j, k)]) / VOXEL_SIZE;
        }
        if (k < Nz - 1)
        {
            m_w[POS_Z(i, j, k + 1)] -= DT * (m_pressure[POS(i, j, k + 1)] - m_pressure[POS(i, j, k)]) / VOXEL_SIZE;
        }
    }
    std::copy(m_u.begin(), m_u.end(), m_u0.begin());
    std::copy(m_v.begin(), m_v.end(), m_v0.begin());
    std::copy(m_w.begin(), m_w.end(), m_w0.begin());
}

double inline getVelocityX(Vec3 pos, double *u_x, const std::vector<int> &dims)
{
    switch (INTERPOLATION_METHOD)
    {
    case E_LINEAR:
        return linearInterpolation<double>(
            pos - 0.5 * Vec3(0.0, VOXEL_SIZE, VOXEL_SIZE),
            u_x,
            {dims[0] + 1, dims[1], dims[2]},
            VOXEL_SIZE);

    case E_MONOTONIC_CUBIC:
        return monotonicCubicInterpolation<double>(
            pos - 0.5 * Vec3(0.0, VOXEL_SIZE, VOXEL_SIZE),
            u_x,
            {dims[0] + 1, dims[1], dims[2]},
            VOXEL_SIZE);
    }
}

double inline getVelocityY(Vec3 pos, double *u_y, const std::vector<int> &dims)
{
    switch (INTERPOLATION_METHOD)
    {
    case E_LINEAR:
        return linearInterpolation<double>(
            pos - 0.5 * Vec3(VOXEL_SIZE, 0.0, VOXEL_SIZE),
            u_y,
            {dims[0], dims[1] + 1, dims[2]},
            VOXEL_SIZE);

    case E_MONOTONIC_CUBIC:
        return monotonicCubicInterpolation<double>(
            pos - 0.5 * Vec3(VOXEL_SIZE, 0.0, VOXEL_SIZE),
            u_y,
            {dims[0], dims[1] + 1, dims[2]},
            VOXEL_SIZE);
    }
}

double inline getVelocityZ(Vec3 pos, double *u_z, const std::vector<int> &dims)
{
    switch (INTERPOLATION_METHOD)
    {
    case E_LINEAR:
        return linearInterpolation<double>(
            pos - 0.5 * Vec3(VOXEL_SIZE, VOXEL_SIZE, 0.0),
            u_z,
            {dims[0], dims[1], dims[2] + 1},
            VOXEL_SIZE);

    case E_MONOTONIC_CUBIC:
        return monotonicCubicInterpolation<double>(
            pos - 0.5 * Vec3(VOXEL_SIZE, VOXEL_SIZE, 0.0),
            u_z,
            {dims[0], dims[1], dims[2] + 1},
            VOXEL_SIZE);
    }
}

Vec3 inline getVelocity(const Vec3 &pos, const std::vector<double *> src, const std::vector<int> &dims)
{
    double vel_u = getVelocityX(pos, src[0], dims);
    double vel_v = getVelocityY(pos, src[1], dims);
    double vel_w = getVelocityZ(pos, src[2], dims);
    return Vec3(vel_u, vel_v, vel_w);
}

Vec3 inline getCenter(int i, int j, int k)
{
    double half_dx = 0.5 * VOXEL_SIZE;

    double x = half_dx + i * VOXEL_SIZE;
    double y = half_dx + j * VOXEL_SIZE;
    double z = half_dx + k * VOXEL_SIZE;
    return Vec3(x, y, z);
}

void Simulator::advectVelocity()
{
    FOR_EACH_CELL
    {
        Vec3 pos_u = getCenter(i, j, k) - 0.5 * Vec3(VOXEL_SIZE, 0, 0);
        Vec3 vel_u =
            getVelocity(
                pos_u,
                std::vector<double *>{m_u0.data(), m_v0.data(), m_w0.data()},
                std::vector<int>{Nx, Ny, Nz});

        pos_u -= DT * vel_u;
        m_u[POS_X(i, j, k)] =
            getVelocityX(pos_u, m_u0.data(), std::vector<int>{Nx, Ny, Nz});
    }
    FOR_EACH_CELL
    {
        Vec3 pos_v = getCenter(i, j, k) - 0.5 * Vec3(0, VOXEL_SIZE, 0);
        // Vec3 vel_v = m_grids->getVelocity(pos_v);
        Vec3 vel_v =
            getVelocity(
                pos_v,
                std::vector<double *>{m_u0.data(), m_v0.data(), m_w0.data()},
                std::vector<int>{Nx, Ny, Nz});
        pos_v -= DT * vel_v;
        // m_v[POS_Y(i, j, k)] = m_grids->getVelocityY(pos_v);
        m_v[POS_Y(i, j, k)] =
            getVelocityY(pos_v, m_v0.data(), std::vector<int>{Nx, Ny, Nz});
    }
    FOR_EACH_CELL
    {
        Vec3 pos_w = getCenter(i, j, k) - 0.5 * Vec3(0, 0, VOXEL_SIZE);
        // Vec3 vel_w = m_grids->getVelocity(pos_w);
        Vec3 vel_w =
            getVelocity(
                pos_w,
                std::vector<double *>{m_u0.data(), m_v0.data(), m_w0.data()},
                std::vector<int>{Nx, Ny, Nz});
        pos_w -= DT * vel_w;
        // m_w[POS_Z(i, j, k)] = m_grids->getVelocityZ(pos_w);
        m_w[POS_Z(i, j, k)] =
            getVelocityZ(pos_w, m_w0.data(), std::vector<int>{Nx, Ny, Nz});
    }
}

void Simulator::advectScalar()
{
    std::copy(m_density.begin(), m_density.end(), m_density0.begin());
    std::copy(m_temperature.begin(), m_temperature.end(), m_temperature0.begin());

    FOR_EACH_CELL
    {
        Vec3 pos_cell = getCenter(i, j, k);
        // Vec3 vel_cell = m_grids->getVelocity(pos_cell);
        Vec3 vel_cell = getVelocity(
            pos_cell,
            std::vector<double *>{m_u0.data(), m_v0.data(), m_w0.data()},
            std::vector<int>{Nx, Ny, Nz});

        pos_cell -= DT * vel_cell;

        m_density[POS(i, j, k)] =
            linearInterpolation<double>(
                pos_cell - 0.5 * Vec3(VOXEL_SIZE, VOXEL_SIZE, VOXEL_SIZE),
                m_density0.data(),
                std::vector<int>{Nx, Ny, Nz},
                VOXEL_SIZE);

        m_temperature[POS(i, j, k)] =
            linearInterpolation<double>(
                pos_cell - 0.5 * Vec3(VOXEL_SIZE, VOXEL_SIZE, VOXEL_SIZE),
                m_temperature0.data(),
                std::vector<int>{Nx, Ny, Nz},
                VOXEL_SIZE);
    }
}