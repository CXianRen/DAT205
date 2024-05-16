#include "simulator.h"
#include <cmath>
#include <random>
#include <algorithm>
#include "common/globalvar.h"
#include <chrono>

#include "common/debug.h"

void Simulator::addSource()
{
    // @todo visualize the source
    for (int k = (gZ - SOURCE_SIZE_Z) / 2; k < (gZ + SOURCE_SIZE_Z) / 2; ++k)
    {
        // 64-3-3 = 58, 64-3 = 61
        for (int j = gY - SOURCE_Y_MERGIN - SOURCE_SIZE_Y; j < gY - SOURCE_Y_MERGIN; ++j)
        {
            // (32-8) / 2 = 12, (32+8) / 2 = 20
            for (int i = (gX - SOURCE_SIZE_X) / 2; i < (gX + SOURCE_SIZE_X) / 2; ++i)
            {
                density[ACCESS3D(i, j, k)] = INIT_DENSITY;
            }
        }
    }
}

void Simulator::setEmitterVelocity()
{
    for (int k = (gZ - SOURCE_SIZE_Z) / 2; k < (gZ + SOURCE_SIZE_Z) / 2; ++k)
    {
        // @todo pluse 1: why?
        for (int j = gY - SOURCE_Y_MERGIN - SOURCE_SIZE_Y; j < gY - SOURCE_Y_MERGIN + 1; ++j)
        {
            for (int i = (gX - SOURCE_SIZE_X) / 2; i < (gX + SOURCE_SIZE_X) / 2; ++i)
            {
                u[ACCESS3D(i, j, k)] = -INIT_VELOCITY;
                u0[ACCESS3D(i, j, k)] = u[ACCESS3D(i, j, k)];
            }
        }
    }
}

void Simulator::init()
{
    DEBUG_PRINT("Simulator init");
    // init solver
    tripletList.reserve(7 * gSIZE);
    solver.setTolerance(1e-6);

    // init temperature field
    std::random_device rnd;
    std::mt19937 engine(rnd());
    std::uniform_real_distribution<float> dist(0, T_AMP);

    for (int i = 0; i < gX; i++)
    {
        for (int j = 0; j < gY; j++)
        {
            for (int k = 0; k < gZ; k++)
            {
                //@todo: why temprature is higher at the top?
                temperature[ACCESS3D(i, j, k)] =
                    (j / (float)gY) * T_AMP + dist(engine) + T_AMBIENT;
            }
        }
    }
    // init density field
    DEBUG_PRINT("Simulator init density field");
    addSource();
    // init velocity field
    DEBUG_PRINT("Simulator init velocity field");
    setEmitterVelocity();
}

void Simulator::resetForce()
{
    for (int i = 0; i < gX; i++)
        for (int j = 0; j < gY; j++)
            for (int k = 0; k < gZ; k++)
            {
                f[ACCESS3D(i, j, k)] =
                    glm::vec3(
                        0.0,
                        ALPHA * density[ACCESS3D(i, j, k)] - BETA * (temperature[ACCESS3D(i, j, k)] - T_AMBIENT),
                        0.0);
            }
    
}

void Simulator::calculateVorticity()
{
    // average velocity
    for (int i = 1; i < gX; i++)
        for (int j = 1; j < gY; j++)
            for (int k = 1; k < gZ; k++)
            {
                auto &p_u = u[ACCESS3D(i, j, k)];
                auto &p_i1 = u[ACCESS3D(i + 1, j, k)].x;
                auto &p_j1 = u[ACCESS3D(i, j + 1, k)].y;
                auto &p_k1 = u[ACCESS3D(i, j, k + 1)].z;

                avg_u[ACCESS3D(i, j, k)] =
                    glm::vec3(p_u.x + p_i1, p_u.y + p_j1, p_u.z + p_k1) * 0.5f;
            }

    // compute omega (vorticity field)
    for (int i = 1; i < gX; i++)
        for (int j = 1; j < gY; j++)
            for (int k = 1; k < gZ; k++)
            {
                // ignore boundary cells
                if (i == 0 || j == 0 || k == 0)
                {
                    continue;
                }
                if (i == gX - 1 || j == gY - 1 || k == gZ - 1)
                {
                    continue;
                }
                auto &p_omg = omg[ACCESS3D(i, j, k)];

                auto omg_x = (avg_u[ACCESS3D(i, j + 1, k)].z - avg_u[ACCESS3D(i, j - 1, k)].z - (avg_u[ACCESS3D(i, j, k + 1)].y - avg_u[ACCESS3D(i, j, k - 1)].y)) * 0.5f / gCONST_h;

                auto omg_y = (avg_u[ACCESS3D(i, j, k + 1)].x - avg_u[ACCESS3D(i, j, k - 1)].x - (avg_u[ACCESS3D(i + 1, j, k)].z - avg_u[ACCESS3D(i - 1, j, k)].z)) * 0.5f / gCONST_h;

                auto omg_z = (avg_u[ACCESS3D(i + 1, j, k)].y - avg_u[ACCESS3D(i - 1, j, k)].y - (avg_u[ACCESS3D(i, j + 1, k)].x - avg_u[ACCESS3D(i, j - 1, k)].x)) * 0.5f / gCONST_h;
                p_omg = glm::vec3(omg_x, omg_y, omg_z);
            }
    // compute vorticity force
    for (int i = 1; i < gX; i++)
        for (int j = 1; j < gY; j++)
            for (int k = 1; k < gZ; k++)
            {
                // ignore boundary cells
                if (i == 0 || j == 0 || k == 0)
                    continue;

                if (i == gX - 1 || j == gY - 1 || k == gZ - 1)
                    continue;

                // compute gradient of vorticity
                float p, q;
                p = omg[ACCESS3D(i + 1, j, k)].length();
                q = omg[ACCESS3D(i - 1, j, k)].length();
                float grad1 = (p - q) * 0.5 / gCONST_h;

                p = omg[ACCESS3D(i, j + 1, k)].length();
                q = omg[ACCESS3D(i, j - 1, k)].length();
                float grad2 = (p - q) * 0.5 / gCONST_h;

                p = omg[ACCESS3D(i, j, k + 1)].length();
                q = omg[ACCESS3D(i, j, k - 1)].length();
                float grad3 = (p - q) * 0.5 / gCONST_h;

                glm::vec3 gradVort(grad1, grad2, grad3);
                // compute N vector
                glm::vec3 N_ijk(0, 0, 0);
                float norm = gradVort.length();
                if (norm != 0)
                {
                    N_ijk = gradVort / norm;
                }

                glm::vec3 vorticity = omg[ACCESS3D(i, j, k)];

                glm::vec3 ft = glm::cross(vorticity, N_ijk) * VORT_EPS * (float)gCONST_h;
                vort[ACCESS3D(i, j, k)] = ft.length();
                f[ACCESS3D(i, j, k)] = ft;
            }
   
}

void Simulator::addForce()
{
    // u = dt * df * 0.5
    // F=ma, a = F/m, dv = a * dt, m=1 => dv = F * dt
    for (int i = 0; i < gX; i++)
        for (int j = 0; j < gY; j++)
            for (int k = 0; k < gZ; k++)
            {
                auto &ft = f[ACCESS3D(i, j, k)];

                if (i < gX - 1)
                {
                    u[ACCESS3D(i + 1, j, k)].x += DT * (ft.x + f[ACCESS3D(i + 1, j, k)].x) * 0.5;
                }
                if (j < gY - 1)
                {
                    u[ACCESS3D(i, j + 1, k)].y += DT * (ft.y + f[ACCESS3D(i, j + 1, k)].y) * 0.5;
                }
                if (k < gZ - 1)
                {
                    u[ACCESS3D(i, j, k + 1)].z += DT * (ft.z + f[ACCESS3D(i, j, k + 1)].z) * 0.5;
                }
            }
    
}

// from paper, pressure term is solved by linear solver
void Simulator::calPressure()
{
    DEBUG_PRINT("Simulator calPressure");
    tripletList.clear();
    A.setZero();
    b.setZero();
    x.setZero();

    float coeff = gCONST_h / DT;
    // build laplacian matrix
    for (int i = 0; i < gX; i++)
        for (int j = 0; j < gY; j++)
            for (int k = 0; k < gZ; k++)
            {
                float F[6] = {static_cast<float>(k > 0), static_cast<float>(j > 0), static_cast<float>(i > 0),
                              static_cast<float>(i < gX - 1), static_cast<float>(j < gY - 1), static_cast<float>(k < gZ - 1)};
                float D[6] = {-1.0, -1.0, -1.0, 1.0, 1.0, 1.0};
                float U[6];

                U[0] = u[ACCESS3D(i, j, k)].z;
                U[1] = u[ACCESS3D(i, j, k)].y;
                U[2] = u[ACCESS3D(i, j, k)].x;

                U[3] = u[ACCESS3D(i + 1, j, k)].x;
                U[4] = u[ACCESS3D(i, j + 1, k)].y;
                U[5] = u[ACCESS3D(i, j, k + 1)].z;
                float sum_F = 0.0;

                for (int n = 0; n < 6; ++n)
                {
                    sum_F += F[n];
                    b(ACCESS3D(i, j, k)) += D[n] * F[n] * U[n];
                }
                b(ACCESS3D(i, j, k)) *= coeff;

                if (k > 0)
                {
                    tripletList.push_back(Eigen::Triplet<float>(ACCESS3D(i, j, k), ACCESS3D(i, j, k - 1), F[0]));
                }
                if (j > 0)
                {
                    tripletList.push_back(Eigen::Triplet<float>(ACCESS3D(i, j, k), ACCESS3D(i, j - 1, k), F[1]));
                }
                if (i > 0)
                {
                    tripletList.push_back(Eigen::Triplet<float>(ACCESS3D(i, j, k), ACCESS3D(i - 1, j, k), F[2]));
                }

                tripletList.push_back(Eigen::Triplet<float>(ACCESS3D(i, j, k), ACCESS3D(i, j, k), -sum_F));

                if (i < gX - 1)
                {
                    tripletList.push_back(Eigen::Triplet<float>(ACCESS3D(i, j, k), ACCESS3D(i + 1, j, k), F[3]));
                }
                if (j < gY - 1)
                {
                    tripletList.push_back(Eigen::Triplet<float>(ACCESS3D(i, j, k), ACCESS3D(i, j + 1, k), F[4]));
                }
                if (k < gZ - 1)
                {
                    tripletList.push_back(Eigen::Triplet<float>(ACCESS3D(i, j, k), ACCESS3D(i, j, k + 1), F[5]));
                }
            }

    DEBUG_PRINT("Solver build matrix");
    A.setFromTriplets(tripletList.begin(), tripletList.end());

    DEBUG_PRINT("Solver start");
    /* solve sparse lenear system by ICCG solver */
    solver.compute(A);
    if (solver.info() == Eigen::Success)
    {
        printf("SUCCESS: Convergence\n");
    }
    else
    {
        fprintf(stderr, "FAILED: No Convergence\n");
    }
    DEBUG_PRINT("Solver solve");
    x = solver.solve(b);
    printf("#iterations:     %d \n", static_cast<int>(solver.iterations()));
    printf("estimated error: %e \n", solver.error());

    // asign x to m_grids->pressure
    // Eigen::Map<Eigen::VectorXf>(p.begin(), gSIZE) = x;
}

void Simulator::applyPressureTerm()
{
    for (int i = 0; i < gX; i++)
        for (int j = 0; j < gY; j++)
            for (int k = 0; k < gZ; k++)
            {
                // compute gradient of pressure
                if (i < gX - 1)
                {
                    u[ACCESS3D(i + 1, j, k)].x -= DT * (p[ACCESS3D(i + 1, j, k)] - p[ACCESS3D(i, j, k)]) / gCONST_h;
                }
                if (j < gY - 1)
                {
                    u[ACCESS3D(i, j + 1, k)].y -= DT * (p[ACCESS3D(i, j + 1, k)] - p[ACCESS3D(i, j, k)]) / gCONST_h;
                }
                if (k < gZ - 1)
                {
                    u[ACCESS3D(i, j, k + 1)].z -= DT * (p[ACCESS3D(i, j, k + 1)] - p[ACCESS3D(i, j, k)]) / gCONST_h;
                }
            }
    std::copy(u.begin(), u.end(), u0.begin());

}

float Simulator::linearInterpolation(int dim, const std::array<glm::vec3, gSIZE> &data, const glm::vec3 &pt)
{
    glm::vec3 pos;
    // clamp position
    pos[0] = std::min(std::max(0.0f, pt[0]), (float)(gX - 1) * gCONST_h - 1e-6f);
    pos[1] = std::min(std::max(0.0f, pt[1]), (float)(gY - 1) * gCONST_h - 1e-6f);
    pos[2] = std::min(std::max(0.0f, pt[2]), (float)(gZ - 1) * gCONST_h - 1e-6f);

    int i = (int)(pos[0] / gCONST_h);
    int j = (int)(pos[1] / gCONST_h);
    int k = (int)(pos[2] / gCONST_h);

    float scale = 1.0 / gCONST_h;
    float fractx = scale * (pos[0] - i * gCONST_h);
    float fracty = scale * (pos[1] - j * gCONST_h);
    float fractz = scale * (pos[2] - k * gCONST_h);

    assert(fractx < 1.0 && fractx >= 0);
    assert(fracty < 1.0 && fracty >= 0);
    assert(fractz < 1.0 && fractz >= 0);

    // Y @ low X, low Z:
    float tmp1 = data[ACCESS3D(i, j, k)][dim];
    float tmp2 = data[ACCESS3D(i, j + 1, k)][dim];
    // Y @ high X, low Z:
    float tmp3 = data[ACCESS3D(i + 1, j, k)][dim];
    float tmp4 = data[ACCESS3D(i + 1, j + 1, k)][dim];

    // Y @ low X, high Z:
    float tmp5 = data[ACCESS3D(i, j, k + 1)][dim];
    float tmp6 = data[ACCESS3D(i, j + 1, k + 1)][dim];
    // Y @ high X, high Z:
    float tmp7 = data[ACCESS3D(i + 1, j, k + 1)][dim];
    float tmp8 = data[ACCESS3D(i + 1, j + 1, k + 1)][dim];

    // Y @ low X, low Z
    float tmp12 = (1 - fracty) * tmp1 + fracty * tmp2;
    // Y @ high X, low Z
    float tmp34 = (1 - fracty) * tmp3 + fracty * tmp4;

    // Y @ low X, high Z
    float tmp56 = (1 - fracty) * tmp5 + fracty * tmp6;
    // Y @ high X, high Z
    float tmp78 = (1 - fracty) * tmp7 + fracty * tmp8;

    // X @ low Z
    float tmp1234 = (1 - fractx) * tmp12 + fractx * tmp34;
    // X @ high Z
    float tmp5678 = (1 - fractx) * tmp56 + fractx * tmp78;

    // Z
    float tmp = (1 - fractz) * tmp1234 + fractz * tmp5678;


    return tmp;
}

float Simulator::linearInterpolation(const std::array<float, gSIZE> &data, const glm::vec3 &pt)
{
    glm::vec3 pos;
    // clamp position
    pos[0] = std::min(std::max(0.0f, pt[0]), (float)(gX - 1) * gCONST_h - 1e-6f);
    pos[1] = std::min(std::max(0.0f, pt[1]), (float)(gY - 1) * gCONST_h - 1e-6f);
    pos[2] = std::min(std::max(0.0f, pt[2]), (float)(gZ - 1) * gCONST_h - 1e-6f);

    int i = (int)(pos[0] / gCONST_h);
    int j = (int)(pos[1] / gCONST_h);
    int k = (int)(pos[2] / gCONST_h);

    float scale = 1.0 / gCONST_h;
    float fractx = scale * (pos[0] - i * gCONST_h);
    float fracty = scale * (pos[1] - j * gCONST_h);
    float fractz = scale * (pos[2] - k * gCONST_h);

    assert(fractx < 1.0 && fractx >= 0);
    assert(fracty < 1.0 && fracty >= 0);
    assert(fractz < 1.0 && fractz >= 0);

    // Y @ low X, low Z:
    float tmp1 = data[ACCESS3D(i, j, k)];
    float tmp2 = data[ACCESS3D(i, j + 1, k)];
    // Y @ high X, low Z:
    float tmp3 = data[ACCESS3D(i + 1, j, k)];
    float tmp4 = data[ACCESS3D(i + 1, j + 1, k)];

    // Y @ low X, high Z:
    float tmp5 = data[ACCESS3D(i, j, k + 1)];
    float tmp6 = data[ACCESS3D(i, j + 1, k + 1)];
    // Y @ high X, high Z:
    float tmp7 = data[ACCESS3D(i + 1, j, k + 1)];
    float tmp8 = data[ACCESS3D(i + 1, j + 1, k + 1)];

    // Y @ low X, low Z
    float tmp12 = (1 - fracty) * tmp1 + fracty * tmp2;
    // Y @ high X, low Z
    float tmp34 = (1 - fracty) * tmp3 + fracty * tmp4;

    // Y @ low X, high Z
    float tmp56 = (1 - fracty) * tmp5 + fracty * tmp6;
    // Y @ high X, high Z
    float tmp78 = (1 - fracty) * tmp7 + fracty * tmp8;

    // X @ low Z
    float tmp1234 = (1 - fractx) * tmp12 + fractx * tmp34;
    // X @ high Z
    float tmp5678 = (1 - fractx) * tmp56 + fractx * tmp78;

    // Z
    float tmp = (1 - fractz) * tmp1234 + fractz * tmp5678;
    return tmp;
}

glm::vec3 Simulator::getCenter(int i, int j, int k)
{
    static float half_dx = 0.5 * gCONST_h;
    float x = half_dx + i * gCONST_h;
    float y = half_dx + j * gCONST_h;
    float z = half_dx + k * gCONST_h;
    return glm::vec3(x, y, z);
}

#define getVelocityX(pos) linearInterpolation(0, u0, pos - 0.5f * glm::vec3(0.0, gCONST_h, gCONST_h))
#define getVelocityY(pos) linearInterpolation(1, u0, pos - 0.5f * glm::vec3(gCONST_h, 0.0, gCONST_h));
#define getVelocityZ(pos) linearInterpolation(2, u0, pos - 0.5f * glm::vec3(gCONST_h, gCONST_h, 0.0));

glm::vec3 Simulator::getVelocity(const glm::vec3 pos)
{
    glm::vec3 vel_u;
    vel_u.x = getVelocityX(pos);
    vel_u.y = getVelocityY(pos);
    vel_u.z = getVelocityZ(pos);
    return vel_u;
}

void Simulator::advectVelocity()
{
    // semi-langrangian advection
    static float half_dx = 0.5 * gCONST_h;

    for (int i = 0; i < gX; i++)
        for (int j = 0; j < gY; j++)
            for (int k = 0; k < gZ; k++)
            {
                auto center = getCenter(i, j, k);
                // FACE X
                glm::vec3 pos_u = center - 0.5f * glm::vec3(gCONST_h, 0.f, 0.f);
                auto vel_u = getVelocity(pos_u);
                pos_u -= DT * vel_u;
                u[ACCESS3D(i, j, k)].x = getVelocityX(pos_u);

                // FACE Y
                glm::vec3 pos_v = center - 0.5f * glm::vec3(0.f, gCONST_h, 0.f);
                auto vel_v = getVelocity(pos_v);
                pos_v -= vel_v * DT;
                u[ACCESS3D(i, j, k)].y = getVelocityY(pos_v);

                // FACE Z
                glm::vec3 pos_w = center - 0.5f * glm::vec3(0.f, 0.f, gCONST_h);
                auto vel_w = getVelocity(pos_w);
                pos_w -= vel_w * DT;
                u[ACCESS3D(i, j, k)].z = getVelocityZ(pos_w);
            }

}

#define getDensity(pos) linearInterpolation(density0, pos - 0.5f * glm::vec3(gCONST_h, gCONST_h, gCONST_h))
#define getTemperature(pos) linearInterpolation(temperature0, pos - 0.5f * glm::vec3(gCONST_h, gCONST_h, gCONST_h))


void Simulator::advectScalar()
{
    std::copy(density.begin(), density.end(), density0.begin());
    std::copy(temperature.begin(), temperature.end(), temperature0.begin());

    static float half_dx = 0.5 * gCONST_h;

    for (int i = 0; i < gX; i++)
        for (int j = 0; j < gY; j++)
            for (int k = 0; k < gZ; k++)
            {
                auto pos_cell = getCenter(i, j, k);

                auto vel_cell = getVelocity(pos_cell);
                pos_cell -= DT * vel_cell;
                density[ACCESS3D(i, j, k)] = getDensity(pos_cell);
                temperature[ACCESS3D(i, j, k)] = getTemperature(pos_cell);
            }
}

// void Simulator::update()
// {
//     auto start = std::chrono::high_resolution_clock::now();
//     resetForce();
//     auto end = std::chrono::high_resolution_clock::now();

//     auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
   
//     std::cout << "Execution time of resetForce(): " << duration << " milliseconds" << std::endl;

//     start = std::chrono::high_resolution_clock::now();
//     calculateVorticity();
//     end = std::chrono::high_resolution_clock::now();

//     duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

//     std::cout << "Execution time of calculateVorticity(): " << duration << " milliseconds" << std::endl;

//     start = std::chrono::high_resolution_clock::now();
//     addForce();
//     end = std::chrono::high_resolution_clock::now();

//     duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

//     std::cout << "Execution time of addForce(): " << duration << " milliseconds" << std::endl;

//     start = std::chrono::high_resolution_clock::now();
//     calPressure();
//     end = std::chrono::high_resolution_clock::now();

//     duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
 
//     std::cout << "Execution time of calPressure(): " << duration << " milliseconds" << std::endl;

//     start = std::chrono::high_resolution_clock::now();
//     applyPressureTerm();
//     end = std::chrono::high_resolution_clock::now();

//     duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

//     std::cout << "Execution time of applyPressureTerm(): " << duration << " milliseconds" << std::endl;

//     advectVelocity();
    
//     advectScalar();

//     // if (m_time < EMIT_DURATION)
//     // {
//     //     addSource();
//     //     setEmitterVelocity();
//     // }
// }

void Simulator::update()
{

    // resetForce();
 
    // calculateVorticity();

    // addForce();

    // calPressure();

    // applyPressureTerm();

    // advectVelocity();

    // advectScalar();
    // if (m_time < EMIT_DURATION)
    // {
    //     addSource();
    //     setEmitterVelocity();
    // }
}


std::array<float, gSIZE> generateSphereDensity(){
    std::array<float, gSIZE> density;
    int count = 0;
    for (int i = 0; i < gX; i++)
    {
        for (int j = 0; j < gY; j++)
        {
            for (int k = 0; k < gZ; k++)
            {
                float x = (i - gX / 2);
                float y = (j - gY / 2);
                float z = (k - gZ / 2);
                float r = sqrt(x * x + y * y + z * z);
                if (r <=  gX / 4.0)
                {
                    density[ACCESS3D(i, j, k)] = INIT_DENSITY;
                    count ++ ;
                }
                else
                {
                    density[ACCESS3D(i, j, k)] = 0.0;
                }
                // density[ACCESS3D(i, j, k)] = INIT_DENSITY;
            }
        }
    }
    // DEBUG_PRINT("Sphere density count:" << count);
    return density;
}