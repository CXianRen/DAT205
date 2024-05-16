#ifndef __PARTICLE_SIMULATOR_H__
#define __PARTICLE_SIMULATOR_H__

#include "common/debug.h"
#include <array>
#include <glm/glm.hpp>
#include <Eigen/Core>
#include <Eigen/Sparse>

#include "const.h"


class Simulator
{

public:
    Simulator() : A(gSIZE, gSIZE), b(gSIZE), x(gSIZE) {}
    ~Simulator() = default;

    void init();
    void update();
    void render() {};
    float m_time = 0.0f;
    
private:
    void addSource();          // can be seen as a heater, will generate density
    void setEmitterVelocity(); // set the velocity of the emitter

    void resetForce();
    void calculateVorticity(); // calculate the vorticity of the fluid
    void addForce();           // add force to the fluid
    void calPressure();        // calculate the pressure of the fluid
    void applyPressureTerm();  // apply the pressure term to the velocity field
    void advectVelocity();     // advect the velocity field

    void advectScalar(); // advect the scalar field

    std::array<glm::vec3, gSIZE> u;     // velocity field
    std::array<glm::vec3, gSIZE> u0;    // velocity field at previous time step
    std::array<glm::vec3, gSIZE> avg_u; // average velocity field
    std::array<glm::vec3, gSIZE> omg;   // vorticity field  (curl of velocity field)
    std::array<glm::vec3, gSIZE> f;     // force field

    std::array<float, gSIZE> p;            // pressure field
    std::array<float, gSIZE> temperature;  // temperature field
    std::array<float, gSIZE> temperature0; // temperature field at previous time step
    std::array<float, gSIZE> density;      // density field
    std::array<float, gSIZE> density0;     // density field at previous time step
    std::array<float, gSIZE> vort;         // vorticity field, length of f

    glm::vec3 getCenter(int i, int j, int k); // get the center of the voxel

    float linearInterpolation(
        int dim, const std::array<glm::vec3, gSIZE> &data, const glm::vec3 &pt);          // linear interpolation
    float linearInterpolation(const std::array<float, gSIZE> &data, const glm::vec3 &pt); // linear interpolation

    glm::vec3 getVelocity(const glm::vec3 pos); // get the velocity at the point pt

    // solver
    std::vector<Eigen::Triplet<float>> tripletList;
    Eigen::ConjugateGradient<Eigen::SparseMatrix<float>, Eigen::Lower | Eigen::Upper> solver;
    Eigen::SparseMatrix<float, Eigen::RowMajor> A;
    Eigen::VectorXf b;
    Eigen::VectorXf x;
};


// for debug 
std::array<float, gSIZE> generateSphereDensity(); 

#endif