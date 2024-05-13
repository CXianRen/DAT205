#ifndef __PARTICLE_SIMULATOR_H__
#define __PARTICLE_SIMULATOR_H__

#include <glm/glm.hpp>

#include "const.h"
#include "../common/debug.h"

class Simulator {

public:
    Simulator(){}
    ~Simulator() = default;

    void init() {
       DEBUG_PRINT("Simulator init");
    }
    void update(){};
    void render(){};

private:
    glm::vec3 u[SIZE]; // velocity field
    glm::vec3 u0[SIZE]; // velocity field at previous time step

    float p[SIZE]; // pressure field
    float temprarure[SIZE]; // temperature field  
};

#endif