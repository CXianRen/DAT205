#ifndef __VEXELIZATION_H__
#define __VEXELIZATION_H__

#include "constants.h"
#include "math.h"
#include <array>
#include <vector>

#include <glm/glm.hpp>

/**
 * provide function to generate vexelized object
 *
 * the object is represented by a 3D array of bool
 * using 1 to represent the object and 0 to represent the empty space
 */

/* generate a sphere in grid center*/
std::array<bool, SIZE>
generate_vexelized_sphere(int radius);

/* generate a cube in gird center*/
std::array<bool, SIZE>
generate_vexelized_cube(int length);


std::array<bool, SIZE>
generate_vexel(std::vector<glm::vec3> &m_positions, float & max_length);

#endif