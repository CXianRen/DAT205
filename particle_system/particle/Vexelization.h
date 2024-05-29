#ifndef __VEXELIZATION_H__
#define __VEXELIZATION_H__

#include "constants.h"
#include "math.h"
#include <array>   

/**
 * provide function to generate vexelized object
 * 
 * the object is represented by a 3D array of bool
 * using 1 to represent the object and 0 to represent the empty space
*/

/* generate a sphere in grid center*/
std::array<bool,SIZE> 
generate_vexelized_sphere(int radius);


/* generate a cube in gird center*/
std::array<bool,SIZE>
generate_vexelized_cube(int length);


#endif