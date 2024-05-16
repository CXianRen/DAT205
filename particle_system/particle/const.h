#ifndef __PARTICLE_CONST_H__
#define __PARTICLE_CONST_H__

#define ACCESS3D(x, y, z) ((x) + (y) * gX + (z) * gX * gY)

// simulation parameters
#define DT 0.02f // time step

// grid size
// a cubic grid is 32x64x32

// #define gN 16
// #define gX gN
// #define gY (2 * gN)
// #define gZ gN

// for testing
#define gN 4
#define gX gN
#define gY (1 * gN)
#define gZ gN


// total number of voxels
#define gSIZE (gX * gY * gZ)

// voxel size
#define gCONST_h 1

// temperature field
#define T_AMP 5.0
#define T_AMBIENT 50.0
#define ALPHA 9.8
#define BETA 15.0

// the size of te sorce box
#define SOURCE_SIZE_X 8   // 8   (int)(gX / 4)
#define SOURCE_SIZE_Y 3   // 3   (int)(gY / 4)
#define SOURCE_SIZE_Z 8   // 8
#define SOURCE_Y_MERGIN 3 // 3

#define EMIT_DURATION 2.0f // 2.0f

#define INIT_DENSITY 1.0


// velocity field
#define INIT_VELOCITY glm::vec3(0.0, 80.0, 0.0)
// vorticity field
#define VORT_EPS 0.25f // vortex confinement epsilon



#endif // __PARTICLE_CONST_H__