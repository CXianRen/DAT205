#pragma once
#include "Vec3.h"
#include "constants.h"
#include "GridData.h"

#define FOR_EACH_CELL                \
    for (int k = 0; k < Nz; ++k)     \
        for (int j = 0; j < Ny; ++j) \
            for (int i = 0; i < Nx; ++i)

#define FOR_EACH_FACE_X              \
    for (int k = 0; k < Nz; ++k)     \
        for (int j = 0; j < Ny; ++j) \
            for (int i = 0; i < Nx + 1; ++i)

#define FOR_EACH_FACE_Y                  \
    for (int k = 0; k < Nz; ++k)         \
        for (int j = 0; j < Ny + 1; ++j) \
            for (int i = 0; i < Nx; ++i)

#define FOR_EACH_FACE_Z              \
    for (int k = 0; k < Nz + 1; ++k) \
        for (int j = 0; j < Ny; ++j) \
            for (int i = 0; i < Nx; ++i)

class MACGrid
{
public:
    MACGrid();
    ~MACGrid();

    Vec3 getCenter(int i, int j, int k);
    // Vec3 getVelocity(const Vec3 &pos);

    // double getVelocityX(const Vec3 &pos);
    // double getVelocityY(const Vec3 &pos);
    // double getVelocityZ(const Vec3 &pos);
    // double getDensity(const Vec3 &pos);
    
    

    
    

    
   
    
    
};
