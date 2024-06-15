#ifndef __CUDA_RENDER_H__
#define __CUDA_RENDER_H__

#include "CudaWorker.h"

namespace MCUDA
{
    class CudaRender : public CudaWorker
    {
    public:
        CudaRender(
            int workSize,
            int Nx,
            int Ny,
            int Nz) : CudaWorker(workSize, Nx, Ny, Nz)
        {
        }

        void init();

        void cleanup();

        void genTransparencyMap(
            double *density,
            double *result,
            double light_x, double light_y, double light_z,
            double module_scale_factor,
            double factor);

    private:
        double *transparency;
        double *density;
    };
}

#endif // __CUDA_RENDER_H__