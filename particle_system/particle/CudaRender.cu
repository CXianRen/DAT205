#include "CudaRender.h"
#include "common/debug.h"
#include "common/mmath.h"
#include "SimBase.h"

namespace MCUDA
{

    void CudaRender::init()
    {
        DEBUG_PRINT("Initializing CudaRender- allocate memory");
        // allocate memory for transparency map
        cudaMalloc(&transparency, workSize_ * sizeof(double));
        cudaMalloc(&density, workSize_ * sizeof(double));
    }

    void CudaRender::cleanup()
    {
        DEBUG_PRINT("Cleaning up CudaRender- free memory");
        // free memory
        cudaFree(transparency);
        cudaFree(density);
    }

    __global__ void RayMarchingKernel(
        double *result,
        double *density,
        int Nx, int Ny, int Nz,
        double light_x, double light_y, double light_z,
        double module_scale_factor,
        double factor)
    {
        CUDA_FOR_EACH
        if (idx < Nx * Ny * Nz)
        {
            double sample_count = Nx;
            double step = 1.0;

            double pos_cell[3];
            getCenter<double>(i, j, k, pos_cell);

            double dir[3] = {
                light_x - pos_cell[0] / Nx * module_scale_factor,
                light_y - pos_cell[1] / Ny * module_scale_factor,
                light_z - pos_cell[2] / Nz * module_scale_factor};

            // normalize direction
            double norm = VEC3_NORM(dir[0], dir[1], dir[2]);
            dir[0] /= norm;
            dir[1] /= norm;
            dir[2] /= norm;

            double Tl = 1.0;
            for (int s = 0; s < sample_count; s++)
            {
                // update position
                pos_cell[0] += dir[0] * step;
                pos_cell[1] += dir[1] * step;
                pos_cell[2] += dir[2] * step;

                // if pos is out of bound
                if (pos_cell[0] > Nx ||
                    pos_cell[1] > Ny ||
                    pos_cell[2] > Nz)
                {
                    break;
                }

                // interpolate density
                double d = getScalar<double>(
                    pos_cell,
                    density,
                    Nx, Ny, Nz);

                if (d < 0.01)
                {
                    continue;
                }

                Tl *= exp(-factor * d * step / Nx);
                if (Tl < 0.01)
                {
                    break;
                }
            }
            result[idx] = Tl;
        }
    }

    void CudaRender::genTransparencyMap(
        double *density,
        double *result,
        double light_x, double light_y, double light_z,
        double module_scale_factor,
        double factor)
    {
        // copy density to device
        cudaMemcpy(this->density, density, workSize_ * sizeof(double), cudaMemcpyHostToDevice);
        // run kernel
        RayMarchingKernel<<<blocksPerGrid_, threadsPerBlock_>>>(
            transparency,
            this->density,
            Nx_, Ny_, Nz_,
            light_x, light_y, light_z,
            module_scale_factor,
            factor);
        // check error
        cudaError_t error = cudaGetLastError();
        if (error != cudaSuccess)
        {
            DEBUG_PRINT("Error at RayMarchingKernel: " << cudaGetErrorString(error));
            DEBUG_PRINT("\tworkSize: " << workSize_);
            DEBUG_PRINT("\tNx: " << Nx_);
            DEBUG_PRINT("\tNy: " << Ny_);
            DEBUG_PRINT("\tNz: " << Nz_);
            DEBUG_PRINT("\tblocksPerGrid: " << blocksPerGrid_);
            DEBUG_PRINT("\tthreadsPerBlock: " << threadsPerBlock_);
        }

        cudaDeviceSynchronize();
        // copy result back
        cudaMemcpy(result, transparency, workSize_ * sizeof(double), cudaMemcpyDeviceToHost);
    }
}