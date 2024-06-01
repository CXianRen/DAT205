#include "CudaWorker.h"
#include "common/debug.h"
#include "common/mmath.h"

#include <cuda_runtime.h>

__global__ void addKernel()
{
}

namespace MCUDA
{

    __global__ void calculateAverageVelocityKernel(
        double *u, double *v, double *w,
        double *avg_u, double *avg_v, double *avg_w,
        int workSize, int Nx, int Ny, int Nz)
    {
        // int idx = blockIdx.x * blockDim.x + threadIdx.x;

        // // calculate the index of the cell
        // int i = idx % Nx;
        // int j = (idx / Nx) % Ny;
        // int k = idx / (Nx * Ny);

        // if (idx < workSize)
        // {
        //     // calculate average velocity
        //     avg_u[POS(i, j, k)] = (u[POS_X(i, j, k)] + u[POS_X(i + 1, j, k)]) * 0.5;
        //     avg_v[POS(i, j, k)] = (v[POS_Y(i, j, k)] + v[POS_Y(i, j + 1, k)]) * 0.5;
        //     avg_w[POS(i, j, k)] = (w[POS_Z(i, j, k)] + w[POS_Z(i, j, k + 1)]) * 0.5;
        // }
    }

    CudaWorker::CudaWorker(
        int workSize,
        int Nx,
        int Ny,
        int Nz) : workSize_(workSize), Nx_(Nx), Ny_(Ny), Nz_(Nz)
    {
        // check cuda device properties
        int deviceCount;
        cudaGetDeviceCount(&deviceCount);
        if (deviceCount == 0)
        {
            DEBUG_PRINT("No CUDA devices found");
            exit(1);
        }
        DEBUG_PRINT("CUDA Device Count: " << deviceCount);

        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp, 0);
        DEBUG_PRINT("CUDA Device Name: " << deviceProp.name);
        DEBUG_PRINT("CUDA Compute Capability: " << deviceProp.major << "." << deviceProp.minor);
        DEBUG_PRINT("CUDA Device Memory: " << deviceProp.totalGlobalMem / (1024 * 1024) << "MB");

        // SM count
        int SMCount = deviceProp.multiProcessorCount;
        DEBUG_PRINT("CUDA SM Count: " << SMCount);

        // max grid size
        int maxGridSize[3];
        maxGridSize[0] = deviceProp.maxGridSize[0];
        maxGridSize[1] = deviceProp.maxGridSize[1];
        maxGridSize[2] = deviceProp.maxGridSize[2];
        DEBUG_PRINT("CUDA Max Grid Size: " << maxGridSize[0]
                                           << "x" << maxGridSize[1]
                                           << "x" << maxGridSize[2]);
        // max threads per block
        int maxThreadsPerBlock = deviceProp.maxThreadsPerBlock;
        DEBUG_PRINT("CUDA Max Threads Per Block: " << maxThreadsPerBlock);
        threadsPerBlock_ = maxThreadsPerBlock;

        // max warp per block
        int maxWarpsPerBlock = deviceProp.maxThreadsPerBlock / deviceProp.warpSize;
        DEBUG_PRINT("CUDA Max Warps Per Block: " << maxWarpsPerBlock);

        // warp size
        int warpSize = deviceProp.warpSize;
        DEBUG_PRINT("CUDA Warp Size: " << warpSize);

        //
        DEBUG_PRINT("USING threadPerBlock: " << threadsPerBlock_);
        DEBUG_PRINT("USING workSize: " << workSize_);
        blocksPerGrid_ = (workSize_ + threadsPerBlock_ - 1) / threadsPerBlock_;
        DEBUG_PRINT("USING blockCount: " << blocksPerGrid_);
    }

    CudaWorker::~CudaWorker()
    {
        cleanup();
    }

    void CudaWorker::init()
    {
        DEBUG_PRINT("Initializing CudaWorker- allocate memory");
        // allocate memory
        // why (Nx_ + 1) * Ny_ * Nz : because we need to store ,
        // an extra cell to call the boundary condition
        cudaMalloc(&u, (Nx_ + 1) * Ny_ * Nz * sizeof(double));
        cudaMalloc(&v, Nx_ * (Ny_ + 1) * Nz * sizeof(double));
        cudaMalloc(&w, Nx_ * Ny_ * (Nz + 1) * sizeof(double));

        cudaMalloc(&avg_u, workSize_ * sizeof(double));
        cudaMalloc(&avg_v, workSize_ * sizeof(double));
        cudaMalloc(&avg_w, workSize_ * sizeof(double));

        cudaMalloc(&omg_x, workSize_ * sizeof(double));
        cudaMalloc(&omg_y, workSize_ * sizeof(double));
        cudaMalloc(&omg_z, workSize_ * sizeof(double));

        cudaMalloc(&f_x, workSize_ * sizeof(double));
        cudaMalloc(&f_y, workSize_ * sizeof(double));
        cudaMalloc(&f_z, workSize_ * sizeof(double));
    }

    void CudaWorker::cleanup()
    {
        DEBUG_PRINT("Cleaning up CudaWorker- free memory");
        // free memory
        cudaFree(avg_u);
        cudaFree(avg_v);
        cudaFree(avg_w);

        cudaFree(omg_x);
        cudaFree(omg_y);
        cudaFree(omg_z);

        cudaFree(f_x);
        cudaFree(f_y);
        cudaFree(f_z);
    }

    void CudaWorker::copyDataToDevice(
        double *src, double *dst, int size)
    {
        cudaMemcpy(dst, src, size * sizeof(double), cudaMemcpyHostToDevice);
    }

    void CudaWorker::copyDataToHost(double *src, double *dst, int size)
    {
        cudaMemcpy(dst, src, size * sizeof(double), cudaMemcpyDeviceToHost);
    }

    void CudaWorker::caculateVorticity(double *h_u, double *h_v, double *h_w)
    {
        // sync data
        copyDataToDevice(h_u, u, (Nx_ + 1) * Ny_ * Nz_);
        copyDataToDevice(h_v, v, Nx_ * (Ny_ + 1) * Nz_);
        copyDataToDevice(h_w, w, Nx_ * Ny_ * (Nz_ + 1));

        // calculate average velocity
        // calculateAverageVelocityKernel<<<blocksPerGrid_, threadsPerBlock_>>>(
        //     u, v, w,
        //     avg_u, avg_v, avg_w,
        //     workSize_, Nx_, Ny_, Nz_);
     

        DEBUG_PRINT("Running kernel");
        addKernel<<<4, 256>>>();
  
        // wait for kernel to finish
        DEBUG_PRINT("Waiting for kernel to finish");
        cudaDeviceSynchronize();
    }
}