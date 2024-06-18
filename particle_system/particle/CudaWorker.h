#ifndef __CUDA_WORKER_H__
#define __CUDA_WORKER_H__

#include "common/debug.h"
#include <cuda_runtime.h>

#define CUDA_FOR_EACH                                \
    int idx = blockIdx.x * blockDim.x + threadIdx.x; \
    /*compute the index of the cell*/              \
    int i = idx % Nx;                                \
    int j = (idx / Nx) % Ny;                         \
    int k = idx / (Nx * Ny);

namespace MCUDA
{

    class CudaWorker
    {
    protected:
        CudaWorker(
            int workSize,
            int Nx,
            int Ny,
            int Nz) : workSize_(workSize),
                      Nx_(Nx), Ny_(Ny), Nz_(Nz)
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
            cudaSetDevice(0);

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
            threadsPerBlock_ = 256;

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

        ~CudaWorker() = default;
        //
        int workSize_;
        int Nx_, Ny_, Nz_;

        int threadsPerBlock_;
        int blocksPerGrid_;
    };
}
#endif // __CUDA_WORKER_H__