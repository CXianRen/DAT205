#include "CudaWorker.h"

#include <cuda_runtime.h>

#include "common/debug.h"
#include "common/mmath.h"
#include "SimBase.h"

#define CUDA_FOR_EACH                                \
    int idx = blockIdx.x * blockDim.x + threadIdx.x; \
    /*calculate the index of the cell*/              \
    int i = idx % Nx;                                \
    int j = (idx / Nx) % Ny;                         \
    int k = idx / (Nx * Ny);

namespace MCUDA
{
    __global__ void calculateAverageVelocityKernel(
        double *u, double *v, double *w,
        double *avg_u, double *avg_v, double *avg_w,
        int workSize, int Nx, int Ny, int Nz)
    {
        CUDA_FOR_EACH;
        if (idx < workSize)
        {
            calculateAverageVelocity(
                i, j, k,
                Nx, Ny, Nz,
                u, v, w,
                avg_u, avg_v, avg_w);
        }
    }

    __global__ void calculateOmgKernel(
        double *avg_u, double *avg_v, double *avg_w,
        double *omg_x, double *omg_y, double *omg_z,
        int workSize, int Nx, int Ny, int Nz)
    {
        CUDA_FOR_EACH
        if (idx < workSize)
        {
            calculateGradient(
                i, j, k,
                Nx, Ny, Nz,
                avg_u, avg_v, avg_w,
                omg_x, omg_y, omg_z);
        }
    }

    __global__ void calculateVorticityForceKernel(
        double *omg_x, double *omg_y, double *omg_z,
        double *f_x, double *f_y, double *f_z,
        int workSize, int Nx, int Ny, int Nz)
    {

        CUDA_FOR_EACH
        if (idx < workSize)
        {
            calculateVorticityBody<double>(
                i, j, k,
                Nx, Ny, Nz,
                omg_x, omg_y, omg_z,
                f_x, f_y, f_z);
        }
    }

    __global__ void applyExternalForceKernel(
        double *u, double *v, double *w,
        double *f_x, double *f_y, double *f_z,
        int workSize, int Nx, int Ny, int Nz)
    {
        CUDA_FOR_EACH
        if (idx < workSize)
        {
            applyExternalForceBody(
                i, j, k,
                Nx, Ny, Nz,
                f_x, f_y, f_z,
                u, v, w);
        }
    }

    __global__ void advectVelocityFieldKernel(
        double *u, double *v, double *w,
        double *u_0, double *v_0, double *w_0,
        int workSize, int Nx, int Ny, int Nz)
    {
        CUDA_FOR_EACH
        if (idx < workSize)
        {
            advectVelocityBody<double>(
                u, v, w,
                u_0, v_0, w_0,
                i, j, k,
                Nx, Ny, Nz);
        }
    }

    __global__ void advectScalarFieldKernel(
        double *field, double *field_0,
        double *u_0, double *v_0, double *w_0,
        int workSize, int Nx, int Ny, int Nz)
    {
        CUDA_FOR_EACH
        if (idx < workSize)
        {
            advectScalarBody<double>(
                i, j, k,
                Nx, Ny, Nz,
                field, field_0,
                u_0, v_0, w_0);
        }
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

    CudaWorker::~CudaWorker()
    {
        cleanup();
    }

    void CudaWorker::init()
    {
        DEBUG_PRINT("Initializing CudaWorker- allocate memory");
        // allocate memory
        cudaMalloc(&u, (Nx_) * Ny_ * Nz * sizeof(double));
        cudaMalloc(&u_0, (Nx_) * Ny_ * Nz * sizeof(double));
        cudaMalloc(&v, Nx_ * (Ny_) * Nz * sizeof(double));
        cudaMalloc(&v_0, Nx_ * (Ny_) * Nz * sizeof(double));
        cudaMalloc(&w, Nx_ * Ny_ * (Nz) * sizeof(double));
        cudaMalloc(&w_0, Nx_ * Ny_ * (Nz) * sizeof(double));

        cudaMalloc(&avg_u, workSize_ * sizeof(double));
        cudaMalloc(&avg_v, workSize_ * sizeof(double));
        cudaMalloc(&avg_w, workSize_ * sizeof(double));

        cudaMalloc(&omg_x, workSize_ * sizeof(double));
        cudaMalloc(&omg_y, workSize_ * sizeof(double));
        cudaMalloc(&omg_z, workSize_ * sizeof(double));

        cudaMalloc(&f_x, workSize_ * sizeof(double));
        cudaMalloc(&f_y, workSize_ * sizeof(double));
        cudaMalloc(&f_z, workSize_ * sizeof(double));

        // allocate memory for density field
        cudaMalloc(&density, workSize_ * sizeof(double));
        cudaMalloc(&density_0, workSize_ * sizeof(double));

        // allocate memory for temperature field
        cudaMalloc(&temperature, workSize_ * sizeof(double));
        cudaMalloc(&temperature_0, workSize_ * sizeof(double));

        // allocate memory for transparency map
        cudaMalloc(&transparency, workSize_ * sizeof(double));
    }

    void CudaWorker::cleanup()
    {
        DEBUG_PRINT("Cleaning up CudaWorker- free memory");
        // free memory
        cudaFree(avg_u);
        cudaFree(avg_v);
        cudaFree(avg_w);
        cudaFree(u);
        cudaFree(v);
        cudaFree(w);
        cudaFree(u_0);
        cudaFree(v_0);
        cudaFree(w_0);

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

    void CudaWorker::setforceField(
        double *f_x,
        double *f_y,
        double *f_z)
    {
        copyDataToDevice(f_x, this->f_x, workSize_);
        copyDataToDevice(f_y, this->f_y, workSize_);
        copyDataToDevice(f_z, this->f_z, workSize_);
    }

    void CudaWorker::getforceField(
        double *f_x,
        double *f_y,
        double *f_z)
    {
        copyDataToHost(this->f_x, f_x, workSize_);
        copyDataToHost(this->f_y, f_y, workSize_);
        copyDataToHost(this->f_z, f_z, workSize_);
    }

    void CudaWorker::setVelocityField(
        double *u_src,
        double *v_src,
        double *w_src)
    {
        copyDataToDevice(u_src, this->u, (Nx_) * Ny_ * Nz_);
        copyDataToDevice(v_src, this->v, Nx_ * (Ny_) * Nz_);
        copyDataToDevice(w_src, this->w, Nx_ * Ny_ * (Nz_));
    }

    void CudaWorker::getVelocityField(
        double *u_dst,
        double *v_dst,
        double *w_dst)
    {
        copyDataToHost(this->u, u_dst, (Nx_) * Ny_ * Nz_);
        copyDataToHost(this->v, v_dst, Nx_ * (Ny_) * Nz_);
        copyDataToHost(this->w, w_dst, Nx_ * Ny_ * (Nz_));
    }

    void CudaWorker::getPreviosVelocityField(
        double *u_dst,
        double *v_dst,
        double *w_dst)
    {
        copyDataToHost(this->u_0, u_dst, (Nx_) * Ny_ * Nz_);
        copyDataToHost(this->v_0, v_dst, Nx_ * (Ny_) * Nz_);
        copyDataToHost(this->w_0, w_dst, Nx_ * Ny_ * (Nz_));
    }

    void CudaWorker::calculateVorticity()
    {
        // calculate average velocity
        // DEBUG_PRINT("Launching kernel");
        calculateAverageVelocityKernel<<<blocksPerGrid_, threadsPerBlock_>>>(
            u, v, w,
            avg_u, avg_v, avg_w,
            workSize_, Nx_, Ny_, Nz_);
        cudaDeviceSynchronize();

        // calculate omg
        calculateOmgKernel<<<blocksPerGrid_, threadsPerBlock_>>>(
            avg_u, avg_v, avg_w,
            omg_x, omg_y, omg_z,
            workSize_, Nx_, Ny_, Nz_);
        cudaDeviceSynchronize();

        // calculate vorticity force
        calculateVorticityForceKernel<<<blocksPerGrid_, threadsPerBlock_>>>(
            omg_x, omg_y, omg_z,
            f_x, f_y, f_z,
            workSize_, Nx_, Ny_, Nz_);
        cudaDeviceSynchronize();
    }

    void CudaWorker::applyExternalForce()
    {
        applyExternalForceKernel<<<blocksPerGrid_, threadsPerBlock_>>>(
            u, v, w,
            f_x, f_y, f_z,
            workSize_, Nx_, Ny_, Nz_);
        cudaError_t error = cudaGetLastError();
        if (error != cudaSuccess)
        {
            DEBUG_PRINT("CUDA Error: " << cudaGetErrorString(error));
        }
        cudaDeviceSynchronize();
    }

    void CudaWorker::advectVelocityField()
    {
        // copy current velocity field to previous
        cudaMemcpy(u_0, u, (Nx_) * Ny_ * Nz_ * sizeof(double), cudaMemcpyDeviceToDevice);
        cudaMemcpy(v_0, v, Nx_ * (Ny_) * Nz_ * sizeof(double), cudaMemcpyDeviceToDevice);
        cudaMemcpy(w_0, w, Nx_ * Ny_ * (Nz_) * sizeof(double), cudaMemcpyDeviceToDevice);

        advectVelocityFieldKernel<<<blocksPerGrid_, threadsPerBlock_>>>(
            u, v, w,
            u_0, v_0, w_0,
            workSize_, Nx_, Ny_, Nz_);
        // check error
        cudaError_t error = cudaGetLastError();
        if (error != cudaSuccess)
        {
            DEBUG_PRINT("Error at advectVelocityFieldKernel: " << cudaGetErrorString(error));
            DEBUG_PRINT("\tworkSize: " << workSize_);
            DEBUG_PRINT("\tNx: " << Nx_);
            DEBUG_PRINT("\tNy: " << Ny_);
            DEBUG_PRINT("\tNz: " << Nz_);
            DEBUG_PRINT("\tblocksPerGrid: " << blocksPerGrid_);
            DEBUG_PRINT("\tthreadsPerBlock: " << threadsPerBlock_);
        }
        cudaDeviceSynchronize();
    }

    void CudaWorker::setDensityField(
        double *density)
    {
        // copy density field to device
        copyDataToDevice(density, this->density, workSize_);
    }

    void CudaWorker::getDensityField(
        double *density)
    {
        // copy density field to host
        copyDataToHost(this->density, density, workSize_);
    }

    void CudaWorker::setPreviosDensityField(
        double *density_0)
    {
        copyDataToDevice(density_0, this->density_0, workSize_);
    }

    void CudaWorker::getPreviosDensityField(
        double *density_0)
    {
        copyDataToHost(this->density_0,
                       density_0, workSize_);
    }

    void CudaWorker::setTemperatureField(
        double *temperature)
    {
        copyDataToDevice(temperature,
                         this->temperature, workSize_);
    }

    void CudaWorker::getTemperatureField(
        double *temperature)
    {
        copyDataToHost(this->temperature,
                       temperature, workSize_);
    }

    void CudaWorker::setPreviosTemperatureField(
        double *temperature_0)
    {
        copyDataToDevice(temperature_0,
                         this->temperature_0,
                         workSize_);
    }

    void CudaWorker::getPreviosTemperatureField(
        double *temperature_0)
    {
        copyDataToHost(this->temperature_0,
                       temperature_0, workSize_);
    }

    void CudaWorker::advectScalarField()
    {
        // copy velocity field to previous
        cudaMemcpy(u_0, u,
                   (Nx_) * Ny_ * Nz_ * sizeof(double),
                   cudaMemcpyDeviceToDevice);
        // copy current density field to previous
        cudaMemcpy(density_0, density,
                   workSize_ * sizeof(double),
                   cudaMemcpyDeviceToDevice);
        // copy current temperature field to previous
        cudaMemcpy(temperature_0, temperature,
                   workSize_ * sizeof(double),
                   cudaMemcpyDeviceToDevice);

        // advect density field
        advectScalarFieldKernel<<<blocksPerGrid_, threadsPerBlock_>>>(
            density, density_0,
            u_0, v_0, w_0,
            workSize_, Nx_, Ny_, Nz_);
        cudaDeviceSynchronize();

        // advect temperature field
        advectScalarFieldKernel<<<blocksPerGrid_, threadsPerBlock_>>>(
            temperature, temperature_0,
            u_0, v_0, w_0,
            workSize_, Nx_, Ny_, Nz_);
        cudaDeviceSynchronize();
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

    void CudaWorker::genTransparencyMap(
        double light_x, double light_y, double light_z,
        double module_scale_factor,
        double factor)
    {
        RayMarchingKernel<<<blocksPerGrid_, threadsPerBlock_>>>(
            transparency,
            density,
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
    }

    void CudaWorker::getTransparencyMap(
        double *transparency)
    {
        copyDataToHost(this->transparency, transparency, workSize_);
    }

}