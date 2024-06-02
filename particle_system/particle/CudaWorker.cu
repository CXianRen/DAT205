#include "CudaWorker.h"
#include "common/debug.h"
#include "common/mmath.h"

#include <cuda_runtime.h>

#define VEC3_NORM(x, y, z) \
    sqrt((x) * (x) + (y) * (y) + (z) * (z))

#define VEC_CROSS(x1, y1, z1, x2, y2, z2, x, y, z) \
    x = y1 * z2 - z1 * y2;                         \
    y = z1 * x2 - x1 * z2;                         \
    z = x1 * y2 - y1 * x2;

namespace MCUDA
{
    //@todo might merge with the one in mmath.h
    template <typename T>
    __device__
        T
        cuda_linearInterpolation(
            const T *pt,
            T *src,
            int *dims,
            int *maxXYZ)
    {
        T pos[3];
        // clamp position
        pos[0] = min(
            max((T)0.0, pt[0]),
            (T)(maxXYZ[0]) * VOXEL_SIZE - (T)1e-6);
        pos[1] = min(
            max((T)0.0, pt[1]),
            (T)(maxXYZ[1]) * VOXEL_SIZE - (T)1e-6);
        pos[2] = min(
            max((T)0.0, pt[2]),
            (T)(maxXYZ[2]) * VOXEL_SIZE - (T)1e-6);

        int i = (int)(pos[0] / VOXEL_SIZE);
        int j = (int)(pos[1] / VOXEL_SIZE);
        int k = (int)(pos[2] / VOXEL_SIZE);

        T scale = 1.0 / VOXEL_SIZE;
        T fractx = scale * (pos[0] - i * VOXEL_SIZE);
        T fracty = scale * (pos[1] - j * VOXEL_SIZE);
        T fractz = scale * (pos[2] - k * VOXEL_SIZE);

        assert(fractx < 1.0 && fractx >= 0);
        assert(fracty < 1.0 && fracty >= 0);
        assert(fractz < 1.0 && fractz >= 0);

        // Y @ low X, low Z:
        T tmp1 = src[ACCESS3D(i, j, k)];
        T tmp2 = src[ACCESS3D(i, j + 1, k)];
        // Y @ high X, low Z:
        T tmp3 = src[ACCESS3D(i + 1, j, k)];
        T tmp4 = src[ACCESS3D(i + 1, j + 1, k)];

        // Y @ low X, high Z:
        T tmp5 = src[ACCESS3D(i, j, k + 1)];
        T tmp6 = src[ACCESS3D(i, j + 1, k + 1)];
        // Y @ high X, high Z:
        T tmp7 = src[ACCESS3D(i + 1, j, k + 1)];
        T tmp8 = src[ACCESS3D(i + 1, j + 1, k + 1)];

        // Y @ low X, low Z
        T tmp12 = ((T)(1) - fracty) * tmp1 + fracty * tmp2;
        // Y @ high X, low Z
        T tmp34 = ((T)(1) - fracty) * tmp3 + fracty * tmp4;

        // Y @ low X, high Z
        T tmp56 = ((T)(1) - fracty) * tmp5 + fracty * tmp6;
        // Y @ high X, high Z
        T tmp78 = ((T)(1) - fracty) * tmp7 + fracty * tmp8;

        // X @ low Z
        T tmp1234 = ((T)(1) - fractx) * tmp12 + fractx * tmp34;
        // X @ high Z
        T tmp5678 = ((T)(1) - fractx) * tmp56 + fractx * tmp78;

        // Z
        T tmp = ((T)(1) - fractz) * tmp1234 + fractz * tmp5678;
        return tmp;
    }

    __device__ double cuda_getVelocityX(
        double *pos_u, double *u, int Nx, int Ny, int Nz)
    {
        int dim[3] = {Nx + 1, Ny, Nz};
        int maxIndex[3] = {Nx, Ny - 1, Nz - 1};
        double pos_t[3];
        pos_t[0] = pos_u[0];
        pos_t[1] = pos_u[1] - 0.5 * VOXEL_SIZE;
        pos_t[2] = pos_u[2] - 0.5 * VOXEL_SIZE;
        return cuda_linearInterpolation<double>(
            pos_t,
            u,
            dim,
            maxIndex);
    }

    __device__ double cuda_getVelocityY(
        double *pos_v, double *v, int Nx, int Ny, int Nz)
    {
        int dim[3] = {Nx, Ny + 1, Nz};
        int maxIndex[3] = {Nx - 1, Ny, Nz - 1};
        double pos_t[3];
        pos_t[0] = pos_v[0] - 0.5 * VOXEL_SIZE;
        pos_t[1] = pos_v[1];
        pos_t[2] = pos_v[2] - 0.5 * VOXEL_SIZE;
        return cuda_linearInterpolation<double>(
            pos_t,
            v,
            dim,
            maxIndex);
    }

    __device__ double cuda_getVelocityZ(
        double *pos_w, double *w, int Nx, int Ny, int Nz)
    {
        int dim[3] = {Nx, Ny, Nz + 1};
        int maxIndex[3] = {Nx - 1, Ny - 1, Nz};
        double pos_t[3];
        pos_t[0] = pos_w[0] - 0.5 * VOXEL_SIZE;
        pos_t[1] = pos_w[1] - 0.5 * VOXEL_SIZE;
        pos_t[2] = pos_w[2];
        return cuda_linearInterpolation<double>(
            pos_t,
            w,
            dim,
            maxIndex);
    }

    __device__ double cuda_getVelocity(
        double *pos,
        double *vel,
        double *u,
        double *v,
        double *w,
        int Nx, int Ny, int Nz)
    {
        vel[0] = cuda_getVelocityX(pos, u, Nx, Ny, Nz);
        vel[1] = cuda_getVelocityY(pos, v, Nx, Ny, Nz);
        vel[2] = cuda_getVelocityZ(pos, w, Nx, Ny, Nz);
    }

    __global__ void calculateAverageVelocityKernel(
        double *u, double *v, double *w,
        double *avg_u, double *avg_v, double *avg_w,
        int workSize, int Nx, int Ny, int Nz)
    {
        int idx = blockIdx.x * blockDim.x + threadIdx.x;

        // calculate the index of the cell
        int i = idx % Nx;
        int j = (idx / Nx) % Ny;
        int k = idx / (Nx * Ny);

        if (idx < workSize)
        {
            // calculate average velocity
            avg_u[POS(i, j, k)] = (u[POS_X(i, j, k)] + u[POS_X(i + 1, j, k)]) * 0.5;
            avg_v[POS(i, j, k)] = (v[POS_Y(i, j, k)] + v[POS_Y(i, j + 1, k)]) * 0.5;
            avg_w[POS(i, j, k)] = (w[POS_Z(i, j, k)] + w[POS_Z(i, j, k + 1)]) * 0.5;
        }
    }

    __global__ void calculateOmgKernel(
        double *avg_u, double *avg_v, double *avg_w,
        double *omg_x, double *omg_y, double *omg_z,
        int workSize, int Nx, int Ny, int Nz)
    {
        int idx = blockIdx.x * blockDim.x + threadIdx.x;

        int i = idx % Nx;
        int j = (idx / Nx) % Ny;
        int k = idx / (Nx * Ny);

        if (idx < workSize)
        {
            // ignore boundary cells
            if (i == 0 || j == 0 || k == 0)
            {
                return;
            }
            if (i == Nx - 1 || j == Ny - 1 || k == Nz - 1)
            {
                return;
            }
            // calculate vorticity
            omg_x[POS(i, j, k)] = (avg_w[POS(i, j + 1, k)] - avg_w[POS(i, j - 1, k)] - avg_v[POS(i, j, k + 1)] + avg_v[POS(i, j, k - 1)]) * 0.5 / VOXEL_SIZE;
            omg_y[POS(i, j, k)] = (avg_u[POS(i, j, k + 1)] - avg_u[POS(i, j, k - 1)] - avg_w[POS(i + 1, j, k)] + avg_w[POS(i - 1, j, k)]) * 0.5 / VOXEL_SIZE;
            omg_z[POS(i, j, k)] = (avg_v[POS(i + 1, j, k)] - avg_v[POS(i - 1, j, k)] - avg_u[POS(i, j + 1, k)] + avg_u[POS(i, j - 1, k)]) * 0.5 / VOXEL_SIZE;
        }
    }

    __global__ void calculateVorticityForceKernel(
        double *omg_x, double *omg_y, double *omg_z,
        double *f_x, double *f_y, double *f_z,
        int workSize, int Nx, int Ny, int Nz)
    {
        int idx = blockIdx.x * blockDim.x + threadIdx.x;

        int i = idx % Nx;
        int j = (idx / Nx) % Ny;
        int k = idx / (Nx * Ny);

        if (idx < workSize)
        {
            // ignore boundary cells
            if (i == 0 || j == 0 || k == 0)
            {
                return;
            }
            if (i == Nx - 1 || j == Ny - 1 || k == Nz - 1)
            {
                return;
            }
            // compute gradient of vorticity
            double p, q;
            p = VEC3_NORM(omg_x[POS(i + 1, j, k)], omg_y[POS(i + 1, j, k)], omg_z[POS(i + 1, j, k)]);
            q = VEC3_NORM(omg_x[POS(i - 1, j, k)], omg_y[POS(i - 1, j, k)], omg_z[POS(i - 1, j, k)]);

            double grad1 = (p - q) / (2.0 * VOXEL_SIZE);

            p = VEC3_NORM(omg_x[POS(i, j + 1, k)], omg_y[POS(i, j + 1, k)], omg_z[POS(i, j + 1, k)]);
            q = VEC3_NORM(omg_x[POS(i, j - 1, k)], omg_y[POS(i, j - 1, k)], omg_z[POS(i, j - 1, k)]);
            double grad2 = (p - q) / (2.0 * VOXEL_SIZE);

            p = VEC3_NORM(omg_x[POS(i, j, k + 1)], omg_y[POS(i, j, k + 1)], omg_z[POS(i, j, k + 1)]);
            q = VEC3_NORM(omg_x[POS(i, j, k - 1)], omg_y[POS(i, j, k - 1)], omg_z[POS(i, j, k - 1)]);
            double grad3 = (p - q) / (2.0 * VOXEL_SIZE);

            double norm = VEC3_NORM(grad1, grad2, grad3);

            double ni = 0.0, nj = 0.0, nk = 0.0;
            if (norm != 0)
            {
                ni = grad1 / norm;
                nj = grad2 / norm;
                nk = grad3 / norm;
            }

            double f1, f2, f3;

            VEC_CROSS(
                omg_x[POS(i, j, k)],
                omg_y[POS(i, j, k)],
                omg_z[POS(i, j, k)],
                ni, nj, nk,
                f1, f2, f3);

            f_x[POS(i, j, k)] += VORT_EPS * VOXEL_SIZE * f1;
            f_y[POS(i, j, k)] += VORT_EPS * VOXEL_SIZE * f2;
            f_z[POS(i, j, k)] += VORT_EPS * VOXEL_SIZE * f3;
        }
    }

    __global__ void applyExternalForceKernel(
        double *u, double *v, double *w,
        double *f_x, double *f_y, double *f_z,
        int workSize, int Nx, int Ny, int Nz)
    {
        int idx = blockIdx.x * blockDim.x + threadIdx.x;

        int i = idx % Nx;
        int j = (idx / Nx) % Ny;
        int k = idx / (Nx * Ny);

        if (idx < workSize)
        {
            if (i < Nx - 1)
            {
                u[POS_X(i + 1, j, k)] += DT * (f_x[POS(i, j, k)] + f_x[POS(i + 1, j, k)]) * 0.5;
            }
            if (j < Ny - 1)
            {
                v[POS_Y(i, j + 1, k)] += DT * (f_y[POS(i, j, k)] + f_x[POS(i, j + 1, k)]) * 0.5;
            }
            if (k < Nz - 1)
            {
                w[POS_Z(i, j, k + 1)] += DT * (f_z[POS(i, j, k)] + f_x[POS(i, j, k + 1)]) * 0.5;
            }
        }
    }

    __global__ void advectVelocityFieldKernel(
        double *u, double *v, double *w,
        double *u_0, double *v_0, double *w_0,
        int workSize, int Nx, int Ny, int Nz)
    {
        int idx = blockIdx.x * blockDim.x + threadIdx.x;
        int i = idx % Nx;
        int j = (idx / Nx) % Ny;
        int k = idx / (Nx * Ny);

        if (idx < workSize)
        {
            double half_dx = 0.5 * VOXEL_SIZE;

            double center[3];
            center[0] = half_dx + i * VOXEL_SIZE;
            center[1] = half_dx + j * VOXEL_SIZE;
            center[2] = half_dx + k * VOXEL_SIZE;

            double pos_u[3];
            pos_u[0] = center[0] - half_dx;
            pos_u[1] = center[1];
            pos_u[2] = center[2];

            double pos_v[3];
            pos_v[0] = center[0];
            pos_v[1] = center[1] - half_dx;
            pos_v[2] = center[2];

            double pos_w[3];
            pos_w[0] = center[0];
            pos_w[1] = center[1];
            pos_w[2] = center[2] - half_dx;

            // advect u

            double vel_u[3];
            cuda_getVelocity(
                pos_u,
                vel_u,
                u_0, v_0, w_0, Nx, Ny, Nz);

            double vel_v[3];
            cuda_getVelocity(
                pos_v,
                vel_v,
                u_0, v_0, w_0, Nx, Ny, Nz);

            double vel_w[3];
            cuda_getVelocity(
                pos_w,
                vel_w,
                u_0, v_0, w_0, Nx, Ny, Nz);

            pos_u[0] -= DT * vel_u[0];
            pos_u[1] -= DT * vel_u[1];
            pos_u[2] -= DT * vel_u[2];

            pos_v[0] -= DT * vel_v[0];
            pos_v[1] -= DT * vel_v[1];
            pos_v[2] -= DT * vel_v[2];

            pos_w[0] -= DT * vel_w[0];
            pos_w[1] -= DT * vel_w[1];
            pos_w[2] -= DT * vel_w[2];

            u[POS_X(i, j, k)] =
                cuda_getVelocityX(pos_u, u_0, Nx, Ny, Nz);

            v[POS_Y(i, j, k)] =
                cuda_getVelocityY(pos_v, v_0, Nx, Ny, Nz);

            w[POS_Z(i, j, k)] =
                cuda_getVelocityZ(pos_w, w_0, Nx, Ny, Nz);
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
        cudaMalloc(&u_0, (Nx_ + 1) * Ny_ * Nz * sizeof(double));
        cudaMalloc(&v, Nx_ * (Ny_ + 1) * Nz * sizeof(double));
        cudaMalloc(&v_0, Nx_ * (Ny_ + 1) * Nz * sizeof(double));
        cudaMalloc(&w, Nx_ * Ny_ * (Nz + 1) * sizeof(double));
        cudaMalloc(&w_0, Nx_ * Ny_ * (Nz + 1) * sizeof(double));

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
        double *u,
        double *v,
        double *w)
    {
        copyDataToDevice(u, this->u, (Nx_ + 1) * Ny_ * Nz_);
        copyDataToDevice(v, this->v, Nx_ * (Ny_ + 1) * Nz_);
        copyDataToDevice(w, this->w, Nx_ * Ny_ * (Nz_ + 1));
    }

    void CudaWorker::getVelocityField(
        double *u,
        double *v,
        double *w)
    {
        copyDataToHost(this->u, u, (Nx_ + 1) * Ny_ * Nz_);
        copyDataToHost(this->v, v, Nx_ * (Ny_ + 1) * Nz_);
        copyDataToHost(this->w, w, Nx_ * Ny_ * (Nz_ + 1));
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
        cudaDeviceSynchronize();
    }

    void CudaWorker::advectVelocityField()
    {
        // swap u and u_0
        double *temp = u;
        u = u_0;
        u_0 = temp;

        // swap v and v_0
        temp = v;
        v = v_0;
        v_0 = temp;

        // swap w and w_0
        temp = w;
        w = w_0;
        w_0 = temp;

        advectVelocityFieldKernel<<<blocksPerGrid_, threadsPerBlock_>>>(
            u, v, w,
            u_0, v_0, w_0,
            workSize_, Nx_, Ny_, Nz_);
        cudaDeviceSynchronize();
    }

}

__global__ void addKernel()
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
}

void test_cuda()
{
    // launch kernel
    DEBUG_PRINT("Launching kernel\n");
    addKernel<<<4, 256>>>();
    cudaDeviceSynchronize();
    DEBUG_PRINT("Kernel finished\n");
}