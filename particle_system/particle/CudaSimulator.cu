#include "CudaSimulator.h"
#include "common/debug.h"
#include "common/mmath.h"
#include "SimBase.h"

namespace MCUDA
{
    __global__ void calculateExternalForceKernel(
        double *f_x, double *f_y, double *f_z,
        double *density, double *temperature,
        int workSize, int Nx, int Ny, int Nz)
    {
        CUDA_FOR_EACH
        if (idx < workSize)
        {
            calculateBuyancyForceBody(
                i, j, k,
                Nx, Ny, Nz,
                density, temperature,
                f_x, f_y, f_z);
        }
    }

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

    CudaSimulator::CudaSimulator(
        int workSize,
        int Nx,
        int Ny,
        int Nz) : CudaWorker(workSize, Nx, Ny, Nz)
    {
    }

    CudaSimulator::~CudaSimulator()
    {
        cleanup();
    }

    void CudaSimulator::init()
    {
        DEBUG_PRINT("Initializing CudaSimulator- allocate memory");
        // allocate memory
        cudaMalloc(&u, (Nx_)*Ny_ * Nz * sizeof(double));
        cudaMalloc(&u_0, (Nx_)*Ny_ * Nz * sizeof(double));
        cudaMalloc(&v, Nx_ * (Ny_)*Nz * sizeof(double));
        cudaMalloc(&v_0, Nx_ * (Ny_)*Nz * sizeof(double));
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

        // allocate memory for pressure
        cudaMalloc(&pressure, workSize_ * sizeof(double));
        // fill pressure with zeros
        cudaMemset(pressure, 0, workSize_ * sizeof(double));

        // allocate memory for rhs
        cudaMalloc(&rhs_, workSize_ * sizeof(double));
        // init solver
        static auto L = build_3d_laplace<double>(Nx_, Ny_, Nz_);
        DEBUG_PRINT("Solver initialized");
        solver_.compute(L);
    }

    void CudaSimulator::cleanup()
    {
        DEBUG_PRINT("Cleaning up CudaSimulator- free memory");
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

        cudaFree(density);
        cudaFree(density_0);

        cudaFree(temperature);
        cudaFree(temperature_0);

        cudaFree(pressure);
        cudaFree(rhs_);
    }

    void CudaSimulator::copyDataToDevice(
        double *src, double *dst, int size)
    {
        cudaMemcpy(dst, src, size * sizeof(double), cudaMemcpyHostToDevice);
    }

    void CudaSimulator::copyDataToHost(double *src, double *dst, int size)
    {
        cudaMemcpy(dst, src, size * sizeof(double), cudaMemcpyDeviceToHost);
    }

    void CudaSimulator::setforceField(
        double *f_x,
        double *f_y,
        double *f_z)
    {
        copyDataToDevice(f_x, this->f_x, workSize_);
        copyDataToDevice(f_y, this->f_y, workSize_);
        copyDataToDevice(f_z, this->f_z, workSize_);
    }

    void CudaSimulator::getforceField(
        double *f_x,
        double *f_y,
        double *f_z)
    {
        copyDataToHost(this->f_x, f_x, workSize_);
        copyDataToHost(this->f_y, f_y, workSize_);
        copyDataToHost(this->f_z, f_z, workSize_);
    }

    void CudaSimulator::setVelocityField(
        double *u_src,
        double *v_src,
        double *w_src)
    {
        copyDataToDevice(u_src, this->u, (Nx_)*Ny_ * Nz_);
        copyDataToDevice(v_src, this->v, Nx_ * (Ny_)*Nz_);
        copyDataToDevice(w_src, this->w, Nx_ * Ny_ * (Nz_));
    }

    void CudaSimulator::getVelocityField(
        double *u_dst,
        double *v_dst,
        double *w_dst)
    {
        copyDataToHost(this->u, u_dst, (Nx_)*Ny_ * Nz_);
        copyDataToHost(this->v, v_dst, Nx_ * (Ny_)*Nz_);
        copyDataToHost(this->w, w_dst, Nx_ * Ny_ * (Nz_));
    }

    void CudaSimulator::setPreviosVelocityField(
        double *u,
        double *v,
        double *w)
    {
        copyDataToDevice(u, this->u_0, (Nx_)*Ny_ * Nz_);
        copyDataToDevice(v, this->v_0, Nx_ * (Ny_)*Nz_);
        copyDataToDevice(w, this->w_0, Nx_ * Ny_ * (Nz_));
    }

    void CudaSimulator::getPreviosVelocityField(
        double *u_dst,
        double *v_dst,
        double *w_dst)
    {
        copyDataToHost(this->u_0, u_dst, (Nx_)*Ny_ * Nz_);
        copyDataToHost(this->v_0, v_dst, Nx_ * (Ny_)*Nz_);
        copyDataToHost(this->w_0, w_dst, Nx_ * Ny_ * (Nz_));
    }

    void CudaSimulator::calculateExternalForce()
    {
        calculateExternalForceKernel<<<blocksPerGrid_, threadsPerBlock_>>>(
            f_x, f_y, f_z,
            density, temperature,
            workSize_, Nx_, Ny_, Nz_);
        cudaDeviceSynchronize();
    }

    void CudaSimulator::calculateVorticity()
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

    void CudaSimulator::applyExternalForce()
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

    void CudaSimulator::advectVelocityField()
    {
        // copy current velocity field to previous
        cudaMemcpy(u_0, u, (Nx_)*Ny_ * Nz_ * sizeof(double), cudaMemcpyDeviceToDevice);
        cudaMemcpy(v_0, v, Nx_ * (Ny_)*Nz_ * sizeof(double), cudaMemcpyDeviceToDevice);
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

    void CudaSimulator::setDensityField(
        double *density)
    {
        // copy density field to device
        copyDataToDevice(density, this->density, workSize_);
    }

    void CudaSimulator::getDensityField(
        double *density)
    {
        // copy density field to host
        copyDataToHost(this->density, density, workSize_);
    }

    void CudaSimulator::setPreviosDensityField(
        double *density_0)
    {
        copyDataToDevice(density_0, this->density_0, workSize_);
    }

    void CudaSimulator::getPreviosDensityField(
        double *density_0)
    {
        copyDataToHost(this->density_0,
                       density_0, workSize_);
    }

    void CudaSimulator::setTemperatureField(
        double *temperature)
    {
        copyDataToDevice(temperature,
                         this->temperature, workSize_);
    }

    void CudaSimulator::getTemperatureField(
        double *temperature)
    {
        copyDataToHost(this->temperature,
                       temperature, workSize_);
    }

    void CudaSimulator::setPreviosTemperatureField(
        double *temperature_0)
    {
        copyDataToDevice(temperature_0,
                         this->temperature_0,
                         workSize_);
    }

    void CudaSimulator::getPreviosTemperatureField(
        double *temperature_0)
    {
        copyDataToHost(this->temperature_0,
                       temperature_0, workSize_);
    }

    void CudaSimulator::advectScalarField()
    {
        // copy velocity field to previous
        cudaMemcpy(u_0, u,
                   (Nx_)*Ny_ * Nz_ * sizeof(double),
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

    __global__ void buildRhsKernel(
        double *rhs,
        double *u, double *v, double *w,
        int workSize, int Nx, int Ny, int Nz)
    {
        CUDA_FOR_EACH
        if (idx < workSize)
        {
            rhs[idx] = 0.0;
            buildRhsBody<double>(
                i, j, k,
                Nx, Ny, Nz,
                u, v, w,
                rhs);
        }
    }

    void CudaSimulator::calculatePressure()
    {
        // build rhs
        buildRhsKernel<<<blocksPerGrid_, threadsPerBlock_>>>(
            rhs_, u, v, w,
            workSize_, Nx_, Ny_, Nz_);
        // check error
        cudaError_t error = cudaGetLastError();
        if (error != cudaSuccess)
        {
            DEBUG_PRINT("Error at buildRhsKernel: " << cudaGetErrorString(error));
        }
        cudaDeviceSynchronize();

        // solve pressure
        solver_.solve_from_gpu(pressure, rhs_);
    }

    void CudaSimulator::getPressureField(
        double *pressure)
    {
        copyDataToHost(this->pressure, pressure, workSize_);
    }

    void CudaSimulator::getRhsField(
        double *rhs)
    {
        copyDataToHost(rhs_, rhs, workSize_);
    }

    __global__ void applyPressureKernel(
        double *u, double *v, double *w,
        double *pressure,
        int workSize, int Nx, int Ny, int Nz)
    {
        CUDA_FOR_EACH
        if (idx < workSize)
        {
            applyPressureBody<double>(
                i, j, k,
                Nx, Ny, Nz,
                pressure,
                u, v, w);
        }
    }

    void CudaSimulator::applyPressure()
    {
        applyPressureKernel<<<blocksPerGrid_, threadsPerBlock_>>>(
            u, v, w,
            pressure,
            workSize_, Nx_, Ny_, Nz_);
        // check error
        cudaError_t error = cudaGetLastError();
        if (error != cudaSuccess)
        {
            DEBUG_PRINT("Error at applyPressureKernel: " << cudaGetErrorString(error));
        }

        cudaDeviceSynchronize();
    }
}