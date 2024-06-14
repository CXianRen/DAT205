#ifndef __CUDA_WORKER_H__
#define __CUDA_WORKER_H__

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
                      Nx_(Nx), Ny_(Ny), Nz_(Nz) {}

        ~CudaWorker() = default;
        //
        int workSize_;
        int Nx_, Ny_, Nz_;

        int threadsPerBlock_;
        int blocksPerGrid_;
    };

    class CudaSimulator : CudaWorker
    {

    public:
        CudaSimulator(
            int workSize,
            int Nx,
            int Ny,
            int Nz);

        ~CudaSimulator();
        void init();
        void cleanup();

        void setforceField(
            double *f_x,
            double *f_y,
            double *f_z);

        void getforceField(
            double *f_x,
            double *f_y,
            double *f_z);

        void setVelocityField(
            double *u,
            double *v,
            double *w);

        void getVelocityField(
            double *u,
            double *v,
            double *w);

        void getPreviosVelocityField(
            double *u,
            double *v,
            double *w);

        void setDensityField(
            double *density);

        void getDensityField(
            double *density);

        void setPreviosDensityField(
            double *density_0);

        void getPreviosDensityField(
            double *density_0);

        void setTemperatureField(
            double *temperature);

        void getTemperatureField(
            double *temperature);

        void setPreviosTemperatureField(
            double *temperature_0);

        void getPreviosTemperatureField(
            double *temperature_0);

        void calculateVorticity();

        void applyExternalForce();

        void advectVelocityField();

        void advectScalarField();

        void genTransparencyMap(
            double light_x, double light_y, double light_z,
            double module_scale_factor,
            double factor);

        void getTransparencyMap(
            double *transparency);

    private:
        void copyDataToDevice(double *src, double *dst, int size);
        void copyDataToHost(double *src, double *dst, int size);

        // data
        double *u;
        double *v;
        double *w;

        double *u_0;
        double *v_0;
        double *w_0;

        double *avg_u;
        double *avg_v;
        double *avg_w;

        double *omg_x;
        double *omg_y;
        double *omg_z;

        // forece  field
        double *f_x;
        double *f_y;
        double *f_z;

        // scalar field
        double *temperature;
        double *temperature_0;

        double *density;
        double *density_0;

        double *transparency;
    };

} // namespace MCUDA

#endif // __CUDA_WORKER_H__