#ifndef __CUDA_WORKER_H__
#define __CUDA_WORKER_H__

namespace MCUDA
{
    class CudaWorker
    {

    public:
        CudaWorker(
            int workSize,
            int Nx,
            int Ny,
            int Nz);

        ~CudaWorker();
        void init();
        void cleanup();

        void caculateVorticity(
            double *u,
            double *v,
            double *w);

    private:
        void copyDataToDevice(double *src, double *dst, int size);
        void copyDataToHost(double *src, double *dst, int size);

        //
        int workSize_;
        int Nx_, Ny_, Nz_;

        int threadsPerBlock_;
        int blocksPerGrid_;

        // data
        double *u;
        double *v;
        double *w;

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
    };

} // namespace MCUDA
#endif // __CUDA_WORKER_H__