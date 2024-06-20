#ifndef __CPU_Simulator_H__
#define __CPU_Simulator_H__

#include "SimBase.h"
#include "Solver.h"

#include "mmath.h"

namespace CPUSIM
{
    class CpuSimulator
    {
    public:
        CpuSimulator(
            int workSize,
            int Nx, int Ny, int Nz)
            : workSize_(workSize),
              Nx_(Nx), Ny_(Ny), Nz_(Nz),
              b(workSize), x(workSize)
        {
        }

        ~CpuSimulator() = default;
        void init()
        {
            auto L =
                build_3d_laplace<double>(Nx_, Ny_, Nz_);

            solver_.compute(L);
        }
        void cleanup() {}

    private:
        int workSize_;
        int Nx_;
        int Ny_;
        int Nz_;

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

        // presure field
        double *pressure;

        // presure solver
        Eigen::VectorXd b;
        Eigen::VectorXd x;

    public:
        EigenSolver solver_;
        void computeExternalForce(
            double alpha, double beta, double t_ambient);
        void computeVorticity();
        void applyExternalForce();
        void computePressure();
        void applyPressure();
        void advectVelocityField();
        void advectScalarField();

        void setforceField(
            double *f_x,
            double *f_y,
            double *f_z)
        {
            this->f_x = f_x;
            this->f_y = f_y;
            this->f_z = f_z;
        }

        void setVelocityField(
            double *u,
            double *v,
            double *w)
        {
            this->u = u;
            this->v = v;
            this->w = w;
        }

        void setPreviosVelocityField(
            double *u,
            double *v,
            double *w)
        {
            this->u_0 = u;
            this->v_0 = v;
            this->w_0 = w;
        }

        void setDensityField(
            double *density)
        {
            this->density = density;
        }

        void setPreviosDensityField(
            double *density_0)
        {
            this->density_0 = density_0;
        }

        void setTemperatureField(
            double *temperature)
        {
            this->temperature = temperature;
        }

        void setPreviosTemperatureField(
            double *temperature_0)
        {
            this->temperature_0 = temperature_0;
        }

        void setAvgVelocityField(
            double *avg_u,
            double *avg_v,
            double *avg_w)
        {
            this->avg_u = avg_u;
            this->avg_v = avg_v;
            this->avg_w = avg_w;
        }

        void setVorticityField(
            double *omg_x,
            double *omg_y,
            double *omg_z)
        {
            this->omg_x = omg_x;
            this->omg_y = omg_y;
            this->omg_z = omg_z;
        }

        void setPressureField(
            double *pressure)
        {
            this->pressure = pressure;
        }
    };

}

#endif // __CPU_Simulator_H__