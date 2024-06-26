#ifndef __M_EXAMPLE_H__
#define __M_EXAMPLE_H__

#include <string>
#include <map>
#include <cmath>
#include <random>

#include "labhelper.h"
#include "common/debug.h"
#include "common/globalvar.h"
#include "particle/Vexelization.h"
#include "particle/Simulator.h"

class Example
{
public:
    Example(
        const std::string &name)
        : name(name), model(nullptr)
    {
    }

    ~Example()
    {
        if (model)
        {
            labhelper::freeModel(model);
        }
    }

    std::string name;
    glm::mat4 model_matrix;
    labhelper::Model *model;

    std::array<bool, SIZE> obj_vexels;
    Emitter emitter;
};

std::map<Demo, std::shared_ptr<Example>>
    examples;

void bottom_cubic_emitter(
    double *u, double *v, double *w,
    double *u0, double *v0, double *w0,
    double *density,
    double *density0,
    double *temperature,
    double *temperature0,
    double *pressure)
{
    constexpr int SOURCE_SIZE_X = (int)(16);
    constexpr int SOURCE_SIZE_Y = (int)(3);
    constexpr int SOURCE_SIZE_Z = (int)(16);
    constexpr int SOURCE_Y_MERGIN = (int)(3);

    std::random_device rnd;
    std::mt19937 engine(rnd());
    std::uniform_real_distribution<double> dist(800, 1000);

    for (int k = (Nz - SOURCE_SIZE_Z) / 2; k < (Nz + SOURCE_SIZE_Z) / 2; ++k)
    {
        for (int j = SOURCE_Y_MERGIN; j < SOURCE_Y_MERGIN + SOURCE_SIZE_Y; ++j)
        {
            for (int i = (Nx - SOURCE_SIZE_X) / 2; i < (Nx + SOURCE_SIZE_X) / 2; ++i)
            {
                density[ACC3D(i, j, k, Ny, Nx)] = INIT_DENSITY;
                temperature[ACC3D(i, j, k, Ny, Nx)] = dist(engine);
            }
        }
    }
}

void center_sphere_emiiter(
    double *u, double *v, double *w,
    double *u0, double *v0, double *w0,
    double *density,
    double *density0,
    double *temperature,
    double *temperature0,
    double *pressure)
{   
    int radius = 5;
    int x0 = Nx / 2 - radius;
    int x1 = Nx / 2 + radius;
    int y0 = Ny / 2 - radius;
    int y1 = Ny / 2 + radius;
    int z0 = Nz / 2 - radius;
    int z1 = Nz / 2 + radius;

    FOR_EACH_CELL
    {
        if (
            (i - Nx / 2) * (i - Nx / 2) +
                (j - Ny / 2) * (j - Ny / 2) +
                (k - Nz / 2) * (k - Nz / 2) <=
            radius * radius)
        {
            density[ACC3D(i, j, k, Ny, Nx)] = 1.0;
            temperature[ACC3D(i, j, k, Ny, Nx)] = 1000.0;
        }
    }
}

void init_examples()
{
    // example empty
    {
        examples[SMOKE_EMPTY] =
            std::make_shared<Example>("empty");
        examples[SMOKE_EMPTY]->obj_vexels.fill(false);
        examples[SMOKE_EMPTY]->emitter = center_sphere_emiiter;
    }

    // example cube
    {
        examples[SMOKE_CUBE] =
            std::make_shared<Example>("cube");
        examples[SMOKE_CUBE]->obj_vexels =
            generate_vexelized_cube((int)(12));
        examples[SMOKE_CUBE]->emitter = bottom_cubic_emitter;
    }

    // example sphere
    {
        examples[SMOKE_SPHERE] =
            std::make_shared<Example>("sphere");
        examples[SMOKE_SPHERE]->obj_vexels =
            generate_vexelized_sphere((int)(10));
        examples[SMOKE_SPHERE]->emitter = bottom_cubic_emitter;
    }

    // example tree
    {
        examples[SMOKE_TREE] =
            std::make_shared<Example>("tree");

        auto &tree_exp = examples[Demo::SMOKE_TREE];

        tree_exp->model =
            labhelper::loadModelFromOBJ(
                "../scenes/Lowpoly_tree_sample.obj");
        float tree_max_length = 0.0f;
        tree_exp->obj_vexels =
            generate_vexel(
                tree_exp->model->m_positions,
                tree_max_length);
        tree_exp->emitter = bottom_cubic_emitter;
    }
}

#endif // __M_EXAMPLE_H__