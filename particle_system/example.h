#ifndef __M_EXAMPLE_H__
#define __M_EXAMPLE_H__

#include <string>
#include <map>

#include "common/debug.h"

#include "common/globalvar.h"
#include "particle/Vexelization.h"

#include "labhelper.h"

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
};

std::map<Demo, std::shared_ptr<Example>>
    examples;

void init_examples()
{
    // example empty
    {
        examples[SMOKE_EMPTY] =
            std::make_shared<Example>("empty");
        examples[SMOKE_EMPTY]->obj_vexels.fill(false);
    }

    // example cube
    {
        examples[SMOKE_CUBE] =
            std::make_shared<Example>("cube");
        examples[SMOKE_CUBE]->obj_vexels =
            generate_vexelized_cube((int)(12));
    }

    // example sphere
    {
        examples[SMOKE_SPHERE] =
            std::make_shared<Example>("sphere");
        examples[SMOKE_SPHERE]->obj_vexels =
            generate_vexelized_sphere((int)(10));
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
    }
}

#endif // __M_EXAMPLE_H__