#ifndef __M_RENDER_H__
#define __M_RENDER_H__
#include <array>

#include <glm/glm.hpp>

#include "const.h"

class SmokeRenderer{

public:
    SmokeRenderer();
    ~SmokeRenderer();

    void render(const std::array<float, gSIZE>& density);

private:
    uint32_t vaoID; // vertex array object id
    uint32_t vboID; // vertex buffer object id
    uint32_t indexID; // index buffer object id (face)
    uint32_t textureID; // texture id

};   
#endif