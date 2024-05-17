#ifndef __M_RENDER_H__
#define __M_RENDER_H__
#include <array>

#include <glm/glm.hpp>

#include "const.h"

#define gX 16
#define gY 32
#define gZ 16
#define gSIZE (gX * gY * gZ)

class SmokeRenderer{

public:
    SmokeRenderer();
    ~SmokeRenderer();

    void render(const std::array<double, gSIZE>& density);
    void render_frame(const glm::mat4& projectionViewMatrix);
private:
    uint32_t vaoID; // vertex array object id
    uint32_t vboID; // vertex buffer object id
    uint32_t indexID; // index buffer object id (face)
    uint32_t textureID; // texture id

};   
#endif