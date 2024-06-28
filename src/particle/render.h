#ifndef __M_RENDER_H__
#define __M_RENDER_H__
#include <array>

#include <glm/glm.hpp>

#include "mmath.h"
#include "constants.h"
#include <GL/glew.h>
class SmokeRenderer
{

public:
    SmokeRenderer();
    ~SmokeRenderer();

    void render(std::array<double, SIZE> &density, double *transparency);
    void set_occupied_texture(std::array<bool, SIZE> &vexels);

private:
    uint32_t vaoID;                 // vertex array object id
    uint32_t vboID;                 // vertex buffer object id
    uint32_t indexID;               // index buffer object id (face)
    uint32_t textureID;             // texture id
    uint32_t occupiedTextureID;     // texture id
    uint32_t transparencyTextureID; // texture id
};

template <typename T>
void setTexture(T *t_data, uint32_t id)
{
    glBindTexture(GL_TEXTURE_3D, id);

    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_BORDER);

    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

    static GLubyte *data = new GLubyte[SIZE];
    GLubyte *ptr = data;

    for (int z = 0; z < Nz; ++z)
    {
        for (int y = 0; y < Ny; ++y)
        {
            for (int x = 0; x < Nx; ++x)
            {
                auto f = t_data[ACC3D(x, y, z, Ny, Nx)];
                *ptr++ = std::max(0, std::min(255, (int)std::floor((float)f * 256.0)));
            }
        }
    }

    glTexImage3D(GL_TEXTURE_3D,
                 0,                // mip map level, 0 means no mip map
                 GL_RED,           // internal format, single channel, 8-bit data, red
                 Nx,               // width
                 Ny,               // height
                 Nz,               // depth
                 0,                // border size
                 GL_RED,           // format of the pixel data
                 GL_UNSIGNED_BYTE, // data type of the pixel data, each pixel is a byte
                 data);

    glBindTexture(GL_TEXTURE_3D, 0);
}

#endif
