#include "particle/render.h"
#include "common/debug.h"
#include <glm/glm.hpp>
#include <GL/glew.h>
#include <glm/gtx/transform.hpp>

#include "common/debug.h"

SmokeRenderer::SmokeRenderer()
{
    DEBUG_PRINT("SmokeRenderer::SmokeRenderer()\n");
    // Create a obj, init the VAO, VBO
    /* domain cube */
    const float vertices[8][3] = {
        {0.0f, 0.0f, 0.0f},
        {1.0f, 0.0f, 0.0f},
        {0.0f, 1.0f, 0.0f},
        {0.0f, 0.0f, 1.0f},
        {1.0f, 1.0f, 0.0f},
        {0.0f, 1.0f, 1.0f},
        {1.0f, 0.0f, 1.0f},
        {1.0f, 1.0f, 1.0f}};

    const unsigned int indices[12][3] = {
        // z->o
        {3, 7, 5},
        {3, 6, 7},
        // x->o
        {1, 7, 6},
        {1, 4, 7},
        // o->z
        {0, 4, 1},
        {0, 2, 4},
        // o->x
        {0, 5, 2},
        {0, 3, 5},
        // y->o
        {2, 5, 7},
        {2, 7, 4},
        // o->y
        {0, 1, 6},
        {0, 6, 3}};

    // generate VAO
    DEBUG_PRINT("SmokeRenderer::SmokeRenderer() vaoID:" << vaoID);
    glGenVertexArrays(1, &vaoID);

    // set current VAO
    glBindVertexArray(vaoID);

    // generate VBO
    glGenBuffers(1, &vboID);
    glBindBuffer(GL_ARRAY_BUFFER, vboID);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(
        0,        // attribute
        3,        // size
        GL_FLOAT, // type
        GL_FALSE, // normalized?
        0,        // stride
        (void *)0 // array buffer offset
    );

    // generate index buffer
    glGenBuffers(1, &indexID);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexID);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

    glBindVertexArray(0);

    // Create a texture
    glGenTextures(1, &textureID);
    glBindTexture(GL_TEXTURE_3D, textureID);
    glBindTexture(GL_TEXTURE_3D, 0);
    DEBUG_PRINT("SmokeRenderer::SmokeRenderer() end\n");
}

SmokeRenderer::~SmokeRenderer()
{
    glDeleteTextures(1, &textureID);
}

void SmokeRenderer::render(const std::array<double, SIZE> &density)
{
    // update the texture
    // Bind the texture
    glBindTexture(GL_TEXTURE_3D, textureID);
    glBindTexture(GL_TEXTURE_3D, textureID);

    // @todo: set the texture parameters GL_CLAMP_TO_BORDER is not supported ?
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_BORDER);

    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

    GLubyte *data = new GLubyte[SIZE];
    GLubyte *ptr = data;

    for (int z = 0; z < Nz; ++z)
    {
        for (int y = 0; y < Ny; ++y)
        {
            for (int x = 0; x < Nx; ++x)
            {
                auto f = density[POS(x, y, z)];
                *ptr++ = std::max(0, std::min(255, (int)std::floor(f * 256.0)));
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

    // draw the cube
    glBindVertexArray(vaoID);

    // bind the texture
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_3D, textureID);

    // 36 : the number of indices (12 triangles * 3 vertices per triangle)
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
}

void SmokeRenderer::render_frame(const glm::mat4 &projectionViewMatrix)
{
    // render the frame
    float ratio_x = 1.f * 0.5f;
    float ratio_y = Ny / Nx * 0.5f;
    float ratio_z = Nz / Nx * 0.5f;

    float vertices[8][3] = {
        {-1.0f, -1.0f, -1.0f},
        {1.0f, -1.0f, -1.0f},
        {-1.0f, 1.0f, -1.0f},
        {-1.0f, -1.0f, 1.0f},
        {1.0f, 1.0f, -1.0f},
        {-1.0f, 1.0f, 1.0f},
        {1.0f, -1.0f, 1.0f},
        {1.0f, 1.0f, 1.0f}};

    const unsigned int indices[12][2] = {
        {0, 1}, {1, 6}, {6, 3}, {3, 0}, {0, 2}, {1, 4}, {6, 7}, {3, 5}, {2, 4}, {4, 7}, {7, 5}, {5, 2}};

    for (int i = 0; i < 12; i++)
    {
        drawLine(glm::vec3(vertices[indices[i][0]][0] * ratio_x, vertices[indices[i][0]][1] * ratio_y, vertices[indices[i][0]][2] * ratio_z),
                 glm::vec3(vertices[indices[i][1]][0] * ratio_x, vertices[indices[i][1]][1] * ratio_y, vertices[indices[i][1]][2] * ratio_z),
                 projectionViewMatrix,
                 glm::vec3(0.0f, 0.0f, 0.0f));
    }
}