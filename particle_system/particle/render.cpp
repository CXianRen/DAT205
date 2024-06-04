#include "particle/render.h"
#include "common/debug.h"
#include <glm/glm.hpp>
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

    glGenTextures(1, &occupiedTextureID);
    glBindTexture(GL_TEXTURE_3D, occupiedTextureID);
    glBindTexture(GL_TEXTURE_3D, 0);

    // for transparency
    glGenTextures(1, &transparencyTextureID);
    glBindTexture(GL_TEXTURE_3D, transparencyTextureID);
    glBindTexture(GL_TEXTURE_3D, 0);

    DEBUG_PRINT("SmokeRenderer::SmokeRenderer() end\n");
}

SmokeRenderer::~SmokeRenderer()
{
    glDeleteTextures(1, &textureID);
}

void SmokeRenderer::set_occupied_texture(std::array<bool, SIZE> &vexels)
{
    setTexture<bool>(vexels.data(), occupiedTextureID);
}

void SmokeRenderer::render(std::array<double, SIZE> &density, double *transparency)
{
    // update the texture
    // Bind the texture
    setTexture<double>(density.data(), textureID);
    setTexture<double>(transparency, transparencyTextureID);

    // draw the cube
    glBindVertexArray(vaoID);

    // bind the texture
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_3D, textureID);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_3D, occupiedTextureID);
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_3D, transparencyTextureID);

    // 36 : the number of indices (12 triangles * 3 vertices per triangle)
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
}

void SmokeRenderer::render_frame(const glm::mat4 &projectionViewMatrix)
{
}