#include "Vexelization.h"

#include <cmath>

#include "mmath.h"

std::array<bool, SIZE>
generate_vexelized_sphere(int radius)
{
    std::array<bool, SIZE> sphere;

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
            sphere[POS(i, j, k)] = true;
        }
    }
    return sphere;
}

std::array<bool, SIZE>
generate_vexelized_cube(int length)
{
    std::array<bool, SIZE> cube;

    int x0 = Nx / 2 - length / 2;
    int x1 = Nx / 2 + length / 2;
    int y0 = Ny / 2 - length / 2;
    int y1 = Ny / 2 + length / 2;
    int z0 = Nz / 2 - length / 2;
    int z1 = Nz / 2 + length / 2;

    FOR_EACH_CELL
    {
        if (
            i >= x0 && i <= x1 &&
            j >= y0 && j <= y1 &&
            k >= z0 && k <= z1)
        {
            cube[POS(i, j, k)] = true;
        }
    }
    return cube;
}

// check if a ray intersect with a triangle
bool ray_triangle_intersect(
    glm::vec3 &ray_origin,
    glm::vec3 &ray_direction,
    glm::vec3 &v0, glm::vec3 &v1, glm::vec3 &v2,
    float &t)
{
    glm::vec3 edge1 = v1 - v0;
    glm::vec3 edge2 = v2 - v0;
    glm::vec3 h = cross(ray_direction, edge2);
    float a = dot(edge1, h);
    if (a > -0.00001 && a < 0.00001)
        return false;
    float f = 1.0 / a;
    glm::vec3 s = ray_origin - v0;
    float u = f * dot(s, h);
    if (u < 0.0 || u > 1.0)
        return false;
    glm::vec3 q = cross(s, edge1);
    float v = f * dot(ray_direction, q);
    if (v < 0.0 || u + v > 1.0)
        return false;
    t = f * dot(edge2, q);
    return t > 0;
}

/* general function to generate a vexelized object*/
// m_positions: the position of the vertices of the object
std::array<bool, SIZE>
generate_vexel(std::vector<glm::vec3> &m_positions, float &max_length)
{
    // find the bounding box of the object
    glm::vec3 min = m_positions[0];
    glm::vec3 max = m_positions[0];

    for (int i = 1; i < m_positions.size(); i++)
    {
        min = glm::min(min, m_positions[i]);
        max = glm::max(max, m_positions[i]);
    }
    std::cout << "m_positions.size() = " << m_positions.size() << std::endl;
    std::cout << "max = " << max.x << " " << max.y << " " << max.z << std::endl;
    std::cout << "min = " << min.x << " " << min.y << " " << min.z << std::endl;

    // x direction length
    float x_length = (int)(max.x - min.x);
    float y_length = (int)(max.y - min.y);
    float z_length = (int)(max.z - min.z);
    // max length
    max_length = std::max(x_length, std::max(y_length, z_length));
    std::cout << "max_length = " << max_length << std::endl;
    // step
    float step = max_length / Nx;
    std::cout << "step = " << step << std::endl;

    std::array<bool, SIZE> object;
    object.fill(false);

    for (int k = 0; k < Nz; k++)
    {
        for (int j = 0; j < Ny; j++)
        {
            glm::vec3 ray_origin = glm::vec3(min.x, min.y + j * step, min.z + k * step);
            glm::vec3 ray_direction = glm::vec3(1, 0, 0);
            // iterate through all the triangles
            std::vector<glm::vec3> intersections;
            for (int i = 0; i < m_positions.size(); i += 3)
            {
                glm::vec3 v0 = m_positions[i];
                glm::vec3 v1 = m_positions[i + 1];
                glm::vec3 v2 = m_positions[i + 2];
                float t;
                if (ray_triangle_intersect(ray_origin, ray_direction, v0, v1, v2, t))
                {
                    intersections.push_back(ray_origin + t * ray_direction);
                }
            }
            // sort the intersections by x value
            std::sort(
                intersections.begin(), intersections.end(),
                [](glm::vec3 a, glm::vec3 b)
                { return a.x < b.x; });
            // iterate through the intersections
            // std::cout << "intersections.size() = " << intersections.size() << std::endl;
            // if (intersections.size() > 0)
            // {
            //     std::cout << "intersections[0].x = " << intersections[0].x << std::endl;
            // }

            for (int i = 0; intersections.size() > 0 && i < intersections.size() - 1; i += 2)
            {
                int x0 = (intersections[i].x - min.x) / step;
                int x1 = (intersections[i + 1].x - min.x) / step;
                // std::cout << x0 << " " << x1 << std::endl;
                for (int x = x0; x < x1; x++)
                {
                    object[POS(x, j, k)] = true;
                }
            }
        }
    }
    return object;
}
