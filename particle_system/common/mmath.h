#ifndef __MMATH_H__
#define __MMATH_H__

#include <vector>

#define ACCESS3D(x, y, z) ((x) + (y) * dims[0] + (z) * dims[0] * dims[1])
#define ACCESS3D_X(x, y, z) ((x) + (y) * (dims[0] + 1) + (z) * (dims[0] + 1) * dims[1])
#define ACCESS3D_Y(x, y, z) ((x) + (y) * dims[0] + (z) * dims[0] * (dims[1] + 1))
#define ACCESS3D_Z(x, y, z) ((x) + (y) * dims[0] + (z) * dims[0] * dims[1])

template <typename T>
void inline calculate_3D_Field_average(
    const std::vector<T *> &src_list,
    const std::vector<T *> &dst_list,
    const std::vector<int> &dims)
{
    for (int k = 0; k < dims[2]; ++k)
        for (int j = 0; j < dims[1]; ++j)
            for (int i = 0; i < dims[0]; ++i)
            {
                // x component
                dst_list[0][ACCESS3D(i, j, k)] =
                    (src_list[0][ACCESS3D_X(i, j, k)] + src_list[0][ACCESS3D_X(i + 1, j, k)]) * 0.5;
                // y component
                dst_list[1][ACCESS3D(i, j, k)] =
                    (src_list[1][ACCESS3D_Y(i, j, k)] + src_list[1][ACCESS3D_Y(i, j + 1, k)]) * 0.5;
                // z component
                dst_list[2][ACCESS3D(i, j, k)] =
                    (src_list[2][ACCESS3D_Z(i, j, k)] + src_list[2][ACCESS3D_Z(i, j, k + 1)]) * 0.5;
            }
}

template <typename T>
void inline calculate_3D_Field_vorticity(
    const std::vector<T *> &src_list,
    const std::vector<T *> &dst_list,
    const std::vector<int> &dims, T h)
{
    for (int k = 1; k < dims[2] - 1; ++k)
        for (int j = 1; j < dims[1] - 1; ++j)
            for (int i = 1; i < dims[0] - 1; ++i)
            {
                dst_list[0][ACCESS3D(i, j, k)] =
                    (src_list[2][ACCESS3D(i, j + 1, k)] - src_list[2][ACCESS3D(i, j - 1, k)] - (src_list[1][ACCESS3D(i, j, k + 1)] - src_list[1][ACCESS3D(i, j, k - 1)])) * 0.5 / h;

                dst_list[1][ACCESS3D(i, j, k)] =
                    (src_list[0][ACCESS3D(i, j, k + 1)] - src_list[0][ACCESS3D(i, j, k - 1)] - (src_list[2][ACCESS3D(i + 1, j, k)] - src_list[2][ACCESS3D(i - 1, j, k)])) * 0.5 / h;

                dst_list[2][ACCESS3D(i, j, k)] =
                    (src_list[1][ACCESS3D(i + 1, j, k)] - src_list[1][ACCESS3D(i - 1, j, k)] - (src_list[0][ACCESS3D(i, j + 1, k)] - src_list[0][ACCESS3D(i, j - 1, k)])) * 0.5 / h;
            }
}

template <typename T>
void inline calculate_Scalar_Field_gradient(
    T* src_list,
    const std::vector<T *> &dst_list,
    const std::vector<int> &dims, T h)
{
    double p, q;
    for (int k = 1; k < dims[2] - 1; ++k)
        for (int j = 1; j < dims[1] - 1; ++j)
            for (int i = 1; i < dims[0] - 1; ++i)
            {

                q = src_list[ACCESS3D(i, j, k)];
                p = src_list[ACCESS3D(i + 1, j, k)];
                dst_list[0][ACCESS3D(i, j, k)] = (p - q) * 0.5 / h;

                p = src_list[ACCESS3D(i, j + 1, k)];
                dst_list[1][ACCESS3D(i, j, k)] = (p - q) * 0.5 / h;

                p = src_list[ACCESS3D(i, j, k + 1)];
                dst_list[1][ACCESS3D(i, j, k)] = (p - q) * 0.5 / h;
            }
}

#endif // __MMATH_H__