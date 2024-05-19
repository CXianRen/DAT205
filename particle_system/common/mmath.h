#ifndef __MMATH_H__
#define __MMATH_H__

#include <vector>
#include "Vec3.h"

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
    T *src_list,
    const std::vector<T *> &dst_list,
    const std::vector<int> &dims, T h)
{
    T p, q;
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

// interpolate 3D scalar field
template <typename T>
T inline linearInterpolation(const Vec3 &pt, T *src, const std::vector<int> &dims, T h)
{
    Vec3 pos;
    // clamp position
    pos[0] = std::min(std::max((T)0.0, pt[0]), (T)(dims[0] - 1) * h - (T)1e-6);
    pos[1] = std::min(std::max((T)0.0, pt[1]), (T)(dims[1] - 1) * h - (T)1e-6);
    pos[2] = std::min(std::max((T)0.0, pt[2]), (T)(dims[2] - 1) * h - (T)1e-6);

    int i = (int)(pos[0] / h);
    int j = (int)(pos[1] / h);
    int k = (int)(pos[2] / h);

    T scale = 1.0 / h;
    T fractx = scale * (pos[0] - i * h);
    T fracty = scale * (pos[1] - j * h);
    T fractz = scale * (pos[2] - k * h);

    assert(fractx < 1.0 && fractx >= 0);
    assert(fracty < 1.0 && fracty >= 0);
    assert(fractz < 1.0 && fractz >= 0);

    // Y @ low X, low Z:
    T tmp1 = src[ACCESS3D(i, j, k)];
    T tmp2 = src[ACCESS3D(i, j + 1, k)];
    // Y @ high X, low Z:
    T tmp3 = src[ACCESS3D(i + 1, j, k)];
    T tmp4 = src[ACCESS3D(i + 1, j + 1, k)];

    // Y @ low X, high Z:
    T tmp5 = src[ACCESS3D(i, j, k + 1)];
    T tmp6 = src[ACCESS3D(i, j + 1, k + 1)];
    // Y @ high X, high Z:
    T tmp7 = src[ACCESS3D(i + 1, j, k + 1)];
    T tmp8 = src[ACCESS3D(i + 1, j + 1, k + 1)];

    // Y @ low X, low Z
    T tmp12 = ((T)(1) - fracty) * tmp1 + fracty * tmp2;
    // Y @ high X, low Z
    T tmp34 = ((T)(1) - fracty) * tmp3 + fracty * tmp4;

    // Y @ low X, high Z
    T tmp56 = ((T)(1) - fracty) * tmp5 + fracty * tmp6;
    // Y @ high X, high Z
    T tmp78 = ((T)(1) - fracty) * tmp7 + fracty * tmp8;

    // X @ low Z
    T tmp1234 = ((T)(1) - fractx) * tmp12 + fractx * tmp34;
    // X @ high Z
    T tmp5678 = ((T)(1) - fractx) * tmp56 + fractx * tmp78;

    // Z
    T tmp = ((T)(1) - fractz) * tmp1234 + fractz * tmp5678;
    return tmp;
}

// monotonicCubicInterpolation
int inline constrainIndex(const int idx, const int N)
{
    if (idx == -1)
    {
        return 0;
    }
    if (idx > N)
    {
        return N;
    }
    return idx;
}

template <typename T>
int inline sign(const T a)
{
    // if a is positive, return 1
    // if a is negative, return -1
    // if a is zero, return 0
    return (a > (T)(0)) - (a < (T)(0));
}

template <typename T> 
T inline axis_monotonicCubicInterpolation(const T *f, const T t)
{
    // f[0]:f_k-1, f[1]:f_k, f[2]:f_k+1, f[3]:f_k+2
    T delta = f[2] - f[1];
    T d0 = (T)(0.5) * (f[2] - f[0]);
    T d1 = (T)(0.5) * (f[3] - f[1]);

    // neccessary condition for monotonic
    d0 = (T)sign(delta) * std::abs(d0);
    d1 = (T)sign(delta) * std::abs(d1);

    T a0 = f[1];
    T a1 = d0;
    T a2 = (T)(3) * delta - T(2) * d0 - d1;
    T a3 = d0 + d1 - T(2) * delta;
    return a3 * t * t * t + a2 * t * t + a1 * t + a0;
}

template <typename T>
T inline monotonicCubicInterpolation(const Vec3 &pt, T *src, const std::vector<int> &dims, T h)
{
    Vec3 pos;
    // clamp position
    pos[0] = std::min(std::max((T)(0.0), pt[0]), (T)(dims[0]-1) * h - (T)(1e-6));
    pos[1] = std::min(std::max((T)(0.0), pt[1]), (T)(dims[1]-1) * h - (T)(1e-6));
    pos[2] = std::min(std::max((T)(0.0), pt[2]), (T)(dims[2]-1) * h - (T)(1e-6));

    int i = (int)(pos[0] / h);
    int j = (int)(pos[1] / h);
    int k = (int)(pos[2] / h);

    T scale = (T)(1.0) / h;
    T fractx = scale * (pos[0] - i * h);
    T fracty = scale * (pos[1] - j * h);
    T fractz = scale * (pos[2] - k * h);

    assert(fractx < 1.0 && fractx >= 0);
    assert(fracty < 1.0 && fracty >= 0);
    assert(fractz < 1.0 && fractz >= 0);

    // Z:
    T arr_z[4];
    for (int z = 0; z < 4; ++z)
    {
        // X @ Z_(z-1):
        T arr_x[4];
        for (int x = 0; x < 4; ++x)
        {
            // Y @ X_(x-1), Z_(z-1):
            int i1 = constrainIndex(i + x - 1, dims[0]-1);
            int j1 = constrainIndex(j - 1, dims[1]-1);
            int j2 = constrainIndex(j + 1, dims[1]-1);
            int j3 = constrainIndex(j + 2, dims[1]-1);
            int k1 = constrainIndex(k + z - 1, dims[2]-1);

            T arr_y[4] = {src[ACCESS3D(i1, j1, k1)], src[ACCESS3D(i1, j, k1)], src[ACCESS3D(i1, j2, k1)], src[ACCESS3D(i1, j3, k1)]};
            arr_x[x] = axis_monotonicCubicInterpolation(arr_y, fracty);
        }
        arr_z[z] = axis_monotonicCubicInterpolation(arr_x, fractx);
    }

    return axis_monotonicCubicInterpolation(arr_z, fractz);
}

#endif // __MMATH_H__