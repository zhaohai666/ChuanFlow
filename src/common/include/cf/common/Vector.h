#pragma once

#include <cmath>
#include <ostream>

#include "cf/common/Labels.h"

namespace cf
{

/// 3D vector; component access via .x/.y/.z or operator[].
template<typename T>
struct Vector3
{
    T x = T(0);
    T y = T(0);
    T z = T(0);

    constexpr Vector3() = default;
    constexpr Vector3(T xx, T yy, T zz) : x(xx), y(yy), z(zz) {}

    static constexpr Vector3 zero() { return {}; }

    T& operator[](int i) { return i == 0 ? x : (i == 1 ? y : z); }
    const T& operator[](int i) const { return i == 0 ? x : (i == 1 ? y : z); }

    Vector3 operator+(const Vector3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vector3 operator-(const Vector3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vector3 operator-() const { return {-x, -y, -z}; }
    Vector3 operator*(T s) const { return {x * s, y * s, z * s}; }
    Vector3 operator/(T s) const { return {x / s, y / s, z / s}; }

    Vector3& operator+=(const Vector3& o) { x += o.x; y += o.y; z += o.z; return *this; }
    Vector3& operator-=(const Vector3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
    Vector3& operator*=(T s) { x *= s; y *= s; z *= s; return *this; }

    bool operator==(const Vector3& o) const { return x == o.x && y == o.y && z == o.z; }
};

template<typename T>
Vector3<T> operator*(T s, const Vector3<T>& v) { return v * s; }

template<typename T>
T dot(const Vector3<T>& a, const Vector3<T>& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

template<typename T>
Vector3<T> cross(const Vector3<T>& a, const Vector3<T>& b)
{
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

template<typename T>
T mag2(const Vector3<T>& v) { return dot(v, v); }

template<typename T>
T mag(const Vector3<T>& v) { return std::sqrt(dot(v, v)); }

/// Dyadic product a_i b_j (used for non-orthogonal correction weights).
template<typename T>
struct Tensor3
{
    T xx{}, xy{}, xz{};
    T yx{}, yy{}, yz{};
    T zx{}, zy{}, zz{};

    static Tensor3 outer(const Vector3<T>& a, const Vector3<T>& b)
    {
        return {a.x * b.x, a.x * b.y, a.x * b.z,
                a.y * b.x, a.y * b.y, a.y * b.z,
                a.z * b.x, a.z * b.y, a.z * b.z};
    }

    Vector3<T> operator&(const Vector3<T>& v) const
    {
        return {xx * v.x + xy * v.y + xz * v.z,
                yx * v.x + yy * v.y + yz * v.z,
                zx * v.x + zy * v.y + zz * v.z};
    }

    T trace() const { return xx + yy + zz; }
};

/// Double contraction A:B of two tensors (== dot(A&v-ish sum)); returns sum_ij A_ij B_ij.
template<typename T>
T doubleContract(const Tensor3<T>& a, const Tensor3<T>& b)
{
    return a.xx * b.xx + a.xy * b.xy + a.xz * b.xz
         + a.yx * b.yx + a.yy * b.yy + a.yz * b.yz
         + a.zx * b.zx + a.zy * b.zy + a.zz * b.zz;
}

using Vector = Vector3<scalar>;
using Tensor = Tensor3<scalar>;

inline std::ostream& operator<<(std::ostream& os, const Vector& v)
{
    return os << '(' << v.x << ' ' << v.y << ' ' << v.z << ')';
}

} // namespace cf
