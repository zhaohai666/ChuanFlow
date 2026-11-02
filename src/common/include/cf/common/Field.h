#pragma once

#include <algorithm>
#include <numeric>
#include <span>
#include <vector>

#include "cf/common/Vector.h"

namespace cf
{

/// Contiguous value field; thin alias over std::vector with numeric helpers.
template<typename T>
class Field : public std::vector<T>
{
public:
    using base = std::vector<T>;
    using base::base;

    Field(std::size_t n, T v = T{}) : base(n, v) {}

    std::span<T> span_() { return {this->data(), this->size()}; }
};

using ScalarField = Field<scalar>;

template<typename T>
T sum(const Field<T>& f)
{
    return std::accumulate(f.begin(), f.end(), T{});
}

/// L2 norm of the field values (vector components summed).
inline scalar l2Norm(const Field<scalar>& f)
{
    return std::sqrt(std::inner_product(f.begin(), f.end(), f.begin(), 0.0));
}

/// Max abs value.
inline scalar lInfNorm(const Field<scalar>& f)
{
    scalar m = 0;
    for (const auto v : f) m = std::max(m, std::abs(v));
    return m;
}

} // namespace cf
