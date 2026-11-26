#pragma once

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

    std::span<T> span_() { return {this->data(), this->size()}; }
};

using ScalarField = Field<scalar>;

template<typename T>
T sum(const Field<T>& f)
{
    return std::accumulate(f.begin(), f.end(), T{});
}

} // namespace cf
