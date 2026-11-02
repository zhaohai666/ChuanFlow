#pragma once

// Shared internal helper for the 1D kernels: dense linear solve with
// partial pivoting (network Newton systems are small: a few hundred nodes).
#include <algorithm>
#include <cmath>
#include <vector>

#include "cf/common/Error.h"

namespace cf
{
namespace detail
{

inline void solveDenseLinear(std::vector<std::vector<double>>& a, std::vector<double>& b)
{
    const std::size_t n = b.size();
    for (std::size_t k = 0; k < n; ++k)
    {
        std::size_t piv = k;
        double best = std::abs(a[k][k]);
        for (std::size_t i = k + 1; i < n; ++i)
        {
            if (std::abs(a[i][k]) > best) { best = std::abs(a[i][k]); piv = i; }
        }
        if (best < 1e-12 * std::max(1.0, *std::max_element(
                a[piv].begin(), a[piv].end())))
        {
            cfFatal("1d: singular Newton matrix - check that every node has a"
                    " pressure reference (fixed pressure or compliance)");
        }
        if (piv != k) std::swap(a[piv], a[k]), std::swap(b[piv], b[k]);
        for (std::size_t i = k + 1; i < n; ++i)
        {
            const double m = a[i][k] / a[k][k];
            if (m == 0) continue;
            for (std::size_t j = k; j < n; ++j) a[i][j] -= m * a[k][j];
            b[i] -= m * b[k];
        }
    }
    for (std::size_t i = n; i-- > 0;)
    {
        double s = b[i];
        for (std::size_t j = i + 1; j < n; ++j) s -= a[i][j] * b[j];
        b[i] = s / a[i][i];
    }
}

} // namespace detail
} // namespace cf
