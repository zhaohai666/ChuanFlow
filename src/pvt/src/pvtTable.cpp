#include <algorithm>
#include <cmath>

#include "cf/common/Error.h"
#include "cf/pvt/Flash.h"
#include "cf/pvt/PvtTable.h"

namespace cf::pvt
{

BilinearTable::BilinearTable(std::vector<scalar> xGrid, std::vector<scalar> yGrid,
                             std::vector<scalar> values)
    : x_(std::move(xGrid)), y_(std::move(yGrid)), values_(std::move(values))
{
    if (x_.size() < 2 || y_.size() < 2)
        cfFatal("pvt: a bilinear table needs at least 2 points per axis");
    if (values_.size() != x_.size() * y_.size())
        cfFatal("pvt: table has " + std::to_string(values_.size())
                + " values for a " + std::to_string(x_.size()) + " x "
                + std::to_string(y_.size()) + " grid");
    for (std::size_t i = 1; i < x_.size(); ++i)
        if (!(x_[i] > x_[i - 1])) cfFatal("pvt: table x grid must be strictly ascending");
    for (std::size_t j = 1; j < y_.size(); ++j)
        if (!(y_[j] > y_[j - 1])) cfFatal("pvt: table y grid must be strictly ascending");
}

std::size_t BilinearTable::locate(const std::vector<scalar>& grid, scalar v,
                                  std::size_t& lo, scalar& frac)
{
    if (v <= grid.front())
    {
        lo = 0;
        frac = 0.0;
        return 1;
    }
    if (v >= grid.back())
    {
        lo = grid.size() - 2;
        frac = 1.0;
        return grid.size() - 1;
    }
    const auto it = std::upper_bound(grid.begin(), grid.end(), v);
    const std::size_t hi = static_cast<std::size_t>(it - grid.begin());
    lo = hi - 1;
    frac = (v - grid[lo]) / (grid[hi] - grid[lo]);
    return hi;
}

scalar BilinearTable::sample(scalar x, scalar y) const
{
    std::size_t i0 = 0, i1 = 0, j0 = 0, j1 = 0;
    scalar fx = 0.0, fy = 0.0;
    i1 = locate(x_, x, i0, fx);
    j1 = locate(y_, y, j0, fy);
    const std::size_t ny = y_.size();
    const scalar v00 = values_[i0 * ny + j0];
    const scalar v01 = values_[i0 * ny + j1];
    const scalar v10 = values_[i1 * ny + j0];
    const scalar v11 = values_[i1 * ny + j1];
    const scalar a = v00 + fx * (v10 - v00);
    const scalar b = v01 + fx * (v11 - v01);
    return a + fy * (b - a);
}

BilinearTable buildTable(const std::vector<scalar>& pGrid,
                         const std::vector<scalar>& tGrid,
                         const PropertyFunction& f)
{
    std::vector<scalar> values(pGrid.size() * tGrid.size());
    for (std::size_t i = 0; i < pGrid.size(); ++i)
        for (std::size_t j = 0; j < tGrid.size(); ++j)
            values[i * tGrid.size() + j] = f(pGrid[i], tGrid[j]);
    return BilinearTable(pGrid, tGrid, std::move(values));
}

namespace
{

FlashResult flashAt(const Mixture& feed, scalar p, scalar T)
{
    return tpFlash(feed, T, p);
}

} // namespace

BilinearTable liquidDensityTable(const Mixture& feed,
                                 const std::vector<scalar>& pGrid,
                                 const std::vector<scalar>& tGrid)
{
    return buildTable(pGrid, tGrid, [&](scalar p, scalar T)
                      { return flashAt(feed, p, T).liquidMassDensity; });
}

BilinearTable vapourDensityTable(const Mixture& feed,
                                 const std::vector<scalar>& pGrid,
                                 const std::vector<scalar>& tGrid)
{
    return buildTable(pGrid, tGrid, [&](scalar p, scalar T)
                      { return flashAt(feed, p, T).vapourMassDensity; });
}

BilinearTable vapourFractionTable(const Mixture& feed,
                                  const std::vector<scalar>& pGrid,
                                  const std::vector<scalar>& tGrid)
{
    return buildTable(pGrid, tGrid, [&](scalar p, scalar T)
                      { return flashAt(feed, p, T).vaporFraction; });
}

} // namespace cf::pvt
