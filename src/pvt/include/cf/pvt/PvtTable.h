#pragma once

#include <functional>
#include <vector>

#include "cf/pvt/Components.h"

namespace cf::pvt
{

/// Tabulated PVT properties. Solvers that need a property at every cell
/// (or every time step) can pre-compute it on a (p, T) grid once and pay a
/// bilinear lookup instead of an EOS iteration; the tables are generated from
/// the same Peng-Robinson flash, so 1D and 3D stay thermodynamically
/// consistent with the direct evaluation.

class BilinearTable
{
public:
    /// xGrid / yGrid strictly ascending; values[i * yGrid.size() + j].
    BilinearTable(std::vector<scalar> xGrid, std::vector<scalar> yGrid,
                  std::vector<scalar> values);

    scalar sample(scalar x, scalar y) const;

    const std::vector<scalar>& xGrid() const { return x_; }
    const std::vector<scalar>& yGrid() const { return y_; }
    std::size_t size() const { return values_.size(); }

private:
    static std::size_t locate(const std::vector<scalar>& grid, scalar v,
                              std::size_t& lo, scalar& frac);

    std::vector<scalar> x_, y_, values_;
};

using PropertyFunction = std::function<scalar(scalar p, scalar T)>;

/// Evaluate f over the (p, T) grid. Pressures in Pa, temperatures in K.
BilinearTable buildTable(const std::vector<scalar>& pGrid,
                         const std::vector<scalar>& tGrid,
                         const PropertyFunction& f);

/// Phase densities (kg/m3) and vapour fraction (mole basis) of a feed fluid.
BilinearTable liquidDensityTable(const Mixture& feed,
                                 const std::vector<scalar>& pGrid,
                                 const std::vector<scalar>& tGrid);
BilinearTable vapourDensityTable(const Mixture& feed,
                                 const std::vector<scalar>& pGrid,
                                 const std::vector<scalar>& tGrid);
BilinearTable vapourFractionTable(const Mixture& feed,
                                  const std::vector<scalar>& pGrid,
                                  const std::vector<scalar>& tGrid);

} // namespace cf::pvt
