#include "cf/fvm/Reconstruction.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "cf/common/Error.h"

namespace cf
{
namespace reconstruction
{

Limiter parseLimiter(const std::string& name)
{
    if (name == "minmod") return Limiter::Minmod;
    if (name == "vanLeer") return Limiter::VanLeer;
    if (name == "MC" || name == "mc") return Limiter::MC;
    if (name == "superbee") return Limiter::Superbee;
    cfFatal("reconstruction: unknown limiter '" + name + "'");
    return Limiter::Minmod; // unreachable; keeps -Wreturn-type quiet
}

const char* limiterName(Limiter limiter)
{
    switch (limiter)
    {
        case Limiter::Minmod:   return "minmod";
        case Limiter::VanLeer:  return "vanLeer";
        case Limiter::MC:       return "MC";
        case Limiter::Superbee: return "superbee";
    }
    return "unknown";
}

scalar limiterPhi(scalar r, Limiter limiter)
{
    if (r <= 0.0) return 0.0; // opposite-signed gradients: no extrapolation
    switch (limiter)
    {
        case Limiter::Minmod:
            return std::min(r, 1.0);
        case Limiter::VanLeer:
            return 2.0 * r / (1.0 + r);
        case Limiter::MC:
            return std::min(std::min(2.0 * r, 0.5 * (1.0 + r)), 2.0);
        case Limiter::Superbee:
            return std::max(std::min(2.0 * r, 1.0), std::min(r, 2.0));
    }
    return 0.0;
}

scalar limitedSlope(scalar a, scalar b, scalar c, Limiter limiter)
{
    const scalar s1 = b - a; // upwind-biased slope ~ gP*dx
    if (s1 == 0.0) return 0.0;
    const scalar s2 = c - b; // face slope ~ gC*L
    const scalar r = s2 / s1; // audit's rP = gC/gP
    return limiterPhi(r, limiter) * s1; // audit's phi(rP)*gP, scaled by dx
}

Field<scalar> upwindFaceValues(const polyMesh& mesh, const Field<scalar>& cellData,
                               const Field<scalar>& phi)
{
    if (cellData.size() < static_cast<std::size_t>(mesh.nCells()) ||
        phi.size() < static_cast<std::size_t>(mesh.nFaces()))
        cfFatal("reconstruction: field/flux size does not match the mesh");

    Field<scalar> faceData(static_cast<std::size_t>(mesh.nFaces()), 0.0);
    const Field<label>& own = mesh.owner();
    const Field<label>& nb = mesh.neighbour();
    for (label f = 0; f < mesh.nFaces(); ++f)
    {
        const std::size_t u = static_cast<std::size_t>(f);
        const label o = own[u];
        if (f >= mesh.nInternalFaces() || phi[u] >= 0.0)
            faceData[u] = cellData[static_cast<std::size_t>(o)];
        else
            faceData[u] = cellData[static_cast<std::size_t>(nb[u])];
    }
    return faceData;
}

Field<scalar> musclFaceValues(const polyMesh& mesh, const Field<scalar>& cellData,
                              const Field<scalar>& phi, Limiter limiter)
{
    if (cellData.size() < static_cast<std::size_t>(mesh.nCells()) ||
        phi.size() < static_cast<std::size_t>(mesh.nFaces()))
        cfFatal("reconstruction: field/flux size does not match the mesh");

    const label nInternal = mesh.nInternalFaces();
    const Field<Vector>& centre = mesh.cellCentres();
    const Field<Vector>& fcentre = mesh.faceCentres();
    // Audit: use the raw owner->neighbour centre vector d = cN - cO. Never the
    // clamped makeInternalOwnerWeights, which would break linear-field exactness.
    const Field<Vector>& delta = mesh.deltaCoeffVec(); // d for internal faces
    const Field<label>& own = mesh.owner();
    const Field<label>& nb = mesh.neighbour();
    const std::vector<std::vector<label>>& cellFaces = mesh.cells();

    Field<scalar> faceData(static_cast<std::size_t>(mesh.nFaces()), 0.0);
    for (label f = 0; f < mesh.nFaces(); ++f)
    {
        const std::size_t u = static_cast<std::size_t>(f);
        const label o = own[u];

        if (f >= nInternal)
        {
            faceData[u] = cellData[static_cast<std::size_t>(o)]; // zeroGradient
            continue;
        }

        const label n = nb[u];
        const label up   = phi[u] >= 0.0 ? o : n;
        const label down = phi[u] >= 0.0 ? n : o;
        const Vector& cu = centre[static_cast<std::size_t>(up)];
        const Vector& d = delta[u]; // d = cN - cO, L = |d|, e = d/L
        const scalar dUD = mag(d);
        const Vector e = d / dUD;
        // flow direction along the centre line (upwind -> downwind)
        const Vector flowE = phi[u] >= 0.0 ? e : -e;

        // The stencil needs the cell just upstream of the upwind cell. Scan
        // the upwind cell's own faces and keep the neighbour that lies farthest
        // upstream along the owner->neighbour centre line.
        label far = -1;
        scalar best = -std::numeric_limits<scalar>::max();
        for (const label g : cellFaces[static_cast<std::size_t>(up)])
        {
            if (g == f) continue;
            const std::size_t v = static_cast<std::size_t>(g);
            label other = -1;
            if (own[v] == up)
            {
                if (g < nInternal) other = nb[v]; // boundary face has no cell
            }
            else if (nb[v] == up)
            {
                other = own[v];
            }
            if (other < 0) continue;
            const scalar score = dot(cu - centre[static_cast<std::size_t>(other)], flowE);
            // the neighbour must actually lie upstream (a transverse neighbour
            // scores ~0 and must not be mistaken for an upstream cell)
            if (score > 1e-12 && score > best)
            {
                best = score;
                far = other;
            }
        }

        const scalar Tu = cellData[static_cast<std::size_t>(up)];
        const scalar Td = cellData[static_cast<std::size_t>(down)];
        scalar value = Tu;
        if (far >= 0)
        {
            const scalar L = limitedSlope(cellData[static_cast<std::size_t>(far)], Tu, Td,
                                          limiter);
            const scalar r = mag(fcentre[u] - cu) / dUD; // 0.5 on a uniform mesh
            value = Tu + L * r;
        }
        else
        {
            // No upstream cell (first cell along the flow): linear ghost cell
            // T_ghost = 2*Tu - Td makes s1 == s2, which keeps the slope exact
            // on linear fields and first order at the domain inlet.
            const scalar L = limitedSlope(2.0 * Tu - Td, Tu, Td, limiter);
            const scalar r = mag(fcentre[u] - cu) / dUD;
            value = Tu + L * r;
        }
        // Monotonicity clamp: a face value may never leave the range of the
        // two cells that share the face (TVD on discontinuities).
        value = std::min(std::max(value, std::min(Tu, Td)), std::max(Tu, Td));
        faceData[u] = value;
    }
    return faceData;
}

} // namespace reconstruction
} // namespace cf
