#include <algorithm>
#include <cmath>
#include <cstddef>

#include "cf/common/Error.h"
#include "cf/parallel/DomainDecomposition.h"

namespace cf
{
namespace parallel
{
namespace
{

/// Local slots of the global owner and neighbour behind a stencil face. The
/// domain stores the cell it updates in `cell`, so the owner is either that
/// cell or its peer.
void ownerNeighbour(const DomainFace& df, label& ownerSlot, label& nbrSlot)
{
    if (df.sign > 0)
    {
        ownerSlot = df.cell;
        nbrSlot = df.peer;
    }
    else
    {
        ownerSlot = df.peer;
        nbrSlot = df.cell;
    }
}

scalar faceValueAt(const DomainMesh& dm, std::size_t k, const Field<scalar>& localT,
                   const FaceValue& faceValue)
{
    const DomainFace& df = dm.faces()[k];
    label o = -1, n = -1;
    ownerNeighbour(df, o, n);
    if (n < 0) return faceValue(df.face, dm.faceCentres()[k]);
    const scalar w = dm.faceOwnerWeight()[k];
    return w * localT[static_cast<std::size_t>(o)] + (1.0 - w) * localT[static_cast<std::size_t>(n)];
}

Vector faceGradAt(const DomainMesh& dm, std::size_t k, const Field<Vector>& localGrad)
{
    const DomainFace& df = dm.faces()[k];
    label o = -1, n = -1;
    ownerNeighbour(df, o, n);
    if (n < 0) return localGrad[static_cast<std::size_t>(df.cell)];
    const scalar w = dm.faceOwnerWeight()[k];
    return localGrad[static_cast<std::size_t>(o)] * w
         + localGrad[static_cast<std::size_t>(n)] * (1.0 - w);
}

/// V * grad accumulated over the faces each domain sees. The serial Gauss sum
/// adds +flux to the owner and -flux to the neighbour; a domain holds one of
/// those sides, so it applies the matching sign of the same flux -- and both
/// sides compute that flux from identical inputs, so they agree exactly.
void accumulateFaceSums(const DomainSet& ds,
                        const std::vector<Field<scalar>>& localT,
                        const std::vector<Field<Vector>>& localGrad,
                        const FaceValue& faceValue,
                        bool corrected,
                        std::vector<Field<Vector>>& sum)
{
    for (label d = 0; d < ds.size(); ++d)
    {
        const DomainMesh& dm = ds[d];
        Field<Vector>& dst = sum[static_cast<std::size_t>(d)];
        const Field<scalar>& lt = localT[static_cast<std::size_t>(d)];
        const Field<Vector>& lg = localGrad[static_cast<std::size_t>(d)];

        for (std::size_t k = 0; k < static_cast<std::size_t>(dm.nFaces()); ++k)
        {
            const Vector sfv = dm.faceSf()[k];
            const scalar tf = faceValueAt(dm, k, lt, faceValue);
            Vector flux = sfv * tf;
            if (corrected)
            {
                const Vector dvec = dm.faceDelta()[k];
                const scalar dmag = mag(dvec);
                if (dmag <= 0) cfFatal("assembledGrad: degenerate face-cell distance");
                const Vector dh = dvec / dmag;
                const Vector cfVec = sfv - dh * dot(sfv, dh);
                flux = sfv * tf - cfVec * dot(faceGradAt(dm, k, lg), dh);
            }
            const DomainFace& df = dm.faces()[k];
            if (df.sign > 0) dst[static_cast<std::size_t>(df.cell)] += flux;
            else dst[static_cast<std::size_t>(df.cell)] -= flux;
        }
    }
}

/// sum -> sum/volume, owned slots only: halo cells carry values, never updates.
void divideByVolume(const DomainSet& ds, const std::vector<Field<Vector>>& sum,
                    std::vector<Field<Vector>>& grad)
{
    for (label d = 0; d < ds.size(); ++d)
    {
        const DomainMesh& dm = ds[d];
        for (std::size_t i = 0; i < static_cast<std::size_t>(dm.nOwned()); ++i)
        {
            grad[static_cast<std::size_t>(d)][i] = sum[static_cast<std::size_t>(d)][i] / dm.cellVolumes()[i];
        }
    }
}

} // namespace

Field<Vector> assembledGrad(const DomainSet& ds, const Field<scalar>& cellValues,
                            const FaceValue& boundaryValueOf, const std::string& scheme)
{
    const Decomposition& dec = ds.decomposition();

    std::vector<Field<scalar>> localT = ds.storage<scalar>();
    ds.scatter(cellValues, localT);
    exchange(dec, localT);

    std::vector<Field<Vector>> grad = ds.storage<Vector>();
    std::vector<Field<Vector>> sum = ds.storage<Vector>();

    accumulateFaceSums(ds, localT, grad, boundaryValueOf, false, sum);
    divideByVolume(ds, sum, grad);

    // The correction interpolates the neighbour's gradient, so it is exchanged
    // between passes: the same collective a distributed run issues there.
    if (scheme == "corrected")
    {
        for (int pass = 0; pass < 2; ++pass)
        {
            exchange(dec, grad);
            for (auto& s : sum) std::fill(s.begin(), s.end(), Vector::zero());
            accumulateFaceSums(ds, localT, grad, boundaryValueOf, true, sum);
            divideByVolume(ds, sum, grad);
        }
    }

    return ds.gather(grad);
}

Field<scalar> assembledDivPhi(const DomainSet& ds, const FaceFlux& facePhi)
{
    std::vector<Field<scalar>> sum = ds.storage<scalar>();
    for (label d = 0; d < ds.size(); ++d)
    {
        const DomainMesh& dm = ds[d];
        Field<scalar>& dst = sum[static_cast<std::size_t>(d)];
        for (std::size_t k = 0; k < static_cast<std::size_t>(dm.nFaces()); ++k)
        {
            const DomainFace& df = dm.faces()[k];
            const scalar p = facePhi(df.face, dm.faceCentres()[k], dm.faceSf()[k]);
            if (df.sign > 0) dst[static_cast<std::size_t>(df.cell)] += p;
            else dst[static_cast<std::size_t>(df.cell)] -= p;
        }
    }
    return ds.gather(sum);
}

} // namespace parallel
} // namespace cf
