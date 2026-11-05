#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>

#include "cf/common/Error.h"
#include "cf/fvm/FvMatrix.h"
#include "cf/parallel/DomainDecomposition.h"

namespace cf
{
namespace parallel
{
namespace
{

/// Sort the block by the longest axis of its cell-centre bounding box and cut
/// it at the median, then order each half the same way. Contiguous slices of the
/// result are therefore compact in space, and so are the faces inside one slice.
///
/// Cells in one layer along the sort axis have the same coordinate in exact
/// arithmetic, but a divergence-theorem centroid carries ~1e-16 of rounding
/// noise, and without a tie band the median cut of a layer would be chosen by
/// that noise: the partition would then not be reproducible from the geometry.
/// So equal coordinates are decided by cell id, which is what the checker
/// mirrors.
void orderBlock(const Field<Vector>& centres, std::vector<label>& ids,
                std::size_t from, std::size_t to)
{
    const std::size_t n = to - from;
    if (n <= 1) return;

    Vector lo = centres[static_cast<std::size_t>(ids[from])];
    Vector hi = lo;
    for (std::size_t i = from + 1; i < to; ++i)
    {
        const Vector& c = centres[static_cast<std::size_t>(ids[i])];
        lo.x = std::min(lo.x, c.x); lo.y = std::min(lo.y, c.y); lo.z = std::min(lo.z, c.z);
        hi.x = std::max(hi.x, c.x); hi.y = std::max(hi.y, c.y); hi.z = std::max(hi.z, c.z);
    }
    const Vector e = hi - lo;
    // A symmetric block has mathematically equal extents, and a centroid's
    // rounding noise would then decide the cut direction of the whole
    // decomposition. So the longest axis is chosen within a band, lowest index
    // first: x, then y, then z.
    const scalar emax = std::max(e.x, std::max(e.y, e.z));
    const scalar band = 1e-9 * emax;
    int axis = 0;
    if (e.x < emax - band) axis = (e.y < emax - band) ? 2 : 1;

    const auto coordinate = [&](label c) { return centres[static_cast<std::size_t>(c)][axis]; };
    std::sort(ids.begin() + static_cast<std::ptrdiff_t>(from),
              ids.begin() + static_cast<std::ptrdiff_t>(to),
              [&](label a, label b)
              {
                  const scalar ka = coordinate(a), kb = coordinate(b);
                  if (ka != kb) return ka < kb;
                  return a < b;
              });

    // collapse each noise-tied run and order it by cell id
    const scalar tie = 1e-9 * std::max(e[axis], scalar(1e-30));
    std::size_t g = from;
    while (g < to)
    {
        const scalar anchor = coordinate(ids[g]);
        std::size_t h = g + 1;
        while (h < to && std::abs(coordinate(ids[h]) - anchor) <= tie) ++h;
        std::sort(ids.begin() + static_cast<std::ptrdiff_t>(g),
                  ids.begin() + static_cast<std::ptrdiff_t>(h));
        g = h;
    }

    const std::size_t mid = from + n / 2;
    orderBlock(centres, ids, from, mid);
    orderBlock(centres, ids, mid, to);
}

Decomposition::Message* messageFor(std::vector<Decomposition::Message>& plan, label peer)
{
    for (auto& m : plan)
    {
        if (m.peer == peer) return &m;
    }
    plan.push_back(Decomposition::Message{peer, {}});
    return &plan.back();
}

} // namespace

Decomposition::Decomposition(const polyMesh& mesh, label nDomains)
    : mesh_(mesh), nDomains_(nDomains)
{
    const label nC = mesh.nCells();
    if (nDomains <= 0)
    {
        cfFatal("decomposition: " + std::to_string(nDomains) + " domains is not a decomposition");
    }
    if (nDomains > nC)
    {
        cfFatal("decomposition: " + std::to_string(nDomains) + " domains for "
                + std::to_string(nC) + " cells; an empty domain owns no face and its stencil degenerates");
    }

    std::vector<label> order(static_cast<std::size_t>(nC));
    std::iota(order.begin(), order.end(), 0);
    orderBlock(mesh.cellCentres(), order, 0, order.size());

    cells_.resize(static_cast<std::size_t>(nDomains));
    haloCells_.resize(static_cast<std::size_t>(nDomains));
    faces_.resize(static_cast<std::size_t>(nDomains));
    send_.resize(static_cast<std::size_t>(nDomains));
    recv_.resize(static_cast<std::size_t>(nDomains));
    slot_.resize(static_cast<std::size_t>(nDomains));
    cellDomain_.assign(static_cast<std::size_t>(nC), -1);

    for (label i = 0; i < nC; ++i)
    {
        const label d = (i * nDomains) / nC;
        const label c = order[static_cast<std::size_t>(i)];
        cells_[static_cast<std::size_t>(d)].push_back(c);
        slot_[static_cast<std::size_t>(d)][c] = static_cast<label>(cells_[static_cast<std::size_t>(d)].size()) - 1;
        cellDomain_[static_cast<std::size_t>(c)] = d;
    }

    const Field<label>& owner = mesh.owner();
    const Field<label>& nbr = mesh.neighbour();
    const label nInt = mesh.nInternalFaces();

    // Sweep the global faces in ascending order: a cell then receives its face
    // contributions in the order the serial sweep gives them, which is what lets
    // a decomposed stencil reproduce the serial one to the last bit.
    for (label f = 0; f < mesh.nFaces(); ++f)
    {
        const std::size_t u = static_cast<std::size_t>(f);
        const label o = owner[u];
        const label n = (f < nInt) ? nbr[u] : -1;
        const label od = cellDomain_[static_cast<std::size_t>(o)];

        addSide(od, f, o, n, 1.0);
        if (n >= 0)
        {
            const label nd = cellDomain_[static_cast<std::size_t>(n)];
            if (nd == od)
            {
                addSide(od, f, n, o, -1.0);
            }
            else
            {
                addSide(nd, f, n, o, -1.0);
                ++nCoupled_;

                // Both sides need the cell on the other side, so the face
                // produces two messages. Each is appended to the sender's send
                // list and the receiver's matching receive list in this same
                // iteration, which keeps the two index-matched.
                messageFor(send_[static_cast<std::size_t>(od)], nd)->slots.push_back(localSlot(od, o));
                messageFor(recv_[static_cast<std::size_t>(nd)], od)->slots.push_back(localSlot(nd, o));
                messageFor(send_[static_cast<std::size_t>(nd)], od)->slots.push_back(localSlot(nd, n));
                messageFor(recv_[static_cast<std::size_t>(od)], nd)->slots.push_back(localSlot(od, n));
            }
        }
    }
}

void Decomposition::addSide(label d, label face, label cell, label peer, scalar sign)
{
    DomainFace df;
    df.face = face;
    df.cell = slotOf(d, cell);
    df.peer = (peer >= 0) ? slotOf(d, peer) : -1;
    df.sign = sign;
    faces_[static_cast<std::size_t>(d)].push_back(df);
}

label Decomposition::slotOf(label d, label cell)
{
    auto& m = slot_[static_cast<std::size_t>(d)];
    const auto it = m.find(cell);
    if (it != m.end()) return it->second;

    const label slot = nCells(d) + static_cast<label>(haloCells_[static_cast<std::size_t>(d)].size());
    haloCells_[static_cast<std::size_t>(d)].push_back(cell);
    m[cell] = slot;
    return slot;
}

label Decomposition::localSlot(label d, label cell) const
{
    const auto& m = slot_[static_cast<std::size_t>(d)];
    const auto it = m.find(cell);
    return it == m.end() ? -1 : it->second;
}

const Decomposition::Message* Decomposition::recvFrom(label d, label peer) const
{
    for (const Message& m : recv_[static_cast<std::size_t>(d)])
    {
        if (m.peer == peer) return &m;
    }
    return nullptr;
}

label Decomposition::nMessages() const
{
    label n = 0;
    for (const auto& plan : send_) n += static_cast<label>(plan.size());
    return n;
}

label Decomposition::minCellsPerDomain() const
{
    label m = labelMax;
    for (label d = 0; d < nDomains_; ++d) m = std::min(m, nCells(d));
    return m;
}

label Decomposition::maxCellsPerDomain() const
{
    label m = 0;
    for (label d = 0; d < nDomains_; ++d) m = std::max(m, nCells(d));
    return m;
}

// ------------------------------------------------------------------
// DomainMesh
// ------------------------------------------------------------------

DomainMesh::DomainMesh(const Decomposition& dec, label domain)
    : domain_(domain),
      cells_(dec.cells()[static_cast<std::size_t>(domain)]),
      haloCells_(dec.haloCells()[static_cast<std::size_t>(domain)]),
      faces_(dec.faces(domain))
{
    nSlots_ = dec.nLocal(domain);

    const polyMesh& mesh = dec.mesh();
    // Only owned cells are ever updated, so their geometry is what a local
    // stencil divides by; halo slots carry values, not volumes.
    cellVolumes_.assign(static_cast<std::size_t>(cells_.size()), 0.0);
    cellCentres_.assign(static_cast<std::size_t>(cells_.size()), Vector::zero());
    for (std::size_t i = 0; i < cells_.size(); ++i)
    {
        cellVolumes_[i] = mesh.cellVolumes()[static_cast<std::size_t>(cells_[i])];
        cellCentres_[i] = mesh.cellCentres()[static_cast<std::size_t>(cells_[i])];
    }

    const label nF = static_cast<label>(faces_.size());
    sf_.assign(static_cast<std::size_t>(nF), Vector::zero());
    faceCentres_.assign(static_cast<std::size_t>(nF), Vector::zero());
    delta_.assign(static_cast<std::size_t>(nF), Vector::zero());
    ownerWeight_.assign(static_cast<std::size_t>(nF), 0.0);

    const Field<scalar> ow = makeInternalOwnerWeights(mesh);
    const label nInt = mesh.nInternalFaces();
    for (std::size_t k = 0; k < faces_.size(); ++k)
    {
        const label f = faces_[k].face;
        const std::size_t u = static_cast<std::size_t>(f);
        sf_[k] = mesh.faceSf()[u];
        faceCentres_[k] = mesh.faceCentres()[u];
        // The non-orthogonal lever is measured from the *global owner* centre,
        // so both sides of a face see the same one and compute the same flux;
        // a rank could only do that if the lever came with the face.
        delta_[k] = faceCentres_[k] - mesh.cellCentres()[static_cast<std::size_t>(mesh.owner()[u])];
        ownerWeight_[k] = (f < nInt) ? ow[u] : 1.0;
    }
}

// ------------------------------------------------------------------
// DomainSet
// ------------------------------------------------------------------

DomainSet::DomainSet(const Decomposition& dec) : dec_(dec)
{
    domains_.reserve(static_cast<std::size_t>(dec.nDomains()));
    for (label d = 0; d < dec.nDomains(); ++d) domains_.emplace_back(dec, d);
}

void DomainSet::scatter(const Field<scalar>& global, std::vector<Field<scalar>>& local) const
{
    for (label d = 0; d < dec_.nDomains(); ++d)
    {
        const std::vector<label>& cells = dec_.cells()[static_cast<std::size_t>(d)];
        Field<scalar>& dst = local[static_cast<std::size_t>(d)];
        for (std::size_t i = 0; i < cells.size(); ++i) dst[i] = global[static_cast<std::size_t>(cells[i])];
    }
}

Field<scalar> DomainSet::gather(const std::vector<Field<scalar>>& local) const
{
    Field<scalar> out(static_cast<std::size_t>(dec_.nCells()), 0.0);
    for (label d = 0; d < dec_.nDomains(); ++d)
    {
        const std::vector<label>& cells = dec_.cells()[static_cast<std::size_t>(d)];
        const Field<scalar>& src = local[static_cast<std::size_t>(d)];
        for (std::size_t i = 0; i < cells.size(); ++i) out[static_cast<std::size_t>(cells[i])] = src[i];
    }
    return out;
}

Field<Vector> DomainSet::gather(const std::vector<Field<Vector>>& local) const
{
    Field<Vector> out(static_cast<std::size_t>(dec_.nCells()), Vector::zero());
    for (label d = 0; d < dec_.nDomains(); ++d)
    {
        const std::vector<label>& cells = dec_.cells()[static_cast<std::size_t>(d)];
        const Field<Vector>& src = local[static_cast<std::size_t>(d)];
        for (std::size_t i = 0; i < cells.size(); ++i) out[static_cast<std::size_t>(cells[i])] = src[i];
    }
    return out;
}

// ------------------------------------------------------------------
// reductions
// ------------------------------------------------------------------

scalar reduceSum(const std::vector<scalar>& partials)
{
    if (partials.empty()) cfFatal("reduce: no domain contribution to sum");
    return std::accumulate(partials.begin(), partials.end(), 0.0);
}

scalar reduceMin(const std::vector<scalar>& partials)
{
    if (partials.empty()) cfFatal("reduce: no domain contribution to minimise");
    return *std::min_element(partials.begin(), partials.end());
}

scalar reduceMax(const std::vector<scalar>& partials)
{
    if (partials.empty()) cfFatal("reduce: no domain contribution to maximise");
    return *std::max_element(partials.begin(), partials.end());
}

} // namespace parallel
} // namespace cf
