#pragma once

// Real MPI transport for the decomposed stencils (A11). This header is used only
// by the chflow-3d-mpi verification app; the cf3d library itself stays MPI-free.
// The buffers it posts are exactly the ones the serial copy-loop exchange uses:
// one message per peer, the sender's owned slots and the receiver's halo slots
// index-matched by Decomposition, so a distributed halo holds the same doubles
// the serial one does and the per-domain kernels reproduce the serial sweep bit
// for bit. The one place a real collective differs from the serial reduction is
// the SUM: MPI_Allreduce adds the partials in a tree order, so it matches the
// serial left-to-right std::accumulate to round-off, not to the last bit. MIN and
// MAX are order-independent and stay exact.

#ifdef CF_HAVE_MPI

#include <mpi.h>

#include <type_traits>
#include <vector>

#include "cf/common/Error.h"
#include "cf/common/Field.h"
#include "cf/common/Vector.h"
#include "cf/parallel/DomainDecomposition.h"

namespace cf
{
namespace parallel
{
namespace mpi
{

template<typename T>
inline constexpr int nComp = std::is_same_v<T, scalar> ? 1 : 3;

template<typename T>
inline void packSlot(const Field<T>& local, label slot, scalar* out)
{
    if constexpr (std::is_same_v<T, scalar>)
    {
        out[0] = local[static_cast<std::size_t>(slot)];
    }
    else
    {
        const Vector& v = local[static_cast<std::size_t>(slot)];
        out[0] = v.x; out[1] = v.y; out[2] = v.z;
    }
}

template<typename T>
inline void unpackSlot(Field<T>& local, label slot, const scalar* in)
{
    if constexpr (std::is_same_v<T, scalar>)
    {
        local[static_cast<std::size_t>(slot)] = in[0];
    }
    else
    {
        local[static_cast<std::size_t>(slot)] = Vector{in[0], in[1], in[2]};
    }
}

/// Fill this rank's halo slots from its peers' owned slots over MPI_Sendrecv.
/// `rank` is this process's domain. Returns the number of values this rank sent
/// (one per coupled-face side it owns), the same count the serial exchange
/// moves for that domain. Owned slots are only read, halo slots only written.
template<typename T>
label exchange(const Decomposition& dec, label rank, Field<T>& local, MPI_Comm comm)
{
    constexpr int nc = nComp<T>;
    label moved = 0;
    for (const Decomposition::Message& m : dec.sendPlan(rank))
    {
        const Decomposition::Message* r = dec.recvFrom(rank, m.peer);
        if (!r || r->slots.size() != m.slots.size())
        {
            cfFatal("mpi exchange: rank " + std::to_string(rank) + " sends "
                    + std::to_string(m.slots.size()) + " values to "
                    + std::to_string(m.peer) + ", which did not ask for that many");
        }
        const int ns = static_cast<int>(m.slots.size());
        std::vector<scalar> sbuf(static_cast<std::size_t>(ns * nc));
        std::vector<scalar> rbuf(static_cast<std::size_t>(ns * nc));
        for (int k = 0; k < ns; ++k) packSlot<T>(local, m.slots[static_cast<std::size_t>(k)], &sbuf[static_cast<std::size_t>(k * nc)]);
        MPI_Sendrecv(sbuf.data(), ns * nc, MPI_DOUBLE, static_cast<int>(m.peer), 0,
                     rbuf.data(), ns * nc, MPI_DOUBLE, static_cast<int>(m.peer), 0,
                     comm, MPI_STATUS_IGNORE);
        for (int k = 0; k < ns; ++k) unpackSlot<T>(local, r->slots[static_cast<std::size_t>(k)], &rbuf[static_cast<std::size_t>(k * nc)]);
        moved += ns;
    }
    return moved;
}

inline scalar reduceSum(scalar partial, MPI_Comm comm)
{
    scalar g = 0;
    MPI_Allreduce(&partial, &g, 1, MPI_DOUBLE, MPI_SUM, comm);
    return g;
}

inline scalar reduceMin(scalar partial, MPI_Comm comm)
{
    scalar g = 0;
    MPI_Allreduce(&partial, &g, 1, MPI_DOUBLE, MPI_MIN, comm);
    return g;
}

inline scalar reduceMax(scalar partial, MPI_Comm comm)
{
    scalar g = 0;
    MPI_Allreduce(&partial, &g, 1, MPI_DOUBLE, MPI_MAX, comm);
    return g;
}

inline label reduceSumInt(label partial, MPI_Comm comm)
{
    int p = static_cast<int>(partial), g = 0;
    MPI_Allreduce(&p, &g, 1, MPI_INT, MPI_SUM, comm);
    return static_cast<label>(g);
}

/// Gather the owned slots of every rank into global cell order on rank 0.
/// Other ranks get an empty field. The owned slots are the first nCells(rank)
/// entries of `local`; rank 0 places them at the global ids the decomposition
/// assigned, so the result is comparable cell by cell with a serial sweep.
inline Field<scalar> gatherScalar(const Decomposition& dec, label rank,
                                  const Field<scalar>& local, MPI_Comm comm)
{
    const int n = static_cast<int>(dec.nDomains());
    std::vector<int> counts(static_cast<std::size_t>(n)), displs(static_cast<std::size_t>(n));
    int off = 0;
    for (int d = 0; d < n; ++d)
    {
        counts[static_cast<std::size_t>(d)] = static_cast<int>(dec.nCells(d));
        displs[static_cast<std::size_t>(d)] = off;
        off += counts[static_cast<std::size_t>(d)];
    }
    const int nOwn = static_cast<int>(dec.nCells(rank));
    std::vector<scalar> sendbuf(local.begin(), local.begin() + nOwn);
    std::vector<scalar> recvbuf(static_cast<std::size_t>(off));
    MPI_Gatherv(sendbuf.data(), nOwn, MPI_DOUBLE,
                recvbuf.data(), counts.data(), displs.data(), MPI_DOUBLE, 0, comm);

    Field<scalar> out(static_cast<std::size_t>(dec.nCells()), 0.0);
    if (rank == 0)
    {
        for (int d = 0; d < n; ++d)
        {
            for (int i = 0; i < counts[static_cast<std::size_t>(d)]; ++i)
            {
                const label c = dec.cells()[static_cast<std::size_t>(d)][static_cast<std::size_t>(i)];
                out[static_cast<std::size_t>(c)] = recvbuf[static_cast<std::size_t>(displs[static_cast<std::size_t>(d)] + i)];
            }
        }
    }
    return out;
}

inline Field<Vector> gatherVector(const Decomposition& dec, label rank,
                                  const Field<Vector>& local, MPI_Comm comm)
{
    const int n = static_cast<int>(dec.nDomains());
    std::vector<int> counts(static_cast<std::size_t>(n)), displs(static_cast<std::size_t>(n));
    int off = 0;
    for (int d = 0; d < n; ++d)
    {
        counts[static_cast<std::size_t>(d)] = 3 * static_cast<int>(dec.nCells(d));
        displs[static_cast<std::size_t>(d)] = off;
        off += counts[static_cast<std::size_t>(d)];
    }
    const int nOwn = static_cast<int>(dec.nCells(rank));
    std::vector<scalar> sendbuf(static_cast<std::size_t>(3 * nOwn));
    for (int i = 0; i < nOwn; ++i) packSlot<Vector>(local, static_cast<label>(i), &sendbuf[static_cast<std::size_t>(3 * i)]);
    std::vector<scalar> recvbuf(static_cast<std::size_t>(off));
    MPI_Gatherv(sendbuf.data(), 3 * nOwn, MPI_DOUBLE,
                recvbuf.data(), counts.data(), displs.data(), MPI_DOUBLE, 0, comm);

    Field<Vector> out(static_cast<std::size_t>(dec.nCells()), Vector::zero());
    if (rank == 0)
    {
        for (int d = 0; d < n; ++d)
        {
            const int nd = counts[static_cast<std::size_t>(d)] / 3;
            for (int i = 0; i < nd; ++i)
            {
                const label c = dec.cells()[static_cast<std::size_t>(d)][static_cast<std::size_t>(i)];
                const scalar* p = &recvbuf[static_cast<std::size_t>(displs[static_cast<std::size_t>(d)] + 3 * i)];
                out[static_cast<std::size_t>(c)] = Vector{p[0], p[1], p[2]};
            }
        }
    }
    return out;
}

} // namespace mpi
} // namespace parallel
} // namespace cf

#endif // CF_HAVE_MPI