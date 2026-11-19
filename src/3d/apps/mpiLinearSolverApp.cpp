// chflow-3d-mpi-linear-solver: distributed sparse linear-system solve over a
// REAL MPI transport (A11 / §5.122). One domain per rank.
//
// A global graph-Laplacian CSR is assembled from the polyMesh cell adjacency
// (diagonal = face count + shift, off-diagonal -1 per internal face; SPD for
// shift > 0). Its rows are distributed by the coordinate-median Decomposition:
// rank r owns the global rows dec.cells(r). Each rank assembles only its owned
// rows into a local CSR whose column indices are domain-local slots (owned
// cells first, then haloed cells).
//
// The iterative shell is a diagonal-preconditioned CG with the same algorithm
// as cf::solveCG, but every matvec is a real MPI_Sendrecv halo exchange of the
// iterate followed by the rank's local CSR product, and every dot product and
// residual norm is an MPI_Allreduce. The gathered owned solution is compared
// with the serial solveCG on the SAME global CSR, and with a known analytic
// solution whose exact right-hand side b = A*xExact was used to drive solve.
// cf3d / cf_common stay MPI-free; MPI lives only here and in mpiTransport.h.
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include <mpi.h>

#include "cf/common/Dictionary.h"
#include "cf/common/Error.h"
#include "cf/common/Field.h"
#include "cf/common/Log.h"
#include "cf/common/Matrix.h"
#include "cf/fvm/MeshLoader.h"
#include "cf/fvm/PolyMesh.h"
#include "cf/parallel/DomainDecomposition.h"

#include "mpiTransport.h"

using namespace cf;

namespace
{

std::string num(double v)
{
    std::ostringstream os;
    os << std::setprecision(17) << v;
    return os.str();
}

/// c0 + c1 x + c2 y + c3 z + c4 x^2 + c5 y^2 + c6 z^2 + c7 xy + c8 yz + c9 zx
double polyValue(const std::vector<double>& c, const Vector& p)
{
    return c[0] + c[1] * p.x + c[2] * p.y + c[3] * p.z
         + c[4] * p.x * p.x + c[5] * p.y * p.y + c[6] * p.z * p.z
         + c[7] * p.x * p.y + c[8] * p.y * p.z + c[9] * p.z * p.x;
}

double maxDev(const std::vector<double>& a, const std::vector<double>& b)
{
    double m = 0;
    for (std::size_t i = 0; i < a.size(); ++i) m = std::max(m, std::abs(a[i] - b[i]));
    return m;
}

using Vec = Field<double>;

/// One rank's owned rows in CSR form. Columns are domain-local slots so a
/// product reads only owned + halo slots of a domain-local iterate vector.
struct LocalCSR
{
    label nRows = 0;
    std::vector<label> rowStart, cols;
    std::vector<double> vals, invD;
};

LocalCSR extractLocal(const SparseMatrix& A, const parallel::Decomposition& dec, label rank)
{
    const label nOwn = dec.nCells(rank);
    LocalCSR L;
    L.nRows = nOwn;
    L.rowStart.assign(static_cast<std::size_t>(nOwn) + 1, 0);
    for (label i = 0; i < nOwn; ++i)
    {
        const label g = dec.cells()[static_cast<std::size_t>(rank)][static_cast<std::size_t>(i)];
        L.rowStart[static_cast<std::size_t>(i)] = static_cast<label>(L.cols.size());
        for (label k = A.rowStart()[static_cast<std::size_t>(g)];
             k < A.rowStart()[static_cast<std::size_t>(g) + 1]; ++k)
        {
            const label gc = A.cols()[static_cast<std::size_t>(k)];
            const label s = dec.localSlot(rank, gc);
            if (s < 0)
                cfFatal("mpiLinearSolver: row " + std::to_string(g) + " needs column "
                        + std::to_string(gc) + " which the decomposition neither owns nor halos");
            L.cols.push_back(s);
            L.vals.push_back(A.vals()[static_cast<std::size_t>(k)]);
        }
    }
    L.rowStart[static_cast<std::size_t>(nOwn)] = static_cast<label>(L.cols.size());

    // Jacobi preconditioner from the owned rows' diagonal (no halo needed).
    L.invD.assign(static_cast<std::size_t>(nOwn), 1.0);
    for (label i = 0; i < nOwn; ++i)
    {
        for (label k = L.rowStart[static_cast<std::size_t>(i)];
             k < L.rowStart[static_cast<std::size_t>(i) + 1]; ++k)
        {
            if (L.cols[static_cast<std::size_t>(k)] == i)
            {
                const double d = L.vals[static_cast<std::size_t>(k)];
                L.invD[static_cast<std::size_t>(i)] = (std::abs(d) > 1e-300) ? 1.0 / d : 1.0;
                break;
            }
        }
    }
    return L;
}

/// y_owned = A_local * x. The halo slots of x are filled from peers first by
/// MPI_Sendrecv; a copy is exchanged so the caller's iterate is never mutated.
Vec localMatvec(const parallel::Decomposition& dec, label rank, const LocalCSR& L,
                const Vec& x, MPI_Comm comm)
{
    Vec xl(x);
    parallel::mpi::exchange(dec, rank, xl, comm);
    Vec y(static_cast<std::size_t>(L.nRows), 0.0);
    for (label i = 0; i < L.nRows; ++i)
    {
        double s = 0;
        for (label k = L.rowStart[static_cast<std::size_t>(i)];
             k < L.rowStart[static_cast<std::size_t>(i) + 1]; ++k)
        {
            s += L.vals[static_cast<std::size_t>(k)]
               * xl[static_cast<std::size_t>(L.cols[static_cast<std::size_t>(k)])];
        }
        y[static_cast<std::size_t>(i)] = s;
    }
    return y;
}

/// Owned partial dot, summed over all ranks. Owned slots are 0..nOwn-1 by
/// construction of Decomposition, so the global dot is exactly the serial one.
double globalDot(const parallel::Decomposition& dec, label rank,
                 const Vec& a, const Vec& b, MPI_Comm comm)
{
    double s = 0;
    const label nOwn = dec.nCells(rank);
    for (label i = 0; i < nOwn; ++i)
        s += a[static_cast<std::size_t>(i)] * b[static_cast<std::size_t>(i)];
    return parallel::mpi::reduceSum(s, comm);
}

/// Owned partial inf-norm. MPI_MAX is order independent, hence exact.
double globalInfNorm(const parallel::Decomposition& dec, label rank,
                     const Vec& a, MPI_Comm comm)
{
    double m = 0;
    const label nOwn = dec.nCells(rank);
    for (label i = 0; i < nOwn; ++i)
        m = std::max(m, std::abs(a[static_cast<std::size_t>(i)]));
    return parallel::mpi::reduceMax(m, comm);
}

struct DistReport
{
    int iterations = 0;
    double finalResidual = 0;
    double bnorm = 0;
    bool converged = false;
    Vec x; // domain-local; owned slots hold this rank's solution
};

/// Diagonal-preconditioned CG, distributed. This mirrors cf::solveCG step for
/// step: the only changes are that A.product becomes localMatvec (halo + local
/// CSR) and dot/infNorm become MPI_Allreduce collectives.
DistReport solveDistCG(const parallel::Decomposition& dec, label rank, const LocalCSR& L,
                       const std::vector<double>& bOwned, double tol, int maxIter,
                       MPI_Comm comm)
{
    const label nOwn = dec.nCells(rank);
    const label nLocal = dec.nLocal(rank);
    const std::size_t no = static_cast<std::size_t>(nOwn);

    DistReport rep;
    Vec x(static_cast<std::size_t>(nLocal), 0.0);
    Vec r(no, 0.0), z(no, 0.0), p(no, 0.0), Ap(no, 0.0);

    const Vec Ax = localMatvec(dec, rank, L, x, comm);
    for (label i = 0; i < nOwn; ++i)
        r[static_cast<std::size_t>(i)] =
            bOwned[static_cast<std::size_t>(i)] - Ax[static_cast<std::size_t>(i)];

    double bInf = 0;
    for (label i = 0; i < nOwn; ++i)
        bInf = std::max(bInf, std::abs(bOwned[static_cast<std::size_t>(i)]));
    rep.bnorm = std::max(parallel::mpi::reduceMax(bInf, comm), 1e-300);

    rep.finalResidual = globalInfNorm(dec, rank, r, comm);
    if (rep.finalResidual / rep.bnorm < tol)
    {
        rep.converged = true;
        rep.x = x;
        return rep;
    }

    for (label i = 0; i < nOwn; ++i)
        z[static_cast<std::size_t>(i)] =
            L.invD[static_cast<std::size_t>(i)] * r[static_cast<std::size_t>(i)];

    Vec pLocal(static_cast<std::size_t>(nLocal), 0.0);
    for (label i = 0; i < nOwn; ++i)
        pLocal[static_cast<std::size_t>(i)] = z[static_cast<std::size_t>(i)];
    double rz = globalDot(dec, rank, r, z, comm);

    for (int it = 1; it <= maxIter; ++it)
    {
        Ap = localMatvec(dec, rank, L, pLocal, comm);
        const double pAp = globalDot(dec, rank, pLocal, Ap, comm);
        if (pAp <= 0 || std::abs(rz) < 1e-300) break;

        const double alpha = rz / pAp;
        for (label i = 0; i < nOwn; ++i)
        {
            const std::size_t u = static_cast<std::size_t>(i);
            x[u] += alpha * pLocal[u];
            r[u] -= alpha * Ap[u];
        }
        rep.iterations = it;
        rep.finalResidual = globalInfNorm(dec, rank, r, comm);
        if (rep.finalResidual / rep.bnorm < tol)
        {
            rep.converged = true;
            break;
        }

        for (label i = 0; i < nOwn; ++i)
            z[static_cast<std::size_t>(i)] =
                L.invD[static_cast<std::size_t>(i)] * r[static_cast<std::size_t>(i)];
        const double rzNew = globalDot(dec, rank, r, z, comm);
        const double beta = rzNew / rz;
        rz = rzNew;
        for (label i = 0; i < nOwn; ++i)
        {
            const std::size_t u = static_cast<std::size_t>(i);
            pLocal[u] = z[u] + beta * pLocal[u];
        }
    }

    rep.x = x;
    return rep;
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank = 0, size = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    int exitCode = 0;
    try
    {
        if (argc < 2)
            cfFatal("usage: chflow-3d-mpi-linear-solver <case.dict> (run under mpiexec -np N)");
        const std::string casePath = argv[1];
        const Dictionary dict = Dictionary::parseFile(casePath);
        const polyMesh mesh = MeshLoader::fromCase(dict, casePath);
        const label nC = mesh.nCells();
        const label nDomains = static_cast<label>(size);

        const Dictionary& ls = dict.getDict("linearSolver");
        const double tol = ls.getScalar("tol");
        const int maxIter = ls.getInt("maxIter");
        const double shift = ls.getScalar("shift");
        const std::vector<double> coeffs = ls.getScalarList("exactCoeffs");
        if (coeffs.size() != 10)
            cfFatal("mpiLinearSolver: exactCoeffs needs ten coefficients c0..c9");

        std::string dumpName = "solution.dat";
        if (dict.has("output"))
            dumpName = dict.getDict("output").getWordOr("dump", dumpName);

        const parallel::Decomposition dec(mesh, nDomains);

        // Global graph Laplacian: off-diagonal -1 for each internal face,
        // diagonal = number of incident faces + shift. Symmetric positive
        // definite for shift > 0, so CG (and the checker) is well posed.
        SparseMatrix::Builder gb(nC);
        const Field<label>& owner = mesh.owner();
        const Field<label>& neigh = mesh.neighbour();
        for (label f = 0; f < mesh.nFaces(); ++f)
        {
            const label o = owner[static_cast<std::size_t>(f)];
            gb.add(o, o, 1.0);
            if (f < mesh.nInternalFaces())
            {
                const label nb = neigh[static_cast<std::size_t>(f)];
                gb.add(nb, nb, 1.0);
                gb.add(o, nb, -1.0);
                gb.add(nb, o, -1.0);
            }
        }
        for (label c = 0; c < nC; ++c) gb.add(c, c, shift);
        const SparseMatrix A = gb.finish();
        const label nnz = static_cast<label>(A.cols().size());

        // Known analytic solution and the exact right-hand side b = A*xExact.
        Field<double> xExact(static_cast<std::size_t>(nC), 0.0);
        for (label c = 0; c < nC; ++c)
            xExact[static_cast<std::size_t>(c)] =
                polyValue(coeffs, mesh.cellCentres()[static_cast<std::size_t>(c)]);
        const std::vector<double> xev(xExact.begin(), xExact.end());
        std::vector<double> bvec;
        A.product(xev, bvec);
        const Field<double> bGlobal(bvec.begin(), bvec.end());

        if (rank == 0)
            log::info("mpiLinearSolver start ranks ", nDomains, " cells ", nC,
                      " nnz ", nnz, " shift ", num(shift), " tol ", num(tol),
                      " maxIter ", maxIter);

        const LocalCSR L = extractLocal(A, dec, static_cast<label>(rank));
        const label nOwn = dec.nCells(static_cast<label>(rank));
        std::vector<double> bOwned(static_cast<std::size_t>(nOwn));
        for (label i = 0; i < nOwn; ++i)
            bOwned[static_cast<std::size_t>(i)] = bGlobal[static_cast<std::size_t>(
                dec.cells()[static_cast<std::size_t>(rank)][static_cast<std::size_t>(i)])];

        const DistReport rep =
            solveDistCG(dec, static_cast<label>(rank), L, bOwned, tol, maxIter, comm);
        int myConv = rep.converged ? 1 : 0, allConv = 0;
        MPI_Allreduce(&myConv, &allConv, 1, MPI_INT, MPI_SUM, comm);

        // Gather this rank's owned solution into global cell order on rank 0.
        const Field<double> xDist =
            parallel::mpi::gatherScalar(dec, static_cast<label>(rank), rep.x, comm);

        // Verify the distributed matvec itself: the gathered y = A_local*x
        // must equal the serial CSR product A*x on the very same vector.
        const Vec yLocal = localMatvec(dec, static_cast<label>(rank), L, rep.x, comm);
        const Field<double> yGather =
            parallel::mpi::gatherScalar(dec, static_cast<label>(rank), yLocal, comm);

        if (rank == 0)
        {
            Field<double> xSerial(static_cast<std::size_t>(nC), 0.0);
            const SolverReport sr = solveCG(A, bGlobal, xSerial, tol, maxIter);
            const std::vector<double> xdv(xDist.begin(), xDist.end());
            std::vector<double> ySerial;
            A.product(xdv, ySerial);

            const double matvecDev =
                maxDev(std::vector<double>(yGather.begin(), yGather.end()), ySerial);
            const double maxDiff =
                maxDev(xdv, std::vector<double>(xSerial.begin(), xSerial.end()));
            const double analyticDiff =
                maxDev(xdv, std::vector<double>(xExact.begin(), xExact.end()));
            const double serAnalyticDiff =
                maxDev(std::vector<double>(xSerial.begin(), xSerial.end()), xev);

            log::info("mpiLinearSolver solve ranks ", nDomains, " iters ", rep.iterations,
                      " serIters ", sr.iterations, " converged ", (rep.converged ? 1 : 0),
                      " convRanks ", allConv, "/", size,
                      " resRel ", num(rep.finalResidual / std::max(rep.bnorm, 1e-300)),
                      " matvecDev ", num(matvecDev));
            log::info("mpiLinearSolver verify ranks ", nDomains, " maxDiff ", num(maxDiff),
                      " analyticDiff ", num(analyticDiff),
                      " serAnalyticDiff ", num(serAnalyticDiff),
                      " bnorm ", num(rep.bnorm));

            const std::string::size_type dot = dumpName.find_last_of('.');
            const std::string base =
                (dot == std::string::npos) ? dumpName : dumpName.substr(0, dot);
            const std::string path = base + "_" + std::to_string(nDomains) + ".dat";
            std::ofstream out(path);
            out << std::setprecision(17);
            for (label c = 0; c < nC; ++c)
            {
                const Vector& p = mesh.cellCentres()[static_cast<std::size_t>(c)];
                out << c << ' ' << xDist[static_cast<std::size_t>(c)] << ' '
                    << xSerial[static_cast<std::size_t>(c)] << ' '
                    << xExact[static_cast<std::size_t>(c)] << ' '
                    << p.x << ' ' << p.y << ' ' << p.z << '\n';
            }
        }
    }
    catch (const std::exception& e)
    {
        log::error("mpiLinearSolver (rank ", rank, "): ", e.what());
        exitCode = 1;
    }

    // A failure on any rank must fail the whole run.
    int anyBad = 0;
    MPI_Allreduce(&exitCode, &anyBad, 1, MPI_INT, MPI_MAX, comm);
    MPI_Finalize();
    return anyBad ? 1 : 0;
}
