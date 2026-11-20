// chflow-3d-petsc-distributed: solve the SAME distributed system as
// chflow-mpi-linear-solver (§5.122) with PETSc's own MPI machinery, so the
// two independent distributed solvers cross-check each other (A11 / the
// "PETSc 大规模分布式对照" boundary item).
//
// The system is the same graph-Laplacian CSR from the polyMesh cell adjacency
// (diagonal = face count + shift, off-diagonal -1 per internal face; SPD for
// shift > 0) with the exact right-hand side b = A*xExact for a known
// polynomial. Every rank builds the same global mesh and matrix (the case is
// small enough for that), then PETSc takes over: rows are owned by PETSc's
// own contiguous block partition -- deliberately NOT the coordinate-median
// Decomposition the homegrown solver uses, so the cross-check compares two
// genuinely different domain decompositions. Assembly goes through
// MatCreateAIJ with exact per-row diagonal/off-diagonal preallocation counts
// taken from the CSR (the sequential branch's zero-preallocation malloc bug
// does not get a chance to repeat), the solve is KSPCG with Jacobi
// preconditioning, and rank 0 dumps cell-by-cell (cell, xPetsc, xSerial,
// xExact, centre) in the same six-column format the mpiLinearSolver dump
// uses, so one checker can compare both dumps directly.
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include <mpi.h>
#include <petsc.h>

#include "cf/common/Dictionary.h"
#include "cf/common/Error.h"
#include "cf/common/Field.h"
#include "cf/common/Log.h"
#include "cf/common/Matrix.h"
#include "cf/fvm/MeshLoader.h"
#include "cf/fvm/PolyMesh.h"

using namespace cf;

// CHKERRABORT wrapper: a PETSc failure anywhere aborts the whole MPI run
// loudly instead of silently returning an error code nobody reads.
#define PCHK(call)                                                              \
    do                                                                          \
    {                                                                           \
        PetscErrorCode ierr_ = (call);                                          \
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);                                   \
    } while (0)


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

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    MPI_Comm comm = PETSC_COMM_WORLD;
    int rank = 0, size = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    int exitCode = 0;
    try
    {
        if (argc < 2)
            cfFatal("usage: chflow-3d-petsc-distributed <case.dict> (run under mpiexec -np N)");
        const std::string casePath = argv[1];
        const Dictionary dict = Dictionary::parseFile(casePath);
        const polyMesh mesh = MeshLoader::fromCase(dict, casePath);
        const label nC = mesh.nCells();

        const Dictionary& ls = dict.getDict("linearSolver");
        const double tol = ls.getScalar("tol");
        const int maxIter = ls.getInt("maxIter");
        const double shift = ls.getScalar("shift");
        const std::vector<double> coeffs = ls.getScalarList("exactCoeffs");
        if (coeffs.size() != 10)
            cfFatal("petscDist: exactCoeffs needs ten coefficients c0..c9");

        std::string dumpName = "petscSolution.dat";
        if (dict.has("output"))
            dumpName = dict.getDict("output").getWordOr("dump", dumpName);

        // Same global graph Laplacian as the homegrown distributed solver.
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

        // Known analytic solution and the exact right-hand side b = A*xExact.
        Field<double> xExact(static_cast<std::size_t>(nC), 0.0);
        for (label c = 0; c < nC; ++c)
            xExact[static_cast<std::size_t>(c)] =
                polyValue(coeffs, mesh.cellCentres()[static_cast<std::size_t>(c)]);
        const std::vector<double> xev(xExact.begin(), xExact.end());
        std::vector<double> bvec;
        A.product(xev, bvec);

        if (rank == 0)
            log::info("petscDist start ranks ", size, " cells ", nC,
                      " nnz ", A.cols().size(), " shift ", num(shift),
                      " tol ", num(tol), " maxIter ", maxIter);

        // PETSc's own contiguous block ownership. This is deliberately a
        // different decomposition from the coordinate-median one the
        // homegrown CG uses -- the cross-check must compare two genuinely
        // independent distributed paths, not one partition twice.
        Mat M;
        PetscInt n, N;
        MatCreateAIJ(comm, PETSC_DECIDE, PETSC_DECIDE,
                     static_cast<PetscInt>(nC), static_cast<PetscInt>(nC),
                     0, nullptr, 0, nullptr, &M);
        MatGetLocalSize(M, &n, nullptr);
        MatGetSize(M, &N, nullptr);
        PetscInt lo = 0, hi = 0;
        MatGetOwnershipRange(M, &lo, &hi);

        // Exact per-row preallocation from the CSR row census (dnnz = columns
        // inside this rank's row block, onnz = the rest).
        std::vector<PetscInt> dnnz, onnz;
        for (PetscInt r = lo; r < hi; ++r)
        {
            PetscInt d = 0, o = 0;
            for (label k = A.rowStart()[static_cast<std::size_t>(r)];
                 k < A.rowStart()[static_cast<std::size_t>(r) + 1]; ++k)
            {
                const PetscInt col = static_cast<PetscInt>(A.cols()[static_cast<std::size_t>(k)]);
                (col >= lo && col < hi) ? ++d : ++o;
            }
            dnnz.push_back(d);
            onnz.push_back(o);
        }
        // MatXAIJSetPreallocation is the type-agnostic form: with np == 1
        // PETSc may materialise the matrix as SEQAIJ (diag block only), and
        // MatMPIAIJSetPreallocation does not transfer the counts on that
        // path -- MatXAIJSetPreallocation routes to the right variant.
        PCHK(MatXAIJSetPreallocation(M, 1, dnnz.data(), onnz.data(),
                                     nullptr, nullptr));

        for (PetscInt r = lo; r < hi; ++r)
        {
            for (label k = A.rowStart()[static_cast<std::size_t>(r)];
                 k < A.rowStart()[static_cast<std::size_t>(r) + 1]; ++k)
            {
                const PetscInt col = static_cast<PetscInt>(A.cols()[static_cast<std::size_t>(k)]);
                const PetscScalar val = A.vals()[static_cast<std::size_t>(k)];
                PCHK(MatSetValues(M, 1, &r, 1, &col, &val, ADD_VALUES));
            }
        }
        PCHK(MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY));
        PCHK(MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY));

        Vec b, x, xe;
        PCHK(VecCreateMPI(comm, n, N, &b));
        PCHK(VecDuplicate(b, &x));
        PCHK(VecDuplicate(b, &xe));
        std::vector<PetscInt> rows;
        std::vector<PetscScalar> bLocal, xeLocal;
        for (PetscInt r = lo; r < hi; ++r)
        {
            rows.push_back(r);
            bLocal.push_back(bvec[static_cast<std::size_t>(r)]);
            xeLocal.push_back(xExact[static_cast<std::size_t>(r)]);
        }
        PCHK(VecSetValues(b, static_cast<PetscInt>(rows.size()), rows.data(),
                          bLocal.data(), INSERT_VALUES));
        PCHK(VecSetValues(xe, static_cast<PetscInt>(rows.size()), rows.data(),
                          xeLocal.data(), INSERT_VALUES));
        PCHK(VecAssemblyBegin(b));
        PCHK(VecAssemblyEnd(b));
        PCHK(VecAssemblyBegin(xe));
        PCHK(VecAssemblyEnd(xe));

        KSP ksp;
        PCHK(KSPCreate(comm, &ksp));
        PCHK(KSPSetOperators(ksp, M, M));
        PCHK(KSPSetType(ksp, KSPCG));
        PC pc;
        PCHK(KSPGetPC(ksp, &pc));
        PCHK(PCSetType(pc, PCJACOBI));
        PCHK(KSPSetTolerances(ksp, static_cast<PetscReal>(tol),
                                      PETSC_DEFAULT, PETSC_DEFAULT,
                                      static_cast<PetscInt>(maxIter)));
        PCHK(KSPSolve(ksp, b, x));

        PetscReal resNorm = 0;
        PetscInt its = 0;
        KSPConvergedReason reason = KSP_CONVERGED_ITERATING;
        PCHK(KSPGetResidualNorm(ksp, &resNorm));
        PCHK(KSPGetIterationNumber(ksp, &its));
        PCHK(KSPGetConvergedReason(ksp, &reason));
        const bool converged = reason > 0;

        // Local analytic deviation, then MPI_MAX across ranks.
        const PetscScalar* xp = nullptr;
        PCHK(VecGetArrayRead(x, &xp));
        double localDev = 0;
        for (PetscInt i = 0; i < n; ++i)
        {
            const double d = std::abs(static_cast<double>(xp[i])
                                      - static_cast<double>(xeLocal[static_cast<std::size_t>(i)]));
            localDev = std::max(localDev, d);
        }

        double bnormLocal = 0;
        for (const PetscScalar v : bLocal)
            bnormLocal = std::max(bnormLocal, std::abs(static_cast<double>(v)));

        double analyticDiff = 0, bnorm = 0;
        MPI_Reduce(&localDev, &analyticDiff, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
        MPI_Reduce(&bnormLocal, &bnorm, 1, MPI_DOUBLE, MPI_MAX, 0, comm);

        // Gather the PETSc solution in global cell order on rank 0. Ownership
        // is contiguous blocks, so a Gatherv over rank order is cell order.
        const int myN = static_cast<int>(n);
        std::vector<int> counts(size, 0), displs(size, 0);
        MPI_Gather(&myN, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, comm);
        if (rank == 0)
        {
            int acc = 0;
            for (int r = 0; r < size; ++r)
            {
                displs[static_cast<std::size_t>(r)] = acc;
                acc += counts[static_cast<std::size_t>(r)];
            }
        }
        std::vector<double> xPetsc(rank == 0 ? static_cast<std::size_t>(nC) : 0, 0.0);
        MPI_Gatherv(xp, myN, MPI_DOUBLE,
                    xPetsc.data(), counts.data(), displs.data(), MPI_DOUBLE, 0, comm);
        // The gather above is the last read through xp; only now release it
        // (restoring before the gather left xp dangling and corrupted the
        // gathered solution).
        PCHK(VecRestoreArrayRead(x, &xp));

        if (rank == 0)
        {
            Field<double> xSerial(static_cast<std::size_t>(nC), 0.0);
            const SolverReport sr = solveCG(A, Field<double>(bvec.begin(), bvec.end()),
                                            xSerial, tol, maxIter);
            const double serAnalyticDiff = maxDev(
                std::vector<double>(xSerial.begin(), xSerial.end()), xev);
            const double petscSerialDiff = maxDev(xPetsc,
                std::vector<double>(xSerial.begin(), xSerial.end()));

            log::info("petscDist solve ranks ", size, " iters ", static_cast<label>(its),
                      " serIters ", sr.iterations, " converged ", (converged ? 1 : 0),
                      " reason ", static_cast<label>(reason),
                      " resNorm ", num(static_cast<double>(resNorm)),
                      " bnorm ", num(bnorm));
            log::info("petscDist verify ranks ", size,
                      " analyticDiff ", num(analyticDiff),
                      " petscSerialDiff ", num(petscSerialDiff),
                      " serAnalyticDiff ", num(serAnalyticDiff));

            const std::string::size_type dot = dumpName.find_last_of('.');
            const std::string base =
                (dot == std::string::npos) ? dumpName : dumpName.substr(0, dot);
            const std::string path = base + "_" + std::to_string(size) + ".dat";
            std::ofstream out(path);
            out << std::setprecision(17);
            for (label c = 0; c < nC; ++c)
            {
                const Vector& p = mesh.cellCentres()[static_cast<std::size_t>(c)];
                out << c << ' ' << xPetsc[static_cast<std::size_t>(c)] << ' '
                    << xSerial[static_cast<std::size_t>(c)] << ' '
                    << xExact[static_cast<std::size_t>(c)] << ' '
                    << p.x << ' ' << p.y << ' ' << p.z << '\n';
            }
        }
        PCHK(MatDestroy(&M));
        PCHK(VecDestroy(&b));
        PCHK(VecDestroy(&x));
        PCHK(VecDestroy(&xe));
        PCHK(KSPDestroy(&ksp));
    }
    catch (const std::exception& e)
    {
        log::error("petscDist (rank ", rank, "): ", e.what());
        exitCode = 1;
    }

    int anyBad = 0;
    MPI_Allreduce(&exitCode, &anyBad, 1, MPI_INT, MPI_MAX, comm);
    PetscFinalize();
    MPI_Finalize();
    return anyBad ? 1 : 0;
}
