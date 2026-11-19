#include "cf/common/Matrix.h"
#include "cf/common/Error.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

#ifdef CHFLOW_WITH_EIGEN
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCore>
#include <Eigen/OrderingMethods>
#endif

#ifdef CHFLOW_WITH_PETSC
#include <petscmat.h>
#include <petscksp.h>
#include <petscvec.h>
#endif

namespace cf
{
namespace
{

double dot(const std::vector<double>& a, const std::vector<double>& b)
{
    double s = 0;
    for (std::size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
    return s;
}

double infNorm(const std::vector<double>& v)
{
    double m = 0;
    for (const double x : v)
    {
        if (std::isnan(x)) return std::numeric_limits<double>::quiet_NaN();
        m = std::max(m, std::abs(x));
    }
    return m;
}

std::vector<double> bVec(const Field<double>& b)
{
    return std::vector<double>(b.begin(), b.end());
}

} // namespace

SparseMatrix SparseMatrix::Builder::finish() const
{
    return finishWithDiagonal({});
}

SparseMatrix SparseMatrix::Builder::finishWithDiagonal(const std::vector<double>& diag) const
{
    std::vector<label> rows = rows_;
    std::vector<label> cols = cols_;
    std::vector<double> vals = vals_;
    const label n = static_cast<label>(diag.size());
    rows.reserve(rows.size() + static_cast<std::size_t>(n));
    cols.reserve(cols.size() + static_cast<std::size_t>(n));
    vals.reserve(vals.size() + static_cast<std::size_t>(n));
    for (label r = 0; r < n; ++r)
    {
        rows.push_back(r);
        cols.push_back(r);
        vals.push_back(diag[static_cast<std::size_t>(r)]);
    }

    const std::size_t nnz = rows.size();
    std::vector<std::size_t> idx(nnz);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) {
        if (rows[a] != rows[b]) return rows[a] < rows[b];
        return cols[a] < cols[b];
    });

    SparseMatrix m;
    m.rowStart_.assign(static_cast<std::size_t>(n_) + 1, 0);
    m.cols_.reserve(nnz);
    m.vals_.reserve(nnz);

    label prevRow = -1;
    for (std::size_t k = 0; k < nnz; ++k)
    {
        const label r = rows[idx[k]];
        const label c = cols[idx[k]];
        const double v = vals[idx[k]];

        if (r != prevRow)
        {
            for (label rr = prevRow + 1; rr <= r; ++rr)
            {
                m.rowStart_[static_cast<std::size_t>(rr)] = static_cast<label>(m.cols_.size());
            }
            m.cols_.push_back(c);
            m.vals_.push_back(v);
            prevRow = r;
        }
        else if (m.cols_.back() == c)
        {
            m.vals_.back() += v; // merge duplicate triplet
        }
        else
        {
            m.cols_.push_back(c);
            m.vals_.push_back(v);
        }
    }
    for (label rr = prevRow + 1; rr <= n_; ++rr)
    {
        m.rowStart_[static_cast<std::size_t>(rr)] = static_cast<label>(m.cols_.size());
    }
    return m;
}

void SparseMatrix::product(const std::vector<double>& x, std::vector<double>& y) const
{
    y.assign(static_cast<std::size_t>(size()), 0.0);
    for (label r = 0; r < size(); ++r)
    {
        double s = 0;
        for (label k = rowStart_[static_cast<std::size_t>(r)];
             k < rowStart_[static_cast<std::size_t>(r) + 1]; ++k)
        {
            s += vals_[static_cast<std::size_t>(k)]
               * x[static_cast<std::size_t>(cols_[static_cast<std::size_t>(k)])];
        }
        y[static_cast<std::size_t>(r)] = s;
    }
}

std::vector<double> SparseMatrix::diagonal() const
{
    std::vector<double> d(static_cast<std::size_t>(size()), 0.0);
    for (label r = 0; r < size(); ++r)
    {
        for (label k = rowStart_[static_cast<std::size_t>(r)];
             k < rowStart_[static_cast<std::size_t>(r) + 1]; ++k)
        {
            if (cols_[static_cast<std::size_t>(k)] == r)
            {
                d[static_cast<std::size_t>(r)] = vals_[static_cast<std::size_t>(k)];
                break;
            }
        }
    }
    return d;
}

SolverReport solveCG(const SparseMatrix& A, const Field<double>& b,
                     Field<double>& x, double tol, int maxIter)
{
    const label n = A.size();
    std::vector<double> xx(x.begin(), x.end());
    const auto diag = A.diagonal();

    std::vector<double> invD(static_cast<std::size_t>(n));
    for (label i = 0; i < n; ++i)
    {
        const double d = diag[static_cast<std::size_t>(i)];
        invD[static_cast<std::size_t>(i)] = (std::abs(d) > 1e-300) ? 1.0 / d : 1.0;
    }

    std::vector<double> Ax(static_cast<std::size_t>(n), 0.0), r(static_cast<std::size_t>(n));
    A.product(xx, Ax);
    const std::vector<double> bb = bVec(b);
    for (label i = 0; i < n; ++i) r[static_cast<std::size_t>(i)] = bb[static_cast<std::size_t>(i)] - Ax[static_cast<std::size_t>(i)];

    SolverReport rep;
    rep.initialResidual = infNorm(r);
    const double bnorm = std::max(infNorm(bb), 1e-300);
    if (rep.initialResidual / bnorm < tol)
    {
        rep.converged = true;
        rep.finalResidual = rep.initialResidual;
        return rep;
    }

    std::vector<double> z(static_cast<std::size_t>(n)), p(static_cast<std::size_t>(n)), Ap(static_cast<std::size_t>(n));
    for (label i = 0; i < n; ++i) z[static_cast<std::size_t>(i)] = invD[static_cast<std::size_t>(i)] * r[static_cast<std::size_t>(i)];
    p = z;
    double rz = dot(r, z);

    for (int it = 1; it <= maxIter; ++it)
    {
        A.product(p, Ap);
        const double pAp = dot(p, Ap);
        if (pAp <= 0 || std::abs(rz) < 1e-300)
        {
            break; // matrix not sufficiently SPD; caller should use BiCGStab
        }
        const double alpha = rz / pAp;
        for (label i = 0; i < n; ++i)
        {
            const std::size_t u = static_cast<std::size_t>(i);
            xx[u] += alpha * p[u];
            r[u] -= alpha * Ap[u];
        }
        rep.iterations = it;
        rep.finalResidual = infNorm(r);
        if (rep.finalResidual / bnorm < tol)
        {
            rep.converged = true;
            break;
        }
        for (label i = 0; i < n; ++i) z[static_cast<std::size_t>(i)] = invD[static_cast<std::size_t>(i)] * r[static_cast<std::size_t>(i)];
        const double rzNew = dot(r, z);
        const double beta = rzNew / rz;
        rz = rzNew;
        for (label i = 0; i < n; ++i)
        {
            const std::size_t u = static_cast<std::size_t>(i);
            p[u] = z[u] + beta * p[u];
        }
    }

    x = Field<double>(xx.begin(), xx.end());
    return rep;
}

SolverReport solveBiCGStab(const SparseMatrix& A, const Field<double>& b,
                           Field<double>& x, double tol, int maxIter)
{
    const label n = A.size();
    std::vector<double> xx(x.begin(), x.end());
    const auto diag = A.diagonal();

    std::vector<double> invD(static_cast<std::size_t>(n));
    for (label i = 0; i < n; ++i)
    {
        const double d = diag[static_cast<std::size_t>(i)];
        invD[static_cast<std::size_t>(i)] = (d > 1e-300) ? 1.0 / d : 1.0;
    }

    std::vector<double> Ax(static_cast<std::size_t>(n), 0.0), r(static_cast<std::size_t>(n));
    A.product(xx, Ax);
    const std::vector<double> bb = bVec(b);
    for (label i = 0; i < n; ++i) r[static_cast<std::size_t>(i)] = bb[static_cast<std::size_t>(i)] - Ax[static_cast<std::size_t>(i)];

    SolverReport rep;
    rep.initialResidual = infNorm(r);
    const double bnorm = std::max(infNorm(bb), 1e-300);
    if (rep.initialResidual / bnorm < tol)
    {
        rep.converged = true;
        rep.finalResidual = rep.initialResidual;
        return rep;
    }

    const std::vector<double> r0 = r;
    // van der Vorst's preconditioned BiCGStab: p stays unconditioned,
    // v = A M^-1 p; z / sh are the transient M^-1 applications.
    std::vector<double> p(static_cast<std::size_t>(n)), v(static_cast<std::size_t>(n));
    std::vector<double> s(static_cast<std::size_t>(n)), z(static_cast<std::size_t>(n));
    std::vector<double> sh(static_cast<std::size_t>(n)), t(static_cast<std::size_t>(n));
    double rho = 1.0, alpha = 1.0, omega = 1.0;

    for (int it = 1; it <= maxIter; ++it)
    {
        const double rhoNew = dot(r0, r);
        if (std::abs(rhoNew) < 1e-300) break;
        const double beta = (rhoNew / rho) * (alpha / omega);
        rho = rhoNew;
        for (label i = 0; i < n; ++i)
        {
            const std::size_t u = static_cast<std::size_t>(i);
            p[u] = r[u] + beta * (p[u] - omega * v[u]);
            z[u] = p[u] * invD[u];
        }
        A.product(z, v);

        const double denom = dot(r0, v);
        if (std::abs(denom) < 1e-300) break;
        alpha = rho / denom;
        for (label i = 0; i < n; ++i)
        {
            const std::size_t u = static_cast<std::size_t>(i);
            xx[u] += alpha * z[u];
            s[u] = r[u] - alpha * v[u];
        }

        for (label i = 0; i < n; ++i)
        {
            const std::size_t u = static_cast<std::size_t>(i);
            sh[u] = s[u] * invD[u];
        }
        A.product(sh, t);
        const double tDenom = dot(t, t);
        if (tDenom < 1e-300)
        {
            r = s;
            rep.iterations = it;
            rep.finalResidual = infNorm(r);
            rep.converged = rep.finalResidual / bnorm < tol;
            break;
        }
        omega = dot(t, s) / tDenom;
        for (label i = 0; i < n; ++i)
        {
            const std::size_t u = static_cast<std::size_t>(i);
            xx[u] += omega * sh[u];
            r[u] = s[u] - omega * t[u];
        }

        rep.iterations = it;
        rep.finalResidual = infNorm(r);
        if (rep.finalResidual / bnorm < tol)
        {
            rep.converged = true;
            break;
        }
    }

    x = Field<double>(xx.begin(), xx.end());
    return rep;
}

const char* eigenBackendName()
{
#ifdef CHFLOW_WITH_EIGEN
    return "eigen";
#else
    return "none";
#endif
}

SolverReport solveWithBackend(const SparseMatrix& A, const Field<double>& b,
                              Field<double>& x, double tol, int maxIter,
                              LinearBackend backend)
{
    if (backend == LinearBackend::Native)
    {
        // Preserve the exact native behaviour: a symmetric-marked matrix is
        // not known here, so match fvMatrix's default (BiCGStab, which handles
        // both symmetric and non-symmetric systems) unless the caller picked
        // CG through the regular entry points.
        return solveBiCGStab(A, b, x, tol, maxIter);
    }

    if (backend == LinearBackend::Eigen)
    {
#ifdef CHFLOW_WITH_EIGEN
    const label n = A.size();
    Eigen::SparseMatrix<double, Eigen::RowMajor> mat(n, n);
    using Triplet = Eigen::Triplet<double>;
    std::vector<Triplet> triplets;
    triplets.reserve(static_cast<std::size_t>(A.vals().size())
                     + static_cast<std::size_t>(n));
    for (label r = 0; r < n; ++r)
    {
        const label b0 = A.rowStart()[static_cast<std::size_t>(r)];
        const label b1 = A.rowStart()[static_cast<std::size_t>(r + 1)];
        for (label k = b0; k < b1; ++k)
            triplets.emplace_back(r, A.cols()[static_cast<std::size_t>(k)],
                                  A.vals()[static_cast<std::size_t>(k)]);
    }
    mat.setFromTriplets(triplets.begin(), triplets.end());
    mat.makeCompressed();

    Eigen::VectorXd rhs(n), sol(n);
    const std::vector<double> bb = bVec(b);
    for (label i = 0; i < n; ++i) rhs[static_cast<Eigen::Index>(i)] = bb[static_cast<std::size_t>(i)];
    for (label i = 0; i < n; ++i)
        sol[static_cast<Eigen::Index>(i)] = (i < static_cast<label>(x.size()))
            ? x[static_cast<std::size_t>(i)] : 0.0;

    Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>, Eigen::IncompleteLUT<double>> solver;
    solver.setTolerance(tol);
    solver.setMaxIterations(maxIter);
    solver.compute(mat);
    sol = solver.solve(rhs);

    std::vector<double> xx(static_cast<std::size_t>(n));
    for (label i = 0; i < n; ++i) xx[static_cast<std::size_t>(i)] = sol[static_cast<Eigen::Index>(i)];
    // Measure the residual of the returned solution against the same assembled
    // CSR the native path solved (Eigen's own `error()` is a relative estimate
    // in its own norm, not the inf-norm the SolverReport contract reports).
    std::vector<double> ax(static_cast<std::size_t>(n), 0.0);
    A.product(xx, ax);
    double residual = 0.0, bnorm = 0.0;
    for (label i = 0; i < n; ++i)
    {
        const std::size_t u = static_cast<std::size_t>(i);
        residual = std::max(residual, std::abs(bb[u] - ax[u]));
        bnorm = std::max(bnorm, std::abs(bb[u]));
    }
    x = Field<double>(xx.begin(), xx.end());
    SolverReport rep;
    rep.iterations = solver.iterations();
    rep.finalResidual = residual;
    rep.initialResidual = bnorm;
    rep.converged = (solver.info() == Eigen::Success)
                    && (residual / std::max(bnorm, 1e-300) < tol);
    return rep;
#else
    cfFatal("solveWithBackend: Eigen backend requested but CHFLOW_WITH_EIGEN "
            "was not compiled in");
#endif
    }

    if (backend == LinearBackend::PETSc)
    {
#ifdef CHFLOW_WITH_PETSC
        // PETSc's sequential AIJ matrix + GMRES with ILU (or Jacobi when ILU
        // is unavailable) solves the identical assembled CSR system. PETSc is
        // initialised once and left alive for the process lifetime; a solver
        // app is responsible for calling PetscInitialize before reaching here.
        const label n = A.size();
        const PetscInt pn = static_cast<PetscInt>(n);
        Mat M;
        PetscErrorCode ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, pn, pn, 0, nullptr, &M);
        CHKERRABORT(PETSC_COMM_SELF, ierr);
        // Row nnz counts, then the CSR values/colindices.
        // Assemble through explicit (row, col, value) triplets, which is
        // portable across PETSc versions and matches the CSR exactly.
        for (label r = 0; r < n; ++r)
        {
            const label b0 = A.rowStart()[static_cast<std::size_t>(r)];
            const label b1 = A.rowStart()[static_cast<std::size_t>(r + 1)];
            for (label k = b0; k < b1; ++k)
            {
                const PetscInt col = static_cast<PetscInt>(A.cols()[static_cast<std::size_t>(k)]);
                const PetscScalar val = A.vals()[static_cast<std::size_t>(k)];
                ierr = MatSetValues(M, 1, &r, 1, &col, &val, ADD_VALUES);
                CHKERRABORT(PETSC_COMM_SELF, ierr);
            }
        }
        ierr = MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_SELF, ierr);
        ierr = MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_SELF, ierr);

        Vec rhs, sol;
        ierr = MatCreateVecs(M, &rhs, &sol); CHKERRABORT(PETSC_COMM_SELF, ierr);
        const std::vector<double> bb = bVec(b);
        for (label i = 0; i < n; ++i)
        {
            const PetscScalar v = bb[static_cast<std::size_t>(i)];
            ierr = VecSetValues(rhs, 1, &i, &v, INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF, ierr);
            ierr = VecSetValue(sol, static_cast<PetscInt>(i),
                               (i < static_cast<label>(x.size())) ? x[static_cast<std::size_t>(i)] : 0.0,
                               INSERT_VALUES);
            CHKERRABORT(PETSC_COMM_SELF, ierr);
        }
        ierr = VecAssemblyBegin(rhs); CHKERRABORT(PETSC_COMM_SELF, ierr);
        ierr = VecAssemblyEnd(rhs); CHKERRABORT(PETSC_COMM_SELF, ierr);
        ierr = VecAssemblyBegin(sol); CHKERRABORT(PETSC_COMM_SELF, ierr);
        ierr = VecAssemblyEnd(sol); CHKERRABORT(PETSC_COMM_SELF, ierr);

        KSP ksp;
        ierr = KSPCreate(PETSC_COMM_SELF, &ksp); CHKERRABORT(PETSC_COMM_SELF, ierr);
        ierr = KSPSetOperators(ksp, M, M); CHKERRABORT(PETSC_COMM_SELF, ierr);
        ierr = KSPSetType(ksp, KSPGMRES); CHKERRABORT(PETSC_COMM_SELF, ierr);
        PC pc;
        ierr = KSPGetPC(ksp, &pc); CHKERRABORT(PETSC_COMM_SELF, ierr);
        ierr = PCSetType(pc, PCILU); CHKERRABORT(PETSC_COMM_SELF, ierr);
        ierr = KSPSetTolerances(ksp, static_cast<PetscReal>(tol), PETSC_DEFAULT,
                                PETSC_DEFAULT, static_cast<PetscInt>(maxIter));
        CHKERRABORT(PETSC_COMM_SELF, ierr);
        ierr = KSPSolve(ksp, rhs, sol); CHKERRABORT(PETSC_COMM_SELF, ierr);

        PetscInt its = 0;
        KSPConvergedReason reason;
        ierr = KSPGetIterationNumber(ksp, &its); CHKERRABORT(PETSC_COMM_SELF, ierr);
        ierr = KSPGetConvergedReason(ksp, &reason); CHKERRABORT(PETSC_COMM_SELF, ierr);

        const PetscScalar* solArr = nullptr;
        ierr = VecGetArrayRead(sol, &solArr); CHKERRABORT(PETSC_COMM_SELF, ierr);
        std::vector<double> xx(static_cast<std::size_t>(n));
        for (label i = 0; i < n; ++i) xx[static_cast<std::size_t>(i)] = solArr[i];
        ierr = VecRestoreArrayRead(sol, &solArr); CHKERRABORT(PETSC_COMM_SELF, ierr);

        std::vector<double> ax(static_cast<std::size_t>(n), 0.0);
        A.product(xx, ax);
        double residual = 0.0, bnorm = 0.0;
        for (label i = 0; i < n; ++i)
        {
            const std::size_t u = static_cast<std::size_t>(i);
            residual = std::max(residual, std::abs(bb[u] - ax[u]));
            bnorm = std::max(bnorm, std::abs(bb[u]));
        }
        x = Field<double>(xx.begin(), xx.end());
        SolverReport rep;
        rep.iterations = static_cast<int>(its);
        rep.finalResidual = residual;
        rep.initialResidual = bnorm;
        rep.converged = (reason > 0)
                        && (residual / std::max(bnorm, 1e-300) < tol);

        ierr = KSPDestroy(&ksp); CHKERRABORT(PETSC_COMM_SELF, ierr);
        ierr = VecDestroy(&rhs); CHKERRABORT(PETSC_COMM_SELF, ierr);
        ierr = VecDestroy(&sol); CHKERRABORT(PETSC_COMM_SELF, ierr);
        ierr = MatDestroy(&M); CHKERRABORT(PETSC_COMM_SELF, ierr);
        return rep;
#else
        cfFatal("solveWithBackend: PETSc backend requested but CHFLOW_WITH_PETSC "
                "was not compiled in");
#endif
    }

    cfFatal("solveWithBackend: unknown backend requested");
}

const char* petscBackendName()
{
#ifdef CHFLOW_WITH_PETSC
    return "petsc";
#else
    return "none";
#endif
}

} // namespace cf
