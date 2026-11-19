// chflow-solver-backend: exercise the selectable sparse linear solver backends
// (README M6 / A11) on the same assembled CSR system with a known analytic
// solution. The 1D Poisson-like tridiagonal A (diag 2, off-diag -1) with
// x_exact[i] = 1 has b = A x_exact, so the exact solution is known to machine
// precision. Two backends solve it from a zero initial guess:
//
//   * Native -- the project's own diagonal-preconditioned BiCGStab
//   * Eigen  -- Eigen's ILU-preconditioned BiCGSTAB (when compiled in)
//
// The app prints each backend's error against the analytic solution, its
// assembled-CSR residual, iteration count and convergence flag. A checker
// verifies both backends land on the same physics; within a build the native
// path stays bit-identical to before, which the rest of the suite witnesses.
//
// usage: chflow-solver-backend <case.dict>
#include "cf/common/Dictionary.h"
#include "cf/common/Error.h"
#include "cf/common/Field.h"
#include "cf/common/Log.h"
#include "cf/common/Matrix.h"

#ifdef CHFLOW_WITH_PETSC
#include <petscsys.h>
#endif

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace
{
std::string num(double v)
{
    std::ostringstream os;
    os << std::setprecision(17) << v;
    return os.str();
}

struct SolveResult
{
    std::vector<double> x;
    cf::SolverReport rep;
    double err = 0;      // inf-norm error vs analytic solution
    double residual = 0; // inf-norm of A x - b
};

SolveResult run(const cf::SparseMatrix& A, const std::vector<double>& b,
                const std::vector<double>& exact, double tol, int maxIter,
                cf::LinearBackend backend)
{
    cf::Field<double> x(static_cast<std::size_t>(A.size()), 0.0);
    const cf::Field<double> rhs(b.begin(), b.end());
    const cf::SolverReport rep =
        cf::solveWithBackend(A, rhs, x, tol, maxIter, backend);
    SolveResult out;
    out.x.assign(x.begin(), x.end());
    out.rep = rep;
    for (cf::label i = 0; i < A.size(); ++i)
        out.err = std::max(out.err, std::abs(out.x[static_cast<std::size_t>(i)]
                                             - exact[static_cast<std::size_t>(i)]));
    std::vector<double> ax(static_cast<std::size_t>(A.size()), 0.0);
    A.product(out.x, ax);
    for (cf::label i = 0; i < A.size(); ++i)
        out.residual = std::max(out.residual,
            std::abs(ax[static_cast<std::size_t>(i)] - b[static_cast<std::size_t>(i)]));
    return out;
}
} // namespace

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        std::cerr << "usage: chflow-solver-backend <case.dict>\n";
        return 2;
    }
    try
    {
        const cf::Dictionary dict = cf::Dictionary::parseFile(argv[1]);
        const cf::label n = static_cast<cf::label>(dict.getScalar("n"));
        const double tol = dict.getScalarOr("tol", 1e-10);
        const int maxIter = static_cast<int>(dict.getScalarOr("maxIter", 500.0));
        if (!(n > 1)) cfFatal("solver-backend: n must exceed 1 (got " + std::to_string(n) + ")");

        // Assemble the tridiagonal SPD system and its analytic right-hand side.
        cf::SparseMatrix::Builder builder(n);
        for (cf::label i = 0; i < n; ++i)
        {
            builder.add(i, i, 2.0);
            if (i > 0) builder.add(i, i - 1, -1.0);
            if (i + 1 < n) builder.add(i, i + 1, -1.0);
        }
        const cf::SparseMatrix A = builder.finish();
        std::vector<double> exact(static_cast<std::size_t>(n), 1.0);
        std::vector<double> b(static_cast<std::size_t>(n), 0.0);
        A.product(exact, b);

        cf::log::info("solver-backend n ", n, " tol ", num(tol),
                      " eigen ", cf::eigenBackendName(),
                      " petsc ", cf::petscBackendName());

        const SolveResult native = run(A, b, exact, tol, maxIter, cf::LinearBackend::Native);
        cf::log::info("solver-backend native converged ", native.rep.converged ? 1 : 0,
                      " iterations ", native.rep.iterations,
                      " err ", num(native.err), " residual ", num(native.residual));

        if (std::string(cf::eigenBackendName()) == "eigen")
        {
            const SolveResult eigen = run(A, b, exact, tol, maxIter, cf::LinearBackend::Eigen);
            cf::log::info("solver-backend eigen converged ", eigen.rep.converged ? 1 : 0,
                          " iterations ", eigen.rep.iterations,
                          " err ", num(eigen.err), " residual ", num(eigen.residual));
            double cross = 0;
            for (cf::label i = 0; i < n; ++i)
                cross = std::max(cross, std::abs(eigen.x[static_cast<std::size_t>(i)]
                                                  - native.x[static_cast<std::size_t>(i)]));
            cf::log::info("solver-backend crossBackendMaxDiff ", num(cross));
            if (!eigen.rep.converged || !(eigen.err < 1e-6) || !(cross < 1e-6))
                cfFatal("solver-backend: Eigen backend did not match the analytic/native solution");
        }

        if (std::string(cf::petscBackendName()) == "petsc")
        {
#ifdef CHFLOW_WITH_PETSC
            PetscInitialize(&argc, &argv, nullptr, nullptr);
#endif
            const SolveResult petsc = run(A, b, exact, tol, maxIter, cf::LinearBackend::PETSc);
            cf::log::info("solver-backend petsc converged ", petsc.rep.converged ? 1 : 0,
                          " iterations ", petsc.rep.iterations,
                          " err ", num(petsc.err), " residual ", num(petsc.residual));
            double cross = 0;
            for (cf::label i = 0; i < n; ++i)
                cross = std::max(cross, std::abs(petsc.x[static_cast<std::size_t>(i)]
                                                  - native.x[static_cast<std::size_t>(i)]));
            cf::log::info("solver-backend petscCrossMaxDiff ", num(cross));
#ifdef CHFLOW_WITH_PETSC
            PetscFinalize();
#endif
            if (!petsc.rep.converged || !(petsc.err < 1e-6) || !(cross < 1e-6))
                cfFatal("solver-backend: PETSc backend did not match the analytic/native solution");
        }

        if (!native.rep.converged || !(native.err < 1e-6))
            cfFatal("solver-backend: native backend did not match the analytic solution");
        return 0;
    }
    catch (const std::exception& e)
    {
        cf::log::error(e.what());
        return 1;
    }
}
