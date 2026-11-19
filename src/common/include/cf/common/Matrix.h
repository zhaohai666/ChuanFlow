#pragma once

#include <vector>

#include "cf/common/Field.h"
#include "cf/common/Labels.h"

namespace cf
{

/// Row-compressed sparse matrix assembled from COO triplets.
/// Serves as the linear system behind fvMatrix (OpenFOAM's lduMatrix role).
class SparseMatrix
{
public:
    class Builder
    {
    public:
        explicit Builder(label n) : n_(n) {}

        void add(label row, label col, double value)
        {
            rows_.push_back(row);
            cols_.push_back(col);
            vals_.push_back(value);
        }

        SparseMatrix finish() const;

        /// finish() plus implicit diagonal entries (r,r)=diag[r].
        SparseMatrix finishWithDiagonal(const std::vector<double>& diag) const;

    private:
        label n_;
        std::vector<label> rows_;
        std::vector<label> cols_;
        std::vector<double> vals_;
    };

    SparseMatrix() = default;

    label size() const { return static_cast<label>(rowStart_.size()) - 1; }
    const std::vector<label>& rowStart() const { return rowStart_; }
    const std::vector<label>& cols() const { return cols_; }
    const std::vector<double>& vals() const { return vals_; }

    void product(const std::vector<double>& x, std::vector<double>& y) const;

    std::vector<double> diagonal() const;

private:
    std::vector<label> rowStart_;
    std::vector<label> cols_;
    std::vector<double> vals_;
};

struct SolverReport
{
    int iterations = 0;
    double initialResidual = 0;
    double finalResidual = 0;
    bool converged = false;
};

/// Symmetric part: diagonal-preconditioned CG (stands in for OpenFOAM's DICG).
/// Falls back to unpreconditioned directions where the diagonal is unusable.
SolverReport solveCG(const SparseMatrix& A, const Field<double>& b,
                     Field<double>& x, double tol, int maxIter);

/// Non-symmetric: diagonal-preconditioned BiCGStab (stands in for PBiCG).
SolverReport solveBiCGStab(const SparseMatrix& A, const Field<double>& b,
                           Field<double>& x, double tol, int maxIter);

/// Selectable sparse linear solver backend (README M6 / A11). `Native` is the
/// project's own CG/BiCGStab above and is the default everywhere, so every
/// existing checker stays bit-identical. `Eigen` dispatches to Eigen's
/// incomplete-LU-preconditioned BiCGSTAB when the library is compiled in
/// (CHFLOW_WITH_EIGEN); requesting an unavailable backend fails loudly rather
/// than silently changing the numerical result.
enum class LinearBackend { Native, Eigen, PETSc };

/// Solve A x = b with an explicitly chosen backend. Both backends solve the
/// same assembled CSR system; they agree to the solver tolerance (they are
/// different iterative algorithms, so bit-identity is not claimed across
/// backends -- only within a backend and against the native path when Eigen
/// is not compiled in).
SolverReport solveWithBackend(const SparseMatrix& A, const Field<double>& b,
                              Field<double>& x, double tol, int maxIter,
                              LinearBackend backend);

/// Name of the compiled-in Eigen backend, or "none" when CHFLOW_WITH_EIGEN
/// was off. Useful for a checker to witness which path actually ran.
const char* eigenBackendName();

/// Name of the compiled-in PETSc backend, or "none" when CHFLOW_WITH_PETSc
/// was off / no PETSc library was found at configure time.
const char* petscBackendName();

} // namespace cf
