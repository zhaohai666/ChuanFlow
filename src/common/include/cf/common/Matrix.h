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

} // namespace cf
