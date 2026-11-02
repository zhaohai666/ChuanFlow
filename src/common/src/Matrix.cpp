#include "cf/common/Matrix.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

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

} // namespace cf
