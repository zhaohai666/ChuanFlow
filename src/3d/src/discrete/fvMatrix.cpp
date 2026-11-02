#include "cf/fvm/FvMatrix.h"

#include <algorithm>
#include <cmath>

#include "cf/common/Error.h"
#include "cf/common/Log.h"

namespace cf
{

FieldView viewOf(volScalarField& f)
{
    FieldView v;
    v.nComp = 1;
    volScalarField* p = &f;
    v.cell = [p](label c, int) { return p->internalField()[static_cast<std::size_t>(c)]; };
    v.lin = [p](label face, int, scalar& coeff, scalar& konst) { p->linearise(face, coeff, konst); };
    return v;
}

FieldView viewOf(volVectorField& f)
{
    FieldView v;
    v.nComp = 3;
    volVectorField* p = &f;
    v.cell = [p](label c, int comp) { return p->internalField()[static_cast<std::size_t>(c)][comp]; };
    v.lin = [p](label face, int comp, scalar& coeff, scalar& konst)
    { p->linearise(face, comp, coeff, konst); };
    return v;
}

fvMatrix::fvMatrix(const polyMesh& mesh, FieldView unknown, std::string name)
    : mesh_(mesh), unknown_(std::move(unknown)), name_(std::move(name))
{
    const std::size_t nc = static_cast<std::size_t>(mesh.nCells());
    for (int c = 0; c < unknown_.nComp; ++c)
    {
        builders_.emplace_back(mesh.nCells());
        diag_.emplace_back(nc, 0.0);
        rhs_.emplace_back(nc, 0.0);
        x_.emplace_back(nc, 0.0);
        matrices_.emplace_back();
    }
}

void fvMatrix::addDiag(label cell, int comp, scalar v)
{
    diag_[static_cast<std::size_t>(comp)][static_cast<std::size_t>(cell)] += v;
}

void fvMatrix::addRhs(label cell, int comp, scalar v)
{
    rhs_[static_cast<std::size_t>(comp)][static_cast<std::size_t>(cell)] += v;
}

void fvMatrix::addCorrection(label cell, int comp, scalar v)
{
    rhs_[static_cast<std::size_t>(comp)][static_cast<std::size_t>(cell)] += v;
}

void fvMatrix::addOffDiag(label row, label col, int comp, scalar v)
{
    builders_[static_cast<std::size_t>(comp)].add(row, col, v);
}

void fvMatrix::registerCorrector(
    std::function<void(fvMatrix&, const std::vector<Field<scalar>>&)> fn)
{
    correctors_.push_back(std::move(fn));
}

SolverReport fvMatrix::solve(const Controls& ctl)
{
    baseRhs_ = rhs_;

    std::vector<double> diagVec;
    for (int comp = 0; comp < unknown_.nComp; ++comp)
    {
        const auto& d = diag_[static_cast<std::size_t>(comp)];
        diagVec.assign(d.begin(), d.end());
        matrices_[static_cast<std::size_t>(comp)] =
            builders_[static_cast<std::size_t>(comp)].finishWithDiagonal(diagVec);
    }

    const int passes = correctors_.empty() ? 1 : std::max(1, ctl.nonOrth + 1);
    SolverReport last;
    for (int pass = 0; pass < passes; ++pass)
    {
        rhs_ = baseRhs_;
        for (auto& corr : correctors_) corr(*this, x_);

        const bool finalPass = (pass + 1 == passes);
        for (int comp = 0; comp < unknown_.nComp; ++comp)
        {
            Field<scalar> xOld = x_[static_cast<std::size_t>(comp)];
            SolverReport rep;
            if (preferCG_)
            {
                rep = solveCG(matrices_[static_cast<std::size_t>(comp)],
                              rhs_[static_cast<std::size_t>(comp)],
                              x_[static_cast<std::size_t>(comp)], ctl.tol, ctl.maxIter);
            }
            else
            {
                rep = solveBiCGStab(matrices_[static_cast<std::size_t>(comp)],
                                    rhs_[static_cast<std::size_t>(comp)],
                                    x_[static_cast<std::size_t>(comp)], ctl.tol, ctl.maxIter);
            }
            if (finalPass && schemes_.relaxation < 1.0)
            {
                const scalar rl = schemes_.relaxation;
                auto& x = x_[static_cast<std::size_t>(comp)];
                for (std::size_t i = 0; i < x.size(); ++i)
                {
                    x[i] = xOld[i] + rl * (x[i] - xOld[i]);
                }
            }
            last = rep;
        }
    }

    if (!last.converged)
    {
        log::warn("fvMatrix '", name_, "': linear solver did not converge (",
                  last.iterations, " iters, residual ", last.finalResidual, ")");
    }
    return last;
}

void fvMatrix::explicitOperator(const Field<scalar>& xIn, int comp, Field<scalar>& h) const
{
    // OpenFOAM's H(): source minus off-diagonal product. The diagonal part
    // is subtracted from A*x so that x = rAU * H is the fixed-point solution
    // of A x = rhs.
    h = rhs_[static_cast<std::size_t>(comp)];
    const auto& A = matrices_[static_cast<std::size_t>(comp)];
    std::vector<double> ax;
    A.product(std::vector<double>(xIn.begin(), xIn.end()), ax);
    const auto& d = diag_[static_cast<std::size_t>(comp)];
    for (label c = 0; c < mesh_.nCells(); ++c)
    {
        const std::size_t u = static_cast<std::size_t>(c);
        h[u] -= ax[u] - d[u] * xIn[u];
    }
}

Schemes loadSchemes(const Dictionary& caseDict)
{
    Schemes s;
    if (!caseDict.has("schemes")) return s;
    const Dictionary& d = caseDict.getDict("schemes");
    s.ddt = d.getWordOr("ddt", s.ddt);
    s.div = d.getWordOr("div", s.div);
    s.laplacian = d.getWordOr("laplacian", s.laplacian);
    s.grad = d.getWordOr("grad", s.grad);
    s.blendCoeff = d.getScalarOr("blendCoeff", s.blendCoeff);
    s.relaxation = d.getScalarOr("relaxation", s.relaxation);
    return s;
}

} // namespace cf
