#include "cf/turbulence/RsmModel.h"

#include <algorithm>
#include <cmath>

#include "cf/common/Dictionary.h"
#include "cf/common/Error.h"
#include "gradU.h"

namespace cf
{
namespace
{

// stress order: 0=11 1=22 2=33 3=12 4=13 5=23
constexpr int ci[6] = {0, 1, 2, 0, 0, 1};
constexpr int cj[6] = {0, 1, 2, 1, 2, 2};

} // namespace

RsmModel::Config RsmModel::Config::fromCase(const Dictionary& caseDict)
{
    Config c;
    if (!caseDict.has("rsm")) return c;
    const Dictionary& d = caseDict.getDict("rsm");
    c.C1 = d.getScalarOr("C1", c.C1);
    c.C2 = d.getScalarOr("C2", c.C2);
    c.C1eps = d.getScalarOr("C1eps", c.C1eps);
    c.C2eps = d.getScalarOr("C2eps", c.C2eps);
    c.sigmaR = d.getScalarOr("sigmaR", c.sigmaR);
    c.sigmaEps = d.getScalarOr("sigmaEps", c.sigmaEps);
    c.relaxation = d.getScalarOr("relaxation", c.relaxation);
    c.k0 = d.getScalarOr("k0", c.k0);
    c.eps0 = d.getScalarOr("eps0", c.eps0);
    c.validate();
    return c;
}

void RsmModel::Config::validate() const
{
    if (!(C1 > 0)) cfFatal("rsm: C1 must be positive");
    if (C2 != 0) cfFatal("rsm: only the slow LRR regression is implemented (C2 must be 0)");
    if (!(C1eps > 0) || !(C2eps > 0)) cfFatal("rsm: eps constants must be positive");
    if (!(sigmaR > 0) || !(sigmaEps > 0)) cfFatal("rsm: sigma must be positive");
    if (!(relaxation > 0 && relaxation <= 1)) cfFatal("rsm: relaxation must be in (0, 1]");
    if (!(k0 > 0) || !(eps0 > 0)) cfFatal("rsm: k0 and eps0 must be positive");
}

RsmModel::RsmModel(const polyMesh& mesh, volVectorField& U, const surfaceScalarField& phi,
                   scalar nu, const Schemes& schemes, const Config& cfg)
    : mesh_(mesh), U_(U), phi_(phi), nu_(nu), schemes_(schemes), cfg_(cfg),
      R_{volScalarField(mesh, "R11", 0.0), volScalarField(mesh, "R22", 0.0),
         volScalarField(mesh, "R33", 0.0), volScalarField(mesh, "R12", 0.0),
         volScalarField(mesh, "R13", 0.0), volScalarField(mesh, "R23", 0.0)},
      eps_(mesh, "epsilon", cfg.eps0),
      k_(static_cast<std::size_t>(mesh.nCells()), 0.0),
      nut_(static_cast<std::size_t>(mesh.nCells()), 0.0),
      nuEffFace_(static_cast<std::size_t>(mesh.nFaces()), nu),
      pk_(static_cast<std::size_t>(mesh.nCells()), 0.0)
{
    if (!(nu_ > 0)) cfFatal("rsm: molecular viscosity must be positive");
    cfg_.validate();

    // No wall functions in this closure: every patch keeps a zero normal
    // gradient for every stress, and the dissipation too. The cases that use
    // the model are homogeneous-shear columns and decaying turbulence, where
    // the wall layer of the two-equation models would only add a second,
    // undeclared closure.
    for (const polyMesh::Patch& p : mesh_.patches())
    {
        const label patch = mesh_.patchIdByName(p.name);
        for (int c = 0; c < 6; ++c) R_[c].setPatchValue(patch, "zeroGradient", 0.0);
        eps_.setPatchValue(patch, "zeroGradient", 0.0);
    }

    // isotropic start: R_ii = 2 k0 / 3, off-diagonals zero
    for (int c = 0; c < 6; ++c)
    {
        const scalar v = c < 3 ? 2.0 * cfg_.k0 / 3.0 : 0.0;
        for (scalar& x : R_[c].internalField()) x = v;
    }
    for (scalar& x : eps_.internalField()) x = cfg_.eps0;

    updateNut();
    refreshFaceNut();
}

void RsmModel::correctTurbulence(scalar dt)
{
    Field<Vector> gradU[3];
    turb_detail::computeGradU(mesh_, U_, schemes_, gradU);

    const std::size_t nC = static_cast<std::size_t>(mesh_.nCells());

    // Production of every stress from the current Reynolds stress and mean
    // gradient: P_ij = -(R_jk dU_i/dx_k + R_ik dU_j/dx_k).
    Field<scalar> prod[6];
    for (int c = 0; c < 6; ++c) prod[c].assign(nC, 0.0);
    for (std::size_t c = 0; c < nC; ++c)
    {
        scalar R[3][3];
        for (int m = 0; m < 6; ++m)
        {
            const scalar v = R_[m].internalField()[c];
            R[ci[m]][cj[m]] = v;
            R[cj[m]][ci[m]] = v;
        }
        for (int m = 0; m < 6; ++m)
        {
            const int i = ci[m], j = cj[m];
            scalar p = 0;
            for (int k = 0; k < 3; ++k)
            {
                p += -(R[j][k] * gradU[i][c][k] + R[i][k] * gradU[j][c][k]);
            }
            prod[m][c] = p;
        }
        scalar S[3][3];
        scalar magS = 0;
        turb_detail::strainAt(gradU, c, S, magS);
        scalar rs = 0;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) rs += R[i][j] * S[i][j];
        pk_[c] = -rs; // P_k = -R_ij S_ij
    }

    // k from the trace of the current stresses
    for (std::size_t c = 0; c < nC; ++c)
    {
        k_[c] = 0.5 * (R_[0].internalField()[c] + R_[1].internalField()[c]
                     + R_[2].internalField()[c]);
        k_[c] = std::max(k_[c], kFloor);
    }

    for (int m = 0; m < 6; ++m) solveStress(m, dt, prod, k_);
    solveEpsilon(dt);

    // rebuild k after the transport, clamp realizability, then nut
    for (std::size_t c = 0; c < nC; ++c)
    {
        for (int m = 0; m < 3; ++m)
        {
            scalar& v = R_[m].internalField()[c];
            v = std::max(v, 0.0);
        }
        for (int m = 3; m < 6; ++m)
        {
            const scalar b = std::sqrt(R_[ci[m]].internalField()[c]
                                     * R_[cj[m]].internalField()[c]);
            scalar& v = R_[m].internalField()[c];
            v = std::min(std::max(v, -b), b);
        }
        k_[c] = std::max(0.5 * (R_[0].internalField()[c] + R_[1].internalField()[c]
                              + R_[2].internalField()[c]), kFloor);
    }

    updateNut();
    refreshFaceNut();
}

void RsmModel::solveStress(int comp, scalar dt, const Field<scalar>* prod,
                           const Field<scalar>& kCell)
{
    volScalarField& R = R_[comp];
    const Field<scalar> old = R.internalField();
    const auto& vol = mesh_.cellVolumes();
    const int i = ci[comp], j = cj[comp];
    const scalar delta = i == j ? 1.0 : 0.0;

    fvMatrix eqn(mesh_, viewOf(R), "REqn");
    eqn.schemes() = schemes_;
    eqn.schemes().relaxation = 1.0;
    if (dt > 0) fvm::ddt(eqn, old, dt);
    fvm::div(eqn, phi_.values(), schemes_.div);

    const std::size_t nC = static_cast<std::size_t>(mesh_.nCells());
    Field<scalar> gamma(nC);
    Field<scalar> sp(nC, 0.0);
    Field<scalar> rhs(nC, 0.0);
    for (std::size_t c = 0; c < nC; ++c)
    {
        gamma[c] = nu_ + nut_[c] / cfg_.sigmaR;
        const scalar e = std::max(eps_.internalField()[c], dissFloor);
        const scalar coef = cfg_.C1 * e / kCell[c]; // (eps/k) C1
        sp[c] = vol[c] * coef;                      // implicit on R_ij
        // explicit remainder: P_ij + (2/3) eps (C1 - 1) delta_ij
        rhs[c] = vol[c] * (prod[comp][c] + (2.0 / 3.0) * e * (cfg_.C1 - 1.0) * delta);
    }
    fvm::laplacianFace(eqn, harmonicFaceGamma(mesh_, gamma), schemes_.laplacian, -1);
    fvm::sp(eqn, sp, 1);
    for (label c = 0; c < mesh_.nCells(); ++c) eqn.addRhs(c, 0, rhs[static_cast<std::size_t>(c)]);

    eqn.solve(fvMatrix::Controls{linearTol, linearMaxIter, 2});

    const Field<scalar>& x = eqn.solution()[0];
    const scalar rl = cfg_.relaxation;
    for (std::size_t c = 0; c < nC; ++c)
    {
        R.internalField()[c] = old[c] + rl * (x[c] - old[c]);
    }
}

void RsmModel::solveEpsilon(scalar dt)
{
    const Field<scalar> old = eps_.internalField();
    const auto& vol = mesh_.cellVolumes();
    const std::size_t nC = static_cast<std::size_t>(mesh_.nCells());

    fvMatrix eqn(mesh_, viewOf(eps_), "epsilonEqn");
    eqn.schemes() = schemes_;
    eqn.schemes().relaxation = 1.0;
    if (dt > 0) fvm::ddt(eqn, old, dt);
    fvm::div(eqn, phi_.values(), schemes_.div);

    Field<scalar> gamma(nC), sp(nC, 0.0), rhs(nC, 0.0);
    for (std::size_t c = 0; c < nC; ++c)
    {
        gamma[c] = nu_ + nut_[c] / cfg_.sigmaEps;
        const scalar e = std::max(old[c], dissFloor);
        const scalar k = k_[c];
        sp[c] = vol[c] * cfg_.C2eps * e / k;
        rhs[c] = vol[c] * cfg_.C1eps * e * pk_[c] / k;
    }
    fvm::laplacianFace(eqn, harmonicFaceGamma(mesh_, gamma), schemes_.laplacian, -1);
    fvm::sp(eqn, sp, 1);
    for (label c = 0; c < mesh_.nCells(); ++c) eqn.addRhs(c, 0, rhs[static_cast<std::size_t>(c)]);

    eqn.solve(fvMatrix::Controls{linearTol, linearMaxIter, 2});

    const Field<scalar>& x = eqn.solution()[0];
    const scalar rl = cfg_.relaxation;
    for (std::size_t c = 0; c < nC; ++c)
    {
        eps_.internalField()[c] = std::max(old[c] + rl * (x[c] - old[c]), dissFloor);
    }
}

void RsmModel::updateNut()
{
    Field<Vector> gradU[3];
    turb_detail::computeGradU(mesh_, U_, schemes_, gradU);
    for (std::size_t c = 0; c < nut_.size(); ++c)
    {
        scalar S[3][3];
        scalar magS = 0;
        turb_detail::strainAt(gradU, c, S, magS);
        scalar s2 = 0;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) s2 += S[i][j] * S[i][j];
        // nut |S|^2 = -R:S = P_k keeps the modelled shear stress exact
        nut_[c] = s2 > 0 ? std::max(pk_[c], 0.0) / s2 : 0.0;
    }
}

void RsmModel::refreshFaceNut()
{
    const Field<scalar> nutFace = harmonicFaceGamma(mesh_, nut_);
    for (std::size_t f = 0; f < nuEffFace_.size(); ++f)
    {
        nuEffFace_[f] = nu_ + nutFace[f];
    }
}

} // namespace cf
