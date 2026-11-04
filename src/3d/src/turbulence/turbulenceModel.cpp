#include "cf/turbulence/TurbulenceModel.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <queue>
#include <utility>

#include "cf/common/Error.h"
#include "cf/common/Log.h"
#include "cf/turbulence/KEpsilon.h"
#include "cf/turbulence/KOmegaSST.h"

namespace cf
{
namespace
{

Field<scalar> component(const Field<Vector>& u, int comp)
{
    Field<scalar> out(static_cast<std::size_t>(u.size()));
    for (std::size_t i = 0; i < u.size(); ++i) out[i] = u[i][comp];
    return out;
}

bool isWallPatch(const polyMesh::Patch& p) { return p.type == "wall"; }
bool isInletPatch(const polyMesh::Patch& p) { return p.type == "inlet"; }

} // namespace

void TurbulenceConfig::validate() const
{
    if (!(intensity >= 0)) cfFatal("turbulence: intensity must be non-negative");
    if (!(relaxation > 0 && relaxation <= 1))
    {
        cfFatal("turbulence: relaxation must be in (0, 1]");
    }
}

TurbulenceConfig TurbulenceConfig::fromCase(const Dictionary& caseDict)
{
    TurbulenceConfig cfg;
    if (!caseDict.has("turbulence")) return cfg;
    const Dictionary& td = caseDict.getDict("turbulence");
    cfg.model = td.getWordOr("model", cfg.model);
    cfg.intensity = td.getScalarOr("intensity", cfg.intensity);
    cfg.lengthScale = td.getScalarOr("lengthScale", cfg.lengthScale);
    cfg.relaxation = td.getScalarOr("relaxation", cfg.relaxation);
    cfg.validate();
    return cfg;
}

std::unique_ptr<RansModel> RansModel::create(const polyMesh& mesh,
                                            volVectorField& U,
                                            const surfaceScalarField& phi,
                                            scalar nu,
                                            const Schemes& schemes,
                                            const TurbulenceConfig& cfg)
{
    cfg.validate();
    if (cfg.model == "laminar")
    {
        auto m = std::make_unique<LaminarModel>(mesh, U, phi, nu, schemes, cfg);
        return m;
    }
    if (cfg.model == "kEpsilon")
    {
        auto m = std::make_unique<KEpsilon>(mesh, U, phi, nu, schemes, cfg);
        m->init();
        return m;
    }
    if (cfg.model == "kOmegaSST")
    {
        auto m = std::make_unique<KOmegaSST>(mesh, U, phi, nu, schemes, cfg);
        m->init();
        return m;
    }
    cfFatal("turbulence: unknown model '" + cfg.model
            + "' (laminar | kEpsilon | kOmegaSST)");
}

// ------------------------------------------------------------------
// base
// ------------------------------------------------------------------

RansModel::RansModel(const polyMesh& mesh, volVectorField& U, const surfaceScalarField& phi,
                     scalar nu, const Schemes& schemes, const TurbulenceConfig& cfg,
                     std::string dissName)
    : mesh_(mesh), U_(U), phi_(phi), nu_(nu), schemes_(schemes), cfg_(cfg),
      k_(mesh, "k", 0.0), diss_(mesh, dissName, 0.0),
      nut_(static_cast<std::size_t>(mesh.nCells()), 0.0),
      nutFace_(static_cast<std::size_t>(mesh.nFaces()), 0.0),
      nuEffFace_(static_cast<std::size_t>(mesh.nFaces()), nu),
      strain_(static_cast<std::size_t>(mesh.nCells()), 0.0),
      prod_(static_cast<std::size_t>(mesh.nCells()), 0.0),
      eps_(static_cast<std::size_t>(mesh.nCells()), 0.0)
{
    if (!(nu_ > 0)) cfFatal("turbulence: molecular viscosity must be positive");

    for (const polyMesh::Patch& p : mesh_.patches())
    {
        const label patch = mesh_.patchIdByName(p.name);
        k_.setPatchValue(patch, isWallPatch(p) ? "zeroGradient"
                                               : (isInletPatch(p) ? "fixedValue" : "zeroGradient"), 0.0);
        diss_.setPatchValue(patch, isWallPatch(p) || isInletPatch(p) ? "fixedValue" : "zeroGradient", 0.0);
    }

    for (const polyMesh::Patch& p : mesh_.patches())
    {
        if (!isWallPatch(p)) continue;
        for (label f = p.start; f < p.start + p.size; ++f) wallFaces_.push_back(f);
    }

    const std::size_t nw = wallFaces_.size();
    uTau_.assign(nw, 0.0);
    yPlus_.assign(nw, 0.0);
    nutWall_.assign(nw, 0.0);
    wallDist_.assign(nw, 0.0);
    for (std::size_t i = 0; i < nw; ++i)
    {
        const std::size_t f = static_cast<std::size_t>(wallFaces_[i]);
        const Vector& sf = mesh_.faceSf()[f];
        const scalar a = mag(sf);
        if (a <= 0) cfFatal("turbulence: zero-area wall face");
        // normal distance from the owner cell centre to the wall face
        wallDist_[i] = dot(mesh_.deltaCoeffVec()[f], sf / a);
        if (!(wallDist_[i] > 0)) cfFatal("turbulence: wall face not outside its owner cell");
    }
}

scalar RansModel::uReference() const
{
    scalar u = 0;
    for (const polyMesh::Patch& p : mesh_.patches())
    {
        if (!isInletPatch(p)) continue;
        const label patch = mesh_.patchIdByName(p.name);
        if (U_.patchBcType(patch) != "fixedValue") continue;
        u = std::max(u, mag(U_.patchValue(patch)));
    }
    if (u > 0) return u;
    for (const Vector& v : U_.internalField()) u = std::max(u, mag(v));
    return u;
}

scalar RansModel::lengthScale() const
{
    if (cfg_.lengthScale > 0) return cfg_.lengthScale;
    const Vector d = mesh_.maxX() - mesh_.minX();
    return 0.07 * std::sqrt(dot(d, d));
}

void RansModel::initializeUniform(scalar k0, scalar diss0)
{
    k0 = std::max(k0, turb::kFloor);
    diss0 = std::max(diss0, turb::dissFloor);

    for (const polyMesh::Patch& p : mesh_.patches())
    {
        if (!isInletPatch(p)) continue;
        const label patch = mesh_.patchIdByName(p.name);
        k_.setPatchValue(patch, "fixedValue", k0);
        diss_.setPatchValue(patch, "fixedValue", diss0);
    }

    for (label c = 0; c < mesh_.nCells(); ++c)
    {
        k_.internalField()[static_cast<std::size_t>(c)] = k0;
        diss_.internalField()[static_cast<std::size_t>(c)] = diss0;
    }

    updateNut();
    updateTurbDissipation();
    refreshFaceNut();
}

Field<scalar> RansModel::diffusionFace(const Field<scalar>& gammaCell,
                                      const Field<scalar>& sigmaCell) const
{
    Field<scalar> g = harmonicFaceGamma(mesh_, gammaCell);
    for (std::size_t i = 0; i < wallFaces_.size(); ++i)
    {
        const label f = wallFaces_[i];
        const label o = mesh_.owner()[static_cast<std::size_t>(f)];
        g[static_cast<std::size_t>(f)] = nu_ + nutWall_[i] * sigmaCell[static_cast<std::size_t>(o)];
    }
    return g;
}

void RansModel::computeStrainRate()
{
    const FieldView view = viewOf(U_);
    Field<scalar> g[3];
    for (int comp = 0; comp < 3; ++comp)
    {
        g[comp] = component(U_.internalField(), comp);
    }

    // gradU[c] = gradient of the c-th velocity component, so
    // gradU[c][s] = d U_c / d x_s
    Field<Vector> gradU[3];
    for (int comp = 0; comp < 3; ++comp)
    {
        gradU[comp] = fvc::grad(mesh_, g[comp], makeFaceValueFn(mesh_, view, comp), schemes_.grad);
    }

    const std::size_t nC = static_cast<std::size_t>(mesh_.nCells());
    for (std::size_t c = 0; c < nC; ++c)
    {
        scalar s2 = 0;
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                const scalar sij = 0.5 * (gradU[i][c][j] + gradU[j][c][i]);
                s2 += sij * sij;
            }
        }
        // |S| = sqrt(2 S:S)
        strain_[c] = std::sqrt(2.0 * s2);
        prod_[c] = nut_[c] * 2.0 * s2;
    }
}

void RansModel::computeWallLayer()
{
    const auto& sf = mesh_.faceSf();
    const auto& own = mesh_.owner();

    for (std::size_t i = 0; i < wallFaces_.size(); ++i)
    {
        const label f = wallFaces_[i];
        const std::size_t u = static_cast<std::size_t>(f);
        const Vector n = sf[u] / mag(sf[u]);
        const Vector& up = U_.internalField()[static_cast<std::size_t>(own[u])];
        const Vector un = n * dot(up, n);
        const scalar parallel = mag(up - un);

        const WallState ws = wallState(parallel, wallDist_[i], nu_);
        uTau_[i] = ws.uTau;
        yPlus_[i] = ws.yPlus;
        // The laminar model reports the y+ its mesh produces but must leave the
        // momentum equation at the molecular viscosity, so its wall nut stays 0.
        nutWall_[i] = hasDissipation() ? ws.nut : 0.0;

        if (hasDissipation())
        {
            const scalar kCell = k_.internalField()[static_cast<std::size_t>(own[u])];
            diss_.setFaceValue(f, wallDissipation(ws, kCell, wallDist_[i]));
        }
    }
}

void RansModel::solveK(scalar dt)
{
    const Field<scalar> kOld = k_.internalField();
    const auto& vol = mesh_.cellVolumes();

    fvMatrix eqn(mesh_, viewOf(k_), "kEqn");
    eqn.schemes() = schemes_;
    eqn.schemes().relaxation = 1.0;
    if (dt > 0) fvm::ddt(eqn, kOld, dt);
    fvm::div(eqn, phi_.values(), schemes_.div);

    const Field<scalar> sigma = kDiffusivity();
    Field<scalar> gamma(static_cast<std::size_t>(mesh_.nCells()));
    for (std::size_t c = 0; c < gamma.size(); ++c)
    {
        gamma[c] = nu_ + nut_[c] * sigma[c];
    }
    fvm::laplacianFace(eqn, diffusionFace(gamma, sigma), schemes_.laplacian, -1);

    fvm::sp(eqn, destructionSp(), 1);
    const Field<scalar> pk = kProduction();
    for (label c = 0; c < mesh_.nCells(); ++c)
    {
        const std::size_t u = static_cast<std::size_t>(c);
        eqn.addRhs(c, 0, vol[u] * pk[u]);
    }

    eqn.solve(fvMatrix::Controls{turb::linearTol, turb::linearMaxIter, 2});

    const Field<scalar>& x = eqn.solution()[0];
    const scalar rl = cfg_.relaxation;
    for (std::size_t c = 0; c < k_.internalField().size(); ++c)
    {
        const scalar kNew = kOld[c] + rl * (x[c] - kOld[c]);
        k_.internalField()[c] = std::max(kNew, turb::kFloor);
    }
}

void RansModel::refreshFaceNut()
{
    nutFace_ = harmonicFaceGamma(mesh_, nut_);
    for (std::size_t i = 0; i < wallFaces_.size(); ++i)
    {
        nutFace_[static_cast<std::size_t>(wallFaces_[i])] = nutWall_[i];
    }
    for (std::size_t f = 0; f < nuEffFace_.size(); ++f)
    {
        nuEffFace_[f] = nu_ + nutFace_[f];
    }
}

void RansModel::imposeWallEquilibrium()
{
    if (!hasDissipation()) return;

    const std::size_t nCells = static_cast<std::size_t>(mesh_.nCells());
    Field<scalar> kEq(nCells, 0.0), dissEq(nCells, 0.0);
    for (std::size_t i = 0; i < wallFaces_.size(); ++i)
    {
        if (!(uTau_[i] > 0)) continue;
        scalar k = 0, d = 0;
        const WallState ws{yPlus_[i], uTau_[i], nutWall_[i]};
        wallEquilibrium(ws, wallDist_[i], k, d);
        if (!(k > 0) || !(d > 0)) continue;

        const std::size_t o =
            static_cast<std::size_t>(mesh_.owner()[static_cast<std::size_t>(wallFaces_[i])]);
        // the wall with the strongest shear governs a cell that touches two
        // layers, and it sets both values so that the pair stays consistent and
        // the eddy viscosity it implies is nut = kappa u_tau y
        if (k > kEq[o])
        {
            kEq[o] = k;
            dissEq[o] = d;
        }
    }
    for (std::size_t c = 0; c < nCells; ++c)
    {
        if (!(kEq[c] > 0)) continue;
        k_.internalField()[c] = kEq[c];
        diss_.internalField()[c] = dissEq[c];
    }
}

void RansModel::correctTurbulence(scalar dt)
{
    computeStrainRate();
    computeWallLayer();
    if (hasDissipation())
    {
        solveK(dt);
        solveDissipation(dt);
        imposeWallEquilibrium();
    }
    updateNut();
    updateTurbDissipation();
    refreshFaceNut();
}

scalar RansModel::maxWallYPlus() const
{
    scalar m = 0;
    for (const scalar v : yPlus_) m = std::max(m, v);
    return m;
}

scalar RansModel::meanWallYPlus() const
{
    scalar s = 0, n = 0;
    for (const scalar v : yPlus_)
    {
        if (v <= 0) continue;
        s += v;
        ++n;
    }
    return n > 0 ? s / n : 0;
}

scalar RansModel::meanFrictionVelocity() const
{
    const auto& area = mesh_.faceAreas();
    scalar s = 0, a = 0;
    for (std::size_t i = 0; i < wallFaces_.size(); ++i)
    {
        const scalar af = area[static_cast<std::size_t>(wallFaces_[i])];
        s += uTau_[i] * af;
        a += af;
    }
    return a > 0 ? s / a : 0;
}

const Field<scalar>& RansModel::cellWallDistance() const
{
    if (cellDistCalculated_) return cellWallDist_;
    cellDistCalculated_ = true;

    const label nC = mesh_.nCells();
    cellWallDist_.assign(static_cast<std::size_t>(nC), -1.0);

    using Node = std::pair<scalar, label>;
    std::priority_queue<Node, std::vector<Node>, std::greater<Node>> pq;

    for (std::size_t i = 0; i < wallFaces_.size(); ++i)
    {
        const label o = mesh_.owner()[static_cast<std::size_t>(wallFaces_[i])];
        const std::size_t u = static_cast<std::size_t>(o);
        if (cellWallDist_[u] < 0 || wallDist_[i] < cellWallDist_[u])
        {
            cellWallDist_[u] = wallDist_[i];
        }
    }
    for (label c = 0; c < nC; ++c)
    {
        const std::size_t u = static_cast<std::size_t>(c);
        if (cellWallDist_[u] >= 0) pq.push({cellWallDist_[u], c});
    }

    // Dijkstra over the cell graph: crossing a face costs the owner-neighbour
    // centre distance, so the shortest path is the wall distance of an
    // orthogonal column measured through cell centres.
    while (!pq.empty())
    {
        const auto [d, c] = pq.top();
        pq.pop();
        const std::size_t u = static_cast<std::size_t>(c);
        if (d > cellWallDist_[u] + 1e-12 * (d + 1.0)) continue;
        for (const label f : mesh_.cells()[u])
        {
            if (f >= mesh_.nInternalFaces()) continue;
            const std::size_t fu = static_cast<std::size_t>(f);
            // faces are stored with owner < neighbour, so a cell is reached from
            // either side of the face
            const label o = mesh_.owner()[fu];
            const label n = (o == c) ? mesh_.neighbour()[fu] : o;
            const scalar nd = d + mag(mesh_.deltaCoeffVec()[fu]);
            const std::size_t nu = static_cast<std::size_t>(n);
            if (cellWallDist_[nu] < 0 || nd < cellWallDist_[nu] - 1e-12 * (nd + 1.0))
            {
                cellWallDist_[nu] = nd;
                pq.push({nd, n});
            }
        }
    }

    // A domain without walls has no near-wall region at all: the blend
    // functions then read a distance of the order of the whole domain, which is
    // the free-stream limit of F1. A cell left unreached in a mesh that does
    // have walls is a disconnected component, and the same far-field value is
    // the safest reading for it.
    scalar reached = 0;
    for (const scalar v : cellWallDist_) reached = std::max(reached, v);
    if (reached <= 0)
    {
        const Vector d = mesh_.maxX() - mesh_.minX();
        reached = std::sqrt(dot(d, d));
    }
    for (scalar& v : cellWallDist_)
    {
        if (v < 0) v = reached;
    }
    return cellWallDist_;
}

// ------------------------------------------------------------------
// laminar
// ------------------------------------------------------------------

LaminarModel::LaminarModel(const polyMesh& mesh, volVectorField& U, const surfaceScalarField& phi,
                           scalar nu, const Schemes& schemes, const TurbulenceConfig& cfg)
    : RansModel(mesh, U, phi, nu, schemes, cfg, "none")
{
    refreshFaceNut();
}

void LaminarModel::initializeFields() {}

void LaminarModel::updateNut()
{
    nut_.assign(static_cast<std::size_t>(mesh_.nCells()), 0.0);
}

Field<scalar> LaminarModel::kDiffusivity() const
{
    return Field<scalar>(static_cast<std::size_t>(mesh_.nCells()), 0.0);
}

Field<scalar> LaminarModel::destructionSp() const
{
    return Field<scalar>(static_cast<std::size_t>(mesh_.nCells()), 0.0);
}

scalar LaminarModel::wallDissipation(const WallState&, scalar, scalar) const
{
    return 0;
}

} // namespace cf
