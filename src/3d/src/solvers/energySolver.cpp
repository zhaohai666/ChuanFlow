#include "cf/fvm/EnergySolver.h"

#include <algorithm>
#include <cmath>

#include "cf/common/Error.h"
#include "cf/common/Log.h"
#include "cf/turbulence/WallFunction.h"

namespace cf
{

namespace
{

int axisIndex(const std::string& axis)
{
    if (axis == "x") return 0;
    if (axis == "y") return 1;
    if (axis == "z") return 2;
    cfFatal("energy: initialAxis must be x, y or z");
}

bool finiteAll(const EnergyProps& p)
{
    if (!std::isfinite(p.rho) || !std::isfinite(p.cp) || !std::isfinite(p.kappa)) return false;
    if (!std::isfinite(p.turbulentPrandtl) || !std::isfinite(p.volumetricSource)) return false;
    if (!std::isfinite(p.beta) || !std::isfinite(p.tRef) || !std::isfinite(p.relaxation)) return false;
    if (!std::isfinite(p.initialT) || !std::isfinite(p.initialAmplitude)) return false;
    return std::isfinite(p.gravity.x) && std::isfinite(p.gravity.y) && std::isfinite(p.gravity.z);
}

} // namespace

EnergyProps EnergyProps::fromCase(const Dictionary& caseDict)
{
    if (!caseDict.has("energy")) cfFatal("energy: case has no energy{} block");
    const Dictionary& e = caseDict.getDict("energy");

    EnergyProps p;
    p.rho = e.getScalarOr("rho", p.rho);
    p.cp = e.getScalarOr("cp", p.cp);
    p.kappa = e.getScalarOr("kappa", p.kappa);
    p.turbulentPrandtl = e.getScalarOr("turbulentPrandtl", p.turbulentPrandtl);
    p.volumetricSource = e.getScalarOr("volumetricSource", p.volumetricSource);
    p.viscousDissipation = e.has("viscousDissipation") ? e.getBool("viscousDissipation") : false;
    p.beta = e.getScalarOr("expansion", p.beta);
    p.tRef = e.getScalarOr("referenceT", p.tRef);
    if (e.has("gravity"))
    {
        const auto g = e.getScalarList("gravity");
        p.gravity = Vector{g[0], g[1], g[2]};
    }
    p.relaxation = e.getScalarOr("relaxation", p.relaxation);
    p.initialT = e.getScalarOr("initialT", p.initialT);
    p.initialMode = e.has("initialMode") ? e.getInt("initialMode") : 0;
    p.initialAmplitude = e.getScalarOr("initialAmplitude", p.initialAmplitude);
    p.initialAxis = e.getWordOr("initialAxis", p.initialAxis);
    p.validate();
    return p;
}

void EnergyProps::validate() const
{
    if (!finiteAll(*this)) cfFatal("energy: non-finite property");
    if (!(rho > 0) || !(cp > 0) || !(kappa > 0))
        cfFatal("energy: rho, cp and kappa must be positive");
    if (!(turbulentPrandtl > 0))
        cfFatal("energy: turbulentPrandtl must be positive");
    if (!(relaxation > 0 && relaxation <= 1))
        cfFatal("energy: relaxation must be in (0, 1]");
    if (initialMode < 0) cfFatal("energy: negative initialMode");
}

EnergyModel::EnergyModel(const polyMesh& mesh, const surfaceScalarField& phi,
                         const Schemes& schemes, const EnergyProps& props, scalar nu,
                         const Field<scalar>* nut, const Field<scalar>* strainRate,
                         const EnergyWallLayer* wallLayer)
    : mesh_(mesh), phi_(phi), schemes_(schemes), props_(props), nu_(nu), nut_(nut),
      strainRate_(strainRate), t_(mesh, "T", props.initialT),
      buoyancy_(static_cast<std::size_t>(mesh.nCells()), Vector::zero())
{
    if (static_cast<label>(phi_.values().size()) != mesh_.nFaces())
        cfFatal("energy: flux size does not match the mesh");
    if (nut_ && static_cast<label>(nut_->size()) != mesh_.nCells())
        cfFatal("energy: eddy viscosity size does not match the mesh");
    if (strainRate_ && static_cast<label>(strainRate_->size()) != mesh_.nCells())
        cfFatal("energy: strain rate size does not match the mesh");
    if (props_.viscousDissipation && !strainRate_)
        cfFatal("energy: viscousDissipation needs the strain rate the turbulence "
                "layer evaluates; add turbulence { model laminar; }");
    if (wallLayer) wallLayer_ = *wallLayer;

    hasBuoyancy_ = props_.beta != 0 && mag2(props_.gravity) > 0;
    refreshBuoyancy();
}

void EnergyModel::initialize()
{
    const auto& cc = mesh_.cellCentres();
    const Vector lo = mesh_.minX();
    const Vector hi = mesh_.maxX();
    const int axis = axisIndex(props_.initialAxis);
    const scalar span = hi[axis] - lo[axis];
    if (props_.initialMode > 0 && !(span > 0))
        cfFatal("energy: initialAxis has no extent");

    auto& t = t_.internalField();
    for (label c = 0; c < mesh_.nCells(); ++c)
    {
        const std::size_t u = static_cast<std::size_t>(c);
        scalar v = props_.initialT;
        if (props_.initialMode > 0)
        {
            const scalar xi = (cc[u][axis] - lo[axis]) / span;
            v += props_.initialAmplitude * std::cos(props_.initialMode * M_PI * xi);
        }
        t[u] = v;
    }
    refreshBuoyancy();
}

Field<scalar> EnergyModel::alphaCell() const
{
    Field<scalar> a(static_cast<std::size_t>(mesh_.nCells()), props_.alphaMolecular());
    if (nut_)
    {
        const scalar invPr = 1.0 / props_.turbulentPrandtl;
        for (std::size_t u = 0; u < a.size(); ++u) a[u] += (*nut_)[u] * invPr;
    }
    return a;
}

Field<scalar> EnergyModel::sourceRate() const
{
    const std::size_t n = static_cast<std::size_t>(mesh_.nCells());
    Field<scalar> s(n, props_.volumetricSource / (props_.rho * props_.cp));
    if (props_.viscousDissipation)
    {
        // molecular stress only: the turbulent work goes into k and its heat
        // release would come from the k transport, which this increment does not
        // read. cp converts m2/s3 to K/s.
        for (std::size_t u = 0; u < n; ++u)
        {
            const scalar g = (*strainRate_)[u];
            s[u] += nu_ * g * g / props_.cp;
        }
    }
    return s;
}

Field<scalar> EnergyModel::faceGamma() const
{
    Field<scalar> g = harmonicFaceGamma(mesh_, alphaCell());
    if (!wallLayer_.faces) return g;

    // the wall function only owns faces whose patch prescribes the wall
    // temperature; a heatFlux or adiabatic wall keeps its own flux definition
    bool anyFixedValue = false;
    for (label f : *wallLayer_.faces)
    {
        if (t_.patchBcType(t_.patchIdOfFace(f)) == "fixedValue")
        {
            anyFixedValue = true;
            break;
        }
    }
    if (!anyFixedValue) return g;

    const scalar Pr = nu_ / props_.alphaMolecular();
    if (Pr < 0.5)
        cfFatal("energy: the thermal wall function needs a molecular Prandtl of "
                "at least 0.5 (got " + std::to_string(Pr) +
                "); a liquid-metal temperature law is a separate increment");
    const scalar yT = thermalSublayerYPlus(Pr, props_.turbulentPrandtl);

    for (std::size_t i = 0; i < wallLayer_.faces->size(); ++i)
    {
        const label f = (*wallLayer_.faces)[i];
        if (t_.patchBcType(t_.patchIdOfFace(f)) != "fixedValue") continue;
        const scalar tPlus = thermalTPlus((*wallLayer_.yPlus)[i], Pr,
                                          props_.turbulentPrandtl, yT);
        g[static_cast<std::size_t>(f)] = (*wallLayer_.uTau)[i]
                                       * (*wallLayer_.wallDist)[i] / tPlus;
    }
    return g;
}

void EnergyModel::correctEnergy(scalar dt)
{
    const label nC = mesh_.nCells();
    const auto& vol = mesh_.cellVolumes();
    const Field<scalar> oldT = t_.internalField();
    const Field<scalar> src = sourceRate();

    fvMatrix eqn(mesh_, viewOf(t_), "TEqn");
    eqn.schemes() = schemes_;
    // the matrix' own relaxation would blend against its zero initial guess;
    // the steady blend below is against the field the equation started from
    eqn.schemes().relaxation = 1.0;
    if (dt > 0) fvm::ddt(eqn, oldT, dt);
    fvm::div(eqn, phi_.values(), schemes_.div);
    fvm::laplacianFace(eqn, faceGamma(), schemes_.laplacian, -1);
    for (label c = 0; c < nC; ++c)
    {
        eqn.addRhs(c, 0, vol[static_cast<std::size_t>(c)] * src[static_cast<std::size_t>(c)]);
    }

    const SolverReport rep = eqn.solve(fvMatrix::Controls{linearTol, linearMaxIter, 2});
    if (!rep.converged)
    {
        log::warn("energy: temperature solve did not converge (", rep.iterations,
                  " iters, residual ", rep.finalResidual, ")");
    }

    const Field<scalar>& x = eqn.solution()[0];
    auto& t = t_.internalField();
    if (dt > 0)
    {
        // implicit Euler: the step is the solution, no relaxation
        t = x;
    }
    else
    {
        const scalar rl = props_.relaxation;
        for (std::size_t u = 0; u < t.size(); ++u) t[u] = oldT[u] + rl * (x[u] - oldT[u]);
    }
    refreshBuoyancy();
}

const Field<Vector>* EnergyModel::buoyancyForce() const
{
    return hasBuoyancy_ ? &buoyancy_ : nullptr;
}

void EnergyModel::refreshBuoyancy()
{
    if (!hasBuoyancy_) return;
    const auto& t = t_.internalField();
    for (std::size_t u = 0; u < buoyancy_.size(); ++u)
    {
        buoyancy_[u] = -props_.beta * (t[u] - props_.tRef) * props_.gravity;
    }
}

scalar EnergyModel::patchHeatRate(label patch) const
{
    const std::size_t np = static_cast<std::size_t>(patch);
    if (patch < 0 || patch >= mesh_.nPatches())
        cfFatal("energy: no such patch for the heat rate");
    const Field<scalar> gamma = faceGamma();
    const auto& sf = mesh_.faceSf();
    const auto& d = mesh_.deltaCoeffVec();
    const auto& patches = mesh_.patches();
    const auto& t = t_.internalField();

    scalar rate = 0;
    for (label f = patches[np].start; f < patches[np].start + patches[np].size; ++f)
    {
        const std::size_t u = static_cast<std::size_t>(f);
        const scalar d2 = dot(d[u], d[u]);
        if (d2 <= 0) continue;
        const scalar coef = gamma[u] * dot(sf[u], d[u]) / d2;
        scalar c = 0, k = 0;
        t_.linearise(f, c, k);
        const scalar tp = t[static_cast<std::size_t>(mesh_.owner()[u])];
        // heat entering the fluid = rho cp * the diffusive flux the equation
        // carries, so a fixedGradient face reports its prescribed flux exactly
        rate += coef * (c * tp + k - tp);
    }
    return rate * props_.rho * props_.cp;
}

scalar EnergyModel::volumetricHeatRate() const
{
    return props_.volumetricSource * sum(mesh_.cellVolumes());
}

scalar EnergyModel::netEnthalpyRate() const
{
    const auto& ph = phi_.values();
    scalar out = 0;
    for (label f = mesh_.nInternalFaces(); f < mesh_.nFaces(); ++f)
    {
        out += ph[static_cast<std::size_t>(f)] * t_.faceValue(f);
    }
    return out * props_.rho * props_.cp;
}

EnergyModel::ThermalWallReport EnergyModel::wallLayerReport(label patch) const
{
    ThermalWallReport rep;
    if (patch < 0 || patch >= mesh_.nPatches() || !wallLayer_.faces) return rep;

    const auto& sf = mesh_.faceSf();
    const auto& patches = mesh_.patches();
    const auto& t = t_.internalField();
    const scalar Pr = nu_ / props_.alphaMolecular();
    const scalar yT = thermalSublayerYPlus(Pr, props_.turbulentPrandtl);

    scalar area = 0;
    for (std::size_t i = 0; i < wallLayer_.faces->size(); ++i)
    {
        const label f = (*wallLayer_.faces)[i];
        const std::size_t u = static_cast<std::size_t>(f);
        if (f < patches[static_cast<std::size_t>(patch)].start
            || f >= patches[static_cast<std::size_t>(patch)].start
                       + patches[static_cast<std::size_t>(patch)].size)
            continue;
        if (t_.patchBcType(patch) != "fixedValue") continue;
        const scalar a = mag(sf[u]);
        const scalar tPlus = thermalTPlus((*wallLayer_.yPlus)[i], Pr,
                                          props_.turbulentPrandtl, yT);
        rep.uTau += a * (*wallLayer_.uTau)[i];
        rep.yPlus += a * (*wallLayer_.yPlus)[i];
        rep.tPlus += a * tPlus;
        const scalar tw = t_.patchValue(patch);
        const scalar tp = t[static_cast<std::size_t>(mesh_.owner()[u])];
        // q = rho cp u_tau (T_w - T_p)/T+ into the fluid, per the same law the
        // face diffusivity came from
        rep.heatRate += props_.rho * props_.cp * (*wallLayer_.uTau)[i] / tPlus
                      * (tw - tp) * a;
        area += a;
    }
    if (area > 0)
    {
        rep.uTau /= area;
        rep.yPlus /= area;
        rep.tPlus /= area;
        rep.active = true;
    }
    return rep;
}

} // namespace cf
