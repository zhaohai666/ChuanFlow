#include "cf/1d/TwoFluid.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "cf/common/Error.h"

#include "cf/1d/Network.h" // churchillFriction
#include "DenseSolve.h"

namespace cf
{
namespace
{

constexpr double kGravity = 9.80665;

double mixtureMu(const TwoFluidFluid& f, double a)
{
    return a * f.muL + (1.0 - a) * f.muG;
}

/// Godunov flux of the quadratic F(a) = qa a^2 + qb a + qc between the left
/// and right states. Convex F (qa > 0) shocks on a falling jump and
/// rarefies on a rising one, concave F mirrors both, and the linear case
/// (qa = 0) is plain upwind on the constant slope. The rarefaction branch
/// evaluates F at the sonic state -qb/(2 qa) clamped into the jump.
double godunovQuadratic(double qa, double qb, double qc, double aL, double aR)
{
    const double fL = qc + aL * (qb + qa * aL);
    if (aL == aR) return fL;
    const double fR = qc + aR * (qb + qa * aR);
    const bool shock = qa > 0 ? aL > aR : aL < aR;
    if (shock)
    {
        const double s = (fL - fR) / (aL - aR);
        return s >= 0 ? fL : fR;
    }
    if (qa == 0) return qb >= 0 ? fL : fR;
    const double sonic = -qb / (2.0 * qa);
    const double lo = std::min(aL, aR), hi = std::max(aL, aR);
    const double a = std::min(std::max(sonic, lo), hi);
    return qc + a * (qb + qa * a);
}

} // namespace

TwoFluidPipe::TwoFluidPipe(const TwoFluidCase& c)
    : cfg_(c)
{
    const TwoFluidFluid& f = cfg_.fluid;
    if (!(f.rhoL > 0) || !(f.rhoG > 0) || !(f.muL > 0) || !(f.muG > 0))
        cfFatal("twoFluid: fluid needs positive densities and viscosities");
    if (f.gasR < 0)
        cfFatal("twoFluid: negative gasR");
    if (cfg_.pipes.empty())
        cfFatal("twoFluid: case has no pipe segment");
    diameter_ = cfg_.pipes.front().diameter;
    length_ = 0;
    for (const TwoFluidGeo& g : cfg_.pipes)
    {
        if (!(g.length > 0) || !(g.diameter > 0) || g.sections < 1)
            cfFatal("twoFluid: pipe needs positive length/diameter and >= 1"
                    " sections");
        if (g.roughness < 0)
            cfFatal("twoFluid: negative roughness");
        // the mixture momentum is lumped on one velocity, which needs a
        // constant bore along the chain
        if (g.diameter != diameter_)
            cfFatal("twoFluid: a serial chain needs a uniform diameter"
                    " (the lumped mixture momentum carries one velocity)");
        length_ += g.length;
    }

    const TwoFluidBoundary& in = cfg_.inlet;
    const TwoFluidBoundary& out = cfg_.outlet;
    if (!cfg_.networkMode)
    {
        // serial-case boundary combinations; network members are checked
        // by the network loader instead (the junction end boundary is
        // owned by the TwoFluidNetwork and rewritten every step)
        const bool inSealed = in.type == TwoFluidBoundary::Type::sealed;
        const bool outSealed = out.type == TwoFluidBoundary::Type::sealed;
        if (inSealed && !outSealed)
            cfFatal("twoFluid: a sealed inlet needs a sealed outlet"
                    " (incompressible phases cannot draw through-flow)");
        if (outSealed && in.type == TwoFluidBoundary::Type::fixedPressure)
            cfFatal("twoFluid: fixed-pressure inlet into a sealed pipe has no"
                    " steady state with incompressible phases");
        if (in.type == TwoFluidBoundary::Type::phaseRates
            && in.vsl + in.vsg < 0)
            cfFatal("twoFluid: phaseRates injection with negative total rate"
                    " is not supported");
    }
    if (cfg_.closure.model == 2 && !(cfg_.closure.interfacialDrag > 0))
        cfFatal("twoFluid: the dynamic closure needs a positive"
                " interfacialDrag");

    n_ = 0;
    for (const TwoFluidGeo& g : cfg_.pipes) n_ += g.sections;
    dxs_.reserve(n_);
    sinT_.reserve(n_);
    relrough_.reserve(n_);
    for (const TwoFluidGeo& g : cfg_.pipes)
    {
        const double dx = g.length / g.sections;
        const double st = std::sin(g.angleDeg * M_PI / 180.0);
        for (int k = 0; k < g.sections; ++k)
        {
            dxs_.push_back(dx);
            sinT_.push_back(st);
            relrough_.push_back(g.roughness / g.diameter);
        }
    }
    alpha_.assign(n_, cfg_.alpha0);
    water_.assign(n_, cfg_.water0);
    slip_.assign(n_, 0.0);
    sDotOld_.assign(n_, 0.0);
    std::vector<double> centers(n_);
    double arc = 0;
    for (int i = 0; i < n_; ++i)
    {
        centers[i] = arc + 0.5 * dxs_[i];
        arc += dxs_[i];
    }
    for (const TwoFluidFill& fill : cfg_.fills)
    {
        if (fill.from > fill.to || fill.alpha < 0 || fill.alpha > 1)
            cfFatal("twoFluid: fill block needs from <= to and alpha in"
                    " [0,1]");
        for (int i = 0; i < n_; ++i)
        {
            if (centers[i] >= fill.from && centers[i] <= fill.to)
                alpha_[i] = fill.alpha;
        }
    }
    // incompressible phases: an imposed rate is felt along the whole pipe
    // from t = 0; pressure-driven and sealed pipes start from rest
    if (cfg_.inlet.type == TwoFluidBoundary::Type::phaseRates)
        vM_ = vMOld_ = cfg_.inlet.vsl + cfg_.inlet.vsg;

    // vertical drop of each cell centre from the chain inlet (geothermal
    // depth), accumulated along the per-segment slopes
    centerDrop_.assign(n_, 0.0);
    double drop = 0;
    for (int i = 0; i < n_; ++i)
    {
        drop += 0.5 * dxs_[i] * sinT_[i];
        centerDrop_[i] = drop;
        drop += 0.5 * dxs_[i] * sinT_[i];
    }

    // gas inventory per unit area: rhoG (1 - alpha) dx; with the EOS on it
    // is the transported state, otherwise it stays at its initial value
    // and every density lookup falls back to the constant fluid rhoG
    gasMass_.resize(n_);
    for (int i = 0; i < n_; ++i)
        gasMass_[i] = f.rhoG * (1.0 - alpha_[i]) * dxs_[i];

    if (cfg_.thermalOn)
    {
        const TwoFluidThermal& th = cfg_.thermal;
        if (!(f.cL > 0) || !(f.cG > 0))
            cfFatal("twoFluid: thermal on needs positive fluid cL/cG");
        if (th.wallModel != 0 && th.wallModel != 1)
            cfFatal("twoFluid: thermal wallModel is constant or geothermal");
        if (th.wallModel == 1)
        {
            if (!(th.geothermalGradient > 0))
                cfFatal("twoFluid: geothermal wall needs a positive gradient");
            double dropTotal = 0;
            for (int i = 0; i < n_; ++i) dropTotal += dxs_[i] * sinT_[i];
            const double dExit = th.inletDepth - dropTotal;
            if (dExit < 0)
                cfFatal("twoFluid: geothermal pipe leaves the ground "
                        "(depth at the outlet is " + std::to_string(dExit)
                        + " m)");
        }
        // both ghost temperatures are required: the outlet one only feeds
        // backflow, but backflow is a normal fixed-pressure outcome
        if (cfg_.inlet.type != TwoFluidBoundary::Type::sealed
            && !(cfg_.inlet.temperature > 0))
            cfFatal("twoFluid: thermal on needs an inlet boundary temperature");
        if (cfg_.outlet.type != TwoFluidBoundary::Type::sealed
            && !(cfg_.outlet.temperature > 0))
            cfFatal("twoFluid: thermal on needs an outlet boundary temperature");
        T_.assign(n_, cfg_.temperature0);
    }

    if (cfg_.hydrateOn)
    {
        if (!eosOn())
            cfFatal("twoFluid: hydrate kinetics needs the gas EOS"
                    " (fluid gasR > 0)");
        if (!cfg_.thermalOn)
            cfFatal("twoFluid: hydrate kinetics needs a thermal{} block");
        const TwoFluidHydrate& hy = cfg_.hydrate;
        if (!(hy.growthRate >= 0) || !(hy.dissociationRate >= 0))
            cfFatal("twoFluid: negative hydrate rate");
        if (!(hy.gasContent > 0))
            cfFatal("twoFluid: hydrate needs a positive gasContent");
        if (!(hy.latentHeat > 0))
            cfFatal("twoFluid: hydrate needs a positive latentHeat");
        if (hy.initial < 0 || hy.initial > 1)
            cfFatal("twoFluid: initial hydrate fraction outside [0,1]");
        if (cfg_.water0 + hy.initial > 1.0)
            cfFatal("twoFluid: water cut plus initial hydrate fraction"
                    " exceeds 1");
        hydrate_.assign(n_, hy.initial);
        hydrateDelta_.assign(n_, 0.0);
    }
}

void TwoFluidPipe::fluxCoeffs(double& qa, double& qb, double& qc) const
{
    // F(a) = a_l * vL(a): liquid volume flux per unit area (algebraic
    // closures; the dynamic closure builds face-local coefficients from
    // the slip state)
    if (cfg_.closure.model == 0) // wallis: vG = vM + alpha*us
    {
        qa = cfg_.closure.slipVelocity;
        qb = vM_ - cfg_.closure.slipVelocity;
        qc = 0.0;
    }
    else // zuber: vG = C0 vM + vD
    {
        qa = 0.0;
        qb = cfg_.closure.c0 * vM_ + cfg_.closure.driftVelocity;
        qc = vM_ * (1.0 - cfg_.closure.c0) - cfg_.closure.driftVelocity;
    }
}

double TwoFluidPipe::liquidFlux(double alpha, double s) const
{
    if (cfg_.closure.model == 2) // dynamic: vL = vM - alpha_g s
        return s * alpha * alpha + (vM_ - s) * alpha;
    double qa, qb, qc;
    fluxCoeffs(qa, qb, qc);
    return qc + alpha * (qb + qa * alpha);
}

double TwoFluidPipe::gasVelocity(double alpha, double s) const
{
    if (cfg_.closure.model == 0)
        return vM_ + alpha * cfg_.closure.slipVelocity;
    if (cfg_.closure.model == 2) // dynamic: vG = vM + alpha s
        return vM_ + alpha * s;
    return cfg_.closure.c0 * vM_ + cfg_.closure.driftVelocity;
}

void TwoFluidPipe::gasFluxCoeffs(double& qa, double& qb, double& qc) const
{
    // G(a) = alpha_g * vG(a): gas volume flux per unit area (algebraic
    // closures; the dynamic closure builds face-local coefficients from
    // the slip state)
    if (cfg_.closure.model == 0) // wallis: vG = vM + alpha*us
    {
        qa = -cfg_.closure.slipVelocity;
        qb = cfg_.closure.slipVelocity - vM_;
        qc = vM_;
    }
    else // zuber: vG = C0 vM + vD (constant)
    {
        qa = 0.0;
        qb = -(cfg_.closure.c0 * vM_ + cfg_.closure.driftVelocity);
        qc = cfg_.closure.c0 * vM_ + cfg_.closure.driftVelocity;
    }
}

double TwoFluidPipe::gasVolumeFlux(double alpha, double s) const
{
    if (cfg_.closure.model == 2) // dynamic: G = (1 - alpha)(vM + alpha s)
        return vM_ + (s - vM_) * alpha - s * alpha * alpha;
    double qa, qb, qc;
    gasFluxCoeffs(qa, qb, qc);
    return qc + alpha * (qb + qa * alpha);
}

double TwoFluidPipe::mixtureRhoCell(int i) const
{
    double rhoG = cfg_.fluid.rhoG;
    if (eosOn())
    {
        const double ag = 1.0 - alpha_[i];
        rhoG = ag > 1e-12 ? gasMass_[i] / (ag * dxs_[i]) : cfg_.fluid.rhoG;
    }
    return alpha_[i] * cfg_.fluid.rhoL + (1.0 - alpha_[i]) * rhoG;
}

std::vector<double> TwoFluidPipe::cellPressure() const
{
    std::vector<double> p(n_, 0.0);
    if (!eosOn()) return p;
    for (int i = 0; i < n_; ++i)
    {
        const double ag = 1.0 - alpha_[i];
        // a section with no gas volume has no meaningful pressure
        p[i] = ag > 1e-12 ? gasMass_[i] * cfg_.fluid.gasR
                            * sectionTemperature(i) / (ag * dxs_[i])
                          : 0.0;
    }
    return p;
}

double TwoFluidPipe::cellGasDensity(int i) const
{
    const double ag = 1.0 - alpha_[i];
    return eosOn() && ag > 1e-12 ? gasMass_[i] / (ag * dxs_[i])
                                 : cfg_.fluid.rhoG;
}

double TwoFluidPipe::ghostGasDensity(const TwoFluidBoundary& b) const
{
    if (!eosOn()) return cfg_.fluid.rhoG;
    return gasDensityOf(b.pressure,
                        b.temperature > 0 ? b.temperature : cfg_.temperature0);
}

double TwoFluidPipe::cellFriction(int i, double vM) const
{
    // homogeneous Churchill law on mixture properties; |dp| per cell
    const TwoFluidFluid& f = cfg_.fluid;
    const double a = alpha_[i];
    const double rho = mixtureRhoCell(i);
    const double mu = mixtureMu(f, a);
    const double v = std::abs(vM);
    if (v == 0) return 0.0;
    const double re = rho * v * diameter_ / mu;
    const double fr = churchillFriction(re, relrough_[i])
                    * (dxs_[i] / diameter_) * 0.5 * rho * v * v;
    return fr;
}

std::vector<double> TwoFluidPipe::cellPressureLosses(double vM) const
{
    // signed pressure loss per cell: friction (odd in vM) + gravity head
    std::vector<double> loss(n_);
    for (int i = 0; i < n_; ++i)
    {
        const double grav = mixtureRhoCell(i) * kGravity * sinT_[i] * dxs_[i];
        loss[i] = (vM >= 0 ? cellFriction(i, vM) : -cellFriction(i, vM)) + grav;
    }
    return loss;
}

double TwoFluidPipe::momentumScale(double dt) const
{
    double rhoBar = 0;
    for (int i = 0; i < n_; ++i) rhoBar += mixtureRhoCell(i);
    rhoBar /= n_;
    double gravSum = 0;
    for (int i = 0; i < n_; ++i)
        gravSum += mixtureRhoCell(i) * kGravity * sinT_[i] * dxs_[i];
    return rhoBar * length_ * (std::abs(vMOld_) + 1.0) / dt
         + std::abs(gravSum) + 1.0;
}

double TwoFluidPipe::momentumResidual(double dt, double u,
                                      double dpDrive) const
{
    // implicit Euler on the mixture momentum with the holdup frozen at
    // the old time level; the dynamic closure adds the lagged slip
    // inertia alpha alpha_g (rhoG - rhoL) ds/dt of the two phase
    // momenta (it vanishes at steady state). vM is a velocity, so the
    // inertia of the column is its mass per unit area rho_m L (not
    // divided by the area - that would be the volumetric-flow-rate form)
    double rhoBar = 0;
    for (int i = 0; i < n_; ++i) rhoBar += mixtureRhoCell(i);
    rhoBar /= n_;
    const double inertia = rhoBar * length_;
    double gravSum = 0;
    double slipInertia = 0;
    for (int i = 0; i < n_; ++i)
    {
        gravSum += mixtureRhoCell(i) * kGravity * sinT_[i] * dxs_[i];
        slipInertia += alpha_[i] * (1.0 - alpha_[i])
                     * (cellGasDensity(i) - cfg_.fluid.rhoL) * sDotOld_[i]
                     * dxs_[i];
    }
    double fr = 0;
    for (int i = 0; i < n_; ++i)
        fr += (u >= 0 ? cellFriction(i, u) : -cellFriction(i, u));
    return inertia * (u - vMOld_) / dt - slipInertia + fr + gravSum
         - dpDrive;
}

void TwoFluidPipe::momentumSolve(double dt)
{
    const TwoFluidBoundary& in = cfg_.inlet;
    const TwoFluidBoundary& out = cfg_.outlet;

    if (in.type == TwoFluidBoundary::Type::phaseRates)
    {
        vM_ = in.vsl + in.vsg;
        return;
    }
    if (in.type == TwoFluidBoundary::Type::sealed)
    {
        // sealed both ends: the mixture volume flux vanishes identically
        vM_ = 0.0;
        return;
    }

    const double dpDrive = in.pressure - out.pressure;
    auto residual = [&](double u)
    {
        return momentumResidual(dt, u, dpDrive);
    };

    double u = vMOld_;
    const double scale = momentumScale(dt) + std::abs(dpDrive);
    for (int it = 0; it < 60; ++it)
    {
        const double r = residual(u);
        if (std::abs(r) <= 1e-12 * scale) break;
        const double h = std::max(1e-6 * std::max(std::abs(u), 1e-3), 1e-12);
        const double slope = (residual(u + h) - residual(u - h)) / (2.0 * h);
        if (!(slope > 0))
            cfFatal("twoFluid: non-monotone momentum residual");
        u -= r / slope;
    }
    if (!std::isfinite(u)) cfFatal("twoFluid: momentum solve diverged");
    vM_ = u;
}

void TwoFluidPipe::slipStep(double dt)
{
    if (cfg_.closure.model != 2) return;

    // reduced two-momentum slip ODE per section (holdup frozen):
    //   K_s ds/dt = (1/rhoG - 1/rhoL)(rho_m du/dt + loss/dx) - D sigma,
    // sigma = C s |s| the interfacial drag. The drag is semi-implicit as
    // the secant through the origin, which is unconditionally stable and
    // lands exactly on the terminal slip s* = sign(R) sqrt(|R|/(C D)).
    const TwoFluidFluid& f = cfg_.fluid;
    const double aM = (vM_ - vMOld_) / dt;
    std::vector<double> old = slip_;
    for (int i = 0; i < n_; ++i)
    {
        const double a = alpha_[i];
        const double ag = 1.0 - a;
        if (a <= 1e-12 || ag <= 1e-12)
        {
            // single-phase cell: no interface, no drag, no slip
            slip_[i] = 0.0;
            sDotOld_[i] = -old[i] / dt;
            continue;
        }
        // per-section gas density: with the EOS the slip inertia factor,
        // drive and drag denominator are section-local
        const double rhoG = cellGasDensity(i);
        const double invRhoDiff = 1.0 / rhoG - 1.0 / f.rhoL;
        const double drho = rhoG - f.rhoL;
        const double ks = 1.0 + a * ag * drho * drho / (f.rhoL * rhoG);
        const double loss = (vM_ >= 0 ? cellFriction(i, vM_)
                                      : -cellFriction(i, vM_))
                          + mixtureRhoCell(i) * kGravity * sinT_[i] * dxs_[i];
        const double drive = invRhoDiff
                           * (mixtureRhoCell(i) * aM + loss / dxs_[i]);
        const double d = 1.0 / (f.rhoL * a) + 1.0 / (rhoG * ag);
        const double rate = cfg_.closure.interfacialDrag
                          * std::abs(old[i]) * d / ks;
        slip_[i] = (old[i] / dt + drive / ks) / (1.0 / dt + rate);
        sDotOld_[i] = (slip_[i] - old[i]) / dt;
    }
}

void TwoFluidPipe::transportStep(double dt)
{
    double qa = 0, qb = 0, qc = 0;
    if (cfg_.closure.model != 2) fluxCoeffs(qa, qb, qc);

    // per-face slip of the dynamic closure (ghost slip at the open
    // boundaries, cell pair average inside); unused by the algebraic
    // closures
    auto faceSlip = [&](int j)
    {
        return j == 0 ? slip_[0]
             : (j == n_ ? slip_[n_ - 1]
                        : 0.5 * (slip_[j - 1] + slip_[j]));
    };

    // CFL guard for the explicit Godunov transport (monotone under
    // dt * max|F'| / dx <= 1; F' is monotone so the extremes sit at 0 and 1;
    // the dynamic closure has a face-local flux, so take the max over faces)
    double cmax = 0;
    if (cfg_.closure.model == 2)
    {
        for (int j = 0; j <= n_; ++j)
        {
            const double s = faceSlip(j);
            cmax = std::max(cmax, std::max(std::abs(vM_ - s),
                                           std::abs(vM_ + s)));
        }
    }
    else
    {
        cmax = std::max(std::abs(qb), std::abs(2.0 * qa + qb));
    }
    const double dxMin = *std::min_element(dxs_.begin(), dxs_.end());
    if (dt * cmax / dxMin > 1.0 + 1e-12)
        cfFatal("twoFluid: CFL violated (dt * |dF/dalpha| / dx = "
                + std::to_string(dt * cmax / dxMin)
                + ") - reduce dt");

    const TwoFluidBoundary& in = cfg_.inlet;
    const TwoFluidBoundary& out = cfg_.outlet;
    std::vector<double> flux(n_ + 1), wflux(n_ + 1);
    std::vector<double> hflux(n_ + 1, 0.0);
    std::vector<double> gvol(n_ + 1, 0.0), gmflux(n_ + 1, 0.0);

    for (int j = 0; j <= n_; ++j)
    {
        if (j == 0)
        {
            if (in.type == TwoFluidBoundary::Type::sealed) flux[0] = 0;
            else if (in.type == TwoFluidBoundary::Type::phaseRates)
                flux[0] = in.vsl; // injected liquid superficial velocity
            else if (cfg_.closure.model == 2)
            {
                const double s = faceSlip(0);
                flux[0] = godunovQuadratic(s, vM_ - s, 0.0, in.alpha,
                                           alpha_[0]);
            }
            else
                flux[0] = godunovQuadratic(qa, qb, qc, in.alpha, alpha_[0]);
        }
        else if (j == n_)
        {
            if (out.type == TwoFluidBoundary::Type::sealed) flux[n_] = 0;
            else if (out.type == TwoFluidBoundary::Type::phaseRates)
                flux[n_] = out.vsl; // junction-imposed split (supplied
                                    // branch fed through its outlet)
            else
                flux[n_] = liquidFlux(alpha_[n_ - 1], faceSlip(n_));
        }
        else if (cfg_.closure.model == 2)
        {
            const double s = faceSlip(j);
            flux[j] = godunovQuadratic(s, vM_ - s, 0.0, alpha_[j - 1],
                                       alpha_[j]);
        }
        else
        {
            flux[j] = godunovQuadratic(qa, qb, qc, alpha_[j - 1], alpha_[j]);
        }

        // water cut rides the liquid flux, upwinded on the flux direction
        if (flux[j] == 0) wflux[j] = 0;
        else if (j == 0)
        {
            if (in.type == TwoFluidBoundary::Type::phaseRates)
                wflux[0] = in.waterFraction * flux[0];
            else
                wflux[0] = (flux[0] >= 0 ? in.waterFraction : water_[0])
                         * flux[0];
        }
        else if (j == n_)
        {
            if (out.type == TwoFluidBoundary::Type::phaseRates)
                wflux[n_] = out.waterFraction * flux[n_];
            else
                wflux[n_] = water_[n_ - 1] * flux[n_];
        }
        else
        {
            wflux[j] = (flux[j] >= 0 ? water_[j - 1] : water_[j]) * flux[j];
        }

        // hydrate fraction rides the liquid flux beside the water cut,
        // same upwind rule (the inflow composition carries hydrateFraction)
        if (!cfg_.hydrateOn) continue;
        if (flux[j] == 0) hflux[j] = 0;
        else if (j == 0)
        {
            if (in.type == TwoFluidBoundary::Type::phaseRates)
                hflux[0] = in.hydrateFraction * flux[0];
            else
                hflux[0] = (flux[0] >= 0 ? in.hydrateFraction : hydrate_[0])
                         * flux[0];
        }
        else if (j == n_)
        {
            if (out.type == TwoFluidBoundary::Type::phaseRates)
                hflux[n_] = out.hydrateFraction * flux[n_];
            else
                hflux[n_] = hydrate_[n_ - 1] * flux[n_];
        }
        else
        {
            hflux[j] = (flux[j] >= 0 ? hydrate_[j - 1] : hydrate_[j])
                     * flux[j];
        }
    }

    // gas volume flux per face: the gas counterpart of the liquid flux,
    // quadratic/linear in alpha for every closure with G(1) = 0, so the
    // same Godunov machinery and the same CFL bound apply. The volume flux
    // exists regardless of the EOS (the thermal step advects gas enthalpy
    // with it); with the EOS on the mass flux upwinds the density from the
    // side the volume flux comes from: the injected reference density at a
    // rate boundary, the ghost density at a pressure boundary, the cell
    // density inside.
    for (int j = 0; j <= n_; ++j)
    {
        if (j == 0)
        {
            if (in.type == TwoFluidBoundary::Type::phaseRates)
                gvol[0] = in.vsg;
            else if (in.type == TwoFluidBoundary::Type::fixedPressure)
            {
                if (cfg_.closure.model == 2)
                {
                    const double s = faceSlip(0);
                    gvol[0] = godunovQuadratic(-s, s - vM_, vM_, in.alpha,
                                               alpha_[0]);
                }
                else
                {
                    double gqa, gqb, gqc;
                    gasFluxCoeffs(gqa, gqb, gqc);
                    gvol[0] = godunovQuadratic(gqa, gqb, gqc, in.alpha,
                                               alpha_[0]);
                }
            }
            // sealed inlet: gvol[0] stays 0
        }
        else if (j == n_)
        {
            if (out.type == TwoFluidBoundary::Type::phaseRates)
                gvol[n_] = out.vsg; // junction-imposed split
            else if (out.type != TwoFluidBoundary::Type::sealed)
                gvol[n_] = gasVolumeFlux(alpha_[n_ - 1], faceSlip(n_));
        }
        else if (cfg_.closure.model == 2)
        {
            const double s = faceSlip(j);
            gvol[j] = godunovQuadratic(-s, s - vM_, vM_, alpha_[j - 1],
                                       alpha_[j]);
        }
        else
        {
            double gqa, gqb, gqc;
            gasFluxCoeffs(gqa, gqb, gqc);
            gvol[j] = godunovQuadratic(gqa, gqb, gqc, alpha_[j - 1],
                                       alpha_[j]);
        }

        if (eosOn())
        {
            // rate-type faces inject the boundary's gas density: the
            // fluid reference rhoG, or the junction mix (rhoGInjected)
            const double rhoIn0 = cfg_.inlet.rhoGInjected > 0
                                      ? cfg_.inlet.rhoGInjected
                                      : cfg_.fluid.rhoG;
            const double rhoInN = cfg_.outlet.rhoGInjected > 0
                                      ? cfg_.outlet.rhoGInjected
                                      : cfg_.fluid.rhoG;
            const double rhoUp =
                gvol[j] >= 0
                    ? (j == 0
                           ? (in.type == TwoFluidBoundary::Type::phaseRates
                                  ? rhoIn0
                                  : ghostGasDensity(in))
                           : j == n_ ? cellGasDensity(n_ - 1)
                                     : cellGasDensity(j - 1))
                    : (j == n_ ? (out.type == TwoFluidBoundary::Type::phaseRates
                                      ? rhoInN
                                      : ghostGasDensity(out))
                               : cellGasDensity(j));
            gmflux[j] = rhoUp * gvol[j];
        }
    }
    liquidFaceFlux_ = flux;
    gasFaceFlux_ = gvol;

    static constexpr double kTol = 1e-9;
    for (int i = 0; i < n_; ++i)
    {
        const double aNew = alpha_[i] + dt / dxs_[i] * (flux[i] - flux[i + 1]);
        if (aNew < -kTol || aNew > 1.0 + kTol)
            cfFatal("twoFluid: holdup left [0,1] - reduce dt");
        const double afOld = alpha_[i] * water_[i];
        const double ahOld = cfg_.hydrateOn ? alpha_[i] * hydrate_[i] : 0.0;
        double afNew = afOld + dt / dxs_[i] * (wflux[i] - wflux[i + 1]);
        const double aCap = std::max(aNew, 0.0);
        if (afNew < -kTol || afNew > aCap + kTol)
            cfFatal("twoFluid: water cut left [0, holdup] - reduce dt");
        alpha_[i] = std::min(std::max(aNew, 0.0), 1.0);
        afNew = std::min(std::max(afNew, 0.0), aCap);
        water_[i] = aCap > kTol ? afNew / aCap : 0.0;
        if (cfg_.hydrateOn)
        {
            // hydrate rides the liquid with the same upwind rule and the
            // same clamp semantics as the water cut (ahOld was taken
            // against the pre-transport holdup, beside afOld)
            double ahNew = ahOld
                         + dt / dxs_[i] * (hflux[i] - hflux[i + 1]);
            ahNew = std::min(std::max(ahNew, 0.0), aCap);
            hydrate_[i] = aCap > kTol ? ahNew / aCap : 0.0;
        }
        if (eosOn())
        {
            // gas mass continuity on the per-unit-area inventory
            // rhoG (1 - alpha) dx: the mass-flux divergence changes it
            // directly (no dx factor); upwinding keeps it non-negative
            const double agNew = 1.0 - alpha_[i];
            const double numer
                = gasMass_[i] - dt * (gmflux[i + 1] - gmflux[i]);
            const double trace = kTol * cfg_.fluid.rhoG * dxs_[i];
            if (agNew <= 1e-12 && numer > trace)
                cfFatal("twoFluid: gas sealed in a closing cell - reduce dt");
            if (numer < -trace)
                cfFatal("twoFluid: gas mass went negative - reduce dt");
            gasMass_[i] = std::max(numer, 0.0);
        }
    }
}

void TwoFluidPipe::hydrateStep(double dt)
{
    if (!cfg_.hydrateOn) return;
    static constexpr double kTol = 1e-9;
    const TwoFluidHydrate& hy = cfg_.hydrate;
    // the subcooling uses the just-transported state: the EOS pressure of
    // the transported gas inventory and the current section temperature
    const std::vector<double> p = cellPressure();
    for (int i = 0; i < n_; ++i)
    {
        const double ag = 1.0 - alpha_[i];
        if (ag <= 1e-12)
            cfFatal("twoFluid: hydrate kinetics needs gas volume in every"
                    " cell - keep a gas-bearing holdup");
        const double dT = equilibriumTemperature(p[i]) - sectionTemperature(i);
        const double hOld = hydrate_[i];
        // hydrate-able slurry material: free water plus existing hydrate
        // (oil is inert); the growth law converges to the full conversion
        const double total = water_[i] + hOld;
        double hNew;
        if (dT >= 0)
        {
            // implicit Euler on dh/dt = kg dT (total - h)
            hNew = (hOld / dt + hy.growthRate * dT * total)
                 / (1.0 / dt + hy.growthRate * dT);
            // gas availability: gasContent kg per m3 of hydrate off the EOS
            // inventory, capped so the mass stays non-negative
            const double trace = kTol * cfg_.fluid.rhoG * dxs_[i];
            const double dhMax = (gasMass_[i] - trace)
                               / (hy.gasContent * alpha_[i] * dxs_[i]);
            hNew = std::min(hNew, hOld + std::max(dhMax, 0.0));
        }
        else
        {
            // implicit Euler on dh/dt = -kd |dT| h
            hNew = hOld / (1.0 + dt * hy.dissociationRate * (-dT));
        }
        hNew = std::min(std::max(hNew, 0.0), total);
        const double dh = hNew - hOld;
        // 1:1 volume move between free water and hydrate (the
        // volume-preserving convention of the network kernel)
        water_[i] = std::min(std::max(total - hNew, 0.0), 1.0);
        hydrate_[i] = hNew;
        // stoichiometric gas mass sink on the per-area EOS inventory
        // (reappears on dissociation)
        gasMass_[i] -= hy.gasContent * alpha_[i] * dh * dxs_[i];
        // hydrate volume change per bulk m3: the latent-heat source the
        // thermal step consumes (positive = release)
        hydrateDelta_[i] = alpha_[i] * dh;
    }
}

void TwoFluidPipe::step(double dt)
{
    if (!(dt > 0)) cfFatal("twoFluid: dt must be positive");
    // junction branches have their mixture velocity imposed by the
    // network's coupled Newton (externalMomentum_ is a one-shot flag)
    if (externalMomentum_) externalMomentum_ = false;
    else momentumSolve(dt);
    slipStep(dt);
    transportStep(dt);
    hydrateStep(dt);
    thermalStep(dt);
    vMOld_ = vM_;
    t_ += dt;
}

double TwoFluidPipe::wallTemperature(int i) const
{
    const TwoFluidThermal& th = cfg_.thermal;
    if (th.wallModel == 0) return th.wallTemperature;
    const double depth = th.inletDepth - centerDrop_[i];
    return th.surfaceTemperature + th.geothermalGradient * depth;
}

void TwoFluidPipe::thermalStep(double dt)
{
    if (!cfg_.thermalOn) return;
    const TwoFluidFluid& f = cfg_.fluid;
    const TwoFluidThermal& th = cfg_.thermal;
    const double area = M_PI * diameter_ * diameter_ / 4.0;

    // implicit Euler on the section temperatures. The phase enthalpies ride
    // the phase volume fluxes of this step's transport - with slip the
    // liquid no longer moves at the mixture velocity, so a mixture-rhoc * vM
    // weighting would over-drag the liquid capacity against the slip - and
    // each phase upwinds its temperature on its own flux direction. The gas
    // density upwinds like the mass flux does (injected reference density at
    // a rate boundary, ghost density at a pressure boundary, cell density
    // inside).
    std::vector<std::vector<double>> mat(n_, std::vector<double>(n_, 0.0));
    std::vector<double> rhs(n_, 0.0);
    for (int i = 0; i < n_; ++i)
    {
        const double cap = area * dxs_[i]
                         * (alpha_[i] * f.rhoL * f.cL
                            + (1.0 - alpha_[i]) * cellGasDensity(i) * f.cG);
        const double kw = th.wallHeatTransfer * dxs_[i];
        mat[i][i] = cap / dt + kw;
        rhs[i] = cap / dt * T_[i] + kw * wallTemperature(i);
        if (cfg_.hydrateOn)
        {
            // hydrate formation releases its latent heat into the section
            // (dissociation absorbs it), J/s over the cell bulk volume
            rhs[i] += cfg_.hydrate.latentHeat * hydrateDelta_[i] * area
                    * dxs_[i] / dt;
        }
    }

    // rate-type faces inject the boundary's gas density: the fluid
    // reference rhoG, or the junction mix (rhoGInjected)
    auto injectedRhoG = [&](const TwoFluidBoundary& b) -> double
    {
        return b.rhoGInjected > 0 ? b.rhoGInjected : f.rhoG;
    };
    auto gasRho = [&](int j) -> double
    {
        if (!eosOn()) return f.rhoG;
        const double q = gasFaceFlux_[j];
        if (q > 0)
            return j == 0
                ? (cfg_.inlet.type == TwoFluidBoundary::Type::phaseRates
                       ? injectedRhoG(cfg_.inlet)
                       : ghostGasDensity(cfg_.inlet))
                : j == n_ ? cellGasDensity(n_ - 1) : cellGasDensity(j - 1);
        if (q < 0)
            return j == n_
                ? (cfg_.outlet.type == TwoFluidBoundary::Type::phaseRates
                       ? injectedRhoG(cfg_.outlet)
                       : ghostGasDensity(cfg_.outlet))
                : cellGasDensity(j);
        return 0.0;
    };
    // one face-phase contribution: the upstream side drains at its own
    // temperature (diagonal), the downstream side receives the upwind
    // temperature - an implicit coupling to the upstream cell, or the
    // boundary temperature of the ghost reservoir
    auto faceTerm = [&](int j, double q, double k)
    {
        if (q == 0) return;
        const int up = q > 0 ? j - 1 : j;
        const int dn = q > 0 ? j : j - 1;
        if (up >= 0 && up < n_) mat[up][up] += k;
        if (dn >= 0 && dn < n_)
        {
            if (up >= 0 && up < n_)
                mat[dn][up] -= k;
            else
                rhs[dn] += k * (j == 0 ? cfg_.inlet.temperature
                                       : cfg_.outlet.temperature);
        }
    };

    inletHeatRate_ = 0;
    for (int j = 0; j <= n_; ++j)
    {
        const double kl = area * f.rhoL * f.cL
                        * std::abs(liquidFaceFlux_[j]);
        const double kg = area * gasRho(j) * f.cG
                        * std::abs(gasFaceFlux_[j]);
        faceTerm(j, liquidFaceFlux_[j], kl);
        faceTerm(j, gasFaceFlux_[j], kg);
        if (j == 0)
        {
            if (liquidFaceFlux_[j] > 0)
                inletHeatRate_ += kl * cfg_.inlet.temperature;
            if (gasFaceFlux_[j] > 0)
                inletHeatRate_ += kg * cfg_.inlet.temperature;
        }
        else if (j == n_)
        {
            if (liquidFaceFlux_[j] < 0)
                inletHeatRate_ += kl * cfg_.outlet.temperature;
            if (gasFaceFlux_[j] < 0)
                inletHeatRate_ += kg * cfg_.outlet.temperature;
        }
    }
    detail::solveDenseLinear(mat, rhs);

    wallHeatRate_ = 0;
    for (int i = 0; i < n_; ++i)
    {
        wallHeatRate_ += th.wallHeatTransfer * dxs_[i]
                       * (wallTemperature(i) - rhs[i]);
        T_[i] = rhs[i];
    }
}

double TwoFluidPipe::pressureDrop() const
{
    const std::vector<double> p = nodePressure();
    return p.front() - p.back();
}

std::vector<double> TwoFluidPipe::nodePressure() const
{
    const std::vector<double> loss = cellPressureLosses(vM_);
    std::vector<double> p(n_ + 1, 0.0);
    if (cfg_.inlet.type == TwoFluidBoundary::Type::fixedPressure)
    {
        p[0] = cfg_.inlet.pressure;
        for (int i = 0; i < n_; ++i) p[i + 1] = p[i] - loss[i];
    }
    else
    {
        // rate-driven or sealed: gauge pressure relative to the outlet
        p[n_] = cfg_.outlet.type == TwoFluidBoundary::Type::fixedPressure
                    ? cfg_.outlet.pressure
                    : 0.0;
        for (int i = n_ - 1; i >= 0; --i) p[i] = p[i + 1] + loss[i];
    }
    return p;
}

std::vector<double> TwoFluidPipe::cellWater() const { return water_; }

std::vector<double> TwoFluidPipe::cellVsl() const
{
    // superficial liquid velocity = liquid volume flux per unit area
    std::vector<double> v(n_);
    for (int i = 0; i < n_; ++i)
        v[i] = liquidFlux(alpha_[i], slip_[i]);
    return v;
}

std::vector<double> TwoFluidPipe::cellVsg() const
{
    std::vector<double> v(n_);
    for (int i = 0; i < n_; ++i)
        v[i] = (1.0 - alpha_[i]) * gasVelocity(alpha_[i], slip_[i]);
    return v;
}

namespace
{

TwoFluidBoundary loadTwoFluidBoundary(const Dictionary& d)
{
    TwoFluidBoundary b;
    const std::string type = d.getWord("type");
    if (type == "sealed") b.type = TwoFluidBoundary::Type::sealed;
    else if (type == "phaseRates")
        b.type = TwoFluidBoundary::Type::phaseRates;
    else if (type == "fixedPressure")
        b.type = TwoFluidBoundary::Type::fixedPressure;
    else cfFatal("twoFluid: unknown boundary type '" + type + "'");
    b.vsl = d.getScalarOr("vsl", 0.0);
    b.vsg = d.getScalarOr("vsg", 0.0);
    b.pressure = d.getScalarOr("p", 0.0);
    b.alpha = d.getScalarOr("alpha", 0.0);
    b.waterFraction = d.getScalarOr("waterFraction", 0.0);
    b.hydrateFraction = d.getScalarOr("hydrateFraction", 0.0);
    b.temperature = d.getScalarOr("temperature", 0.0);
    return b;
}

// the physics blocks shared by the serial case and the branch network:
// fluid / closure / initial / thermal / hydrate
void loadTwoFluidCommon(const Dictionary& caseDict, TwoFluidCase& c)
{
    if (!caseDict.has("fluid")) cfFatal("twoFluid: case has no fluid{} block");
    const Dictionary& fd = caseDict.getDict("fluid");
    c.fluid.rhoL = fd.getScalar("rhoL");
    c.fluid.rhoG = fd.getScalar("rhoG");
    c.fluid.muL = fd.getScalar("muL");
    c.fluid.muG = fd.getScalar("muG");
    c.fluid.cL = fd.getScalarOr("cL", 0.0);
    c.fluid.cG = fd.getScalarOr("cG", 0.0);
    c.fluid.gasR = fd.getScalarOr("gasR", 0.0);

    if (!caseDict.has("closure"))
        cfFatal("twoFluid: case has no closure{} block");
    const Dictionary& cd = caseDict.getDict("closure");
    const std::string model = cd.getWordOr("model", "wallis");
    if (model == "wallis") c.closure.model = 0;
    else if (model == "zuber") c.closure.model = 1;
    else if (model == "dynamic") c.closure.model = 2;
    else cfFatal("twoFluid: unknown closure model '" + model + "'");
    c.closure.slipVelocity = cd.getScalarOr("slipVelocity", 0.0);
    c.closure.c0 = cd.getScalarOr("C0", 1.0);
    c.closure.driftVelocity = cd.getScalarOr("driftVelocity", 0.0);
    c.closure.interfacialDrag = cd.getScalarOr("interfacialDrag", 0.0);

    if (!caseDict.has("initial"))
        cfFatal("twoFluid: case has no initial{} block");
    const Dictionary& id = caseDict.getDict("initial");
    c.alpha0 = id.getScalarOr("alpha", 0.0);
    c.water0 = id.getScalarOr("water", 0.0);
    c.temperature0 = id.getScalarOr("temperature", 300.0);
    if (c.alpha0 < 0 || c.alpha0 > 1 || c.water0 < 0 || c.water0 > 1)
        cfFatal("twoFluid: initial alpha/water must be in [0,1]");
    if (id.has("fills"))
    {
        const Dictionary& fills = id.getDict("fills");
        for (const std::string& name : fills.keys())
        {
            const Dictionary& f = fills.getDict(name);
            TwoFluidFill fill;
            fill.from = f.getScalar("from");
            fill.to = f.getScalar("to");
            fill.alpha = f.getScalar("alpha");
            c.fills.push_back(fill);
        }
    }

    if (caseDict.has("thermal"))
    {
        c.thermalOn = true;
        const Dictionary& th = caseDict.getDict("thermal");
        c.thermal.wallHeatTransfer = th.getScalar("wallHeatTransfer");
        const std::string wm = th.getWordOr("wallModel", "constant");
        if (wm == "constant") c.thermal.wallModel = 0;
        else if (wm == "geothermal") c.thermal.wallModel = 1;
        else cfFatal("twoFluid: unknown thermal wallModel '" + wm + "'");
        c.thermal.wallTemperature = th.getScalarOr("wallTemperature", 0.0);
        c.thermal.surfaceTemperature
            = th.getScalarOr("surfaceTemperature", 0.0);
        c.thermal.geothermalGradient
            = th.getScalarOr("gradient", 0.0);
        c.thermal.inletDepth = th.getScalarOr("inletDepth", 0.0);
        if (c.thermal.wallHeatTransfer < 0)
            cfFatal("twoFluid: negative wallHeatTransfer");
    }

    if (caseDict.has("hydrate"))
    {
        c.hydrateOn = true;
        const Dictionary& hd = caseDict.getDict("hydrate");
        c.hydrate.curveA = hd.getScalar("curveA");
        c.hydrate.curveB = hd.getScalar("curveB");
        c.hydrate.growthRate = hd.getScalar("growthRate");
        c.hydrate.dissociationRate = hd.getScalar("dissociationRate");
        c.hydrate.gasContent = hd.getScalar("gasContent");
        c.hydrate.latentHeat = hd.getScalar("latentHeat");
        c.hydrate.initial = hd.getScalarOr("initial", 0.0);
    }
}

} // namespace

TwoFluidCase loadTwoFluidCase(const Dictionary& caseDict)
{
    TwoFluidCase c;
    loadTwoFluidCommon(caseDict, c);

    if (caseDict.has("pipe"))
    {
        const Dictionary& pd = caseDict.getDict("pipe");
        TwoFluidGeo g;
        g.name = pd.getWordOr("name", "pipe");
        g.length = pd.getScalar("length");
        g.diameter = pd.getScalar("diameter");
        g.roughness = pd.getScalarOr("roughness", 0.0);
        g.angleDeg = pd.getScalarOr("angleDeg", 0.0);
        g.sections = static_cast<int>(pd.getScalarOr("sections", 50.0));
        c.pipes.push_back(g);
    }
    else if (caseDict.has("pipes"))
    {
        // named sub-blocks in first-occurrence order form the serial chain
        const Dictionary& ps = caseDict.getDict("pipes");
        for (const std::string& k : ps.keys())
        {
            const Dictionary& pd = ps.getDict(k);
            TwoFluidGeo g;
            g.name = k;
            g.length = pd.getScalar("length");
            g.diameter = pd.getScalar("diameter");
            g.roughness = pd.getScalarOr("roughness", 0.0);
            g.angleDeg = pd.getScalarOr("angleDeg", 0.0);
            g.sections = static_cast<int>(pd.getScalarOr("sections", 50.0));
            c.pipes.push_back(g);
        }
    }
    else
    {
        cfFatal("twoFluid: case has no pipe{} or pipes{} block");
    }

    if (!caseDict.has("boundaries"))
        cfFatal("twoFluid: case has no boundaries{} block");
    const Dictionary& bd = caseDict.getDict("boundaries");
    if (!bd.has("inlet") || !bd.has("outlet"))
        cfFatal("twoFluid: boundaries{} needs inlet{} and outlet{}");
    c.inlet = loadTwoFluidBoundary(bd.getDict("inlet"));
    c.outlet = loadTwoFluidBoundary(bd.getDict("outlet"));

    return c;
}

TwoFluidNetworkCase loadTwoFluidNetworkCase(const Dictionary& caseDict)
{
    TwoFluidNetworkCase nc;
    nc.perPipe.networkMode = true;
    loadTwoFluidCommon(caseDict, nc.perPipe);

    if (!caseDict.has("network"))
        cfFatal("twoFluid: network case has no network{} block");
    const Dictionary& nd = caseDict.getDict("network");
    if (!nd.has("pipes")) cfFatal("twoFluid: network{} needs pipes{}");
    const Dictionary& ps = nd.getDict("pipes");
    for (const std::string& k : ps.keys())
    {
        const Dictionary& pd = ps.getDict(k);
        TwoFluidGeo g;
        g.name = k;
        g.length = pd.getScalar("length");
        g.diameter = pd.getScalar("diameter");
        g.roughness = pd.getScalarOr("roughness", 0.0);
        g.angleDeg = pd.getScalarOr("angleDeg", 0.0);
        g.sections = static_cast<int>(pd.getScalarOr("sections", 50.0));
        nc.pipes.push_back(g);
    }
    if (!nd.has("junctions")) cfFatal("twoFluid: network{} needs junctions{}");
    const Dictionary& js = nd.getDict("junctions");
    for (const std::string& k : js.keys())
    {
        const Dictionary& jd = js.getDict(k);
        TwoFluidJunctionSpec s;
        s.name = k;
        s.ends = jd.getWordList("ends");
        if (s.ends.size() < 2)
            cfFatal("twoFluid: junction '" + k + "' needs at least two ends");
        nc.junctions.push_back(s);
    }
    if (nd.has("boundaries"))
    {
        const Dictionary& bd = nd.getDict("boundaries");
        for (const std::string& k : bd.keys())
        {
            if (k.rfind('.') == std::string::npos)
                cfFatal("twoFluid: network boundary key '" + k
                        + "' must be pipe.inlet or pipe.outlet");
            nc.boundaries[k] = loadTwoFluidBoundary(bd.getDict(k));
        }
    }

    // topology validation: every pipe exactly one junction end, that end
    // boundary-free, the other end carrying a boundary; a rate-driven far
    // end must inject toward the junction; the dynamic closure does not
    // support junctions yet
    if (nc.perPipe.closure.model == 2)
        cfFatal("twoFluid: the dynamic closure does not support junctions"
                " yet");
    std::map<std::string, char> junctionEnd; // "pipe.inlet" -> i/o
    for (const TwoFluidJunctionSpec& j : nc.junctions)
    {
        for (const std::string& ref : j.ends)
        {
            const auto dot = ref.rfind('.');
            if (dot == std::string::npos)
                cfFatal("twoFluid: junction end '" + ref
                        + "' must be pipe.inlet or pipe.outlet");
            const std::string p = ref.substr(0, dot);
            const std::string e = ref.substr(dot + 1);
            if (e != "inlet" && e != "outlet")
                cfFatal("twoFluid: junction end '" + ref + "' must end in"
                        " .inlet or .outlet");
            bool found = false;
            for (const TwoFluidGeo& g : nc.pipes) found |= g.name == p;
            if (!found)
                cfFatal("twoFluid: junction '" + j.name + "' references"
                        " unknown pipe '" + p + "'");
            const std::string key = p + "." + e;
            if (junctionEnd.count(key))
                cfFatal("twoFluid: pipe end '" + key + "' joins two"
                        " junctions");
            junctionEnd[key] = e == "outlet" ? 'o' : 'i';
        }
    }
    for (const TwoFluidGeo& g : nc.pipes)
    {
        const bool inJ = junctionEnd.count(g.name + ".inlet") > 0;
        const bool outJ = junctionEnd.count(g.name + ".outlet") > 0;
        if (inJ == outJ)
            cfFatal("twoFluid: pipe '" + g.name
                    + "' needs exactly one junction end");
        const std::string farKey = g.name + (inJ ? ".outlet" : ".inlet");
        if (!nc.boundaries.count(farKey))
            cfFatal("twoFluid: pipe '" + g.name + "' end '" + farKey
                    + "' has no boundary and no junction");
        const std::string jKey = g.name + (inJ ? ".inlet" : ".outlet");
        if (nc.boundaries.count(jKey))
            cfFatal("twoFluid: junction end '" + jKey
                    + "' must not also carry a boundary");
        const TwoFluidBoundary& far = nc.boundaries.at(farKey);
        if (far.type == TwoFluidBoundary::Type::phaseRates)
        {
            if (far.vsl + far.vsg < 0)
                cfFatal("twoFluid: phaseRates injection with negative total"
                        " rate is not supported");
            if (!outJ)
                cfFatal("twoFluid: pipe '" + g.name + "': a rate-driven far"
                        " end needs the junction at its outlet");
        }
        if (nc.perPipe.thermalOn
            && far.type != TwoFluidBoundary::Type::sealed
            && !(far.temperature > 0))
            cfFatal("twoFluid: pipe '" + g.name
                    + "': thermal on needs a far-end temperature");
    }

    return nc;
}

TwoFluidNetwork::TwoFluidNetwork(const TwoFluidNetworkCase& c)
    : cfg_(c)
{
    for (int i = 0; i < static_cast<int>(cfg_.pipes.size()); ++i)
        pipeIndex_[cfg_.pipes[i].name] = i;

    for (const TwoFluidJunctionSpec& js : cfg_.junctions)
    {
        Junction jn;
        jn.name = js.name;
        jn.alphaJ = cfg_.perPipe.alpha0;
        jn.waterJ = cfg_.perPipe.water0;
        jn.hydrateJ = cfg_.perPipe.hydrateOn ? cfg_.perPipe.hydrate.initial
                                             : 0.0;
        jn.tJ = cfg_.perPipe.temperature0;
        jn.rhoJ = cfg_.perPipe.fluid.rhoG;
        bool pressureRef = false;
        for (const std::string& ref : js.ends)
        {
            const auto dot = ref.rfind('.');
            EndRef e;
            e.pipe = pipeIndex_.at(ref.substr(0, dot));
            e.atOutlet = ref.substr(dot + 1) == "outlet";
            const std::string farKey = cfg_.pipes[e.pipe].name
                                     + (e.atOutlet ? ".inlet" : ".outlet");
            e.far = cfg_.boundaries.at(farKey);
            pressureRef |= e.far.type
                        == TwoFluidBoundary::Type::fixedPressure;
            jn.ends.push_back(e);
        }
        if (!pressureRef)
            cfFatal("twoFluid: junction '" + jn.name + "' has no pressure"
                    " reference (needs a fixed-pressure far end)");
        // initial junction pressure: the steady conductance-weighted mean
        // of the far fixed pressures (the t = 0 report and the Newton
        // warm start)
        const TwoFluidFluid& f = cfg_.perPipe.fluid;
        const double a0 = cfg_.perPipe.alpha0;
        const double muBar = a0 * f.muL + (1.0 - a0) * f.muG;
        double num = 0, den = 0;
        for (const EndRef& e : jn.ends)
        {
            if (e.far.type != TwoFluidBoundary::Type::fixedPressure)
                continue;
            const TwoFluidGeo& g = cfg_.pipes[e.pipe];
            const double cond = M_PI * g.diameter * g.diameter / 4.0
                              * g.diameter * g.diameter
                              / (32.0 * muBar * g.length);
            num += cond * e.far.pressure;
            den += cond;
        }
        jn.p = num / den;
        junctions_.push_back(jn);
    }

    for (int i = 0; i < static_cast<int>(cfg_.pipes.size()); ++i)
    {
        const TwoFluidGeo& g = cfg_.pipes[i];
        TwoFluidCase pc = cfg_.perPipe;
        pc.pipes.assign(1, g);
        // the junction end gets a placeholder ghost (fixed pressure at the
        // initial state) that the per-step junction solve overwrites
        TwoFluidBoundary junctionGhost;
        junctionGhost.type = TwoFluidBoundary::Type::fixedPressure;
        junctionGhost.pressure = 0;
        junctionGhost.alpha = pc.alpha0;
        junctionGhost.waterFraction = pc.water0;
        junctionGhost.hydrateFraction
            = pc.hydrateOn ? pc.hydrate.initial : 0.0;
        junctionGhost.temperature = pc.temperature0;
        bool jAtOutlet = false;
        for (const Junction& jn : junctions_)
            for (const EndRef& e : jn.ends)
                if (e.pipe == i) jAtOutlet = e.atOutlet;
        const std::string farKey = g.name + (jAtOutlet ? ".inlet" : ".outlet");
        if (jAtOutlet)
        {
            pc.inlet = cfg_.boundaries.at(farKey);
            pc.outlet = junctionGhost;
        }
        else
        {
            pc.inlet = junctionGhost;
            pc.outlet = cfg_.boundaries.at(farKey);
        }
        pipes_.push_back(std::make_unique<TwoFluidPipe>(pc));
    }
}

void TwoFluidNetwork::step(double dt)
{
    if (!(dt > 0)) cfFatal("twoFluid: dt must be positive");
    for (Junction& jn : junctions_) solveJunction(dt, jn);
    for (auto& p : pipes_) p->step(dt);
    t_ += dt;
}

void TwoFluidNetwork::solveJunction(double dt, Junction& jn)
{
    struct Branch
    {
        int pipe = 0;
        double s = 0;      // +1 junction at the outlet, -1 at the inlet
        double area = 0;
        bool free = false; // pressure-driven far end: velocity unknown
        double vM = 0;     // known velocity (rate far end / dead leg)
        double pFar = 0;   // far fixed pressure (free branches)
        int pos = -1;      // row/column in the Newton system (free branches)
    };
    std::vector<Branch> br;
    for (const EndRef& e : jn.ends)
    {
        Branch b;
        b.pipe = e.pipe;
        b.s = e.atOutlet ? 1.0 : -1.0;
        b.area = pipes_[e.pipe]->area();
        b.free = e.far.type == TwoFluidBoundary::Type::fixedPressure;
        if (b.free) b.pFar = e.far.pressure;
        else if (e.far.type == TwoFluidBoundary::Type::phaseRates)
            b.vM = e.far.vsl + e.far.vsg;
        // sealed far end: dead leg, vM = 0
        br.push_back(b);
    }
    const int m = static_cast<int>(br.size());
    int nf = 0;
    for (int k = 0; k < m; ++k)
        if (br[k].free) br[k].pos = nf++;

    // unknowns: the free-branch velocities, then the shared junction
    // pressure; the momentum rows are branch-diagonal (each free branch's
    // residual sees only its own velocity), the balance row couples them
    std::vector<double> x(nf + 1);
    for (int k = 0; k < m; ++k)
        if (br[k].free) x[br[k].pos] = pipes_[br[k].pipe]->mixtureVelocity();
    x[nf] = jn.p;

    auto vMOf = [&](const std::vector<double>& xv, int k)
    {
        return br[k].free ? xv[br[k].pos] : br[k].vM;
    };
    auto balance = [&](const std::vector<double>& xv)
    {
        double q = 0;
        for (int k = 0; k < m; ++k)
            q += br[k].s * br[k].area * vMOf(xv, k);
        return q;
    };
    auto dpDriveOf = [&](int k, double pJ)
    {
        // the driving difference toward the junction, signed like the
        // serial fixed-pressure form: positive drives flow into the
        // junction for s = +1 branches
        return br[k].s * (br[k].pFar - pJ);
    };

    std::vector<double> r(nf + 1), scale(nf + 1);
    auto residuals = [&](const std::vector<double>& xv)
    {
        for (int k = 0; k < m; ++k)
        {
            if (!br[k].free) continue;
            r[br[k].pos] = pipes_[br[k].pipe]->momentumResidual(
                dt, vMOf(xv, k), dpDriveOf(k, xv[nf]));
        }
        r[nf] = balance(xv);
    };
    for (int k = 0; k < m; ++k)
    {
        if (!br[k].free) continue;
        scale[br[k].pos] = pipes_[br[k].pipe]->momentumScale(dt)
                         + std::abs(br[k].pFar) + std::abs(x[nf]);
    }
    double qScale = 1e-12;
    for (int k = 0; k < m; ++k)
        qScale += br[k].area * (std::abs(vMOf(x, k)) + 1.0);
    scale[nf] = qScale;

    bool ok = false;
    for (int it = 0; it < 60; ++it)
    {
        residuals(x);
        ok = true;
        for (int q = 0; q <= nf; ++q)
            if (!(std::abs(r[q]) <= 1e-12 * scale[q])) { ok = false; break; }
        if (ok) break;

        std::vector<std::vector<double>> J(
            nf + 1, std::vector<double>(nf + 1, 0.0));
        for (int k = 0; k < m; ++k)
        {
            if (!br[k].free) continue;
            TwoFluidPipe& p = *pipes_[br[k].pipe];
            const double u = x[br[k].pos];
            const double h
                = std::max(1e-6 * std::max(std::abs(u), 1e-3), 1e-12);
            const double slope
                = (p.momentumResidual(dt, u + h, dpDriveOf(k, x[nf]))
                   - p.momentumResidual(dt, u - h, dpDriveOf(k, x[nf])))
                  / (2.0 * h);
            if (!(slope > 0))
                cfFatal("twoFluid: non-monotone branch momentum at junction '"
                        + jn.name + "'");
            J[br[k].pos][br[k].pos] = slope;
            J[br[k].pos][nf] = br[k].s; // dR/dpJ = s (R = ... - s(pFar-pJ))
        }
        for (int k = 0; k < m; ++k)
            if (br[k].free) J[nf][br[k].pos] = br[k].s * br[k].area;
        std::vector<double> dx = r;
        for (double& v : dx) v = -v;
        detail::solveDenseLinear(J, dx);
        for (int q = 0; q <= nf; ++q) x[q] += dx[q];
    }
    if (!ok)
        cfFatal("twoFluid: junction '" + jn.name
                + "' Newton did not converge");
    jn.p = x[nf];
    for (int k = 0; k < m; ++k)
        if (br[k].free) br[k].vM = x[br[k].pos];

    // every branch takes its velocity (free branches the converged one),
    // then the inflow streams mix: holdup/water/hydrate volume-flow
    // weighted, temperature heat-capacity-flow weighted, gas density
    // gas-mass-flux weighted (with the EOS this closes the junction gas
    // mass balance exactly)
    for (int k = 0; k < m; ++k)
        pipes_[br[k].pipe]->setExternalMomentum(br[k].vM);

    const TwoFluidFluid& f = cfg_.perPipe.fluid;
    double qL = 0, qG = 0, wNum = 0, hNum = 0, eNum = 0, eDen = 0, gNum = 0;
    for (int k = 0; k < m; ++k)
    {
        if (br[k].s * br[k].vM <= 0) continue; // fed or static: no inflow
        const TwoFluidPipe& p = *pipes_[br[k].pipe];
        const bool atInlet = !jn.ends[k].atOutlet;
        const double fl = p.endLiquidFlux(atInlet);
        const double fg = p.endGasVolumeFlux(atInlet);
        if (fl * br[k].s < 0 || fg * br[k].s < 0)
            cfFatal("twoFluid: countercurrent phase flux at junction '"
                    + jn.name + "' not supported");
        const double aL = br[k].area * std::abs(fl);
        const double aG = br[k].area * std::abs(fg);
        const double rhoG = p.endGasDensity(atInlet);
        qL += aL;
        qG += aG;
        wNum += aL * p.endWater(atInlet);
        hNum += aL * p.endHydrate(atInlet);
        eNum += (aL * f.rhoL * f.cL + aG * rhoG * f.cG)
              * p.endTemperature(atInlet);
        eDen += aL * f.rhoL * f.cL + aG * rhoG * f.cG;
        gNum += aG * rhoG;
    }
    const double qM = qL + qG;
    if (qM > 0)
    {
        jn.alphaJ = qL / qM;
        if (qL > 0)
        {
            jn.waterJ = wNum / qL;
            jn.hydrateJ = hNum / qL;
        }
        if (eDen > 0) jn.tJ = eNum / eDen;
        if (qG > 0) jn.rhoJ = gNum / qG;
    }

    // junction faces: feeding branches get a fixed-pressure ghost at the
    // mixed state (their face flux upwinds from the boundary cell);
    // supplied branches get the junction split imposed as phaseRates at
    // the pipe's own flux functions at the mixed holdup (F + G = vM keeps
    // the volume balance exact for every closure; the phase split is
    // exact for the homogeneous closure, closure-approximate with slip)
    for (int k = 0; k < m; ++k)
    {
        TwoFluidPipe& p = *pipes_[br[k].pipe];
        TwoFluidBoundary b;
        b.pressure = jn.p;
        b.alpha = jn.alphaJ;
        b.waterFraction = jn.waterJ;
        b.hydrateFraction = jn.hydrateJ;
        b.temperature = jn.tJ;
        b.rhoGInjected = jn.rhoJ;
        if (br[k].s * br[k].vM > 0)
        {
            b.type = TwoFluidBoundary::Type::fixedPressure;
        }
        else
        {
            b.type = TwoFluidBoundary::Type::phaseRates;
            b.vsl = p.liquidFluxAt(jn.alphaJ);
            b.vsg = p.gasVolumeFluxAt(jn.alphaJ);
        }
        if (jn.ends[k].atOutlet) p.setOutletBoundary(b);
        else p.setInletBoundary(b);
    }
}

std::vector<std::string> TwoFluidNetwork::pipeNames() const
{
    std::vector<std::string> names;
    for (const TwoFluidGeo& g : cfg_.pipes) names.push_back(g.name);
    return names;
}

const TwoFluidPipe& TwoFluidNetwork::pipe(const std::string& name) const
{
    return *pipes_.at(pipeIndex_.at(name));
}

std::vector<std::string> TwoFluidNetwork::junctionNames() const
{
    std::vector<std::string> names;
    for (const Junction& jn : junctions_) names.push_back(jn.name);
    return names;
}

double TwoFluidNetwork::junctionPressure(const std::string& name) const
{
    for (const Junction& jn : junctions_)
        if (jn.name == name) return jn.p;
    cfFatal("twoFluid: unknown junction '" + name + "'");
    return 0;
}

double TwoFluidNetwork::junctionNetRate(const std::string& name) const
{
    double net = 0;
    for (double q : junctionPipeRates(name)) net += q;
    return net;
}

std::vector<double> TwoFluidNetwork::junctionPipeRates(
    const std::string& name) const
{
    for (const Junction& jn : junctions_)
        if (jn.name == name)
        {
            std::vector<double> rates;
            for (const EndRef& e : jn.ends)
                rates.push_back((e.atOutlet ? 1.0 : -1.0)
                                * pipes_[e.pipe]->area()
                                * pipes_[e.pipe]->mixtureVelocity());
            return rates;
        }
    cfFatal("twoFluid: unknown junction '" + name + "'");
    return {};
}

} // namespace cf
