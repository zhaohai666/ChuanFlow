#include <algorithm>
#include <cmath>
#include <limits>

#include "cf/common/Error.h"
#include "cf/fvm/CombustionDetailed.h"

namespace cf
{

scalar detailedArrhenius(scalar A, scalar Ea, scalar R, scalar T)
{
    if (!(T > 0) || !(R > 0))
        cfFatal("detailedChemistry: Arrhenius needs positive T and gas constant");
    if (Ea == 0) return A;
    return A * std::exp(-Ea / (R * T));
}

void DetailedChemistryConfig::validate() const
{
    const auto nonNeg = [](scalar v) { return std::isfinite(v) && v >= 0; };
    if (!nonNeg(A1)) cfFatal("detailedChemistry: A1 branching pre-exponential must be finite and non-negative");
    if (!nonNeg(E1)) cfFatal("detailedChemistry: E1 branching activation energy must be finite and non-negative");
    if (!nonNeg(A2)) cfFatal("detailedChemistry: A2 termination pre-exponential must be finite and non-negative");
    if (!nonNeg(E2)) cfFatal("detailedChemistry: E2 termination activation energy must be finite and non-negative");
    if (!nonNeg(q)) cfFatal("detailedChemistry: q heat release must be finite and non-negative");
    if (!std::isfinite(T0) || !(T0 > 0))
        cfFatal("detailedChemistry: T0 initial temperature must be finite and positive (K)");
    if (!std::isfinite(gasConstant) || !(gasConstant > 0))
        cfFatal("detailedChemistry: gasConstant R must be finite and positive");
    if (!std::isfinite(Yf0) || Yf0 < 0 || Yf0 > 1)
        cfFatal("detailedChemistry: Yf0 initial fuel fraction must be in [0,1]");
    if (!std::isfinite(Yr0) || Yr0 < 0 || Yr0 > 1)
        cfFatal("detailedChemistry: Yr0 initial radical fraction must be in [0,1]");
    if (!nonNeg(cooling)) cfFatal("detailedChemistry: cooling coefficient must be finite and non-negative");
    if (!std::isfinite(Tenv) || !(Tenv > 0))
        cfFatal("detailedChemistry: Tenv ambient temperature must be finite and positive (K)");
    if (!std::isfinite(dt) || !(dt > 0))
        cfFatal("detailedChemistry: dt must be finite and positive (s)");
    if (nSteps <= 0) cfFatal("detailedChemistry: nSteps must be positive");
    if (nSub <= 0) cfFatal("detailedChemistry: nSub must be positive");
    if (!std::isfinite(ignitionRise) || !(ignitionRise > 0) || !(ignitionRise < 1))
        cfFatal("detailedChemistry: ignitionRise must be in (0,1)");
}

DetailedChemistryConfig readDetailedChemistryConfig(const Dictionary& caseDict)
{
    DetailedChemistryConfig cfg;
    if (!caseDict.has("detailedChemistry")) return cfg;
    const Dictionary& cd = caseDict.getDict("detailedChemistry");
    cfg.enabled = true;
    cfg.A1 = cd.getScalarOr("A1", 0.0);
    cfg.E1 = cd.getScalarOr("E1", 0.0);
    cfg.A2 = cd.getScalarOr("A2", 0.0);
    cfg.E2 = cd.getScalarOr("E2", 0.0);
    cfg.q = cd.getScalarOr("q", 0.0);
    cfg.T0 = cd.getScalarOr("T0", 0.0);
    cfg.gasConstant = cd.getScalarOr("gasConstant", 8.314462618);
    cfg.Yf0 = cd.getScalarOr("Yf0", 1.0);
    cfg.Yr0 = cd.getScalarOr("Yr0", 1e-6);
    cfg.cooling = cd.getScalarOr("cooling", 0.0);
    cfg.Tenv = cd.getScalarOr("Tenv", cfg.T0);
    cfg.dt = cd.getScalarOr("dt", 0.01);
    cfg.nSteps = cd.getInt("nSteps");
    cfg.nSub = cd.getInt("nSub");
    cfg.ignitionRise = cd.getScalarOr("ignitionRise", 0.3);
    cfg.validate();
    return cfg;
}

namespace
{

struct Deriv
{
    scalar dT, dYf, dYr;
};

/// Right-hand side of the two-step reactor. `beta` is a cooling override so
/// the same routine can bisect for the critical heat-loss coefficient.
Deriv rhs(scalar T, scalar Yf, scalar Yr, const DetailedChemistryConfig& c, scalar beta)
{
    const scalar a = detailedArrhenius(c.A1, c.E1, c.gasConstant, T);
    const scalar b = detailedArrhenius(c.A2, c.E2, c.gasConstant, T);
    const scalar r1 = a * Yf * Yr; // branching  F + R -> P + 2R
    const scalar r2 = b * Yr;      // termination R -> P
    return {c.q * r1 - beta * (T - c.Tenv), -r1, r1 - r2};
}

void rk4Step(scalar& T, scalar& Yf, scalar& Yr, scalar h,
             const DetailedChemistryConfig& c, scalar beta)
{
    const Deriv k1 = rhs(T, Yf, Yr, c, beta);
    const Deriv k2 = rhs(T + 0.5 * h * k1.dT, Yf + 0.5 * h * k1.dYf,
                         Yr + 0.5 * h * k1.dYr, c, beta);
    const Deriv k3 = rhs(T + 0.5 * h * k2.dT, Yf + 0.5 * h * k2.dYf,
                         Yr + 0.5 * h * k2.dYr, c, beta);
    const Deriv k4 = rhs(T + h * k3.dT, Yf + h * k3.dYf,
                         Yr + h * k3.dYr, c, beta);
    const scalar s = h / 6;
    T += s * (k1.dT + 2 * k2.dT + 2 * k3.dT + k4.dT);
    Yf += s * (k1.dYf + 2 * k2.dYf + 2 * k3.dYf + k4.dYf);
    Yr += s * (k1.dYr + 2 * k2.dYr + 2 * k3.dYr + k4.dYr);
}

struct RunOut
{
    std::vector<DetailedState> hist;
    scalar peakT = 0;
    scalar tIgn = -1;
};

RunOut runReactor(const DetailedChemistryConfig& c, scalar beta, bool collect)
{
    RunOut o;
    o.peakT = c.T0;
    scalar T = c.T0, Yf = c.Yf0, Yr = c.Yr0;
    const scalar h = c.dt / static_cast<scalar>(c.nSub);
    const scalar trigger = c.T0 + c.ignitionRise * c.q * c.Yf0;
    if (collect) o.hist.push_back({0.0, T, Yf, Yr});
    for (int i = 0; i < c.nSteps; ++i)
    {
        for (int j = 0; j < c.nSub; ++j)
            rk4Step(T, Yf, Yr, h, c, beta);
        const scalar t = static_cast<scalar>(i + 1) * c.dt;
        if (T > o.peakT) o.peakT = T;
        if (o.tIgn < 0 && T >= trigger) o.tIgn = t;
        if (collect) o.hist.push_back({t, T, Yf, Yr});
    }
    return o;
}

bool ignitesAt(const DetailedChemistryConfig& c, scalar beta)
{
    const scalar trigger = c.T0 + c.ignitionRise * c.q * c.Yf0;
    return runReactor(c, beta, false).peakT >= trigger;
}

} // namespace

std::vector<DetailedState> integrateDetailed(const DetailedChemistryConfig& cfg)
{
    cfg.validate();
    return runReactor(cfg, cfg.cooling, true).hist;
}

DetailedDiagnostics diagnoseDetailed(const DetailedChemistryConfig& cfg)
{
    cfg.validate();
    DetailedDiagnostics d;
    const scalar a0 = detailedArrhenius(cfg.A1, cfg.E1, cfg.gasConstant, cfg.T0);
    const scalar b0 = detailedArrhenius(cfg.A2, cfg.E2, cfg.gasConstant, cfg.T0);
    d.branchDa = (b0 > 0) ? a0 * cfg.Yf0 / b0
                          : std::numeric_limits<scalar>::infinity();

    // Extinction Damkoehler threshold: the smallest heat-loss coefficient that
    // suppresses thermal runaway, found by bisection on a full reactor run.
    // betaCrit = 0 means even the adiabatic reactor never ignites (chemistry
    // is termination-dominated).
    if (!ignitesAt(cfg, 0))
    {
        d.betaCrit = 0;
    }
    else
    {
        scalar lo = 0, hi = 1;
        while (ignitesAt(cfg, hi) && hi < 1e6) hi *= 2;
        for (int it = 0; it < 36; ++it)
        {
            const scalar mid = 0.5 * (lo + hi);
            if (ignitesAt(cfg, mid)) lo = mid; else hi = mid;
        }
        d.betaCrit = 0.5 * (lo + hi);
    }
    const scalar loss = b0 + d.betaCrit;
    d.DaCrit = (loss > 0) ? a0 * cfg.Yf0 / loss
                          : std::numeric_limits<scalar>::infinity();

    const RunOut mainRun = runReactor(cfg, cfg.cooling, false);
    d.ignitionDelay = mainRun.tIgn;
    d.peakT = mainRun.peakT;
    d.ignited = d.ignitionDelay >= 0;
    return d;
}

} // namespace cf
