#pragma once

#include <vector>

#include "cf/common/Dictionary.h"
#include "cf/common/Labels.h"

namespace cf
{

/// Batch D / F7: opt-in algebraic mixture slip (drift-flux) plus a
/// monodisperse three-moment population balance. This is a compact, spatially
/// uniform moment kernel (the source term a DPM/species framework would call),
/// not a multiphase CFD solver. Nothing is registered unless a case declares a
/// mixturePbm{} block, so existing apps and their cases are left untouched.
///
/// Monodisperse closure: every particle (nucleus and parent) has the single
/// representative volume v = pi*d0^3/6. Nucleation creates particles of that
/// volume at rate J [1/(m3 s)]; breakage at frequency b [1/s] splits a parent
/// into nFrag equal fragments (volume conserved). The four transported moments
/// M_k = integral n(v) v^k dv (k = 0..3) obey, with f2 = 1 - 1/nFrag and
/// f3 = 1 - 1/nFrag^2:
///
///   dM0/dt = J + (nFrag-1) b M0          number (nucleation + binary gain)
///   dM1/dt = J v                         volume (breakage conserves volume)
///   dM2/dt = J v^2 - f2 b v^2 M0         surface-like moment
///   dM3/dt = J v^3 - f3 b v^3 M0         volume-squared moment
///
/// Because the parent size v is fixed by the closure (it is NOT M1/M0), the
/// source has no moment ratios and the ODEs have a closed form: M1 is linear,
/// M0 is exponential (or linear when b=0), and M2/M3 are driven by the
/// integral of M0. The checker replays the SAME explicit-Euler rule from the
/// case constants alone and also grades the closed forms.
///
/// Drift-flux: with mixture superficial velocity j and slip (drift) velocity
/// s = uG-uL, the holdup is redistributed by k = j/(j+s), so the in-situ gas
/// fraction is alpha(t) = M1(t) k. The phase velocities are uG = j+(1-alpha)s
/// and uL = j-alpha*s, hence uG-uL = s identically. Changing slip therefore
/// changes the reported phase distribution without touching the PBM source.
struct MixturePbmConfig
{
    bool enabled = false;
    scalar alpha0 = 0;   // inlet gas volume fraction (0 < alpha0 < 1)
    scalar d0 = 0;       // representative particle diameter [m]
    scalar rhoL = 0;     // liquid density [kg/m3], finite positive
    scalar rhoG = 0;     // gas density [kg/m3], finite positive, < rhoL
    scalar slip = 0;     // drift velocity uG-uL [m/s]
    scalar j = 1;        // mixture superficial velocity [m/s]
    scalar J = 0;        // nucleation rate [1/(m3 s)], non-negative
    scalar breakage = 0; // breakage frequency [1/s], non-negative
    scalar nFrag = 2;    // daughter fragments per breakage event (>= 2)
    scalar dt = 0;       // s
    int steps = 0;
    int reportEvery = 1;

    void validate() const;
};

MixturePbmConfig readMixturePbmConfig(const Dictionary& caseDict);

/// Closure-derived constants, re-derivable by a checker from the case text.
struct MixturePbmProps
{
    scalar vPart = 0;    // pi*d0^3/6, representative particle volume
    scalar k = 1;        // drift-flux holdup factor j/(j+slip)
    scalar m0_0 = 0;     // alpha0 / vPart
    scalar m1_0 = 0;     // alpha0
    scalar m2_0 = 0;     // alpha0 * vPart
    scalar m3_0 = 0;     // alpha0 * vPart^2
};

MixturePbmProps mixtureProperties(const MixturePbmConfig& cfg);

/// Documented monodisperse moment source, exposed so a checker can replay the
/// identical law without trusting the solver's trajectory.
void mixturePbmSource(const MixturePbmConfig& cfg, const scalar m[4], scalar dm[4]);

struct MixturePbmStep
{
    int step = 0;
    scalar t = 0;
    scalar m[4] = {0, 0, 0, 0};
    scalar alpha = 0;    // in-situ gas holdup = M1 * k
    scalar uG = 0;
    scalar uL = 0;
};

struct MixturePbmResult
{
    std::vector<MixturePbmStep> history;
    MixturePbmProps props;
    bool valid = false;
};

MixturePbmResult runMixturePbm(const MixturePbmConfig& cfg);

} // namespace cf
