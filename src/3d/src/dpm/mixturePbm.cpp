#include <cmath>

#include "cf/common/Error.h"
#include "cf/fvm/MixturePbm.h"

namespace cf
{
namespace
{
constexpr scalar kPi = 3.14159265358979323846;

scalar particleVolume(const MixturePbmConfig& cfg)
{
    return kPi * cfg.d0 * cfg.d0 * cfg.d0 / 6.0;
}
} // namespace

void MixturePbmConfig::validate() const
{
    if (!enabled) return; // the app reports a missing declaration itself

    if (!std::isfinite(rhoL) || !(rhoL > 0))
        cfFatal("mixturePbm: rhoL must be finite and positive");
    if (!std::isfinite(rhoG) || !(rhoG > 0))
        cfFatal("mixturePbm: rhoG must be finite and positive");
    if (!(rhoL > rhoG))
        cfFatal("mixturePbm: rhoL must exceed rhoG (invalid density ordering)");

    if (!std::isfinite(J) || J < 0)
        cfFatal("mixturePbm: J (nucleation rate) must be finite and non-negative");
    if (!std::isfinite(breakage) || breakage < 0)
        cfFatal("mixturePbm: breakage frequency must be finite and non-negative");

    if (!std::isfinite(alpha0) || !(alpha0 > 0))
        cfFatal("mixturePbm: alpha0 must be positive");
    if (!(alpha0 < 1))
        cfFatal("mixturePbm: alpha0 too large (must be < 1)");

    if (!std::isfinite(d0) || !(d0 > 0))
        cfFatal("mixturePbm: d0 must be finite and positive");
    if (!std::isfinite(slip))
        cfFatal("mixturePbm: slip must be finite");
    if (!std::isfinite(j) || !(j > 0))
        cfFatal("mixturePbm: j must be finite and positive");
    if (!(j + slip > 0))
        cfFatal("mixturePbm: j+slip must be positive (invalid slip)");
    if (!(nFrag >= 2))
        cfFatal("mixturePbm: nFrag must be at least 2");

    if (!std::isfinite(dt) || !(dt > 0))
        cfFatal("mixturePbm: dt must be finite and positive");
    if (steps < 1)
        cfFatal("mixturePbm: steps must be at least 1");
    if (reportEvery < 1)
        cfFatal("mixturePbm: reportEvery must be at least 1");
}

MixturePbmConfig readMixturePbmConfig(const Dictionary& caseDict)
{
    MixturePbmConfig cfg;
    if (!caseDict.has("mixturePbm")) return cfg;
    const Dictionary& d = caseDict.getDict("mixturePbm");
    cfg.enabled = true;
    cfg.alpha0 = d.getScalarOr("alpha0", 0.0);
    cfg.d0 = d.getScalarOr("d0", 0.0);
    cfg.rhoL = d.getScalarOr("rhoL", 0.0);
    cfg.rhoG = d.getScalarOr("rhoG", 0.0);
    cfg.slip = d.getScalarOr("slip", 0.0);
    cfg.j = d.getScalarOr("j", 1.0);
    cfg.J = d.getScalarOr("J", 0.0);
    cfg.breakage = d.getScalarOr("breakage", 0.0);
    cfg.nFrag = d.getScalarOr("nFrag", 2.0);
    cfg.dt = d.getScalarOr("dt", 0.0);
    cfg.steps = static_cast<int>(d.getScalarOr("steps", 0.0));
    cfg.reportEvery = static_cast<int>(d.getScalarOr("reportEvery", 1.0));
    cfg.validate();
    return cfg;
}

MixturePbmProps mixtureProperties(const MixturePbmConfig& cfg)
{
    cfg.validate();
    MixturePbmProps p;
    p.vPart = particleVolume(cfg);
    p.k = cfg.j / (cfg.j + cfg.slip);
    p.m1_0 = cfg.alpha0;
    p.m0_0 = cfg.alpha0 / p.vPart;
    p.m2_0 = p.m0_0 * p.vPart * p.vPart;
    p.m3_0 = p.m0_0 * p.vPart * p.vPart * p.vPart;
    return p;
}

void mixturePbmSource(const MixturePbmConfig& cfg, const scalar m[4], scalar dm[4])
{
    const scalar v = particleVolume(cfg);
    const scalar f2 = 1.0 - 1.0 / cfg.nFrag;
    const scalar f3 = 1.0 - 1.0 / (cfg.nFrag * cfg.nFrag);

    dm[0] = cfg.J + (cfg.nFrag - 1.0) * cfg.breakage * m[0];
    dm[1] = cfg.J * v;
    dm[2] = cfg.J * v * v - f2 * cfg.breakage * v * v * m[0];
    dm[3] = cfg.J * v * v * v - f3 * cfg.breakage * v * v * v * m[0];
}

MixturePbmResult runMixturePbm(const MixturePbmConfig& cfg)
{
    cfg.validate();
    MixturePbmResult out;
    out.props = mixtureProperties(cfg);
    const MixturePbmProps& p = out.props;

    scalar m[4] = {p.m0_0, p.m1_0, p.m2_0, p.m3_0};

    for (int step = 0; step <= cfg.steps; ++step)
    {
        if (step % cfg.reportEvery == 0 || step == cfg.steps)
        {
            MixturePbmStep s;
            s.step = step;
            s.t = step * cfg.dt;
            for (int k = 0; k < 4; ++k) s.m[k] = m[k];
            s.alpha = m[1] * p.k;
            if (!(s.alpha > 0) || !(s.alpha < 1))
                cfFatal("mixturePbm: alpha left (0,1) at step " + std::to_string(step));
            s.uG = cfg.j + (1.0 - s.alpha) * cfg.slip;
            s.uL = cfg.j - s.alpha * cfg.slip;
            out.history.push_back(s);
        }
        if (step == cfg.steps) break;

        scalar dm[4];
        mixturePbmSource(cfg, m, dm);
        scalar next[4];
        for (int k = 0; k < 4; ++k) next[k] = m[k] + cfg.dt * dm[k];
        for (int k = 0; k < 4; ++k)
            if (!std::isfinite(next[k]) || !(next[k] > 0))
                cfFatal("mixturePbm: moment M" + std::to_string(k)
                        + " became non-physical at step " + std::to_string(step + 1));
        for (int k = 0; k < 4; ++k) m[k] = next[k];
    }
    out.valid = true;
    return out;
}

} // namespace cf
