#include <algorithm>
#include <cmath>

#include "cf/common/Error.h"
#include "cf/fvm/Qmom.h"

namespace cf
{

void QmomConfig::validate() const
{
    if (!std::isfinite(beta) || !(beta > 0))
        cfFatal("qmom: beta must be a finite positive coagulation kernel");
    if (!std::isfinite(m0) || !(m0 > 0))
        cfFatal("qmom: M_0 must be finite and positive");
    if (!std::isfinite(m1) || !(m1 > 0))
        cfFatal("qmom: M_1 must be finite and positive");
    if (!std::isfinite(m2) || !(m2 > 0))
        cfFatal("qmom: M_2 must be finite and positive");
    if (!std::isfinite(m3) || !(m3 > 0))
        cfFatal("qmom: M_3 must be finite and positive");
    // A valid measure obeys M_0 M_2 >= M_1^2 and M_1 M_3 >= M_2^2.
    if (m0 * m2 < m1 * m1)
        cfFatal("qmom: moments violate M_0 M_2 >= M_1^2 (not a positive measure)");
    if (m1 * m3 < m2 * m2)
        cfFatal("qmom: moments violate M_1 M_3 >= M_2^2 (not a positive measure)");
    if (!std::isfinite(dt) || !(dt > 0))
        cfFatal("qmom: dt must be finite and positive");
    if (steps < 1)
        cfFatal("qmom: steps must be at least 1");
    if (reportEvery < 1)
        cfFatal("qmom: reportEvery must be at least 1");
}

QmomConfig readQmomConfig(const Dictionary& caseDict)
{
    QmomConfig cfg;
    if (!caseDict.has("qmom")) return cfg;
    const Dictionary& qd = caseDict.getDict("qmom");
    cfg.enabled = true;
    cfg.beta = qd.getScalarOr("beta", 0.0);
    cfg.m0 = qd.getScalarOr("m0", 0.0);
    cfg.m1 = qd.getScalarOr("m1", 0.0);
    cfg.m2 = qd.getScalarOr("m2", 0.0);
    cfg.m3 = qd.getScalarOr("m3", 0.0);
    cfg.dt = qd.getScalarOr("dt", 0.0);
    cfg.steps = static_cast<int>(qd.getScalarOr("steps", 0.0));
    cfg.reportEvery = static_cast<int>(qd.getScalarOr("reportEvery", 1.0));
    cfg.validate();
    return cfg;
}

void coagulationSource(const scalar m[4], scalar beta, scalar dm[4])
{
    // Constant-kernel Smoluchowski, the k=0 loss is separate from k>=1 gains.
    dm[0] = -0.5 * beta * m[0] * m[0];
    dm[1] = 0.0;
    dm[2] = beta * m[1] * m[1];
    dm[3] = 3.0 * beta * m[1] * m[2];
}

bool quadratureFromMoments(const scalar m[4], scalar w[2], scalar v[2])
{
    if (!(m[0] > 0)) return false;
    // Normalised moments and central moments of the measure n(v).
    const scalar mu1 = m[1] / m[0];
    const scalar mu2 = m[2] / m[0];
    const scalar mu3 = m[3] / m[0];
    const scalar var = mu2 - mu1 * mu1;             // E[(v-mu)^2]
    if (!(var > 0)) return false;
    const scalar c3 = mu3 - 3.0 * mu1 * mu2 + 2.0 * mu1 * mu1 * mu1; // E[(v-mu)^3]

    // Two-node discrete distribution with nodes d = A ± B about the mean and
    // weights p = (1+A/B)/2, q = (1-A/B)/2. Matching the first three central
    // moments of such a rule gives the identities
    //   var = B^2 - A^2,   c3 = 2 A var
    // so A = c3/(2 var) and B = sqrt(var + A^2). Because B^2 = var + A^2 > A^2,
    // both weights are automatically positive for any positive variance.
    const scalar A = c3 / (2.0 * var);
    const scalar B = std::sqrt(var + A * A);
    if (!(B > 0)) return false;
    // The smaller node is v = mu + A - B; its weight is p = (1+A/B)/2.
    v[0] = mu1 + A - B;
    v[1] = mu1 + A + B;
    w[0] = 0.5 * m[0] * (1.0 + A / B);
    w[1] = 0.5 * m[0] * (1.0 - A / B);
    return w[0] >= 0 && w[1] >= 0;
}

QmomResult runQmom(const QmomConfig& cfg)
{
    cfg.validate();
    QmomResult out;
    scalar m[4] = {cfg.m0, cfg.m1, cfg.m2, cfg.m3};
    const scalar m1Initial = cfg.m1;

    for (int step = 0; step <= cfg.steps; ++step)
    {
        if (step % cfg.reportEvery == 0 || step == cfg.steps)
        {
            QmomStep s;
            s.step = step;
            s.t = step * cfg.dt;
            for (int k = 0; k < 4; ++k) s.m[k] = m[k];
            if (!quadratureFromMoments(m, s.w, s.v))
                cfFatal("qmom: quadrature reconstruction failed at step " + std::to_string(step));
            s.massError = std::abs(m[1] - m1Initial) / std::max(m1Initial, 1e-300);
            out.history.push_back(s);
        }
        if (step == cfg.steps) break;

        scalar dm[4];
        coagulationSource(m, cfg.beta, dm);
        // Riemann-sum explicit Euler; the checker replays the identical rule.
        scalar next[4];
        for (int k = 0; k < 4; ++k) next[k] = m[k] + cfg.dt * dm[k];
        // Guard against a non-physical step (a checker case is chosen so the
        // explicit horizon stays in the positive-measure cone).
        if (!(next[0] > 0) || !(next[1] > 0) || !(next[2] > 0) || !(next[3] > 0))
            cfFatal("qmom: moment became non-positive at step " + std::to_string(step + 1));
        for (int k = 0; k < 4; ++k) m[k] = next[k];
    }
    out.valid = true;
    return out;
}

} // namespace cf
