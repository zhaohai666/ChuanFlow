#include <cmath>
#include <limits>

#include "cf/common/Error.h"
#include "cf/fvm/Combustion.h"

namespace cf
{

scalar arrheniusRate(scalar preExp, scalar activationE, scalar gasConstant, scalar T)
{
    if (!(T > 0) || !(gasConstant > 0))
        cfFatal("combustion: Arrhenius needs positive T and gas constant");
    if (activationE == 0) return preExp;
    return preExp * std::exp(-activationE / (gasConstant * T));
}

scalar CombustionConfig::rateConstant() const
{
    return arrheniusRate(preExp, activationE, gasConstant, temperature);
}

void CombustionConfig::validate() const
{
    if (!std::isfinite(diffusivity) || !(diffusivity > 0))
        cfFatal("combustion: diffusivity must be finite and positive (m2/s)");
    if (!std::isfinite(preExp) || !(preExp > 0))
        cfFatal("combustion: preExp A must be finite and positive (1/s)");
    if (!std::isfinite(activationE) || activationE < 0)
        cfFatal("combustion: activationE Ea cannot be negative (J/mol)");
    if (!std::isfinite(gasConstant) || !(gasConstant > 0))
        cfFatal("combustion: gasConstant R must be finite and positive");
    if (!std::isfinite(temperature) || !(temperature > 0))
        cfFatal("combustion: temperature must be a finite positive field (K)");
    if (!std::isfinite(stoich) || !(stoich > 0))
        cfFatal("combustion: stoich s must be finite and positive");
    if (fixedOxidiser >= 0 && (!std::isfinite(fixedOxidiser) || fixedOxidiser > 1))
        cfFatal("combustion: fixedOxidiser must be in [0,1]");
    if (!std::isfinite(fuelInlet) || fuelInlet < 0 || fuelInlet > 1)
        cfFatal("combustion: fuelInlet must be in [0,1]");
    if (!std::isfinite(oxInlet) || oxInlet < 0 || oxInlet > 1)
        cfFatal("combustion: oxInlet must be in [0,1]");
}

CombustionConfig readCombustionConfig(const Dictionary& caseDict)
{
    CombustionConfig cfg;
    if (!caseDict.has("combustion")) return cfg;
    const Dictionary& cd = caseDict.getDict("combustion");
    cfg.enabled = true;
    cfg.diffusivity = cd.getScalarOr("diffusivity", 0.0);
    cfg.preExp = cd.getScalarOr("preExp", 0.0);
    cfg.activationE = cd.getScalarOr("activationE", 0.0);
    cfg.gasConstant = cd.getScalarOr("gasConstant", 8.314462618);
    cfg.temperature = cd.getScalarOr("temperature", 0.0);
    cfg.stoich = cd.getScalarOr("stoich", 1.0);
    cfg.fixedOxidiser = cd.getScalarOr("fixedOxidiser", -1.0);
    cfg.fuelInlet = cd.getScalarOr("fuelInlet", 0.0);
    cfg.oxInlet = cd.getScalarOr("oxInlet", 0.0);
    cfg.validate();
    return cfg;
}

namespace
{

/// Assemble and solve one linear species equation:
///   div(phi Y) - div(D grad Y) + Sp Y = rhs_explicit   (all on lhs = 0 form).
/// `Sp` is per-cell reaction coefficient (already includes volume scaling),
/// `rhs` is the explicit part with volume. Returns the iteration count.
int solveSpecies(fvMatrix& eqn, const polyMesh& mesh, volScalarField& Y,
                 const surfaceScalarField& phi, const Schemes& schemes,
                 scalar D, const Field<scalar>& Sp,
                 const Field<scalar>& rhsExplicit)
{
    eqn.schemes() = schemes;
    eqn.schemes().relaxation = 1.0;
    fvm::div(eqn, phi.values(), schemes.div);
    fvm::laplacian(eqn, D, schemes.laplacian, -1);
    fvm::sp(eqn, Sp, +1);
    const label nC = mesh.nCells();
    for (label c = 0; c < nC; ++c)
        eqn.addRhs(c, 0, rhsExplicit[static_cast<std::size_t>(c)]);
    const SolverReport rep = eqn.solve(fvMatrix::Controls{1e-12, 20000, 2});
    Y.internalField() = eqn.solution()[0];
    if (!rep.converged)
        cfFatal("combustion: species solve did not converge ("
                + std::to_string(rep.iterations) + " iters, residual "
                + std::to_string(rep.finalResidual) + ")");
    return rep.iterations;
}

} // namespace

CombustionResult solveCombustion(const polyMesh& mesh, const surfaceScalarField& phi,
                                 volScalarField& Yf, volScalarField& Yo,
                                 const CombustionConfig& cfg, const Schemes& schemes,
                                 bool fixedOxidiser)
{
    cfg.validate();
    const label nC = mesh.nCells();
    const auto& vol = mesh.cellVolumes();
    const scalar D = cfg.diffusivity;
    const scalar k = cfg.rateConstant();
    const scalar s = cfg.stoich;

    CombustionResult out;
    out.Yf.assign(nC, 0.0);
    out.Yo.assign(nC, 0.0);
    out.reactionRate.assign(nC, 0.0);

    // Yf / Yo arrive with their patch BCs set by the app.
    const Field<scalar> zero(static_cast<std::size_t>(nC), 0.0);
    Field<scalar> spF(static_cast<std::size_t>(nC), 0.0);
    Field<scalar> rhsF(static_cast<std::size_t>(nC), 0.0);
    Field<scalar> spO(static_cast<std::size_t>(nC), 0.0);
    Field<scalar> rhsO(static_cast<std::size_t>(nC), 0.0);

    if (fixedOxidiser)
    {
        const scalar yo = cfg.fixedOxidiser;
        for (label c = 0; c < nC; ++c)
            spF[static_cast<std::size_t>(c)] = k * yo * vol[static_cast<std::size_t>(c)];
        fvMatrix eqnF(mesh, viewOf(Yf), "YfEqn");
        out.fuelIters = solveSpecies(eqnF, mesh, Yf, phi, schemes, D, spF, zero);
        out.Yf = Yf.internalField();
        for (label c = 0; c < nC; ++c)
        {
            const std::size_t u = static_cast<std::size_t>(c);
            out.Yo[u] = yo;
            out.reactionRate[u] = k * out.Yf[u] * yo;
        }
    }
    else
    {
        // Picard: freeze Yo in the fuel equation and Yf in the oxidiser
        // equation, solving each linear system in turn until both settle.
        Field<scalar> yf(static_cast<std::size_t>(nC), 0.5 * cfg.fuelInlet);
        Field<scalar> yo(static_cast<std::size_t>(nC), 0.5 * cfg.oxInlet);
        scalar lastResidual = std::numeric_limits<scalar>::infinity();
        const int maxOuter = 100;
        int used = 0;
        for (int outer = 1; outer <= maxOuter; ++outer)
        {
            for (label c = 0; c < nC; ++c)
            {
                const std::size_t u = static_cast<std::size_t>(c);
                spF[u] = k * yo[u] * vol[u];
                spO[u] = s * k * yf[u] * vol[u];
            }
            fvMatrix eqnF(mesh, viewOf(Yf), "YfEqn");
            fvMatrix eqnO(mesh, viewOf(Yo), "YoEqn");
            out.fuelIters = solveSpecies(eqnF, mesh, Yf, phi, schemes, D, spF, zero);
            out.oxIters = solveSpecies(eqnO, mesh, Yo, phi, schemes, D, spO, zero);
            const Field<scalar>& nf = Yf.internalField();
            const Field<scalar>& no = Yo.internalField();
            scalar delta = 0, scale = 0;
            for (std::size_t u = 0; u < nf.size(); ++u)
            {
                delta = std::max(delta, std::abs(nf[u] - yf[u]));
                delta = std::max(delta, std::abs(no[u] - yo[u]));
                scale = std::max(scale, std::abs(nf[u]));
                scale = std::max(scale, std::abs(no[u]));
            }
            lastResidual = delta / std::max(scale, 1e-30);
            yf = nf;
            yo = no;
            used = outer;
            if (lastResidual <= 1e-9) break;
        }
        if (used >= maxOuter)
            cfFatal("combustion: coupled species Picard did not converge in "
                    + std::to_string(maxOuter) + " sweeps (residual "
                    + std::to_string(lastResidual) + ")");
        out.Yf = yf;
        out.Yo = yo;
    }

    // Conservation diagnostic: integrate the reaction sink over the domain.
    out.fuelConsumed = 0;
    out.oxConsumed = 0;
    for (label c = 0; c < nC; ++c)
    {
        const std::size_t u = static_cast<std::size_t>(c);
        out.reactionRate[u] = k * out.Yf[u] * out.Yo[u];
        out.fuelConsumed += vol[u] * out.reactionRate[u];
        out.oxConsumed += s * vol[u] * out.reactionRate[u];
    }
    out.solved = true;
    return out;
}

} // namespace cf
