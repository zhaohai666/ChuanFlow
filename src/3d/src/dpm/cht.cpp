#include <algorithm>
#include <cmath>
#include <limits>

#include "cf/common/Error.h"
#include "cf/fvm/Cht.h"

namespace cf
{

void ChtConfig::validate() const
{
    if (!std::isfinite(solid.k) || !(solid.k > 0))
        cfFatal("cht: solid k must be finite and positive (W/m K)");
    if (!std::isfinite(fluid.k) || !(fluid.k > 0))
        cfFatal("cht: fluid k must be finite and positive (W/m K)");
    if (!std::isfinite(solid.outerT) || !std::isfinite(fluid.outerT))
        cfFatal("cht: outer temperatures must be finite");
    if (!(solid.outerT != fluid.outerT))
        cfFatal("cht: the two outer temperatures must differ to drive a flux");
    if (!std::isfinite(interfaceTol) || !(interfaceTol > 0))
        cfFatal("cht: interfaceTol must be finite and positive");
    if (maxIterations < 1)
        cfFatal("cht: maxIterations must be at least 1");
    if (!std::isfinite(relax) || !(relax > 0) || relax > 1)
        cfFatal("cht: relax must be in (0, 1]");
}

ChtConfig readChtConfig(const Dictionary& caseDict)
{
    ChtConfig cfg;
    if (!caseDict.has("cht")) return cfg;
    const Dictionary& cd = caseDict.getDict("cht");
    cfg.enabled = true;
    {
        const Dictionary& sd = cd.getDict("solid");
        cfg.solid.k = sd.getScalarOr("k", 0.0);
        cfg.solid.outerT = sd.getScalarOr("outerT", 0.0);
        cfg.solid.initialT = sd.getScalarOr("initialT", cfg.solid.outerT);
    }
    {
        const Dictionary& fd = cd.getDict("fluid");
        cfg.fluid.k = fd.getScalarOr("k", 0.0);
        cfg.fluid.outerT = fd.getScalarOr("outerT", 0.0);
        cfg.fluid.initialT = fd.getScalarOr("initialT", cfg.fluid.outerT);
    }
    cfg.interfaceTol = cd.getScalarOr("interfaceTol", 1e-9);
    cfg.maxIterations = static_cast<int>(cd.getScalarOr("maxIterations", 50.0));
    cfg.relax = cd.getScalarOr("relax", 1.0);
    cfg.validate();
    return cfg;
}

RegionResult solveRegion(const polyMesh& mesh, const Schemes& schemes,
                         scalar k, scalar outerT, scalar interfaceT,
                         const std::string& outerPatch,
                         const std::string& interfacePatch)
{
    const label nC = mesh.nCells();
    RegionResult out;
    out.T.assign(nC, 0.0);

    volScalarField T(mesh, "T", 0.5 * (outerT + interfaceT));
    for (label pid = 0; pid < mesh.nPatches(); ++pid)
    {
        const std::string& pname = mesh.patches()[static_cast<std::size_t>(pid)].name;
        if (pname == outerPatch)
            T.setPatchValue(pid, "fixedValue", outerT);
        else if (pname == interfacePatch)
            T.setPatchValue(pid, "fixedValue", interfaceT);
        else
            T.setPatchValue(pid, "zeroGradient", 0.0);
    }

    fvMatrix eqn(mesh, viewOf(T), "TEqn");
    eqn.schemes() = schemes;
    eqn.schemes().relaxation = 1.0;
    fvm::laplacian(eqn, k, schemes.laplacian, -1);
    const SolverReport rep = eqn.solve(fvMatrix::Controls{1e-13, 20000, 2});
    out.iterations = rep.iterations;
    out.residual = rep.finalResidual;
    out.solved = rep.converged;
    if (!rep.converged)
        cfFatal("cht: region temperature solve did not converge ("
                + std::to_string(rep.iterations) + " iters)");
    T.internalField() = eqn.solution()[0];
    for (label c = 0; c < nC; ++c)
        out.T[static_cast<std::size_t>(c)] = T.internalField()[static_cast<std::size_t>(c)];

    // Conductive heat flux on a named patch: average of -k * oc * (T_face - T_owner)
    // over the patch faces, per unit face area (the sign convention is the flux
    // leaving the region through that face along its outward normal Sf).
    auto patchFlux = [&](const std::string& name) {
        const label pid = mesh.patchIdByName(name);
        const polyMesh::Patch& pa = mesh.patches()[static_cast<std::size_t>(pid)];
        const auto& sf = mesh.faceSf();
        const auto& dc = mesh.deltaCoeffVec();
        scalar flux = 0, area = 0;
        for (label f = pa.start; f < pa.start + pa.size; ++f)
        {
            const std::size_t u = static_cast<std::size_t>(f);
            const label owner = mesh.owner()[u];
            const Vector& d = dc[u];
            const scalar d2 = dot(d, d);
            if (!(d2 > 0)) continue;
            const scalar oc = dot(sf[u], d) / d2;
            scalar c = 0, kk = 0;
            T.linearise(f, c, kk);
            const scalar tOwner = T.internalField()[static_cast<std::size_t>(owner)];
            const scalar dT = (c - 1.0) * tOwner + kk; // T_face - T_owner
            // Flux vector q = -k grad T; q . Sf = -k oc (T_face - T_owner).
            flux += -k * oc * dT;
            area += mag(sf[u]);
        }
        return area > 0 ? flux / area : 0.0;
    };
    out.interfaceFlux = patchFlux(interfacePatch);
    out.outerFlux = patchFlux(outerPatch);
    return out;
}

ChtResult solveCht(const polyMesh& solidMesh, const polyMesh& fluidMesh,
                   const Schemes& schemes, const ChtConfig& cfg,
                   const std::string& solidOuter, const std::string& solidIface,
                   const std::string& fluidOuter, const std::string& fluidIface)
{
    cfg.validate();
    ChtResult out;
    // Secant on T_i: define F(T_i) = q_s(T_i) - q_f(T_i). q_s increases with
    // T_i (hotter interface reduces solid->interface flux), and q_f decreases
    // with T_i, so F is monotone decreasing -- a bracketed secant is robust.
    scalar ti = 0.5 * (cfg.solid.outerT + cfg.fluid.outerT);
    scalar tiPrev = ti;
    scalar fPrev = std::numeric_limits<scalar>::quiet_NaN();
    scalar fNow = std::numeric_limits<scalar>::quiet_NaN();

    for (int it = 1; it <= cfg.maxIterations; ++it)
    {
        const RegionResult rs = solveRegion(solidMesh, schemes, cfg.solid.k,
                                            cfg.solid.outerT, ti, solidOuter, solidIface);
        const RegionResult rf = solveRegion(fluidMesh, schemes, cfg.fluid.k,
                                            cfg.fluid.outerT, ti, fluidOuter, fluidIface);
        // On the solid interface face Sf points +x (out of solid): q_s . n.
        // On the fluid interface face Sf points -x (out of fluid): q_f . n.
        // Continuity of the physical +x flux requires q_s.n + q_f.n = 0, i.e.
        // the two outward fluxes are equal and opposite. We compare their sum.
        fNow = rs.interfaceFlux + rf.interfaceFlux;

        out.solid = rs;
        out.fluid = rf;
        out.interfaceTemperature = ti;
        out.outerIterations = it;
        out.fluxMismatch = std::abs(fNow);

        const scalar scale = std::max(std::abs(rs.outerFlux), std::abs(rf.outerFlux));
        const scalar rel = out.fluxMismatch / std::max(scale, 1e-30);
        if (rel <= cfg.interfaceTol && it > 1)
        {
            out.solved = true;
            break;
        }

        // Secant update using the previous two (T_i, F) pairs. On the first
        // iteration use a fixed small perturbation to get a second point.
        scalar tiNext;
        if (it == 1 || !std::isfinite(fPrev) || std::abs(fNow - fPrev) < 1e-30)
        {
            const scalar dT = 0.1 * std::abs(cfg.solid.outerT - cfg.fluid.outerT);
            tiNext = ti - 0.05 * dT * (fNow >= 0 ? 1.0 : -1.0);
        }
        else
        {
            tiNext = ti - fNow * (ti - tiPrev) / (fNow - fPrev);
        }
        // Clamp to the physical bracket and relax.
        const scalar lo = std::min(cfg.solid.outerT, cfg.fluid.outerT);
        const scalar hi = std::max(cfg.solid.outerT, cfg.fluid.outerT);
        tiNext = std::max(lo, std::min(hi, tiNext));
        tiPrev = ti;
        fPrev = fNow;
        ti = ti + cfg.relax * (tiNext - ti);
    }
    return out;
}

} // namespace cf
