// chflow-3d-compressible: density-based Euler solver, first increment.
// First-order Godunov with the HLLC approximate Riemann flux (Toro) on the
// polyhedral mesh, ideal-gas EOS. Every boundary face is an impermeable slip
// wall (wall/symmetry patches alike), so only mass and total energy are
// conserved exactly -- the walls exchange momentum, which the app accumulates
// as an impulse time series the momentum balance needs.
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "cf/common/Dictionary.h"
#include "cf/common/Log.h"
#include "cf/fvm/MeshLoader.h"
#include "cf/fvm/VtkWriter.h"

using namespace cf;

namespace
{

std::string num(double v)
{
    std::ostringstream os;
    os << std::setprecision(12) << v;
    return os.str();
}

/// Per-volume conserved state [rho, rho*u, rho*E]; E is the specific total
/// energy. Fluxes and star states are built from this triple.
struct State
{
    scalar rho = 0;
    Vector mom = Vector::zero();
    scalar rhoE = 0;
};

struct Flux
{
    scalar mass = 0;
    Vector mom = Vector::zero();
    scalar energy = 0;
};

/// Physical flux through a face with unit normal nHat, per unit area:
/// [rho*un, rho*un*u + p*nHat, (rhoE + p)*un].
Flux physicalFlux(const State& u, scalar p, scalar un, const Vector& nHat)
{
    return {u.rho * un, u.mom * un + p * nHat, (u.rhoE + p) * un};
}

/// The pressure the two-rarefaction Riemann solver would report: exact when
/// both sides release a rarefaction, and the star-pressure estimate that tells
/// the wave-speed formula below which side actually carries a shock.
scalar twoRarefactionPressure(scalar pL, scalar unL, scalar cL,
                              scalar pR, scalar unR, scalar cR, scalar gamma)
{
    const scalar k = (gamma - 1.0) / (2.0 * gamma);
    const scalar num = 0.5 * (gamma - 1.0) * (unL - unR) + cL + cR;
    const scalar den = cL * std::pow(pL, -k) + cR * std::pow(pR, -k);
    return std::pow(num / den, 1.0 / k);
}

/// HLLC approximate Riemann flux per unit area across a face whose unit
/// normal nHat points from the left (owner) to the right (neighbour) state.
Flux hllc(const State& uL, scalar pL, const Vector& vL, scalar cL,
          const State& uR, scalar pR, const Vector& vR, scalar cR,
          const Vector& nHat, scalar gamma)
{
    const scalar unL = dot(vL, nHat);
    const scalar unR = dot(vR, nHat);
    // Toro's shock-adaptive speeds. The acoustic bounds min(un - c) / max(un + c)
    // do not bracket a strong shock: at the Sod diaphragm the true shock runs at
    // 1.75 while they say 1.18, and a flux whose fan is narrower than the wave it
    // carries cannot preserve the Rankine-Hugoniot jump across it. q_K widens the
    // side the estimated star pressure compresses -- exactly the Rankine-Hugoniot
    // shock speed when p* is right, and 1 (the rarefaction head) when it does not.
    const scalar pStar = twoRarefactionPressure(pL, unL, cL, pR, unR, cR, gamma);
    const auto q = [pStar, gamma](scalar pK)
    {
        if (pStar <= pK) return scalar{1};
        return std::sqrt((gamma + 1.0) / (2.0 * gamma) * (pStar / pK)
                       + (gamma - 1.0) / (2.0 * gamma));
    };
    const scalar sL = unL - cL * q(pL);
    const scalar sR = unR + cR * q(pR);
    const scalar sM = (pR - pL + uL.rho * unL * (sL - unL) - uR.rho * unR * (sR - unR))
                    / (uL.rho * (sL - unL) - uR.rho * (sR - unR));

    const Flux fL = physicalFlux(uL, pL, unL, nHat);
    const Flux fR = physicalFlux(uR, pR, unR, nHat);

    // star states; the tangential velocity rides through the contact unchanged
    const Vector tL = vL - unL * nHat;
    const Vector tR = vR - unR * nHat;
    const scalar rhoLs = uL.rho * (sL - unL) / (sL - sM);
    const scalar rhoRs = uR.rho * (sR - unR) / (sR - sM);
    const State usL{rhoLs, rhoLs * (sM * nHat + tL),
                    rhoLs * (uL.rhoE / uL.rho + (sM - unL) * (sM + pL / (uL.rho * (sL - unL))))};
    const State usR{rhoRs, rhoRs * (sM * nHat + tR),
                    rhoRs * (uR.rhoE / uR.rho + (sM - unR) * (sM + pR / (uR.rho * (sR - unR))))};

    if (sL >= 0.0) return fL;
    if (sR <= 0.0) return fR;

    // flux* + s*(U* - U) on the supersonic side of the contact
    const State uSel = sM >= 0.0 ? usL : usR;
    const Flux fSel = sM >= 0.0 ? fL : fR;
    const scalar sSel = sM >= 0.0 ? sL : sR;
    const State uOwn = sM >= 0.0 ? uL : uR;
    return {fSel.mass + sSel * (uSel.rho - uOwn.rho),
            fSel.mom + sSel * (uSel.mom - uOwn.mom),
            fSel.energy + sSel * (uSel.rhoE - uOwn.rhoE)};
}

} // namespace

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        log::error("usage: chflow-3d-compressible <case.dict>");
        return 2;
    }
    try
    {
        const Dictionary caseDict = Dictionary::parseFile(argv[1]);
        polyMesh mesh = MeshLoader::fromCase(caseDict, argv[1]);

        const Dictionary& fluid = caseDict.getDict("fluid");
        const scalar gamma = fluid.getScalar("gamma");

        const Dictionary& ini = caseDict.getDict("initial");
        const scalar xDiaphragm = ini.getScalar("diaphragm");
        const Dictionary& iniL = ini.getDict("left");
        const Dictionary& iniR = ini.getDict("right");
        const scalar rhoLi = iniL.getScalar("rho");
        const scalar pLi = iniL.getScalar("p");
        const scalar uLi = iniL.getScalar("u");
        const scalar rhoRi = iniR.getScalar("rho");
        const scalar pRi = iniR.getScalar("p");
        const scalar uRi = iniR.getScalar("u");

        const Dictionary& timeD = caseDict.getDict("time");
        const scalar dt = timeD.getScalar("dt");
        const scalar endTime = timeD.getScalar("endTime");
        const int reportEvery = static_cast<int>(timeD.getScalar("reportEvery"));

        const std::string vtk = caseDict.has("output")
                                  ? caseDict.getDict("output").getWordOr("vtk", "")
                                  : "";

        const label n = mesh.nCells();
        const Field<scalar>& vol = mesh.cellVolumes();
        const Field<Vector>& sf = mesh.faceSf();
        const Field<scalar>& magSf = mesh.faceAreas();
        const Field<label>& owner = mesh.owner();
        const Field<label>& neighbour = mesh.neighbour();
        const label nInt = mesh.nInternalFaces();

        // characteristic cell size V/maxFaceArea, for the CFL guard
        scalar dxMin = std::numeric_limits<scalar>::max();
        const auto& cellFaces = mesh.cells();
        for (label i = 0; i < n; ++i)
        {
            scalar aMax = 0;
            for (label f : cellFaces[i]) aMax = std::max(aMax, magSf[f]);
            dxMin = std::min(dxMin, vol[i] / aMax);
        }

        std::vector<scalar> massV(n), energyV(n);
        std::vector<Vector> momV(n);
        for (label i = 0; i < n; ++i)
        {
            const Vector& xc = mesh.cellCentres()[i];
            const bool leftSide = xc.x < xDiaphragm;
            const scalar rho = leftSide ? rhoLi : rhoRi;
            const scalar p = leftSide ? pLi : pRi;
            const Vector u(leftSide ? uLi : uRi, 0.0, 0.0);
            massV[i] = rho * vol[i];
            momV[i] = rho * u * vol[i];
            energyV[i] = rho * (p / ((gamma - 1.0) * rho) + 0.5 * dot(u, u)) * vol[i];
        }

        scalar mass0 = 0, energy0 = 0, mom0x = 0;
        for (label i = 0; i < n; ++i)
        {
            mass0 += massV[i];
            energy0 += energyV[i];
            mom0x += momV[i].x;
        }

        const int nSteps = static_cast<int>(std::lround(endTime / dt));
        log::info("euler start: cells ", n, " mass ", num(mass0), " mom ", num(mom0x),
                  " energy ", num(energy0), " dt ", num(dt), " endTime ", num(endTime),
                  " gamma ", num(gamma));

        std::vector<scalar> rho(n), p(n), c(n);
        std::vector<Vector> vel(n);

        auto primitives = [&]()
        {
            for (label i = 0; i < n; ++i)
            {
                rho[i] = massV[i] / vol[i];
                vel[i] = momV[i] / massV[i];
                const scalar e = energyV[i] / massV[i] - 0.5 * dot(vel[i], vel[i]);
                p[i] = (gamma - 1.0) * rho[i] * e;
                c[i] = std::sqrt(gamma * p[i] / rho[i]);
            }
        };

        bool ok = true;
        scalar t = 0;
        scalar impulseX = 0;
        for (int step = 1; step <= nSteps && ok; ++step)
        {
            primitives();
            scalar amax = 0;
            for (label i = 0; i < n; ++i) amax = std::max(amax, mag(vel[i]) + c[i]);
            if (!std::isfinite(amax) || amax <= 0)
            {
                log::error("euler: non-finite state at step ", step);
                ok = false;
                break;
            }
            const scalar cfl = amax * dt / dxMin;
            if (cfl > 1.0 + 1e-12)
            {
                log::error("euler: CFL violated (dt * amax / dxMin = ", num(cfl), ")");
                ok = false;
                break;
            }

            std::vector<scalar> dm(n, 0.0), dE(n, 0.0);
            std::vector<Vector> dmom(n, Vector::zero());

            for (label f = 0; f < nInt; ++f)
            {
                const label o = owner[f], nn = neighbour[f];
                const Vector nHat = sf[f] / magSf[f];
                const Flux phi = hllc({rho[o], rho[o] * vel[o], energyV[o] / vol[o]},
                                      p[o], vel[o], c[o],
                                      {rho[nn], rho[nn] * vel[nn], energyV[nn] / vol[nn]},
                                      p[nn], vel[nn], c[nn], nHat, gamma);
                const Flux face{phi.mass * magSf[f], phi.mom * magSf[f], phi.energy * magSf[f]};
                dm[o] -= face.mass;
                dmom[o] -= face.mom;
                dE[o] -= face.energy;
                dm[nn] += face.mass;
                dmom[nn] += face.mom;
                dE[nn] += face.energy;
            }
            for (label f = nInt; f < mesh.nFaces(); ++f)
            {
                const label o = owner[f];
                // impermeable slip wall: no mass or energy crosses, the wall
                // takes the pressure force the fluid loses
                dmom[o] -= p[o] * sf[f];
                impulseX += dt * (p[o] * sf[f]).x;
            }
            for (label i = 0; i < n; ++i)
            {
                // the totals are extensive: dU/dt = -div Phi, so dt * flux
                massV[i] += dt * dm[i];
                momV[i] += dt * dmom[i];
                energyV[i] += dt * dE[i];
            }
            t += dt;

            if (step % reportEvery != 0 && step != nSteps) continue;
            primitives();
            scalar pmin = p[0], pmax = p[0];
            for (label i = 0; i < n; ++i)
            {
                pmin = std::min(pmin, p[i]);
                pmax = std::max(pmax, p[i]);
            }
            scalar mass = 0, energy = 0, momx = 0;
            for (label i = 0; i < n; ++i)
            {
                mass += massV[i];
                energy += energyV[i];
                momx += momV[i].x;
            }
            log::info("euler t ", num(t), " mass ", num(mass), " mom ", num(momx),
                      " energy ", num(energy), " impulse ", num(impulseX),
                      " pmin ", num(pmin), " pmax ", num(pmax), " cfl ", num(cfl));
            if (!std::isfinite(mass) || !std::isfinite(energy) ||
                !std::isfinite(pmin) || !std::isfinite(pmax))
            {
                log::error("euler: non-finite field at step ", step);
                ok = false;
                break;
            }
        }

        if (ok && !vtk.empty())
        {
            Field<scalar> rhoF(n), pF(n);
            Field<Vector> uF(n);
            primitives();
            for (label i = 0; i < n; ++i)
            {
                rhoF[i] = rho[i];
                pF[i] = p[i];
                uF[i] = vel[i];
            }
            VtkWriter w(mesh);
            w.addCellScalar("rho", rhoF);
            w.addCellScalar("p", pF);
            w.addCellVector("U", uF);
            w.write(vtk);
        }

        return ok ? 0 : 1;
    }
    catch (const std::exception& e)
    {
        log::error("euler: ", e.what());
        return 1;
    }
}
