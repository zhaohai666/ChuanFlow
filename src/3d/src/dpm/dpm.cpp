#include "cf/fvm/Dpm.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

#include "cf/common/Error.h"

namespace cf
{

namespace
{

Vector vectorFromList(const std::vector<double>& v, const std::string& what)
{
    if (v.size() != 3)
        cfFatal("dpm: " + what + " needs exactly three components, got "
                + std::to_string(v.size()));
    return Vector{v[0], v[1], v[2]};
}

bool finiteVector(const Vector& v)
{
    return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

} // namespace

void DpmConfig::validate() const
{
    if (streams.empty())
        cfFatal("dpm: at least one stream0{ point; velocity; diameter; density; count; }"
                " block is required");
    if (!(dt > 0) || !std::isfinite(dt))
        cfFatal("dpm: dt must be finite and positive (s)");
    if (maxSteps < 1)
        cfFatal("dpm: maxSteps must be at least 1");
    if (!(carrierRho > 0) || !std::isfinite(carrierRho))
        cfFatal("dpm: carrierRho must be finite and positive (kg/m3)");
    if (!finiteVector(gravity))
        cfFatal("dpm: gravity must be finite");
    if (drag != "schillerNaumann")
        cfFatal("dpm: drag law '" + drag + "' is not implemented; this increment"
                " only owns the published Schiller-Naumann sphere law");
    if (reportEvery < 1 || traceSteps < 0 || traceId < 0)
        cfFatal("dpm: reportEvery must be >= 1 and traceId/traceSteps non-negative");
    for (std::size_t i = 0; i < streams.size(); ++i)
    {
        const DpmStream& s = streams[i];
        if (!finiteVector(s.point) || !finiteVector(s.velocity))
            cfFatal("dpm: stream " + std::to_string(i)
                    + " has a non-finite point or velocity");
        if (!(s.diameter > 0) || !std::isfinite(s.diameter))
            cfFatal("dpm: stream " + std::to_string(i)
                    + " diameter must be finite and positive (m)");
        if (!(s.density > 0) || !std::isfinite(s.density))
            cfFatal("dpm: stream " + std::to_string(i)
                    + " density must be finite and positive (kg/m3)");
        if (s.count < 1)
            cfFatal("dpm: stream " + std::to_string(i) + " count must be >= 1");
    }
}

DpmConfig readDpmConfig(const Dictionary& caseDict)
{
    if (!caseDict.has("dpm"))
        cfFatal("dpm: the chflow-3d-dpm app needs a dpm{} declaration");
    const Dictionary& d = caseDict.getDict("dpm");
    DpmConfig cfg;
    cfg.dt = d.getScalarOr("dt", cfg.dt);
    cfg.maxSteps = static_cast<int>(d.getScalarOr("maxSteps", cfg.maxSteps));
    cfg.carrierRho = d.getScalarOr("carrierRho", cfg.carrierRho);
    cfg.drag = d.getWordOr("drag", cfg.drag);
    cfg.reportEvery = static_cast<int>(d.getScalarOr("reportEvery", cfg.reportEvery));
    cfg.traceId = static_cast<int>(d.getScalarOr("traceId", cfg.traceId));
    cfg.traceSteps = static_cast<int>(d.getScalarOr("traceSteps", cfg.traceSteps));
    if (d.has("gravity")) cfg.gravity = vectorFromList(d.getScalarList("gravity"), "gravity");

    for (const std::string& key : d.keys())
    {
        if (key.rfind("stream", 0) != 0) continue;
        const Dictionary& sd = d.getDict(key);
        DpmStream s;
        s.point = vectorFromList(sd.getScalarList("point"), key + " point");
        s.velocity = vectorFromList(sd.getScalarList("velocity"), key + " velocity");
        s.diameter = sd.getScalar("diameter");
        s.density = sd.getScalar("density");
        s.count = static_cast<int>(sd.getScalar("count"));
        cfg.streams.push_back(s);
    }
    cfg.validate();
    return cfg;
}

scalar sphereDragCoefficient(scalar reynolds)
{
    if (!(reynolds > 0))
        cfFatal("dpm: particle Reynolds number must be positive");
    if (reynolds <= 1000.0)
        return (24.0 / reynolds) * (1.0 + 0.15 * std::pow(reynolds, 0.687));
    return 0.44;
}

Vector particleDragAcceleration(const Vector& uFluid, const Vector& uParticle,
                                scalar carrierRho, scalar rhoParticle,
                                scalar diameter, scalar mu)
{
    const Vector rel = uFluid - uParticle;
    const scalar relMag = mag(rel);
    if (!(relMag > 0)) return Vector::zero();
    const scalar re = carrierRho * relMag * diameter / mu;
    // Below Re = 1e-6 the published law is its Stokes limit; evaluating
    // 24/Re directly would multiply two huge/opposite numbers.
    if (re < 1.0e-6)
        return rel * (18.0 * mu / (rhoParticle * diameter * diameter));
    const scalar cd = sphereDragCoefficient(re);
    return rel * (3.0 * carrierRho * cd / (4.0 * rhoParticle * diameter));
}

const char* dpmStatusName(DpmStatus status)
{
    switch (status)
    {
    case DpmStatus::active: return "active";
    case DpmStatus::outlet: return "outlet";
    case DpmStatus::wall: return "wall";
    case DpmStatus::expired: return "expired";
    }
    return "unknown";
}

label locateDpmCell(const polyMesh& mesh, const Vector& point)
{
    const auto& centres = mesh.cellCentres();
    const Vector lo = mesh.minX();
    const Vector hi = mesh.maxX();
    const scalar scale = std::max({hi.x - lo.x, hi.y - lo.y, hi.z - lo.z, 1.0});
    const scalar tol = 1.0e-9 * scale;
    if (point.x < lo.x - tol || point.x > hi.x + tol
        || point.y < lo.y - tol || point.y > hi.y + tol
        || point.z < lo.z - tol || point.z > hi.z + tol)
        cfFatal("dpm: point (" + std::to_string(point.x) + " "
                + std::to_string(point.y) + " " + std::to_string(point.z)
                + ") is outside the mesh; particles must be injected inside it");

    // Zero-order sampling: the nearest CELL CENTRE owns the particle. For
    // the block/gmsh cells this is the Voronoi location and needs no face
    // half-space test; it also stays defined at the outlet face, where a
    // strict inside test can reject a particle the mesh itself still owns.
    label best = -1;
    scalar bestD2 = std::numeric_limits<scalar>::max();
    for (label c = 0; c < mesh.nCells(); ++c)
    {
        const scalar d2 = mag2(point - centres[static_cast<std::size_t>(c)]);
        if (d2 < bestD2)
        {
            bestD2 = d2;
            best = c;
        }
    }
    return best;
}

DpmResult runDpm(const polyMesh& mesh, const volVectorField& U,
                 const DpmConfig& cfg, scalar mu)
{
    if (!(mu > 0) || !std::isfinite(mu))
        cfFatal("dpm: carrier dynamic viscosity must be finite and positive");
    const auto& uCells = U.internalField();
    const Vector lo = mesh.minX();
    const Vector hi = mesh.maxX();
    DpmResult out;
    int id = 0;

    for (std::size_t s = 0; s < cfg.streams.size(); ++s)
    {
        const DpmStream& stream = cfg.streams[s];
        const label startCell = locateDpmCell(mesh, stream.point);
        for (int n = 0; n < stream.count; ++n, ++id)
        {
            DpmParticle p;
            p.id = id;
            p.stream = static_cast<int>(s);
            p.x0 = stream.point;
            p.u0 = stream.velocity;
            p.x = stream.point;
            p.u = stream.velocity;
            p.cell = startCell;
            p.diameter = stream.diameter;
            p.rhoParticle = stream.density;

            for (int step = 0; step <= cfg.maxSteps; ++step)
            {
                p.cell = locateDpmCell(mesh, p.x);
                const Vector uFluid = uCells[static_cast<std::size_t>(p.cell)];
                const scalar relMag = mag(uFluid - p.u);
                const scalar re = cfg.carrierRho * relMag * p.diameter / mu;
                if (re > 0 && (p.id == cfg.traceId)
                    && (step < cfg.traceSteps || step % cfg.reportEvery == 0))
                {
                    DpmSample sample;
                    sample.id = p.id;
                    sample.step = step;
                    sample.t = step * cfg.dt;
                    sample.x = p.x;
                    sample.u = p.u;
                    sample.uFluid = uFluid;
                    sample.reynolds = re;
                    sample.cd = sphereDragCoefficient(re);
                    sample.dragAcceleration = particleDragAcceleration(
                        uFluid, p.u, cfg.carrierRho, p.rhoParticle, p.diameter, mu);
                    out.samples.push_back(sample);
                }
                if (step == cfg.maxSteps)
                {
                    p.status = DpmStatus::expired;
                    break;
                }
                const Vector drag = particleDragAcceleration(
                    uFluid, p.u, cfg.carrierRho, p.rhoParticle, p.diameter, mu);
                const Vector body = (1.0 - cfg.carrierRho / p.rhoParticle) * cfg.gravity;
                p.u = p.u + cfg.dt * (drag + body);
                p.x = p.x + cfg.dt * p.u;
                p.age = (step + 1) * cfg.dt;

                const scalar r = 0.5 * p.diameter;
                if (p.x.x >= hi.x - r)
                {
                    p.status = DpmStatus::outlet;
                    break;
                }
                if (p.x.y <= lo.y + r || p.x.y >= hi.y - r
                    || p.x.z <= lo.z + r || p.x.z >= hi.z - r)
                {
                    p.status = DpmStatus::wall;
                    break;
                }
            }
            if (p.status == DpmStatus::outlet) ++out.outletCount;
            else if (p.status == DpmStatus::wall) ++out.wallCount;
            else ++out.expiredCount;
            out.particles.push_back(p);
        }
    }
    return out;
}

} // namespace cf
