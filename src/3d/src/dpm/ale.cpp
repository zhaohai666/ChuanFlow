#include <algorithm>
#include <cmath>

#include "cf/common/Error.h"
#include "cf/fvm/Ale.h"

namespace cf
{

void AleConfig::validate() const
{
    if (motion != "translate" && motion != "expand")
        cfFatal("ale: motion must be 'translate' or 'expand' (got '" + motion + "')");
    if (!std::isfinite(dt) || !(dt > 0))
        cfFatal("ale: dt must be finite and positive");
    if (steps < 1)
        cfFatal("ale: steps must be at least 1");
    if (reportEvery < 1)
        cfFatal("ale: reportEvery must be at least 1");
    if (motion == "expand" && (!std::isfinite(expansionRate) || expansionRate <= 0))
        cfFatal("ale: expansionRate must be finite and positive for motion=expand");
}

AleConfig readAleConfig(const Dictionary& caseDict)
{
    AleConfig cfg;
    if (!caseDict.has("ale")) return cfg;
    const Dictionary& ad = caseDict.getDict("ale");
    cfg.enabled = true;
    cfg.motion = ad.getWordOr("motion", "translate");
    if (ad.has("translation"))
    {
        const auto v = ad.getScalarList("translation");
        if (v.size() != 3) cfFatal("ale: translation needs three components");
        cfg.translation = Vector{v[0], v[1], v[2]};
    }
    cfg.expansionRate = ad.getScalarOr("expansionRate", 0.0);
    cfg.dt = ad.getScalarOr("dt", 0.0);
    cfg.steps = static_cast<int>(ad.getScalarOr("steps", 0.0));
    cfg.reportEvery = static_cast<int>(ad.getScalarOr("reportEvery", 1.0));
    cfg.validate();
    return cfg;
}

Field<Vector> meshVertexVelocity(const polyMesh& mesh, const AleConfig& cfg)
{
    cfg.validate();
    const Field<Vector>& pts = mesh.points();
    Field<Vector> v(pts.size(), Vector::zero());
    for (std::size_t i = 0; i < pts.size(); ++i)
    {
        if (cfg.motion == "translate")
            v[i] = cfg.translation;
        else
            v[i] = pts[i] * cfg.expansionRate; // U_m = rate * x
    }
    return v;
}

Field<Vector> advanceVertices(const polyMesh& mesh, const AleConfig& cfg)
{
    cfg.validate();
    const Field<Vector>& pts = mesh.points();
    Field<Vector> np(pts.size(), Vector::zero());
    for (std::size_t i = 0; i < pts.size(); ++i)
    {
        if (cfg.motion == "translate")
            np[i] = pts[i] + cfg.translation * cfg.dt;
        else
            np[i] = pts[i] * (1.0 + cfg.expansionRate * cfg.dt);
    }
    return np;
}

Field<scalar> gclResidual(const polyMesh& mesh, const Field<scalar>& oldVolumes,
                          const Field<Vector>& vertexVelocity, scalar dt)
{
    if (!(dt > 0)) cfFatal("ale: GCL needs a positive dt");
    const label nC = mesh.nCells();
    Field<scalar> res(static_cast<std::size_t>(nC), 0.0);

    const auto& faces = mesh.faces();
    const auto& sf = mesh.faceSf();

    // Mesh-face velocity = mean of its vertices' velocities (exact for a
    // linearly-moving face).
    Field<Vector> faceVel(static_cast<std::size_t>(mesh.nFaces()), Vector::zero());
    for (label f = 0; f < mesh.nFaces(); ++f)
    {
        const auto& verts = faces[static_cast<std::size_t>(f)];
        Vector sum = Vector::zero();
        for (label vid : verts) sum = sum + vertexVelocity[static_cast<std::size_t>(vid)];
        if (!verts.empty())
            faceVel[static_cast<std::size_t>(f)] = sum * (1.0 / static_cast<scalar>(verts.size()));
    }

    // Sum signed U_m . Sf per cell with the NEW geometry (owner +, neighbour -).
    for (label f = 0; f < mesh.nFaces(); ++f)
    {
        const std::size_t u = static_cast<std::size_t>(f);
        const scalar flux = dot(faceVel[u], sf[u]);
        const label owner = mesh.owner()[u];
        res[static_cast<std::size_t>(owner)] += flux;
        if (f < mesh.nInternalFaces())
        {
            const label nb = mesh.neighbour()[u];
            res[static_cast<std::size_t>(nb)] -= flux;
        }
    }

    const auto& vol = mesh.cellVolumes(); // triggers recompute after movePoints
    for (label c = 0; c < nC; ++c)
    {
        const std::size_t u = static_cast<std::size_t>(c);
        const scalar dVdt = (vol[u] - oldVolumes[u]) / dt;
        res[u] = dVdt - res[u];
    }
    return res;
}

} // namespace cf
