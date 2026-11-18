#include <cmath>
#include <vector>

#include "cf/common/Error.h"
#include "cf/pvt/ComponentData.h"

namespace cf::pvt
{

// Weinaug & Katz (1943) Table 18.1, cgs parachor [cm3/mol (dyn/cm)^0.25].
// See ComponentData.h for the full provenance. Exact one-decimal values as
// published; no interpolation, no extrapolation.
scalar parachorFromTableCgs(const std::string& name, bool& found)
{
    static const struct
    {
        const char* name;
        scalar value;
    } table[] = {
        {"carbonDioxide", 78.0},   // CO2
        {"nitrogen", 41.0},        // N2
        {"methane", 77.0},         // C1
        {"ethane", 108.0},         // C2
        {"propane", 150.3},        // C3
        {"iButane", 181.5},        // iC4
        {"nButane", 189.9},        // nC4
        {"iPentane", 225.0},       // iC5 (not in the registry; kept for record)
        {"nPentane", 231.5},       // nC5
        {"nHexane", 271.0},        // nC6
        {"nHeptane", 312.5},       // nC7 (not in the registry; kept for record)
        {"nOctane", 351.5},        // nC8 (not in the registry; kept for record)
    };
    for (const auto& row : table)
    {
        if (name == row.name)
        {
            found = true;
            return row.value;
        }
    }
    found = false;
    return 0.0;
}

scalar requireParachorCgs(const std::string& name)
{
    bool found = false;
    const scalar p = parachorFromTableCgs(name, found);
    if (!found)
        cfFatal("pvt: component '" + name
                + "' has no Weinaug-Katz (1943) tabulated parachor; declare one explicitly");
    return p;
}

std::vector<scalar> parachorsFromTableSi(const std::vector<Component>& comps)
{
    // P_SI = P_cgs * 1e-6.75 (SurfaceTension.h documents the conversion).
    static const scalar kCgsToSi = std::pow(10.0, -6.75);
    std::vector<scalar> out;
    out.reserve(comps.size());
    for (const Component& c : comps)
        out.push_back(requireParachorCgs(c.name) * kCgsToSi);
    return out;
}

} // namespace cf::pvt
