#include <cmath>
#include <string>

#include "cf/common/Error.h"
#include "cf/pvt/PvtCase.h"

namespace cf::pvt
{

Mixture loadComposition(const Dictionary& compDict)
{
    Mixture mix;
    for (const std::string& key : compDict.keys())
    {
        const Dictionary::Node& node = compDict.find(key);
        // 5.238: a block entry declares a CUSTOM component -- the built-in
        // tables stop at water (omega 0.344), and the PR1978 alpha branch is
        // for omega > 0.491, so the heavier fractions a field model carries
        // (C7+, H2S-rich runs) must be declarable. fraction is the mole
        // fraction; Tc/Pc/omega/molarMass are the four EOS constants. A
        // scalar entry is the historical built-in lookup, byte-identical.
        if (node.kind == Dictionary::Kind::block)
        {
            const Dictionary& cd = *node.block;
            Component c;
            c.name = key;
            c.Tc = cd.getScalar("Tc");
            c.Pc = cd.getScalar("Pc");
            c.acentric = cd.getScalar("omega");
            c.molarMass = cd.getScalar("molarMass");
            const double fraction = cd.getScalar("fraction");
            if (!(c.Tc > 0.0) || !(c.Pc > 0.0) || !(c.molarMass > 0.0))
                cfFatal("pvt: custom component '" + key
                        + "' needs positive Tc, Pc and molarMass");
            if (!std::isfinite(c.acentric))
                cfFatal("pvt: custom component '" + key
                        + "' needs a finite omega");
            if (!(fraction > 0.0))
                cfFatal("pvt: custom component '" + key
                        + "' needs a positive fraction");
            mix.add(c, fraction);
            continue;
        }
        if (node.kind != Dictionary::Kind::scalar)
            cfFatal("pvt: composition entry '" + key + "' must be a number"
                    " or a custom-component block");
        mix.add(requireComponent(key), node.scalar);
    }
    if (mix.size() == 0) cfFatal("pvt: empty composition block");
    mix.normalise();
    return mix;
}

FlashCase loadFlashCase(const Dictionary& caseDict)
{
    FlashCase c;
    if (!caseDict.has("composition"))
        cfFatal("pvt: case needs a composition { } block");
    c.feed = loadComposition(caseDict.getDict("composition"));

    // Optional declared binary-interaction parameters: a name-based symmetric
    // matrix, binaryInteraction { methane { carbonDioxide 0.10; water 0.50; } }.
    // Absent == kij = 0 everywhere, which is bit-identical to the legacy rule.
    if (caseDict.has("binaryInteraction"))
    {
        const Dictionary& bi = caseDict.getDict("binaryInteraction");
        for (const std::string& ni : bi.keys())
        {
            const std::size_t i = c.feed.index(ni);
            const Dictionary::Node& nodeI = bi.find(ni);
            if (nodeI.kind != Dictionary::Kind::block)
                cfFatal("pvt: binaryInteraction entry '" + ni
                        + "' must be a block of { component value; } pairs");
            const Dictionary& row = *nodeI.block;
            for (const std::string& nj : row.keys())
            {
                const std::size_t j = c.feed.index(nj);
                c.feed.setBinaryInteraction(i, j, row.getScalar(nj));
            }
        }
    }

    c.temperature = caseDict.getScalarOr("temperature", 300.0);
    if (caseDict.has("pressures"))
        c.pressures = caseDict.getScalarList("pressures");
    else if (caseDict.has("pressure"))
        c.pressures = {caseDict.getScalar("pressure")};
    else
        cfFatal("pvt: case needs pressure or pressures ( ... )");

    if (caseDict.has("solve"))
    {
        const Dictionary& solve = caseDict.getDict("solve");
        c.controls.tol = solve.getScalarOr("tol", c.controls.tol);
        if (solve.has("maxIter")) c.controls.maxIter = solve.getInt("maxIter");
        if (solve.has("bubblePoint")) c.computeBubblePoint = solve.getBool("bubblePoint");
    }
    return c;
}

} // namespace cf::pvt
