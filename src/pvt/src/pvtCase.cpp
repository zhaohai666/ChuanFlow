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
        if (node.kind != Dictionary::Kind::scalar)
            cfFatal("pvt: composition entry '" + key + "' must be a number");
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
