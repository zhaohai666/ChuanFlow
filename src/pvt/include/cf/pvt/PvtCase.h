#pragma once

#include <vector>

#include "cf/common/Dictionary.h"
#include "cf/pvt/Components.h"
#include "cf/pvt/Flash.h"

namespace cf::pvt
{

/// Case grammar (see tests/cases/pvtFlash/case.dict):
///
///   temperature 320.0;                       // K
///   pressures (1.0e6 4.0e6 8.0e6);           // Pa, flash each
///   composition { methane 0.34; ethane 0.18; ... }   // mole fractions
///   solve { tol 1.0e-12; maxIter 500; bubblePoint yes; }
struct FlashCase
{
    Mixture feed;
    scalar temperature = 300.0;
    std::vector<scalar> pressures;
    FlashControls controls;
    bool computeBubblePoint = false;
};

FlashCase loadFlashCase(const Dictionary& caseDict);

/// Load a mixture from a `composition { name fraction; ... }` block.
Mixture loadComposition(const Dictionary& compDict);

} // namespace cf::pvt
