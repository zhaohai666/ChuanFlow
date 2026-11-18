#pragma once

#include <map>
#include <string>

#include "cf/common/Field.h"
#include "cf/fvm/PolyMesh.h"

namespace cf
{

/// Optional HDF5 result writer (skeleton).
///
/// Compiled only when the build defines CF_HAVE_HDF5 and links the HDF5 C
/// library. Without it the function still exists (so callers link unchanged)
/// but throws a clear error instead of silently dropping results.
///
/// Layout (when enabled):
///   /points                 nPoints x 3 double
///   /owner, /neighbour      face addressing (int64)
///   /faceOwnerStart? -- not written; topology is fully given by owner/neighbour
///   /fields/<name>          one dataset per cell scalar, length nCells
///
/// `fields` values must have one entry per cell; a size mismatch is fatal.
void writeFieldsHdf5(const std::string& path,
                     const polyMesh& mesh,
                     const std::map<std::string, Field<scalar>>& fields);

} // namespace cf
