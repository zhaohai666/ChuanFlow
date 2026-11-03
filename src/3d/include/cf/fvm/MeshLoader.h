#pragma once

#include <string>

#include "cf/common/Dictionary.h"
#include "cf/fvm/PolyMesh.h"

namespace cf
{

/// Case-file mesh selection.
///
///   block { min (...); max (...); cells (...); xmin { ... } ... }
///
///   mesh {
///       format gmsh;                       # only gmsh for now
///       file "meshes/channel.msh";         # ASCII .msh 2.2 or 4.1
///       patches {
///           inlet  { type inlet; }         # default type: wall
///           outlet { type outlet; }
///       }
///   }
///
/// A relative `file` is resolved against the case file's directory first, so
/// a case can ship its mesh next to its .dict and still run from any cwd.
class MeshLoader
{
public:
    /// `casePath` is the .dict the dictionary was parsed from (may be empty
    /// when it came from a string); only its directory is used.
    static polyMesh fromCase(const Dictionary& caseDict, const std::string& casePath);
};

} // namespace cf
