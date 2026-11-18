#pragma once

#include <map>
#include <string>

#include "cf/fvm/PolyMesh.h"

namespace cf
{

/// Reader for the ASCII representation of an OpenFOAM polyMesh.
///
/// Reads the five classic files from a `constant/polyMesh` directory:
///
///   points      list of (x y z)
///   faces       list of N(v0 v1 ...)
///   owner       one cell id per face
///   neighbour   one cell id per internal face (prefix of `faces`)
///   boundary    list of patches, each with type/nFaces/startFace
///
/// The reader is deliberately conservative: it accepts both the `FoamFile {
/// ... }` header and a bare `N ( ... )` list, tolerates arbitrary whitespace,
/// extra blanks and `//` or `/* */` comments, and matches keywords
/// case-insensitively. It never silently truncates a malformed list: declared
/// counts must match the entries actually read, face vertices must exist and
/// owner/neighbour must address a cell that has at least one face. Binary
/// OpenFOAM meshes are rejected with a clear message.
///
/// Boundary patch types are mapped onto the kernel's vocabulary:
///
///   wall / "67" / patch        -> wall
///   empty / symmetryPlane      -> symmetry
///   anything else              -> kept verbatim (no guessing)
///
/// A caller may override any patch type via `patchTypes` (keyed by patch name).
class OpenFoamMesh
{
public:
    /// `directory` is the polyMesh directory itself. If it does not contain
    /// `points` directly but does contain `constant/polyMesh/points`, that
    /// subdirectory is used instead, so a case root can be passed as well.
    static polyMesh fromDirectory(
        const std::string& directory,
        const std::map<std::string, std::string>& patchTypes = {});

    /// Read from an explicit set of five file paths (mainly for testing).
    static polyMesh fromFiles(
        const std::string& pointsFile,
        const std::string& facesFile,
        const std::string& ownerFile,
        const std::string& neighbourFile,
        const std::string& boundaryFile,
        const std::map<std::string, std::string>& patchTypes = {});
};

} // namespace cf
