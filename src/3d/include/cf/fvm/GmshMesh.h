#pragma once

#include <map>
#include <string>

#include "cf/fvm/PolyMesh.h"

namespace cf
{

/// Reader for Gmsh ASCII .msh files (formats 2.2 and 4.1).
///
/// First-order volume elements are supported: tetrahedra (4), hexahedra (5)
/// and prisms/wedges (6). Faces are recovered from the volume connectivity,
/// deduplicated by their vertex sets and addressed as owner/neighbour pairs,
/// so any mesh that polyMesh can represent can be imported - the reader does
/// not assume a structured block.
///
/// Patch names come from the Gmsh physical groups of the surface elements
/// (matching each boundary face by its vertices; a quad face also matches the
/// two triangles a triangulated surface splits it into). When a file carries
/// no surface information at all, boundary faces are named after the
/// bounding-box plane they lie on (xmin/xmax/ymin/ymax/zmin/zmax) and a
/// warning is logged. A boundary face that matches neither is an error: the
/// kernel refuses to guess a boundary condition.
///
/// Pyramids (7) and all second-order elements are rejected with a message
/// rather than silently truncated to their corner nodes.
class GmshMesh
{
public:
    /// `patchTypes` maps a patch name onto its boundary type
    /// (wall | symmetry | inlet | outlet | empty); names absent from the map
    /// default to "wall".
    static polyMesh fromFile(const std::string& path,
                             const std::map<std::string, std::string>& patchTypes = {});
};

} // namespace cf
