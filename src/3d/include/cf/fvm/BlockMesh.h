#pragma once

#include <string>
#include <vector>

#include "cf/common/Dictionary.h"
#include "cf/fvm/PolyMesh.h"

namespace cf
{

/// Side patch descriptor for BlockMesh: (name, type) per box side.
struct BlockSides
{
    std::pair<std::string, std::string> xmin{"xmin", "wall"};
    std::pair<std::string, std::string> xmax{"xmax", "wall"};
    std::pair<std::string, std::string> ymin{"ymin", "wall"};
    std::pair<std::string, std::string> ymax{"ymax", "wall"};
    std::pair<std::string, std::string> zmin{"zmin", "symmetry"};
    std::pair<std::string, std::string> zmax{"zmax", "symmetry"};
};

/// Structured single-block hex mesher (the role of OpenFOAM's blockMesh for
/// our kernel's test cases): a box [x0,x1]x[y0,y1]x[z0,z1] split into
/// (nx, ny, nz) cells, with one named patch per outer side.
///
/// Sides are named by the coordinate they sit on:
///   xmin/xmax/ymin/ymax/zmin/zmax; sides with a single cell layer may be
///   merged into a named patch group via patchNames/patchTypes entries in
///   the case dictionary.
class BlockMesh
{
public:
    using Sides = BlockSides;

    static polyMesh create(const Vector& min, const Vector& max,
                           label nx, label ny, label nz,
                           const Sides& sides = {});

    /// Build from a case dictionary block:
    ///   block
    ///   {
    ///       min      (0 0 0);
    ///       max      (1 1 0.1);
    ///       cells    (10 10 1);
    ///       xmin     { name left;  type inlet; }
    ///       xmax     { name right; type outlet; }
    ///       zmin     { name front; type symmetry; }
    ///       zmax     { name back;  type symmetry; }
    ///   }
    /// Sides omitted default to wall (x,y) / symmetry (z).
    static polyMesh fromDict(const Dictionary& block);
};

} // namespace cf
