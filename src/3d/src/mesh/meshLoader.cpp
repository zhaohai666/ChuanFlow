#include "cf/fvm/MeshLoader.h"

#include <filesystem>
#include <map>
#include <string>

#include "cf/common/Error.h"
#include "cf/fvm/BlockMesh.h"
#include "cf/fvm/GmshMesh.h"

namespace cf
{

polyMesh MeshLoader::fromCase(const Dictionary& caseDict, const std::string& casePath)
{
    const bool hasFile = caseDict.has("mesh");
    const bool hasBlock = caseDict.has("block");
    if (hasFile && hasBlock)
    {
        cfFatal("case: block {} and mesh {} are mutually exclusive");
    }
    if (hasBlock)
    {
        return BlockMesh::fromDict(caseDict.getDict("block"));
    }
    if (!hasFile)
    {
        cfFatal("case: no mesh; add a block {} or a mesh {} section");
    }

    const Dictionary& md = caseDict.getDict("mesh");
    const std::string format = md.getWordOr("format", "gmsh");
    if (format != "gmsh")
    {
        cfFatal("mesh: format '" + format + "' is unknown (gmsh is supported)");
    }

    std::filesystem::path file(md.getWord("file"));
    if (file.is_relative() && !std::filesystem::exists(file) && !casePath.empty())
    {
        // a case ships its mesh next to the .dict, but the process cwd may be
        // elsewhere
        const std::filesystem::path besideCase =
            std::filesystem::path(casePath).parent_path() / file;
        if (std::filesystem::exists(besideCase)) file = besideCase;
    }

    std::map<std::string, std::string> patchTypes;
    if (md.has("patches"))
    {
        const Dictionary& ps = md.getDict("patches");
        for (const auto& name : ps.keys())
        {
            patchTypes[name] = ps.getDict(name).getWordOr("type", "wall");
        }
    }

    return GmshMesh::fromFile(file.string(), patchTypes);
}

} // namespace cf
