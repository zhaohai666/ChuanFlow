#include <cctype>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

#include "cf/common/Dictionary.h"
#include "cf/common/Error.h"
#include "cf/common/Log.h"
#include "cf/fvm/OpenFoamMesh.h"
#include "cf/fvm/PolyMesh.h"

using namespace cf;

namespace
{

std::string lower(std::string s)
{
    for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

std::string num(scalar v)
{
    std::ostringstream os;
    os << std::setprecision(15) << v;
    return os.str();
}

std::string vec(const Vector& v)
{
    return num(v.x) + " " + num(v.y) + " " + num(v.z);
}

void usage()
{
    std::cerr << "usage: chflow-openfoam-mesh <case.dict>\n"
              << "  mesh { format openfoam; directory <polyMesh dir>; }\n";
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        usage();
        return 2;
    }

    try
    {
        const Dictionary caseDict = Dictionary::parseFile(argv[1]);
        if (!caseDict.has("mesh"))
            cfFatal("openfoam-mesh: case has no mesh {} block");

        const Dictionary& md = caseDict.getDict("mesh");
        const std::string format = md.getWordOr("format", "openfoam");
        if (lower(format) != "openfoam")
        {
            cfFatal("openfoam-mesh: mesh format '" + format
                    + "' is not supported (this binary reads OpenFOAM ASCII polyMesh)");
        }

        std::string directory;
        if (md.has("directory")) directory = md.getWord("directory");
        else if (md.has("file")) directory = md.getWord("file");
        else cfFatal("openfoam-mesh: mesh {} needs 'directory <polyMesh dir>'");

        std::filesystem::path path(directory);
        if (path.is_relative() && !std::filesystem::exists(path))
        {
            const std::filesystem::path beside =
                std::filesystem::path(argv[1]).parent_path() / path;
            if (std::filesystem::exists(beside)) path = beside;
        }
        if (std::filesystem::is_regular_file(path)) path = path.parent_path();

        std::map<std::string, std::string> patchTypes;
        if (md.has("patches"))
        {
            const Dictionary& ps = md.getDict("patches");
            for (const auto& name : ps.keys())
                patchTypes[name] = ps.getDict(name).getWordOr("type", "wall");
        }

        const polyMesh mesh = OpenFoamMesh::fromDirectory(path.string(), patchTypes);

        log::info("openfoam-mesh: nPoints ", mesh.nPoints(), ";");
        log::info("openfoam-mesh: nCells ", mesh.nCells(), ";");
        log::info("openfoam-mesh: nFaces ", mesh.nFaces(), ";");
        log::info("openfoam-mesh: nInternalFaces ", mesh.nInternalFaces(), ";");
        log::info("openfoam-mesh: nBoundaryFaces ", mesh.nBoundaryFaces(), ";");

        for (const auto& p : mesh.patches())
        {
            log::info("openfoam-mesh: patch ", p.name, " type ", p.type,
                      " start ", p.start, " size ", p.size, ";");
        }

        scalar volume = 0;
        for (label c = 0; c < mesh.nCells(); ++c)
        {
            const scalar v = mesh.cellVolumes()[static_cast<std::size_t>(c)];
            volume += v;
            log::info("openfoam-mesh: cellVolume ", c, " ", num(v), ";");
        }
        log::info("openfoam-mesh: cellVolumeSum ", num(volume), ";");

        log::info("openfoam-mesh: minX ", vec(mesh.minX()), ";");
        log::info("openfoam-mesh: maxX ", vec(mesh.maxX()), ";");

        const auto& areas = mesh.faceAreas();
        for (label f = 0; f < mesh.nFaces(); ++f)
        {
            log::info("openfoam-mesh: faceArea ", f, " ",
                      num(areas[static_cast<std::size_t>(f)]), ";");
        }

        return 0;
    }
    catch (const std::exception& e)
    {
        log::error("openfoam-mesh: ", e.what());
        return 1;
    }
}
