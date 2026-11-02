#include "cf/fvm/VtkWriter.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "cf/common/Error.h"
#include "cf/common/Log.h"

namespace cf
{

void VtkWriter::addCellScalar(const std::string& name, const Field<scalar>& values)
{
    if (static_cast<label>(values.size()) != mesh_.nCells())
    {
        cfFatal("vtk: scalar '" + name + "' size mismatch");
    }
    scalars_.emplace_back(name, &values);
}

void VtkWriter::addCellVector(const std::string& name, const Field<Vector>& values)
{
    if (static_cast<label>(values.size()) != mesh_.nCells())
    {
        cfFatal("vtk: vector '" + name + "' size mismatch");
    }
    vectors_.emplace_back(name, &values);
}

void VtkWriter::write(const std::string& path) const
{
    std::ofstream out(path);
    if (!out)
    {
        cfFatal("vtk: cannot open '" + path + "' for writing");
    }
    out << std::setprecision(12);

    out << "<?xml version=\"1.0\"?>\n"
        << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
        << "<UnstructuredGrid>\n"
        << "<Piece NumberOfPoints=\"" << mesh_.nPoints()
        << "\" NumberOfCells=\"" << mesh_.nCells() << "\">\n";

    // points
    out << "<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const Vector& p : mesh_.points())
    {
        out << p.x << ' ' << p.y << ' ' << p.z << '\n';
    }
    out << "</DataArray>\n</Points>\n";

    // cells: connectivity / offsets / types / faces (polyhedron -> type 42,
    // which requires the per-cell face list in VTK)
    out << "<Cells>\n";
    out << "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
    const auto& faces = mesh_.faces();
    const auto& cells = mesh_.cells();
    label offset = 0;
    std::vector<label> offsets;
    offsets.reserve(cells.size());
    std::ostringstream facesStream;
    for (const auto& cell : cells)
    {
        std::vector<label> verts;
        for (const label f : cell)
        {
            for (const label v : faces[static_cast<std::size_t>(f)])
            {
                if (std::find(verts.begin(), verts.end(), v) == verts.end())
                {
                    verts.push_back(v);
                }
            }
        }
        for (const label v : verts) out << v << ' ';
        out << '\n';
        offset += static_cast<label>(verts.size());
        offsets.push_back(offset);

        // face stream: nFaces n0 v.. n1 v..
        facesStream << cell.size() << ' ';
        for (const label f : cell)
        {
            const auto& fv = faces[static_cast<std::size_t>(f)];
            facesStream << fv.size() << ' ';
            for (const label v : fv) facesStream << v << ' ';
        }
        facesStream << '\n';
    }
    out << "</DataArray>\n";

    out << "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
    for (const label o : offsets) out << o << ' ';
    out << "</DataArray>\n";

    out << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (std::size_t c = 0; c < cells.size(); ++c) out << "42 ";
    out << "</DataArray>\n";

    out << "<DataArray type=\"Int64\" Name=\"faces\" format=\"ascii\">\n"
        << facesStream.str() << "</DataArray>\n</Cells>\n";

    // cell data
    if (!scalars_.empty() || !vectors_.empty())
    {
        out << "<CellData>\n";
        for (const auto& [name, vals] : scalars_)
        {
            out << "<DataArray type=\"Float64\" Name=\"" << name << "\" format=\"ascii\">\n";
            for (const scalar v : *vals) out << v << '\n';
            out << "</DataArray>\n";
        }
        for (const auto& [name, vals] : vectors_)
        {
            out << "<DataArray type=\"Float64\" Name=\"" << name
                << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
            for (const Vector& v : *vals) out << v.x << ' ' << v.y << ' ' << v.z << '\n';
            out << "</DataArray>\n";
        }
        out << "</CellData>\n";
    }

    out << "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
    log::info("vtk: wrote ", path);
}

} // namespace cf
