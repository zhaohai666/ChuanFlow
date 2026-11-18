#include "cf/fvm/Hdf5Writer.h"

#include <vector>

#include "cf/common/Error.h"

#ifdef CF_HAVE_HDF5
#include <hdf5.h>
#endif

namespace cf
{

#ifdef CF_HAVE_HDF5
namespace
{

void writeDataset1D(hid_t loc, const char* name, hid_t dtype,
                    const void* data, std::size_t n)
{
    const hsize_t dim = static_cast<hsize_t>(n);
    const hid_t space = H5Screate_simple(1, &dim, nullptr);
    const hid_t ds = H5Dcreate2(loc, name, dtype, space,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (ds < 0) cfFatal(std::string("hdf5: cannot create dataset '") + name + "'");
    if (H5Dwrite(ds, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data) < 0)
        cfFatal(std::string("hdf5: cannot write dataset '") + name + "'");
    H5Dclose(ds);
    H5Sclose(space);
}

} // namespace
#endif

void writeFieldsHdf5(const std::string& path,
                     const polyMesh& mesh,
                     const std::map<std::string, Field<scalar>>& fields)
{
#ifndef CF_HAVE_HDF5
    (void)path;
    (void)mesh;
    (void)fields;
    cfFatal("HDF5 support not compiled: rebuild with CF_HAVE_HDF5 and the HDF5 C library");
#else
    for (const auto& [name, values] : fields)
    {
        if (static_cast<label>(values.size()) != mesh.nCells())
            cfFatal("hdf5: field '" + name + "' size mismatch (expected "
                    + std::to_string(mesh.nCells()) + " cells)");
    }

    const hid_t file = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) cfFatal("hdf5: cannot open '" + path + "' for writing");

    // points, flattened as (x0 y0 z0 x1 y1 z1 ...)
    std::vector<double> coords(static_cast<std::size_t>(mesh.nPoints()) * 3);
    for (label p = 0; p < mesh.nPoints(); ++p)
    {
        const Vector& v = mesh.points()[static_cast<std::size_t>(p)];
        coords[static_cast<std::size_t>(p) * 3 + 0] = v.x;
        coords[static_cast<std::size_t>(p) * 3 + 1] = v.y;
        coords[static_cast<std::size_t>(p) * 3 + 2] = v.z;
    }
    const hsize_t pointDims[2] = {static_cast<hsize_t>(mesh.nPoints()), 3};
    const hid_t pointSpace = H5Screate_simple(2, pointDims, nullptr);
    const hid_t pointDs = H5Dcreate2(file, "points", H5T_NATIVE_DOUBLE, pointSpace,
                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (pointDs >= 0)
    {
        H5Dwrite(pointDs, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, coords.data());
        H5Dclose(pointDs);
    }
    H5Sclose(pointSpace);

    writeDataset1D(file, "owner", H5T_NATIVE_INT64,
                   mesh.owner().data(), mesh.owner().size());
    writeDataset1D(file, "neighbour", H5T_NATIVE_INT64,
                   mesh.neighbour().data(), mesh.neighbour().size());

    const hid_t group = H5Gcreate2(file, "fields", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (group < 0) cfFatal("hdf5: cannot create /fields group");
    for (const auto& [name, values] : fields)
    {
        writeDataset1D(group, name.c_str(), H5T_NATIVE_DOUBLE,
                       values.data(), values.size());
    }
    H5Gclose(group);
    H5Fclose(file);
#endif
}

} // namespace cf
