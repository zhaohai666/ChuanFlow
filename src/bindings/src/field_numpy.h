// Conversion layer contract between the kernel Field<T> storage and NumPy
// ndarrays. Shared by the coupling and 3D binding families.
//
// Contract (pinned by tests/python/test_bindings_3d.py):
//   * dtype: always float64 (cf::scalar is double). A float32 / integer /
//     complex input is refused with TypeError, never silently cast.
//   * layout: returned arrays are C-order (row-major) COPIES owned by Python;
//     they are not views into kernel storage, so mutating the array cannot
//     corrupt a field the solver still owns.
//   * scalar field  -> ndarray shape (n,), C-order.
//   * vector field  -> ndarray shape (n, 3), C-order, components x, y, z.
//   * write-back (numpy -> field): 1-D length n for a scalar field, 2-D
//     shape (n, 3) for a vector field; any other rank or an element count
//     that does not match the field raises ValueError.
#pragma once

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <stdexcept>

#include "cf/common/Field.h"
#include "cf/common/Vector.h"

namespace py = pybind11;

namespace cfbind
{

inline py::array_t<double> scalarFieldToNumpy(const cf::Field<cf::scalar>& f)
{
    py::array_t<double> out(static_cast<py::ssize_t>(f.size()));
    std::copy(f.begin(), f.end(), out.mutable_data());
    return out;
}

inline py::array_t<double> vectorFieldToNumpy(const cf::Field<cf::Vector>& f)
{
    py::array_t<double> out({static_cast<py::ssize_t>(f.size()),
                             static_cast<py::ssize_t>(3)});
    auto r = out.mutable_unchecked<2>();
    for (std::size_t i = 0; i < f.size(); ++i)
    {
        r(i, 0) = f[i].x;
        r(i, 1) = f[i].y;
        r(i, 2) = f[i].z;
    }
    return out;
}

namespace detail
{

// Validate dtype / layout up front; pybind11's array_t caster does not reject
// a mismatched dtype on its own, so the float64 / C-order contract is enforced
// here explicitly.
inline const double* checkArray(const py::array& arr, const char* what)
{
    if (!arr.dtype().equal(py::dtype::of<double>()))
        throw py::type_error(std::string(what) + ": dtype must be float64");
    if ((arr.flags() & py::array::c_style) == 0)
        throw py::type_error(std::string(what) + ": array must be C-order");
    return static_cast<const double*>(arr.data());
}

} // namespace detail

inline void scalarFieldFromNumpy(cf::Field<cf::scalar>& f, const py::array& arr)
{
    const double* data = detail::checkArray(arr, "scalar field write-back");
    if (arr.ndim() != 1)
        throw std::invalid_argument("scalar field expects a 1-D array");
    if (static_cast<std::size_t>(arr.size()) != f.size())
        throw std::invalid_argument(
            "scalar field size mismatch: field has " + std::to_string(f.size())
            + " cells, array has " + std::to_string(arr.size()));
    std::copy(data, data + arr.size(), f.begin());
}

inline void vectorFieldFromNumpy(cf::Field<cf::Vector>& f, const py::array& arr)
{
    const double* data = detail::checkArray(arr, "vector field write-back");
    if (arr.ndim() != 2 || arr.shape(1) != 3)
        throw std::invalid_argument(
            "vector field expects a 2-D array of shape (n, 3)");
    if (static_cast<std::size_t>(arr.shape(0)) != f.size())
        throw std::invalid_argument(
            "vector field size mismatch: field has " + std::to_string(f.size())
            + " cells, array has " + std::to_string(arr.shape(0)));
    for (std::size_t i = 0; i < f.size(); ++i)
    {
        f[i].x = data[3 * i + 0];
        f[i].y = data[3 * i + 1];
        f[i].z = data[3 * i + 2];
    }
}

} // namespace cfbind
