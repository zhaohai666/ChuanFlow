cmake_minimum_required(VERSION 3.24)

project(ChuanFlow
    VERSION 0.1.0
    DESCRIPTION "ChuanFlow: integrated multiphase flow simulation platform (1D network + 1D transient + 3D FVM)"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(NOT MSVC)
    add_compile_options(-Wall -Wextra -Wno-unused-parameter)
endif()

enable_testing()

# Eigen sparse linear solver backend (README M6 / A11). Eigen is header-only
# and vendored under third_party/eigen (3.4.0), so no system package is needed.
# Declared before src/common so cf_common sees the option and include dir.
# Default OFF: the default build keeps the project's own CG/BiCGStab and every
# checker stays bit-identical. Turn ON to dispatch solveWithBackend(...,Eigen)
# to Eigen's ILU-preconditioned BiCGSTAB.
option(CHFLOW_WITH_EIGEN "Use the vendored Eigen sparse solver backend" OFF)
if(CHFLOW_WITH_EIGEN)
    set(CF_EIGEN_DIR "${PROJECT_SOURCE_DIR}/third_party/eigen")
    if(NOT EXISTS "${CF_EIGEN_DIR}/Eigen/Dense")
        message(WARNING "CHFLOW_WITH_EIGEN=ON but third_party/eigen is missing; "
                        "keeping the native solver backend")
        set(CHFLOW_WITH_EIGEN OFF)
    endif()
endif()

# PETSc sparse linear solver backend (README M6 / A11). PETSc is a compiled
# library; discovered via find_package and only enabled when found. A missing
# library downgrades to the native solver with a warning (same QUIET philosophy
# as HDF5 / Kokkos). This keeps the default build free of heavy MPI deps.
option(CHFLOW_WITH_PETSC "Use the PETSc sparse solver backend when available" OFF)
if(CHFLOW_WITH_PETSC)
    find_package(PETSc QUIET)
    if(PETSc_FOUND)
        message(STATUS "PETSc found: ${PETSc_VERSION}")
    else()
        message(WARNING "CHFLOW_WITH_PETSC=ON but no PETSc library found; "
                        "keeping the native solver backend")
        set(CHFLOW_WITH_PETSC OFF)
    endif()
endif()

add_subdirectory(src/common)
add_subdirectory(src/pvt)
# Optional parallel backends (README M6). Both default OFF: the default build
# is the serial CPU backend, which every checker requires to be bit-identical.
# A requested backend that cannot be found downgrades to CPU with a warning
# (same QUIET philosophy as the optional MPI / Python-bindings blocks).
option(CHFLOW_WITH_OPENMP "Use OpenMP in the cf::backend parallelFor" OFF)
option(CHFLOW_WITH_KOKKOS "Use a Kokkos backend (needs the Kokkos SDK)" OFF)
if(CHFLOW_WITH_KOKKOS)
    find_package(Kokkos CONFIG QUIET)
    if(Kokkos_FOUND)
        set(CF_BACKEND_KOKKOS ON)
    else()
        message(WARNING "CHFLOW_WITH_KOKKOS=ON but no Kokkos SDK found; keeping the CPU backend")
    endif()
elseif(CHFLOW_WITH_OPENMP)
    find_package(OpenMP QUIET)
    if(OpenMP_CXX_FOUND)
        set(CF_BACKEND_OPENMP ON)
    else()
        message(WARNING "CHFLOW_WITH_OPENMP=ON but no OpenMP found; keeping the CPU backend")
    endif()
endif()

add_subdirectory(src/1d)
add_subdirectory(src/3d)
add_subdirectory(src/coupling)

# Python bindings (README M1): pybind11 via FetchContent. Skipped when no
# interpreter development headers are available so the C++ build never breaks.
option(CHFLOW_WITH_PYTHON "Build the pybind11 Python bindings" ON)
if(CHFLOW_WITH_PYTHON)
    find_package(Python3 COMPONENTS Interpreter Development.Module QUIET)
    if(Python3_FOUND)
        # use the modern FindPython3 result above (the same interpreter the
        # tests run) instead of pybind11's classic FindPythonLibs fallback,
        # which can resolve a different Python and build an ABI-mismatched .pyd
        set(PYBIND11_FINDPYTHON ON CACHE BOOL "" FORCE)
        set(PYTHON_EXECUTABLE "${Python3_EXECUTABLE}" CACHE FILEPATH "" FORCE)
        set(Python3_EXECUTABLE "${Python3_EXECUTABLE}" CACHE FILEPATH "" FORCE)
        include(FetchContent)
        FetchContent_Declare(pybind11
            GIT_REPOSITORY https://github.com/pybind/pybind11
            GIT_TAG v2.13.6
            GIT_SHALLOW TRUE)
        FetchContent_MakeAvailable(pybind11)
        add_subdirectory(src/bindings)
        set(CHFLOW_PYTHON_BINDINGS ON)
    else()
        message(STATUS "ChuanFlow: Python development headers not found - skipping bindings")
    endif()
endif()

add_subdirectory(tests)
