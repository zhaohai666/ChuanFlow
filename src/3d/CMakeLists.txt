# Per-batch source hooks (parallel development): each hook may set its own
# <NAME>_SRCS list; missing hooks are no-ops.
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake_src_3dturb.cmake OPTIONAL)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake_src_3dphys.cmake OPTIONAL)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake_src_3dmesh.cmake OPTIONAL)

add_library(cf3d STATIC
    src/mesh/polyMesh.cpp
    src/mesh/blockMesh.cpp
    src/mesh/gmshMesh.cpp
    src/mesh/meshLoader.cpp
    src/mesh/openFoamMesh.cpp
    src/parallel/domainDecomposition.cpp
    src/parallel/domainOperators.cpp
    src/parallel/backend.cpp
    src/fields/volFields.cpp
    src/discrete/fvMatrix.cpp
    src/discrete/interpolate.cpp
    src/discrete/divScheme.cpp
    src/discrete/laplacianScheme.cpp
    src/discrete/gradScheme.cpp
    src/discrete/reconstruction.cpp
    src/io/vtkWrite.cpp
    src/io/hdf5Write.cpp
    src/turbulence/wallFunction.cpp
    src/turbulence/turbulenceModel.cpp
    src/turbulence/kEpsilon.cpp
    src/turbulence/kOmegaSST.cpp
    src/solvers/pressureCorrection.cpp
    src/solvers/compressibleFlow.cpp
    src/solvers/simpleSolver.cpp
    src/solvers/generalizedNewtonian.cpp
    src/solvers/pisoSolver.cpp
    src/solvers/energySolver.cpp
    src/dpm/dpm.cpp
    src/dpm/radiation.cpp
    src/dpm/combustion.cpp
    src/dpm/qmom.cpp
    src/dpm/cht.cpp
    src/dpm/ale.cpp
    src/multiphase/vof.cpp
    src/multiphase/plic.cpp
    src/multiphase/plicAdvect.cpp
    ${CF3D_TURB_SRCS}
    ${CF3D_PHYS_SRCS}
    ${CF3D_MESH_SRCS}
)

target_include_directories(cf3d PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(cf3d PUBLIC cf_common)

# Optional parallel backends (M6): both default OFF, so the default build is
# the bit-identical serial CPU backend. Kokkos needs the external SDK; a
# missing SDK downgrades to the CPU backend rather than failing the build.
if(CF_BACKEND_KOKKOS)
    target_compile_definitions(cf3d PUBLIC CF_WITH_KOKKOS)
    target_link_libraries(cf3d PUBLIC Kokkos::kokkos)
elseif(CF_BACKEND_OPENMP)
    target_compile_definitions(cf3d PUBLIC CF_WITH_OPENMP)
    target_link_libraries(cf3d PUBLIC OpenMP::OpenMP_CXX)
endif()

# Optional HDF5 result writer (README IO): silently skipped without the
# system library; hdf5Write.cpp always compiles (fallback branch cfFatals).
find_package(HDF5 QUIET COMPONENTS C)
if(HDF5_FOUND)
    target_compile_definitions(cf3d PUBLIC CF_HAVE_HDF5)
    target_include_directories(cf3d PRIVATE ${HDF5_INCLUDE_DIRS})
    target_link_libraries(cf3d PUBLIC ${HDF5_C_LIBRARIES})
endif()

add_subdirectory(apps)
