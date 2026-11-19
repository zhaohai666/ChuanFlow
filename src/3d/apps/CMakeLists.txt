add_executable(chflow-scalar-transport scalarTransportApp.cpp)
target_link_libraries(chflow-scalar-transport PRIVATE cf3d)

add_executable(chflow-simple-foam simpleFoamApp.cpp)
target_link_libraries(chflow-simple-foam PRIVATE cf3d)

add_executable(chflow-3d-dpm dpmApp.cpp)
target_link_libraries(chflow-3d-dpm PRIVATE cf3d)

add_executable(chflow-3d-p1 p1App.cpp)
target_link_libraries(chflow-3d-p1 PRIVATE cf3d)

add_executable(chflow-3d-combustion combustionApp.cpp)
target_link_libraries(chflow-3d-combustion PRIVATE cf3d)

add_executable(chflow-3d-qmom qmomApp.cpp)
target_link_libraries(chflow-3d-qmom PRIVATE cf3d)

add_executable(chflow-3d-cht chtApp.cpp)
target_link_libraries(chflow-3d-cht PRIVATE cf3d)

add_executable(chflow-3d-ale aleApp.cpp)
target_link_libraries(chflow-3d-ale PRIVATE cf3d)

add_executable(chflow-piso-foam pisoFoamApp.cpp)
target_link_libraries(chflow-piso-foam PRIVATE cf3d)

add_executable(chflow-vof vofApp.cpp)
target_link_libraries(chflow-vof PRIVATE cf3d cf_pvt)

add_executable(chflow-3d-compressible compressibleApp.cpp)
target_link_libraries(chflow-3d-compressible PRIVATE cf3d)

# MUSCL limited reconstruction standalone verifier (5.130): replays the
# reconstruction operator face by face on analytical profiles so the checker
# can compare against closed forms and the limiter functions independently.
add_executable(chflow-reconstruction reconstructionApp.cpp)
target_link_libraries(chflow-reconstruction PRIVATE cf3d)

add_executable(chflow-3d-parallel parallelApp.cpp)
target_link_libraries(chflow-3d-parallel PRIVATE cf3d)


add_executable(chflow-openfoam-mesh openFoamMeshApp.cpp)
target_link_libraries(chflow-openfoam-mesh PRIVATE cf3d)

add_executable(chflow-3d-backend backendApp.cpp)
target_link_libraries(chflow-3d-backend PRIVATE cf3d)

# README M3 / A11: real MPI distributed linear solve (halo transport + CSR
# CG with Allreduce dot products). Built only when an MPI C++ toolchain is
# available; cf3d / cf_common stay MPI-free and the default (serial) build is
# unchanged. The test runs ranks 1/2/4 and compares the distributed solution
# with the serial solveCG on the same global CSR.
find_package(MPI QUIET COMPONENTS CXX)
if(MPI_CXX_FOUND)
    add_executable(chflow-mpi-linear-solver mpiLinearSolverApp.cpp)
    target_compile_definitions(chflow-mpi-linear-solver PRIVATE CF_HAVE_MPI)
    target_include_directories(chflow-mpi-linear-solver PRIVATE
        ${MPI_CXX_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})
    target_link_libraries(chflow-mpi-linear-solver PRIVATE cf3d cf_common MPI::MPI_CXX)
endif()
