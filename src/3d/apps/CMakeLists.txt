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

# 5.224: the two-fluid Euler-Euler kernel (per-phase velocities, shared
# pressure, interphase drag) and its case-driven app.
add_executable(chflow-3d-euler-euler eulerEulerApp.cpp)
target_link_libraries(chflow-3d-euler-euler PRIVATE cf3d)

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

# Per-batch app hooks (parallel development); missing hooks are no-ops.
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake_apps_3dturb.cmake OPTIONAL)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake_apps_3dphys.cmake OPTIONAL)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake_apps_3dmesh.cmake OPTIONAL)

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

    # 5.225: the full SIMPLE pressure-velocity coupling on Decomposition-split
    # domains -- rank-local per-domain stencils, Sendrecv halos before every
    # stencil application, Allreduce collectives for the global reductions.
    add_executable(chflow-3d-mpi-flow mpiFlowSolverApp.cpp)
    target_compile_definitions(chflow-3d-mpi-flow PRIVATE CF_HAVE_MPI)
    target_include_directories(chflow-3d-mpi-flow PRIVATE
        ${MPI_CXX_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})
    target_link_libraries(chflow-3d-mpi-flow PRIVATE cf3d cf_common MPI::MPI_CXX)

    # A11 cross-check: the SAME graph-Laplacian system solved by PETSc's own
    # MPI machinery (MatCreateAIJ + KSPCG) under a deliberately different
    # partition, dumped in the same six-column format so one checker compares
    # both distributed solvers. Only built when a real PETSc was found.
    if(CHFLOW_WITH_PETSC AND PETSc_LIBRARIES)
        add_executable(chflow-3d-petsc-distributed petscDistributedApp.cpp)
        target_compile_definitions(chflow-3d-petsc-distributed PRIVATE CF_HAVE_MPI)
        target_include_directories(chflow-3d-petsc-distributed PRIVATE
            ${MPI_CXX_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})
        if(PETSc_INCLUDE_DIRS)
            target_include_directories(chflow-3d-petsc-distributed PRIVATE ${PETSc_INCLUDE_DIRS})
        elseif(PETSC_INCLUDE_DIRS)
            target_include_directories(chflow-3d-petsc-distributed PRIVATE ${PETSC_INCLUDE_DIRS})
        endif()
        target_link_libraries(chflow-3d-petsc-distributed PRIVATE
            cf3d cf_common MPI::MPI_CXX ${PETSc_LIBRARIES})
    endif()
endif()
