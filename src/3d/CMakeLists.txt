add_library(cf3d STATIC
    src/mesh/polyMesh.cpp
    src/mesh/blockMesh.cpp
    src/mesh/gmshMesh.cpp
    src/mesh/meshLoader.cpp
    src/fields/volFields.cpp
    src/discrete/fvMatrix.cpp
    src/discrete/interpolate.cpp
    src/discrete/divScheme.cpp
    src/discrete/laplacianScheme.cpp
    src/discrete/gradScheme.cpp
    src/io/vtkWrite.cpp
    src/turbulence/wallFunction.cpp
    src/turbulence/turbulenceModel.cpp
    src/turbulence/kEpsilon.cpp
    src/turbulence/kOmegaSST.cpp
    src/solvers/pressureCorrection.cpp
    src/solvers/simpleSolver.cpp
    src/solvers/pisoSolver.cpp
    src/solvers/energySolver.cpp
    src/multiphase/vof.cpp
    src/multiphase/plic.cpp
)

target_include_directories(cf3d PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(cf3d PUBLIC cf_common)

add_subdirectory(apps)
