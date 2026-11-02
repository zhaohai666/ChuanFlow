add_library(cf3d STATIC
    src/mesh/polyMesh.cpp
    src/mesh/blockMesh.cpp
    src/fields/volFields.cpp
    src/discrete/fvMatrix.cpp
    src/discrete/interpolate.cpp
    src/discrete/divScheme.cpp
    src/discrete/laplacianScheme.cpp
    src/discrete/gradScheme.cpp
    src/io/vtkWrite.cpp
    src/solvers/pressureCorrection.cpp
    src/solvers/simpleSolver.cpp
    src/solvers/pisoSolver.cpp
)

target_include_directories(cf3d PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(cf3d PUBLIC cf_common)

add_subdirectory(apps)
