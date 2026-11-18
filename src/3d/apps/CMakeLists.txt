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
target_link_libraries(chflow-vof PRIVATE cf3d)

add_executable(chflow-3d-compressible compressibleApp.cpp)
target_link_libraries(chflow-3d-compressible PRIVATE cf3d)

add_executable(chflow-3d-parallel parallelApp.cpp)
target_link_libraries(chflow-3d-parallel PRIVATE cf3d)


add_executable(chflow-openfoam-mesh openFoamMeshApp.cpp)
target_link_libraries(chflow-openfoam-mesh PRIVATE cf3d)

add_executable(chflow-3d-backend backendApp.cpp)
target_link_libraries(chflow-3d-backend PRIVATE cf3d)
