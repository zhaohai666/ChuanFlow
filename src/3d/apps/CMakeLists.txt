add_executable(chflow-scalar-transport scalarTransportApp.cpp)
target_link_libraries(chflow-scalar-transport PRIVATE cf3d)

add_executable(chflow-simple-foam simpleFoamApp.cpp)
target_link_libraries(chflow-simple-foam PRIVATE cf3d)

add_executable(chflow-piso-foam pisoFoamApp.cpp)
target_link_libraries(chflow-piso-foam PRIVATE cf3d)

add_executable(chflow-vof vofApp.cpp)
target_link_libraries(chflow-vof PRIVATE cf3d)

add_executable(chflow-3d-compressible compressibleApp.cpp)
target_link_libraries(chflow-3d-compressible PRIVATE cf3d)

add_executable(chflow-3d-parallel parallelApp.cpp)
target_link_libraries(chflow-3d-parallel PRIVATE cf3d)
