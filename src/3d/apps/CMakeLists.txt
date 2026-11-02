add_executable(chflow-scalar-transport scalarTransportApp.cpp)
target_link_libraries(chflow-scalar-transport PRIVATE cf3d)

add_executable(chflow-simple-foam simpleFoamApp.cpp)
target_link_libraries(chflow-simple-foam PRIVATE cf3d)

add_executable(chflow-piso-foam pisoFoamApp.cpp)
target_link_libraries(chflow-piso-foam PRIVATE cf3d)
