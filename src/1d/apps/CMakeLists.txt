add_executable(chflow-1d-steady steadyNetworkApp.cpp)
target_link_libraries(chflow-1d-steady PRIVATE cf1d)

add_executable(chflow-1d-transient transientApp.cpp)
target_link_libraries(chflow-1d-transient PRIVATE cf1d)

add_executable(chflow-1d-twofluid twoFluidApp.cpp)
target_link_libraries(chflow-1d-twofluid PRIVATE cf1d)

add_executable(chflow-1d-nodal nodalApp.cpp)
target_link_libraries(chflow-1d-nodal PRIVATE cf1d)
