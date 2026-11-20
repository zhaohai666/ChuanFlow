add_executable(chflow-1d-steady steadyNetworkApp.cpp)
target_link_libraries(chflow-1d-steady PRIVATE cf1d)

add_executable(chflow-1d-transient transientApp.cpp)
target_link_libraries(chflow-1d-transient PRIVATE cf1d)

add_executable(chflow-1d-twofluid twoFluidApp.cpp)
target_link_libraries(chflow-1d-twofluid PRIVATE cf1d)

add_executable(chflow-1d-nodal nodalApp.cpp)
target_link_libraries(chflow-1d-nodal PRIVATE cf1d)

add_executable(chflow-1d-point-elements pointElementsApp.cpp)
target_link_libraries(chflow-1d-point-elements PRIVATE cf1d)

# Per-batch app hooks (parallel development); missing hooks are no-ops.
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake_apps_1dsteady.cmake OPTIONAL)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake_apps_1dtransient.cmake OPTIONAL)
