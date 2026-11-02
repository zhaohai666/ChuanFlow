add_executable(chflow-1d-steady steadyNetworkApp.cpp)
target_link_libraries(chflow-1d-steady PRIVATE cf1d)

add_executable(chflow-1d-transient transientApp.cpp)
target_link_libraries(chflow-1d-transient PRIVATE cf1d)
