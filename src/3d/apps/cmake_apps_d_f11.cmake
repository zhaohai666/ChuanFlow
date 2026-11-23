# Batch D item F11 app: moving-mesh relative flux + rigid 6DOF driver.
# Included by the aggregate hook cmake_apps_3dphys.cmake.
add_executable(chflow-3d-dynamic-mesh dynamicMesh6DofApp.cpp)
target_link_libraries(chflow-3d-dynamic-mesh PRIVATE cf3d)
