# F9 fragment: one binary, model do|s2s|p1. Included by cmake_apps_3dphys.cmake.
add_executable(chflow-3d-radiation-do radiationDoS2sApp.cpp)
target_link_libraries(chflow-3d-radiation-do PRIVATE cf3d)
