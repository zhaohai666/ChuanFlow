# Batch D item F7 app hook: opt-in executable. Included by the aggregate
# cmake_apps_3dphys.cmake; a no-op when this fragment is absent.
add_executable(chflow-3d-mixture-pbm mixturePbmApp.cpp)
target_link_libraries(chflow-3d-mixture-pbm PRIVATE cf3d)
