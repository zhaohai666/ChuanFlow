# Batch wt-3dturb apps (F4 + F12). Each app reads its own opt-in case block,
# so the existing piso/simple apps stay untouched.
add_executable(chflow-3d-les lesApp.cpp)
target_link_libraries(chflow-3d-les PRIVATE cf3d)

add_executable(chflow-3d-rsm rsmApp.cpp)
target_link_libraries(chflow-3d-rsm PRIVATE cf3d)

add_executable(chflow-3d-melting meltingApp.cpp)
target_link_libraries(chflow-3d-melting PRIVATE cf3d)

add_executable(chflow-3d-cavitation cavitationApp.cpp)
target_link_libraries(chflow-3d-cavitation PRIVATE cf3d)
