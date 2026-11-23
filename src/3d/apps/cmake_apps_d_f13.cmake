# Batch D / F13: standalone adaptive-refinement verifier.
add_executable(chflow-3d-adaptive-refine adaptiveRefineApp.cpp)
target_link_libraries(chflow-3d-adaptive-refine PRIVATE cf3d)
