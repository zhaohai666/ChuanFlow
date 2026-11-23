# Batch D / F10: standalone 0D detailed-chemistry reactor driver.
add_executable(chflow-3d-combustion-detailed combustionDetailedApp.cpp)
target_link_libraries(chflow-3d-combustion-detailed PRIVATE cf3d)
