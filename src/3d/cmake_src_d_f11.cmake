# Batch D item F11: moving-mesh relative-flux transport + rigid 6DOF kernel.
# Included by the aggregate hook cmake_src_3dphys.cmake; CF3D_PHYS_SRCS is
# consumed by src/3d/CMakeLists.txt after the hook.
list(APPEND CF3D_PHYS_SRCS src/dpm/dynamicMesh6Dof.cpp)
