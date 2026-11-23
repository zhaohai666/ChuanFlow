# Batch wt-3dturb (F4 high-order turbulence + F12 phase change): sources added
# to the cf3d library. Missing from the main list on purpose - the main
# CMakeLists includes this hook OPTIONAL, so the default build is unchanged
# until this file exists.
set(CF3D_TURB_SRCS
    src/turbulence/lesModel.cpp
    src/turbulence/rsmModel.cpp
    src/turbulence/desModel.cpp
    src/turbulence/transitionModel.cpp
    src/solvers/cavitationModel.cpp
)
