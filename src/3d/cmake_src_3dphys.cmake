# Batch D aggregate hook: each physics item owns a uniquely named fragment, so
# parallel batches never edit the same file. Missing fragments are no-ops.
# The main CMakeLists consumes CF3D_PHYS_SRCS after this include.
set(CF3D_PHYS_SRCS "")
foreach(_d_fragment d_f7 d_f9 d_f10 d_f11 d_f13)
    include(${CMAKE_CURRENT_SOURCE_DIR}/cmake_src_${_d_fragment}.cmake OPTIONAL)
endforeach()
