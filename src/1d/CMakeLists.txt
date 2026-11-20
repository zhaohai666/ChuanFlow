# Per-batch source hooks (parallel development): each hook may set its own
# <NAME>_SRCS list; missing hooks are no-ops.
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake_src_1dsteady.cmake OPTIONAL)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake_src_1dtransient.cmake OPTIONAL)

add_library(cf1d STATIC
    src/network.cpp
    src/transient.cpp
    src/beggsBrill.cpp
    src/gasLiftValve.cpp
    src/rodPump.cpp
    src/jetPump.cpp
    src/progressiveCavityPump.cpp
    src/plungerGasLift.cpp
    src/separator.cpp
    src/compressor.cpp
    src/twoFluid.cpp
    src/nodal.cpp
    ${CF1D_STEADY_SRCS}
    ${CF1D_TRANSIENT_SRCS}
)

target_include_directories(cf1d PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(cf1d PUBLIC cf_common cf_pvt)

add_subdirectory(apps)
