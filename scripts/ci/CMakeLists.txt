# Standalone Kokkos backend probe (GitHub Actions compile gate + cloud GPU
# execution). Deliberately NOT part of the main build: cf3d's parallelFor
# kernels are host lambdas over host memory and only compile under Kokkos
# host execution spaces; this probe is the device-callable slice of the
# cf::backend::parallelFor contract documented in src/3d/include/cf/parallel/Backend.h.
#
# configure:
#   cmake -S scripts/ci -B <build>                        # host/cpu-serial smoke
#   cmake -S scripts/ci -B <build> -DPROBE_WITH_KOKKOS=ON \
#         [-DKokkos_DIR=<install>/lib/cmake/Kokkos]       \
#         [-DCMAKE_CUDA_ARCHITECTURES=<arch>]             # compile / GPU run

cmake_minimum_required(VERSION 3.24)
project(cf_kokkos_probe LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# The probe TU is compiled by nvcc on CUDA builds; it only needs C++17, and
# CUDA C++20 needs nvcc >= 12 which cloud images do not all carry.
set(CMAKE_CUDA_STANDARD 17)

option(PROBE_WITH_KOKKOS "Drive the probe through a real Kokkos install" OFF)
if(PROBE_WITH_KOKKOS)
    find_package(Kokkos CONFIG REQUIRED)
endif()

add_executable(cf-kokkos-probe
    kokkos_cuda_probe.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../../src/3d/src/parallel/backend.cpp
)
target_include_directories(cf-kokkos-probe PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../src/3d/include)

if(PROBE_WITH_KOKKOS)
    target_compile_definitions(cf-kokkos-probe PRIVATE CF_WITH_KOKKOS)
    target_link_libraries(cf-kokkos-probe PRIVATE Kokkos::kokkos)
    if(Kokkos_ENABLE_CUDA)
        enable_language(CUDA)
        if(NOT CMAKE_CUDA_ARCHITECTURES)
            message(FATAL_ERROR
                "Kokkos was built with CUDA but CMAKE_CUDA_ARCHITECTURES is unset; "
                "native detection fails on GPU-less compile-gate runners")
        endif()
        set_source_files_properties(kokkos_cuda_probe.cpp PROPERTIES LANGUAGE CUDA)
    endif()
endif()

enable_testing()
add_test(NAME verify.kokkosProbeHost
         COMMAND cf-kokkos-probe 1048576)
set_tests_properties(verify.kokkosProbeHost PROPERTIES
    PASS_REGULAR_EXPRESSION "PROBE PASS")
