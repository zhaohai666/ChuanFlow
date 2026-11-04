find_package(Python3 COMPONENTS Interpreter REQUIRED)

# ---- C++ unit tests ----
add_executable(kernel_tests unit/kernelTests.cpp)
target_link_libraries(kernel_tests PRIVATE cf3d)
add_test(NAME kernel.unit COMMAND kernel_tests)

add_executable(network_tests unit/networkTests.cpp)
target_link_libraries(network_tests PRIVATE cf1d)
target_compile_definitions(network_tests PRIVATE
    CF_CASE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/cases")
add_test(NAME network.unit COMMAND network_tests)

add_executable(transient_tests unit/transientTests.cpp)
target_link_libraries(transient_tests PRIVATE cf1d)
target_compile_definitions(transient_tests PRIVATE
    CF_CASE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/cases")
add_test(NAME transient.unit COMMAND transient_tests)

add_executable(pvt_tests unit/pvtTests.cpp)
target_link_libraries(pvt_tests PRIVATE cf_pvt)
target_compile_definitions(pvt_tests PRIVATE
    CF_CASE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/cases")
add_test(NAME pvt.unit COMMAND pvt_tests)

add_executable(piso_tests unit/pisoTests.cpp)
target_link_libraries(piso_tests PRIVATE cf3d)
add_test(NAME piso.unit COMMAND piso_tests)

add_executable(gmsh_mesh_tests unit/gmshMeshTests.cpp)
target_link_libraries(gmsh_mesh_tests PRIVATE cf3d)
target_compile_definitions(gmsh_mesh_tests PRIVATE
    CF_MESH_DIR="${CMAKE_CURRENT_SOURCE_DIR}/meshes")
add_test(NAME gmshMesh.unit COMMAND gmsh_mesh_tests)

add_executable(turbulence_tests unit/turbulenceTests.cpp)
target_link_libraries(turbulence_tests PRIVATE cf3d)
add_test(NAME turbulence.unit COMMAND turbulence_tests)

add_executable(vof_tests unit/vofTests.cpp)
target_link_libraries(vof_tests PRIVATE cf3d)
add_test(NAME vof.unit COMMAND vof_tests)

add_executable(plic_tests unit/plicTests.cpp)
target_link_libraries(plic_tests PRIVATE cf3d)
target_compile_definitions(plic_tests PRIVATE
    CF_MESH_DIR="${CMAKE_CURRENT_SOURCE_DIR}/meshes")
add_test(NAME plic.unit COMMAND plic_tests)

# ---- analytic verification cases ----
add_test(
    NAME verify.advectionDiffusion
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_scalar_profile.py
        $<TARGET_FILE:chflow-scalar-transport>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/advectionDiffusion/case.dict
        advection 5e-2)

add_test(
    NAME verify.pureDiffusion
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_scalar_profile.py
        $<TARGET_FILE:chflow-scalar-transport>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/pureDiffusion/case.dict
        linearZ 1e-8)

add_test(
    NAME verify.gmshDiffusion
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_scalar_profile.py
        $<TARGET_FILE:chflow-scalar-transport>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/gmshDiffusion/case.dict
        linearZ 1e-10 out_gmsh_poisson.vtu)

add_test(
    NAME verify.poiseuille
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_poiseuille.py
        $<TARGET_FILE:chflow-simple-foam>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/poiseuilleChannel/case.dict)

add_test(
    NAME verify.startedChannel
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_piso.py
        $<TARGET_FILE:chflow-piso-foam>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/startedChannel/case.dict)

add_test(
    NAME verify.turbulentChannel
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_turbulent_channel.py
        $<TARGET_FILE:chflow-piso-foam>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/turbulentChannel/case.dict)

# the same friction physics through the steady coupling; the case carries both a
# simple{} and a piso{} block and each app reads its own
add_test(
    NAME verify.turbulentChannelSimple
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_turbulent_channel.py
        $<TARGET_FILE:chflow-simple-foam>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/turbulentChannel/case.dict)

add_test(
    NAME verify.vofTranslation
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_vof.py
        $<TARGET_FILE:chflow-vof>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/vofTranslation/case.dict)

# the interface geometry the transport has to preserve: solid rotation is rigid,
# so the droplet keeps its radius, its centre turns at omega, and its interface
# keeps reading -1/R
add_test(
    NAME verify.vofRotation
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_vof.py
        $<TARGET_FILE:chflow-vof>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/vofRotation/case.dict)

add_test(
    NAME verify.vofDroplet
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_vof.py
        $<TARGET_FILE:chflow-vof>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/vofDroplet/case.dict)

# the same droplet, cell-sharp, with the surface force reading its curvature from
# the reconstructed planes: the jump the pressure solves to is still sigma/(rho R),
# on a field whose own slope cannot say what the curvature is
add_test(
    NAME verify.vofDropletPlic
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_vof.py
        $<TARGET_FILE:chflow-vof>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/vofDropletPlic/case.dict)

# the two-way loop: the solved flux carries alpha and the force is rebuilt from
# the shape it leaves, so an ellipse oscillates about the circle; the crossing
# clock decides the density-matched Lamb period against the free-drop one
add_test(
    NAME verify.vofOscillating
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_vof.py
        $<TARGET_FILE:chflow-vof>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/vofOscillating/case.dict)

# the same droplet cell-sharp, with the interface reconstructed as a plane in every
# band cell: the band, the fraction each plane reads back, the area of a one-cell
# band and the volume of the exact filling all have closed forms
add_test(
    NAME verify.vofPlic
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_vof.py
        $<TARGET_FILE:chflow-vof>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/vofPlic/case.dict)

# the same droplet carried by its own planes instead of by a flux of alpha: the
# phase volume is conserved by construction, the front stays one cell thick over a
# full turn, and the field is compared cell by cell with the analytic rotated disk
add_test(
    NAME verify.vofZalesak
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_vof.py
        $<TARGET_FILE:chflow-vof>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/vofZalesak/case.dict)

# the same turn and the same planes, with each cell's plane found by the least-squares
# fit of its neighbourhood: the direction converges with the mesh, and the turn leaves a
# narrower band, a shorter interface and a shape nearer the disk than the cell-local one
add_test(
    NAME verify.vofFit
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_vof.py
        $<TARGET_FILE:chflow-vof>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/vofFit/case.dict)

# 3D energy equation. The three conduction cases have no flow at all, so they
# isolate the diffusion operator, the temperature conditions (fixedValue,
# heatFlux, adiabatic) and the heat-balance report; the decay case runs the same
# physics through the transient coupling against a discrete eigenmode; the rest
# case is the Boussinesq force a uniform temperature offset implies, and the
# channel case is where convection, wall flux and the continuity flux have to
# agree - its Nusselt number is the number the kernel is not given.
add_test(
    NAME verify.energyConductionFixed
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_energy.py
        $<TARGET_FILE:chflow-simple-foam>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/energyConductionFixed/case.dict
        fixed)

add_test(
    NAME verify.energyConductionFlux
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_energy.py
        $<TARGET_FILE:chflow-simple-foam>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/energyConductionFlux/case.dict
        flux)

add_test(
    NAME verify.energyConductionSource
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_energy.py
        $<TARGET_FILE:chflow-simple-foam>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/energyConductionSource/case.dict
        source)

add_test(
    NAME verify.energyDecay
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_energy.py
        $<TARGET_FILE:chflow-piso-foam>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/energyDecay/case.dict
        decay)

add_test(
    NAME verify.energyBuoyantRest
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_energy.py
        $<TARGET_FILE:chflow-piso-foam>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/energyBuoyantRest/case.dict
        rest)

add_test(
    NAME verify.energyChannelFlux
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_energy.py
        $<TARGET_FILE:chflow-simple-foam>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/energyChannelFlux/case.dict
        channel)

add_test(
    NAME verify.steadyNetwork
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_steady_network.py
        $<TARGET_FILE:chflow-1d-steady>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/steadyNetwork/case.dict)

# black-oil Beggs-Brill gathering network: balance, holdup/regime sanity and
# a dry-oil GOR sensitivity twin
add_test(
    NAME verify.steadyNetworkMultiphase
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_steady_network_multiphase.py
        $<TARGET_FILE:chflow-1d-steady>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/steadyNetworkMultiphase/case.dict)

foreach(tmode "blowdown;BLOWDOWN" "rlc;RLC" "slug;Slug")
    list(GET tmode 0 mode)
    list(GET tmode 1 Mode)
    add_test(
        NAME verify.transient${Mode}
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/check_transient.py
            $<TARGET_FILE:chflow-1d-transient>
            ${CMAKE_CURRENT_SOURCE_DIR}/cases/transient${Mode}/case.dict
            ${mode})
endforeach()

# distributed two-fluid kernel: Godunov shock/fan kinematics, steady
# drift-flux holdup (Wallis and Zuber closures), sealed-pipe gravity
# segregation, laminar pressure-driven start-up
foreach(tcase twoFluidFront twoFluidSteadyHoldupWallis
              twoFluidSteadyHoldupZuber twoFluidSettling
              twoFluidPressureLaminar)
    string(REGEX REPLACE "^twoFluid" "" mode ${tcase})
    string(TOLOWER ${mode} mode)
    add_test(
        NAME verify.${tcase}
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/check_two_fluid.py
            $<TARGET_FILE:chflow-1d-twofluid>
            ${CMAKE_CURRENT_SOURCE_DIR}/cases/${tcase}/case.dict
            ${mode})
endforeach()

# The PVT case sweeps a pressure envelope through the flash; the checker
# asserts thermodynamic identities (material balance, fugacity equality,
# envelope collapse, table vs direct EOS) rather than stored numbers.
add_test(
    NAME verify.pvtFlash
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_pvt.py
        $<TARGET_FILE:chflow-pvt>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/pvtFlash/case.dict)

if(CHFLOW_PYTHON_BINDINGS)
    add_test(
        NAME python.bindings
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/python/test_bindings.py)
    # ':' not ';': CMake would split an ENVIRONMENT value containing ';'
    set_tests_properties(python.bindings PROPERTIES
        ENVIRONMENT "PYTHONPATH=${PROJECT_SOURCE_DIR}/python:$<TARGET_FILE_DIR:chflow_core>")
endif()

# ---- CLI (Rust) ----
find_program(CARGO_EXECUTABLE cargo)
if(CARGO_EXECUTABLE)
    add_test(NAME cli.unit
        COMMAND ${CARGO_EXECUTABLE} test --offline --quiet
            --manifest-path ${PROJECT_SOURCE_DIR}/cli/Cargo.toml
            --target-dir ${CMAKE_BINARY_DIR}/cli)
    add_test(NAME cli.runSteady
        COMMAND ${CARGO_EXECUTABLE} run --offline --quiet
            --manifest-path ${PROJECT_SOURCE_DIR}/cli/Cargo.toml
            --target-dir ${CMAKE_BINARY_DIR}/cli
            -- run steady ${CMAKE_CURRENT_SOURCE_DIR}/cases/steadyNetwork/case.dict)
    set_tests_properties(cli.runSteady PROPERTIES
        ENVIRONMENT "CHUANFLOW_BIN_DIR=$<TARGET_FILE_DIR:chflow-1d-steady>")
    add_test(NAME cli.runTransient
        COMMAND ${CARGO_EXECUTABLE} run --offline --quiet
            --manifest-path ${PROJECT_SOURCE_DIR}/cli/Cargo.toml
            --target-dir ${CMAKE_BINARY_DIR}/cli
            -- run transient ${CMAKE_CURRENT_SOURCE_DIR}/cases/transientBlowdown/case.dict)
    set_tests_properties(cli.runTransient PROPERTIES
        ENVIRONMENT "CHUANFLOW_BIN_DIR=$<TARGET_FILE_DIR:chflow-1d-transient>")
    add_test(NAME cli.runPvt
        COMMAND ${CARGO_EXECUTABLE} run --offline --quiet
            --manifest-path ${PROJECT_SOURCE_DIR}/cli/Cargo.toml
            --target-dir ${CMAKE_BINARY_DIR}/cli
            -- run pvt ${CMAKE_CURRENT_SOURCE_DIR}/cases/pvtFlash/case.dict)
    set_tests_properties(cli.runPvt PROPERTIES
        ENVIRONMENT "CHUANFLOW_BIN_DIR=$<TARGET_FILE_DIR:chflow-pvt>")
endif()
