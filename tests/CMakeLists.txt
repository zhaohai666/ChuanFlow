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

add_executable(nodal_tests unit/nodalTests.cpp)
target_link_libraries(nodal_tests PRIVATE cf1d)
target_compile_definitions(nodal_tests PRIVATE
    CF_CASE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/cases")
add_test(NAME nodal.unit COMMAND nodal_tests)

# nodal analysis: a Vogel inflow curve crossed with the Beggs-Brill lift curve;
# the living gas well, its dead-oil twin (gas lifts the well) and the dead well
add_test(
    NAME verify.nodalAnalysis
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_nodal.py
        $<TARGET_FILE:chflow-1d-nodal>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/nodalVogel/case.dict
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/nodalDeadOil/case.dict
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/nodalDeadWell/case.dict)

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

# the partitioned 1D<->3D interface loop against analytic partner models: the
# linear case must need one coefficient refresh, the nonlinear one must land on
# the exact root of dp(Q) = dp_network(Q), and every number a caller reads must
# describe one state (rate the partner was asked for == rate reported)
add_executable(coupling_tests unit/couplingTests.cpp)
target_link_libraries(coupling_tests PRIVATE cfcoupling)
add_test(NAME coupling.unit COMMAND coupling_tests)

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

# the thermal wall function, on the same force-driven column with 310 K / 300 K
# walls: the checker restates the discrete conduction problem (harmonic
# diffusivity inside, the two-branch T+ law at the walls) in Python and compares
# the field and the reported fluxes against it, for both turbulence models
add_test(
    NAME verify.thermalWallChannel
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_thermal_wall.py
        $<TARGET_FILE:chflow-simple-foam>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/thermalWallChannel/case.dict)

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

# density-based Euler, first increment: the Sod shock tube on the polyhedral
# kernel. The checker re-runs the HLLC recursion cell by cell (1e-9), holds the
# mass/energy totals and the wall-impulse momentum balance exact, and compares
# the plateaus with the exact Riemann solution whose star state it solves itself
add_test(
    NAME verify.sodShockTube
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_euler.py
        $<TARGET_FILE:chflow-3d-compressible>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/sodShockTube/case.dict)

# 3D domain decomposition, first increment: the partition, the halo exchange and
# the global reductions, verified with no MPI in sight. The checker re-derives
# the coordinate-median ordering and every domain's halo/face classes from the
# mesh geometry in the .vtu, holds the assembled gradient and net out-flux to the
# serial operators exactly, and requires the written fields to be identical for
# 1, 2, 3, 5, 7, 8 and 16 domains (hex) and 1, 3, 7 and 13 (tets).
add_test(
    NAME verify.parallelDecomposition
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_parallel.py
        $<TARGET_FILE:chflow-3d-parallel>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/parallelDecomposition/case.dict)

add_test(
    NAME verify.parallelDecompositionTet
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_parallel.py
        $<TARGET_FILE:chflow-3d-parallel>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/parallelDecomposition/caseTet.dict)

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

# closed box + uniform gravity: the patch-uniform p BCs cannot hold the linear
# hydrostatic field (the same box without pRgh creeps to uMax 1.83 m/s and
# grows); p_rgh absorbs the force exactly and the box stays at rest
add_test(
    NAME verify.pRghClosedBox
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_prgh.py
        $<TARGET_FILE:chflow-piso-foam>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/closedBoxGravity/case.dict)

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

# inline equipment on laminar single-phase pipes: the pipe law reduces to an
# explicit quadratic (valve loss K, quadratic pump curve), so the checker
# never iterates; the network mode is also a pump-head sign test
foreach(tcase steadyEquipmentValve steadyEquipmentPump steadyEquipmentNetwork)
    string(REGEX REPLACE "^steadyEquipment" "" mode ${tcase})
    string(TOLOWER ${mode} mode)
    add_test(
        NAME verify.${tcase}
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/check_steady_equipment.py
            $<TARGET_FILE:chflow-1d-steady>
            ${CMAKE_CURRENT_SOURCE_DIR}/cases/${tcase}/case.dict
            ${mode})
endforeach()

# non-Newtonian pipe friction, one substitution deep: the Metzner-Reed
# generalized Reynolds number goes where rho*u*d/mu used to go, and because the
# fluid is then a power law the laminar branch of the same Churchill call stops
# being a correlation and becomes the exact discharge
# Q = (pi d^3/32) [n dp d/((3n+1) K L)]^(1/n). Nine pipes with both ends fixed
# therefore assert against closed forms -- the dp exponent 1/n, the bore
# exponent (3n+1)/n, an unmoved length exponent, the wall stress that no K
# changes, the bit-identity of the declared-but-Newtonian pair, and the two
# parse-time refusals (a flow index without a consistency, a rheology on a
# coupled element)
add_test(
    NAME verify.steadyPowerLawRheology
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_steady_power_law_rheology.py
        $<TARGET_FILE:chflow-1d-steady>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/steadyPowerLawRheology/case.dict)

# black-oil Beggs-Brill gathering network: balance, holdup/regime sanity and
# a dry-oil GOR sensitivity twin
add_test(
    NAME verify.steadyNetworkMultiphase
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_steady_network_multiphase.py
        $<TARGET_FILE:chflow-1d-steady>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/steadyNetworkMultiphase/case.dict)

# tabulated-property single-phase law (fluid { pvt ... }): the checker
# re-parses the case table and re-marches both pipes from the reported
# pressures — compressibility closure + expansion + dp visibility
add_test(
    NAME verify.steadyPvtColumn
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_steady_pvt.py
        $<TARGET_FILE:chflow-1d-steady>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/steadyPvtColumn/case.dict)

# partitioned steady 1D<->3D coupling: the same duct solved in 3D replaces the
# analytic law of edge X, so the checker compares the interface against the
# closed-form network (branch split, pressure chain), the block against planar
# Poiseuille, and the fixed point against three independent starting coefficients
add_test(
    NAME verify.coupled1d3d
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_coupled_1d3d.py
        $<TARGET_FILE:chflow-coupled-1d3d>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/coupled1d3d/case.dict)

# the same duct and the same network, put in motion: the interface is iterated
# inside every time step, so the checker closes each step's momentum row against
# the block's own drop, and regresses the ramp the two solvers produce together
# against the discrete closed form built from the region's inertance
add_test(
    NAME verify.transientCoupled1d3d
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_transient_coupled_1d3d.py
        $<TARGET_FILE:chflow-coupled-1d3d-transient>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/transientCoupled1d3d/case.dict)

# two detailed regions on one network: a parallel duct whose rate the split
# decides and a series duct whose rate topology pins, so the checker reads the
# interface per edge, the network laws off the 1D side alone, each block's field
# against its own duct, and the whole thing against a series-parallel closed form
add_test(
    NAME verify.coupled1d3dMulti
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_coupled_1d3d_multi.py
        $<TARGET_FILE:chflow-coupled-1d3d-multi>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/coupled1d3dMulti/case.dict)

# the same two regions again, but the interface now hands each inlet a SHAPE:
# the checker integrates the profile it asked for, reads the momentum-flux factor
# back off the boundary faces, follows it down the region to the plateau, and
# prices the plug inlet against the Poiseuille one -- three shapes, two edges,
# one run where they differ
add_test(
    NAME verify.coupled1d3dProfileMap
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_coupled_1d3d_profile_map.py
        $<TARGET_FILE:chflow-coupled-1d3d-multi>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/coupled1d3dProfileMap/case.dict)

# the same interface against a partner that is not a duct law: two kOmegaSST
# wall-function ducts, so the checker can watch the region develop to its OWN
# momentum-flux factor rather than the continuum 6/5, read the sign of
# (beta_in - beta_out) off the pressure stations, and price each shape against
# the wall balance 2*rho*u_tau^2/H instead of against Poiseuille
add_test(
    NAME verify.coupled1d3dTurbulentProfile
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_coupled_1d3d_turbulent_profile.py
        $<TARGET_FILE:chflow-coupled-1d3d-multi>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/coupled1d3dTurbulentProfile/case.dict)

# the last of the three axes put together: a SHAPE handed across the interface, a
# wall-function partner that answers with its own beta, and a time level. The
# steady sign law -- beta_in above beta_out means the entry region GIVES back
# static pressure -- is now required at every accepted level, so the recovery is a
# curve in time (pRisePa per level) and not only a station list; and the shape is
# charged against the ramp itself, where a parabolic inlet doubles the inertial
# time constant while the laminar partner of the previous case moved it by 3%
add_test(
    NAME verify.transientCoupledMultiTurbulentProfile
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_transient_coupled_multi_turbulent_profile.py
        $<TARGET_FILE:chflow-coupled-1d3d-transient-multi>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/transientCoupledMultiTurbulentProfile/case.dict)

# the 2x2 the previous case leaves open: the same shaped inlet, but now the block
# itself answers du/dt (schemes{ ddt Euler; } honoured inside the region). What is
# settled here is which half of the interface a shape lives in -- the measured
# inertia is the same to 0.4 pp for plug, parabola and 1/7 while the same three
# shapes move the friction half by 15-51%, so the row's inertance claim is blind to
# the shape and the chord it negotiates is not
add_test(
    NAME verify.transientCoupledMultiTurbulentProfileUnsteady
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_transient_coupled_multi_turbulent_profile_unsteady.py
        $<TARGET_FILE:chflow-coupled-1d3d-transient-multi>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/transientCoupledMultiTurbulentProfileUnsteady/case.dict)

# the boundary the previous case's own per-arm table leaves open: `ddts` answers
# for every coupled edge at once, so within one run both edges are asked to step
# or to be quasi-steady together. A per-edge table crosses the formulation axis
# against the edge - one time level holds two different interface readings, one
# edge's resistance charged with a friction half and the other's with its whole
# drop - and the two mirror arms say which edge the mix acts on. What is settled
# here is separability: the edge left alone is the control's computation to
# exactly 0.0 (flat and shaped alike) and the ramp follows its free edge, while
# what the mix changes is only which half the interface stops on, where the
# double count lives, and the price the negotiation ends on
add_test(
    NAME verify.transientCoupledMultiTurbulentEdgeDdt
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_transient_coupled_multi_turbulent_edge_ddt.py
        $<TARGET_FILE:chflow-coupled-1d3d-transient-multi>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/transientCoupledMultiTurbulentEdgeDdt/case.dict)

# the same two regions with BOTH per-edge tables declared at once, on bores typed
# short of their sections so the three inertance sources are three different
# numbers: the row's error then splits into the chord gap the interface never
# iterates plus the bill for whose inertance was installed, edge by edge, and the
# crossing says which of the two the source can remove -- and the grade cannot see
# the second bracket at all on a quasi edge, where it is paid anyway
add_test(
    NAME verify.transientCoupledMultiTurbulentEdgeDdtInertia
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_transient_coupled_multi_turbulent_edge_ddt_inertia.py
        $<TARGET_FILE:chflow-coupled-1d3d-transient-multi>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/transientCoupledMultiTurbulentEdgeDdtInertia/case.dict)

# those two brackets again, now crossed with the inlet shape the interface writes
# on the partner's faces: the bill is arithmetic in the installed M and in dq/dt,
# so a shaped inlet reaches it only through the trajectory -- on the pinned edge,
# where topology gives dq, shaping moves the bill 0.02% and the block's own half
# 21%, while on the free edge the bill and the rate move together by 5% -- and the
# one bracket shape does move is the chord, which stays tolerance-sized either way
add_test(
    NAME verify.transientCoupledMultiTurbulentProfileInertia
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_transient_coupled_multi_turbulent_profile_inertia.py
        $<TARGET_FILE:chflow-coupled-1d3d-transient-multi>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/transientCoupledMultiTurbulentProfileInertia/case.dict)

# the same two regions, inside every time step: the parallel edge ramps on the
# inertial split while the series edge is pinned by topology and only ever
# negotiates its drop, so the checker closes each edge's momentum row with that
# edge's own inertance, counts the block solves against passes times edges, and
# reads the discrete ramp back out of the parallel branch
add_test(
    NAME verify.transientCoupledMulti1d3d
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_transient_coupled_multi.py
        $<TARGET_FILE:chflow-coupled-1d3d-transient-multi>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/transientCoupledMulti1d3d/case.dict)

# the shape question again, now that it is asked inside a time level: five arms
# hand the two edges the same rate in a different shape, and the checker counts
# 60 levels x 2 edges of shape statements, separates what the inlet was given
# from what the region's own planes read, and shows the destination the parallel
# branch relaxes to -- and the tau it relaxes with -- both move with the shape
add_test(
    NAME verify.transientCoupledProfileMap
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_transient_coupled_profile_map.py
        $<TARGET_FILE:chflow-coupled-1d3d-transient-multi>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/transientCoupledProfileMap/case.dict)

# the same network and the same two blocks, with only the source of the region's
# inertia switched: one arm pays rho*L/A from the proxy bore, the other the
# rho*V/Area^2 the block measures on its own mesh. The two edges are mis-stated
# by different amounts, so the checker closes each edge's momentum row with that
# arm's number, shows the interface coefficient is unchanged by the attribution,
# and reads the different first-step split and tau back out of the parallel branch
add_test(
    NAME verify.transientCoupledMultiInertia
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_transient_coupled_multi_inertia.py
        $<TARGET_FILE:chflow-coupled-1d3d-transient-multi>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/transientCoupledMultiInertia/case.dict)

# both edges handed to a nonlinear partner at once, on the schedule that stops on
# the worst of them, one arm plain Gauss-Seidel and the next quasi-Newton. The
# proxy bores are the regions' exact equivalent sections, so §5.26's attribution
# degree of freedom cancels and the only thing left to measure is partner
# nonlinearity: the checker closes each edge's momentum row with that edge's own
# chord and inertance, shows the chord identity survives a RANS partner on both
# edges, and reads the ramp's restoring term off the partners' SLOPES, not chords
add_test(
    NAME verify.transientCoupledMultiTurbulent1d3d
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_transient_coupled_multi_turbulent.py
        $<TARGET_FILE:chflow-coupled-1d3d-transient-multi>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/transientCoupledMultiTurbulent1d3d/case.dict)

# the same two RANS edges, with the attribution degree of freedom the case above
# neutralised put back on the table: both proxy bores are mis-stated, by
# different amounts, so the bore arm and the region arm install genuinely
# different inertances under a partner whose resistance moves with the rate. The
# checker keeps the two accounts apart - the interface stays a chord at either
# inertance, the ramp's time constant is (M_X1+M_By) over the sum of SLOPES so M
# moves tau without touching m, the resistive destination holds no M at all, and
# the pinned series edge is bit-identical across arms while its first-step kick
# is that arm's M*Q/dt
add_test(
    NAME verify.transientCoupledMultiTurbulentInertia
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_transient_coupled_multi_turbulent_inertia.py
        $<TARGET_FILE:chflow-coupled-1d3d-transient-multi>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/transientCoupledMultiTurbulentInertia/case.dict)

# the same two RANS edges and the same two mis-stated bores, with the choice of
# inertance source moved from the interface to the edge: `couple.inertiaFrom`
# gives one arm X1 from its block's mesh and X2 from its drawing's bore inside a
# single run. The checker requires that decision to be exactly separable - the
# mixed arm's X1 states are bit-identical to the region arm's and its X2 pays the
# bore arm's first-step kick - so a per-edge table is a sum of two whole-interface
# runs, not a third reading of the interface
add_test(
    NAME verify.transientCoupledMultiMixedInertia
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_transient_coupled_multi_mixed_inertia.py
        $<TARGET_FILE:chflow-coupled-1d3d-transient-multi>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/transientCoupledMultiMixedInertia/case.dict)

# the two accounts this chain has kept separate - where a region's inertia comes
# from, and which operator negotiates the interface - crossed: a complete 2x2 of
# arms on one topology, both partners RANS, both bores mis-stated. The checker
# requires the axes to be independent in both directions: the operator never
# enters the inertance table, the chord identity and the momentum rows survive
# quasi-Newton unchanged, and within a source pair both operators land on one
# rate, one asymptote and one time constant - so what acceleration buys is the
# schedule, not the physics
add_test(
    NAME verify.transientCoupledMultiTurbulentInertiaNewton
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_transient_coupled_multi_turbulent_inertia_newton.py
        $<TARGET_FILE:chflow-coupled-1d3d-transient-multi>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/transientCoupledMultiTurbulentInertiaNewton/case.dict)

# the second half of that account: the region's changing profile stops being a
# number the block reports and becomes a term the 1D row solves for. The case
# declares it per edge and per arm, so the same run holds the control (nothing
# charged, which must reproduce 5.41's p_net state for state and pass for pass),
# the relabeling control (charged on the pinned edge, where the rate cannot move
# and only the split of one bill changes), and the real statement (charged on the
# free edge, where the row now carries a term quadratic in the rate). The chord
# gives up exactly what the row takes on, so the interface's own residual moves
# with the declaration - the check that nothing is billed twice.
add_test(
    NAME verify.transientCoupledMultiTurbulentConvectiveRow
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_transient_coupled_multi_turbulent_convective_row.py
        $<TARGET_FILE:chflow-coupled-1d3d-transient-multi>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/transientCoupledMultiTurbulentConvectiveRow/case.dict)

# both who-pays tables at once, on one edge: 5.43 moved the region's changing
# profile out of the block's chord and into the row, and recorded that the row's
# other debt - the inertial bill 5.40 decomposed out - was still unpaid because
# the attribution table was held at the drawing everywhere. Here the two tables
# cross on the pinned edge, three rungs by two placements, with the free edge held
# at the drawing and reporting in every arm. Two of the six cells are 5.43's own
# arms, so the factorial is anchored to shipped numbers, and the question is
# whether the two levers have different addresses: the charge must not bend the
# attribution staircase, the rung must not bend the shape account, and the row's
# error must fall by both halves at once only when both are declared.
add_test(
    NAME verify.transientCoupledMultiTurbulentShapeAndInertiaSplit
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_transient_coupled_multi_turbulent_shape_and_inertia_split.py
        $<TARGET_FILE:chflow-coupled-1d3d-transient-multi>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/transientCoupledMultiTurbulentShapeAndInertiaSplit/case.dict)

# the three degrees of freedom this chain kept separate, all three present at once:
# where each edge's region inertia comes from (per arm AND per edge) crossed with
# which operator negotiates the interface - four tables x two maps, eight arms on
# one topology, and the separation measured rather than inferred
add_test(
    NAME verify.transientCoupledMultiThreeAxis
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_transient_coupled_multi_three_axis.py
        $<TARGET_FILE:chflow-coupled-1d3d-transient-multi>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/transientCoupledMultiThreeAxis/case.dict)

# the topology every earlier coupled run (architecture.md 5.25 through 5.31)
# could not put two negotiating rates into: both coupled edges sit between the
# same node pair, so the interface has no passenger.
# The worst-residual index hands over between them, one node drop has to close
# against two momentum rows carrying two different inertances, and the four
# attribution tables now set the ramp itself - a strictly ordered ladder of four
# time constants whose ratios are the sums M_X1 + M_X2
add_test(
    NAME verify.transientCoupledMultiBothParallel
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_transient_coupled_multi_both_parallel.py
        $<TARGET_FILE:chflow-coupled-1d3d-transient-multi>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/transientCoupledMultiBothParallel/case.dict)

# the one axis every coupled run in architecture.md 5.25 through 5.32 left in the
# 1D row: here the partner solves du/dt on its own mesh inside the time level, so
# the drop it reports is no longer a function of the rate alone. The interface
# still iterates on the friction part of it - refreshing on the whole drop is the
# a - c/R map that loses its crossing when dt < M/R - and grades the
# acceleration against the inertance the row installed, which is the first
# measurable error an attribution choice can carry
add_test(
    NAME verify.transientCoupledUnsteadyPartner
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_transient_coupled_unsteady_partner.py
        $<TARGET_FILE:chflow-coupled-1d3d-transient-multi>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/transientCoupledUnsteadyPartner/case.dict)

# the same partitioned interface with a nonlinear partner: edge X is a RANS
# channel, so its drop-per-rate is a function of the rate the interface is
# hunting for. One seed runs under both operators, which is what lets the
# checker read the saving as passes rather than as luck
add_test(
    NAME verify.coupledTurbulent
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_coupled_turbulent.py
        $<TARGET_FILE:chflow-coupled-1d3d>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/coupledTurbulent/case.dict)

# one region, two sources for its inertia: the same case runs once with the
# coupled pipe's own rho*L/A and once with the rho*V/A^2 the block measures on
# the region's mesh, so the checker closes every momentum row against the
# inertance that arm installed and reads the bore's shape error back out of the
# ramp the two arms produce
add_test(
    NAME verify.transientCoupledInertia
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_transient_coupled_inertia.py
        $<TARGET_FILE:chflow-coupled-1d3d-transient>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/transientCoupledInertia/case.dict)

# the RANS partner inside a time-stepping network: every time level agrees the
# interface again, and the checker reads the ramp back out of the series - the
# restoring term is the two partners' slopes, so the chord formula that closes
# the laminar case misreads this one by the exponents it never had to carry
add_test(
    NAME verify.transientCoupledTurbulent
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_transient_coupled_turbulent.py
        $<TARGET_FILE:chflow-coupled-1d3d-transient>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/transientCoupledTurbulent/case.dict)

foreach(tmode "blowdown;BLOWDOWN" "rlc;RLC" "slug;Slug" "valve;Valve")
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
# segregation, laminar pressure-driven start-up, dynamic-slip segregation
# and holdup (transported slip state + interfacial drag)
foreach(tcase twoFluidFront twoFluidSteadyHoldupWallis
              twoFluidSteadyHoldupZuber twoFluidSettling
              twoFluidPressureLaminar twoFluidDynamicSlipSettling
              twoFluidDynamicSlipHoldup twoFluidCompressibleHoldup
              twoFluidWellChain twoFluidHydrateSealed twoFluidHydrateFlow
              twoFluidJunction twoFluidNetworkSpan
              twoFluidEntrained twoFluidFullMomentum
              twoFluidSlugPulse twoFluidPressureFull)
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

# temperature equation: lumped-capacity relaxation onto the discrete closed
# form, the steady heat-exchanger recurrence + wall/stream enthalpy identity,
# and a shut-in well relaxing onto the geothermal line
foreach(tcase twoFluidThermalRelax twoFluidThermalExchanger
              twoFluidThermalGeothermal)
    string(REGEX REPLACE "^twoFluidThermal" "" mode ${tcase})
    string(TOLOWER ${mode} mode)
    add_test(
        NAME verify.${tcase}
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/check_two_fluid_thermal.py
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
