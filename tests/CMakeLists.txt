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

# 5.46 Bingham yield stress in the steady network. The drop no longer vanishes
# with the rate -- at q = 0 the wall holds dpStatic = 4 tau_0 L / d -- so this is
# not another substitution into Re: the rate law needs an explicit dead branch
# before the Illinois bracket, and the flow side is Buckingham's 1921 laminar
# solution inverted for the wall stress by a monotone bracket. The checker
# therefore asserts the closed form at a known dp (eight fixed-end pipes) and at
# a known Q (case_driven, where the drop is the unknown), the exact static
# threshold (Q, v, Re and f all print 0 and dp stays at the asked drop), the
# identity f Re g(xi) = 64 that proves the reported factor is the solution rather
# than a Churchill number, Q/Q_Hagen-Poiseuille == g(xi), the wall stress that no
# tau_0 moves, the dead branch through the network Newton (a yielded series pair
# lands the free node exactly on its symmetric pressure, a stuck one empties the
# line), and the parse-time refusals -- of which the one that refused a yield
# stress together with a flow index is lifted again by 5.47, and is now asserted
# as a solved state rather than a refusal
add_test(
    NAME verify.steadyBinghamYield
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_steady_bingham_yield.py
        $<TARGET_FILE:chflow-1d-steady>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/steadyBinghamYield/case.dict)

# 5.64 / B6: close the "every edge below its own threshold" cell of 5.49. Zero
# kernel changes -- two checked-in dicts (5 pipes / 4 pipes) plus four arms
# derived in a temp dir. Thresholds 4 tau0 L/d, offered drops, the free-node
# initial guess (mean of fixed pressures) and Buckingham's Q are all re-derived
# from the dict by the checker: iterations=0 residual=0 converged on the dead
# network; p(A) shifted moves E 25000 -> 25500; dp == dpStatic still gives Q=0;
# one released parallel pipe flows at Buckingham's exact Q while iterations is
# STILL 0 (the zero is the free node's residual, not a dead-network shortcut);
# a 1e-6 source at the free node forces Newton to run 19 steps and raises E to
# p(B) + dp, where dp is Buckingham's inverse at that rate -- the HP cap on a
# dead branch is what keeps the Jacobian invertible there.
add_test(
    NAME verify.steadyAllBelowThreshold
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_steady_all_below_threshold.py
        $<TARGET_FILE:chflow-1d-steady>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/steadyAllBelowThreshold/case.dict)

# 5.47 Herschel-Bulkley in the same network -- tau = tau_0 + K gdot^n, the branch
# 5.46 refused. One closed form is the laminar profile integral,
# Q = pi R^3 tw^-3 K^(-1/n) n D^(1+1/n) [tau0^2/(n+1) + 2 tau0 D/(2n+1) +
# D^2/(3n+1)] with D = tw - tau0, which is Buckingham at n = 1 and the exact
# power law of 5.45 at tau0 = 0. That second limit is also what fixes 5.45's
# laminar constant: it carried the turbulent-collapse K (3n+1)/(4n) where the
# profile needs K ((3n+1)/(4n))^n, 25% of rate at n = 0.5 and exactly nothing at
# n = 1, which is why its own PowN1-vs-Newt control never saw it. The checker
# reads the closed form at two drops and on a series pair, re-derives it by
# quadrature rather than term by term, checks both limits against the branches
# that shipped them, keeps the threshold geometric (dp_static does not know n),
# and accounts for the 1e-9 slack on the two softest landings as the solver's own
# rate bracket rather than as a looser law
add_test(
    NAME verify.steadyHerschelBulkley
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_steady_herschel_bulkley.py
        $<TARGET_FILE:chflow-1d-steady>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/steadyHerschelBulkley/case.dict)

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

# a coupled edge declares the fluid of its 3D partner region. Five arms of one
# case (undeclared / declared-Newtonian / shear-thinning / shear-thickening /
# yielded) check the per-face apparent viscosity against K*gdot^(n-1)+tau0/gdot,
# the plane-channel closed form of the same law, the interface residual, and a
# bit-identical Newtonian limit; a sixth arm forwards it through the transient
# app. The 1D row's Re and parallel drop must not move at all.
add_test(
    NAME verify.coupledRheology
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_coupled_rheology.py
        $<TARGET_FILE:chflow-coupled-1d3d>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/coupledRheology/case.dict
        $<TARGET_FILE:chflow-coupled-1d3d-transient>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/transientCoupled1d3d/case.dict)

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

# A10 / 5.72: the double charge becomes a criterion. Same case as the split
# factorial above, zero kernel changes -- what is new is the policy read off
# the rowError identity: the bill is (M_row - M_dyn)*dq/dt with M_dyn a
# property of the block (all three rungs read the same 2.1435e6 kg/m4), the
# level-1 ladder is a per-rung bound (0.6635 / 0.5602 / 0.0046, strictly
# ordered), the rung is bit-invisible where dq = 0 (levels 2..7 spread 0.0),
# the collapse schedule is ONE settling level above the floor (level 2 =
# 1.07e-2, the block's own dpInertia) and then the tolerance floor itself
# (max 8.07e-4 < 1e-3 from level 3 on), and the free edge closes its identity
# to 2.3e-13 Pa with a positive, capped bill.
add_test(
    NAME verify.transientRowErrorCriterion
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_transient_row_error_criterion.py
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

# 5.71 / A9 the opposite corner of 5.32's topology: BOTH coupled edges in the
# SERIES position (A -> X1 -> J1 -> X2 -> Snk, one free node, no third
# branch). The network carries ONE rate -- q_X1 == q_X2 bit for bit in all six
# arms x 12 levels -- the rows ADD into a weight-free numerator (rowDrop ==
# R*q + M*dq/dt bit for bit, the free node closes both rows, the merged row is
# the plain sum), the constant inertial total splits EXACTLY as M1/(M1+M2) and
# moves the middle node without touching the rate, the four tables give four
# slopes inversely proportional to their M sums (the linear part of the
# relaxation, tau ~ 0.52 s against dt 5e-5), quasi-Newton keeps the rung, and
# both edges lead the residual schedule.
add_test(
    NAME verify.transientCoupledMultiBothSeries
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_transient_coupled_multi_both_series.py
        $<TARGET_FILE:chflow-coupled-1d3d-transient-multi>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/transientCoupledMultiBothSeries/case.dict)

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

# 5.48 the wall law of 5.45 - 5.47 in the transient momentum row. The row used to
# build `rho u d / mu` itself, so a declared `flowIndex` in a transient case was
# not refused but ignored: it ran Newtonian friction and printed a Newtonian
# rate. Both ends fixed-pressure makes each pipe one row,
# rho L/A dQ/dt = dp - dp_wall(Q), so the checker solves that recursion itself
# with the law written out from the case file and walks all 2001 reported states:
# the terminal rate must be the profile integral 5.47 corrected 5.45 with (one
# law, two kernels), the printed slope of the approach must be the law's own
# n dp/Q rather than the chord dp/Q or the Newtonian floor, n = 1 with K = mu
# must change no string, the odd continuation must show up as exactly one turn
# around from a negative start. Its ninth family used to end in a threshold
# refusal; 5.49 below turns that refusal into a solved case, so the refusal
# file is now `case_yieldNowSolves.dict` and the checker measures the yielded
# terminal rates instead
add_test(
    NAME verify.transientPowerLawRheology
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_transient_power_law.py
        $<TARGET_FILE:chflow-1d-transient>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/transientPowerLawRheology/case.dict)

# 5.49 the static branch: a yielded wall in the transient momentum row. The odd
# continuation 5.48 leans on has no root where |drive| <= dp_static = 4 tau_0 L/d,
# so a marching row there would chatter around zero with amplitude
# (dp_static - drive) dt / (rho L/A) instead of standing still. The row now
# latches at step entry (the committed rate is zero and the wall can hold the
# drive), clamps the crossing inside the solve, and reads an excess secant rather
# than a chord, because the yielded law's slope at rest is infinite. Verified by
# walking all 2001 states of 10 pipes against that recursion including `Q = 0`,
# by the four exactly-zero pipes each zero for a different reason (threshold,
# HB above its own threshold, gravity in the drive, a drive under the wall), by
# the stop happening *through* zero at a step the recursion also predicts, by the
# node storage balance closing to the solver's own exit residual when the pipe
# feeding it is held, and by the released pair ending on the rate the *steady*
# kernel of 5.46 computes for the same topology
add_test(
    NAME verify.transientYieldStaticHold
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_transient_yield_static_hold.py
        $<TARGET_FILE:chflow-1d-transient>
        $<TARGET_FILE:chflow-1d-steady>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/transientYieldStaticHold/case.dict)

# 5.50 the yield wall and the local loss on the same row: the cell 5.46 and 5.49
# each put in their own boundary list. The wall still asks for 4 tau_0 L/d at zero
# rate while the valve asks for nothing there, so the latch's drive contains no K
# (three Stuck pipes five decades apart in K print one shared trace, and the
# steady kernel reports no dpValve field for a pipe at rest), while every pipe the
# wall is *not* holding lands on the root of both losses together -- asserted
# against the checker's own bisection in both kernels, monotone in K across one
# threshold, and met between the kernels at t = 4. The stop is the row's
# arithmetic rather than the latch's: the same 2000 Pa margin under a K of 5e3
# stops at the same step as no valve (where the valve is worth 4.7e-02 Pa) and at
# step 5 under a K of 1e8 (worth 1.7e10 Pa at the starting rate), and the
# recursion with both terms reproduces all three traces and all three steps
add_test(
    NAME verify.transientYieldValve
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_transient_yield_valve.py
        $<TARGET_FILE:chflow-1d-transient>
        $<TARGET_FILE:chflow-1d-steady>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/transientYieldValve/case.dict)

# 5.65 / C2 (gap_analysis): the missing cell of the (valve, schedule, index)
# table -- the HB wall (5.47's flowIndex + consistency + yieldStress) and the
# BARE valve K (no schedule, 5.50's other half) on one row. Zero kernel changes:
# the parser has taken all four keys on one pipe since 5.50; nothing was refused,
# it was only never asked. What the index adds is that the valve's share of the
# asked drop is set by the wall law: the checker re-derives all 2001 states of
# all 16 pipes from the HB profile integral + K rho|Q|Q/2A^2 (worst 4.3e-10),
# asserts five decades of K on a held wall print one trace, one tau_0 prints one
# threshold at n = 0.5 and n = 1.5, the four-K family lands on its own roots
# monotone in K, the same K = 500 takes 2.33 Pa of the n = 0.5 pipe's drop but
# 308.76 Pa of the n = 1.5 pipe's (133x), the three stop paths are reproduced
# step for step, and the steady kernel meets the marched row on the same
# two-term root with the printed split read back through the inversion
add_test(
    NAME verify.transientYieldIndexValveK
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_transient_yield_index_valve_k.py
        $<TARGET_FILE:chflow-1d-transient>
        $<TARGET_FILE:chflow-1d-steady>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/transientYieldIndexValveK/case.dict
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/transientYieldIndexValveK/case_steady.dict)

# 5.66 endTime and reportEvery are decoupled: the loop takes
# ceil(endTime/dt - 1e-9) steps whose last one is clamped to the remainder and
# always reports the final state, so a case whose endTime is NOT a multiple of
# dt (0.1 with dt 0.03 -> 0.03/0.03/0.03/0.01) lands on endTime exactly, a
# near-multiple endTime (0.12/0.03 = 4.000000000000001) is not miscounted, and
# a multiple-of-dt case prints the same sample set the old fixed-step loop did.
# Every reported state is re-derived from the implicit-Euler row with the
# per-step h -- closed form for the Newtonian pipes, bisection in wall-stress
# space for the yielded one.
add_test(
    NAME verify.transientTimeDecoupled
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_transient_time_decoupled.py
        $<TARGET_FILE:chflow-1d-transient>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/transientTimeDecoupled/case.dict)

# B5 the full coupling grammar in ONE case: ddtFrom + inertiaFrom +
# convectiveToRow and the per-edge profile have never met in a single file
# before (5.40/5.41 carried two tables + profile, 5.43/5.44 the other two).
# Five arms on 5.41's topology, meshes, dt and seeds; the new cell is
# mix_charge -- quasi AND slab AND charged on one edge at once -- whose row
# pays M*dq/dt with the slab's number and hands back dpConv exactly where the
# charge table says. The base arm reproduces 5.41's published p_net states and
# pass sequence field for field, the row closes with zero slack in all 80
# states, the charge moves exactly one term of one row (the neighbour is
# bit-identical), and the chord stays inside the interface tolerance.
add_test(
    NAME verify.transientCoupledMultiTurbulentTripleTable
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_transient_coupled_multi_turbulent_triple_table.py
        $<TARGET_FILE:chflow-coupled-1d3d-transient-multi>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/transientCoupledMultiTurbulentTripleTable/case.dict)

# 5.68 / D5 CFL-adaptive time stepping: `cfl` in time{} switches the loop to
# dt_next = clamp(cfl * min_k(L_k A_k / |q_k|), dtMin, dtMax) read off the
# flows the step just committed -- first step is the case's dt, a resting
# network takes dtMax, t accumulates, the last step is clamped to the
# remainder and t snaps to endTime. The case's schedule visits every regime
# (initial / ceiling / rule / floor / remainder) and stays laminar so every
# state re-derives in closed form. The checker re-derives the whole schedule
# from the states -- the schedule is a property of the physics, not a
# remembered list -- plus the fixed march (cfl absent) untouched and three
# refusals.
add_test(
    NAME verify.transientAdaptiveCfl
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_transient_adaptive_cfl.py
        $<TARGET_FILE:chflow-1d-transient>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/transientAdaptiveCfl/case.dict)

# 5.69 / A7 the solver exit-gate policy as a published contract: the five
# named constants in SolverGates.h equal the documented values, the
# architecture document states each gate by name AND value, the call sites
# reference the constants (no second source), and the gates' measured
# footprints hold -- steadyBinghamYield within 1e-11 of Buckingham's exact
# solution (bracket gate far below its nominal 1e-9 width), tol 1e-16
# bit-identical to the default (inner gates are floors), maxIter 0 fails
# loudly, and the transient march reproduces the exact implicit row with a
# real nonzero drift under 1e-9 (the Newton exits ON its gate, not on a
# fixed point).
add_test(
    NAME verify.solverGates
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_solver_gates.py
        $<TARGET_FILE:chflow-1d-steady>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/steadyBinghamYield/case.dict)

# 5.73 / A6 the annulus cross-section: FB 1958 laminar anchor (f.Re = C(kappa)
# closed to 1e-9, Q to 1e-9), a monotone kappa ladder falling to the circular
# 64, the power-law FB profile (lambda by wall matching) with its own C_pl
# ladder toward 64((3n+1)/4n)^n, Dh-consistent reports at 1e-9 on every tier,
# and loud rejections outside the declaration scope. The circular branch is
# bit-still: the whole pre-existing suite is the witness.
add_test(
    NAME verify.steadyAnnulus
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_steady_annulus.py
        $<TARGET_FILE:chflow-1d-steady>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/steadyAnnulus/case.dict)

# 5.74 / A2 gel strength and sub-threshold creep in the transient row: the
# gel arm opens exactly on the latch recursion's predicted layer (22 layers
# after the plain twin -- the restart pressure as calendar), the creep arms
# answer pi d^4 drive / (128 mu_eff L) at every step of the 600-step march
# and stay linear at 1.25x, the never-opening wall stays 5.49's Q = 0, the
# opened walls park on 5.46's Buckingham law, and the five declaration
# scope violations are refused loudly. A zero declaration is the old row
# bit for bit -- the standing suite is the witness.
add_test(
    NAME verify.transientGelCreep
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_transient_gel_creep.py
        $<TARGET_FILE:chflow-1d-transient>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/transientGelCreep/case.dict)

# 5.70 / A3 the momentum-flux factors betaIn/betaOut, declarable on a
# non-coupled transient pipe: the row charges rho*(betaOut - betaIn)/A^2 *
# q|q| -- the same account a coupled interface installs (5.43), now owned by
# the case. Three geometrically identical pipes separate by the declaration
# alone (billed slower, refunded faster, every state on the implicit row's
# quadratic in closed form), the term flips sign with q|q| (driven backwards
# bit-identical in |q|), yield + beta share one row (bisection in wall
# stress), a held pipe reports NO account (Q = 0 kills it, 5.50's statement),
# and a factor below 1 or a coupled pipe declaring beta is refused. No
# declaration means the old row bit for bit.
add_test(
    NAME verify.transientBetaMomentum
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_transient_beta_momentum.py
        $<TARGET_FILE:chflow-1d-transient>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/transientBetaMomentum/case.dict)


# 5.51 the published laminar-turbulent locus of a yielded pipe, read off the
# pipe's own wall law. The turbulent-yield FRICTION closure could not be
# transcribed with its constants intact, so nothing is re-priced here: what
# lands is He, Re_B, the three published bands of Re_B,t(He), and the one number
# that is the pipe's own physics -- the drop dp_t its 5.46/5.47 law asks for at
# the locus rate. Asserted as geometry + fluid rather than state (twins at
# different drops bit-identical), against the bands' own fits, through the
# identity Trel = Re_B,t / He^0.5 (which forces Trel = 25 across the He^0.5
# band and 80 He^-0.10 in the mid band), and inert on the answer (every rate
# still the exact root of the asked drop, f Re_B g = 64 past the locus, the
# tau_0 -> 0 pipe on its Newtonian twin times the bracket's own factor). The two
# band joins are measured, not smoothed: 1e5 steps down 1.17%, 1700 steps down
# 25.3% because the paper's bridging Eq (17) is unread
add_test(
    NAME verify.steadyYieldTransition
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_steady_yield_transition.py
        $<TARGET_FILE:chflow-1d-steady>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/steadyYieldTransition/case.dict)

# 5.52 the valve *timetable* on the same row as the yield wall's static branch --
# which is the cell 5.50 put in its own boundary 3 and 4. The two switches on
# that row are read from different places: the latch compares the drive (which
# carries no K) with 4 tau_0 L/d, while the schedule only prices the term the
# moving row pays. So a closing valve can move the stop step (54 / 72 / 184 for a
# ramp ending at 0.15 / 0.6 / never, and monotone in it) but cannot cause it, a
# ramp that arrives after the wall already stopped the pipe is invisible in all
# 2001 printed tokens and in all 2001 latch flags, no schedule wakes what the
# wall holds while its Newtonian twin returns to its own rate, and a pipe already
# held at rest prints the schedule-free trace string for string. Above the
# threshold the valve is a discount and not a refusal: K = 1e6 is worth 9.8x the
# whole drive there and still holds nothing. Sampling convention is settled the
# cheap way -- the app stops on the t + dt step, the t variant is still moving one
# step later, and for a ramp that never overlaps motion the two are one trace
# (0.000e+00 apart). All 12 x 2001 states reproduce this file's own recursion of
# both laws to 2.14e-09 relative on the moving states and 2.67e-15 m3/s in
# absolute units, and the settled row meets the steady kernel with every valve
# frozen at its end-of-schedule K (2e-9), with no dpValve field on a pipe
# reporting Q 0
add_test(
    NAME verify.transientYieldValveSchedule
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_transient_yield_valve_schedule.py
        $<TARGET_FILE:chflow-1d-transient>
        $<TARGET_FILE:chflow-1d-steady>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/transientYieldValveSchedule/case.dict)

# 5.53 the free node between them: the yield wall's static branch, a valve
# timetable and a node pressure the solve has to find -- the cell 5.50 put in its
# boundary 4 and 5.52 in its own boundary 2, which needed no kernel change, only
# a case file that stops pinning both ends. Nine pairs, each a held yielded pipe
# into a floating node and a Newtonian bleed out of it. The drive the latch
# compares against is now whatever the network says, so a *neighbour's* schedule
# owns the release: the step goes 253 / 299 / 1830 for K = 0 / 5e3 / 1e6 on the
# bleed and a 1e8 valve pushes it past the end of the run, a ramp that lands
# before the drain gets there cancels it and one that lands after is invisible
# until the sampled t + dt steps inside it, and a neighbour that seals and then
# reopens wakes the wall at step 1290 -- 0.578 s after the schedule first moved
# anything, because the network, not the valve, is what reaches the node. On the
# held pipe's own row the same two prescriptions release on K0's step with an
# identical node history and cost rate only from the first moving state. The wall
# answers one number: all 7 releases happen at 21990.9 - 21999.2 Pa against
# pA - 4 tau_0 L/d = 22000 Pa across three decades of K and two schedules, and
# none re-latches -- the three pairs whose neighbour seals after the release peak
# at 21980.9 - 21983.0 Pa and the steady kernel settles those same nodes at
# 20450.6 - 21634.2 Pa, all below it. All 9 x 2001
# states of the 4 s run and both pairs of a 30 s run (15000 marched steps)
# reproduce this file's own copy of the two rows, the node row and the sampled K
# (worst 1.51e-09 relative, 1.10e-16 absolute), and every floating node closes to
# the solver's own 1e-11 exit criterion -- a held pipe enters continuity as the
# forced zero it prints. Against the steady kernel the nine settled rates are this
# file's bisection to 2.0e-8 with dp_h + dp_b = 1e4 Pa to 4e-9, dpValve = K rho
# Q^2/(2 A^2) to 2.4e-12, the three all-open placements on one rate and the four
# placements of one 1e8 on another within 1.1e-8 -- while the node the valve sits
# opposite splits by exactly that valve's 1603.08 Pa drop. Held for 4 s is the
# calendar: the same pairs trickle at 6.1015e-08 m3/s with no latch in the way,
# the 30 s march hands R1 back at t = 24.25 s (entry drive 7999.48 over 8000) and
# K8 is still 304.35 Pa short at 30 s with its node falling monotonically at
# 102.55 Pa/s, which is that drain's own rate over the node's compliance
add_test(
    NAME verify.transientYieldValveNodeRelease
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_transient_yield_valve_node_release.py
        $<TARGET_FILE:chflow-1d-transient>
        $<TARGET_FILE:chflow-1d-steady>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/transientYieldValveNodeRelease/case.dict)

# 5.54 the branching node under the latch: two yielded walls on ONE free node, a
# scheduled drain, and the two cells 5.53 had to declare unreachable -- its
# boundary 2 (a node with more than one pair) and its boundary 6 (a latch that
# closes again after it opened). Seven arms, 13 yielded walls, 8 Newtonian pipes,
# 9 nodes, and no kernel change: the node row was already degree-generic, the
# entry latch was already a test on Q_prev == 0, and the in-step clamp was
# already there. The second wall gives the node something the first never had, a
# ceiling: while the weaker wall flows the limit is *that* wall's release
# pressure, so every stronger wall crosses it from below and can be put back to
# sleep. The latch reads one number in both directions -- 14 release events on 11
# walls caught 1.00 - 16.54 Pa over 4 tau_0 L/d and 3 re-latches 3.10 - 3.55 Pa
# under the same number, 17 flips in the file, so there is no hysteresis to fit.
# Who re-latches is an ordering, not a
# calendar: Sa 8000 under Wa 6000, Wt/Ws 6000 under Ts/Tt 4000, and arm `i`'s two
# identical 8000 walls re-latch neither, as arm `u`'s single wall re-latches
# nothing. The inequality predicts two walls the file refuses, each for its own
# reason: arm `t`'s 8000 wall is never asked (its node bottoms at 22593.53 Pa,
# 593.53 Pa high) and arm `e`'s never comes back under (the largest node after
# its release is 21979.46 Pa, 20.54 Pa short of 22000, and its steady destination
# 620.27 Pa short). The ceiling is a wall pressure only when the weakest inflow is
# a wall: arm `p`'s Newtonian neighbour makes it pA, and its wall gets no closer
# to the drive it needs than 5101.52 Pa of 8000. The held row is silent in the
# strongest sense the print allows -- arm `t` is arm `s` plus a wall that never
# wakes, and the two agree to 0.0 Pa on the node and token for token on every
# flowing pipe, in the transient and in the steady kernel (Ft = Fs = 22657.9280).
# 23000 states of the 4 s and 18 s runs reproduce this file's own three rows, node
# row and sampled K (1.36e-09 relative on rates, 7.83e-12 on pressures, no latch
# flag differs), all 19658 held states print the single token 0, and the node
# balance closes to 4.951e-14 m3/s on the printed tokens alone. The calendar buys
# 7001 steps of seal: arm `a`'s weak wall is 1048.65 Pa under its release pressure
# after 2 s and still 394.36 Pa under it after 14 s, the per-thousand-step shrink
# fallen to 0.99392, while the same window hands the strong wall 7295 steps of
# sleep -- `not yet` and `not ever` are the same asymptote read at two times
add_test(
    NAME verify.transientYieldValveBranchRelease
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_transient_yield_valve_branch_release.py
        $<TARGET_FILE:chflow-1d-transient>
        $<TARGET_FILE:chflow-1d-steady>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/transientYieldValveBranchRelease/case.dict)

# the same latch seen from the OTHER side of the node: 7 arms, 14 yielded walls
# and 19 pipes marched 2001 + 9001 states, where a drain wall's drive is
# `pF - pB` and the free node is the upstream end. 8 flips (5 wakes, 3
# re-latches) and the re-latching law is measured, not assumed -- Sa by Wa
# (8000 by 6000, inflow side), Wz by Tz and Wb by the NEWTONIAN drain Db (both
# outflow side): a wall is re-latched iff a row on the same side of its node has
# a strictly smaller threshold, and the 11 walls with no such carrier never
# re-latch. The floors are approached and never reached: Ff 303.12 Pa above
# pB + 6000 at 4 s and still 82.54 Pa at 18 s, Fz 292.22 -> 60.44 Pa, arm z's
# per-thousand-step shrink 0.43708 -> 0.73507. The dead band is a fixed point of
# the WHOLE solve -- Fd prints 27000.0 in all 2001 states -- and steady says the
# same arm has an EXACT solution (iterations=0 residual=0, converged) at four
# nodes that are the solver's own guess 25000; the seven rows and four nodes are
# token-identical in the 7-arm and 3-arm files, so what keeps Newton from
# converging there (iterations=200, residual 0.0321364485 = arm a's node
# imbalance at the sources' 6e-5 scale) is the two arms that flow, and at
# maxIter 50000 Fa creeps to 21379.7317 Pa with residual 2.6336e-11 -- a rate,
# not a refusal. Every node row closes on the printed tokens alone to
# 4.961e-14 m3/s, under one unit of the last printed pressure digit
add_test(
    NAME verify.transientYieldValveSplitRelease
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_transient_yield_valve_split_release.py
        $<TARGET_FILE:chflow-1d-transient>
        $<TARGET_FILE:chflow-1d-steady>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/transientYieldValveSplitRelease/case.dict)

# the dead band itself, swept from the inside: 12 arms, 24 yielded walls, the
# free node's OWN initial pressure walked across [pA - 4 tau_0 L/d, pB + 4 tau_0
# L/d] = [22000, 28000] with the kernel untouched.  6 arms parked (both edges
# included -- `<=`, not `<`: Fd prints 22000.0 in all 2001 states and Fh 28000.0)
# and 6 sliding (3 below, 3 above, mirrored token-for-token: Fa + Fk = 50000 Pa
# at every one of 2001 states, max 0 Pa off, identical rates).  Boundary bracketed
# to 1 Pa on each side: Fc 21999 -> 21999.00196 moves, Fd 22000 (drive 8000 =
# threshold 8000) holds; Fi 28001 moves, Fh holds.  Nothing overshoots -- 20001
# states put the closest slider 0.9807 Pa above its floor.  The approach is a 1/t
# tail whose coefficient comes out of the documented Bingham near-threshold
# quadratic (both the constant and linear terms of 1 - 4t/3 + t^4/3 cancel at
# t=1), k_Q = 4.9087e-13 m3/(s Pa2) so 1/gap gains 4.9087e-4 Pa-1 per second:
# measured slopes land within 1.45%, the integrated 40 s arrival within 3.62%,
# d ln Q / d ln gap = 1.9647 - 1.9974, and the per-2 s shrink climbs 0.22988 ->
# 0.89994 toward 1 without ever reaching it.  A wall that woke never re-arms
# (smallest printed waking rate 4.78156e-13 m3/s), the calendar cannot park a
# slider (arm m == arm f token-for-token while K hits 1e8), and the band is local:
# a 3-arm file reproduces the parked arm's tokens exactly and the sliders to one
# printed quantum.  Steady pays zero on all 22 rows and leaves every swept node at
# the solver's own guess 25000 Pa -- which is the band centre by arithmetic,
# (pA - thr + pB + thr)/2 = (pA + pB)/2, so the blind default is always parked
add_test(
    NAME verify.transientYieldParkingBand
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_transient_yield_parking_band.py
        $<TARGET_FILE:chflow-1d-transient>
        $<TARGET_FILE:chflow-1d-steady>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/transientYieldParkingBand/case.dict)


# 5.57: does the dead band know the flow index?  Herschel-Bulkley walls
# (yieldStress 2 Pa with flowIndex / consistency set per pipe) parked at the
# same 8000 Pa threshold for n = 0.5, 1, 2: yieldThresholdPa is 4 tau0 L/d and
# reads no rheology, so the desk -- 25000 / 22000 / 28000, one token each over
# 2001 and 20001 states, and arms b/c HELD at i=1 with drive exactly 8000 -- is
# index-blind while everything downstream of the latch is not.  The HB section
# integral reproduces the tokens from the FIRST decided step on every sliding
# arm (1312 sampled rates per index, worst relative error 4.698e-07 / 3.612e-07 /
# 2.874e-07) and at n=1 the general bracket is 5.56's Bingham bracket to 6.3e-16.
# The first row is the ramp: inertia takes 539.5 / 1635.8 / 1821.0 Pa of the
# 8000 Pa drive (6.7 / 20.4 / 22.8%), monotone in n, so the wall starts at
# 2.365 / 2.091 / 2.045 instead of the rest-state 2.500.  The decay is a
# geography, not a number: the bracket's own d ln Q / d ln gap reads 3.000 /
# 2.000 / 1.500 one nanoPa above yield -- exactly (n+1)/n -- and 2.000 / 1.000 /
# 0.500 five thousand yield stresses above it -- exactly 1/n -- while the real
# trajectories only ever cross between, with 18/26/34 chords strictly inside
# (1/n, m) and each one within 0.027% / 3.922% / 2.023% of the bracket's own
# chord over the same two gaps.  The near-edge arms read m off their own kernel:
# 1.99976-1.99983 at n=1, 1.49985-1.49995 at n=2 (0.012% / 0.010% of m), and at
# n=0.5 they contribute NO chord at all because a 1 Pa start moves less than one
# printed quantum in 40 s.  The straight line moves from 1/gap to gap^(-1/n):
# measured slopes 1.66666e-08 / 4.90772e-04 / 1.46335e-02 per Pa per second
# (1.859% / 0.021% / 0.011% off the leading term A_n/n, and at n=1 that is 5.56's
# 4.9087e-04) -- the (n, K) pair moves that rate by 2.9e4 from n=0.5 to n=1 and
# by another 30 to n=2 while the band does not move at all.  Same 1 Pa start, 40
# s later: gap 1.000 / 0.981 / 0.398 Pa.  Mirror arms d/f and e/g are
# token-for-token identical (pressures sum to 50000 within 0.0e+00), the latch
# never re-arms (smallest waking rate 8.18e-18 / 4.72e-13 / 7.35e-12 m3/s), all
# 14 static-hold flags at t=0 read 0 (the scan must start at i=1), and the
# steady kernel is as index-blind as the transient one: six rows print Q 0, three
# nodes print 25000, iterations=0
add_test(
    NAME verify.transientYieldIndexParkingBand
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_transient_yield_index_parking_band.py
        $<TARGET_FILE:chflow-1d-transient>
        $<TARGET_FILE:chflow-1d-steady>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/transientYieldIndexParkingBand/case_5.dict)


# 5.58: once the band is index-blind, what can a calendar touch?  Ten HB arms per
# index (n = 0.5 / 1 / 2) share one 4 s clock carrying valveClose(0.5 -> 1.0,
# closedK 1e8) + valveOpen(2.0 -> 2.5) on some of their rows.  The desk and the
# latch never see the calendar: 25000 / 25000 / 22000 / 28000 each hold ONE token
# over 2001 states at every index, arm h (schedule on BOTH rows) is arm a
# token-for-token in 2001 node + 2001+2001 rate tokens because a held row becomes
# Q = 0 before kv is read, and arm b is still HELD at i = 1 with drive exactly
# 8000 while K(t) sweeps to 1e8 behind it.  The armed row's own arithmetic closes
# only with the valve term: the HB section integral reproduces its printed rate
# from the FIRST decided step to 3.757e-10 / 1.001e-09 / 8.235e-10 with
# K(t)rho|Q|Q/2A^2 moved off the drive and to 6.4e-01 / 1.4e+01 / 1.9e+01 without
# it.  One calendar, three latencies: armed and unarmed rows agree through steps
# 0-250 (K = 0 there, `valve 0`), the far slider leaves its twin at step 251 --
# the first sample past closeStart, K = 4.000e+05 -- at every index, while the
# 1 Pa armed arm leaves it at step None / 757 / 268 (n = 0.5 / 1 / 2).  What the
# closed valve costs, per index: the rate at t = 1.5 s is throttled 1.59x / 7.22x
# / 8.20x and the approach loses 24.1797 / 279.2209 / 253.1707 Pa of 167.8146 /
# 1540.2417 / 1822.4064 Pa (14.4% / 18.1% / 13.9% -- the fraction is the
# calendar's, the distance is the index's), and the money is real: the valve term
# peaks at 323.0 / 2483.1 / 3701.8 Pa, i.e. 4.0% / 31.0% / 46.3% of the 8000 Pa
# threshold the latch reads, monotone in n, while the wall stress the armed row
# implies stays within 3.1822% / 10.2906% / 11.1758% of its unarmed twin's.  With
# K back to exactly 0 the debt is repaid, strictly and monotonically, at an
# index-dependent rate -- 8.3% / 59.6% / 70.2% by t = 4 s -- and in every one of
# the 750 post-open states the armed arm flows FASTER than its twin (1.03..1.04 /
# 1.59..3.57 / 1.93..5.31x) because the same law is handed more drive.  The band
# itself does not move at any index: no arm reaches 22000 or 28000 (nearest
# 0.89252 Pa at n = 2, a full 1.00000 Pa at n = 0.5) and no throttled row ever
# re-arms -- smallest rate 8.1797e-18 / 4.8957e-13 / 2.6280e-11 m3/s.  Mirror e/f
# is token-for-token exact (nodes sum to 50000 within 0.0e+00, 2001/2001 rates
# identical) because |Q|Q is odd, all 20 t = 0 flags read 0, and every node row
# closes on the printed tokens alone to 4.9163e-14 / 4.9640e-14 / 4.9382e-14 m3/s.
add_test(
    NAME verify.transientYieldIndexValveBand
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_transient_yield_index_valve_band.py
        $<TARGET_FILE:chflow-1d-transient>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/transientYieldIndexValveBand/case_5.dict)


# 5.59: the re-latch moved onto an HB wall -- 5.54's two-way latch rewritten with
# (n, K) = (0.5, 1.0) / (1.0, 0.1) / (2.0, 0.02), 4 arms x 10 pipes x 6 nodes,
# 2001 states, and NO kernel change (threshold :28-:31, holdsStatic :40-:43,
# in-step clamp :353-:355, degree-generic node row :395-:416).  The desk is as
# index-blind as 5.57 said: Wa's flag changes at 1 and 180, arm u's lone wall at
# 1 and 318, arm p's at 1 only -- the same three integers in all three tiers --
# and Fa / Fu print ONE token for their first 180 / 318 states in all three files.
# The arrival is not: arm a's 8000 Pa wall wakes at step 318 / 333 / 343
# (0.636 / 0.666 / 0.686 s, a 25-step spread) on entry drives 8001.7883 /
# 8006.1182 / 8000.6927 Pa, each the first above the threshold and each within
# 80 Pa of a number no tier can see the inside of.  A sleep needs BOTH halves of
# :340-:360, and the coast between them is the law's: at n = 1 the drive falls
# under 8000 at step 507 and the flag flips there (lag 0); at n = 2 the drive is
# already under at 465 while the row still prints 7.28766e-09 m3/s and the flip
# only comes 7 steps (1.4e-02 s) later.  At n = 0.5 the desk is never asked -- the
# smallest entry drive after the wake is 8001.7883, and the sealed window's best
# node is 21400.49321 Pa against the 22000 Pa release pressure, 599.507 Pa short,
# while n = 1 / n = 2 overshoot by 1192.446 / 1544.178 Pa.  The race is supply
# vs leak at the first sealed state: 1.30083e-07 / 1.83744e-06 / 2.0945e-06 m3/s
# against 5.43469e-08 / 6.73837e-08 / 7.0236e-08, ratio 2.394 / 27.268 / 29.821,
# and the node climbs 108.2611 / 1210.4094 / 1391.8863 Pa before the reopened
# drain takes over.  Every unheld step with a root closes on the HB section
# integral (worst 1.957e-08 / 4.561e-07 / 4.429e-08 relative on Sa) and in the
# whole run EXACTLY ONE unheld step implies a sub-yield stress -- step 506 / 471,
# one before the clamp -- so holding is a lost root, not a residual.  n = 1 with
# K = mu is 5.54's wall written twice: identical flags, ladders [1, 333, 507, 1530]
# and [1, 318] in both files, nodes within 1.0e-07 Pa, and bit-identical tokens on
# the arm whose wall never wakes.  87 judgments, three case files, 0 kernel lines
add_test(
    NAME verify.transientYieldIndexValveRelatch
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_transient_yield_index_valve_relatch.py
        $<TARGET_FILE:chflow-1d-transient>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/transientYieldIndexValveRelatch/case_5.dict)

# 5.60: the tau_0 x flowIndex cross 5.57 and 5.59 both left empty -- 4 arms x 3
# tiers = 12 cells, strong wall tau_0 = 1.75 / 2.0 / 2.25 / 2.5 against a 1.5
# neighbour, 12 pipes / 6 nodes / 2001 states per file, NO kernel change
# (threshold :28-:31, desk :40-:43, clamp :353-:355, node row :395-:416).
# The desk is blind to BOTH axes: all twelve weak walls release at step 180 in
# all three tiers, and the four arms of a tier are the same strings until the
# earliest strong wake (239 / 240 / 241).  The re-latch is a WINDOW in tau_0, not
# a half-line: the census is 4 SLEPT (n = 1 / 2 at tau_0 1.75 and 2.0), 2 AWAKE
# never asked (n = 0.5 at 1.75 / 2.0, smallest entry drive 7005.6287 / 8001.7883),
# 1 LATE (n = 0.5 at 2.25 wakes at 1479, 129 steps AFTER the seal broke) and
# 5 NEVER (tau_0 2.25 / 2.5 at n = 1 / 2, 2.5 at n = 0.5).  The two axes buy
# different things: tau_0 moves the release pressure 3000 Pa and the sealed
# window's best node only 51.36260 / 45.20672 / 44.18058 Pa, while n moves the
# best node 2138.60589 / 2143.68512 / 2145.78792 / 2145.78792 Pa at fixed release
# -- a ratio of 41.6, and the four n-spreads agree within 7.18202 Pa.  A positive
# margin is NOT a sleep: n = 0.5 arm c tops out 397.040 Pa ABOVE its release
# pressure and still never sleeps, because it was held for all 851 sealed states.
# The coast is the law's AND the threshold's: lag 0 / 0 at n = 1, 5 / 7 at n = 2.
# Every unheld step with a root is the HB section integral (worst 1.876e-06
# relative, n = 1 arm a at 961) and the only sub-yield unheld step in a sleeping
# cell is the one before the clamp (962 / 506 / 708 / 471).  Arm b is 5.59's arm a
# in a new file: identical flags and ladders, nodes within 1.0e-07 Pa, rates worst
# 1.780e-08 relative; all 12 node rows close on printed tokens alone to
# 4.9661e-14 m3/s, all 24 t = 0 flags read 0, smallest decision gap 129 steps.
# 146 judgments, three case files, 0 kernel lines
add_test(
    NAME verify.transientYieldTauIndexRelatch
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_transient_yield_tau_index_relatch.py
        $<TARGET_FILE:chflow-1d-transient>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/transientYieldTauIndexRelatch/case_5.dict)

# 5.61 -- source schedule x HB yield wall x drain valve. New kernel hook:
# `sources{ <name> { node N; schedule { openRate; close { start; end; }
# open { start; end; } } } }`, mirror of valveClose / valveOpen on a pipe,
# evaluated at the same implicit new time level t + dt (transient.cpp:125-:145
# shared ramp helper, sourceRateAt; :329-:337 per-step reset+sum into
# nd.source; :400 node balance; :480 liquid inventory with sourceFraction).
# Parser rejects `sources` entries with neither `rate` nor `schedule`. The
# desk (transient.cpp:40-:43) and clamp (:353-:355) are schedule-blind, so
# the release step (318) is identical in case_off and case_on; the schedule
# only moves the re-latch and re-release steps.
# case_off transitions: [1, 318] (events(), 2 entries).
# case_on  transitions: [1, 318, 531, 1466] (events(), 4 entries).
# Release drive 8006.59 Pa > threshold 8000, one step earlier 7999.38 Pa.
# HB section integral worst 1.4e-8 rel, node closure worst 1.2e-7 Pa/step,
# 600-state prefix bit-for-bit identical across the two files.
# 6423 judgments, two case files, 0 kernel lines (counted at the next hook).
add_test(
    NAME verify.transientSourceYieldWake
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_transient_source_yield_wake.py
        $<TARGET_FILE:chflow-1d-transient>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/transientSourceYieldWake/case_off.dict)

# 5.62 -- band edge: pB lowered from 20000 to 19500 so arm d (n=0.5, tau_0=2.5)
# wakes for the first time. Zero kernel change -- the same holdsStatic, in-step
# clamp, and degree-generic node row as 5.49/5.54/5.60. The only move is the
# sink pressure, which lets the node drop below the release floor pA - 4 tau_0
# L/d = 20000 Pa that was asymptotically unreachable at pB = 20000.
# pB scan at n=0.5: 19800->NEVER, 19780->4 states, 19760->61, 19740->103,
# 19720->136, 19700->163, 19500->313 flowing states out of 2000.
# Arm c in n=0.5 now has re-latch cycle: release 388, re-latch 929, re-release
# 1390 (was release 1479, no re-latch in 5.60).
# Node Fd at step 2000: 19752.84 Pa (below 20000 threshold floor, only possible
# because pB = 19500 < 20000).
# Weak walls release at same step [1,1,1,1] across all tiers -- the desk is
# blind to the far-end boundary.
# 84 judgments, three case files, 0 kernel lines
add_test(
    NAME verify.transientYieldTauIndexBandEdge
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_transient_yield_tau_index_band_edge.py
        $<TARGET_FILE:chflow-1d-transient>)


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
