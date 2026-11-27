# Batch C hook: Fluent-side F4 (LES SGS) + F12 (cavitation / melting) regressions.

add_test(
    NAME verify.lesSgs
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_les_sgs.py
        $<TARGET_FILE:chflow-3d-les>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/lesSgs/case.dict)

# RSM (LRR slow pressure-strain) second-moment closure on a homogeneous plane
# shear: the stresses are spatially uniform, so the checker re-derives the
# homogeneous-shear ODE system from the case text and replays the app's implicit
# Euler + relaxation + realizability scheme step by step (driving the mean shear
# by the reported shear trajectory), then checks k=tr/2, realizability, the
# developing anisotropy, the energy bound and loud refusals.
add_test(
    NAME verify.rsmShear
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_rsm_shear.py
        $<TARGET_FILE:chflow-3d-rsm>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/rsmShear/case.dict)

add_test(
    NAME verify.cavitationTransport
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_cavitation_venturi.py
        $<TARGET_FILE:chflow-3d-cavitation>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/cavitationVenturi/case.dict)

add_test(
    NAME verify.meltingFront
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_melting_front.py
        $<TARGET_FILE:chflow-3d-melting>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/meltingFront/case.dict)

# 5.196 (F12 fix): the latent heat of a step that CROSSES the phase interval.
# Three single-cell cases (all-symmetry faces, source only, so the exact answer
# is arithmetic) land the enthalpy below-then-inside (a), deeper inside (b) and
# past the liquidus (c); the checker inverts H/cp = T + (L/cp) f_l in closed
# form and asserts the printed T, the printed f_l, the enthalpy gain and that
# the sensible-only answer differs by more than 1 K (so a regression to the
# pre-fix assembly, which only priced cells whose OLD T was inside, cannot
# pass silently).
add_test(
    NAME verify.meltingLatentStep
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_melting_latent_step.py
        $<TARGET_FILE:chflow-3d-melting>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/meltingLatentStep/a/case.dict
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/meltingLatentStep/b/case.dict
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/meltingLatentStep/c/case.dict)
