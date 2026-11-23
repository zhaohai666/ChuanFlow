# Batch C hook: Fluent-side F4 (LES SGS) + F12 (cavitation / melting) regressions.

add_test(
    NAME verify.lesSgs
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_les_sgs.py
        $<TARGET_FILE:chflow-3d-les>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/lesSgs/case.dict)

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
