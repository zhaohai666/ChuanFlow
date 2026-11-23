# F9 fragment: included by verify_3dphys.cmake. The checker rebuilds the
# analytic S_N slab solution (cell-averaged exponential attenuation, closed
# form wall flux), the exact discrete balance, the S2S view-factor limit, the
# P1/DO trend and five loud refusals from the case constants alone.
add_test(
    NAME verify.radiationDoS2s
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_radiation_do.py
        $<TARGET_FILE:chflow-3d-radiation-do>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/radiationDo/case.dict
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/radiationDo/caseP1.dict
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/radiationDo/caseReject.dict)
