# Batch D item F7 verification hook, included by verify_3dphys.cmake.
add_test(
    NAME verify.mixturePbm
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_mixture_pbm.py
        $<TARGET_FILE:chflow-3d-mixture-pbm>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/mixturePbm/case.dict)
