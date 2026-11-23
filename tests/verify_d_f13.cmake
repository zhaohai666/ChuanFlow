# Batch D / F13: adaptive mesh refinement regression.
add_test(
    NAME verify.adaptiveRefine
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_adaptive_refine.py
        $<TARGET_FILE:chflow-3d-adaptive-refine>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/adaptiveRefine/case.dict
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/adaptiveRefine/caseUniform.dict
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/adaptiveRefine/caseReject.dict)
