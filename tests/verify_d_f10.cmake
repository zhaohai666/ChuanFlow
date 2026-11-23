# Batch D / F10: independent RK4 grading of the two-step chain-branching
# 0D reactor (history, invariant, ignition delay, extinction threshold).
add_test(
    NAME verify.combustionDetailed
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_combustion_detailed.py
        $<TARGET_FILE:chflow-3d-combustion-detailed>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/combustionDetailed/case.dict
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/combustionDetailed/caseIgnite.dict
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/combustionDetailed/caseExtinguish.dict
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/combustionDetailed/caseReject.dict)
