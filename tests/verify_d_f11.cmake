# Batch D item F11 CTest fragment. Included by verify_3dphys.cmake.
# The checker independently integrates the rigid-body ODE (closed forms) and
# the uniform-translation scalar invariance law, plus a static twin and loud
# rejections of negative mass / drag / inertia.
add_test(
    NAME verify.dynamicMesh6Dof
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_dynamic_mesh_6dof.py
        $<TARGET_FILE:chflow-3d-dynamic-mesh>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/dynamicMesh6Dof/case.dict
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/dynamicMesh6Dof/caseStatic.dict
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/dynamicMesh6Dof/caseReject.dict)
