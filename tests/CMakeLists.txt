find_package(Python3 COMPONENTS Interpreter REQUIRED)

# ---- C++ unit tests ----
add_executable(kernel_tests unit/kernelTests.cpp)
target_link_libraries(kernel_tests PRIVATE cf3d)
add_test(NAME kernel.unit COMMAND kernel_tests)

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
    NAME verify.poiseuille
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_poiseuille.py
        $<TARGET_FILE:chflow-simple-foam>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/poiseuilleChannel/case.dict)
