add_library(cf_common STATIC
    src/Dictionary.cpp
    src/Matrix.cpp
)

target_include_directories(cf_common PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Optional Eigen sparse solver backend (README M6 / A11). Header-only; the
# vendored 3.4.0 tree is used when the top-level switch is on. Off by default
# so the native CG/BiCGStab path (and every bit-identity checker) is unchanged.
if(CHFLOW_WITH_EIGEN)
    target_compile_definitions(cf_common PUBLIC CHFLOW_WITH_EIGEN)
    target_include_directories(cf_common PUBLIC ${CF_EIGEN_DIR})
endif()

# Optional PETSc sparse solver backend (README M6 / A11). Only compiled when
# find_package found a PETSc library; otherwise the native solver stays.
if(CHFLOW_WITH_PETSC)
    target_compile_definitions(cf_common PUBLIC CHFLOW_WITH_PETSC)
    if(PETSc_INCLUDE_DIRS)
        target_include_directories(cf_common PUBLIC ${PETSc_INCLUDE_DIRS})
    elseif(PETSC_INCLUDE_DIRS)
        target_include_directories(cf_common PUBLIC ${PETSC_INCLUDE_DIRS})
    endif()
    target_link_libraries(cf_common PUBLIC ${PETSc_LIBRARIES})
endif()

add_subdirectory(apps)
