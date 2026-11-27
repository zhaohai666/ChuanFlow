pybind11_add_module(chflow_core
    src/module.cpp
    src/coupling_bind.cpp
    src/coupling_transient_bind.cpp
    src/scalar_transport_bind.cpp
    src/melting_bind.cpp
    src/combustion_bind.cpp
    src/vof_bind.cpp
    src/dpm_bind.cpp
    src/pbm_bind.cpp
    src/qmom_bind.cpp
    src/les_bind.cpp
    src/rsm_bind.cpp
    src/radiation_do_bind.cpp
    src/ale_bind.cpp
    src/dynmesh_bind.cpp
    src/ns_bind.cpp
    src/fields3d_bind.cpp)
# cfcoupling pulls cf1d + cf3d transitively (PUBLIC link), which is what the
# coupling and 3D binding families need.
target_link_libraries(chflow_core PRIVATE cf1d cf_pvt cfcoupling)

# NumPy headers for the Field <-> ndarray conversion layer. The interpreter
# here is the same one CTest runs (PYTHON_EXECUTABLE forced by the parent).
execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -c
    "import numpy, sys; sys.stdout.write(numpy.get_include())"
    OUTPUT_VARIABLE CF_NUMPY_INC
    ERROR_QUIET)
if(CF_NUMPY_INC)
    target_include_directories(chflow_core PRIVATE ${CF_NUMPY_INC})
endif()
