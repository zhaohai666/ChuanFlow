pybind11_add_module(chflow_core src/module.cpp)
target_link_libraries(chflow_core PRIVATE cf1d cf_pvt)
