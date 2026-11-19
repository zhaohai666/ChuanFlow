# Sparse linear solver backend verification app (README M6 / A11).
add_executable(chflow-solver-backend solverBackendApp.cpp)
target_link_libraries(chflow-solver-backend PRIVATE cf_common)
