add_library(cfcoupling STATIC
    src/coupledNetwork.cpp
)

target_include_directories(cfcoupling PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(cfcoupling PUBLIC cf1d cf3d)

add_subdirectory(apps)
