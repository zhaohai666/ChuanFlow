add_library(cf1d STATIC
    src/network.cpp
    src/transient.cpp
    src/beggsBrill.cpp
    src/twoFluid.cpp
    src/nodal.cpp
)

target_include_directories(cf1d PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(cf1d PUBLIC cf_common cf_pvt)

add_subdirectory(apps)
