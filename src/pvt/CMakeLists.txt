add_library(cf_pvt STATIC
    src/components.cpp
    src/pengRobinson.cpp
    src/flash.cpp
    src/pvtTable.cpp
    src/pvtCase.cpp
)

target_include_directories(cf_pvt PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(cf_pvt PUBLIC cf_common)

add_subdirectory(apps)
