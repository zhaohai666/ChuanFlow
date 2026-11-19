add_library(cf_pvt STATIC
    src/components.cpp
    src/componentData.cpp
    src/pengRobinson.cpp
    src/flash.cpp
    src/pvtTable.cpp
    src/pvtCase.cpp
    src/blackOil.cpp
    src/separator.cpp
    src/viscosity.cpp
    src/surfaceTension.cpp
    src/waxSolubility.cpp
    src/thermalConductivity.cpp
    src/hydrateEquilibrium.cpp
)

target_include_directories(cf_pvt PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(cf_pvt PUBLIC cf_common)

add_subdirectory(apps)
