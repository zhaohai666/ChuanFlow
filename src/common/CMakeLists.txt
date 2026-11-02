add_library(cf_common STATIC
    src/Dictionary.cpp
    src/Matrix.cpp
)

target_include_directories(cf_common PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
