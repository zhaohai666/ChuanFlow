cmake_minimum_required(VERSION 3.24)

project(ChuanFlow
    VERSION 0.1.0
    DESCRIPTION "ChuanFlow: integrated multiphase flow simulation platform (1D network + 1D transient + 3D FVM)"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_compile_options(-Wall -Wextra -Wno-unused-parameter)

enable_testing()

add_subdirectory(src/common)
add_subdirectory(src/3d)
add_subdirectory(tests)
