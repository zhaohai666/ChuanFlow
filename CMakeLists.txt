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
add_subdirectory(src/pvt)
add_subdirectory(src/1d)
add_subdirectory(src/3d)

# Python bindings (README M1): pybind11 via FetchContent. Skipped when no
# interpreter development headers are available so the C++ build never breaks.
option(CHFLOW_WITH_PYTHON "Build the pybind11 Python bindings" ON)
if(CHFLOW_WITH_PYTHON)
    find_package(Python3 COMPONENTS Interpreter Development.Module QUIET)
    if(Python3_FOUND)
        include(FetchContent)
        FetchContent_Declare(pybind11
            GIT_REPOSITORY https://github.com/pybind/pybind11
            GIT_TAG v2.13.6
            GIT_SHALLOW TRUE)
        FetchContent_MakeAvailable(pybind11)
        add_subdirectory(src/bindings)
        set(CHFLOW_PYTHON_BINDINGS ON)
    else()
        message(STATUS "ChuanFlow: Python development headers not found - skipping bindings")
    endif()
endif()

add_subdirectory(tests)
