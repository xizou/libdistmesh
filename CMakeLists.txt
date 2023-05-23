cmake_minimum_required(VERSION 3.0)
project(libdistmesh LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# set(CMAKE_COLOR_MAKEFILE ON)

if(WIN32)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/build") # .exe and .dll
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BUILD_DIR}/lib") # .so and .dylib
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BUILD_DIR}/lib") # .lib and .a
endif()

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

if(MSVC)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /Wall")
else()
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wextra")
endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_subdirectory(src)
add_subdirectory(examples)

message(STATUS "Using CMAKE_CXX_COMPILER = " ${CMAKE_CXX_COMPILER})
