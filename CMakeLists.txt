cmake_minimum_required(VERSION 3.8)
project(camera_calibration)


set(CMAKE_CXX_STANDARD 14)

IF(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Debug Mode")
    set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/pipeline/bin_debug)
    #"Valgrind erf_inv_initializer crashes" workaround
    set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-DBOOST_MATH_PROMOTE_DOUBLE_POLICY=false -g3")
    set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} "-g3")
ENDIF()

IF(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/pipeline/bin)
ENDIF()

add_subdirectory(subroutines)
add_subdirectory(core)

