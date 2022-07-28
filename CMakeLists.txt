cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

project(ac_gcam_pose)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif ()
set( CMAKE_CXX_FLAGS "-std=c++11" )

find_package(Eigen3 REQUIRED)

include_directories(
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/src
        ${EIGEN3_INCLUDE_DIR}
)

add_library(solver_lib SHARED
            src/solver_ac_pose3d_1AC_plane.cpp
            src/solver_ac_pose3d_2AC_plane.cpp
            src/solver_ac_pose4d.cpp
            )
target_link_libraries(solver_lib)