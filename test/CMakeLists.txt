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

add_library(solver_lib1 SHARED
            src/solver_ac_pose3d_1AC_plane.cpp
            )

add_library(solver_lib2 SHARED
            src/solver_ac_pose3d_2AC_plane.cpp
            )

add_library(solver_lib3 SHARED
            src/solver_ac_pose4d.cpp
            )

add_executable(test3d_inter_1AC_plane src/test3d_inter_1AC_plane.cpp)
add_executable(test3d_inter_2AC_plane src/test3d_inter_2AC_plane.cpp )
add_executable(test3d_intra_2AC_plane src/test3d_intra_2AC_plane.cpp)
add_executable(test4d_inter_2AC_vertical src/test4d_inter_2AC_vertical.cpp)
add_executable(test4d_intra_2AC_vertical src/test4d_intra_2AC_vertical.cpp)
target_link_libraries(test3d_inter_1AC_plane solver_lib1)
target_link_libraries(test3d_inter_2AC_plane solver_lib2)
target_link_libraries(test3d_intra_2AC_plane solver_lib2)
target_link_libraries(test4d_inter_2AC_vertical solver_lib3)
target_link_libraries(test4d_intra_2AC_vertical solver_lib3)