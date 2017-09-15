cmake_minimum_required(VERSION 3.7)
project(arc)

# The version number.
set(ARC_MAJOR_VERSION 1)
set(ARC_MINOR_VERSION 1)
set(ARC_PATCH_VERSION 1)
set(ARC_VERSION
        ${ARC_MAJOR_VERSION}.${ARC_MINOR_VERSION}.${ARC_PATCH_VERSION})

# Offer the user the choice of overriding the installation directories
set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_BUILD_TYPE "Debug")
set(CMAKE_CXX_FLAGS "-std=c++11")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -g")

MESSAGE("Build type: " ${CMAKE_BUILD_TYPE})

set(QT_USE_QTOPENGL TRUE)

find_package(Qt4 REQUIRED)

include(${QT_USE_FILE})
add_definitions(${QT_DEFINITIONS})

# Check C++11 or C++0x support
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    add_definitions(-DCOMPILEDWITHC11)
    message(STATUS "Using flag -std=c++11.")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
    add_definitions(-DCOMPILEDWITHC0X)
    message(STATUS "Using flag -std=c++0x.")
else()
    message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no
                         C++11 support. Please use a different C++ compiler.")
endif()
add_definitions(-w)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

if (WITH_GFLAGS)
    find_package (gflags 2.2.0 REQUIRED)
    if (gflags_FOUND)
        set (HAVE_LIB_GFLAGS 1)
    endif (gflags_FOUND)
endif (WITH_GFLAGS)

include_directories(
    ${PROJECT_SOURCE_DIR}/arc_log/include
    ${PROJECT_SOURCE_DIR}/arc/include
    ${GFLAGS_INCLUDE_DIR}
    ${PROJECT_SOURCE_DIR}/arc_pf/lib/include
    /usr/include/eigen3
    /usr/include/qt4/Qt
    /usr/include/qt4
    /usr/include/qt4/QtGui
    /usr/include/qt4/QtCore
    ${PROJECT_SOURCE_DIR}/arc_plot/include
    ${PROJECT_SOURCE_DIR}/arc_cmn/src/newmat
)

add_subdirectory(arc_log)
add_subdirectory(arc_pf)
add_subdirectory(arc_plot)
add_subdirectory(arc_cmn/src/newmat)

add_library(${PROJECT_NAME}_rtk STATIC
    arc_cmn/src/arc_cmn.cc
    arc_eph/src/arc_eph.cc
    arc_rinex/src/arc_rinex.cc
    arc_srtk/src/arc_srtk.cc
    arc_srtk/src/arc_pnt.cc
    arc_cmn/src/arc_lambda.cc
    arc/src/arc_core.cc
    arc/src/arc_States.cc
    arc/src/arc_MovementModel.cc
    arc/src/arc_ObservationModel.cpp
    arc/src/arc_StatesDistribution.cpp
    arc/src/arc_ukf.cc
    arc/src/arc_ffratio_table.cc
    arc_cmn/src/arc_opt.cc
    arc_cmn/src/arc_solution.cc
    arc_srtk/src/arc_srtk_dd.cc
    arc_srtk/src/arc_srtkpos.cc)

target_link_libraries(${PROJECT_NAME}_rtk glog PF ceres)

add_executable(arc_test1
               arc_test/src/test1.cpp arc_test/src/huace_test.cpp)
add_executable(arc_test2
               arc_test/src/test2.cpp)

add_executable(ARC-SRTK
               arc/src/arc_app.cc)

add_executable(huace_test
              arc_test/src/huace_test.cpp)

target_link_libraries(arc_test1 ${PROJECT_NAME}_rtk )
target_link_libraries(arc_test2 ${PROJECT_NAME}_rtk )
target_link_libraries(ARC-SRTK  ${PROJECT_NAME}_rtk QtGui QtCore)
target_link_libraries(huace_test ${PROJECT_NAME}_rtk)


