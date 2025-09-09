cmake_minimum_required(VERSION 3.9)
project(StiefelManifoldExample)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release)

if(POLICY CMP0167)
    cmake_policy(SET CMP0167 OLD)
endif()

#set(GTSAM_DIR /home/aryaman/gtsam/cmake-build-release)
#set(GTSAM_UNSTABLE_DIR /home/aryaman/gtsam/cmake-build-release)
#find_package(SPQR REQUIRED)
find_package(BLAS REQUIRED)
include_directories(/usr/include/suitesparse)
#set(Eigen3_DIR "/home/jason/SLAM_LIB/eigen-3.4.0/install/share/eigen3/cmake")
# Require *at least* 3.4.0 and force use of the Eigen3::Eigen imported target.
#find_package(Eigen3 3.4.0 REQUIRED NO_MODULE)
#set(Eigen3_DIR "/home/jason/SLAM_LIB/eigen-3.4.0/install/share/eigen3/cmake")
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
# Directory for built libraries
set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/lib CACHE PATH "The directory in which to place the SE-Sync library built by this project")
# Directory for built executables
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/bin CACHE PATH "The directory in which to place executables built by this project")
find_package(Eigen3 3.3.7 REQUIRED)
if(EIGEN3_FOUND)
    message(STATUS "Found Eigen3 library: ${EIGEN3_INCLUDE_DIR}")
endif()


# Find GTSAM components
set(GTSAM_DIR /home/nikolas/gtsam/build)
set(GTSAM_UNSTABLE_DIR /home/nikolas/gtsam/build)
set(GTSAM_MIN_VERSION "4.3.0")

# Make sure CMake searches these directories first
list(APPEND CMAKE_PREFIX_PATH "${GTSAM_DIR}" "${GTSAM_UNSTABLE_DIR}")

# require GTSAM 4.3.0 at least and use GTSAM_DIR to find the package
find_package(GTSAM REQUIRED ${GTSAM_MIN_VERSION} PATHS ${GTSAM_DIR}) # Uses installed package
if(GTSAM_FOUND)
    message(STATUS "Found GTSAM library: ${GTSAM_INCLUDE_DIR}\n")
endif()
include_directories(${GTSAM_INCLUDE_DIR})

find_package(GTSAM_UNSTABLE REQUIRED ${GTSAM_MIN_VERSION} PATHS ${GTSAM_UNSTABLE_DIR}) # Uses installed package
if(GTSAM_UNSTABLE_FOUND)
    message(STATUS "Found GTSAM_UNSTABLE library: ${GTSAM_INCLUDE_DIR}\n")
endif()

set(GTSAM_LIBRARIES gtsam)   # TODO: automatic search libs
set(GTSAM_UNSTABLE_LIBRARIES gtsam_unstable)   # TODO: automatic search libs

find_package(GTSAMCMakeTools)
if(GTSAMCMakeTools_FOUND)
    message(STATUS "Found GTSAMCMakeTools library: ${GTSAMCMakeTools_DIR}\n")
endif()
find_package(OpenMP REQUIRED)
if(OpenMP_FOUND)
    message(STATUS "Found OpenMP library: ${OpenMP_DIR}\n")
endif()
set(EIGEN_USE_BLAS ON)
find_library(PROFILER_LIB profiler REQUIRED)
include(GtsamMakeConfigFile)
include(GtsamBuildTypes)
include(GtsamTesting)

# #–– gflags
# find_package(gflags REQUIRED)
# if (gflags_FOUND)
#     message(STATUS "gflags: found (version ${gflags_VERSION})")
# else()
#     message(WARNING "gflags: NOT found")
# endif()

# #–– glog
# find_package(glog REQUIRED)
# if (glog_FOUND)
#     message(STATUS "glog: found (version ${glog_VERSION})")
# else()
#     message(WARNING "glog: NOT found")
# endif()

# #–– GoogleTest
# find_package(GTest REQUIRED)
# if (GTest_FOUND)
#     message(STATUS "GTest: found (version ${GTest_VERSION})")
# else()
#     message(WARNING "GTest: NOT found")
# endif()

# for unittest scripts
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${GTSAM_DIR}/../GTSAMCMakeTools")

# Process source subdirs
add_subdirectory(cpp)


#From fast verifications
# Find Optimization library
set(OPTIMIZATION_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Optimization/ CACHE PATH "Path to top-level directory of Optimization library (the one containing CMakeLists.txt)")

# Run the Optimization package's CMake file
add_subdirectory(${OPTIMIZATION_DIR})

# Find Preconditioners library
set(PRECONDITIONERS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Preconditioners/ CACHE PATH "Path to top-level directory of Preconditioners library (the one containing CMakeLists.txt)")

# Run the Preconditioners package's CMake file
add_subdirectory(${PRECONDITIONERS_DIR})

