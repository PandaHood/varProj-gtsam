# Install script for directory: /home/nikolas/varProj-gtsam/cpp

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/StiefelManifoldExample" TYPE FILE FILES
    "/home/nikolas/varProj-gtsam/cpp/CertifiableLandmark.h"
    "/home/nikolas/varProj-gtsam/cpp/CertifiablePGO.h"
    "/home/nikolas/varProj-gtsam/cpp/CertifiableProblemOpts.h"
    "/home/nikolas/varProj-gtsam/cpp/CertifiableRA.h"
    "/home/nikolas/varProj-gtsam/cpp/CertifiableRangeSLAM.h"
    "/home/nikolas/varProj-gtsam/cpp/CertifiableResults.h"
    "/home/nikolas/varProj-gtsam/cpp/Certifiable_Incremental_PGO.h"
    "/home/nikolas/varProj-gtsam/cpp/Certifiable_problem.h"
    "/home/nikolas/varProj-gtsam/cpp/ConcurrentCertifiableBatchSmoother.h"
    "/home/nikolas/varProj-gtsam/cpp/LandmarkFactor.h"
    "/home/nikolas/varProj-gtsam/cpp/LiftedPose.h"
    "/home/nikolas/varProj-gtsam/cpp/LiftedPosePriorFactor.h"
    "/home/nikolas/varProj-gtsam/cpp/LiftedRangeFactor.h"
    "/home/nikolas/varProj-gtsam/cpp/NormalRangeFactor.h"
    "/home/nikolas/varProj-gtsam/cpp/RaFactor.h"
    "/home/nikolas/varProj-gtsam/cpp/RelativePoseMeasurement.h"
    "/home/nikolas/varProj-gtsam/cpp/RunSAM.h"
    "/home/nikolas/varProj-gtsam/cpp/SEsyncFactor.h"
    "/home/nikolas/varProj-gtsam/cpp/StiefelManifold-inl.h"
    "/home/nikolas/varProj-gtsam/cpp/StiefelManifold.h"
    "/home/nikolas/varProj-gtsam/cpp/StiefelManifoldPriorFactor.h"
    "/home/nikolas/varProj-gtsam/cpp/UnitSphere.h"
    "/home/nikolas/varProj-gtsam/cpp/types.h"
    "/home/nikolas/varProj-gtsam/cpp/utils.h"
    )
endif()

