
# Copyright Disney Enterprises, Inc.  All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License
# and the following modification to it: Section 6 Trademarks.
# deleted and replaced with:
#
# 6. Trademarks. This License does not grant permission to use the
# trade names, trademarks, service marks, or product names of the
# Licensor and its affiliates, except as required for reproducing
# the content of the NOTICE file.
#
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0


## CMake compatibility issues: don't modify this, please!
CMAKE_MINIMUM_REQUIRED( VERSION 2.4.6 )
MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)
## allow more human readable "if then else" constructs
SET( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE )
## Use verbose make file
SET( CMAKE_VERBOSE_MAKEFILE TRUE )
set(CMAKE_SKIP_BUILD_RPATH TRUE)

## project name & version
PROJECT( ptexutils )
SET( ${PROJECT_NAME}_MAJOR_VERSION 0 )
SET( ${PROJECT_NAME}_MINOR_VERSION 1 )
SET( ${PROJECT_NAME}_PATCH_LEVEL 0 )

set(BIN_INSTALL_DIR bin CACHE STRING "Executable install directory")

# output executable and library path
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

## Setup platform specific helper defines build variants
IF(WIN32)
  ADD_DEFINITIONS (-DPTEXUTILS_WIN32 -D_USE_MATH_DEFINES)
ENDIF(WIN32)

# Adjust OpenMP in Clang on Apple platforms
IF(APPLE)
    IF(CMAKE_C_COMPILER_ID MATCHES "Clang")
        SET(OpenMP_C "${CMAKE_C_COMPILER}" CACHE STRING "" FORCE)
        SET(OpenMP_C_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument" CACHE STRING "" FORCE)
        SET(OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5" CACHE STRING "" FORCE)
        SET(OpenMP_libomp_LIBRARY ${OpenMP_C_LIB_NAMES} CACHE STRING "" FORCE)
        SET(OpenMP_libgomp_LIBRARY ${OpenMP_C_LIB_NAMES} CACHE STRING "" FORCE)
        SET(OpenMP_libiomp5_LIBRARY ${OpenMP_C_LIB_NAMES} CACHE STRING "" FORCE)
    ENDIF()

    IF(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        SET(OpenMP_CXX "${CMAKE_CXX_COMPILER}" CACHE STRING "" FORCE)
        SET(OpenMP_CXX_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument" CACHE STRING "" FORCE)
        SET(OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5" CACHE STRING "" FORCE)
        SET(OpenMP_libomp_LIBRARY ${OpenMP_CXX_LIB_NAMES} CACHE STRING "" FORCE)
        SET(OpenMP_libgomp_LIBRARY ${OpenMP_CXX_LIB_NAMES} CACHE STRING "" FORCE)
        SET(OpenMP_libiomp5_LIBRARY ${OpenMP_CXX_LIB_NAMES} CACHE STRING "" FORCE)
    ENDIF()

    # Default (Command-line tools & Xcode toolchain) compiler won't provide proper OpenMP functionality on macOS.
    # Compiler obtained from `brew install llvm` will suffice. 
    # Create a symbolic link for ease of use.
    SET(CMAKE_CXX_COMPILER clang-omp++ CACHE STRING "C++ compiler" FORCE)
ENDIF()


## Choose build options
# Disney specific method of choosing variant
IF("$ENV{FLAVOR}" STREQUAL "optimize")
    SET(CMAKE_BUILD_TYPE "Release" CACHE STRING "type of build" FORCE)
ENDIF("$ENV{FLAVOR}" STREQUAL "optimize")
IF("$ENV{FLAVOR}" STREQUAL "debug")
    SET(CMAKE_BUILD_TYPE "Debug" CACHE STRING "type of build" FORCE)
ENDIF("$ENV{FLAVOR}" STREQUAL "debug")
# Set to release if nothing else defined
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

## Set install location
IF (NOT DEFINED CMAKE_INSTALL_PREFIX)
    EXECUTE_PROCESS(COMMAND sh -c "echo `uname`-`uname -r | cut -d'-' -f1`-`uname -m`"
                    OUTPUT_VARIABLE VARIANT_DIRECTORY
		    OUTPUT_STRIP_TRAILING_WHITESPACE)
    SET(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/${VARIANT_DIRECTORY}")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)

## Set the cmake module path
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

## Search for useful libraries
find_package(PTex REQUIRED)
find_package(OpenMP REQUIRED)
find_package(GLUT REQUIRED)
find_package(OpenGL REQUIRED)
find_package(OIIO REQUIRED)

## Trasitive dependency, needed when linking with static libPTex
find_package(ZLIB REQUIRED)

#INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/src/lib )

## Traverse subdirectories
ADD_SUBDIRECTORY (src)
