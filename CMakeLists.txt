
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

cmake_minimum_required(VERSION 3.15.0)
project(ptexutils LANGUAGES CXX)

set(${PROJECT_NAME}_MAJOR_VERSION 0)
set(${PROJECT_NAME}_MINOR_VERSION 1)
set(${PROJECT_NAME}_PATCH_LEVEL 0)

set(CMAKE_CXX_EXTENSIONS OFF)
SET(CMAKE_VERBOSE_MAKEFILE TRUE)
set(CMAKE_SKIP_BUILD_RPATH TRUE)
set(CMAKE_FIND_PACKAGE_RESOLVE_SYMLINKS TRUE)
set(OpenGL_GL_PREFERENCE GLVND)

include(GNUInstallDirs)

## Setup platform specific helper defines build variants
if (WIN32)
  add_compile_definitions(-DPTEXUTILS_WIN32 -D_USE_MATH_DEFINES)
endif()

## Choose build options
# Disney specific method of choosing variant
if ("$ENV{FLAVOR}" STREQUAL "optimize")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "type of build" FORCE)
endif()

if ("$ENV{FLAVOR}" STREQUAL "debug")
    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "type of build" FORCE)
endif()

# Set to release if nothing else defined
if(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif()

if (DEFINED ENV{CXXFLAGS_STD})
    string(REPLACE "c++" "" cxx_std $ENV{CXXFLAGS_STD})
    set(CMAKE_CXX_STANDARD ${cxx_std})
endif()

## Set install location
if (NOT DEFINED CMAKE_INSTALL_PREFIX)
    execute_process(COMMAND sh -c "echo `uname`-`uname -r | cut -d'-' -f1`-`uname -m`"
                    OUTPUT_VARIABLE VARIANT_DIRECTORY
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/${VARIANT_DIRECTORY}")
endif()

if (DEFINED ENV{RP_glew})
    set(GLEW_ROOT $ENV{RP_glew} CACHE PATH "glew location" FORCE)
endif()

if (DEFINED ENV{RP_freeglut})
    set(GLUT_ROOT $ENV{RP_freeglut} CACHE PATH "glut location" FORCE)
endif()


if (DEFINED ENV{RP_mesa})
    set(OPENGL_INCLUDE_DIR "$ENV{RP_mesa}/include" CACHE PATH "OpenGL headers" FORCE)
    set(OPENGL_EGL_INCLUDE_DIR "$ENV{RP_mesa}/include" CACHE PATH "OpenGL EGL headers" FORCE)
    set(OPENGL_GL_INCLUDE_DIR "$ENV{RP_mesa}/include" CACHE PATH "OpenGL GL headers" FORCE)
    set(OPENGL_GLX_INCLUDE_DIR "$ENV{RP_mesa}/include" CACHE PATH "OpenGL GLX headers" FORCE)
endif()
if (DEFINED ENV{RP_zlib})
    set(ZLIB_ROOT $ENV{RP_zlib} CACHE PATH "zlib location" FORCE)
endif()

## Set the cmake module path
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

## Search for useful libraries
find_package(Ptex REQUIRED)
find_package(OpenMP REQUIRED)
find_package(GLUT REQUIRED)
find_package(GLEW REQUIRED)
find_package(OpenGL REQUIRED)
find_package(OIIO REQUIRED)

## Trasitive dependency, needed when linking with static libPTex
find_package(ZLIB REQUIRED)

## Use static ptex library
option(PTEX_USE_STATIC_LIBS "Link against Ptex static libraries" ON)
if (PTEX_USE_STATIC_LIBS)
    set(PTEX_LIBRARIES Ptex::Ptex_static)
else()
    set(PTEX_LIBRARIES Ptex::Ptex_dynamic)
endif()

## Traverse subdirectories
add_subdirectory(src)
