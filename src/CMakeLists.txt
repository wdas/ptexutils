
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

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -Wall -pedantic")

add_executable(ptxtransfer ptxtransfer.cpp)
target_link_libraries(ptxtransfer ${PTEX_LIBRARIES} gomp)

add_executable(ptxview ptxview.cpp)
target_link_libraries(ptxview ${PTEX_LIBRARIES} glut GLU)

install(TARGETS ptxtransfer DESTINATION bin)

install(TARGETS ptxview DESTINATION bin)
