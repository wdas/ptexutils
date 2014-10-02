
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

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS
      "${CMAKE_CXX_FLAGS} -pedantic -Wall -Wextra -Wno-unused-parameter -std=c++0x -Wl,--allow-shlib-undefined")
endif()

include_directories(${PTEX_INCLUDE_DIR}
                    ${OPENGL_INCLUDE_DIR}
                    ${GLUT_INCLUDE_DIR}
                    ${OIIO_INCLUDE_DIR})

add_executable(ptxtransfer ptxtransfer.cpp)
target_link_libraries(ptxtransfer ${PTEX_LIBRARIES} ${ZLIB_LIBRARIES})

add_executable(ptxview ptxview.cpp)
target_link_libraries(ptxview ${PTEX_LIBRARIES}
                              ${GLUT_LIBRARIES}
                              ${OPENGL_LIBRARIES}
                              ${ZLIB_LIBRARIES})

add_executable(ptxconvert ptxconvert.cpp)
target_link_libraries(ptxconvert ${PTEX_LIBRARIES}
                                 ${OIIO_LIBRARIES}
                                 ${ZLIB_LIBRARIES})

install(TARGETS ptxtransfer DESTINATION ${BIN_INSTALL_DIR})

install(TARGETS ptxview DESTINATION ${BIN_INSTALL_DIR})

install(TARGETS ptxconvert DESTINATION ${BIN_INSTALL_DIR})
