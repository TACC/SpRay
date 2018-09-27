# ========================================================================== #
# Copyright (c) 2017-2018 The University of Texas at Austin.                 #
# All rights reserved.                                                       #
#                                                                            #
# Licensed under the Apache License, Version 2.0 (the "License");            #
# you may not use this file except in compliance with the License.           #
# A copy of the License is included with this software in the file LICENSE.  #
# If your copy does not contain the License, you may obtain a copy of the    #
# License at:                                                                #
#                                                                            #
#     https://www.apache.org/licenses/LICENSE-2.0                            #
#                                                                            #
# Unless required by applicable law or agreed to in writing, software        #
# distributed under the License is distributed on an "AS IS" BASIS, WITHOUT  #
# WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.           #
# See the License for the specific language governing permissions and        #
# limitations under the License.                                             #
#                                                                            #
# ========================================================================== #

cmake_minimum_required(VERSION 3.5.1)

########################################
# project info
########################################
project(spray)

set(spray_VERSION_MAJOR   1   CACHE STRING "Project major version number.")
set(spray_VERSION_MINOR   0   CACHE STRING "Project minor version number.")
set(spray_VERSION_PATCH   0   CACHE STRING "Project patch version number.")

MESSAGE(STATUS "PROJECT_SOURCE_DIR: " ${PROJECT_SOURCE_DIR})

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

SET(HOST_NAME $ENV{HOSTNAME})
message(STATUS "HOST_NAME: " ${HOST_NAME})

########################################
# isa
########################################
set(SPRAY_ISA_TYPE "AVX2" CACHE STRING "Select ISA type.")
set_property(CACHE SPRAY_ISA_TYPE PROPERTY STRINGS AVX2)

if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
  if(${SPRAY_ISA_TYPE} STREQUAL "AVX2")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xCORE-AVX2")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -xCORE-AVX2")

  elseif(${SPRAY_ISA_TYPE} STREQUAL "AVX512KNL")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xMIC-AVX512")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -xMIC-AVX512")

  elseif(${SPRAY_ISA_TYPE} STREQUAL "AVX512SKX")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xCORE-AVX512")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -xCORE-AVX512")

  else()
    message(FATAL_ERROR "Unsupported ISA type: " ${SPRAY_ISA_TYPE})
  endif()
endif()

########################################
# dependencies
########################################
add_subdirectory(${CMAKE_SOURCE_DIR}/deps/icet)
include_directories(${CMAKE_SOURCE_DIR}/deps/icet/src/include)
include_directories(${CMAKE_BINARY_DIR}/deps/icet/src/include)

include_directories(${CMAKE_SOURCE_DIR}/deps)

add_subdirectory(${CMAKE_SOURCE_DIR}/deps/glfw)

include(pbrt_compile)
include(pbrt_flags)
add_subdirectory(${CMAKE_SOURCE_DIR}/deps/pbrt)

set(WITH_GFLAGS OFF CACHE BOOL "Use gflags")
set(BUILD_SHARED_LIBS OFF CACHE BOOL " " FORCE)
add_subdirectory(${CMAKE_SOURCE_DIR}/deps/glog)

add_subdirectory(${CMAKE_SOURCE_DIR}/deps/googletest)

########################################
# main source files
########################################
add_subdirectory(src)

