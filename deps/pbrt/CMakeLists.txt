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

project(pbrt-v3)
cmake_minimum_required(VERSION 2.8.6)

# include(pbrt_flags)
# include(pbrt_compile)
set(PBRT_SOURCE_LIST ${CMAKE_CURRENT_SOURCE_DIR}/memory.cpp)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

include_directories (${CMAKE_SOURCE_DIR}/deps/glog/src
                     ${CMAKE_BINARY_DIR}/deps/glog/src
                     ${CMAKE_BINARY_DIR}/deps/glog)

add_library(pbrt ${PBRT_SOURCE_LIST})
