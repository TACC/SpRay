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

########################################
# paths
########################################
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

set(SPRAY_SRC_DIR ${CMAKE_SOURCE_DIR}/src)
include_directories(${CMAKE_SOURCE_DIR}/src ${CMAKE_BINARY_DIR})

########################################
# configure
########################################
set(SPRAY_CFG_RAY_DOMAIN_LIST_SIZE "16" CACHE STRING "Maximum number of domains a ray can intersect. A small value may incur a program crash.")
set(SPRAY_CFG_SPECU_HISTORY_SIZE "4" CACHE STRING "Speculation history size. A larage value may incur a memory explosion.")
set(SPRAY_CFG_SHADOW_TRANSPARENCY "0.3f" CACHE STRING "Shadow blending factor.")
set(SPRAY_CFG_L1_CACHE_LINE_SIZE "64" CACHE STRING "L1 cache size of target architecture.")

configure_file (
  "${SPRAY_SRC_DIR}/cmake_config.h.in"
  "${CMAKE_BINARY_DIR}/cmake_config.h" @ONLY
)

option(SPRAY_FILTER_ON "Enable filtered ray processing." ON)
if(SPRAY_FILTER_ON)
  add_definitions(-DSPRAY_FILTER_ON)
endif(SPRAY_FILTER_ON)

########################################
# compiler
########################################
if (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
  add_definitions(-DSPRAY_USING_CLANG_COMPILER)

  find_package(OpenMP REQUIRED)

  if (APPLE)
    include_directories(${OpenMP_INCLUDE_DIR})
    link_directories(${OpenMP_LIBRARY_DIR})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Xpreprocessor -fopenmp")
  else (APPLE)
    message(FATAL_ERROR "Clang not tested. Comment out this message at your discretion.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")
  endif (APPLE)

elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
  add_definitions(-DSPRAY_USING_INTEL_COMPILER)

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -qopt-report=5 -qopenmp")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -qopt-report=5 -qopenmp")

elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  add_definitions(-DSPRAY_USING_GNU_COMPILER)

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11 -fopenmp")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")

else ()
  MESSAGE(FATAL_ERROR "unsupported compiler: " ${CMAKE_CXX_COMPILER_ID})
endif ()

add_definitions(-D "__forceinline=inline")


########################################
# isa
########################################
if(${SPRAY_ISA_TYPE} STREQUAL "AVX2")
  add_definitions(-DSPRAY_AVX2)

elseif(${SPRAY_ISA_TYPE} STREQUAL "AVX512KNL")
  add_definitions(-DSPRAY_AVX512)

elseif(${SPRAY_ISA_TYPE} STREQUAL "AVX512SKX")
  add_definitions(-DSPRAY_AVX512)

else()
  message(FATAL_ERROR "Unsupported ISA type: " ${SPRAY_ISA_TYPE})
endif()

########################################
# mpi
########################################
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

########################################
# glog
########################################
set_property(TARGET glog PROPERTY FOLDER "dependencies")
include_directories (
  ${CMAKE_SOURCE_DIR}/deps/glog/src
  ${CMAKE_BINARY_DIR}/deps/glog/src
  ${CMAKE_BINARY_DIR}/deps/glog
)

########################################
# googletest
########################################
set_property(TARGET gtest PROPERTY FOLDER "dependencies")
set_property(TARGET gtest_main PROPERTY FOLDER "dependencies")
include_directories (${CMAKE_SOURCE_DIR}/deps/googletest/googletest/include)

########################################
# glm
########################################
include_directories(${CMAKE_SOURCE_DIR}/deps/glm)

########################################
# open gl
########################################
find_package(OpenGL REQUIRED)
if(APPLE)
  include_directories(${OPENGL_INCLUDE_DIR}/OpenGL)
else()
  include_directories(${OPENGL_INCLUDE_DIR}/GL)
endif()
message(STATUS "OPENGL_INCLUDE_DIR: " ${OPENGL_INCLUDE_DIR})

# WARNING! DO NOT EVER USE GLU
# message(STATUS "OPENGL_glu_LIBRARY: " ${OPENGL_glu_LIBRARY})

########################################
#glfw3
########################################
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL " " FORCE)

# add_subdirectory(${CMAKE_SOURCE_DIR}/deps/glfw)
set_property(TARGET glfw PROPERTY FOLDER "dependencies")

include_directories(${CMAKE_SOURCE_DIR}/deps/glfw/include)

########################################
# glut
########################################
find_package(GLUT REQUIRED)

########################################
# Embree
########################################
find_package(Embree REQUIRED)
include_directories(${EMBREE_INCLUDE_DIR})
link_directories(${EMBREE_LIBRARY_DIR})
message(STATUS "EMBREE_INSTALL_DIR: " ${EMBREE_INSTALL_DIR})

########################################
# SPRAY options
########################################

# glog
option(SPRAY_GLOG_CHECK "Enable glog check." OFF)
if(SPRAY_GLOG_CHECK)
  add_definitions(-DSPRAY_GLOG_CHECK)
endif(SPRAY_GLOG_CHECK)

# timing
option(SPRAY_TIMING "Enable timing for overall performance." ON)
if(SPRAY_TIMING)
  add_definitions(-DSPRAY_TIMING)
endif(SPRAY_TIMING)

########################################
# intersection mode (based on embree)
########################################
set(SPRAY_ISECT_MODE "PACKET1" CACHE STRING "Select mode for intersection tests.")
set_property(CACHE SPRAY_ISECT_MODE PROPERTY STRINGS PACKET1)

if (${SPRAY_ISECT_MODE} STREQUAL "PACKET1")
  add_definitions(-DSPRAY_ISECT_PACKET1)
else()
  message(FATAL_ERROR "Unsupported mode ${SPRAY_ISECT_MODE}")

endif()

########################################
# libraries
########################################
set(DEP_LIBS glfw glog gtest_main gtest pbrt
    ${OpenMP_LIBRARY}
    ${OPENGL_LIBRARIES}
    ${GLUT_LIBRARIES}
    ${TBB_LIBRARIES}
    ${EMBREE_LIBRARY}
    ${MPI_LIBRARIES}
    ${ICET_LIBRARIES})

########################################
# source files
########################################
set(SPRAY_SRC_LIST
    # utils
    utils/profiler.cc
    utils/comm.cc

    # scene
    scene/camera.cc

    # materials
    materials/reflection.cc

    # display
    display/composite.cc
    
    # caches
    caches/infinite_cache.cc
    caches/lru_cache.cc

    # accelerators
    accelerators/wbvh_embree.cc

    # partition
    partition/aabb.cc
    partition/block_buffer.cc
    partition/domain.cc
    partition/ply_loader.cc
    partition/tile.cc
    partition/trimesh_buffer.cc
    
    # renderers
    renderers/config.cc
    # renderers/rays.cc
    renderers/spray.cc
    renderers/sampler.cc
    # renderers/visibility_buffer.cc
    # renderers/scheduler.cc
    # renderers/async_comm_shared.cc
    # renderers/async_comm_insitu.cc
    # renderers/tile_context.cc

    # insitu
    insitu/insitu_comm.cc
    insitu/insitu_tiler.cc
    insitu/insitu_vbuf.cc
    insitu/insitu_work_stats.cc
    
    # ooc
    ooc/ooc_tiler.cc
    ooc/ooc_vbuf.cc
    ooc/ooc_domain_stats.cc
    ooc/ooc_pcontext.cc
    ooc/ooc_tcontext.cc

    # baseline
    baseline/baseline_ray.cc
    baseline/baseline_schedulers.cc
)

add_library(spray STATIC ${SPRAY_SRC_LIST})

add_executable(baseline_ooc apps/baseline_ooc.cc)
target_link_libraries(baseline_ooc spray ${DEP_LIBS})

add_executable(baseline_insitu apps/baseline_insitu.cc)
target_link_libraries(baseline_insitu spray ${DEP_LIBS})

add_executable(spray_insitu_singlethread apps/spray_insitu_singlethread.cc)
target_link_libraries(spray_insitu_singlethread spray ${DEP_LIBS})

add_executable(spray_insitu_multithread apps/spray_insitu_multithread.cc)
target_link_libraries(spray_insitu_multithread spray ${DEP_LIBS})

add_executable(spray_ooc apps/spray_ooc.cc)
target_link_libraries(spray_ooc spray ${DEP_LIBS})

# intallation
install (TARGETS baseline_ooc DESTINATION bin)
install (TARGETS spray_insitu_singlethread DESTINATION bin)
install (TARGETS spray_insitu_multithread DESTINATION bin)
install (TARGETS spray_ooc DESTINATION bin)
install (TARGETS spray DESTINATION lib)

