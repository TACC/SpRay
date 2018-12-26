// ========================================================================== //
// Copyright (c) 2017-2018 The University of Texas at Austin.                 //
// All rights reserved.                                                       //
//                                                                            //
// Licensed under the Apache License, Version 2.0 (the "License");            //
// you may not use this file except in compliance with the License.           //
// A copy of the License is included with this software in the file LICENSE.  //
// If your copy does not contain the License, you may obtain a copy of the    //
// License at:                                                                //
//                                                                            //
//     https://www.apache.org/licenses/LICENSE-2.0                            //
//                                                                            //
// Unless required by applicable law or agreed to in writing, software        //
// distributed under the License is distributed on an "AS IS" BASIS, WITHOUT  //
// WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.           //
// See the License for the specific language governing permissions and        //
// limitations under the License.                                             //
//                                                                            //
// ========================================================================== //

#pragma once

#include <cstdint>

#include <embree2/rtcore_ray.h>
#include <embree2/rtcore_scene.h>

#include "glm/glm.hpp"
#include "pbrt/memory.h"

#include "io/ply_loader.h"

#define NUM_VERTICES_PER_FACE 3  // triangle

namespace spray {

struct RTCRayIntersection;

class TriMeshBuffer {
 public:
  TriMeshBuffer();
  ~TriMeshBuffer();

 public:
  void init(int max_cache_size_ndomains, std::size_t max_nvertices,
            std::size_t max_nfaces, bool compute_normals);

  RTCScene load(const std::string& filename, int cache_block,
                const glm::mat4& transform, bool apply_transform);
  RTCScene get(int cache_block) { return scenes_[cache_block]; }

  void updateIntersection(int cache_block, RTCRayIntersection* isect) const;

 private:
  void getColorTuple(int cache_block, uint32_t primID,
                     uint32_t colors[3]) const;

  void getNormalTuple(int cache_block, uint32_t primID,
                      float normals_out[9]) const;

  void computeNormals(int cache_block);

  std::size_t vertexBaseIndex(int cache_block) const {
    return (cache_block * max_nvertices_ * 3);
  }

  std::size_t normalBaseIndex(int cache_block) const {
    return (cache_block * max_nvertices_ * 3);
  }

  std::size_t faceBaseIndex(int cache_block) const {
    return (cache_block * max_nfaces_ * NUM_VERTICES_PER_FACE);
  }

  std::size_t colorBaseIndex(int cache_block) const {
    return (cache_block * max_nvertices_);
  }

  void cleanup();
  void mapEmbreeBuffer(int cache_block, float* vertices,
                       std::size_t num_vertices, uint32_t* faces,
                       std::size_t num_faces);

 private:
  enum MeshStatus { CREATED = -1, DESTROYED = 0 };

 private:
  int max_cache_size_;  // in number of domains
  std::size_t max_nvertices_;
  std::size_t max_nfaces_;

  float* vertices_;  //!< per-cache-block vertices. 1d array as 2d.
  float* normals_;   //!< per-cache-block normals. unnormalized. 1d array as 2d.
  uint32_t* faces_;  //!< per-cache-block faces. 1d array as 2d.
  uint32_t* colors_;  //!< per-cache-block packed rgb colors. 1d array as 2d.

  std::size_t* num_vertices_;
  std::size_t* num_faces_;

  RTCDevice device_;
  RTCScene* scenes_;  //!< 1D array of per-cache-block Embree scenes.

  int* embree_mesh_created_;  // -1: initialized, 0: not initialized

  MemoryArena arena_;
  PlyLoader loader_;

  bool compute_normals_;
};

}  // namespace spray

