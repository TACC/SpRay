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

#include "render/shape.h"
#include "render/rays.h"

namespace spray {

struct RTCRayIntersection;
struct Domain;
class Material;

class ShapeBuffer {
 public:
  ShapeBuffer();
  ~ShapeBuffer();

 public:
  void init(int max_cache_size_ndomains, std::size_t max_nvertices,
            std::size_t max_nfaces);

  RTCScene load(int cache_block, Domain& domain);

 private:
  void loadShapes(std::vector<Shape*>& shapes, int cache_block);

 public:
  RTCScene get(int cache_block) { return scenes_[cache_block]; }

  void updateIntersection(int cache_block, RTCRayIntersection* isect) const {
    isect->material = getMaterial(cache_block, isect->primID);
  }

  static void intersect(void* shape_ptr, RTCRay& ray_i);

 private:
  // void getColorTuple(int cache_block, uint32_t primID,
  //                    uint32_t colors[3]) const;

  void cleanup();

  static void sphereBoundsCallback(void* shape_ptr, std::size_t item,
                                   RTCBounds& bounds_o);
  static void sphereIntersect1Callback(void* shape_ptr, RTCRay& ray_i,
                                       std::size_t item);
  static void sphereOccluded1Callback(void* shape_ptr, RTCRay& ray_i,
                                      std::size_t item);

  Material* getMaterial(int cache_block, int prim_id) const {
    Shape* shape = shapes_[cache_block]->at(prim_id);
    return shape->material;
  }

 private:
  enum MeshStatus { CREATED = -1, DESTROYED = 0 };

 private:
  int max_cache_size_;  // in number of domains

  // uint32_t* colors_;  //!< per-cache-block packed rgb colors. 2d array
  // organized as 1d

  RTCDevice device_;

  std::vector<RTCScene>
      scenes_;  //!< 1D array of per-cache-block Embree scenes.

  // MemoryArena arena_;

  bool loaded_;
  std::vector<const std::vector<Shape*>*> shapes_;
};

}  // namespace spray

