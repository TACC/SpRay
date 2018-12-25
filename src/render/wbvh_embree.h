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

#include <embree2/rtcore.h>

#include "glog/logging.h"

#include "partition/aabb.h"
#include "partition/domain.h"
#include "renderers/rays.h"

namespace spray {

class WbvhNode;

struct WbvhEmbreePrimitive {
  glm::vec3 bounds_min;
  unsigned id;
  glm::vec3 bounds_max;
  unsigned nprims;
};

class WbvhEmbree {
 public:
  enum BuildMode { NORMAL1, NORMAL8, NORMAL16, STREAM_1M };

  WbvhEmbree() : device_(nullptr), scene_(nullptr) {}
  ~WbvhEmbree() { cleanup(); }

  void initialize(const Aabb &bound, const std::vector<Domain> &domains);

  void runBuilder() { LOG(FATAL) << "unsupported method"; }

  void build(BuildMode mode);

  // APIs for intersection tests

  void intersect(RTCRayExt &ray) { rtcIntersect(scene_, (RTCRay &)ray); }

  void intersect8(const unsigned valid[8], RTCRayExt8 &ray) {
    rtcIntersect8((const void *)valid, scene_, (RTCRay8 &)ray);
  }

  void intersect16(const void *valid, RTCRay16 &ray) {
    rtcIntersect16(valid, scene_, ray);
  }

  WbvhNode *getRoot() { return nullptr; }

 private:
  void cleanup() {
    if (scene_) {
      rtcDeleteScene(scene_);
    }
    if (device_) {
      rtcDeleteDevice(device_);
    }
  }

  static void cbBounds(void *ptr,   /*!< pointer to user data */
                       size_t item, /*!< item to calculate bounds for */
                       RTCBounds &bounds_o /*!< returns calculated bounds */);

  static void cbIntersect1(void *ptr,     /*!< pointer to user data */
                           RTCRay &ray_i, /*!< ray to intersect */
                           size_t item /*!< item to intersect */);

  static void cbIntersect8(const void *valid, /*!< pointer to valid mask */
                           void *ptr,         /*!< pointer to user data */
                           RTCRay8 &ray,      /*!< ray packet to intersect */
                           size_t item /*!< item to intersect */);

  static void cbIntersect16(const void *valid, /*!< pointer to valid mask */
                            void *ptr,         /*!< pointer to user data */
                            RTCRay16 &ray,     /*!< ray packet to intersect */
                            size_t item /*!< item to intersect */);

  static void cbIntersectStream1M(
      const int *valid,                   /*!< pointer to valid mask */
      void *ptr,                          /*!< pointer to geometry user data */
      const RTCIntersectContext *context, /*!< intersection context as passed to
                                             rtcIntersect/rtcOccluded */
      RTCRayN *rays,                      /*!< ray packet to intersect */
      size_t N,                           /*!< number of rays in packet */
      size_t item /*!< item to intersect */);

 private:
  Aabb bound_;
  std::vector<WbvhEmbreePrimitive> prims_;

  RTCDevice device_;
  RTCScene scene_;
};

}  // namespace spray
