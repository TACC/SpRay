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

#include "render/wbvh_embree.h"

#include "render/aabb.h"
#include "render/rays.h"

#define PRINT_WBVH_BOUNDS
#undef PRINT_WBVH_BOUNDS

namespace spray {

void WbvhEmbree::init(const Aabb& bound, const std::vector<Domain>& domains) {
  // populate wbvh primitives
  prims_.resize(domains.size());

  Aabb aabb;
  for (std::size_t i = 0; i < domains.size(); ++i) {
    prims_[i].id = i;
    prims_[i].nprims = domains[i].num_faces;
    aabb = domains[i].world_aabb;

    prims_[i].bounds_min = aabb.bounds[0];
    prims_[i].bounds_max = aabb.bounds[1];

#ifdef PRINT_WBVH_BOUNDS
    std::cout << "[WbvhDmerge::init] [domain " << i
              << "] [num_prims: " << prims_[i].nprims << "] [bound: " << aabb
              << "]" << std::endl;
#endif
  }

  bound_ = bound;

  cleanup();
  device_ = nullptr;
  scene_ = nullptr;
}

void WbvhEmbree::build(BuildMode mode) {
  if (!device_) {
    device_ = rtcNewDevice(nullptr);
  }

  RTCSceneFlags sflags = RTC_SCENE_STATIC | RTC_SCENE_HIGH_QUALITY;

  RTCAlgorithmFlags aflags;

  CHECK(mode == NORMAL1 || mode == NORMAL8 || mode == NORMAL16)
      << "unsupported build mode: " << mode;

  if (mode == NORMAL1) {
    aflags = RTC_INTERSECT1;
  } else if (mode == NORMAL8) {
    aflags = RTC_INTERSECT8;
  } else if (mode == NORMAL16) {
    aflags = RTC_INTERSECT16;
  } else if (mode == STREAM_1M) {
    aflags = RTC_INTERSECT_STREAM;
  }

  // scene
  if (!scene_) {
    scene_ = rtcDeviceNewScene(device_, sflags, aflags);
  }

  // register domains
  unsigned geom_id = rtcNewUserGeometry(scene_, prims_.size());
  void* prims_ptr = (void*)&prims_[0];
  rtcSetUserData(scene_, geom_id, prims_ptr);

  // set callbacks for bounds
  rtcSetBoundsFunction(scene_, geom_id, cbBounds);

  // set callbacks for intersection tests
  if (mode == NORMAL1) {
    rtcSetIntersectFunction(scene_, geom_id, cbIntersect1);
  } else if (mode == NORMAL8) {
    rtcSetIntersectFunction8(scene_, geom_id, cbIntersect8);
  } else {
    LOG(FATAL) << "unsupported intersection mode " << mode;
  }
  // } else if (mode == NORMAL16) {
  //   rtcSetIntersectFunction16(scene_, geom_id, cbIntersect16);
  // }
  //   } else if (mode == STREAM_1M) {
  // #error unsupported
  //     // rtcSetIntersectFunctionN(scene_, geom_id, cbIntersectStream1M);
  //   }

  // commit scene
  rtcCommit(scene_);
}

void WbvhEmbree::cbBounds(void* ptr, std::size_t item, RTCBounds& bounds_o) {
  const WbvhEmbreePrimitive* prims = (const WbvhEmbreePrimitive*)ptr;
  const WbvhEmbreePrimitive& prim = prims[item];

  bounds_o.lower_x = prim.bounds_min[0];
  bounds_o.lower_y = prim.bounds_min[1];
  bounds_o.lower_z = prim.bounds_min[2];
  //
  bounds_o.upper_x = prim.bounds_max[0];
  bounds_o.upper_y = prim.bounds_max[1];
  bounds_o.upper_z = prim.bounds_max[2];
}

void WbvhEmbree::cbIntersect1(void* ptr, RTCRay& ray_i, std::size_t item) {
  //
  const WbvhEmbreePrimitive* prims = (const WbvhEmbreePrimitive*)ptr;
  RTCRayExt& ray = (RTCRayExt&)ray_i;

  const WbvhEmbreePrimitive& prim = prims[item];
  const Aabb aabb(prim.bounds_min, prim.bounds_max);

  float tmin, tmax;
  bool hit =
      intersectAabb(aabb, ray.org, ray.dir, ray.tnear, ray.tfar, &tmin, &tmax);

  if (hit) {
    DomainList* domains = ray.domains;
    unsigned offset = domains->count;

    CHECK_LT(offset, SPRAY_RAY_DOMAIN_LIST_SIZE);

    domains->count = offset + 1;
    domains->ids[offset] = prims[item].id;
    domains->ts[offset] = tmin;
  }
}

void WbvhEmbree::cbIntersect8(const void* valid, void* ptr, RTCRay8& ray,
                              std::size_t item) {
  // primitive
  const WbvhEmbreePrimitive* prims = (const WbvhEmbreePrimitive*)ptr;

  const WbvhEmbreePrimitive& prim = prims[item];
  const Aabb aabb(prim.bounds_min, prim.bounds_max);

  unsigned id = prim.id;

  // ray packet
  RTCRayExt8& packet = (RTCRayExt8&)ray;

  DomainList8* domains = packet.domains;

  // flag
  const unsigned* active = (const unsigned*)valid;

  // ray-box intersection tests
#pragma omp simd
  for (unsigned i = 0; i < 8; ++i) {
    if (active[i] != -1) continue;

    float org[3];
    org[0] = packet.orgx[i];
    org[1] = packet.orgy[i];
    org[2] = packet.orgz[i];
    float dir[3];
    dir[0] = packet.dirx[i];
    dir[1] = packet.diry[i];
    dir[2] = packet.dirz[i];

    float ray_tnear = packet.tnear[i];
    float ray_tfar = packet.tfar[i];

    float tmin, tmax;
    bool hit = intersectAabb(aabb, org, dir, ray_tnear, ray_tfar, &tmin, &tmax);

    unsigned count = domains->count[i];
    unsigned offset = count * 8 + i;

#ifdef SPRAY_GLOG_CHECK
    CHECK(offset < RAY8_DOMAIN_LIST_SIZE)
        << "domain list overflow. offset " << offset
        << " domains.count: " << count;
#endif

    domains->count[i] = count + (unsigned)hit;
    domains->ids[offset] = id;
    domains->ts[offset] = tmin;
  }
}

void WbvhEmbree::cbIntersect16(const void* valid, void* ptr, RTCRay16& ray,
                               std::size_t item) {
  // primitive
  const WbvhEmbreePrimitive* prims = (const WbvhEmbreePrimitive*)ptr;

  const WbvhEmbreePrimitive& prim = prims[item];
  const Aabb aabb(prim.bounds_min, prim.bounds_max);

  unsigned id = prim.id;

  // ray packet
  RTCRayExt16& packet = (RTCRayExt16&)ray;

  // flag
  const unsigned* active = (const unsigned*)valid;

  // ray-box intersection tests
#pragma omp simd
  for (unsigned i = 0; i < 16; ++i) {
    if (active[i] != -1) continue;

    const glm::vec3 ray_org =
        glm::vec3(packet.orgx[i], packet.orgy[i], packet.orgz[i]);

    const glm::vec3 ray_dir =
        glm::vec3(packet.dirx[i], packet.diry[i], packet.dirz[i]);

    float ray_tnear = packet.tnear[i];
    float ray_tfar = packet.tfar[i];

    float tmin, tmax;
    bool hit = intersectAabb(aabb, ray_org, ray_dir, ray_tnear, ray_tfar, &tmin,
                             &tmax);

    unsigned sp = packet.dom_count[i];
    unsigned offset = SPRAY_RAY_DOMAIN_LIST_SIZE * i + sp;

#ifdef SPRAY_GLOG_CHECK
    CHECK(offset < PACKET16_STACK_SIZE)
        << "domain stack overflow, packet.dom_count: " << sp;
#endif

    packet.dom_ids[offset] = id;
    packet.dom_ts[offset] = tmin;

    packet.dom_count[i] = sp + (unsigned)hit;
  }
}

}  // namespace spray
