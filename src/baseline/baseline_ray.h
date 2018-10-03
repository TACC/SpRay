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

#include <queue>

#include "glm/glm.hpp"
#include "glog/logging.h"

#include "renderers/rays.h"
#include "renderers/spray.h"

namespace spray {
namespace baseline {

struct SPRAY_ALIGN(16) DRay {
 public:
  float org[3];
  int pixid;

  float dir[3];
  int samid;

 public:
  float w[3];  //!< Radiance weight.
  int depth;   //!< Bounce number starting from 0.
  float t;
  float u;          //!< Barycentric u coordinate of hit
  float v;          //!< Barycentric v coordinate of hit
  unsigned geomID;  //!< geometry ID
  unsigned primID;
  unsigned flag;
  int domid;  // closest domain ID

  int domain_pos;   // current domain position
  float next_tdom;  // distance to next domain

#ifdef SPRAY_GLOG_CHECK
  friend std::ostream& operator<<(std::ostream& os, const DRay& r);
#endif
};

struct SPRAY_ALIGN(16) DRayQItem {
  int dummy;
  DRay* ray;
};

typedef std::queue<DRay*> DRayQ;

struct RayUtil {
  // inline static void makeEyeRay(const DRay& r, DomainList* domains,
  //                               RTCRayExt* rout) {
  //   rout->org[0] = r.org[0];
  //   rout->org[1] = r.org[1];
  //   rout->org[2] = r.org[2];

  //   rout->dir[0] = r.dir[0];
  //   rout->dir[1] = r.dir[1];
  //   rout->dir[2] = r.dir[2];

  //   rout->tnear = SPRAY_RAY_EPSILON;
  //   rout->tfar = SPRAY_FLOAT_INF;
  //   rout->geomID = RTC_INVALID_GEOMETRY_ID;
  //   rout->primID = RTC_INVALID_GEOMETRY_ID;

  //   rout->domains = domains;

  //   domains->count = 0;
  // }

  inline static void makeIntersection(const DRay& r,
                                      RTCRayIntersection* isect) {
    isect->org[0] = r.org[0];
    isect->org[1] = r.org[1];
    isect->org[2] = r.org[2];
    isect->dir[0] = r.dir[0];
    isect->dir[1] = r.dir[1];
    isect->dir[2] = r.dir[2];
    isect->tfar = r.t;
    isect->u = r.u;
    isect->v = r.v;
    isect->geomID = r.geomID;
    isect->primID = r.primID;
  }
};

//////////////////////////////////////////////////////////////////////////

enum DRayFlag {
  DRAY_FLAG_OCCLUDED = 0,
  DRAY_FLAG_SHADE = 1,
  DRAY_FLAG_SHADOW = 2
};

struct DRayUtil {
  /////////////////////////////// DRayUtil //////////////////////////////////
  inline static bool hasNextDomain(const DRay* ray) {
    return (std::isinf(ray->next_tdom) == false);
  }
  inline static bool hasCurrentDomain(const DRay* ray) {
    return (ray->domain_pos < INT_MAX);
  }

  /////////////////////////////// DRayUtil //////////////////////////////////

  // [0]: occluded (don't care for radiance)
  // [1]: shade
  // [2]: shadow
  //
  inline static uint32_t initShadowFlag() {
    return (0x00000001 << DRAY_FLAG_SHADOW);
  }

  inline static uint32_t initRadianceFlag() { return 0x00000000; }

  /////////////////////////////// DRayUtil //////////////////////////////////

  inline static void setOccluded(DRay* r) {
    r->flag |= (0x00000001 << DRAY_FLAG_OCCLUDED);
  }

  inline static uint32_t getOccluded(const DRay& r) {
    return ((r.flag >> DRAY_FLAG_OCCLUDED) & 0x00000001);
  }

  /////////////////////////////// DRayUtil //////////////////////////////////

  inline static void setShade(DRay* r) {
    r->flag |= (0x00000001 << DRAY_FLAG_SHADE);
  }
  inline static uint32_t getShade(const DRay& r) {
    return ((r.flag >> DRAY_FLAG_SHADE) & 0x00000001);
  }

  /////////////////////////////// DRayUtil //////////////////////////////////

  inline static void setShadow(DRay* r) {
    r->flag |= (0x00000001 << DRAY_FLAG_SHADOW);
  }
  inline static uint32_t getShadow(const DRay& r) {
    return ((r.flag >> DRAY_FLAG_SHADOW) & 0x00000001);
  }

  /////////////////////////////// DRayUtil //////////////////////////////////

  inline static void makeEyeRay(int pixid, int samid, DRay* r) {
    // assume org, dir already initialized
    r->pixid = pixid;
    r->samid = samid;
    r->depth = 0;

    r->w[0] = 1.0f;
    r->w[1] = 1.0f;
    r->w[2] = 1.0f;

    r->t = SPRAY_FLOAT_INF;
    r->geomID = RTC_INVALID_GEOMETRY_ID;
    r->primID = RTC_INVALID_GEOMETRY_ID;
    r->flag = DRayUtil::initRadianceFlag();

    r->domid = INT_MAX;
    r->domain_pos = INT_MAX;
    r->next_tdom = SPRAY_FLOAT_INF;
  }

  /////////////////////////////// DRayUtil //////////////////////////////////

  inline static void makeShadowRay(const DRay& rin, const glm::vec3& org,
                                   const glm::vec3& dir, const glm::vec3& w,
                                   float t_hit, DRay* r) {
    // assume org, dir already initialized
    r->org[0] = org[0];
    r->org[1] = org[1];
    r->org[2] = org[2];
    r->pixid = rin.pixid;

    r->dir[0] = dir[0];
    r->dir[1] = dir[1];
    r->dir[2] = dir[2];
    r->samid = rin.samid;

    r->depth = rin.depth + 1;

    r->w[0] = w[0];
    r->w[1] = w[1];
    r->w[2] = w[2];

    // r->t
    // r->u
    // r->v
    // r->geomID
    // r->primID

    r->flag = DRayUtil::initShadowFlag();

    r->domain_pos = INT_MAX;
    r->next_tdom = SPRAY_FLOAT_INF;
  }

  /////////////////////////////// DRayUtil //////////////////////////////////

  inline static void makeRadianceRay(int depth, const DRay& rin,
                                     const glm::vec3& org, const glm::vec3& dir,
                                     const glm::vec3& w, DRay* r) {
    // assume org, dir already initialized
    r->org[0] = org[0];
    r->org[1] = org[1];
    r->org[2] = org[2];
    r->pixid = rin.pixid;

    r->dir[0] = dir[0];
    r->dir[1] = dir[1];
    r->dir[2] = dir[2];
    r->samid = rin.samid;

    r->depth = depth;

    r->w[0] = w[0];
    r->w[1] = w[1];
    r->w[2] = w[2];

    r->t = SPRAY_FLOAT_INF;
    r->geomID = RTC_INVALID_GEOMETRY_ID;
    r->primID = RTC_INVALID_GEOMETRY_ID;
    r->flag = DRayUtil::initRadianceFlag();

    r->domid = INT_MAX;
    r->domain_pos = INT_MAX;
    r->next_tdom = SPRAY_FLOAT_INF;
  }

  /////////////////////////////// DRayUtil //////////////////////////////////

  inline static void updateIntersection(int id, const RTCRayIntersection& isect,
                                        DRay* r) {
    r->t = isect.tfar;
    r->u = isect.u;
    r->v = isect.v;
    r->geomID = isect.geomID;
    r->primID = isect.primID;

    r->domid = id;
  }
};  // end of DRayUtil

//////////////////////////////////////////////////////////////////////////

struct RTCRayUtil {
  inline static void makeEyeRay(const DRay& r, DomainList* domains,
                                RTCRayExt* rout) {
    rout->org[0] = r.org[0];
    rout->org[1] = r.org[1];
    rout->org[2] = r.org[2];

    rout->dir[0] = r.dir[0];
    rout->dir[1] = r.dir[1];
    rout->dir[2] = r.dir[2];

    // rout->tnear = 0.0f;
    rout->tnear = SPRAY_RAY_EPSILON;
    rout->tfar = SPRAY_FLOAT_INF;
    // rout->instID = RTC_INVALID_GEOMETRY_ID;
    rout->geomID = RTC_INVALID_GEOMETRY_ID;
    rout->primID = RTC_INVALID_GEOMETRY_ID;
    // rout->mask = 0xFFFFFFFF;
    // rout->time = 0.0f;

    rout->domains = domains;

    domains->count = 0;
  }

  /////////////////////////////// RTCRayUtil //////////////////////////////////

  inline static void makeRayForDomainIntersection(const float org[3],
                                                  const float dir[3],
                                                  DomainList* domains,
                                                  RTCRayExt* rout) {
    rout->org[0] = org[0];
    rout->org[1] = org[1];
    rout->org[2] = org[2];

    rout->dir[0] = dir[0];
    rout->dir[1] = dir[1];
    rout->dir[2] = dir[2];

    // rout->tnear = 0.0f;
    rout->tnear = SPRAY_RAY_EPSILON;
    rout->tfar = SPRAY_FLOAT_INF;
    // rout->instID = RTC_INVALID_GEOMETRY_ID;
    rout->geomID = RTC_INVALID_GEOMETRY_ID;
    rout->primID = RTC_INVALID_GEOMETRY_ID;
    // rout->mask = 0xFFFFFFFF;
    // rout->time = 0.0f;

    rout->domains = domains;

    domains->count = 0;
  }

  /////////////////////////////// RTCRayUtil //////////////////////////////////

  // sort domains front-to-back
  inline static void sortDomains(const DomainList& domains,
                                 DomainHit1 hits[SPRAY_RAY_DOMAIN_LIST_SIZE]) {
    int count = domains.count;

    for (int i = 0; i < count; ++i) {
      hits[i].id = domains.ids[i];
      hits[i].t = domains.ts[i];
    }

    std::sort(hits, hits + count,
              [&](const DomainHit1& a, const DomainHit1& b) {
                // return (a.t < b.t);
                return ((a.t < b.t) || ((a.t == b.t) && (a.id < b.id)));
              });
  }

  /////////////////////////////// RTCRayUtil //////////////////////////////////

  inline static void sortDomains8(int p, const DomainList8& domains,
                                  DomainHit1 hits[SPRAY_RAY_DOMAIN_LIST_SIZE]) {
    int count = domains.count[p];

    for (int i = 0; i < count; ++i) {
      int offset = i * 8 + p;
      hits[i].id = domains.ids[offset];
      hits[i].t = domains.ts[offset];
    }

    std::sort(hits, hits + count,
              [&](const DomainHit1& a, const DomainHit1& b) {
                // return (a.t < b.t);
                return ((a.t < b.t) || ((a.t == b.t) && (a.id < b.id)));
              });
  }

  /////////////////////////////// RTCRayUtil //////////////////////////////////

  // inline static void makeDomainIntersections8(
  //     const RTCRayExt8& packet, DomainHit1
  //     hits[PACKET8_STACK_SIZE]) {
  //   //
  //   for (int p = 0; p < 8; ++p) {
  //     int offset = p * SPRAY_RAY_DOMAIN_LIST_SIZE;

  //     for (int i = 0; i < packet[p].dom_count; ++i) {
  //       int idx = offset + i;
  //       hits[idx].id = packet[p].dom_ids[idx];
  //       hits[idx].t = packet[p].dom_ts[idx];
  //     }

  //     std::sort(&hits[offset], &hits[offset] + packet[p].dom_count,
  //               [&](const DomainHit1& a, const DomainHit1& b)
  //               {
  //                 return (a.t < b.t);
  //               });
  //   }
  // }

  /////////////////////////////// RTCRayUtil //////////////////////////////////

  // inline static void makeRadianceRay(int samid, int depth, const float
  // org[3],
  inline static void makeRadianceRay(const float org[3], const float dir[3],
                                     RTCRayIntersection* r) {
    r->org[0] = org[0];
    r->org[1] = org[1];
    r->org[2] = org[2];
    // r->samID = samid;
    r->dir[0] = dir[0];
    r->dir[1] = dir[1];
    r->dir[2] = dir[2];
    // r->depth = depth;
    // r->tnear = 0.f;
    r->tnear = SPRAY_RAY_EPSILON;
    r->tfar = SPRAY_FLOAT_INF;
    r->instID = RTC_INVALID_GEOMETRY_ID;
    r->geomID = RTC_INVALID_GEOMETRY_ID;
    r->primID = RTC_INVALID_GEOMETRY_ID;
    r->mask = 0xFFFFFFFF;
    r->time = 0.0f;
  }

  /////////////////////////////// RTCRayUtil //////////////////////////////////

  inline static void makeRadianceRay(const glm::vec3& org, const float dir[3],
                                     RTCRayIntersection* r) {
    r->org[0] = org[0];
    r->org[1] = org[1];
    r->org[2] = org[2];
    // r->samID = samid;
    r->dir[0] = dir[0];
    r->dir[1] = dir[1];
    r->dir[2] = dir[2];
    // r->depth = depth;
    // r->tnear = 0.f;
    r->tnear = SPRAY_RAY_EPSILON;
    r->tfar = SPRAY_FLOAT_INF;
    r->instID = RTC_INVALID_GEOMETRY_ID;
    r->geomID = RTC_INVALID_GEOMETRY_ID;
    r->primID = RTC_INVALID_GEOMETRY_ID;
    r->mask = 0xFFFFFFFF;
    r->time = 0.0f;
  }

  /////////////////////////////// RTCRayUtil //////////////////////////////////

  inline static void makeIntersection(const DRay& r,
                                      RTCRayIntersection* isect) {
    isect->org[0] = r.org[0];
    isect->org[1] = r.org[1];
    isect->org[2] = r.org[2];
    // isect->samID = r.samid;
    isect->dir[0] = r.dir[0];
    isect->dir[1] = r.dir[1];
    isect->dir[2] = r.dir[2];
    // isect->depth = r.depth;
    // // isect->tnear = 0.f;
    // isect->tnear = SPRAY_RAY_EPSILON;
    isect->tfar = r.t;
    // isect->instID = RTC_INVALID_GEOMETRY_ID;
    isect->u = r.u;
    isect->v = r.v;
    isect->geomID = r.geomID;
    isect->primID = r.primID;
    // isect->mask = 0xFFFFFFFF;
    // isect->time = 0.0f;
  }

  /////////////////////////////// RTCRayUtil //////////////////////////////////

  inline static void makeShadowRay(const glm::vec3& org, const glm::vec3& dir,
                                   RTCRay* r) {
    r->org[0] = org[0];
    r->org[1] = org[1];
    r->org[2] = org[2];
    r->dir[0] = dir[0];
    r->dir[1] = dir[1];
    r->dir[2] = dir[2];
    // r->tnear = 0.001f;
    r->tnear = SPRAY_RAY_EPSILON;
    r->tfar = SPRAY_FLOAT_INF;
    r->geomID = RTC_INVALID_GEOMETRY_ID;
    r->primID = RTC_INVALID_GEOMETRY_ID;
    r->mask = -1;
    r->time = 0;
  }

  /////////////////////////////// RTCRayUtil //////////////////////////////////

  inline static void makeShadowRay(const float org[3], const float dir[3],
                                   RTCRay* r) {
    r->org[0] = org[0];
    r->org[1] = org[1];
    r->org[2] = org[2];
    r->dir[0] = dir[0];
    r->dir[1] = dir[1];
    r->dir[2] = dir[2];
    // r->tnear = 0.001f;
    r->tnear = SPRAY_RAY_EPSILON;
    r->tfar = SPRAY_FLOAT_INF;
    r->geomID = RTC_INVALID_GEOMETRY_ID;
    r->primID = RTC_INVALID_GEOMETRY_ID;
    r->mask = -1;
    r->time = 0;
  }

  /////////////////////////////// RTCRayUtil //////////////////////////////////

  inline static void hitPosition(const float org[3], const float dir[3],
                                 const float t, float pos[3]) {
    pos[0] = dir[0] * t + org[0];
    pos[1] = dir[1] * t + org[1];
    pos[2] = dir[2] * t + org[2];
  }

  /////////////////////////////// RTCRayUtil //////////////////////////////////

  inline static glm::vec3 hitPosition(const float org[3], const float dir[3],
                                      float t) {
    glm::vec3 pos(dir[0] * t + org[0], dir[1] * t + org[1],
                  dir[2] * t + org[2]);
    return pos;
  }

  /////////////////////////////// RTCRayUtil //////////////////////////////////

  inline static glm::vec3 hitPosition(const glm::vec3& org, const float dir[3],
                                      float t) {
    glm::vec3 pos(dir[0] * t + org[0], dir[1] * t + org[1],
                  dir[2] * t + org[2]);
    return pos;
  }

  /////////////////////////////// RTCRayUtil //////////////////////////////////

  inline static glm::vec3 hitPosition(const RTCRayIntersection& isect) {
    float t = isect.tfar;
    glm::vec3 pos(isect.dir[0] * t + isect.org[0],
                  isect.dir[1] * t + isect.org[1],
                  isect.dir[2] * t + isect.org[2]);
    return pos;
  }
};  // end of struct RTCRayUtil

}  // namespace baseline
}  // namespace spray
