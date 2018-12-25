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

#include "glm/glm.hpp"
#include "glog/logging.h"
#include "pbrt/memory3.h"

#include "cmake_config.h"  // auto generated by cmake

#include "baseline/baseline_ray.h"
#include "partition/arena_queue.h"
#include "partition/data_partition.h"
#include "render/spray.h"
#include "scene/scene.h"

namespace spray {
namespace baseline {

template <typename CacheT>
class DomainIntersector {
 public:
  DomainIntersector(int ndomains, spray::Scene<CacheT>* scene)
      : ndomains_(ndomains), scene_(scene) {}

  void intersect(int current_id, DRay* ray, ArenaQs<DRayQItem>* qs) {
    // setup a ray
    RTCRayUtil::makeRayForDomainIntersection(ray->org, ray->dir, &domains_,
                                             &eray_);
    // ray-domain intersection tests
    scene_->intersectDomains(eray_);

#ifdef SPRAY_GLOG_CHECK
    CHECK_LE(domains_.count, SPRAY_RAY_DOMAIN_LIST_SIZE);
#endif
    if (domains_.count) {
      RTCRayUtil::sortDomains(domains_, hits_);
    }

    // assume that ray generation and processing are done separately
    // therefore, next domain_pos becomes 0 when the ray is newly generated
    // ray->domain_pos = (ray->domain_pos == INT_MAX) ? 0 : ray->domain_pos + 1;

    // newly generated ray, exclude current domain
    int dom_pos = ray->domain_pos;

    if (dom_pos < INT_MAX) {  // ray traversing domains
      ray->domain_pos = dom_pos + 1;
    } else {  // new ray
      if (domains_.count && hits_[0].id == current_id) {
        ray->domain_pos = 1;
      } else {
        ray->domain_pos = 0;
      }
    }

    if (ray->domain_pos < domains_.count) {
      // two domains away from current domain
      int next_pos = ray->domain_pos + 1;
      ray->next_tdom =
          (next_pos < domains_.count) ? hits_[next_pos].t : SPRAY_FLOAT_INF;
      if (qs) {
        ray_data_.ray = ray;
        qs->copy(hits_[ray->domain_pos].id, &ray_data_);
      }
#ifdef SPRAY_GLOG_CHECK
      CHECK_NE(hits_[ray->domain_pos].id, current_id);
#endif
    } else {
      ray->domain_pos = INT_MAX;
      ray->next_tdom = SPRAY_FLOAT_INF;
    }
  }

 private:
  DomainList domains_;
  DomainHit1 hits_[SPRAY_RAY_DOMAIN_LIST_SIZE];

  RTCRayExt eray_;
  DRayQItem ray_data_;

  int ndomains_;
  spray::Scene<CacheT>* scene_;
};

}  // namespace baseline
}  // namespace spray
