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

#include "ooc/ooc_domain_stats.h"
#include "ooc/ooc_ray.h"
#include "render/qvector.h"
#include "render/scene.h"

namespace spray {
namespace ooc {

template <typename SceneT>
class Isector {
 public:
  void init(std::size_t num_domains) { domains_.resize(num_domains); }

 public:
  void intersect(const SceneT* scene, Ray* ray, spray::QVector<RayData>* qs,
                 DomainStats* stats) {
    isectAll(scene, ray, qs, stats);
  }
  //     RTCRayUtil::makeRayForDomainIntersection(ray->org, ray->dir, &domains_,
  //                                              &eray_);
  //
  //     scene->intersectDomains(eray_);
  //
  //     int rank = mpi::rank();
  //
  //     if (domains_.count) {
  // #ifdef SPRAY_GLOG_CHECK
  //       CHECK_LT(domains_.count, SPRAY_RAY_DOMAIN_LIST_SIZE);
  // #endif
  //       RTCRayUtil::sortDomains(domains_, hits_);
  //
  //       for (int d = 0; d < domains_.count; ++d) {
  //         int id = hits_[d].id;
  // #ifdef SPRAY_GLOG_CHECK
  //         CHECK_LT(id, ndomains);
  // #endif
  //         ray_data_.ray = ray;
  //         ray_data_.tdom = hits_[d].t;
  //         ray_data_.dom_depth = d;
  //
  //         qs->push(id, ray_data_);
  //
  //         stats->increment(id, d /*depth*/);
  //       }
  //     }
  //   }

  bool intersect(int exclude_id, const SceneT* scene, Ray* ray,
                 spray::QVector<RayData>* qs, DomainStats* stats) {
    return isectWithoutCurrentDomain(exclude_id, scene, ray, qs, stats);
  }

  //     RTCRayUtil::makeRayForDomainIntersection(ray->org, ray->dir, &domains_,
  //                                              &eray_);
  //
  //     scene->intersectDomains(eray_);
  //
  //     int rank = mpi::rank();
  //
  //     int num_hit_domains = 0;
  //
  //     if (domains_.count) {
  // #ifdef SPRAY_GLOG_CHECK
  //       CHECK_LT(domains_.count, SPRAY_RAY_DOMAIN_LIST_SIZE);
  // #endif
  //       RTCRayUtil::sortDomains(domains_, hits_);
  //
  //       for (int d = 0; d < domains_.count; ++d) {
  //         int id = hits_[d].id;
  // #ifdef SPRAY_GLOG_CHECK
  //         CHECK_LT(id, ndomains);
  // #endif
  //         if (id != exclude_id) {
  //           ++num_hit_domains;
  //
  //           ray_data_.ray = ray;
  //           ray_data_.tdom = hits_[d].t;
  //           ray_data_.dom_depth = d;
  //
  //           qs->push(id, ray_data_);
  //
  //           stats->increment(id, d /*depth*/);
  //         }
  //       }
  //     }
  //     return (num_hit_domains > 0);
  //   }

 private:
  void isectDomains(const SceneT* scene, const Ray* ray) {
#ifdef SPRAY_GLOG_CHECK
    CHECK_GT(domains_.size(), 0);
#endif
    domains_.reset();
    eray_.reset(ray->org, ray->dir, &domains_);
    scene->intersectDomains(eray_);
    domains_.sort();
  }

  void isectAll(const SceneT* scene, Ray* ray, spray::QVector<RayData>* qs,
                DomainStats* stats) {
#ifdef SPRAY_GLOG_CHECK
    CHECK_GT(domains_.size(), 0);
#endif
    isectDomains(scene, ray);

    for (std::size_t i = 0; i < domains_.getNumHits(); ++i) {
#ifdef SPRAY_GLOG_CHECK
      CHECK_LT(domains_.getId(i), domains_.size());
#endif
      ray_data_.ray = ray;
      ray_data_.tdom = domains_.getTnear(i);
      ray_data_.dom_depth = i;

      int id = domains_.getId(i);
      qs->push(id, ray_data_);
      stats->increment(id, i /*depth*/);
    }
  }

  bool isectWithoutCurrentDomain(int exclude_id, const SceneT* scene, Ray* ray,
                                 spray::QVector<RayData>* qs,
                                 DomainStats* stats) {
#ifdef SPRAY_GLOG_CHECK
    CHECK_GT(domains_.size(), 0);
#endif
    isectDomains(scene, ray);

    int num_hit_domains = 0;

    for (std::size_t i = 0; i < domains_.getNumHits(); ++i) {
#ifdef SPRAY_GLOG_CHECK
      CHECK_LT(domains_.getId(i), domains_.size());
#endif
      int id = domains_.getId(i);

      if (id != exclude_id) {
        ++num_hit_domains;
        ray_data_.ray = ray;
        ray_data_.tdom = domains_.getTnear(i);
        ray_data_.dom_depth = i;

        qs->push(id, ray_data_);
        stats->increment(id, i /*depth*/);
      }
    }
    return (num_hit_domains > 0);
  }

 private:
  DomainList domains_;
  RTCRayExt eray_;
  RayData ray_data_;

  // private:
  //  DomainList domains_;
  //  DomainHit1 hits_[SPRAY_RAY_DOMAIN_LIST_SIZE];

  //  RTCRayExt eray_;
  //  RayData ray_data_;
};

}  // namespace ooc
}  // namespace spray
