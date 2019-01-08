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

#include "embree/random_sampler.h"
#include "glm/glm.hpp"
#include "glog/logging.h"

#include "baseline/baseline_isector.h"
#include "baseline/baseline_ray.h"
#include "render/arena_queue.h"
#include "render/config.h"
#include "render/light.h"
#include "render/reflection.h"
#include "utils/util.h"

namespace spray {
namespace baseline {

template <typename CacheT, typename SceneT>
class ShaderAo {
 public:
  void init(const spray::Config &cfg, SceneT *scene) {
    bounces_ = cfg.bounces;
    samples_ = cfg.ao_samples;
    scene_ = scene;
    lights_ = scene->getLights();  // copy lights
  }

  bool isAo() { return true; }

  void operator()(int id, int ndomains, const DRay &rayin, RTCRay *rtc_ray,
                  RTCRayIntersection *isect, spray::ArenaQs<DRayQItem> *qs,
                  spray::ArenaQs<DRayQItem> *sqs, spray::MemoryArena *arena,
                  DomainIntersector<CacheT, SceneT> *domain_isector,
                  CommitBufferB *retire_buf, DRayQ *temp_q);

 private:
  SceneT *scene_;
  std::vector<Light *> lights_;  // copied lights
  int bounces_;
  int samples_;
};

template <typename CacheT, typename SceneT>
void ShaderAo<CacheT, SceneT>::operator()(
    int id, int ndomains, const DRay &rayin, RTCRay *rtc_ray,
    RTCRayIntersection *isect, spray::ArenaQs<DRayQItem> *qs,
    spray::ArenaQs<DRayQItem> *sqs, spray::MemoryArena *arena,
    DomainIntersector<CacheT, SceneT> *domain_isector,
    CommitBufferB *retire_buf, DRayQ *temp_q) {
  // TODO
}

/*
template <typename CacheT, typename SceneT>
void ShaderAo<CacheT, SceneT>::operator()(int id, int ndomains, const
DRay &rayin, RTCRay *rtc_ray, RTCRayIntersection *isect,
spray::ArenaQs<DRayQItem> *qs, spray::ArenaQs<DRayQItem> *sqs,
spray::MemoryArena *arena, DomainIntersector<CacheT, SceneT> *domain_isector,
                  CommitBufferB *retire_buf, DRayQ *temp_q) {
    // temporary variables
    glm::vec3 pos, wi, light_radiance, Lr;

    // weight
    glm::vec3 Lin(rayin.w[0], rayin.w[1], rayin.w[2]);

    // surface radiance
    glm::vec3 surf_radiance;
    util::unpack(isect->color, surf_radiance);

    // direct illumination

    float t_hit = isect->tfar;
    pos = RTCRayUtil::hitPosition(isect->org, isect->dir, t_hit);

    float pdf, costheta;

    glm::vec3 normal_ff = faceForwardFloat(isect->dir, isect->Ns);
    normal_ff = glm::normalize(normal_ff);

    const float ao_weight = 1.0f / static_cast<float>(samples_);
    spray::RandomSampler sampler;
    Bsdf *bsdf = scene_->getBsdf(id);

    for (int i = 0; i < samples_; ++i) {
      RandomSampler_init(sampler, rayin.pixid * (i + 1));
      bsdf->sampleRandom(normal_ff, &sampler, &wi, &pdf);

      costheta = glm::clamp(glm::dot(normal_ff, wi), 0.f, 1.f);
      Lr = surf_radiance * (SPRAY_ONE_OVER_PI * costheta * ao_weight / pdf);

      if (hasPositive(Lr)) {
        if (!scene_->occluded(pos, wi, rtc_ray)) {  // occluded
          DRay *dray = arena->Alloc<DRay>(1, false);
          DRayUtil::makeShadowRay(rayin, pos, wi, Lr, t_hit, dray);

          domain_isector->intersect(id, dray, sqs);
#ifdef SPRAY_GLOG_CHECK
          CHECK(sqs->empty(id));
#endif
          if (!DRayUtil::hasCurrentDomain(dray)) {
            retire_buf->commit(dray->pixid, Lr);
          }
        }
      }
    }
  }
*/

}  // namespace baseline
}  // namespace spray
