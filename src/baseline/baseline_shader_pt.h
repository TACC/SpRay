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

#include "baseline/baseline_commit_buffer.h"
#include "baseline/baseline_isector.h"
#include "baseline/baseline_ray.h"
#include "render/arena_queue.h"
#include "render/light.h"
#include "render/reflection.h"
#include "utils/math.h"
#include "utils/util.h"

namespace spray {
namespace baseline {

template <typename CacheT, typename SceneT>
class ShaderPt {
 public:
  void init(const spray::Config &cfg, SceneT *scene) {
    bounces_ = cfg.bounces;
    samples_ = cfg.ao_samples;
    scene_ = scene;
    lights_ = scene->getLights();  // copy lights
    ks_ = cfg.ks;
    shininess_ = cfg.shininess;
  }

  bool isAo() { return false; }

  void operator()(int id, int ndomains, const DRay &rayin, RTCRay *rtc_ray,
                  RTCRayIntersection *isect, spray::ArenaQs<DRayQItem> *qs,
                  spray::ArenaQs<DRayQItem> *sqs, spray::MemoryArena *arena,
                  DomainIntersector<CacheT, SceneT> *domain_isector,
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

    glm::vec3 normal(isect->Ns[0], isect->Ns[1], isect->Ns[2]);
    glm::vec3 wo(-isect->dir[0], -isect->dir[1], -isect->dir[2]);
    float cos_theta_i = glm::dot(wo, normal);
    bool entering = (cos_theta_i > 0.0f);
    glm::vec3 normal_ff = entering ? normal : -normal;
    normal_ff = glm::normalize(normal_ff);

    // retire_buf->commit(rayin.pixid, surf_radiance);

    Bsdf *bsdf = scene_->getBsdf(id);
    bool delta_dist = bsdf->isDelta();

    int next_ray_depth = rayin.depth + 1;
    spray::RandomSampler light_sampler;
    spray::RandomSampler_init(light_sampler, rayin.samid * next_ray_depth);

    if (!delta_dist) {
      for (std::size_t l = 0; l < lights_.size(); ++l) {
        if (lights_[l]->isAreaLight()) {
          for (int s = 0; s < samples_; ++s) {
            // sample() returns normalized direction
            light_radiance =
                lights_[l]->sampleArea(light_sampler, normal_ff, &wi, &pdf);
            if (pdf > 0.0f) {
              // evaluate direct lighting
              costheta = glm::clamp(glm::dot(normal_ff, wi), 0.0f, 1.0f);
              Lr = Lin *
                   blinnPhong(costheta, surf_radiance, ks_, shininess_,
                              light_radiance, wi, normal_ff, wo) *
                   (1.0f / (pdf * samples_));

              if (hasPositive(Lr)) {
                // intersect domains
                if (!scene_->occluded(pos, wi, rtc_ray)) {  // occluded
                  DRay *dray = arena->Alloc<DRay>(1, false);
                  DRayUtil::makeShadowRay(rayin, pos, wi, Lr, t_hit, dray);

                  domain_isector->intersect(id, dray, sqs);
#ifdef SPRAY_GLOG_CHECK
                  CHECK(sqs->empty(id));
#endif
                  // no more domain
                  if (!DRayUtil::hasCurrentDomain(dray)) {
                    retire_buf->commit(dray->pixid, Lr);
                  }
                }
              }
            }
          }
        } else {  // point light
          // sample() returns normalized direction
          light_radiance = lights_[l]->sample(pos, &wi, &pdf);
          if (pdf > 0.0f) {
            // evaluate direct lighting
            costheta = glm::clamp(glm::dot(normal_ff, wi), 0.0f, 1.0f);
            Lr = Lin *
                 blinnPhong(costheta, surf_radiance, ks_, shininess_,
                            light_radiance, wi, normal_ff, wo) *
                 (1.0f / pdf);

            if (hasPositive(Lr)) {
              // intersect domains
              if (!scene_->occluded(pos, wi, rtc_ray)) {  // occluded
                DRay *dray = arena->Alloc<DRay>(1, false);
                DRayUtil::makeShadowRay(rayin, pos, wi, Lr, t_hit, dray);
                // LOG(FATAL) << "[shadow ray]" << *dray;

                domain_isector->intersect(id, dray, sqs);
#ifdef SPRAY_GLOG_CHECK
                CHECK(sqs->empty(id));
#endif
                // no more domain
                if (!DRayUtil::hasCurrentDomain(dray)) {
                  retire_buf->commit(dray->pixid, Lr);
                }
              }
            }
          }
        }
      }
    }

    if (next_ray_depth < bounces_) {
      wo = glm::normalize(wo);

      if (delta_dist) {
        if (cos_theta_i != 0.0f) {  // rule out 90 degree
          cos_theta_i = glm::clamp(cos_theta_i, -1.0f, 1.0f);
          float abs_cos_theta_i = glm::abs(cos_theta_i);

          if (!entering) {  // i.e cos_theta_i < 0.0f
            cos_theta_i = abs_cos_theta_i;
          }

          uint32_t sample_type;
          float fr;      // prob. of reflection
          glm::vec3 wt;  // direction of transmitted ray
          bsdf->sampleDelta(entering, cos_theta_i, wo, normal_ff, &sample_type,
                            &fr, &wt);

          if (hasReflection(sample_type)) {
            glm::vec3 wr = Reflect(wo, normal_ff);
            wi = glm::normalize(wr);
            // Lr = Lin * fr;
            Lr = Lin * (fr / abs_cos_theta_i);
            if (hasPositive(Lr))
              genRadianceRay(next_ray_depth, rayin, pos, wi, Lr, arena, temp_q);
          }

          if (hasTransmission(sample_type)) {
            wi = glm::normalize(wt);
            Lr = Lin * ((1.0f - fr) / abs_cos_theta_i);
            if (hasPositive(Lr))
              genRadianceRay(next_ray_depth, rayin, pos, wi, Lr, arena, temp_q);
          }
        }
      } else {
        spray::RandomSampler sampler;
        spray::RandomSampler_init(sampler, rayin.samid * next_ray_depth);
        //
        bsdf->sampleRandom(normal_ff, &sampler, &wi, &pdf);
        costheta = glm::clamp(glm::dot(normal_ff, wi), 0.0f, 1.0f);
        Lr = Lin * surf_radiance * SPRAY_ONE_OVER_PI * costheta / pdf;
        if (hasPositive(Lr))
          genRadianceRay(next_ray_depth, rayin, pos, wi, Lr, arena, temp_q);
      }
    }
  }

 private:
  void genRadianceRay(int depth, const DRay &rin, const glm::vec3 &org,
                      const glm::vec3 &dir, const glm::vec3 &w,
                      spray::MemoryArena *arena, DRayQ *temp_q) {
    DRay *r = arena->Alloc<DRay>(1, false);
    DRayUtil::makeRadianceRay(depth, rin, org, dir, w, r);
    temp_q->push(r);
  }

 private:
  SceneT *scene_;
  std::vector<Light *> lights_;  // copied lights
  int bounces_;
  int samples_;
  glm::vec3 ks_;
  float shininess_;
};

}  // namespace baseline
}  // namespace spray
