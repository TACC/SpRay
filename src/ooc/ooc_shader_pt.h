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
#include "pbrt/memory.h"

#include "ooc/ooc_ray.h"
#include "render/config.h"
#include "render/light.h"
#include "render/material.h"
#include "render/rays.h"
#include "render/reflection.h"
#include "utils/util.h"

namespace spray {
namespace ooc {

template<typename SceneT>
class ShaderPt {
 public:
  void init(const spray::Config &cfg, const SceneT &scene) {
    bounces_ = cfg.bounces;
    samples_ = cfg.ao_samples;  // number of samples for area lights
    ks_ = cfg.ks;
    shininess_ = cfg.shininess;
    lights_ = scene.getLights();  // copy lights
  }

 private:
  std::vector<Light *> lights_;
  int bounces_;
  int samples_;
  glm::vec3 ks_;
  float shininess_;

 public:
  bool isAo() { return false; }

 public:
  void operator()(int domain_id, const Ray &rayin,
                  const spray::RTCRayIntersection &isect,
                  spray::MemoryArena *mem, std::queue<Ray *> *sq,
                  std::queue<Ray *> *rq, std::queue<Ray *> *pending_q,
                  int ray_depth);

 private:
  void genR2(const Ray &rayin, const glm::vec3 &org, const glm::vec3 &dir,
             const glm::vec3 &w, float t, spray::MemoryArena *mem,
             std::queue<Ray *> *rq, std::queue<Ray *> *pending_q) {
    Ray *r2 = mem->Alloc<Ray>(1, false);
    CHECK_NOTNULL(r2);

    int next_virtual_depth = rayin.depth + 1;
    if (next_virtual_depth == SPRAY_HISTORY_SIZE) {
      RayUtil::makeRay(rayin, org, dir, w, t, r2, 0);
      pending_q->push(r2);
#ifdef SPRAY_GLOG_CHECK
      CHECK_EQ(r2->depth, 0);
#endif
    } else {
      RayUtil::makeRay(rayin, org, dir, w, t, r2, next_virtual_depth);
      rq->push(r2);
    }
  }
};

template <typename SceneT>
void ShaderPt<SceneT>::operator()(int domain_id, const Ray &rayin,
                                  const spray::RTCRayIntersection &isect,
                                  spray::MemoryArena *mem,
                                  std::queue<Ray *> *sq, std::queue<Ray *> *rq,
                                  std::queue<Ray *> *pending_q, int ray_depth) {
  glm::vec3 pos = RTCRayUtil::hitPosition(rayin.org, rayin.dir, isect.tfar);

  glm::vec3 surf_radiance;

  auto *material = isect.material;

  bool is_shape = (isect.color == SPRAY_INVALID_COLOR);

  glm::vec3 albedo;
  if (is_shape) {
    albedo = material->getAlbedo();
  } else {
    util::unpack(isect.color, albedo);
  }

  glm::vec3 normal(isect.Ns[0], isect.Ns[1], isect.Ns[2]);

  glm::vec3 wo(-rayin.dir[0], -rayin.dir[1], -rayin.dir[2]);
  wo = glm::normalize(wo);

  glm::vec3 Lin(rayin.w[0], rayin.w[1], rayin.w[2]);

#ifdef SPRAY_FACE_FORWARD_OFF
  glm::vec3 normal_ff = glm::normalize(normal);
#else
  glm::vec3 normal_ff;
  if (is_shape) {
    normal_ff = glm::normalize(normal);
  } else {
    float cos_theta_i = glm::dot(wo, normal);
    bool entering = (cos_theta_i > 0.0f);
    normal_ff = glm::normalize(entering ? normal : -normal);
  }
#endif

  glm::vec3 wi, light_color, Lr;
  float pdf, inv_shade_pdf, costheta;
  int nlights = lights_.size();

  std::size_t color_idx = sq->size();

  int next_virtual_depth = rayin.depth + 1;
  int next_ray_depth = ray_depth + next_virtual_depth;

  RandomSampler sampler;
  RandomSampler_init(sampler, rayin.samid * next_ray_depth);

  int light_sample_offset = 0;
  int num_light_samples;
  int light_sample_id;

  glm::vec3 shade_color;

  // direct illumination

  if (material->hasDiffuse()) {
    for (int l = 0; l < nlights; ++l) {
      Light *light = lights_[l];

      num_light_samples = light->getNumSamples();

      // for each light sample
      for (int s = 0; s < num_light_samples; ++s) {
        // light color
        light_color = light->sampleL(pos, sampler, normal_ff, &wi, &pdf);

        if (pdf > 0.0f) {
          // wi, wo, normal_ff: all normalized
          shade_color = material->shade(albedo, wi, wo, normal_ff);

          Lr = Lin * light_color * shade_color *
               (1.0f / (pdf * num_light_samples));

          if (hasPositive(Lr)) {
            // create shadow ray
            Ray *shadow = mem->Alloc<Ray>(1, false);
            CHECK_NOTNULL(shadow);

            light_sample_id = light_sample_offset + s;

            RayUtil::makeShadow(rayin, light_sample_id, pos, wi, Lr, isect.tfar,
                                shadow);

            sq->push(shadow);
          }
        }
      }

      light_sample_offset += num_light_samples;
    }
  }

  // indirect illumination

#ifdef SPRAY_GLOG_CHECK
  CHECK_LT(ray_depth + rayin.depth, bounces_);
#endif

  if (next_ray_depth < bounces_) {
    RandomSampler_init(sampler, rayin.samid * next_ray_depth);
    glm::vec3 weight;
    bool valid =
        material->sample(albedo, wo, normal_ff, sampler, &wi, &weight, &pdf);
    if (valid) {
      Lr = Lin * weight * (1.0f / pdf);
      if (hasPositive(Lr)) {
        genR2(rayin, pos, wi, Lr, isect.tfar, mem, rq, pending_q);
      }
    }
  }
}

/*
template <typename CacheT, typename SceneT>
void ShaderPt<CacheT, SceneT>::operator()(
    int domain_id, const Ray &rayin, const spray::RTCRayIntersection &isect,
    spray::MemoryArena *mem, std::queue<Ray *> *sq, std::queue<Ray *> *rq,
    std::queue<Ray *> *pending_q, int ray_depth) {
  glm::vec3 pos = RTCRayUtil::hitPosition(rayin.org, rayin.dir, isect.tfar);
  glm::vec3 surf_radiance;
  util::unpack(isect.color, surf_radiance);

  glm::vec3 normal(isect.Ns[0], isect.Ns[1], isect.Ns[2]);
  glm::vec3 wo(-rayin.dir[0], -rayin.dir[1], -rayin.dir[2]);
  glm::vec3 Lin(rayin.w[0], rayin.w[1], rayin.w[2]);

  float cos_theta_i = glm::dot(wo, normal);
  bool entering = (cos_theta_i > 0.0f);
  glm::vec3 normal_ff = entering ? normal : -normal;
  normal_ff = glm::normalize(normal_ff);

  glm::vec3 wi, light_radiance, Lr;
  float pdf, costheta;
  int nlights = lights_.size();

  std::size_t color_idx = sq->size();

  Bsdf *bsdf = scene_->getBsdf(domain_id);
  bool delta_dist = bsdf->isDelta();

  int next_virtual_depth = rayin.depth + 1;
  int next_actual_depth = ray_depth + next_virtual_depth;

#ifdef SPRAY_GLOG_CHECK
  CHECK_LE(next_virtual_depth, SPRAY_HISTORY_SIZE);
#endif

  if (!delta_dist) {
    RandomSampler light_sampler;
    RandomSampler_init(light_sampler, rayin.samid * next_actual_depth);

    for (int l = 0; l < nlights; ++l) {
      if (lights_[l]->isAreaLight()) {
        for (int s = 0; s < samples_; ++s) {
          light_radiance =
              lights_[l]->sampleArea(light_sampler, normal_ff, &wi, &pdf);

          if (pdf > 0.0f) {
            costheta = glm::clamp(glm::dot(normal_ff, wi), 0.0f, 1.0f);
            Lr = Lin *
                 blinnPhong(costheta, surf_radiance, ks_, shininess_,
                            light_radiance, wi, normal_ff, wo) *
                 (1.0f / (pdf * samples_));

            if (hasPositive(Lr)) {
              Ray *shadow = mem->Alloc<Ray>(1, false);
              CHECK_NOTNULL(shadow);

              RayUtil::makeShadow(rayin, l, pos, wi, Lr, isect.tfar, shadow);

              sq->push(shadow);
            }
          }
        }
      } else {  // point light
        light_radiance = lights_[l]->sample(pos, &wi, &pdf);
        if (pdf > 0.0f) {
          // evaluate direct lighting
          costheta = glm::clamp(glm::dot(normal_ff, wi), 0.0f, 1.0f);
          Lr = Lin *
               blinnPhong(costheta, surf_radiance, ks_, shininess_,
                          light_radiance, wi, normal_ff, wo) *
               (1.0f / pdf);

          if (hasPositive(Lr)) {
            Ray *shadow = mem->Alloc<Ray>(1, false);
            CHECK_NOTNULL(shadow);

            RayUtil::makeShadow(rayin, l, pos, wi, Lr, isect.tfar, shadow);
            sq->push(shadow);
          }
        }
      }
    }
  }  // if (!delta_dist) {

  if (next_actual_depth < bounces_) {
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
        bool has_reflect = hasReflection(sample_type);

        if (has_reflect) {
          glm::vec3 wr = Reflect(wo, normal_ff);
          wi = glm::normalize(wr);
          Lr = Lin * (fr / abs_cos_theta_i);
          if (hasPositive(Lr)) {
            genR2(rayin, pos, wi, Lr, isect.tfar, mem, rq, pending_q);
          }
        }

        if (hasTransmission(sample_type)) {
          // TODO: support both reflection and refraction
          CHECK_EQ(has_reflect, false);
          wi = glm::normalize(wt);
          Lr = Lin * ((1.0f - fr) / abs_cos_theta_i);
          if (hasPositive(Lr)) {
            genR2(rayin, pos, wi, Lr, isect.tfar, mem, rq, pending_q);
          }
        }
      }
    } else {
      RandomSampler sampler;
      RandomSampler_init(sampler, rayin.samid * next_actual_depth);

      bsdf->sampleRandom(normal_ff, &sampler, &wi, &pdf);
      costheta = glm::clamp(glm::dot(normal_ff, wi), 0.0f, 1.0f);
      Lr = Lin * surf_radiance * SPRAY_ONE_OVER_PI * costheta / pdf;
      if (hasPositive(Lr)) {
        genR2(rayin, pos, wi, Lr, isect.tfar, mem, rq, pending_q);
      }
    }
  }
}
*/

}  // namespace ooc
}  // namespace spray
