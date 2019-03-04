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
#include "utils/util.h"

namespace spray {
namespace ooc {

template <typename SceneT>
class ShaderAo {
 public:
  typedef SceneT SceneType;

  ShaderAo() { matte_material_ = new Matte; }
  ~ShaderAo() { delete matte_material_; }

  void init(const spray::Config &cfg, const SceneT &scene) {
    bounces_ = cfg.bounces;
    samples_ = cfg.ao_samples;  // number of samples for area lights
    ao_weight_ = 1.0f / static_cast<float>(samples_);
#ifdef SPRAY_GLOG_CHECK
    num_pixels_ = cfg.image_w * cfg.image_h;
#endif
    use_spray_color_ = cfg.use_spray_color;
  }

  bool isAo() { return true; }

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

 private:
  const Material *matte_material_;
  int bounces_;
  int samples_;
  float ao_weight_;
#ifdef SPRAY_GLOG_CHECK
  int num_pixels_;
#endif
  bool use_spray_color_;
};

template <typename SceneT>
void ShaderAo<SceneT>::operator()(int domain_id, const Ray &rayin,
                                  const spray::RTCRayIntersection &isect,
                                  spray::MemoryArena *mem,
                                  std::queue<Ray *> *sq, std::queue<Ray *> *rq,
                                  std::queue<Ray *> *pending_q, int ray_depth) {
  glm::vec3 pos = RTCRayUtil::hitPosition(rayin.org, rayin.dir, isect.tfar);

  bool is_shape = (isect.color == SPRAY_INVALID_COLOR);

  glm::vec3 albedo;
#define SPRAY_AO_SHOW_SURFACE_COLOR
#ifdef SPRAY_AO_SHOW_SURFACE_COLOR
  auto *material = isect.material;
  if (material->type() != Material::MATTE) {
    material = matte_material_;
  }
  if (is_shape || use_spray_color_) {
    albedo = material->getAlbedo();
  } else {
    util::unpack(isect.color, albedo);
  }
#else
  auto *material = matte_material_;
  albedo = material->getAlbedo();
#endif

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

  glm::vec3 wi, Lr, weight;
  float pdf;

  int next_virtual_depth = rayin.depth + 1;
  int next_ray_depth = ray_depth + next_virtual_depth;

  RandomSampler sampler;
  RandomSampler_init(sampler, rayin.samid * next_ray_depth);

  // direct illumination

  for (int s = 0; s < samples_; ++s) {
    bool valid =
        material->sample(albedo, wo, normal_ff, sampler, &wi, &weight, &pdf);
    if (valid) {
      Lr = Lin * weight * ao_weight_ * (1.0f / pdf);
      if (hasPositive(Lr)) {
        // create shadow ray
        Ray *shadow = mem->Alloc<Ray>(1, false);
        CHECK_NOTNULL(shadow);

        RayUtil::makeShadow(rayin, s, pos, wi, Lr, isect.tfar, shadow);

        sq->push(shadow);
      }
    }
  }

  // indirect illumination

#ifdef SPRAY_GLOG_CHECK
  CHECK_LT(ray_depth + rayin.depth, bounces_);
#endif

  if (next_ray_depth < bounces_) {
    RandomSampler_init(sampler, rayin.samid * next_ray_depth);
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

}  // namespace ooc
}  // namespace spray
