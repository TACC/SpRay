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

#include "insitu/insitu_ray.h"
#include "render/config.h"
#include "render/light.h"
#include "render/material.h"
#include "render/rays.h"
#include "render/reflection.h"
#include "utils/util.h"

#include "pbrt/memory.h"

namespace spray {
namespace insitu {

template <typename CacheT, typename SceneT>
class ShaderPtShapes {
 public:
  void init(const spray::Config &cfg, SceneT *scene) {
    bounces_ = cfg.bounces;
    samples_ = cfg.ao_samples;  // number of samples for area lights
    ks_ = cfg.ks;
    shininess_ = cfg.shininess;
    scene_ = scene;
    lights_ = scene->getLights();  // copy lights
  }

 private:
  SceneT *scene_;
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
                  std::queue<Ray *> *rq, int ray_depth);
};

template <typename CacheT, typename SceneT>
void ShaderPtShapes<CacheT, SceneT>::operator()(
    int domain_id, const Ray &rayin, const spray::RTCRayIntersection &isect,
    spray::MemoryArena *mem, std::queue<Ray *> *sq, std::queue<Ray *> *rq,
    int ray_depth) {
  glm::vec3 pos = RTCRayUtil::hitPosition(rayin.org, rayin.dir, isect.tfar);

  glm::vec3 surf_radiance;

  Material* material = isect.material;

  // if (isect.material->type() == Material::MATTE) {
  //   surf_radiance = static_cast<Matte *>(isect.material)->albedo;
  // }
  // glm::vec3 surf_radiance;
  // util::unpack(isect.color, surf_radiance);

  glm::vec3 normal(isect.Ng[0], isect.Ng[1], isect.Ng[2]);

  glm::vec3 wo(-rayin.dir[0], -rayin.dir[1], -rayin.dir[2]);
  wo = glm::normalize(wo);

  glm::vec3 Lin(rayin.w[0], rayin.w[1], rayin.w[2]);

  // float cos_theta_i = glm::dot(wo, normal);
  // bool entering = (cos_theta_i > 0.0f);
  // glm::vec3 normal_ff = entering ? normal : -normal;
  // normal_ff = glm::normalize(normal_ff);
  glm::vec3 normal_ff = glm::normalize(normal);

  glm::vec3 wi, light_color, Lr;
  float pdf, inv_shade_pdf, costheta;
  int nlights = lights_.size();

  std::size_t color_idx = sq->size();

  Bsdf *bsdf = scene_->getBsdf(domain_id);
  // bool delta_dist = bsdf->isDelta();

  int next_ray_depth = ray_depth + 1;

  RandomSampler sampler;
  RandomSampler_init(sampler, rayin.samid * next_ray_depth);

  int light_sample_offset = 0;
  int num_light_samples;
  int light_sample_id;

  glm::vec3 shade_color;

  // direct illumination

  // for each light
  if (!material->isMetal()) {
    for (int l = 0; l < nlights; ++l) {
      Light *light = lights_[l];

      num_light_samples = light->getNumSamples();

      // for each light sample
      for (int s = 0; s < num_light_samples; ++s) {
        // light color
        light_color = light->sampleL(pos, sampler, normal_ff, &wi, &pdf);

        if (pdf > 0.0f) {
          // wi, wo, normal_ff: all normalized
          shade_color = material->shade(wi, wo, normal_ff);

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
  CHECK_LT(ray_depth, bounces_);
#endif

  if (next_ray_depth < bounces_) {
    RandomSampler_init(sampler, rayin.samid * next_ray_depth);
    glm::vec3 weight;
    bool valid = material->sample(wo, normal_ff, sampler, &wi, &weight, &pdf);
    if (valid) {
      Lr = Lin * weight * (1.0f / pdf);
      if (hasPositive(Lr)) {
        Ray *r2 = mem->Alloc<Ray>(1, false);
        CHECK_NOTNULL(r2);
        RayUtil::makeRay(rayin, pos, wi, Lr, isect.tfar, r2);
        rq->push(r2);
      }
    }
  }
}

}  // namespace insitu
}  // namespace spray
