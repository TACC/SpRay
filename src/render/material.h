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

#include "render/sampler.h"
#include "render/spray.h"

namespace spray {

class Material {
 public:
  enum { UNDEFINED, MATTE, METAL, DIELECTRIC };

  /** Return a material type. */
  virtual int type() const = 0;

  /**
   * Shade based on the material type.
   *
   * \param wi Normalized incident ray direction. Direction to the light source.
   * \param wo Normalized outgoing ray direction.
   * \param normal Normalized surface normal.
   * \param inv_pdf Reciprocal of a probablity density function value.
   */
  virtual glm::vec3 shade(const glm::vec3& wi, const glm::vec3& wo,
                          const glm::vec3& normal, float* inv_pdf) const = 0;

  virtual glm::vec3 sample(const glm::vec3& wo, const glm::vec3& normal,
                           RandomSampler& sampler, glm::vec3* wi,
                           float* pdf) const = 0;
};

class Matte : public Material {
 public:
  Matte(const glm::vec3& albedo) : albedo_(albedo) {}
  int type() const override { return MATTE; }

  glm::vec3 shade(const glm::vec3& wi, const glm::vec3& wo,
                  const glm::vec3& normal, float* inv_pdf) const override;

  glm::vec3 sample(const glm::vec3& wo, const glm::vec3& normal,
                   RandomSampler& sampler, glm::vec3* wi,
                   float* pdf) const override;

 private:
  const glm::vec3 albedo_;
};

inline glm::vec3 Matte::shade(const glm::vec3& wi, const glm::vec3& wo,
                              const glm::vec3& normal, float* inv_pdf) const {
  *inv_pdf = SPRAY_ONE_OVER_PI;
  return albedo_ * SPRAY_ONE_OVER_PI *
         glm::clamp(glm::dot(wi, normal), 0.0f, 1.0f);
}

inline glm::vec3 Matte::sample(const glm::vec3& wo, const glm::vec3& normal,
                               RandomSampler& sampler, glm::vec3* wi,
                               float* pdf) const {
  glm::vec2 u = RandomSampler_get2D(sampler);
  getCosineHemisphereSample(u.x, u.y, normal, wi, pdf);
  return albedo_ * SPRAY_ONE_OVER_PI *
         glm::clamp(glm::dot(*wi, normal), 0.0f, 1.0f);
}

class Metal : public Material {
 public:
  Metal(const glm::vec3& albedo, float fuzz) : albedo_(albedo), fuzz_(fuzz) {}
  int type() const override { return METAL; }

  glm::vec3 shade(const glm::vec3& wi, const glm::vec3& wo,
                  const glm::vec3& normal, float* inv_pdf) const override {
    *inv_pdf = 0.0f;
    return glm::vec3(0.0f);
  }

  glm::vec3 sample(const glm::vec3& wo, const glm::vec3& normal,
                   RandomSampler& sampler, glm::vec3* wi,
                   float* pdf) const override {
    *pdf = 1.0f;
    glm::vec3 reflect = wo - (2.0f * glm::dot(wo, normal) * normal);
    *wi = glm::normalize(reflect);
    return albedo_;
  }

 private:
  const glm::vec3 albedo_;
  const float fuzz_;
};

class Dielectric : public Material {
 public:
  Dielectric(float index) : index_(index) {}
  int type() const override { return DIELECTRIC; }

  glm::vec3 shade(const glm::vec3& wi, const glm::vec3& wo,
                  const glm::vec3& normal, float* inv_pdf) const override {
    // TODO
    glm::vec3(1.0f);
  }

  glm::vec3 sample(const glm::vec3& wo, const glm::vec3& normal,
                   RandomSampler& sampler, glm::vec3* wi,
                   float* pdf) const override {
    // TODO
    glm::vec3(1.0f);
  }

 private:
  const float index_;  ///< Refraction index.
};


}  // namespace spray
