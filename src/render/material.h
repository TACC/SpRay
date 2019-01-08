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

#include <cmath>
#include <cstdlib>
#include <iostream>

#include "embree/random_sampler.h"
#include "glm/glm.hpp"
#include "glog/logging.h"

#include "render/sampler.h"
#include "render/spray.h"
#include "utils/util.h"

namespace spray {

class Material {
 public:
  enum { UNDEFINED, MATTE, METAL, DIELECTRIC };

  virtual ~Material() {}

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
                          const glm::vec3& normal) const = 0;

  virtual bool sample(const glm::vec3& wo, const glm::vec3& normal,
                      RandomSampler& sampler, glm::vec3* wi, glm::vec3* color,
                      float* pdf) const = 0;

  virtual bool hasDiffuse() const { return false; }
};

class Matte : public Material {
 public:
  Matte() : albedo_(glm::vec3(0.5f, 0.5f, 0.5f)) {}
  Matte(const glm::vec3& albedo) : albedo_(albedo) {}
  int type() const override { return MATTE; }
  bool hasDiffuse() const override { return true; }

  glm::vec3 shade(const glm::vec3& wi, const glm::vec3& wo,
                  const glm::vec3& normal) const override;

  bool sample(const glm::vec3& wo, const glm::vec3& normal,
              RandomSampler& sampler, glm::vec3* wi, glm::vec3* color,
              float* pdf) const override;

 private:
  const glm::vec3 albedo_;
};

inline glm::vec3 Matte::shade(const glm::vec3& wi, const glm::vec3& wo,
                              const glm::vec3& normal) const {
  return albedo_ * glm::clamp(glm::dot(wi, normal), 0.0f, 1.0f);
}

inline bool Matte::sample(const glm::vec3& wo, const glm::vec3& normal,
                          RandomSampler& sampler, glm::vec3* wi,
                          glm::vec3* color, float* pdf) const {
  glm::vec2 u = RandomSampler_get2D(sampler);
  getCosineHemisphereSample(u.x, u.y, normal, wi, pdf);

  if (*pdf > 0.0f) {
    float l_dot_n = glm::dot(*wi, normal);

    if (l_dot_n > 0.0f) {
      *color = albedo_ * SPRAY_ONE_OVER_PI * glm::clamp(l_dot_n, 0.0f, 1.0f);
      return true;
    }
  }
  return false;
}

class Metal : public Material {
 public:
  Metal() : albedo_(glm::vec3(0.5f, 0.5f, 0.5f)), fuzz_(0.0f) {}
  Metal(const glm::vec3& albedo, float fuzz) : albedo_(albedo), fuzz_(fuzz) {}
  int type() const override { return METAL; }

  glm::vec3 shade(const glm::vec3& wi, const glm::vec3& wo,
                  const glm::vec3& normal) const override {
    return glm::vec3(0.0f);
    // glm::vec3 half_hat = glm::normalize(wi + wo);
    // float n_dot_h = glm::clamp(glm::dot(normal, half_hat), 0.0f, 1.0f);

    // float shininess = 20;
    // glm::vec3 cs = albedo_ * glm::pow(n_dot_h, shininess);
    // // glm::vec3 cd = kd * costheta;

    // return cs;
  }

  bool sample(const glm::vec3& wo, const glm::vec3& normal,
              RandomSampler& sampler, glm::vec3* wi, glm::vec3* color,
              float* pdf) const override {
    *pdf = 1.0f;
    glm::vec3 reflect = glm::reflect(-wo, normal);

    if (glm::dot(reflect, normal) > 0) {
      *wi = glm::normalize(reflect);
      *color = albedo_;
      return true;
    }
    return false;
  }

 private:
  const glm::vec3 albedo_;
  const float fuzz_;
};

class Dielectric : public Material {
 public:
  Dielectric() : index_(1.5f) {}
  Dielectric(float index) : index_(index) {}

  float getIndex() const { return index_; }

  int type() const override { return DIELECTRIC; }

  glm::vec3 shade(const glm::vec3& wi, const glm::vec3& wo,
                  const glm::vec3& normal) const override {
    return glm::vec3(0.0f);
  }

  bool sample(const glm::vec3& wo, const glm::vec3& normal,
              RandomSampler& sampler, glm::vec3* wi, glm::vec3* color,
              float* pdf) const override {
    *color = glm::vec3(1.0f);
    *pdf = 1.0f;

    glm::vec3 ray_dir = -wo;
    glm::vec3 outward_normal;
    float ni_over_nt;
    float cosine;

    if (glm::dot(ray_dir, normal) > 0.0f) {
      outward_normal = -normal;
      ni_over_nt = index_;
      // cosine = index_ * glm::dot(ray_dir, normal);
      cosine = glm::dot(ray_dir, normal);
      cosine = glm::sqrt(1.0f - index_ * index_ * (1.0f - cosine * cosine));

    } else {
      outward_normal = normal;
      ni_over_nt = 1.0f / index_;
      cosine = -glm::dot(ray_dir, normal);
    }

    glm::vec3 refracted;
    float reflect_prob;

    if (refract(ray_dir, outward_normal, ni_over_nt, &refracted)) {
      reflect_prob = schlick(cosine, index_);
    } else {
      reflect_prob = 1.0;
    }

    if (drand48() < reflect_prob) {
      glm::vec3 reflected = glm::reflect(ray_dir, normal);
      *wi = glm::normalize(reflected);
    } else {
      *wi = glm::normalize(refracted);
    }

    return true;
  }

 private:
  float schlick(float cosine, float index) const {
    float r0 = (1.0f - index) / (1.0f + index);
    r0 = r0 * r0;
    float c = 1.0f - cosine;
    return r0 + ((1.0f - r0) * (c * c * c * c * c));
  }

  bool refract(const glm::vec3& v, const glm::vec3& n, float ni_over_nt,
               glm::vec3* refracted) const {
    float dt = glm::dot(v, n);
    float discriminant = 1.0f - ni_over_nt * ni_over_nt * (1.0f - dt * dt);

    if (discriminant > 0.0f) {
      *refracted = ni_over_nt * (v - n * dt) - n * std::sqrt(discriminant);
      return true;
    }
    return false;
  }

 private:
  const float index_;  ///< Refraction index.
};

class HybridMaterial {
 public:
  int getType() const { return type_; }

  void set(const Material* material, uint32_t color) {
    type_ = material->type();
    if (type_ == Material::MATTE) {
      spray::util::unpack(color, &albedo_);

    } else if (type_ == Material::METAL) {
      spray::util::unpack(color, &albedo_);

    } else if (type_ == Material::DIELECTRIC) {
      index_ = static_cast<const Dielectric*>(material)->getIndex();
    }
  }

 private:
  int type_;
  glm::vec3 albedo_;
  float index_;
};

}  // namespace spray
