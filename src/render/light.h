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

namespace spray {

class Light {
 public:
  virtual glm::vec3 sample(const glm::vec3& p, glm::vec3* wi,
                           float* pdf) const = 0;
  virtual ~Light() {}
  virtual bool isAreaLight() const = 0;
  virtual glm::vec3 sampleArea(RandomSampler& sampler, const glm::vec3& normal,
                               glm::vec3* wi, float* pdf) const = 0;
  virtual glm::vec3 sampleL(const glm::vec3& hit_point, RandomSampler& sampler,
                            const glm::vec3& normal, glm::vec3* wi,
                            float* pdf) const = 0;
  virtual int getNumSamples() const = 0;
};

class PointLight : public Light {
 public:
  PointLight(const glm::vec3& position, const glm::vec3& radiance)
      : position_(position), radiance_(radiance) {}
  virtual ~PointLight() {}

  glm::vec3 sample(const glm::vec3& p, glm::vec3* wi,
                   float* pdf) const override {
    glm::vec3 dir2light = position_ - p;
    *wi = glm::normalize(dir2light);
    *pdf = 1.0;
    return radiance_;
  }
  bool isAreaLight() const override { return false; }

  glm::vec3 sampleArea(RandomSampler& sampler, const glm::vec3& normal,
                       glm::vec3* wi, float* pdf) const override {
    LOG(FATAL) << "forbidden";
    return radiance_;
  }

  glm::vec3 sampleL(const glm::vec3& hit_point, RandomSampler& sampler,
                    const glm::vec3& normal, glm::vec3* wi,
                    float* pdf) const override {
    glm::vec3 dir2light = position_ - hit_point;
    *wi = glm::normalize(dir2light);
    *pdf = 1.0;
    return radiance_;
  }

  int getNumSamples() const override { return 1; }

 private:
  glm::vec3 position_;
  glm::vec3 radiance_;
};

class DiffuseHemisphereLight : public Light {
 public:
  DiffuseHemisphereLight(const glm::vec3& radiance, int num_samples)
      : radiance_(radiance), num_samples_(num_samples) {}

  virtual ~DiffuseHemisphereLight() {}

  bool isAreaLight() const override { return true; }

  glm::vec3 sample(const glm::vec3& p, glm::vec3* wi,
                   float* pdf) const override {
    CHECK(false);
    return radiance_;
  }
  glm::vec3 sampleArea(RandomSampler& sampler, const glm::vec3& normal,
                       glm::vec3* wi, float* pdf) const override {
    CHECK(false);
    glm::vec2 u = RandomSampler_get2D(sampler);
    getCosineHemisphereSample(u.x, u.y, normal, wi, pdf);
    return radiance_;
  }

  glm::vec3 sampleL(const glm::vec3& hit_point, RandomSampler& sampler,
                    const glm::vec3& normal, glm::vec3* wi,
                    float* pdf) const override {
    glm::vec2 u = RandomSampler_get2D(sampler);
    getCosineHemisphereSample(u.x, u.y, normal, wi, pdf);
    return radiance_;
  }

  int getNumSamples() const override { return num_samples_; }

 private:
  glm::vec3 radiance_;
  int num_samples_;
};

class DiffuseSphereLight : public Light {
 public:
  DiffuseSphereLight(const glm::vec3& radiance, int num_samples)
      : radiance_(radiance), num_samples_(num_samples) {}

  virtual ~DiffuseSphereLight() {}

  bool isAreaLight() const override { return true; }

  glm::vec3 sample(const glm::vec3& p, glm::vec3* wi,
                   float* pdf) const override {
    CHECK(false);
    return radiance_;
  }
  glm::vec3 sampleArea(RandomSampler& sampler, const glm::vec3& normal,
                       glm::vec3* wi, float* pdf) const override {
    CHECK(false);
    return radiance_;
  }

  glm::vec3 sampleL(const glm::vec3& hit_point, RandomSampler& sampler,
                    const glm::vec3& normal, glm::vec3* wi,
                    float* pdf) const override {
    sampleRandomInUnitSphere(hit_point, normal, wi, pdf);
    return radiance_;
  }

  int getNumSamples() const override { return num_samples_; }

 private:
  glm::vec3 radiance_;
  int num_samples_;
};

}  // namespace spray

