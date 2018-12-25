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

#include <math.h>
#include <stdlib.h>

#include "render/spray.h"
#include "utils/math.h"

namespace spray {

void createCoordSystem(const glm::vec3& n, glm::mat3* obj2world);

struct Sample3 {
  glm::vec3 dir;
  float pdf;
};

struct UniformDiskSampling {
  void operator()(float u1, float u2, float* x, float* y) {
    float r = glm::sqrt(u1);
    float theta = 2.0f * SPRAY_PI * u2;
    *x = r * glm::cos(theta);
    *y = r * glm::sin(theta);
  }

  float pdf() { return SPRAY_INV_TWOPI; }
};

struct ConcentricDiskSampling {
  void operator()(float u1, float u2, float* dx, float* dy) {
    float r, theta;
    // Map uniform random numbers to $[-1,1]^2$
    float sx = 2 * u1 - 1;
    float sy = 2 * u2 - 1;

    // Map square to $(r,\theta)$

    // Handle degeneracy at the origin
    if (sx == 0.0 && sy == 0.0) {
      *dx = 0.0;
      *dy = 0.0;
      return;
    }
    if (sx >= -sy) {
      if (sx > sy) {
        // Handle first region of disk
        r = sx;
        if (sy > 0.0)
          theta = sy / r;
        else
          theta = 8.0f + sy / r;
      } else {
        // Handle second region of disk
        r = sy;
        theta = 2.0f - sx / r;
      }
    } else {
      if (sx <= sy) {
        // Handle third region of disk
        r = -sx;
        theta = 4.0f - sy / r;
      } else {
        // Handle fourth region of disk
        r = -sy;
        theta = 6.0f + sx / r;
      }
    }
    theta *= SPRAY_PI / 4.f;
    *dx = r * glm::cos(theta);
    *dy = r * glm::sin(theta);
  }
};

inline float select(bool s, float a, float b) { return s ? a : b; }

inline glm::vec3 select(bool s, const glm::vec3& a, const glm::vec3& b) {
  return glm::vec3(select(s, a.x, b.x), select(s, a.y, b.y),
                   select(s, a.z, b.z));
}

inline glm::vec3 localToWorld(const glm::vec3& N, const glm::vec3& v) {
  glm::vec3 dx0(0, N.z, -N.y);
  glm::vec3 dx1(-N.z, 0, N.x);
  glm::vec3 dx =
      glm::normalize(select(glm::dot(dx0, dx0) > glm::dot(dx1, dx1), dx0, dx1));
  glm::vec3 dy = glm::normalize(glm::cross(N, dx));

  glm::mat3 trans(dx, dy, N);
  return (trans * v);
};

template <class DiskSamplingFunctor = ConcentricDiskSampling>
glm::vec3 cosineHemisphereSample(
    float u1, float u2,
    DiskSamplingFunctor disk_sampler = ConcentricDiskSampling()) {
  glm::vec3 v;
  disk_sampler(u1, u2, &v.x, &v.y);
  v.z = glm::sqrt(glm::max(0.f, 1.f - v.x * v.x - v.y * v.y));
  return v;
}

inline float cosineHemispherePdf(float costheta, float phi) {
  return costheta * SPRAY_INV_PI;
}

inline float cosineHemispherePdf(const glm::vec3& dir) {
  return dir.z * SPRAY_INV_PI;
}

void getCosineHemisphereSample(float u1, float u2, const glm::vec3& N,
                               Sample3* s);

void getCosineHemisphereSample(float u1, float u2, const glm::vec3& N,
                               glm::vec3* wi, float* pdf);

class CosineHemisphereSampler {
 public:
  static Sample3 get(const glm::vec2& random_sample, const glm::vec3& n) {
    glm::vec3 dir = get(random_sample);
    glm::mat3 obj2world;
    createCoordSystem(n, &obj2world);

    Sample3 s;
    s.dir = glm::normalize(obj2world * dir);
    s.pdf = getPdf(s.dir);
    return s;
  }

  static Sample3 get(float u, float v, const glm::vec3& n) {
    return get(glm::vec2(u, v), n);
  }

 private:
  static glm::vec3 cartesian(float phi, float sin_theta, float cos_theta) {
    float sin_phi, cos_phi;
#if defined(__APPLE__)
    __sincosf(phi, &sin_phi, &cos_phi);
#elif defined(SPRAY_USING_INTEL_COMPILER)
    sin_phi = sinf(phi);
    cos_phi = cosf(phi);
#else
    sincosf(phi, &sin_phi, &cos_phi);
#endif
    return glm::vec3(cos_phi * sin_theta, sin_phi * sin_theta, cos_theta);
  }

  static glm::vec3 get(const glm::vec2& s) {
    float phi = float(SPRAY_TWO_PI) * s.x;
    float cos_theta = sqrt(s.y);
    float sin_theta = sqrt(1.0f - s.y);
    return cartesian(phi, sin_theta, cos_theta);
  }

  static float getPdf(const glm::vec3& dir) { return dir.z / float(SPRAY_PI); }
  static float getPdf(float cos_theta) { return cos_theta / float(SPRAY_PI); }
};

}  // namespace spray
