// ======================================================================== //
// Copyright 2009-2018 Intel Corporation                                    //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

#pragma once

/*! \brief utility library containing sampling functions */

// convention is to return the sample (Vec3fa) generated from given Vec2f
// 's'ample as last parameter sampling functions often come in pairs: sample and
// pdf (needed later for MIS) good reference is "Total Compendium" by Philip
// Dutre http://people.cs.kuleuven.be/~philip.dutre/GI/

// adapted to SpRay

// #include "../math/vec.h"

#include <math.h>
#include <algorithm>

#include "glm/glm.hpp"
#include "glm/mat3x3.hpp"

#include "render/spray.h"

namespace spray {

inline float sqr(const float f) { return f * f; }
inline float cos2sin(const float f) {
  return sqrt(std::max(0.f, 1.f - sqr(f)));
}
inline float sin2cos(const float f) { return cos2sin(f); }

inline float rcp(const float x) { 1.0f / x; }

///////////////////////////////////////////////////////////

inline glm::vec3 cartesian(const float phi, const float sinTheta,
                           const float cosTheta) {
  float sinPhi, cosPhi;
  sincosf(phi, &sinPhi, &cosPhi);
  return glm::vec3(cosPhi * sinTheta, sinPhi * sinTheta, cosTheta);
}

inline glm::vec3 cartesian(const float phi, const float cosTheta) {
  return cartesian(phi, cos2sin(cosTheta), cosTheta);
}

/// cosine-weighted sampling of hemisphere oriented along the +z-axis
////////////////////////////////////////////////////////////////////////////////

inline glm::vec3 cosineSampleHemisphere(const glm::vec2 s) {
  const float phi = float(SPRAY_TWO_PI) * s.x;
  const float cosTheta = sqrt(s.y);
  const float sinTheta = sqrt(1.0f - s.y);
  return cartesian(phi, sinTheta, cosTheta);
}

inline float cosineSampleHemispherePDF(const glm::vec3 &dir) {
  return dir.z * SPRAY_ONE_OVER_PI;
}

inline float cosineSampleHemispherePDF(float cosTheta) {
  return cosTheta * SPRAY_ONE_OVER_PI;
}

/// power cosine-weighted sampling of hemisphere oriented along the +z-axis
////////////////////////////////////////////////////////////////////////////////

inline glm::vec3 powerCosineSampleHemisphere(const float n,
                                             const glm::vec2 &s) {
  const float phi = float(SPRAY_TWO_PI) * s.x;
  const float cosTheta = pow(s.y, 1.0f / (n + 1.0f));
  return cartesian(phi, cosTheta);
}

inline float powerCosineSampleHemispherePDF(
    const float cosTheta, const float n)  // TODO: order of arguments
{
  return (n + 1.0f) * (0.5f / float(SPRAY_PI)) * pow(cosTheta, n);
}

inline float powerCosineSampleHemispherePDF(
    const glm::vec3 &dir, const float n)  // TODO: order of arguments
{
  return (n + 1.0f) * (0.5f / float(SPRAY_PI)) * pow(dir.z, n);
}

/// sampling of cone of directions oriented along the +z-axis
////////////////////////////////////////////////////////////////////////////////

inline glm::vec3 uniformSampleCone(const float cosAngle, const glm::vec2 &s) {
  const float phi = float(SPRAY_TWO_PI) * s.x;
  const float cosTheta = 1.0f - s.y * (1.0f - cosAngle);
  return cartesian(phi, cosTheta);
}

inline float uniformSampleConePDF(const float cosAngle) {
  return rcp(float(SPRAY_TWO_PI) * (1.0f - cosAngle));
}

inline float _uniformSampleConePDF(const float cosAngle) {
  return rcp(float(SPRAY_TWO_PI) * (1.0f - cosAngle));
}

/// sampling of disk
////////////////////////////////////////////////////////////////////////////////

inline glm::vec3 uniformSampleDisk(const float radius, const glm::vec2 &s) {
  const float r = sqrtf(s.x) * radius;
  const float phi = float(SPRAY_TWO_PI) * s.y;
  float sinPhi, cosPhi;
  sincosf(phi, &sinPhi, &cosPhi);
  return glm::vec3(r * cosPhi, r * sinPhi, 0.f);
}

inline float uniformSampleDiskPDF(const float radius) {
  return rcp(float(SPRAY_PI) * sqr(radius));
}

inline float _uniformSampleDiskPDF(const float radius) {
  return rcp(float(SPRAY_PI) * sqr(radius));
}

/// sampling of triangle abc
////////////////////////////////////////////////////////////////////////////////

inline glm::vec3 uniformSampleTriangle(const glm::vec3 &a, const glm::vec3 &b,
                                       const glm::vec3 &c, const glm::vec2 &s) {
  const float su = sqrtf(s.x);
  return c + (1.0f - su) * (a - c) + (s.y * su) * (b - c);
}

inline float uniformSampleTrianglePDF(const glm::vec3 &a, const glm::vec3 &b,
                                      const glm::vec3 &c) {
  return 2.0f * rcp(abs(length(cross(a - c, b - c))));
}

struct Sample3f {
  glm::vec3 v;
  float pdf;
};

/* constructs a coordinate frame form a normalized normal */
inline glm::mat3x3 frame(const glm::vec3 &N) {
  const glm::vec3 dx0(0, N.z, -N.y);
  const glm::vec3 dx1(-N.z, 0, N.x);
  const glm::vec3 dx =
      glm::normalize(dot(dx0, dx0) > dot(dx1, dx1) ? dx0 : dx1);
  const glm::vec3 dy = normalize(cross(N, dx));
  return glm::mat3x3(dx, dy, N);
}

/*! Cosine weighted hemisphere sampling. Up direction is provided as argument.
 */
inline Sample3f cosineSampleHemisphere(const float u, const float v,
                                       const glm::vec3 &N) {
  glm::vec3 localDir = cosineSampleHemisphere(glm::vec2(u, v));
  Sample3f s;
  s.v = frame(N) * localDir;
  s.pdf = cosineSampleHemispherePDF(localDir);
  return s;
}

}  // namespace spray
