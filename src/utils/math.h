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

#include <algorithm>

#include "glm/glm.hpp"

namespace spray {

inline glm::vec3 neg(const glm::vec3& v) { return -v; }

inline glm::vec3 neg(float v[3]) { return -glm::vec3(v[0], v[1], v[2]); }

inline float dot(const float v1[3], const float v2[3]) {
  return (v1[0] * v2[0]) + (v1[1] * v2[1]) + (v1[2] * v2[2]);
}

inline float dot(const float v1[3], const glm::vec3& v2) {
  return (v1[0] * v2[0]) + (v1[1] * v2[1]) + (v1[2] * v2[2]);
}

inline float dot(const glm::vec3& v1, const float v2[3]) {
  return (v1[0] * v2[0]) + (v1[1] * v2[1]) + (v1[2] * v2[2]);
}

inline glm::vec3 faceForward(const glm::vec3& dir, const glm::vec3& normal) {
  return glm::dot(dir, normal) < 0.0f ? normal : neg(normal);
}

inline glm::vec3 faceForwardFloat(const float dir[3], const float normal[3]) {
  glm::vec3 n(normal[0], normal[1], normal[2]);
  return dot(dir, normal) < 0.0f ? n : neg(n);
}

inline void faceForward(const float dir[3], float normal[3]) {
  if (dot(dir, normal) >= 0.0f) {
    normal[0] = -normal[0];
    normal[1] = -normal[1];
    normal[2] = -normal[2];
  }
}

// assume v is facing outward from surface
inline glm::vec3 Faceforward(const glm::vec3& v, const glm::vec3& v2) {
  return (glm::dot(v, v2) < 0.f) ? -v : v;
}

template <typename T>
inline T clamp(T x, T lo, T hi) {
  return std::min(std::max(x, lo), hi);
}

inline bool isZero(const glm::vec3& v) {
  return (v[0] == 0.f) && (v[1] == 0.f) && (v[1] == 0.f);
}

inline bool hasPositive(const glm::vec3& v) {
  return (v[0] > 0.0f || v[1] > 0.0f || v[2] > 0.0f);
}

inline bool hasNegative(const glm::vec3& v) {
  return (v[0] < 0.0f || v[1] < 0.0f || v[2] < 0.0f);
}

inline float maxComponent(const glm::vec3 v) {
  return glm::max(glm::max(v.x, v.y), v.z);
}

inline glm::vec3 normalize(const float v[3]) {
  return glm::normalize(glm::vec3(v[0], v[1], v[2]));
}

inline float squaredLength(const glm::vec3& v) {
  return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
}

}  // namespace spray

