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

#include <cstdint>
#include <iostream>
#include <limits>

#include "glm/glm.hpp"

#include "render/spray.h"

namespace spray {

class SPRAY_ALIGN(16) Aabb {
 public:
  Aabb();

  Aabb(const glm::vec3& v0, const glm::vec3& v1);

  // create aabb given triangle vertices
  Aabb(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2);

  friend std::ostream& operator<<(std::ostream& os, const Aabb& a);
  bool operator==(const Aabb& other);
  bool operator!=(const Aabb& other);

  void invalidate();
  void reset() { invalidate(); }

  void merge(const glm::vec3& point);
  void merge(const Aabb& aabb);

  glm::vec3 getVertex(unsigned int i) const;

  std::uint8_t getLongestAxis() const;
  glm::vec3 getExtent() const;
  glm::vec3 getCenter() const;
  float getArea() const;
  float getHalfArea() const;
  bool isValid() const;
  void draw(const glm::vec4& color) const;

  void splitHalf(unsigned int axis, Aabb* a, Aabb* b) const;

  bool contains(const glm::vec3& point) const;
  bool contains(const Aabb& aabb) const;

  const glm::vec3& getMin() const { return bounds[0]; }
  const glm::vec3& getMax() const { return bounds[1]; }
  const glm::vec3& getBound(std::uint8_t i) const { return bounds[i]; }
  float getEdge(std::uint8_t which_bound, std::uint8_t axis) const {
    return bounds[which_bound][axis];
  }
  void setBounds(const float aabb_min[3], const float aabb_max[3]) {
    for (int i = 0; i < 3; ++i) bounds[0][i] = aabb_min[i];
    for (int i = 0; i < 3; ++i) bounds[1][i] = aabb_max[i];
  }

  glm::vec3 bounds[2];
};

inline glm::ivec3 computeSigns(const glm::vec3& dir, glm::vec3* inv_dir) {
  *inv_dir = 1.0f / dir;
  return glm::ivec3(inv_dir->x < 0.0f, inv_dir->y < 0.0f, inv_dir->z < 0.0f);
}

inline bool intersectAabb(const Aabb& box, const glm::vec3& org,
                          const glm::vec3& dir, float* t) {
  const glm::vec3& origin = org;
  glm::vec3 inv_direction;
  glm::ivec3 sign = computeSigns(dir, &inv_direction);

  float tmin, tmax, tymin, tymax, tzmin, tzmax;
  tmin = (box.bounds[sign[0]].x - origin.x) * inv_direction.x;
  tmax = (box.bounds[1 - sign[0]].x - origin.x) * inv_direction.x;
  tymin = (box.bounds[sign[1]].y - origin.y) * inv_direction.y;
  tymax = (box.bounds[1 - sign[1]].y - origin.y) * inv_direction.y;
  if ((tmin > tymax) || (tymin > tmax)) return false;
  if (tymin > tmin) tmin = tymin;
  if (tymax < tmax) tmax = tymax;
  tzmin = (box.bounds[sign[2]].z - origin.z) * inv_direction.z;
  tzmax = (box.bounds[1 - sign[2]].z - origin.z) * inv_direction.z;
  if ((tmin > tzmax) || (tzmin > tmax)) return false;
  if (tzmin > tmin) tmin = tzmin;
  if (tzmax < tmax) tmax = tzmax;

  bool hit = (tmin < SPRAY_REAL_MAX) && (tmax > 0.0);
  *t = tmin;

  return hit;
}

inline bool intersectAabb(const Aabb& box, glm::vec3 org, glm::vec3 dir,
                          float t0, float t1, float* tmin_o, float* tmax_o) {
  glm::vec3 inv_direction;
  inv_direction[0] = 1.0f / dir[0];
  inv_direction[1] = 1.0f / dir[1];
  inv_direction[2] = 1.0f / dir[2];

  glm::ivec3 sign;
  sign[0] = (inv_direction[0] < 0.0f);
  sign[1] = (inv_direction[1] < 0.0f);
  sign[2] = (inv_direction[2] < 0.0f);

  float tmin, tmax, tymin, tymax, tzmin, tzmax;
  tmin = (box.bounds[sign[0]].x - org.x) * inv_direction.x;
  tmax = (box.bounds[1 - sign[0]].x - org.x) * inv_direction.x;
  tymin = (box.bounds[sign[1]].y - org.y) * inv_direction.y;
  tymax = (box.bounds[1 - sign[1]].y - org.y) * inv_direction.y;
  if ((tmin > tymax) || (tymin > tmax)) return false;
  if (tymin > tmin) tmin = tymin;
  if (tymax < tmax) tmax = tymax;
  tzmin = (box.bounds[sign[2]].z - org.z) * inv_direction.z;
  tzmax = (box.bounds[1 - sign[2]].z - org.z) * inv_direction.z;
  if ((tmin > tzmax) || (tzmin > tmax)) return false;
  if (tzmin > tmin) tmin = tzmin;
  if (tzmax < tmax) tmax = tzmax;

  *tmin_o = tmin;
  *tmax_o = tmax;

  return ((tmin < t1) && (tmax > t0));
}

inline bool intersectAabb(const Aabb& box, float org[3], float dir[3], float t0,
                          float t1, float* tmin_o, float* tmax_o) {
  float inv_direction[3];
  inv_direction[0] = 1.0f / dir[0];
  inv_direction[1] = 1.0f / dir[1];
  inv_direction[2] = 1.0f / dir[2];

  int sign[3];
  sign[0] = (inv_direction[0] < 0.0f);
  sign[1] = (inv_direction[1] < 0.0f);
  sign[2] = (inv_direction[2] < 0.0f);

  float tmin, tmax, tymin, tymax, tzmin, tzmax;
  tmin = (box.bounds[sign[0]].x - org[0]) * inv_direction[0];
  tmax = (box.bounds[1 - sign[0]].x - org[0]) * inv_direction[0];
  tymin = (box.bounds[sign[1]].y - org[1]) * inv_direction[1];
  tymax = (box.bounds[1 - sign[1]].y - org[1]) * inv_direction[1];
  if ((tmin > tymax) || (tymin > tmax)) return false;
  if (tymin > tmin) tmin = tymin;
  if (tymax < tmax) tmax = tymax;
  tzmin = (box.bounds[sign[2]].z - org[2]) * inv_direction[2];
  tzmax = (box.bounds[1 - sign[2]].z - org[2]) * inv_direction[2];
  if ((tmin > tzmax) || (tzmin > tmax)) return false;
  if (tzmin > tmin) tmin = tzmin;
  if (tzmax < tmax) tmax = tzmax;

  *tmin_o = tmin;
  *tmax_o = tmax;

  return ((tmin < t1) && (tmax > t0));
}

inline bool intersectAabb(const Aabb& box, const glm::vec3& org,
                          const glm::vec3& dir) {
  const glm::vec3& origin = org;
  glm::vec3 inv_direction;
  glm::ivec3 sign = computeSigns(dir, &inv_direction);

  float tmin, tmax, tymin, tymax, tzmin, tzmax;
  tmin = (box.bounds[sign[0]].x - origin.x) * inv_direction.x;
  tmax = (box.bounds[1 - sign[0]].x - origin.x) * inv_direction.x;
  tymin = (box.bounds[sign[1]].y - origin.y) * inv_direction.y;
  tymax = (box.bounds[1 - sign[1]].y - origin.y) * inv_direction.y;
  if ((tmin > tymax) || (tymin > tmax)) return false;
  if (tymin > tmin) tmin = tymin;
  if (tymax < tmax) tmax = tymax;
  tzmin = (box.bounds[sign[2]].z - origin.z) * inv_direction.z;
  tzmax = (box.bounds[1 - sign[2]].z - origin.z) * inv_direction.z;
  if ((tmin > tzmax) || (tzmin > tmax)) return false;
  if (tzmin > tmin) tmin = tzmin;
  if (tzmax < tmax) tmax = tzmax;

  bool hit = (tmin < SPRAY_REAL_MAX) && (tmax > 0.0);

  return hit;
}

}  // namespace spray

