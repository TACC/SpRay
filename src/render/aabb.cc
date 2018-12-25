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

#include "render/aabb.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <iomanip>

#include "display/opengl.h"
#include "render/spray.h"

namespace spray {

glm::vec3 fminf(const glm::vec3& a, const glm::vec3& b) {
  return glm::vec3(std::fminf(a[0], b[0]), std::fminf(a[1], b[1]),
                   std::fminf(a[2], b[2]));
}

glm::vec3 fmaxf(const glm::vec3& a, const glm::vec3& b) {
  return glm::vec3(std::fmaxf(a[0], b[0]), std::fmaxf(a[1], b[1]),
                   std::fmaxf(a[2], b[2]));
}

Aabb::Aabb() { invalidate(); }

Aabb::Aabb(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2) {
  bounds[0] = fminf(v0, fminf(v1, v2));
  bounds[1] = fmaxf(v0, fmaxf(v1, v2));
}

Aabb::Aabb(const glm::vec3& v0, const glm::vec3& v1) {
  bounds[0] = fminf(v0, v1);
  bounds[1] = fmaxf(v0, v1);
}

std::ostream& operator<<(std::ostream& os, const Aabb& a) {
  os << std::setprecision(15);
  os << "min(" << a.bounds[0].x << "," << a.bounds[0].y << "," << a.bounds[0].z
     << ") max(" << a.bounds[1].x << "," << a.bounds[1].y << ","
     << a.bounds[1].z << ")";
  return os;
}

bool Aabb::operator==(const Aabb& other) {
  return ((bounds[0].x == other.bounds[0].x) &&
          (bounds[0].y == other.bounds[0].y) &&
          (bounds[0].z == other.bounds[0].z) &&
          (bounds[1].x == other.bounds[1].x) &&
          (bounds[1].y == other.bounds[1].y) &&
          (bounds[1].z == other.bounds[1].z));
}

bool Aabb::operator!=(const Aabb& other) {
  return !((bounds[0].x == other.bounds[0].x) &&
           (bounds[0].y == other.bounds[0].y) &&
           (bounds[0].z == other.bounds[0].z) &&
           (bounds[1].x == other.bounds[1].x) &&
           (bounds[1].y == other.bounds[1].y) &&
           (bounds[1].z == other.bounds[1].z));
}

void Aabb::invalidate() {
  bounds[0] = glm::vec3(std::numeric_limits<float>::max());
  bounds[1] = -glm::vec3(std::numeric_limits<float>::max());
}

void Aabb::merge(const glm::vec3& point) {
  bounds[0] = fminf(bounds[0], point);
  bounds[1] = fmaxf(bounds[1], point);
}

void Aabb::merge(const Aabb& aabb) {
  bounds[0] = fminf(aabb.bounds[0], bounds[0]);
  bounds[1] = fmaxf(aabb.bounds[1], bounds[1]);
}

glm::vec3 Aabb::getExtent() const { return (bounds[1] - bounds[0]); }

glm::vec3 Aabb::getVertex(unsigned int i) const {
  // 0,1,2,3: bottom plane
  // 4,5,6,7: upper plane
  glm::vec3 v;
  if (i == 0) {
    v = glm::vec3(bounds[0].x, bounds[0].y, bounds[0].z);
  } else if (i == 1) {
    v = glm::vec3(bounds[1].x, bounds[0].y, bounds[0].z);
  } else if (i == 2) {
    v = glm::vec3(bounds[1].x, bounds[0].y, bounds[1].z);
  } else if (i == 3) {
    v = glm::vec3(bounds[0].x, bounds[0].y, bounds[1].z);
  } else if (i == 4) {
    v = glm::vec3(bounds[0].x, bounds[1].y, bounds[0].z);
  } else if (i == 5) {
    v = glm::vec3(bounds[1].x, bounds[1].y, bounds[0].z);
  } else if (i == 6) {
    v = glm::vec3(bounds[1].x, bounds[1].y, bounds[1].z);
  } else if (i == 7) {
    v = glm::vec3(bounds[0].x, bounds[1].y, bounds[1].z);
  }
  return v;
}

std::uint8_t Aabb::getLongestAxis() const {
  glm::vec3 diag = getExtent();

  if (diag.x > diag.y && diag.x > diag.z)
    return 0;  // x-axis
  else if (diag.y > diag.z)
    return 1;  // y-axis
  else
    return 2;  // z-axis
}

glm::vec3 Aabb::getCenter() const {
  return glm::vec3(0.5 * (bounds[0].x + bounds[1].x),
                   0.5 * (bounds[0].y + bounds[1].y),
                   0.5 * (bounds[0].z + bounds[1].z));
}

float Aabb::getArea() const { return (2.f * getHalfArea()); }

float Aabb::getHalfArea() const {
  glm::vec3 diag = getExtent();
  return ((diag.x * diag.y) + (diag.y * diag.z) + (diag.z * diag.x));
}

bool Aabb::isValid() const {
  return (bounds[0].x <= bounds[1].x && bounds[0].y <= bounds[1].y &&
          bounds[0].z <= bounds[1].z);
}

void Aabb::draw(const glm::vec4& color) const {
  glColor4f(color.r, color.g, color.b, color.a);

  // top
  glBegin(GL_LINE_STRIP);
  glVertex3d(bounds[1].x, bounds[1].y, bounds[1].z);
  glVertex3d(bounds[1].x, bounds[1].y, bounds[0].z);
  glVertex3d(bounds[0].x, bounds[1].y, bounds[0].z);
  glVertex3d(bounds[0].x, bounds[1].y, bounds[1].z);
  glVertex3d(bounds[1].x, bounds[1].y, bounds[1].z);
  glEnd();

  // bottom
  glBegin(GL_LINE_STRIP);
  glVertex3d(bounds[0].x, bounds[0].y, bounds[0].z);
  glVertex3d(bounds[0].x, bounds[0].y, bounds[1].z);
  glVertex3d(bounds[1].x, bounds[0].y, bounds[1].z);
  glVertex3d(bounds[1].x, bounds[0].y, bounds[0].z);
  glVertex3d(bounds[0].x, bounds[0].y, bounds[0].z);
  glEnd();

  // side
  glBegin(GL_LINES);
  glVertex3d(bounds[1].x, bounds[1].y, bounds[1].z);
  glVertex3d(bounds[1].x, bounds[0].y, bounds[1].z);
  glVertex3d(bounds[1].x, bounds[1].y, bounds[0].z);
  glVertex3d(bounds[1].x, bounds[0].y, bounds[0].z);
  glVertex3d(bounds[0].x, bounds[1].y, bounds[0].z);
  glVertex3d(bounds[0].x, bounds[0].y, bounds[0].z);
  glVertex3d(bounds[0].x, bounds[1].y, bounds[1].z);
  glVertex3d(bounds[0].x, bounds[0].y, bounds[1].z);
  glEnd();
}

void Aabb::splitHalf(unsigned int axis, Aabb* a, Aabb* b) const {
  // split point
  float split = (bounds[0][axis] + bounds[1][axis]) / 2.f;
  float x, y, z;

  // a
  x = (axis == 0) ? split : bounds[1][0];
  y = (axis == 1) ? split : bounds[1][1];
  z = (axis == 2) ? split : bounds[1][2];
  a->bounds[0] = bounds[0];
  a->bounds[1] = glm::vec3(x, y, z);

  // b
  x = (axis == 0) ? split : bounds[0][0];
  y = (axis == 1) ? split : bounds[0][1];
  z = (axis == 2) ? split : bounds[0][2];
  b->bounds[0] = glm::vec3(x, y, z);
  b->bounds[1] = bounds[1];
}

bool Aabb::contains(const glm::vec3& point) const {
  return point.x >= bounds[0].x && point.x <= bounds[1].x &&
         point.y >= bounds[0].y && point.y <= bounds[1].y &&
         point.z >= bounds[0].z && point.z <= bounds[1].z;
}

bool Aabb::contains(const Aabb& aabb) const {
  return contains(aabb.bounds[0]) && contains(aabb.bounds[1]);
}

}  // namespace spray
