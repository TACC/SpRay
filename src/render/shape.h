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

#include "glm/glm.hpp"
#include "glog/logging.h"

#include "render/aabb.h"
#include "render/material.h"

namespace spray {

class Shape {
 public:
  enum ShapeType { UNDEFINED, SPHERE };

  Shape(Material* m) : material(m) {}
  virtual ~Shape() { delete material; }

  virtual int type() const {
    LOG(FATAL) << "undefined";
    return UNDEFINED;
  }
  virtual void setGeomId(int id) { LOG(FATAL) << "undefined"; }

  virtual void getBounds(Aabb* aabb) const = 0;

  Material* material;
};

class Sphere : public Shape {
 public:
  Sphere(const glm::vec3& center, float radius, Material* m)
      : Shape(m), center(center), radius(radius) {}

  int type() const override { return Shape::SPHERE; }
  void setGeomId(int id) override { geom_id = id; }

  void getBounds(Aabb* aabb) const override;

  const glm::vec3 center;
  const float radius;
  unsigned int geom_id;
};

void Sphere::getBounds(Aabb* aabb) const {
  aabb->bounds[0].x = center.x - radius;
  aabb->bounds[0].y = center.y - radius;
  aabb->bounds[0].z = center.z - radius;
  aabb->bounds[1].x = center.x + radius;
  aabb->bounds[1].y = center.y + radius;
  aabb->bounds[1].z = center.z + radius;
}

}  // namespace spray
