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

#include "renderers/spray.h"

namespace spray {
namespace insitu {

struct SPRAY_ALIGN(16) Ray {
  float org[3];
  int pixid;
  float dir[3];
  int samid;
  float w[3];
  float t;
  int light;
  int occluded;
};

struct RayData {
  Ray* ray;
  float tdom;
  int dom_depth;
};

struct RayUtil {
  inline static void makeShadow(const Ray& ray, int light, const glm::vec3& pos,
                                const glm::vec3& dir, const glm::vec3& w,
                                float t, Ray* shadow) {
    shadow->org[0] = pos[0];
    shadow->org[1] = pos[1];
    shadow->org[2] = pos[2];

    shadow->pixid = ray.pixid;

    shadow->dir[0] = dir[0];
    shadow->dir[1] = dir[1];
    shadow->dir[2] = dir[2];

    shadow->samid = ray.samid;

    shadow->w[0] = w[0];
    shadow->w[1] = w[1];
    shadow->w[2] = w[2];

    shadow->t = t;

    shadow->light = light;
    shadow->occluded = 0;
  }
  inline static void makeRay(const Ray& rayin, const glm::vec3& pos,
                             const glm::vec3& dir, const glm::vec3& w, float t,
                             Ray* rayout) {
    rayout->org[0] = pos[0];
    rayout->org[1] = pos[1];
    rayout->org[2] = pos[2];

    rayout->pixid = rayin.pixid;

    rayout->dir[0] = dir[0];
    rayout->dir[1] = dir[1];
    rayout->dir[2] = dir[2];

    rayout->samid = rayin.samid;

    rayout->w[0] = w[0];
    rayout->w[1] = w[1];
    rayout->w[2] = w[2];

    rayout->t = t;
  }
};

}  // namespace insitu
}  // namespace spray
