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

#include "render/sampler.h"

namespace spray {

void createCoordSystem(const glm::vec3& n, glm::mat3* obj2world) {
  glm::vec3 z = glm::vec3(n.x, n.y, n.z);
  glm::vec3 h = z;
  if (fabs(h.x) <= fabs(h.y) && fabs(h.x) <= fabs(h.z))
    h.x = 1.0;
  else if (fabs(h.y) <= fabs(h.x) && fabs(h.y) <= fabs(h.z))
    h.y = 1.0;
  else
    h.z = 1.0;

  z = glm::normalize(z);
  glm::vec3 y = glm::cross(h, z);
  y = glm::normalize(y);
  glm::vec3 x = glm::cross(z, y);
  x = glm::normalize(x);

  (*obj2world)[0] = x;
  (*obj2world)[1] = y;
  (*obj2world)[2] = z;
}

void getCosineHemisphereSample(float u1, float u2, const glm::vec3& N,
                               Sample3* s) {
  glm::vec3 v = cosineHemisphereSample(u1, u2);  // v in local space
  v = glm::normalize(v);
  s->dir = glm::normalize(localToWorld(N, v));  // s->v in world space
  s->pdf = cosineHemispherePdf(v);
}

void getCosineHemisphereSample(float u1, float u2, const glm::vec3& N,
                               glm::vec3* wi, float* pdf) {
  glm::vec3 v = cosineHemisphereSample(u1, u2);  // v in local space
  v = glm::normalize(v);
  *wi = glm::normalize(localToWorld(N, v));  // s->v in world space
  *pdf = cosineHemispherePdf(v);
}

}  // namespace spray
