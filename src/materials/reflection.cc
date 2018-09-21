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

#include "materials/reflection.h"

#include <string>
#include <utility>

#include "glm/glm.hpp"
#include "glog/logging.h"

#include "utils/math.h"
#include "renderers/spray.h"

namespace spray {

// I: incident direction, hit position (end)
// R: relfected direction, hit positoin (start)
// T: transmitted direction, hit positoin (start)

bool refract(const glm::vec3& I, const glm::vec3& N, float n1, float n2,
             glm::vec3* T) {
  float n = n1 / n2;
  float cosI = -glm::dot(N, I);
  float sinT2 = n * n * (1.0f - cosI * cosI);

  if (sinT2 > 1.0f) {
    // TIR (Total Interal Reflection)
    return false;
  }

  float cosT = glm::sqrt(1.0f - sinT2);

  *T = n * I + (n * cosI - cosT) * N;

  return true;
}

float reflectanceFresnel(const glm::vec3& I, const glm::vec3& N, float n1,
                         float n2) {
  float n = n1 / n2;
  float cosI = -glm::dot(N, I);
  float sinT2 = n * n * (1.0f - cosI * cosI);

  if (sinT2 > 1.0f) {
    // TIR (Total Interal Reflection)
    return 1.0f;
  }

  float cosT = glm::sqrt(1.0f - sinT2);

  float r_ortho = (n1 * cosI - n2 * cosT) / (n1 * cosI + n2 * cosT);
  float r_para = (n2 * cosI - n1 * cosT) / (n2 * cosI + n1 * cosT);

  return (r_ortho * r_ortho + r_para * r_para) / 2.0f;
}

float reflectanceSchlick(const glm::vec3& I, const glm::vec3& N, float n1,
                         float n2) {
  float r0 = (n1 - n2) / (n1 + n2);
  r0 *= r0;

  float cosI = -glm::dot(N, I);

  if (n1 > n2) {
    float n = n1 / n2;
    float sinT2 = n * n * (1.0f - cosI * cosI);

    if (sinT2 > 1.0f) {  // TIR
      return 1.0f;
    }

    cosI = glm::sqrt(1.0f - sinT2);
  }

  float x = 1.0f - cosI;

  return r0 + (1.0f - r0) * x * x * x * x * x;
}

}  // namespace spray
