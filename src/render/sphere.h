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

namespace spray {

void computeSphereBounds(void* shape_ptr, std::size_t item,
                         RTCBounds& bounds_o) {
  const Sphere** spheres = static_cast<const Sphere**>(shape_ptr);
  const Sphere* sphere = spheres[item];
  bounds_o.lower_x = sphere->center.x - sphere->radius;
  bounds_o.lower_y = sphere->center.y - sphere->radius;
  bounds_o.lower_z = sphere->center.z - sphere->radius;
  bounds_o.upper_x = sphere->center.x + sphere->radius;
  bounds_o.upper_y = sphere->center.y + sphere->radius;
  bounds_o.upper_z = sphere->center.z + sphere->radius;
}

void raySphereIntersectionTest(void* shape_ptr, RTCRay& ray_i,
                               std::size_t item) {
  const Sphere** spheres = static_cast<const Sphere**>(shape_ptr);
  const Sphere* sphere = spheres[item];

  RTCRayIntersection& ray = (RTCRayIntersection&)ray_i;

  const glm::vec3 center_to_origin(ray.org[0] - sphere->center[0],
                                   ray.org[1] - sphere->center[1],
                                   ray.org[2] - sphere->center[2]);

  float a = spray::dot(ray.dir, ray.dir);
  float b = 2.0f * spray::dot(center_to_origin, ray.dir);
  float c = glm::dot(center_to_origin, center_to_origin) -
            (sphere->radius * sphere->radius);

  float discriminant = (b * b) - (4.0f * a * c);

  if (discriminant > 0.0f) {
    float sqrt_d = std::sqrt(discriminant);

    // TODO: use embree's rcp() for 1/a
    float a2 = 2.0f * a;
    float root = (-b - sqrt_d) / a2;

    // root between tnear and tfar
    if (root > ray.tnear && root < ray.tfar) {
      // TODO: update u,v with correct values
      // ray.u = 0.0f;
      // ray.v = 0.0f;
      ray.tfar = root;
      ray.geomID = sphere->geom_id;
      ray.primID = static_cast<unsigned int>(item);
      // NOTE: Ng not normalized
      ray.Ng[0] = ray.org[0] + (root * ray.dir[0]) - sphere->center[0];
      ray.Ng[1] = ray.org[1] + (root * ray.dir[1]) - sphere->center[1];
      ray.Ng[2] = ray.org[2] + (root * ray.dir[2]) - sphere->center[2];
      return;
    }

    root = (-b + sqrt_d) / a2;
    // root between tnear and tfar
    if (root > ray.tnear && root < ray.tfar) {
      // TODO: update u,v with correct values
      // ray.u = 0.0f;
      // ray.v = 0.0f;
      ray.tfar = root;
      ray.geomID = sphere->geom_id;
      ray.primID = static_cast<unsigned int>(item);
      // NOTE: Ng not normalized
      ray.Ng[0] = ray.org[0] + (root * ray.dir[0]) - sphere->center[0];
      ray.Ng[1] = ray.org[1] + (root * ray.dir[1]) - sphere->center[1];
      ray.Ng[2] = ray.org[2] + (root * ray.dir[2]) - sphere->center[2];
    }
  }
}

void raySphereOcclusionTest(void* shape_ptr, RTCRay& ray_i, std::size_t item) {
  const Sphere** spheres = static_cast<const Sphere**>(shape_ptr);
  const Sphere* sphere = spheres[item];

  RTCRayIntersection& ray = (RTCRayIntersection&)ray_i;

  const glm::vec3 center_to_origin(ray.org[0] - sphere->center[0],
                                   ray.org[1] - sphere->center[1],
                                   ray.org[2] - sphere->center[2]);

  float a = spray::dot(ray.dir, ray.dir);
  float b = 2.0f * spray::dot(center_to_origin, ray.dir);
  float c = glm::dot(center_to_origin, center_to_origin) -
            (sphere->radius * sphere->radius);

  float discriminant = (b * b) - (4.0f * a * c);

  if (discriminant > 0.0f) {
    float sqrt_d = std::sqrt(discriminant);

    // TODO: use embree's rcp() for 1/a
    float a2 = 2.0f * a;
    float root = (-b - sqrt_d) / a2;

    // root between tnear and tfar
    if (root > ray.tnear && root < ray.tfar) {
      ray.geomID = 0;  // 0 means occluded
      return;
    }

    root = (-b + sqrt_d) / a2;

    // root between tnear and tfar
    if (root > ray.tnear && root < ray.tfar) {
      ray.geomID = 0;  // 0 means occluded
    }
  }
}

}  // namespace spray
