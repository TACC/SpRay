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

#include "render/shape_buffer.h"

#include <cmath>

#include <embree2/rtcore_geometry.h>
#include "glm/glm.hpp"
#include "glog/logging.h"

#include "render/rays.h"
#include "render/shape.h"
#include "utils/math.h"
#include "utils/util.h"

namespace spray {

ShapeBuffer::ShapeBuffer()
    : max_cache_size_(0),
      // colors_(nullptr),
      device_(nullptr),
      loaded_(false) {}

ShapeBuffer::~ShapeBuffer() { cleanup(); }

void ShapeBuffer::init(int max_cache_size_ndomains, std::size_t max_nvertices,
                       std::size_t max_nfaces, bool compute_normals) {
  // cleanup
  cleanup();

  // sizes
  max_cache_size_ = max_cache_size_ndomains;

  std::size_t cache_size = static_cast<std::size_t>(max_cache_size_ndomains);

  // TODO: colors
  // colors_ = arena_.Alloc<uint32_t>(cache_size * max_nvertices, false);
  // CHECK_NOTNULL(colors_);

  // embree device
  device_ = rtcNewDevice(nullptr);
  CHECK_NOTNULL(device_);

  // embree scenes
  // scenes_ = arena_.Alloc<RTCScene>(cache_size, false);
  // CHECK_NOTNULL(scenes_);
  scenes_.resize(cache_size);

  RTCSceneFlags sflags = RTC_SCENE_STATIC | RTC_SCENE_HIGH_QUALITY;
  RTCAlgorithmFlags aflags = RTC_INTERSECT1;

  for (std::size_t i = 0; i < cache_size; ++i) {
    // scenes_[i] = rtcDeviceNewScene(device_, RTC_SCENE_DYNAMIC,
    // RTC_INTERSECT1);
    scenes_[i] = rtcDeviceNewScene(device_, sflags, aflags);
    CHECK_NOTNULL(scenes_[i]);
  }

  shapes_.resize(cache_size);
  for (std::size_t i = 0; i < cache_size; ++i) {
    shapes_[i] = nullptr;
  }
}

RTCScene ShapeBuffer::load(const std::string& filename, int cache_block,
                           const glm::mat4& transform, bool apply_transform,
                           std::vector<Shape*>& shapes) {
  CHECK(!loaded_);
  loaded_ = true;

  CHECK_LT(cache_block, shapes_.size());
  CHECK(shapes_[cache_block] == nullptr);
  shapes_[cache_block] = &shapes;

  // TODO: transform
  // TODO: can we somehow not delete geometry? i.e. something similar to how we
  // handle triangles using rtcSetBuffer2?
  // TODO: support ooc (ooc not supported at this moment)

  RTCSceneFlags sflags = RTC_SCENE_STATIC | RTC_SCENE_HIGH_QUALITY;
  RTCAlgorithmFlags aflags = RTC_INTERSECT1;
  RTCScene scene = scenes_[cache_block];

  // create geometry
  unsigned int geom_id = rtcNewUserGeometry(scene, shapes.size());

  // set geom id
  void* shape_ptr = (void*)&shapes[0];
  for (std::size_t i = 0; i < shapes.size(); ++i) {
    // TODO: only spheres are supported
    Shape* shape = shapes[i];
    CHECK_EQ(shape->type(), Shape::SPHERE);
    shape->setGeomId(geom_id);
  }

  // set data
  rtcSetUserData(scene, geom_id, shape_ptr);

  // callback functions
  rtcSetBoundsFunction(scene, geom_id, ShapeBuffer::sphereBoundsCallback);
  rtcSetIntersectFunction(scene, geom_id,
                          ShapeBuffer::sphereIntersect1Callback);
  rtcSetOccludedFunction(scene, geom_id, ShapeBuffer::sphereOccluded1Callback);

  // update geometry
  rtcUpdate(scene, geom_id);
  rtcEnable(scene, geom_id);

  rtcCommit(scene);

  return scene;
}

void ShapeBuffer::cleanup() {
  // embree scenes
  for (int i = 0; i < max_cache_size_; ++i) {
    rtcDeleteScene(scenes_[i]);
  }

  // embree device
  rtcDeleteDevice(device_);
  device_ = nullptr;

  // arena
  // arena_.Reset();

  // sizes
  max_cache_size_ = 0;
}

void ShapeBuffer::sphereBoundsCallback(void* shape_ptr, std::size_t item,
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

void ShapeBuffer::sphereIntersect1Callback(void* shape_ptr, RTCRay& ray,
                                           std::size_t item) {
  const Sphere** spheres = static_cast<const Sphere**>(shape_ptr);
  const Sphere* sphere = spheres[item];

  const glm::vec3 center_to_origin(ray.org[0] - sphere->center[0],
                                   ray.org[1] - sphere->center[1],
                                   ray.org[2] - sphere->center[2]);

  float a = spray::dot(ray.dir, ray.dir);
  float b = 2.0f * spray::dot(center_to_origin, ray.dir);
  float c = glm::dot(center_to_origin, center_to_origin) -
            (sphere->radius * sphere->radius);
  float discriminant = (b * b) - (4.0f * a * c);
  if (discriminant < 0.0f) return;

  float sqrt_d = std::sqrt(discriminant);

  // TODO: use embree's rcp() for 1/a
  float a2 = 2.0f * a;
  float root0 = (-b - sqrt_d) / a2;
  float root1 = (-b + sqrt_d) / a2;

  // root0 between tnear and tfar
  if (root0 > ray.tnear && root0 < ray.tfar) {
    // TODO: update u,v with correct values
    ray.u = 0.0f;
    ray.v = 0.0f;
    ray.tfar = root0;
    ray.geomID = sphere->geom_id;
    ray.primID = static_cast<unsigned int>(item);
    // NOTE: Ng not normalized
    ray.Ng[0] = ray.org[0] + (root0 * ray.dir[0]) - sphere->center[0];
    ray.Ng[1] = ray.org[1] + (root0 * ray.dir[1]) - sphere->center[1];
    ray.Ng[2] = ray.org[2] + (root0 * ray.dir[2]) - sphere->center[2];
  }

  // root1 between tnear and tfar
  if (root1 > ray.tnear && root1 < ray.tfar) {
    // TODO: update u,v with correct values
    ray.u = 0.0f;
    ray.v = 0.0f;
    ray.tfar = root1;
    ray.geomID = sphere->geom_id;
    ray.primID = static_cast<unsigned int>(item);
    // NOTE: Ng not normalized
    ray.Ng[0] = ray.org[0] + (root1 * ray.dir[0]) - sphere->center[0];
    ray.Ng[1] = ray.org[1] + (root1 * ray.dir[1]) - sphere->center[1];
    ray.Ng[2] = ray.org[2] + (root1 * ray.dir[2]) - sphere->center[2];
  }
}

void ShapeBuffer::sphereOccluded1Callback(void* shape_ptr, RTCRay& ray,
                                          std::size_t item) {
  const Sphere** spheres = static_cast<const Sphere**>(shape_ptr);
  const Sphere* sphere = spheres[item];

  const glm::vec3 center_to_origin(ray.org[0] - sphere->center[0],
                                   ray.org[1] - sphere->center[1],
                                   ray.org[2] - sphere->center[2]);

  float a = spray::dot(ray.dir, ray.dir);
  float b = 2.0f * spray::dot(center_to_origin, ray.dir);
  float c = glm::dot(center_to_origin, center_to_origin) -
            (sphere->radius * sphere->radius);
  float discriminant = (b * b) - (4.0f * a * c);
  if (discriminant < 0.0f) return;

  float sqrt_d = std::sqrt(discriminant);

  // TODO: use embree's rcp() for 1/a
  float a2 = 2.0f * a;
  float root0 = (-b - sqrt_d) / a2;
  float root1 = (-b + sqrt_d) / a2;

  // root0 between tnear and tfar
  if (root0 > ray.tnear && root0 < ray.tfar) {
    ray.geomID = 0;  // 0 means occluded
  }

  // root1 between tnear and tfar
  if (root1 > ray.tnear && root1 < ray.tfar) {
    ray.geomID = 0;  // 0 means occluded
  }
}

// void ShapeBuffer::getColorTuple(int cache_block, uint32_t primID,
//                                 uint32_t colors[3]) const {
//   uint32_t* faces = &faces_[faceBaseIndex(cache_block)];
// 
//   uint32_t fid = primID * NUM_VERTICES_PER_FACE;
//   uint32_t* c = &colors_[colorBaseIndex(cache_block)];
// 
//   colors[0] = c[faces[fid]];
//   colors[1] = c[faces[fid + 1]];
//   colors[2] = c[faces[fid + 2]];
// }

void ShapeBuffer::updateIntersection(int cache_block,
                                     RTCRayIntersection* isect) const {
  // uint32_t r = 128;
  // uint32_t g = 128;
  // uint32_t b = 128;
  // isect->color = util::pack(r, g, b);
  isect->material = getMaterial(cache_block, isect->primID);

  // shading normal
  // isect->Ns[0] = isect->Ng[0];
  // isect->Ns[1] = isect->Ng[1];
  // isect->Ns[2] = isect->Ng[2];
}

}  // namespace spray

