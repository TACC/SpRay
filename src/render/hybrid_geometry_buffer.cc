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

#include "render/hybrid_geometry_buffer.h"

#include <cmath>

#include <embree2/rtcore_geometry.h>

#include "glm/glm.hpp"
#include "glog/logging.h"

#include "render/rays.h"
#include "render/shape.h"
#include "render/sphere.h"
#include "utils/math.h"
#include "utils/util.h"

#define DEBUG_MESH
#undef DEBUG_MESH

namespace spray {

HybridGeometryBuffer::HybridGeometryBuffer()
    : max_cache_size_(0),
      max_nvertices_(0),
      max_nfaces_(0),
      vertices_(nullptr),
      normals_(nullptr),
      faces_(nullptr),
      colors_(nullptr),
      num_vertices_(nullptr),
      num_faces_(nullptr),
      device_(nullptr),
      scenes_(nullptr),
      embree_mesh_created_(nullptr),
      // embree_mesh_geom_ids_(nullptr),
      shape_created_(nullptr),
      // shape_geom_ids_(nullptr),
      shape_geom_ids_(nullptr),
      compute_normals_(false) {}

HybridGeometryBuffer::~HybridGeometryBuffer() { cleanup(); }

void HybridGeometryBuffer::init(int max_cache_size_ndomains,
                                std::size_t max_nvertices,
                                std::size_t max_nfaces, bool compute_normals) {
  // cleanup
  cleanup();

  compute_normals_ = compute_normals;

  // sizes
  max_cache_size_ = max_cache_size_ndomains;
  max_nvertices_ = max_nvertices;
  max_nfaces_ = max_nfaces;

  std::size_t cache_size = static_cast<std::size_t>(max_cache_size_ndomains);

  // number of vertices for each domain
  num_vertices_ = arena_.Alloc<std::size_t>(cache_size, false);
  CHECK_NOTNULL(num_vertices_);

  // vertices
  vertices_ = arena_.Alloc<float>(cache_size * max_nvertices * 3, false);
  CHECK_NOTNULL(vertices_);

  // per-vertex normals
  if (compute_normals) {
    normals_ = arena_.Alloc<float>(cache_size * max_nvertices * 3, false);
    CHECK_NOTNULL(normals_);
  } else {
    normals_ = nullptr;
  }

  // number of faces for each domain
  num_faces_ = arena_.Alloc<std::size_t>(cache_size, false);
  CHECK_NOTNULL(num_faces_);

  // faces
  faces_ = arena_.Alloc<uint32_t>(
      cache_size * max_nfaces * NUM_VERTICES_PER_FACE, false);
  CHECK_NOTNULL(faces_);

  // colors
  colors_ = arena_.Alloc<uint32_t>(cache_size * max_nvertices, false);
  CHECK_NOTNULL(colors_);

  // embree mesh created
  embree_mesh_created_ = arena_.Alloc<int>(cache_size);
  CHECK_NOTNULL(embree_mesh_created_);

  // embree_mesh_geom_ids_ = arena_.Alloc<unsigned int>(cache_size);
  // CHECK_NOTNULL(embree_mesh_geom_ids_);

  // shape created
  shape_created_ = arena_.Alloc<int>(cache_size);
  CHECK_NOTNULL(shape_created_);

  shape_geom_ids_ = arena_.Alloc<unsigned int>(cache_size);
  CHECK_NOTNULL(shape_geom_ids_);

  // shape_geom_ids_ = arena_.Alloc<int>(cache_size);
  // CHECK_NOTNULL(shape_geom_ids_);

  for (std::size_t i = 0; i < cache_size; ++i) {
    embree_mesh_created_[i] = DESTROYED;
    shape_created_[i] = DESTROYED;
  }

  // embree device
  device_ = rtcNewDevice("tri_accel=bvh4.triangle4v,threads=1");
  CHECK_NOTNULL(device_);

  // embree scenes
  scenes_ = arena_.Alloc<RTCScene>(cache_size, false);
  CHECK_NOTNULL(scenes_);

  // RTCSceneFlags sflags = RTC_SCENE_DYNAMIC;
  RTCSceneFlags sflags = RTC_SCENE_STATIC;
  RTCAlgorithmFlags aflags = RTC_INTERSECT1;

  for (std::size_t i = 0; i < cache_size; ++i) {
    scenes_[i] = rtcDeviceNewScene(device_, sflags, aflags);
    CHECK_NOTNULL(scenes_[i]);
  }

  shapes_.resize(cache_size);
  for (std::size_t i = 0; i < cache_size; ++i) {
    shapes_[i] = nullptr;
  }
}

RTCScene HybridGeometryBuffer::load(const std::string& filename,
                                    int cache_block, const glm::mat4& transform,
                                    bool apply_transform,
                                    std::vector<Shape*>& shapes) {
  CHECK(!(filename.empty() && shapes.empty()));

  bool has_mash = !filename.empty();

  unsigned int shape_geom_id = 0;

  if (!filename.empty()) {
    shape_geom_id = 1;
    loadTriangles(filename, cache_block, transform, apply_transform);
  }

  std::cout << "cache_block " << cache_block << " [shape_geom_id "
            << shape_geom_id << "\n";
  shape_geom_ids_[cache_block] = shape_geom_id;

  if (!shapes.empty()) loadShapes(shapes, cache_block, shape_geom_id);

  RTCScene scene = scenes_[cache_block];
  rtcCommit(scene);

  return scene;
}

void HybridGeometryBuffer::loadTriangles(const std::string& filename,
                                         int cache_block,
                                         const glm::mat4& transform,
                                         bool apply_transform) {
  // setup
  PlyLoader::Data d;
  d.vertices_capacity = max_nvertices_ * 3;                // in
  d.faces_capacity = max_nfaces_ * NUM_VERTICES_PER_FACE;  // in
  // d.colors_capacity = 0;                                         // in
  d.colors_capacity = max_nvertices_;                     // in
  d.vertices = &vertices_[vertexBaseIndex(cache_block)];  // in/out
  // num_vertices;  // out
  d.faces = &faces_[faceBaseIndex(cache_block)];  // in/out
  // num_faces;  // out
  // d.colors = nullptr;  // rgb, in/out
  d.colors = &colors_[colorBaseIndex(cache_block)];  // rgb, in/out

  // load
  loader_.load(filename, &d);

  // update geometry sizes
  num_vertices_[cache_block] = d.num_vertices;
  num_faces_[cache_block] = d.num_faces;

  // glm::vec3 origin(0.0f);

  if (apply_transform) {
    glm::mat4 x = transform;
    glm::vec4 v;
    std::size_t nverts = d.num_vertices * 3;

    for (std::size_t n = 0; n < nverts; n += 3) {
      v = x *
          glm::vec4(d.vertices[n], d.vertices[n + 1], d.vertices[n + 2], 1.0f);
      d.vertices[n] = v.x;
      d.vertices[n + 1] = v.y;
      d.vertices[n + 2] = v.z;
    }
    // v = x * glm::vec4(origin, 1.0f);
    // origin.x = v.x;
    // origin.y = v.y;
    // origin.z = v.z;
  }

  if (compute_normals_) {
    computeNormals(cache_block);
  }

  // map buffers
  mapEmbreeBuffer(cache_block, d.vertices, d.num_vertices, d.faces,
                  d.num_faces);
}

void HybridGeometryBuffer::loadShapes(std::vector<Shape*>& shapes,
                                      int cache_block,
                                      unsigned int shape_geom_id) {
  RTCScene scene = scenes_[cache_block];
  if (shape_created_[cache_block] == DESTROYED) {
    CHECK_LT(cache_block, shapes_.size());
    CHECK(shapes_[cache_block] == nullptr);
    shapes_[cache_block] = &shapes;

    std::cout << "shapes.size: " << shapes.size()
              << " shapes pointer: " << &shapes << " cache block "
              << cache_block << std::endl;

    // std::cout<<"shapes.size: "<<shapes.size()<<std::endl;
    // create geometry
    unsigned int geom_id = rtcNewUserGeometry(scene, shapes.size());
    shape_created_[cache_block] = CREATED;
    // shape_geom_ids_[cache_block] = geom_id;

    CHECK_EQ(geom_id, shape_geom_id);

    // set geom id
    Shape* shape = shapes[shape_geom_id];
    // TODO: only spheres are supported
    CHECK_EQ(shape->type(), Shape::SPHERE);
    shape->setGeomId(shape_geom_id);
  }

  // unsigned int gid = shape_geom_ids_[cache_block];

  // set data
  void* shape_ptr = (void*)&shapes[0];
  for (std::size_t i = 0; i < shapes.size(); ++i) {
    Shape* shape = shapes[i];
    // TODO: only spheres are supported
    CHECK_EQ(shape->type(), Shape::SPHERE);
    shape->setGeomId(shape_geom_id);
  }

  rtcSetUserData(scene, shape_geom_id, shape_ptr);

  // callback functions
  rtcSetBoundsFunction(scene, shape_geom_id, computeSphereBounds);
  rtcSetIntersectFunction(scene, shape_geom_id, raySphereIntersectionTest);
  rtcSetOccludedFunction(scene, shape_geom_id, raySphereOcclusionTest);

  // update geometry
  rtcUpdate(scene, shape_geom_id);
  rtcEnable(scene, shape_geom_id);
}

void HybridGeometryBuffer::cleanup() {
  // embree scenes
  for (int i = 0; i < max_cache_size_; ++i) {
    rtcDeleteScene(scenes_[i]);
  }

  // embree device
  rtcDeleteDevice(device_);
  device_ = nullptr;

  // arena
  arena_.Reset();

  // sizes
  max_cache_size_ = 0;
  max_nvertices_ = 0;
  max_nfaces_ = 0;
}

void HybridGeometryBuffer::mapEmbreeBuffer(int cache_block, float* vertices,
                                           std::size_t num_vertices,
                                           uint32_t* faces,
                                           std::size_t num_faces) {
  // select scene
  RTCScene scene = scenes_[cache_block];

  // create triangle mesh
  if (embree_mesh_created_[cache_block] == DESTROYED) {
    unsigned int geom_id =
        rtcNewTriangleMesh(scene, RTC_GEOMETRY_DYNAMIC, num_faces, num_vertices,
                           1 /*numTimeSteps*/);

// #ifdef DEBUG_MESH
#ifdef DEBUG_MESH
    LOG(INFO) << "created embree triangle mesh geom ID: " << geom_id
              << " cache block " << cache_block;
#endif
    CHECK_EQ(geom_id, 0);
    // CHECK_NE(geom_id, RTC_INVALID_GEOMETRY_ID);

    embree_mesh_created_[cache_block] = CREATED;
    // embree_mesh_geom_ids_[cache_block] = geom_id;
  }

  // map vertices

  // gid = embree_mesh_geom_ids_[cache_block];

  rtcSetBuffer2(scene, 0 /* geomID */, RTC_VERTEX_BUFFER, vertices, 0,
                sizeof(float) * 3, num_vertices);

  // map faces

  rtcSetBuffer2(scene, 0 /* geomID */, RTC_INDEX_BUFFER, faces, 0,
                sizeof(uint32_t) * NUM_VERTICES_PER_FACE, num_faces);

  rtcUpdate(scene, 0 /*geomID*/);
  rtcEnable(scene, 0 /*geomID*/);
}

void HybridGeometryBuffer::getColorTuple(int cache_block, uint32_t primID,
                                         uint32_t colors[3]) const {
  uint32_t* faces = &faces_[faceBaseIndex(cache_block)];

  uint32_t fid = primID * NUM_VERTICES_PER_FACE;
  uint32_t* c = &colors_[colorBaseIndex(cache_block)];

  colors[0] = c[faces[fid]];
  colors[1] = c[faces[fid + 1]];
  colors[2] = c[faces[fid + 2]];
}

void HybridGeometryBuffer::getNormalTuple(int cache_block, uint32_t primID,
                                          float normals_out[9]) const {
  std::size_t vid[3];
  const uint32_t fid = primID * NUM_VERTICES_PER_FACE;
  const uint32_t* faces = &faces_[faceBaseIndex(cache_block)];

  vid[0] = faces[fid] * 3;
  vid[1] = faces[fid + 1] * 3;
  vid[2] = faces[fid + 2] * 3;

  const float* normals = &normals_[normalBaseIndex(cache_block)];

  // vertex 0
  normals_out[0] = normals[vid[0]];
  normals_out[1] = normals[vid[0] + 1];
  normals_out[2] = normals[vid[0] + 2];

  // vertex 1
  normals_out[3] = normals[vid[1]];
  normals_out[4] = normals[vid[1] + 1];
  normals_out[5] = normals[vid[1] + 2];

  // vertex 2
  normals_out[6] = normals[vid[2]];
  normals_out[7] = normals[vid[2] + 1];
  normals_out[8] = normals[vid[2] + 2];
}

void HybridGeometryBuffer::computeNormals(int cache_block) {
  const float* vertices = &vertices_[vertexBaseIndex(cache_block)];
  const uint32_t* faces = &faces_[faceBaseIndex(cache_block)];

  float* normals = &normals_[normalBaseIndex(cache_block)];
  std::size_t normals_size = num_vertices_[cache_block] * 3;

  for (std::size_t i = 0; i < normals_size; i += 3) {
    normals[i] = 0.0f;
    normals[i + 1] = 0.0f;
    normals[i + 2] = 0.0f;
  }

  glm::vec3 v0, v1, v2;  // vertices
  glm::vec3 n, u, v, dir;

  std::size_t fid;
  std::size_t vid[3];

  std::size_t num_faces = num_faces_[cache_block];

  for (std::size_t i = 0; i < num_faces; ++i) {
    fid = i * NUM_VERTICES_PER_FACE;

    vid[0] = faces[fid] * 3;
    vid[1] = faces[fid + 1] * 3;
    vid[2] = faces[fid + 2] * 3;

    v0.x = vertices[vid[0]];
    v0.y = vertices[vid[0] + 1];
    v0.z = vertices[vid[0] + 2];

    v1.x = vertices[vid[1]];
    v1.y = vertices[vid[1] + 1];
    v1.z = vertices[vid[1] + 2];

    v2.x = vertices[vid[2]];
    v2.y = vertices[vid[2] + 1];
    v2.z = vertices[vid[2] + 2];

    u = v1 - v0;
    v = v2 - v0;

    n.x = u.y * v.z - u.z * v.y;
    n.y = u.z * v.x - u.x * v.z;
    n.z = u.x * v.y - u.y * v.x;

    normals[vid[0]] += n.x;
    normals[vid[0] + 1] += n.y;
    normals[vid[0] + 2] += n.z;

    normals[vid[1]] += n.x;
    normals[vid[1] + 1] += n.y;
    normals[vid[1] + 2] += n.z;

    normals[vid[2]] += n.x;
    normals[vid[2] + 1] += n.y;
    normals[vid[2] + 2] += n.z;
  }
}

void HybridGeometryBuffer::updateTriangleIntersection(
    int cache_block, RTCRayIntersection* isect) const {
  // cache_block pointing to the current cache block in the mesh buffer
  // isect->primID, the current primitive intersected
  // colors: per-vertex colors
  uint32_t colors[3];
  getColorTuple(cache_block, isect->primID, colors);

  // interploate color tuple and update isect->color
  float u = isect->u;
  float v = isect->v;

  uint32_t rgb[9];
  util::unpack(colors[0], &rgb[0]);
  util::unpack(colors[1], &rgb[3]);
  util::unpack(colors[2], &rgb[6]);

  float w = 1.f - u - v;
  uint32_t r = (rgb[0] * w) + (rgb[3] * u) + (rgb[6] * v);
  uint32_t g = (rgb[1] * w) + (rgb[4] * u) + (rgb[7] * v);
  uint32_t b = (rgb[2] * w) + (rgb[5] * u) + (rgb[8] * v);

  isect->color = util::pack(r, g, b);

  // shading normal

  float ns[9];
  getNormalTuple(cache_block, isect->primID, ns);

  isect->Ns[0] = (ns[0] * w) + (ns[3] * u) + (ns[6] * v);
  isect->Ns[1] = (ns[1] * w) + (ns[4] * u) + (ns[7] * v);
  isect->Ns[2] = (ns[2] * w) + (ns[5] * u) + (ns[8] * v);
}

}  // namespace spray

