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

#include "render/domain.h"
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
      device_(nullptr),
      scenes_(nullptr),
      embree_mesh_created_(nullptr),
      shape_created_(nullptr),
      shape_geom_ids_(nullptr) {}

HybridGeometryBuffer::~HybridGeometryBuffer() { cleanup(); }

void HybridGeometryBuffer::init(bool use_spray_color,
                                int max_cache_size_ndomains,
                                std::size_t max_nvertices,
                                std::size_t max_nfaces) {
  // cleanup
  cleanup();

  // sizes
  max_cache_size_ = max_cache_size_ndomains;
  max_nvertices_ = max_nvertices;
  max_nfaces_ = max_nfaces;

  std::size_t cache_size = static_cast<std::size_t>(max_cache_size_ndomains);
  CHECK_GT(cache_size, 0);

  // vertices
  vertex_buffer_size_ = 3 * max_nvertices * cache_size;
  if (vertex_buffer_size_) {
    vertices_ = arena_.Alloc<float>(vertex_buffer_size_ * 3, false);
    CHECK_NOTNULL(vertices_);
  }

  // per-vertex normals
  if (vertex_buffer_size_) {
    normals_ = arena_.Alloc<float>(vertex_buffer_size_, false);
    CHECK_NOTNULL(normals_);
  }

  // faces
  face_buffer_size_ = 3 * max_nfaces * cache_size;
  if (face_buffer_size_) {
    faces_ = arena_.Alloc<uint32_t>(face_buffer_size_, false);
    CHECK_NOTNULL(faces_);
  }

  // colors
  if (use_spray_color) {
    color_buffer_size_ = 0;
    colors_ = nullptr;
  } else {
    color_buffer_size_ = max_nvertices * cache_size;
    if (color_buffer_size_) {
      colors_ = arena_.Alloc<uint32_t>(color_buffer_size_, false);
      CHECK_NOTNULL(colors_);
    }
  }

  // domains
  domains_.resize(cache_size);

  // embree mesh created
  if (vertex_buffer_size_ > 0 && face_buffer_size_ > 0) {
    embree_mesh_created_ = arena_.Alloc<int>(cache_size);
    CHECK_NOTNULL(embree_mesh_created_);
    for (std::size_t i = 0; i < cache_size; ++i) {
      embree_mesh_created_[i] = DESTROYED;
    }
  }

  // shape created
  shape_created_ = arena_.Alloc<int>(cache_size);
  CHECK_NOTNULL(shape_created_);

  shape_geom_ids_ = arena_.Alloc<unsigned int>(cache_size);
  CHECK_NOTNULL(shape_geom_ids_);

  for (std::size_t i = 0; i < cache_size; ++i) {
    shape_created_[i] = DESTROYED;
  }

  // embree device
  device_ = rtcNewDevice("tri_accel=bvh4.triangle4v,threads=1");
  CHECK_NOTNULL(device_);

  // embree scenes
  scenes_ = arena_.Alloc<RTCScene>(cache_size, false);
  CHECK_NOTNULL(scenes_);

  // RTCSceneFlags sflags = RTC_SCENE_STATIC; // dont use this
  RTCSceneFlags sflags = RTC_SCENE_DYNAMIC;
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

RTCScene HybridGeometryBuffer::load(int cache_block, Domain& domain) {
  auto& shapes = domain.getShapes();

  CHECK(domain.hasModels() || !shapes.empty());

  unsigned int shape_geom_id = 0;

  if (domain.hasModels()) {
    shape_geom_id = domain.getNumModels();
    loadTriangles(cache_block, domain);
  }

  // std::cout << "cache_block " << cache_block << " [shape_geom_id "
  //           << shape_geom_id << "\n";
  shape_geom_ids_[cache_block] = shape_geom_id;

  if (!shapes.empty()) {
    loadShapes(shapes, cache_block, shape_geom_id);
  }

  RTCScene scene = scenes_[cache_block];
  rtcCommit(scene);

  return scene;
}

void HybridGeometryBuffer::loadTriangles(int cache_block,
                                         const Domain& domain) {
  // set domain
  domains_[cache_block] = &domain;

  // setup
  PlyLoader::Data d;

  const auto& models = domain.getModels();

  std::size_t sum_num_vertices = 0;
  std::size_t sum_num_faces = 0;

  float* vertices_base = &vertices_[vertexBaseIndex(cache_block)];
  uint32_t* faces_base = &faces_[faceBaseIndex(cache_block)];
  uint32_t* colors_base = &colors_[colorBaseIndex(cache_block)];

  // load surface data
  glm::vec4 v;
  for (std::size_t i = 0; i < models.size(); ++i) {
    const auto& model = models[i];

    d.vertices_capacity = model.getNumVertices() * 3;
    d.faces_capacity = model.getNumFaces() * NUM_VERTICES_PER_FACE;
    d.colors_capacity = model.getNumVertices();

    d.vertices = &vertices_base[sum_num_vertices * 3];
    d.faces = &faces_base[sum_num_faces * 3];
    d.colors = &colors_base[sum_num_vertices];

    loader_.load(model.getFilename(), &d);

    // TODO: transform vertices in the loader
    bool apply_transform = (model.getTransform() != glm::mat4(1.0));
    if (apply_transform) {
      float* vertices = &vertices_base[sum_num_vertices * 3];
      // uint32_t* faces = &faces_base[sum_num_faces * 3];

      const glm::mat4& x = model.getTransform();

      std::size_t nverts = model.getNumVertices() * 3;

      for (std::size_t n = 0; n < nverts; n += 3) {
        v = x * glm::vec4(vertices[n], vertices[n + 1], vertices[n + 2], 1.0f);
        vertices[n] = v.x;
        vertices[n + 1] = v.y;
        vertices[n + 2] = v.z;
      }
    }

    mapEmbreeBuffer(cache_block, d.vertices, model.getNumVertices(), d.faces,
                    model.getNumFaces(), models.size(), i);

    float* normals = &normals_[normalBaseIndex(cache_block, sum_num_vertices)];

    computeNormals(cache_block, d.vertices, model.getNumVertices(), d.faces,
                   model.getNumFaces(), normals);

    sum_num_vertices += model.getNumVertices();
    sum_num_faces += model.getNumFaces();

    CHECK_EQ(d.num_vertices, model.getNumVertices());
    CHECK_EQ(d.num_faces, model.getNumFaces());
  }

  CHECK_EQ(sum_num_vertices, domain.getNumVertices());
  CHECK_EQ(sum_num_faces, domain.getNumFaces());
}

void HybridGeometryBuffer::loadShapes(const std::vector<Shape*>& shapes,
                                      int cache_block,
                                      unsigned int shape_geom_id) {
  RTCScene scene = scenes_[cache_block];
  if (shape_created_[cache_block] == DESTROYED) {
    CHECK_LT(cache_block, shapes_.size());
    CHECK(shapes_[cache_block] == nullptr);
    shapes_[cache_block] = &shapes;

    // std::cout << "shapes.size: " << shapes.size()
    //           << " shapes pointer: " << &shapes << " cache block "
    //           << cache_block << std::endl;

    // std::cout<<"shapes.size: "<<shapes.size()<<std::endl;
    // create geometry
    unsigned int geom_id = rtcNewUserGeometry(scene, shapes.size());
    shape_created_[cache_block] = CREATED;

    CHECK_EQ(geom_id, shape_geom_id);
  }

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

void HybridGeometryBuffer::mapEmbreeBuffer(
    int cache_block, float* vertices, std::size_t num_vertices, uint32_t* faces,
    std::size_t num_faces, std::size_t num_models, std::size_t model_id) {
  // select scene
  RTCScene scene = scenes_[cache_block];

  // create triangle mesh
  if (embree_mesh_created_[cache_block] == DESTROYED) {
    unsigned int geom_id =
        rtcNewTriangleMesh(scene, RTC_GEOMETRY_STATIC, num_faces, num_vertices,
                           1 /*numTimeSteps*/);
#ifdef DEBUG_MESH
    std::cout << "created embree triangle mesh geom ID: " << geom_id
              << " cache block " << cache_block;
#endif
    CHECK_EQ(geom_id, model_id);
    CHECK_NE(geom_id, RTC_INVALID_GEOMETRY_ID);

    if (model_id == num_models - 1) {
      embree_mesh_created_[cache_block] = CREATED;
    }
  }

  // map vertices

  rtcSetBuffer2(scene, model_id /* geomID */, RTC_VERTEX_BUFFER, vertices, 0,
                sizeof(float) * 3, num_vertices);

  // map faces

  rtcSetBuffer2(scene, model_id /* geomID */, RTC_INDEX_BUFFER, faces, 0,
                sizeof(uint32_t) * NUM_VERTICES_PER_FACE, num_faces);

  rtcUpdate(scene, model_id /*geomID*/);
  rtcEnable(scene, model_id /*geomID*/);
}

void HybridGeometryBuffer::getColorTuple(const Domain& domain, int cache_block,
                                         uint32_t geomID, uint32_t primID,
                                         uint32_t colors[3]) const {
  uint32_t* faces = &faces_[faceBaseIndex(domain, cache_block, geomID)];

  uint32_t fid = primID * NUM_VERTICES_PER_FACE;
  uint32_t* c = &colors_[colorBaseIndex(domain, cache_block, geomID)];

  colors[0] = c[faces[fid]];
  colors[1] = c[faces[fid + 1]];
  colors[2] = c[faces[fid + 2]];
}

void HybridGeometryBuffer::getNormalTuple(const Domain& domain, int cache_block,
                                          uint32_t geomID, uint32_t primID,
                                          float normals_out[9]) const {
  std::size_t vid[3];
  const uint32_t fid = primID * NUM_VERTICES_PER_FACE;
  const uint32_t* faces = &faces_[faceBaseIndex(domain, cache_block, geomID)];

  vid[0] = faces[fid] * 3;
  vid[1] = faces[fid + 1] * 3;
  vid[2] = faces[fid + 2] * 3;

  const float* normals =
      &normals_[normalBaseIndex(domain, cache_block, geomID)];

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

void HybridGeometryBuffer::computeNormals(
    int cache_block, const float* vertices, std::size_t num_vertices,
    const uint32_t* faces, std::size_t num_faces, float* normals) {
  std::size_t normals_size = num_vertices * 3;

  for (std::size_t i = 0; i < normals_size; i += 3) {
    normals[i] = 0.0f;
    normals[i + 1] = 0.0f;
    normals[i + 2] = 0.0f;
  }

  glm::vec3 v0, v1, v2;  // vertices
  glm::vec3 n, u, v, dir;

  std::size_t fid;
  std::size_t vid[3];

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
  //
  const Domain& domain = *(domains_[cache_block]);

  // cache_block pointing to the current cache block in the mesh buffer
  // isect->primID, the current primitive intersected
  // colors: per-vertex colors

  unsigned int prim_id = isect->primID;
  unsigned int geom_id = isect->geomID;

  uint32_t colors[3];
  getColorTuple(domain, cache_block, geom_id, prim_id, colors);

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

  // material
  isect->material = getTriMeshMaterial(cache_block, geom_id);

  // shading normal

  float ns[9];
  getNormalTuple(domain, cache_block, geom_id, prim_id, ns);

  isect->Ns[0] = (ns[0] * w) + (ns[3] * u) + (ns[6] * v);
  isect->Ns[1] = (ns[1] * w) + (ns[4] * u) + (ns[7] * v);
  isect->Ns[2] = (ns[2] * w) + (ns[5] * u) + (ns[8] * v);
}

void HybridGeometryBuffer::updateTriangleIntersectionNoColor(
    int cache_block, RTCRayIntersection* isect) const {
  //
  const Domain& domain = *(domains_[cache_block]);

  // cache_block pointing to the current cache block in the mesh buffer
  // isect->primID, the current primitive intersected
  // colors: per-vertex colors

  unsigned int prim_id = isect->primID;
  unsigned int geom_id = isect->geomID;

  uint32_t colors[3];
  getColorTuple(domain, cache_block, geom_id, prim_id, colors);

  // interploate color tuple and update isect->color
  float u = isect->u;
  float v = isect->v;

  float w = 1.f - u - v;

  // material
  isect->material = getTriMeshMaterial(cache_block, geom_id);

  // shading normal

  float ns[9];
  getNormalTuple(domain, cache_block, geom_id, prim_id, ns);

  isect->Ns[0] = (ns[0] * w) + (ns[3] * u) + (ns[6] * v);
  isect->Ns[1] = (ns[1] * w) + (ns[4] * u) + (ns[7] * v);
  isect->Ns[2] = (ns[2] * w) + (ns[5] * u) + (ns[8] * v);
}

}  // namespace spray

