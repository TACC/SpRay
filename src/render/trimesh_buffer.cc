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

#include "render/trimesh_buffer.h"

#include <embree2/rtcore_geometry.h>

#include "glog/logging.h"

#include "render/domain.h"
#include "render/rays.h"
#include "utils/util.h"

#define DEBUG_MESH
#undef DEBUG_MESH

namespace spray {

TriMeshBuffer::TriMeshBuffer()
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
      embree_mesh_created_(nullptr) {}

TriMeshBuffer::~TriMeshBuffer() { cleanup(); }

void TriMeshBuffer::init(const std::vector<Domain>& domains,
                         int max_cache_size_ndomains, std::size_t max_nvertices,
                         std::size_t max_nfaces) {
  // cleanup
  cleanup();

  // sizes
  max_cache_size_ = max_cache_size_ndomains;
  max_nvertices_ = max_nvertices;  // max num of vertices per domain
  max_nfaces_ = max_nfaces;        // max num of faces per domain
  
  std::size_t cache_size = static_cast<std::size_t>(max_cache_size_ndomains);

  max_nmodels_ = 0;  // max number of models per domain
  std::size_t model_count;

  for (const Domain& domain : domains) {
    const auto& models = domain.models;
    model_count = models.size();
    if (model_count > max_nmodels_) {
      max_nmodels_ = model_count;
    }
  }

  // resize prefix sums
  prefix_sum_num_vertices_ =
      arena_.Alloc<std::size_t>(cache_size * max_nmodels_, false);
  CHECK_NOTNULL(prefix_sum_num_vertices_);

  prefix_sum_num_faces_ =
      arena_.Alloc<std::size_t>(cache_size * max_nmodels_, false);
  CHECK_NOTNULL(prefix_sum_num_faces_);

  // number of vertices for each domain
  num_vertices_ = arena_.Alloc<std::size_t>(cache_size * max_nmodels_, false);
  CHECK_NOTNULL(num_vertices_);

  // vertices
  vertices_ = arena_.Alloc<float>(cache_size * max_nvertices_ * 3, false);
  CHECK_NOTNULL(vertices_);

  // per-vertex normals
  if (compute_normals) {
    normals_ = arena_.Alloc<float>(cache_size * max_nvertices_ * 3, false);
    CHECK_NOTNULL(normals_);
  } else {
    normals_ = nullptr;
  }

  // number of faces for each domain
  num_faces_ = arena_.Alloc<std::size_t>(cache_size * max_nmodels_, false);
  CHECK_NOTNULL(num_faces_);

  // faces
  faces_ = arena_.Alloc<uint32_t>(
      cache_size * max_nfaces_ * NUM_VERTICES_PER_FACE, false);
  CHECK_NOTNULL(faces_);

  // colors
  colors_ = arena_.Alloc<uint32_t>(cache_size * max_nvertices_, false);
  CHECK_NOTNULL(colors_);

  // materials
  materials_size_ = cache_size * max_nmodels_;
  materials_ = arena_.Alloc<Material*>(materials_size_, false);
  CHECK_NOTNULL(materials_);
  for (std::size_t i = 0; i < material_size_; ++i) {
    materials_[i] = nullptr;
  }

  // embree mesh created
  std::size_t num_meshes = cache_size * max_nmodels_;
  embree_mesh_created_ = arena_.Alloc<int>(num_meshes);
  CHECK_NOTNULL(embree_mesh_created_);

  for (std::size_t i = 0; i < num_meshes; ++i) {
    embree_mesh_created_[i] = DESTROYED;
  }

  // embree device
  device_ = rtcNewDevice("tri_accel=bvh4.triangle4v,threads=1");
  CHECK_NOTNULL(device_);

  // embree scenes
  scenes_ = arena_.Alloc<RTCScene>(cache_size, false);
  CHECK_NOTNULL(scenes_);

  for (std::size_t i = 0; i < cache_size; ++i) {
    scenes_[i] = rtcDeviceNewScene(device_, RTC_SCENE_DYNAMIC, RTC_INTERSECT1);
    CHECK_NOTNULL(scenes_[i]);
  }
}

/*
RTCScene TriMeshBuffer::load(const std::string& filename, int cache_block,
                             const glm::mat4& transform, bool apply_transform,
                             std::vector<Shape*>& shapes) {
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

  // return scene
  return scenes_[cache_block];
}
*/

RTCScene TriMeshBuffer::load(int cache_block, Domain& domain) {
  // init vertex/face counts and prefix sums
  std::size_t sum_num_vertices=0;
  std::size_t sum_num_faces=0;

  const auto& models = domain.models;

  for (std::size_t i = 0; i < models.size(); ++i) {
    setPrefixSumNumVertices(cache_block, i, sum_num_vertices);
    setPrefixSumNumFaces(cache_block, i, sum_num_faces);

    const auto& model = models[i];
    setNumVertices(cache_block, i, model.num_vertices);
    setNumFaces(cache_block, i, model.num_faces);

    sum_num_vertices += model.num_vertices;
    sum_num_faces += model.num_faces;
  }

  // TODO: support other model types
  // setup ply data
  PlyLoader::Data data;
  data.vertices_capacity = max_nvertices_ * 3;                // in
  data.faces_capacity = max_nfaces_ * NUM_VERTICES_PER_FACE;  // in
  data.colors_capacity = max_nvertices_;                      // in

  // std::size_t vertex_offset = vertexBaseIndex(cache_block);
  // std::size_t face_offset = faceBaseIndex(cache_block);
  // std::size_t color_offset = colorBaseIndex(cache_block);

  // float* vertices = &vertices_[vertex_offset];
  // uint32_t* faces = &faces_[face_offset];

  // load models
  for (std::size_t i = 0; i < models.size(); ++i) {
    const auto& model = domain.models[i];

    std::size_t vertex_offset = vertexBaseIndex(cache_block, i);
    std::size_t face_offset = faceBaseIndex(cache_block, i);
    std::size_t color_offset = colorBaseIndex(cache_block, i);

    data.vertices = &vertices_[vertex_offset];  // in/out
    // num_vertices;  // out
    data.faces = &faces_[face_offset];  // in/out
    // num_faces;  // out
    data.colors = &colors_[color_offset];  // rgb, in/out

    // load
    loader_.load(model.filename, &data);

    bool apply_transform = (model.transform != glm::mat4(1.0));

    if (apply_transform) {
      // glm::mat4 x = transform;
      glm::vec4 v;
      std::size_t nverts = model.num_vertices * 3;

      for (std::size_t n = 0; n < nverts; n += 3) {
        v = model.transform * glm::vec4(data.vertices[n], data.vertices[n + 1],
                                        data.vertices[n + 2], 1.0f);
        data.vertices[n] = v.x;
        data.vertices[n + 1] = v.y;
        data.vertices[n + 2] = v.z;
      }
    }
  }

  // update normals
  for (std::size_t i = 0; i < models.size(); ++i) {
    const auto& model = domain.models[i];
    computeNormals(cache_block, i, model);
  }

  // map buffers
  mapEmbreeBuffer(cache_block, vertices, domain.num_vertices, faces,
                  domain.num_faces);

  populateMaterials(cache_block, domain);

  // return scene
  return scenes_[cache_block];
}

void TriMeshBuffer::populateMaterials(int cache_block, Domain& domain) {
  Material** materials = &materials_[materialBaseIndex(cache_block)];
  const auto& models = domain.models;

  std::size_t n = 0;
  for (const auto& m : models) {
    m.
    materials[n] = new Matte();
  }
}

void TriMeshBuffer::cleanup() {
  // embree scenes
  for (int i = 0; i < max_cache_size_; ++i) {
    rtcDeleteScene(scenes_[i]);
  }

  // embree device
  rtcDeleteDevice(device_);

  for (std::size_t i = 0; i < total_num_models_; ++i) {
    delete materials_[i];
    materials_[i] = nullptr;
  }

  // arena
  arena_.Reset();

  // sizes
  max_cache_size_ = 0;
  max_nvertices_ = 0;
  max_nfaces_ = 0;
}

void TriMeshBuffer::mapEmbreeBuffer(int cache_block, float* vertices,
                                    std::size_t num_vertices, uint32_t* faces,
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
    CHECK_NE(geom_id, RTC_INVALID_GEOMETRY_ID);

    embree_mesh_created_[cache_block] = CREATED;
  }

  // map vertices

  rtcSetBuffer2(scene, 0 /*geomID*/, RTC_VERTEX_BUFFER, vertices, 0,
                sizeof(float) * 3, num_vertices);

  // map faces

  rtcSetBuffer2(scene, 0 /*geomID*/, RTC_INDEX_BUFFER, faces, 0,
                sizeof(uint32_t) * NUM_VERTICES_PER_FACE, num_faces);

  rtcUpdate(scene, 0 /*geomID*/);
  rtcEnable(scene, 0 /*geomID*/);

  rtcCommit(scene);
}

void TriMeshBuffer::getColorTuple(int cache_block, uint32_t primID,
                                  uint32_t colors[3]) const {
  uint32_t* faces = &faces_[faceBaseIndex(cache_block)];

  uint32_t fid = primID * NUM_VERTICES_PER_FACE;
  uint32_t* c = &colors_[colorBaseIndex(cache_block)];

  colors[0] = c[faces[fid]];
  colors[1] = c[faces[fid + 1]];
  colors[2] = c[faces[fid + 2]];
}

void TriMeshBuffer::getNormalTuple(int cache_block, uint32_t primID,
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

void TriMeshBuffer::computeNormals(int cache_block, int model_id,
                                   const ModelFile& model) {
  //
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

void TriMeshBuffer::updateIntersection(int cache_block,
                                       RTCRayIntersection* isect) const {
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

