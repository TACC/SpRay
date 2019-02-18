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

#include <string>

#include "glm/glm.hpp"

#include "io/ply_loader.h"
#include "render/aabb.h"
#include "render/material.h"
#include "render/reflection.h"
#include "render/shape.h"
#include "utils/util.h"

namespace spray {

class SurfaceModel {
 public:
  SurfaceModel()
      : material_(nullptr),
        num_vertices_(0),
        num_faces_(0),
        transform_(glm::mat4(1.0f)),
        bitfields_(0) {}
  ~SurfaceModel() { delete material_; }

  void populateModelInfo();
  bool isValid() const { return (num_vertices_ > 0); }

  const std::string& getFilename() const { return filename_; }
  std::size_t getNumVertices() const { return num_vertices_; }
  std::size_t getNumFaces() const { return num_faces_; }
  const glm::mat4& getTransform() const { return transform_; }
  const Aabb& getObjectAabb() const { return object_aabb_; }
  const Material* getMaterial() const { return material_; }
  bool hasMaterial() const { return material_; }
  bool hasFile() const { return !filename_.empty(); }

  bool isNumVerticesSet() const {
    return (BIT_MASK_NUM_VERTICES_FLAG & bitfields_) ==
           BIT_MASK_NUM_VERTICES_FLAG;
  }

  bool isNumFacesSet() const {
    return (BIT_MASK_NUM_FACES_FLAG & bitfields_) == BIT_MASK_NUM_FACES_FLAG;
  }

  bool isObjectAabbSet() const {
    return (BIT_MASK_OBJECT_AABB_FLAG & bitfields_) == BIT_MASK_OBJECT_AABB_FLAG;
  }

  bool isConfigured() const {
    return isNumVerticesSet() && isNumFacesSet() && isObjectAabbSet();
  }

 private:
  friend class SceneLoader;
  friend class Domain;
  void setMaterial(Material* m) { material_ = m; }
  void setFilename(const std::string& filename) { filename_ = filename; }
  void setTransform(const glm::mat4& transform) { transform_ = transform; }

  // vertices

  void setNumVertices(std::size_t n) {
    CHECK_EQ(isNumVerticesSet(), false);
    setNumVerticesFlag();
    num_vertices_ = n;
  }

  void setNumVerticesFlag() { bitfields_ |= BIT_MASK_NUM_VERTICES_FLAG; }

  // faces

  void setNumFaces(std::size_t n) {
    CHECK_EQ(isNumFacesSet(), false);
    setNumFacesFlag();
    num_faces_ = n;
  }

  void setNumFacesFlag() { bitfields_ |= BIT_MASK_NUM_FACES_FLAG; }

  // object aabb

  void setObjectAabb(const float aabb_min[3], const float aabb_max[3]) {
    CHECK_EQ(isObjectAabbSet(), false);
    setObjectAabbFlag();
    object_aabb_.setBounds(aabb_min, aabb_max);
  }

  void setObjectAabbFlag() { bitfields_ |= BIT_MASK_OBJECT_AABB_FLAG; }

 private:
  std::string filename_;

  Material* material_;
  glm::mat4 transform_;
  Aabb object_aabb_;

  std::size_t num_vertices_;
  std::size_t num_faces_;

  /**
   * [31:2] : reserved
   * [2]    : object aabb set
   * [1]    : num_vertices_ set
   * [0]    : num_faces_ set
   */
  uint32_t bitfields_;

  enum BitMask {
    BIT_MASK_OBJECT_AABB_FLAG = 0x00000100,
    BIT_MASK_NUM_VERTICES_FLAG = 0x00000010,
    BIT_MASK_NUM_FACES_FLAG = 0x00000001
  };

  enum BitPosition {
    BIT_POS_OBJECT_AABB_FLAG = 2,
    BIT_POS_NUM_VERTICES_FLAG = 1,
    BIT_POS_NUM_FACES_FLAG = 0
  };
};

inline void SurfaceModel::populateModelInfo() {
  CHECK_EQ(filename_.empty(), false);
  if (isConfigured()) return;

  std::string ext = spray::util::getFileExtension(filename_);

  if (ext == "ply") {
    PlyLoader::LongHeader header;
    PlyLoader::readLongHeader(filename_, &header);

    num_vertices_ = header.num_vertices;
    num_faces_ = header.num_faces;
    object_aabb_ = header.bounds;

    CHECK(object_aabb_.isValid()) << object_aabb_;

  } else {
    CHECK(false) << "unsupported extension" << ext;
  }
}

class Domain {
 public:
  Domain() : num_vertices_(0), num_faces_(0) {}

  ~Domain() {
    for (std::size_t i = 0; i < shapes_.size(); ++i) {
      delete shapes_[i];
    }
  }

  void updateDomainInfo();

  unsigned int getId() const { return id_; }
  std::size_t getNumVertices() const { return num_vertices_; }
  std::size_t getNumVerticesPrefixSum(std::size_t model_id) const {
    return num_vertices_prefix_sums_[model_id];
  }
  std::size_t getNumFacesPrefixSum(std::size_t model_id) const {
    return num_faces_prefix_sums_[model_id];
  }
  std::size_t getNumFaces() const { return num_faces_; }

  const std::vector<SurfaceModel>& getModels() const { return models_; }
  const Aabb& getWorldAabb() const { return world_aabb_; }
  const std::vector<Shape*>& getShapes() const { return shapes_; }

  const Material* getMaterial(std::size_t model_id) const {
    return models_[model_id].getMaterial();
  }

  bool hasModels() const { return !models_.empty(); }
  std::size_t getNumModels() const { return models_.size(); }

  bool hasShapes() const { return !shapes_.empty(); }
  std::size_t getNumShapes() const { return shapes_.size(); }

  const std::string& getFilename(std::size_t i) const {
    return models_[i].getFilename();
  }

  void setNumVertices(std::size_t model_id, std::size_t num_vertices) {
    models_[model_id].setNumVertices(num_vertices);
  }

  void setNumFaces(std::size_t model_id, std::size_t num_faces) {
    models_[model_id].setNumFaces(num_faces);
  }

  void setObjectAabb(std::size_t model_id, const float aabb_min[3],
                     const float aabb_max[3]) {
    models_[model_id].setObjectAabb(aabb_min, aabb_max);
  }

  const SurfaceModel& getModel(std::size_t i) const { return models_[i]; }

  bool isConfigured() const {
    bool domain_configured =
        isNumVerticesSet() && isNumFacesSet() && isWorldAabbSet();
    std::size_t size = 0;
    for (std::size_t i = 0; i < models_.size(); ++i) {
      size += models_[i].isConfigured();
    }
    bool models_configured = (size == models_.size());
    return domain_configured && models_configured;
  }

  // void load(float* vertices, uint32_t* faces, uint32_t* colors);
 private:
  friend class SceneLoader;

  SurfaceModel* getModelPtr(std::size_t i) { return &models_[i]; }

  void setId(int id) { id_ = id; }
  void push(Shape* shape) { shapes_.push_back(shape); }
  void resizeModels(std::size_t n) {
    models_.resize(n);
    num_vertices_prefix_sums_.resize(n);
    num_faces_prefix_sums_.resize(n);
  }

  // world aabb

  void setWorldAabb(const glm::vec3& min, const glm::vec3& max) {
    CHECK_EQ(isWorldAabbSet(), false);
    setWorldAabbFlag();
    world_aabb_.reset(min, max);
  }

  void setWorldAabbFlag() { bitfields_ |= BIT_MASK_WORLD_AABB_FLAG; }

  bool isWorldAabbSet() const {
    return (BIT_MASK_WORLD_AABB_FLAG & bitfields_) == BIT_MASK_WORLD_AABB_FLAG;
  }

  // vertices

  void setNumVertices(std::size_t n) {
    CHECK_EQ(isNumVerticesSet(), false);
    setNumVerticesFlag();
    num_vertices_ = n;
  }

  void setNumVerticesFlag() { bitfields_ |= BIT_MASK_NUM_VERTICES_FLAG; }

  bool isNumVerticesSet() const {
    return (BIT_MASK_NUM_VERTICES_FLAG & bitfields_) ==
           BIT_MASK_NUM_VERTICES_FLAG;
  }

  // faces

  void setNumFaces(std::size_t n) {
    CHECK_EQ(isNumFacesSet(), false);
    setNumFacesFlag();
    num_faces_ = n;
  }

  void setNumFacesFlag() { bitfields_ |= BIT_MASK_NUM_FACES_FLAG; }

  bool isNumFacesSet() const {
    return (BIT_MASK_NUM_FACES_FLAG & bitfields_) == BIT_MASK_NUM_FACES_FLAG;
  }

 private:
  std::vector<SurfaceModel> models_;
  std::vector<std::size_t> num_vertices_prefix_sums_;
  std::vector<std::size_t> num_faces_prefix_sums_;
  Aabb world_aabb_;
  std::vector<Shape*> shapes_;

  std::size_t num_vertices_;
  std::size_t num_faces_;
  int id_;
  
  /**
   * [31:2] : reserved
   * [2]    : world aabb set
   * [1]    : num_vertices_ set
   * [0]    : num_faces_ set
   */
  uint32_t bitfields_;

  enum BitMask {
    BIT_MASK_WORLD_AABB_FLAG = 0x00000100,
    BIT_MASK_NUM_VERTICES_FLAG = 0x00000010,
    BIT_MASK_NUM_FACES_FLAG = 0x00000001
  };

  enum BitPosition {
    BIT_POS_WORLD_AABB_FLAG = 2,
    BIT_POS_NUM_VERTICES_FLAG = 1,
    BIT_POS_NUM_FACES_FLAG = 0
  };
};

inline void Domain::updateDomainInfo() {
  glm::vec3 bounds_min, bounds_max;

  std::size_t num_vertices = 0;
  std::size_t num_faces = 0;

  Aabb temp_aabb, world_aabb;

  for (std::size_t i = 0; i < models_.size(); ++i) {
    num_vertices_prefix_sums_[i] = num_vertices;
    num_faces_prefix_sums_[i] = num_faces;

    SurfaceModel& m = models_[i];

    m.populateModelInfo();

    CHECK(m.isObjectAabbSet());
    bounds_min = m.getObjectAabb().bounds[0];
    bounds_max = m.getObjectAabb().bounds[1];

    temp_aabb.bounds[0] = m.getTransform() * glm::vec4(bounds_min, 1.0f);
    temp_aabb.bounds[1] = m.getTransform() * glm::vec4(bounds_max, 1.0f);
#ifdef PRINT_DOMAIN_BOUNDS
    std::cout << "object: " << m.getObjectAabb() << "\n";
    std::cout << "world: " << temp_aabb << "\n";
#endif
    world_aabb.merge(temp_aabb);

    CHECK(m.isNumVerticesSet());
    num_vertices += m.getNumVertices();

    CHECK(m.isNumFacesSet());
    num_faces += m.getNumFaces();
  }

  for (std::size_t i = 0; i < shapes_.size(); ++i) {
    Shape* shape = shapes_[i];
    shape->getBounds(&temp_aabb);
    world_aabb.merge(temp_aabb);
  }

  if (!isNumVerticesSet()) {
    setNumVertices(num_vertices);
  } else {
    CHECK_EQ(num_vertices, getNumVertices());
  }

  if (!isNumFacesSet()) {
    setNumFaces(num_faces);
  } else {
    CHECK_EQ(num_faces, getNumFaces());
  }

  if (!isWorldAabbSet()) {
    world_aabb_ = world_aabb;
  } else {
    CHECK_EQ(world_aabb.within(world_aabb_), true);
  }

#ifdef PRINT_DOMAIN_BOUNDS
  std::cout << "scene bound (world): " << world_aabb_ << "\n";
#endif
}

}  // namespace spray

