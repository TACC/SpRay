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
        transform_(glm::mat4(1.0f)) {}
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

 private:
  friend class SceneLoader;
  void setMaterial(Material* m) { material_ = m; }
  void setFilename(const std::string& filename) { filename_ = filename; }
  void setTransform(const glm::mat4& transform) { transform_ = transform; }
  void setNumVertices(std::size_t n) { num_vertices_ = n; }
  void setNumFaces(std::size_t n) { num_faces_ = n; }

 private:
  std::string filename_;
  std::size_t num_vertices_;
  std::size_t num_faces_;
  Material* material_;
  glm::mat4 transform_;
  Aabb object_aabb_;
};

inline void SurfaceModel::populateModelInfo() {
  if (!filename_.empty()) {
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
}

class Domain {
 public:
  Domain() : num_vertices_(0), num_faces_(0) {}

  ~Domain() {
    for (std::size_t i = 0; i < shapes_.size(); ++i) {
      delete shapes_[i];
    }
  }

  void populateModelInfo();

  unsigned int getId() const { return id_; }
  std::size_t getNumVertices() const { return num_vertices_; }
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

  // void load(float* vertices, uint32_t* faces, uint32_t* colors);
 private:
  friend class SceneLoader;

  SurfaceModel& getModel(std::size_t i) { return models_[i]; }
  void setId(int id) { id_ = id; }
  void push(Shape* shape) { shapes_.push_back(shape); }
  void setNumVertices(std::size_t n) { num_vertices_ = n; }
  void setNumFaces(std::size_t n) { num_faces_ = n; }
  void resizeModels(std::size_t n) { models_.resize(n); }

 private:
  int id_;
  std::size_t num_vertices_;
  std::size_t num_faces_;

  std::vector<SurfaceModel> models_;
  Aabb world_aabb_;
  std::vector<Shape*> shapes_;
};

inline void Domain::populateModelInfo() {
  glm::vec3 bounds_min, bounds_max, world_bounds_min, world_bounds_max;
  Aabb temp_aabb;
  num_vertices_ = 0;
  num_faces_ = 0;

  for (auto& m : models_) {
    if (!m.isValid()) {
      m.populateModelInfo();

      bounds_min = m.getObjectAabb().bounds[0];
      bounds_max = m.getObjectAabb().bounds[1];

      temp_aabb.bounds[0] = m.getTransform() * glm::vec4(bounds_min, 1.0f);
      temp_aabb.bounds[1] = m.getTransform() * glm::vec4(bounds_max, 1.0f);
#ifdef PRINT_DOMAIN_BOUNDS
      std::cout << "object: " << m.getObjectAabb() << "\n";
      std::cout << "world: " << temp_aabb << "\n";
#endif
      world_aabb_.merge(temp_aabb);
    }
    num_vertices_ += m.getNumVertices();
    num_faces_ += m.getNumFaces();
  }

  for (std::size_t i = 0; i < shapes_.size(); ++i) {
    Shape* shape = shapes_[i];
    shape->getBounds(&temp_aabb);
    world_aabb_.merge(temp_aabb);
  }

#ifdef PRINT_DOMAIN_BOUNDS
  std::cout << "scene bound (world): " << world_aabb_ << "\n";
#endif
}

}  // namespace spray

