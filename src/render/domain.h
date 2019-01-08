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

struct SurfaceModel {
  SurfaceModel()
      : material(nullptr),
        num_vertices(0),
        num_faces(0),
        transform(glm::mat4(1.0f)) {}
  ~SurfaceModel() { delete material; }

  void populateModelInfo();
  bool isValid() const { return (num_vertices > 0); }

  std::string filename;
  std::size_t num_vertices;
  std::size_t num_faces;
  Material* material;
  glm::mat4 transform;
  Aabb object_aabb;
};

struct Domain {
  Domain() : num_vertices(0), num_faces(0) {}

  ~Domain() {
    for (std::size_t i = 0; i < shapes.size(); ++i) delete shapes[i];
  }

  void populateModelInfo();

  unsigned int id;
  std::size_t num_vertices;
  std::size_t num_faces;

  std::vector<SurfaceModel> models;

  Aabb world_aabb;

  std::vector<Shape*> shapes;
};

inline void Domain::populateModelInfo() {
  glm::vec3 bounds_min, bounds_max, world_bounds_min, world_bounds_max;
  Aabb temp_aabb;
  std::size_t num_vertices = 0;
  std::size_t num_faces = 0;

  for (auto& m : models) {
    if (!m.isValid()) {
      m.populateModelInfo();

      bounds_min = m.object_aabb.bounds[0];
      bounds_max = m.object_aabb.bounds[1];

      temp_aabb.bounds[0] = m.transform * glm::vec4(bounds_min, 1.0f);
      temp_aabb.bounds[1] = m.transform * glm::vec4(bounds_max, 1.0f);

      world_aabb.merge(temp_aabb);
    }
    num_vertices += m.num_vertices;
    num_faces += m.num_faces;
  }

  for (std::size_t i = 0; i < shapes.size(); ++i) {
    Shape* shape = shapes[i];
    shape->getBounds(&temp_aabb);
    world_aabb.merge(temp_aabb);
  }
}

inline void SurfaceModel::populateModelInfo() {
  if (!filename.empty()) {
    std::string ext = spray::util::getFileExtension(filename);

    if (ext == "ply") {
      PlyLoader::LongHeader header;
      PlyLoader::readLongHeader(filename, &header);

      num_vertices = header.num_vertices;
      num_faces = header.num_faces;
      object_aabb = header.bounds;

    } else {
      CHECK(false) << "unsupported extension" << ext;
    }
  }
}

}  // namespace spray

