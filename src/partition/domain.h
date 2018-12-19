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
#include <vector>

#include "glm/glm.hpp"

#include "partition/aabb.h"

namespace spray {

class Bsdf;
class Light;

struct Domain {
  Domain() : num_vertices(0), num_faces(0), bsdf(nullptr) {}
  // ~Domain() { delete bsdf; }

  unsigned id;  // single domain id
  std::size_t num_vertices;
  std::size_t num_faces;
  std::string filename;
  Aabb object_aabb;
  Aabb world_aabb;
  glm::mat4 transform;
  Bsdf* bsdf;
};

class SceneParser {
 private:
  int domain_id_;
  int light_id_;
  std::vector<Domain>* domains_;
  std::vector<Light*>* lights_;

 public:
  void parse(const std::string& filename, const std::string& ply_path,
             std::vector<Domain>* domains_out, std::vector<Light*>* lights_out);

 private:
  enum class DomainTokenType {
    kComment,
    kDomain,
    kFile,
    kMaterial,
    kBound,
    kScale,
    kRotate,
    kTranslate,
    kFace,
    kVertex,
    kLight
  };

  void reset(std::vector<Domain>* domains, std::vector<Light*>* lights) {
    domain_id_ = -1;
    light_id_ = 0;

    CHECK_NOTNULL(domains);
    CHECK(domains->size() == 0);

    domains_ = domains;

    if (lights) CHECK(lights->size() == 0);
    lights_ = lights;
  }

  void countAndAllocate(std::ifstream& infile);

  DomainTokenType getTokenType(const std::string& tag);

  void parseDomain(const std::vector<std::string>& tokens);

  void parseFile(const std::string& ply_path,
                 const std::vector<std::string>& tokens);

  void parseMaterial(const std::vector<std::string>& tokens);

  void parseBound(const std::vector<std::string>& tokens);

  void parseScale(const std::vector<std::string>& tokens);

  void parseRotate(const std::vector<std::string>& tokens);

  void parseTranslate(const std::vector<std::string>& tokens);

  void parseFace(const std::vector<std::string>& tokens);

  void parseVertex(const std::vector<std::string>& tokens);

  void parseLight(const std::vector<std::string>& tokens);

  void parseLineTokens(const std::string& ply_path,
                       const std::vector<std::string>& tokens);

  Domain& currentDomain() {
    CHECK_GT(domain_id_, -1);
    CHECK_LT(domain_id_, domains_->size());
    return (*domains_)[domain_id_];
  }

  void addLight(Light* light) {
    CHECK_NOTNULL(lights_);
    CHECK_LT(light_id_, lights_->size());
    (*lights_)[light_id_] = light;
    ++light_id_;
  }

  void nextDomain() { ++domain_id_; }
};

}  // namespace spray

