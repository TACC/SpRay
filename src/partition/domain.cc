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

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>

#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "glog/logging.h"

#include "materials/reflection.h"
#include "partition/aabb.h"
#include "partition/domain.h"
#include "renderers/spray.h"
#include "scene/light.h"

// #define PRINT_LINES
// #define PRINT_TOKENS

namespace spray {

// # domain 0
// domain
// file CornellBox-Original.obj
// mtl 1
// bound -1 -1 -1 1 1 1
// scale 1 1 1
// rotate 0 0 0
// translate 0 0 0
//
// # domain 1
// # tbd

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

inline DomainTokenType getTokenType(const std::string& tag) {
  DomainTokenType type;
  if (tag[0] == '#') {
    type = DomainTokenType::kComment;
  } else if (tag == "domain") {
    type = DomainTokenType::kDomain;
  } else if (tag == "file") {
    type = DomainTokenType::kFile;
  } else if (tag == "mtl") {
    type = DomainTokenType::kMaterial;
  } else if (tag == "bound") {
    type = DomainTokenType::kBound;
  } else if (tag == "scale") {
    type = DomainTokenType::kScale;
  } else if (tag == "rotate") {
    type = DomainTokenType::kRotate;
  } else if (tag == "translate") {
    type = DomainTokenType::kTranslate;
  } else if (tag == "face") {
    type = DomainTokenType::kFace;
  } else if (tag == "vertex") {
    type = DomainTokenType::kVertex;
  } else if (tag == "light") {
    type = DomainTokenType::kLight;
  } else {
    LOG(FATAL) << "unknown tag name " << tag;
  }
  return type;
}

inline void parseDomain(const std::vector<std::string>& tokens,
                        std::vector<spray::Domain>& domain) {
  spray::Domain d;
  d.id = domain.size();
  d.transform = glm::mat4(1.f);
  domain.push_back(d);
}

inline void parseFile(const std::string& ply_path,
                      const std::vector<std::string>& tokens,
                      std::vector<spray::Domain>& domain) {
  CHECK(domain.size());
  spray::Domain& d = domain[domain.size() - 1];

  CHECK_EQ(tokens.size(), 2);

  d.filename = ply_path.empty() ? tokens[1] : ply_path + "/" + tokens[1];
}

inline void parseMaterial(const std::vector<std::string>& tokens,
                          std::vector<spray::Domain>& domain) {
  CHECK(domain.size());
  spray::Domain& d = domain[domain.size() - 1];

  if (tokens[1] == "diffuse") {
    // mtl diffuse albedo<r g b>
    CHECK(tokens.size() == 5);
    glm::vec3 albedo;

    albedo[0] = atof(tokens[2].c_str());
    albedo[1] = atof(tokens[3].c_str());
    albedo[2] = atof(tokens[4].c_str());

    d.bsdf = new spray::DiffuseBsdf(albedo);

  } else if (tokens[1] == "mirror") {
    // mtl mirror reflectance<r g b>
    CHECK(tokens.size() == 5);
    glm::vec3 reflectance;

    reflectance[0] = atof(tokens[2].c_str());
    reflectance[1] = atof(tokens[3].c_str());
    reflectance[2] = atof(tokens[4].c_str());

    d.bsdf = new spray::MirrorBsdf(reflectance);

  } else if (tokens[1] == "glass") {
    // mtl mirror etaA etaB
    CHECK(tokens.size() == 4);

    float eta_a = atof(tokens[2].c_str());
    float eta_b = atof(tokens[3].c_str());

    d.bsdf = new spray::GlassBsdf(eta_a, eta_b);

  } else if (tokens[1] == "transmission") {
    // mtl transmission etaA etaB
    CHECK(tokens.size() == 4);

    float eta_a = atof(tokens[2].c_str());
    float eta_b = atof(tokens[3].c_str());

    d.bsdf = new spray::TransmissionBsdf(eta_a, eta_b);

  } else {
    LOG(FATAL) << "unknown material type " << tokens[1];
  }
}

inline void parseBound(const std::vector<std::string>& tokens,
                       std::vector<spray::Domain>& domain) {
  CHECK(domain.size());
  spray::Domain& d = domain[domain.size() - 1];

  CHECK_EQ(tokens.size(), 7);

  glm::vec3 min(atof(tokens[1].c_str()), atof(tokens[2].c_str()),
                atof(tokens[3].c_str()));
  glm::vec3 max(atof(tokens[4].c_str()), atof(tokens[5].c_str()),
                atof(tokens[6].c_str()));
  d.object_aabb.bounds[0] = min;
  d.object_aabb.bounds[1] = max;
}

inline void parseScale(const std::vector<std::string>& tokens,
                       std::vector<spray::Domain>& domain) {
  CHECK(domain.size());
  spray::Domain& d = domain[domain.size() - 1];

  CHECK_EQ(tokens.size(), 4);

  d.transform = glm::scale(
      d.transform, glm::vec3(atof(tokens[1].c_str()), atof(tokens[2].c_str()),
                             atof(tokens[3].c_str())));
}

inline void parseRotate(const std::vector<std::string>& tokens,
                        std::vector<spray::Domain>& domain) {
  CHECK(domain.size());
  spray::Domain& d = domain[domain.size() - 1];

  CHECK_EQ(tokens.size(), 3);

  glm::vec3 axis;
  if (tokens[1] == "x") {
    axis = glm::vec3(1.f, 0.f, 0.f);
  } else if (tokens[1] == "y") {
    axis = glm::vec3(0.f, 1.f, 0.f);
  } else if (tokens[1] == "z") {
    axis = glm::vec3(0.f, 0.f, 1.f);
  } else {
    LOG(FATAL) << "invalid axis name " << tokens[1];
  }
  d.transform = glm::rotate(d.transform,
                            (float)glm::radians(atof(tokens[2].c_str())), axis);
}

inline void parseTranslate(const std::vector<std::string>& tokens,
                           std::vector<spray::Domain>& domain) {
  CHECK(domain.size());
  spray::Domain& d = domain[domain.size() - 1];

  CHECK(tokens.size() == 4);

  d.transform = glm::translate(
      d.transform, glm::vec3(atof(tokens[1].c_str()), atof(tokens[2].c_str()),
                             atof(tokens[3].c_str())));
}

inline void parseFace(const std::vector<std::string>& tokens,
                      std::vector<spray::Domain>& domain) {
  CHECK(domain.size());
  spray::Domain& d = domain[domain.size() - 1];

  CHECK(tokens.size() == 2);

  d.num_faces = std::stoul(tokens[1]);
}

inline void parseVertex(const std::vector<std::string>& tokens,
                        std::vector<spray::Domain>& domain) {
  CHECK(domain.size());
  spray::Domain& d = domain[domain.size() - 1];

  CHECK(tokens.size() == 2);

  d.num_vertices = std::stoul(tokens[1]);
}

inline void parseLight(const std::vector<std::string>& tokens,
                       std::vector<spray::Light*>& lights) {
  if (tokens[1] == "point") {
    CHECK(tokens.size() == 8);
    glm::vec3 position, radiance;

    position[0] = atof(tokens[2].c_str());
    position[1] = atof(tokens[3].c_str());
    position[2] = atof(tokens[4].c_str());

    radiance[0] = atof(tokens[5].c_str());
    radiance[1] = atof(tokens[6].c_str());
    radiance[2] = atof(tokens[7].c_str());

    spray::Light* light = new spray::PointLight(position, radiance);
    lights.push_back(light);

  } else if (tokens[1] == "diffuse") {
    CHECK(tokens.size() == 5);
    glm::vec3 position, radiance;

    radiance[0] = atof(tokens[2].c_str());
    radiance[1] = atof(tokens[3].c_str());
    radiance[2] = atof(tokens[4].c_str());

    spray::Light* light = new spray::DiffuseHemisphereLight(radiance);
    lights.push_back(light);

  } else {
    LOG(FATAL) << "unknown light source " << tokens[1];
  }
}

inline void parseLineTokens(const std::string& ply_path,
                            const std::vector<std::string>& tokens,
                            std::vector<spray::Domain>& domain,
                            std::vector<spray::Light*>& lights) {
  DomainTokenType type = getTokenType(tokens[0]);
  switch (type) {
    case DomainTokenType::kDomain:
      parseDomain(tokens, domain);
      break;
    case DomainTokenType::kFile:
      parseFile(ply_path, tokens, domain);
      break;
    case DomainTokenType::kMaterial:
      parseMaterial(tokens, domain);
      break;
    case DomainTokenType::kBound:
      parseBound(tokens, domain);
      break;
    case DomainTokenType::kScale:
      parseScale(tokens, domain);
      break;
    case DomainTokenType::kRotate:
      parseRotate(tokens, domain);
      break;
    case DomainTokenType::kTranslate:
      parseTranslate(tokens, domain);
      break;
    case DomainTokenType::kFace:
      parseFace(tokens, domain);
      break;
    case DomainTokenType::kVertex:
      parseVertex(tokens, domain);
      break;
    case DomainTokenType::kLight:
      parseLight(tokens, lights);
      break;
    default:
      break;
  }
}

void loadDescriptor(const std::string& filename, const std::string& ply_path,
                    std::vector<Domain>* domains_out,
                    std::vector<Light*>* lights_out) {
  std::vector<Domain> domains;
  std::vector<Light*> lights;

  std::ifstream infile(filename);
  CHECK(infile.is_open()) << "unable to open input file " << filename;

  char delim[] = " ";

  std::vector<std::string> tokens;
  while (infile.good()) {
    std::string line;
    getline(infile, line);
#ifdef PRINT_LINES
    std::cout << line << "\n";
#endif
    char* token = std::strtok(&line[0], delim);
    tokens.clear();
    while (token != NULL) {
#ifdef PRINT_TOKENS
      std::cout << token << " token len:" << std::string(token).size() << "\n";
#endif
      tokens.push_back(token);
      token = std::strtok(NULL, delim);
    }
    if (tokens.size()) {
      parseLineTokens(ply_path, tokens, domains, lights);
    }
  }

  infile.close();

  // apply transformation matrix
  for (auto& d : domains) {
    glm::vec4 min = d.transform * glm::vec4(d.object_aabb.bounds[0], 1.0f);
    glm::vec4 max = d.transform * glm::vec4(d.object_aabb.bounds[1], 1.0f);
    d.world_aabb.bounds[0] = glm::vec3(min);
    d.world_aabb.bounds[1] = glm::vec3(max);
  }

  *domains_out = std::move(domains);
  *lights_out = std::move(lights);
}

}  // namespace spray
