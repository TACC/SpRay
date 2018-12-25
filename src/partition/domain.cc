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
#include "render/spray.h"
#include "scene/light.h"

// #define SPRAY_PRINT_LINES
// #define SPRAY_PRINT_TOKENS

namespace spray {

SceneParser::DomainTokenType SceneParser::getTokenType(const std::string& tag) {
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
  } else if (tag == "sphere") {
    type = DomainTokenType::kSphere;
  } else {
    LOG(FATAL) << "unknown tag name " << tag;
  }
  return type;
}

void SceneParser::parseDomain(const std::vector<std::string>& tokens) {
  nextDomain();
  Domain& d = currentDomain();
  d.id = domain_id_;
  d.transform = glm::mat4(1.f);
}

void SceneParser::parseFile(const std::string& ply_path,
                            const std::vector<std::string>& tokens) {
  Domain& d = currentDomain();

  CHECK_EQ(tokens.size(), 2);

  d.filename = ply_path.empty() ? tokens[1] : ply_path + "/" + tokens[1];
}

void SceneParser::parseMaterial(const std::vector<std::string>& tokens) {
  Domain& d = currentDomain();

  if (tokens[1] == "diffuse") {
    // mtl diffuse albedo<r g b>
    CHECK_EQ(tokens.size(), 5);
    glm::vec3 albedo;

    albedo[0] = atof(tokens[2].c_str());
    albedo[1] = atof(tokens[3].c_str());
    albedo[2] = atof(tokens[4].c_str());

    d.bsdf = new DiffuseBsdf(albedo);

  } else if (tokens[1] == "mirror") {
    // mtl mirror reflectance<r g b>
    CHECK_EQ(tokens.size(), 5);
    glm::vec3 reflectance;

    reflectance[0] = atof(tokens[2].c_str());
    reflectance[1] = atof(tokens[3].c_str());
    reflectance[2] = atof(tokens[4].c_str());

    d.bsdf = new MirrorBsdf(reflectance);

  } else if (tokens[1] == "glass") {
    // mtl mirror etaA etaB
    CHECK_EQ(tokens.size(), 4);

    float eta_a = atof(tokens[2].c_str());
    float eta_b = atof(tokens[3].c_str());

    d.bsdf = new GlassBsdf(eta_a, eta_b);

  } else if (tokens[1] == "transmission") {
    // mtl transmission etaA etaB
    CHECK_EQ(tokens.size(), 4);

    float eta_a = atof(tokens[2].c_str());
    float eta_b = atof(tokens[3].c_str());

    d.bsdf = new TransmissionBsdf(eta_a, eta_b);

  } else {
    LOG(FATAL) << "unknown material type " << tokens[1];
  }
}

void SceneParser::parseBound(const std::vector<std::string>& tokens) {
  Domain& d = currentDomain();

  CHECK_EQ(tokens.size(), 7);

  glm::vec3 min(atof(tokens[1].c_str()), atof(tokens[2].c_str()),
                atof(tokens[3].c_str()));
  glm::vec3 max(atof(tokens[4].c_str()), atof(tokens[5].c_str()),
                atof(tokens[6].c_str()));
  d.object_aabb.bounds[0] = min;
  d.object_aabb.bounds[1] = max;
}

void SceneParser::parseScale(const std::vector<std::string>& tokens) {
  Domain& d = currentDomain();

  CHECK_EQ(tokens.size(), 4);

  d.transform = glm::scale(
      d.transform, glm::vec3(atof(tokens[1].c_str()), atof(tokens[2].c_str()),
                             atof(tokens[3].c_str())));
}

void SceneParser::parseRotate(const std::vector<std::string>& tokens) {
  Domain& d = currentDomain();

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

void SceneParser::parseTranslate(const std::vector<std::string>& tokens) {
  Domain& d = currentDomain();

  CHECK_EQ(tokens.size(), 4);

  d.transform = glm::translate(
      d.transform, glm::vec3(atof(tokens[1].c_str()), atof(tokens[2].c_str()),
                             atof(tokens[3].c_str())));
}

void SceneParser::parseFace(const std::vector<std::string>& tokens) {
  Domain& d = currentDomain();

  CHECK_EQ(tokens.size(), 2);

  d.num_faces = std::stoul(tokens[1]);
}

void SceneParser::parseVertex(const std::vector<std::string>& tokens) {
  Domain& d = currentDomain();

  CHECK_EQ(tokens.size(), 2);

  d.num_vertices = std::stoul(tokens[1]);
}

void SceneParser::parseLight(const std::vector<std::string>& tokens) {
  if (tokens[1] == "point") {
    CHECK_EQ(tokens.size(), 8);
    glm::vec3 position, radiance;

    position[0] = atof(tokens[2].c_str());
    position[1] = atof(tokens[3].c_str());
    position[2] = atof(tokens[4].c_str());

    radiance[0] = atof(tokens[5].c_str());
    radiance[1] = atof(tokens[6].c_str());
    radiance[2] = atof(tokens[7].c_str());

    addLight(new PointLight(position, radiance));

  } else if (tokens[1] == "diffuse") {
    CHECK_EQ(tokens.size(), 5);
    glm::vec3 position, radiance;

    radiance[0] = atof(tokens[2].c_str());
    radiance[1] = atof(tokens[3].c_str());
    radiance[2] = atof(tokens[4].c_str());

    addLight(new DiffuseHemisphereLight(radiance));

  } else {
    LOG(FATAL) << "unknown light source " << tokens[1];
  }
}

// sphere <center_x> <center_y> <center_z> <radius>
void SceneParser::parseSphere(const std::vector<std::string>& tokens) {
  CHECK_EQ(tokens.size(), 5);
  Domain& d = currentDomain();

  glm::vec3 center;
  float radius;

  center[0] = atof(tokens[1].c_str());
  center[1] = atof(tokens[2].c_str());
  center[2] = atof(tokens[3].c_str());

  radius = atof(tokens[4].c_str());

  d.shapes.push_back(new Sphere(center, radius));

  // vertices, faces
  d.num_vertices = 0;
  d.num_faces = 0;
   
  // object bounds
  glm::vec3 v(radius);
  d.object_aabb.bounds[0] = center - v;
  d.object_aabb.bounds[1] = center + v;
}

void SceneParser::parseLineTokens(const std::string& ply_path,
                                  const std::vector<std::string>& tokens) {
  DomainTokenType type = getTokenType(tokens[0]);
  switch (type) {
    case DomainTokenType::kDomain:
      parseDomain(tokens);
      break;
    case DomainTokenType::kFile:
      parseFile(ply_path, tokens);
      break;
    case DomainTokenType::kMaterial:
      parseMaterial(tokens);
      break;
    case DomainTokenType::kBound:
      parseBound(tokens);
      break;
    case DomainTokenType::kScale:
      parseScale(tokens);
      break;
    case DomainTokenType::kRotate:
      parseRotate(tokens);
      break;
    case DomainTokenType::kTranslate:
      parseTranslate(tokens);
      break;
    case DomainTokenType::kFace:
      parseFace(tokens);
      break;
    case DomainTokenType::kVertex:
      parseVertex(tokens);
      break;
    case DomainTokenType::kLight:
      parseLight(tokens);
      break;
    case DomainTokenType::kSphere:
      parseSphere(tokens);
      break;
    default:
      break;
  }
}

void SceneParser::countAndAllocate(std::ifstream& infile) {
  int ndomains = 0;
  int nlights = 0;

  std::string line;

  while (infile.good()) {
    getline(infile, line);
    char* token = std::strtok(&line[0], " ");
    if (token != NULL) {
      if (strcmp(token, "domain") == 0) {
        ++ndomains;
      } else if (strcmp(token, "light") == 0) {
        ++nlights;
      }
    }
  }

  CHECK_GT(ndomains, 0);
  domains_->resize(ndomains);

  if (nlights)
    lights_->resize(nlights);
  else
    std::cout << "[WARNING] no lights detected\n";

  std::cout << "[INFO] number of domains: " << ndomains << "\n";
  std::cout << "[INFO] number of lights: " << nlights << "\n";
}

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

void SceneParser::parse(const std::string& filename,
                        const std::string& ply_path,
                        std::vector<Domain>* domains_out,
                        std::vector<Light*>* lights_out) {
  reset(domains_out, lights_out);

  std::ifstream infile(filename);
  CHECK(infile.is_open()) << "unable to open input file " << filename;

  countAndAllocate(infile);

  infile.clear();
  infile.seekg(infile.beg);

  char delim[] = " ";

  std::vector<std::string> tokens;
  while (infile.good()) {
    std::string line;
    getline(infile, line);
#ifdef SPRAY_PRINT_LINES
    std::cout << line << "\n";
#endif
    char* token = std::strtok(&line[0], delim);
    tokens.clear();
    while (token != NULL) {
#ifdef SPRAY_PRINT_TOKENS
      std::cout << token << " token len:" << std::string(token).size() << "\n";
#endif
      tokens.push_back(token);
      token = std::strtok(NULL, delim);
    }
    if (tokens.size()) parseLineTokens(ply_path, tokens);
  }

  infile.close();

  // apply transformation matrix
  for (auto& d : (*domains_out)) {
    glm::vec4 min = d.transform * glm::vec4(d.object_aabb.bounds[0], 1.0f);
    glm::vec4 max = d.transform * glm::vec4(d.object_aabb.bounds[1], 1.0f);
    d.world_aabb.bounds[0] = glm::vec3(min);
    d.world_aabb.bounds[1] = glm::vec3(max);
  }
}

}  // namespace spray
