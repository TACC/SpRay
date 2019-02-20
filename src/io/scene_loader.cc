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

#include "io/scene_loader.h"

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>

#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "glog/logging.h"

#include "render/aabb.h"
#include "render/light.h"
#include "render/material.h"
#include "render/reflection.h"
#include "render/spray.h"
#include "utils/comm.h"

// #define SPRAY_PRINT_LINES
// #define SPRAY_PRINT_TOKENS

namespace spray {

SceneLoader::DomainTokenType SceneLoader::getTokenType(const std::string& tag) {
  DomainTokenType type;
  if (tag[0] == '#') {
    type = DomainTokenType::kComment;
  } else if (tag == "DomainBegin") {
    type = DomainTokenType::kDomainBegin;
  } else if (tag == "DomainEnd") {
    type = DomainTokenType::kDomainEnd;
  } else if (tag == "ModelBegin") {
    type = DomainTokenType::kModelBegin;
  } else if (tag == "ModelEnd") {
    type = DomainTokenType::kModelEnd;
  } else if (tag == "file") {
    type = DomainTokenType::kFile;
  } else if (tag == "material") {
    type = DomainTokenType::kMaterial;
  } else if (tag == "DomainBounds") {
    type = DomainTokenType::kDomainBounds;
  } else if (tag == "ModelBounds") {
    type = DomainTokenType::kModelBounds;
  } else if (tag == "scale") {
    type = DomainTokenType::kScale;
  } else if (tag == "rotate") {
    type = DomainTokenType::kRotate;
  } else if (tag == "translate") {
    type = DomainTokenType::kTranslate;
  } else if (tag == "DomainFace") {
    type = DomainTokenType::kDomainFace;
  } else if (tag == "DomainVertex") {
    type = DomainTokenType::kDomainVertex;
  } else if (tag == "ModelFace") {
    type = DomainTokenType::kModelFace;
  } else if (tag == "ModelVertex") {
    type = DomainTokenType::kModelVertex;
  } else if (tag == "light") {
    type = DomainTokenType::kLight;
  } else if (tag == "sphere") {
    type = DomainTokenType::kSphere;
  } else {
    LOG(FATAL) << "unknown tag name " << tag;
  }
  return type;
}

void SceneLoader::parseDomainBegin() {
  CHECK_EQ(num_domain_begins_, 0);
  Domain& d = currentDomain();
  d.setId(domain_id_);
  ++num_domain_begins_;
  num_domain_ends_ = 0;
  // d.transform = glm::mat4(1.f);
}

void SceneLoader::parseDomainEnd() {
  CHECK_EQ(num_domain_ends_, 0);
  nextDomain();
  resetModelId();
  ++num_domain_ends_;
  num_domain_begins_ = 0;
}

void SceneLoader::parseModelBegin() { SurfaceModel* m = currentModel(); }

void SceneLoader::parseModelEnd() {
  SurfaceModel* m = currentModel();
  if (!m->hasMaterial()) {
    m->setMaterial(new Matte());
  }
  CHECK(m->hasFile());
  nextModel();
}

void SceneLoader::parseFile(const std::string& ply_path,
                            const std::vector<std::string>& tokens) {
  SurfaceModel* m = currentModel();
  std::size_t tokens_size = tokens.size();

  CHECK_GE(tokens_size, 2);

  m->setFilename(ply_path.empty() ? tokens[1] : ply_path + "/" + tokens[1]);
}

void SceneLoader::parseMaterial(const std::vector<std::string>& tokens) {
  std::size_t tokens_size = tokens.size();
  CHECK_GT(tokens_size, 1);

  SurfaceModel* model = currentModel();
  CHECK(!model->hasMaterial())
      << "found more than one material [domain " << getDomainId() << "][model "
      << getModelId() << "]";

  if (tokens[1] == "matte") {
    if (tokens_size > 2) {
      CHECK_EQ(tokens_size, 5);

      glm::vec3 albedo;
      albedo[0] = atof(tokens[2].c_str());
      albedo[1] = atof(tokens[3].c_str());
      albedo[2] = atof(tokens[4].c_str());

      model->setMaterial(new Matte(albedo));

    } else {
      model->setMaterial(new Matte());
    }

  } else if (tokens[1] == "metal") {
    if (tokens_size > 2) {
      CHECK_EQ(tokens_size, 6);

      glm::vec3 albedo;
      albedo[0] = atof(tokens[2].c_str());
      albedo[1] = atof(tokens[3].c_str());
      albedo[2] = atof(tokens[4].c_str());

      float fuzz = atof(tokens[5].c_str());

      model->setMaterial(new Metal(albedo, fuzz));

    } else {
      model->setMaterial(new Metal());
    }

  } else if (tokens[1] == "dielectric") {
    if (tokens_size > 2) {
      CHECK_EQ(tokens_size, 3);

      float index = atof(tokens[2].c_str());

      model->setMaterial(new Dielectric(index));

    } else {
      model->setMaterial(new Dielectric());
    }
  } else {
    CHECK(false) << "unsupported material: " << tokens[2];
  }
}

/*
void SceneLoader::parseUnusedMaterial(const std::vector<std::string>& tokens) {
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
*/

void SceneLoader::parseDomainWorldBounds(
    const std::vector<std::string>& tokens) {
  CHECK_EQ(tokens.size(), 7);
  Domain& d = currentDomain();

  glm::vec3 min(atof(tokens[1].c_str()), atof(tokens[2].c_str()),
                atof(tokens[3].c_str()));
  glm::vec3 max(atof(tokens[4].c_str()), atof(tokens[5].c_str()),
                atof(tokens[6].c_str()));
  d.setWorldAabb(min, max);
}

void SceneLoader::parseModelWorldBounds(
    const std::vector<std::string>& tokens) {
  CHECK_EQ(tokens.size(), 7);
  SurfaceModel* m = currentModel();

  glm::vec3 min(atof(tokens[1].c_str()), atof(tokens[2].c_str()),
                atof(tokens[3].c_str()));
  glm::vec3 max(atof(tokens[4].c_str()), atof(tokens[5].c_str()),
                atof(tokens[6].c_str()));
  m->setObjectAabb(&min[0], &max[0]);
}

void SceneLoader::parseScale(const std::vector<std::string>& tokens) {
  CHECK_EQ(tokens.size(), 4);
  SurfaceModel* m = currentModel();

  m->setTransform(
      glm::scale(m->getTransform(),
                 glm::vec3(atof(tokens[1].c_str()), atof(tokens[2].c_str()),
                           atof(tokens[3].c_str()))));
}

void SceneLoader::parseRotate(const std::vector<std::string>& tokens) {
  CHECK_EQ(tokens.size(), 3);
  SurfaceModel* m = currentModel();

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
  m->setTransform(glm::rotate(
      m->getTransform(), (float)glm::radians(atof(tokens[2].c_str())), axis));
}

void SceneLoader::parseTranslate(const std::vector<std::string>& tokens) {
  CHECK_EQ(tokens.size(), 4);
  SurfaceModel* m = currentModel();

  m->setTransform(
      glm::translate(m->getTransform(),
                     glm::vec3(atof(tokens[1].c_str()), atof(tokens[2].c_str()),
                               atof(tokens[3].c_str()))));
}

void SceneLoader::parseDomainFace(const std::vector<std::string>& tokens) {
  CHECK_EQ(tokens.size(), 2);
  Domain& d = currentDomain();
  d.setNumFaces(std::stoul(tokens[1]));
}

void SceneLoader::parseDomainVertex(const std::vector<std::string>& tokens) {
  CHECK_EQ(tokens.size(), 2);
  Domain& d = currentDomain();
  d.setNumVertices(std::stoul(tokens[1]));
}

void SceneLoader::parseModelFace(const std::vector<std::string>& tokens) {
  CHECK_EQ(tokens.size(), 2);
  SurfaceModel* m = currentModel();
  m->setNumFaces(std::stoul(tokens[1]));
}

void SceneLoader::parseModelVertex(const std::vector<std::string>& tokens) {
  CHECK_EQ(tokens.size(), 2);
  SurfaceModel* m = currentModel();
  m->setNumVertices(std::stoul(tokens[1]));
}

void SceneLoader::parseLight(const std::vector<std::string>& tokens) {
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

  } else if (tokens[1] == "diffuse-hemisphere") {
    CHECK_EQ(tokens.size(), 5);
    glm::vec3 position, radiance;

    radiance[0] = atof(tokens[2].c_str());
    radiance[1] = atof(tokens[3].c_str());
    radiance[2] = atof(tokens[4].c_str());

    addLight(new DiffuseHemisphereLight(radiance, num_light_samples_));

  } else if (tokens[1] == "diffuse-sphere") {
    CHECK_EQ(tokens.size(), 5);
    glm::vec3 position, radiance;

    radiance[0] = atof(tokens[2].c_str());
    radiance[1] = atof(tokens[3].c_str());
    radiance[2] = atof(tokens[4].c_str());

    addLight(new DiffuseSphereLight(radiance, num_light_samples_));

  } else {
    LOG(FATAL) << "unknown light source " << tokens[1];
  }
}

// sphere <center_x> <center_y> <center_z> <radius>
void SceneLoader::parseSphere(const std::vector<std::string>& tokens) {
  CHECK_GT(tokens.size(), 6);
  Domain& d = currentDomain();

  glm::vec3 center;
  float radius;

  center[0] = atof(tokens[1].c_str());
  center[1] = atof(tokens[2].c_str());
  center[2] = atof(tokens[3].c_str());

  radius = atof(tokens[4].c_str());

  Material* m = nullptr;

  if (tokens[5] == "matte") {
    CHECK_EQ(tokens.size(), 9);

    glm::vec3 albedo;
    albedo[0] = atof(tokens[6].c_str());
    albedo[1] = atof(tokens[7].c_str());
    albedo[2] = atof(tokens[8].c_str());

    m = new Matte(albedo);

  } else if (tokens[5] == "metal") {
    CHECK_EQ(tokens.size(), 10);

    glm::vec3 albedo;
    albedo[0] = atof(tokens[6].c_str());
    albedo[1] = atof(tokens[7].c_str());
    albedo[2] = atof(tokens[8].c_str());

    float fuzz = atof(tokens[9].c_str());

    m = new Metal(albedo, fuzz);

  } else if (tokens[5] == "dielectric") {
    CHECK_EQ(tokens.size(), 7);

    float index = atof(tokens[6].c_str());

    m = new Dielectric(index);
  } else {
    CHECK(false) << "unknown material " << tokens[5];
  }

  CHECK_NOTNULL(m);

  d.push(new Sphere(center, radius, m));

  // vertices, faces
  // should be set once by Domain::updateDomainInfo
  // d.setNumVertices(0);
  // d.setNumFaces(0);

  // object bounds
  // glm::vec3 v(radius);
  // d.object_aabb.bounds[0] = center - v;
  // d.object_aabb.bounds[1] = center + v;
}

void SceneLoader::parseLineTokens(const std::string& ply_path,
                                  const std::vector<std::string>& tokens) {
  DomainTokenType type = getTokenType(tokens[0]);
  switch (type) {
    case DomainTokenType::kComment:
      // do nothing
      break;
    case DomainTokenType::kDomainBegin:
      parseDomainBegin();
      break;
    case DomainTokenType::kDomainEnd:
      parseDomainEnd();
      break;
    case DomainTokenType::kModelBegin:
      parseModelBegin();
      break;
    case DomainTokenType::kModelEnd:
      parseModelEnd();
      break;
    case DomainTokenType::kFile:
      parseFile(ply_path, tokens);
      break;
    case DomainTokenType::kMaterial:
      parseMaterial(tokens);
      break;
    case DomainTokenType::kDomainBounds:
      parseDomainWorldBounds(tokens);
      break;
    case DomainTokenType::kModelBounds:
      parseModelWorldBounds(tokens);
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
    case DomainTokenType::kDomainFace:
      parseDomainFace(tokens);
      break;
    case DomainTokenType::kDomainVertex:
      parseDomainVertex(tokens);
      break;
    case DomainTokenType::kModelFace:
      parseModelFace(tokens);
      break;
    case DomainTokenType::kModelVertex:
      parseModelVertex(tokens);
      break;
    case DomainTokenType::kLight:
      parseLight(tokens);
      break;
    case DomainTokenType::kSphere:
      parseSphere(tokens);
      break;
    default:
      CHECK(false) "unsupported token: " << static_cast<int>(type);
      break;
  }
}

void SceneLoader::countAndAllocate(std::ifstream& infile) {
  int ndomains = 0;
  int nlights = 0;

  std::string line;
  std::queue<int> num_models_q;  // domain 0 to N-1

  int num_models_per_domain = 0;
  int num_files_per_model = 0;

  while (infile.good()) {
    getline(infile, line);
    char* token = std::strtok(&line[0], " ");
    if (token != NULL) {
      if (strcmp(token, "DomainBegin") == 0) {
        ++ndomains;
        num_models_per_domain = 0;

      } else if (strcmp(token, "DomainEnd") == 0) {
        num_models_q.push(num_models_per_domain);

      } else if (strcmp(token, "ModelBegin") == 0) {
        ++num_models_per_domain;
        num_files_per_model = 0;

      } else if (strcmp(token, "ModelEnd") == 0) {
        CHECK_EQ(num_files_per_model, 1);

      } else if (strcmp(token, "light") == 0) {
        ++nlights;

      } else if (strcmp(token, "file") == 0) {
        ++num_files_per_model;
      }
    }
  }

  CHECK_GT(ndomains, 0);
  domains_->resize(ndomains);

  CHECK_EQ(num_models_q.size(), ndomains);

  for (int i = 0; i < ndomains; ++i) {
    Domain& domain = domains_->at(i);
    domain.resizeModels(num_models_q.front());
    num_models_q.pop();
  }

  if (nlights)
    lights_->resize(nlights);
  else
    std::cout << "[WARNING] no lights detected\n";

  if (mpi::isRootProcess()) {
    std::cout << "[INFO] number of domains: " << ndomains << "\n";
    std::cout << "[INFO] number of lights: " << nlights << "\n";
  }
}

void SceneLoader::load(const std::string& filename, const std::string& ply_path,
                       int num_light_samples, std::vector<Domain>* domains_out,
                       std::vector<Light*>* lights_out) {
  reset(num_light_samples, domains_out, lights_out);

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
  // for (auto& d : (*domains_out)) {
  //   glm::vec4 min = d.transform * glm::vec4(d.object_aabb.bounds[0], 1.0f);
  //   glm::vec4 max = d.transform * glm::vec4(d.object_aabb.bounds[1], 1.0f);
  //   d.world_aabb.bounds[0] = glm::vec3(min);
  //   d.world_aabb.bounds[1] = glm::vec3(max);
  // }
}

}  // namespace spray
