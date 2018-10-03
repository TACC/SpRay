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

#include <glog/logging.h>
#include <cstdlib>
#include <string>
#include <vector>

#include "glm/glm.hpp"

#include "accelerators/wbvh_embree.h"
#include "caches/caches.h"
#include "display/opengl.h"
#include "partition/aabb.h"
#include "partition/data_partition.h"
#include "partition/domain.h"
#include "partition/trimesh_buffer.h"
#include "renderers/rays.h"
#include "renderers/spray.h"
#include "scene/light.h"
#include "utils/math.h"
#include "utils/util.h"
#include "utils/comm.h"

#define PRINT_DOMAIN_BOUNDS
#undef PRINT_DOMAIN_BOUNDS

namespace spray {

class Light;

struct Intersection {
  bool hit;
  float t;
};

struct SceneInfo {
  RTCScene rtc_scene;
  int cache_block;
};

template <class CacheT>
class Scene {
 public:
  Scene() {}
  ~Scene() {
    for (std::size_t i = 0; i < lights_.size(); ++i) {
      delete lights_[i];
    }
    for (std::size_t i = 0; i < domains_.size(); ++i) {
      delete domains_[i].bsdf;
    }
    if (!storage_basepath_.empty()) deleteAllDomainsFromLocalDisk();
  }

  void init(const std::string& desc_filename, const std::string& ply_path,
            const std::string& storage_basepath, int cache_size, int view_mode,
            bool insitu_mode, int num_virtual_ranks);

  const InsituPartition& getInsituPartition() const { return partition_; }
  bool insitu() const { return insitu_; }

  void buildWbvh();

  Aabb getBound() const {
    CHECK(bound_.isValid()) << "invalid scene bound";
    return bound_;
  }

  // drawing domains
  void drawDomains();
  void drawPartitions();

  int nextDomain(int increment) {
    int delta = getNumDomains() - glfw_domain_idx_;
    if (increment < delta) {
      glfw_domain_idx_ += increment;
    } else {
      glfw_domain_idx_ = increment - delta;
    }
#ifdef SPRAY_GLOG_CHECK
    CHECK_GE(glfw_domain_idx_, 0);
    CHECK_LT(glfw_domain_idx_, getNumDomains());
#endif
    return domains_[glfw_domain_idx_].id;
  }

  int prevDomain(int decrement) {
    int delta = glfw_domain_idx_;
    if (decrement < delta) {
      glfw_domain_idx_ -= decrement;
    } else {
      glfw_domain_idx_ = getNumDomains() - (decrement - delta);
    }
#ifdef SPRAY_GLOG_CHECK
    CHECK_GE(glfw_domain_idx_, 0);
    CHECK_LT(glfw_domain_idx_, getNumDomains());
#endif
    return domains_[glfw_domain_idx_].id;
  }

  void setPartitionNumber(int num) { partition_num_ = num; }

  int nextPartition() {
    ++partition_num_;
    if (partition_num_ == num_partitions_) partition_num_ = 0;
    return partition_num_;
  }
  int prevPartition() {
    --partition_num_;
    if (partition_num_ == -1) partition_num_ = num_partitions_ - 1;
    return partition_num_;
  }

  void printPartitionSizes() const {
    std::vector<int> domain_to_partition;
    domain_to_partition = partition_.computePartition(num_partitions_);

    std::vector<int> partition_to_vertices;
    std::vector<int> partition_to_faces;
    partition_to_vertices.resize(num_partitions_, 0);
    partition_to_faces.resize(num_partitions_, 0);

    for (std::size_t i = 0; i < domains_.size(); ++i) {
      const Domain& d = domains_[i];
      int p = domain_to_partition[d.id];
      partition_to_vertices[p] += d.num_vertices;
      partition_to_faces[p] += d.num_faces;
    }
    for (int i = 0; i < num_partitions_; ++i) {
      std::cout << "parti[" << i << "]: v " << partition_to_vertices[i] << " f "
                << partition_to_faces[i] << std::endl;
    }
  }

 public:
  void load(int id);
  void load(int id, SceneInfo* sinfo);

  bool intersect(RTCScene rtc_scene, int cache_block, const float org[3],
                 const float dir[3], RTCRayIntersection* isect);
  bool intersect(RTCScene rtc_scene, int cache_block, const glm::vec3& org,
                 const float dir[3], RTCRayIntersection* isect);

  bool intersect(const float org[3], const float dir[3],
                 RTCRayIntersection* isect);

  bool occluded(const glm::vec3& org, const glm::vec3& dir, RTCRay* ray);
  bool occluded(RTCScene rtc_scene, const glm::vec3& org, const glm::vec3& dir,
                RTCRay* ray);

  bool occluded(const float org[3], const float dir[3], RTCRay* ray);
  bool occluded(RTCScene rtc_scene, const float org[3], const float dir[3],
                RTCRay* ray);

  void intersectDomains(RTCRayExt& ray) { wbvh_.intersect(ray); }

  void intersectDomains8(const unsigned valid[8], RTCRayExt8& ray) {
    wbvh_.intersect8(valid, ray);
  }

 public:
  Bsdf* getBsdf(int id) { return domains_[id].bsdf; }

 public:
  std::size_t getNumDomains() const { return domains_.size(); }
  const std::vector<Domain>& getDomains() const { return domains_; }

  const std::vector<Light*>& getLights() const { return lights_; }
  std::size_t getNumLights() const { return lights_.size(); }

 public:
  void updateIntersection(RTCRayIntersection* isect) const;
  void updateIntersection(int cache_block, RTCRayIntersection* isect) const;

 private:
  void mergeDomainBounds(std::size_t* max_num_vertices,
                         std::size_t* max_num_faces);

  void copyAllDomainsToLocalDisk(const std::string& dest_path,
                                 bool insitu_mode);
  void deleteAllDomainsFromLocalDisk();

 private:
  Aabb bound_;  // bound of entire scene in world space
  std::vector<Domain> domains_;
  std::vector<Light*> lights_;

  std::string storage_basepath_;

  CacheT cache_;
  TriMeshBuffer trimesh_buf_;

  RTCScene scene_;   // current domain's scene
  int cache_block_;  // current cache block

  WbvhEmbree wbvh_;

  InsituPartition partition_;
  bool insitu_;

  int glfw_domain_idx_;
  int partition_num_;
  int num_partitions_;
};

}  // namespace spray

#define SPRAY_SCENE_INL_
#include "scene/scene.inl"
#undef SPRAY_SCENE_INL_

