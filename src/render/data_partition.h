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

#include <algorithm>
#include <list>
#include <vector>

#include "glog/logging.h"

#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"

#include "render/aabb.h"
#include "render/domain.h"
#include "render/morton.h"

#undef GROUP_CLOSE_DOMAINS
#define GROUP_CLOSE_DOMAINS

namespace spray {

class InsituPartition {
 public:
  struct MortonCode {
    int domain;
    uint32_t code;
  };

 public:
  int rank(int id) const {
#ifdef SPRAY_GLOG_CHECK
    CHECK_LT(id, codes_.size());
#endif
    int rank = domain_to_rank_[id];
    return rank;
  }

  const std::vector<MortonCode>& getCodes() const { return codes_; }

  void partition(int ndomains, const std::vector<Domain>& domains,
                 Aabb scene_aabb, int num_ranks) {
#ifdef SPRAY_GLOG_CHECK
    CHECK_EQ(domains.size(), ndomains);
    CHECK_GT(num_ranks, 0);
#endif
    ndomains_ = ndomains;
    codes_.resize(ndomains);
    domain_to_rank_.resize(ndomains);

    glm::vec3 diag = scene_aabb.bounds[1] - scene_aabb.bounds[0];
    glm::vec3 scale_factor =
        glm::vec3(1.0f / diag.x, 1.0f / diag.y, 1.0f / diag.z);

    glm::vec3 min = scene_aabb.bounds[0] * scale_factor;
    glm::vec3 max = scene_aabb.bounds[1] * scale_factor;
    glm::vec3 min_offset = glm::vec3(0.0f, 0.0f, 0.0f) - min;
    glm::vec3 max_offset = glm::vec3(1.0f, 1.0f, 1.0f) - max;

#ifdef SPRAY_GLOG_CHECK
    CHECK_LT(glm::abs(min_offset.x - max_offset.x), 0.0001f);
    CHECK_LT(glm::abs(min_offset.y - max_offset.y), 0.0001f);
    CHECK_LT(glm::abs(min_offset.z - max_offset.z), 0.0001f);
#endif

    glm::mat4 trans = glm::translate(glm::mat4(), min_offset);
    glm::mat4 trasform = trans * glm::scale(glm::mat4(), scale_factor);

    // compute
    MortonCode m;

    for (int i = 0; i < ndomains; ++i) {
      // comptue centroids
      glm::vec3 center = domains[i].getWorldAabb().getCenter();
      glm::vec3 c = glm::vec3(trasform * glm::vec4(center, 1.0f));

      // morton codes
#ifdef SPRAY_GLOG_CHECK
      CHECK_LE(c.x, 1.0f);
      CHECK_LE(c.y, 1.0f);
      CHECK_LE(c.z, 1.0f);
      CHECK_GE(c.x, 0.0f);
      CHECK_GE(c.y, 0.0f);
      CHECK_GE(c.z, 0.0f);
      CHECK_LT(domains[i].getId(), ndomains_);
#endif

      m.domain = domains[i].getId();
      m.code = Morton::compute(c.x, c.y, c.z);
      codes_[i] = m;
    }

    // sort morton codes
    std::sort(codes_.begin(), codes_.end(),
              [](const MortonCode& a, const MortonCode& b) -> bool {
                return a.code < b.code;
              });

    // partition across the cluster
#ifdef GROUP_CLOSE_DOMAINS
    int shares = ndomains / num_ranks;
    int rank = 0;
    int s = 0;

    for (std::size_t i = 0; i < codes_.size(); ++i) {
      int domid = codes_[i].domain;
#ifdef SPRAY_GLOG_CHECK
      CHECK_LT(domid, ndomains_);
      CHECK_LT(domid, domain_to_rank_.size());
#endif
      domain_to_rank_[domid] = rank;
      ++s;
      if (s == shares) {
        s = 0;
        ++rank;
        if (rank == num_ranks) rank = 0;
      }
    }
    mapDomains(num_ranks);

#else
#error unverified
    // distribute domains in round-robin order
    // trying to scatter close domains
    int rank = 0;
    for (std::size_t i = 0; i < codes_.size(); ++i) {
      int domid = codes_[i].domain;
#ifdef SPRAY_GLOG_CHECK
      CHECK_LT(domid, ndomains_);
      CHECK_LT(domid, domain_to_rank_.size());
#endif
      domain_to_rank_[domid] = rank;
      ++rank;
      if (rank == num_ranks) rank = 0;
    }
    mapDomains(num_ranks);
#endif
  }

 public:
  std::vector<int> computePartition(int num_partitions) const {
#ifdef GROUP_CLOSE_DOMAINS
    int shares = ndomains_ / num_partitions;
    int partition = 0;
    int s = 0;

    std::vector<int> partition_numbers;
    partition_numbers.resize(ndomains_);

    CHECK_EQ(ndomains_, codes_.size());

    for (std::size_t i = 0; i < codes_.size(); ++i) {
      int domid = codes_[i].domain;
#ifdef SPRAY_GLOG_CHECK
      CHECK_LT(domid, ndomains_);
      CHECK_LT(domid, domain_to_rank_.size());
#endif
      partition_numbers[domid] = partition;
      ++s;
      if (s == shares) {
        s = 0;
        ++partition;
        if (partition == num_partitions) partition = 0;
      }
    }

    std::vector<int> domain_counts;  // partition to domain count
    domain_counts.resize(num_partitions, 0);
    for (std::size_t id = 0; id < partition_numbers.size(); ++id) {
      int part_id = partition_numbers[id];
      ++domain_counts[part_id];
    }
    for (std::size_t i = 0; i < domain_counts.size(); ++i) {
      CHECK_GT(domain_counts[i], 0)
          << "partition " << i << " number of domains: " << ndomains_;
    }

    return std::move(partition_numbers);
#else
#error unverified
    std::vector<int> partition_numbers;
    partition_numbers.resize(ndomains_);
    // distribute domains in round-robin order
    // trying to scatter close domains
    int partition = 0;
    for (std::size_t i = 0; i < codes_.size(); ++i) {
      int domid = codes_[i].domain;
#ifdef SPRAY_GLOG_CHECK
      CHECK_LT(domid, ndomains_);
      CHECK_LT(domid, domain_to_rank_.size());
#endif
      partition_numbers[domid] = partition;
      ++partition;
      if (partition == num_partitions) partition = 0;
    }
    return std::move(partition_numbers);
#endif
  }

 public:
  const std::list<int>& getDomains(int rank) const {
    return rank_to_domains_[rank];
  }
  std::size_t getNumDomains(int rank) const {
    return rank_to_domains_[rank].size();
  }

  int getNumDomains() const { return ndomains_; }

 private:
  void mapDomains(int num_ranks) {
#ifdef SPRAY_GLOG_CHECK
    CHECK_EQ(domain_to_rank_.size(), ndomains_);
    LOG(INFO) << "<<<< MORTON CODES >>>>";
    for (auto& c : codes_) {
      LOG(INFO) << "domain " << c.domain << " code " << c.code;
    }
#endif
    rank_to_domains_.resize(num_ranks);
    for (std::size_t id = 0; id < domain_to_rank_.size(); ++id) {
      int rank = domain_to_rank_[id];
      rank_to_domains_[rank].push_back(id);
    }
  }

 private:
  std::vector<MortonCode> codes_;
  std::vector<int> domain_to_rank_;
  std::vector<std::list<int>> rank_to_domains_;
  int ndomains_;
};

}  // namespace spray
