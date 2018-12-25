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

#include <mpi.h>
#include <cstdint>
#include <iostream>
#include <vector>

#include "insitu/insitu_tcontext.h"
#include "render/spray.h"
#include "utils/comm.h"

namespace spray {

class InsituPartition;

namespace insitu {

class WorkStats {
 private:
  struct ScatterEntry {
    int world_num_blocks_to_proc;
    int rank_num_blocks_to_proc;
  };

  // number of blocks to recv for each rank
  std::vector<int> reduce_buf_;            // per-rank
  std::vector<int> block_counters_;        // per-domain
  std::vector<ScatterEntry> scatter_buf_;  // per-rank

 public:
  WorkStats() : ndomains_(0), num_blocks_to_recv_(0) {}

  void resize(int nranks, int nthreads, int ndomains) {
    reduce_buf_.resize(nranks, 0);
    scatter_buf_.resize(nranks);

    ndomains_ = ndomains;
    // bug fix: seg faults occuring with the following conditional statement
    // when only 1 thread is launched
    // if (nthreads > 1) {
    block_counters_.resize(ndomains + 1, 0);  // +1 for cached block
    // }
  }

  void reset() {
    num_blocks_to_recv_ = 0;
    for (auto& r : reduce_buf_) r = 0;
  }

 private:
  void resetBlockCounters() {
    for (auto& r : block_counters_) r = 0;
  }

 private:
  int ndomains_;
  int num_blocks_to_recv_;

 public:
  // cluster level reduction
  void reduce();
  bool recvDone(int num_blocks_recved) const {
    return (num_blocks_recved == num_blocks_to_recv_);
  }

 public:
  bool allDone() const {
    return (scatter_buf_[mpi::rank()].world_num_blocks_to_proc == 0);
  }

 public:
  template <typename TContextT>
  void reduceRadianceThreadWorkStats(int this_rank,
                                     const spray::InsituPartition* partition,
                                     std::vector<TContextT>& tcontexts);

  template <typename TContextT>
  void reduceThreadWorkStats(int this_rank,
                             const spray::InsituPartition* partition,
                             std::vector<TContextT>& tcontexts);

  void addNumDomains(int dest_rank, int num_blocks) {
    reduce_buf_[dest_rank] += num_blocks;
  }

 private:
  void reduceRayBlocks(std::queue<int>* blockIds);
  void updateReduceBuffer(const spray::InsituPartition* partition);
};

template <typename TContextT>
void WorkStats::reduceRadianceThreadWorkStats(
    int this_rank, const spray::InsituPartition* partition,
    std::vector<TContextT>& tcontexts) {
  reset();

  // radiance blocks
  resetBlockCounters();

  for (auto& t : tcontexts) {
    reduceRayBlocks(t.getRadianceBlockIds());
  }

  updateReduceBuffer(partition);
}

template <typename TContextT>
void WorkStats::reduceThreadWorkStats(int this_rank,
                                      const spray::InsituPartition* partition,
                                      std::vector<TContextT>& tcontexts) {
  reset();

  // radiance blocks
  resetBlockCounters();

  for (auto& t : tcontexts) {
    reduceRayBlocks(t.getRadianceBlockIds());
  }

  updateReduceBuffer(partition);

  // shadow blocks
  resetBlockCounters();

  for (auto& t : tcontexts) {
    reduceRayBlocks(t.getShadowBlockIds());
  }

  updateReduceBuffer(partition);

  // cached blocks
  for (auto& t : tcontexts) {
    if (t.hasCachedBlock()) {
      addNumDomains(this_rank, 1);
      break;
    }
  }
}

}  // namespace insitu
}  // namespace spray
