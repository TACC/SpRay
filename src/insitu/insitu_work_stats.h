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

#include "partition/data_partition.h"
#include "renderers/spray.h"
#include "utils/comm.h"

namespace spray {
namespace insitu {

class ThreadWorkStats {
 public:
  void resize(int nranks) { reduce_buf_.resize(nranks, 0); }

  void reset() {
    for (auto& r : reduce_buf_) r = 0;
  }

  void addNumDomains(int dest_rank, int num_blocks) {
    reduce_buf_[dest_rank] += num_blocks;
  }

  std::size_t getReduceBufSize() const { return reduce_buf_.size(); }
  int getReduceBuf(std::size_t i) const { return reduce_buf_[i]; }

 private:
  std::vector<int> reduce_buf_;  // per-rank
};

class WorkStats {
 private:
  struct ScatterEntry {
    int world_num_blocks_to_proc;
    int rank_num_blocks_to_proc;
  };

  // organizatin
  //   number of blocks to recv for each rank
  std::vector<int> reduce_buf_;            // per-rank
  std::vector<ScatterEntry> scatter_buf_;  // per-rank

 public:
  void resize() {
    int num_ranks = mpi::size();
    reduce_buf_.resize(num_ranks, 0);
    scatter_buf_.resize(num_ranks);
  }

  void reset() {
    num_blocks_to_recv_ = 0;
    for (auto& r : reduce_buf_) r = 0;
  }

 private:
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
  void addNumDomains(int dest_rank, int num_blocks) {
    reduce_buf_[dest_rank] += num_blocks;
  }

 public:
  void merge(const ThreadWorkStats& stats) {
#ifdef SPRAY_GLOG_CHECK
    CHECK_GT(reduce_buf_.size(), 0);
    CHECK_EQ(reduce_buf_.size(), stats.getReduceBufSize());
#endif
    for (std::size_t i = 0; i < reduce_buf_.size(); ++i) {
      reduce_buf_[i] += stats.getReduceBuf(i);
    }
  }
};

}  // namespace insitu
}  // namespace spray
