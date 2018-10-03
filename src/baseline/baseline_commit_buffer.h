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

#include "glm/glm.hpp"
#include "pbrt/memory3.h"

#include "baseline/baseline_comm.h"
#include "baseline/baseline_ray.h"
#include "display/image.h"
#include "utils/util.h"

namespace spray {
namespace baseline {

struct SPRAY_ALIGN(16) CommitB {
  int pixid;
  float color[3];
};

class CommitBufferB {
 public:
  void commit(int pixid, const glm::vec3& color) {
    CommitB commit;
    commit.pixid = pixid;
    commit.color[0] = color[0];
    commit.color[1] = color[1];
    commit.color[2] = color[2];
    arena_.copy<CommitB>(&commit);
  }
  void reset() { arena_.reset(); }
  bool empty() const { return arena_.empty(); }

  void retire(spray::HdrImage* image, int nsamples) {
    double scale = 1.0 / (double)nsamples;
    const std::list<MemBlock>& blocks = arena_.getBlocks();
    for (auto& blk : blocks) {
      // get block info
      CommitB* buf = (CommitB*)blk.buf;
      std::size_t num = spray::util::getNumOfItems<CommitB>(blk.buf, blk.size);

      for (std::size_t n = 0; n < num; ++n) {
        CommitB& commit = buf[n];
        int pixid = commit.pixid;
        image->add(pixid, commit.color, scale);
      }
    }
  }

 private:
  spray::MemoryArena3 arena_;
};

}  // namespace baseline
}  // namespace spray

