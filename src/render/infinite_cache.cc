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

#include "render/infinite_cache.h"

#include "glog/logging.h"
#include "pbrt/memory.h"

#include "render/trimesh_buffer.h"

namespace spray {

InfiniteCache::InfiniteCache() : capacity_(0), status_(nullptr) {}
InfiniteCache::~InfiniteCache() { FreeAligned(status_); }

// max_aceh_size_ndomains is a don't care
void InfiniteCache::initialize(int num_domains, int cache_size,
                               bool insitu_mode) {
  //
  CHECK(cache_size < 0 || cache_size >= num_domains || insitu_mode);
  capacity_ = num_domains;

  status_ = AllocAligned<int>(capacity_);
  CHECK_NOTNULL(status_);

  for (int i = 0; i < capacity_; ++i) {
    status_[i] = MISS;
  }
}

bool InfiniteCache::load(int domid, int* cache_block_id) {
#ifdef SPRAY_GLOG_CHECK
  CHECK_LT(domid, capacity_);
#endif
  *cache_block_id = domid;
  if (status_[domid] == HIT) {
    return true;
  }
#ifdef SPRAY_GLOG_CHECK
  CHECK_LT(domid, capacity_);
#endif
  status_[domid] = HIT;
  return false;
}

}  // namespace spray
