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

#include <cstddef>

#include "glog/logging.h"

namespace spray {

class InfiniteCache {
 public:
  InfiniteCache();
  ~InfiniteCache();

  // max_aceh_size_ndomains is a don't care
  void initialize(int num_domains, int cache_size, bool insitu_mode = false);

  // returns true if hit, false if miss
  bool load(int domid, int* cache_block_id);

  int getCacheSize() const { return capacity_; }
  int getSize() const { return capacity_; }

 private:
  enum Status { HIT = -1, MISS = 0 };

 private:
  int capacity_;  ///< unit in number of domains
  int* status_;   ///< per-cache-entry status (-1: loaded, 0: not loaded)
};

}  // namespace spray
