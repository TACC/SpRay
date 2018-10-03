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

#include "glog/logging.h"
#include "pbrt/memory.h"

namespace spray {
namespace baseline {

class QStats {
 public:
  QStats() : counts_(nullptr) {}
  ~QStats() { FreeAligned(counts_); }

  void init(int ndomains) {
    ndomains_ = ndomains;
    FreeAligned(counts_);

    counts_ = AllocAligned<int64_t>(ndomains);
  }

  void reset() {
    for (int i = 0; i < ndomains_; ++i) {
      counts_[i] = 0;
    }
  }

  void add(int i, const QStats& stats) {
#ifdef SPRAY_GLOG_CHECK
    CHECK_LT(i, ndomains_);
#endif
    counts_[i] += stats.counts_[i];
  }

  void set(int i, int64_t count) {
#ifdef SPRAY_GLOG_CHECK
    CHECK_LT(i, ndomains_);
#endif
    counts_[i] = count;
  }

  int64_t get(int i) const { return counts_[i]; }

 private:
  int64_t* counts_;
  int ndomains_;
};

}  // namespace baseline
}  // namespace spray
