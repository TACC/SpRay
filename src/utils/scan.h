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

#include <vector>

namespace spray {

template <typename T>
class InclusiveScan {
 public:
  void resize(int size) { buf_.resize(size); }
  std::size_t size() const { return buf_.size(); }

  void reset(T val) {
    for (std::size_t i = 0; i < buf_.size(); ++i) buf_[i] = val;
  }

  void set(int i, T value) {
#ifdef SPRAY_GLOG_CHECK
    CHECK_LT(i, buf_.size());
#endif
    buf_[i] = value;
  }

  T get(int i) const {
#ifdef SPRAY_GLOG_CHECK
    CHECK_LT(i, buf_.size());
#endif
    return buf_[i];
  }

  T sum() const { return buf_[buf_.size() - 1]; }

  void scan() {
#ifdef SPRAY_GLOG_CHECK
    CHECK_GT(buf_.size(), 0);
#endif
    for (std::size_t i = 1; i < buf_.size(); ++i) {
      buf_[i] = buf_[i - 1] + buf_[i];
    }
  }

 private:
  std::vector<T> buf_;
};

}  // namespace spray

