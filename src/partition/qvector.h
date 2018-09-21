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

#include <cstdint>
#include <queue>
#include <vector>

#include "glog/logging.h"

namespace spray {

template <typename T>
class QVector {
 public:
  void resize(std::size_t size) {
    CHECK(empty());
    qs_.resize(size);
  }
  void push(int i, T& data) { qs_[i].push(data); }

  std::queue<T>* getQ(int i) { return &qs_[i]; }
  std::size_t size(int i) const { return qs_[i].size(); }

  bool empty() const {
    for (const auto& q : qs_) {
      if (!q.empty()) return false;
    }
    return true;
  }
  uint64_t size() const {
    uint64_t s = 0;
    for (const auto& q : qs_) {
      s += q.size();
    }
    return s;
  }

  bool empty(int i) const { return qs_[i].empty(); }
  void pop(int i) { return qs_[i].pop(); }
  T& front(int i) { return qs_[i].front(); }

  void flush() {
    for (auto& q : qs_) {
      while (!q.empty()) {
        q.pop();
      }
    }
  }

 private:
  std::vector<std::queue<T>> qs_;
};

}  // namespace spray
