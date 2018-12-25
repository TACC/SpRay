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

#include <vector>

#include "glog/logging.h"
#include "pbrt/memory3.h"

namespace spray {

template <typename T>
class ArenaQ;

template <typename T>
class ArenaQs {
 public:
  ArenaQs() : qs_(nullptr) {}
  ~ArenaQs() { delete[] qs_; }

  void resize(int array_size) {
#ifdef DRT_GLOG_CHECK
    CHECK_GT(array_size, 0);
#endif
    delete[] qs_;
    size_ = array_size;
    qs_ = new MemoryArena3[array_size];
  }

  T* allocate(int i, std::size_t count) {
#ifdef DRT_GLOG_CHECK
    CHECK_LT(i, size_);
#endif
    return qs_[i].template allocate<T>(count, false);
  }

  void copy(int i, const T* src) {
#ifdef DRT_GLOG_CHECK
    CHECK_LT(i, size_);
#endif
    qs_[i].template copy<T>(src);
  }

  void push(int i, const T& src) {
#ifdef DRT_GLOG_CHECK
    CHECK_LT(i, size_);
#endif
    qs_[i].template copy<T>(&src);
  }

  void reset(int i) {
#ifdef DRT_GLOG_CHECK
    CHECK_LT(i, size_);
#endif
    qs_[i].reset();
#ifdef DRT_GLOG_CHECK
    CHECK(qs_[i].empty());
    CHECK_EQ(qs_[i].size<T>(), 0);
#endif
  }

  void reset() {
    for (int i = 0; i < size_; ++i) reset(i);
  }

  bool empty(int i) const {
#ifdef DRT_GLOG_CHECK
    CHECK_LT(i, size_);
#endif
    return qs_[i].empty();
  }

  bool empty() const {
    for (int i = 0; i < size_; ++i) {
      if (!qs_[i].empty()) return false;
    }
    return true;
  }

  std::size_t size(int i) const {
#ifdef DRT_GLOG_CHECK
    CHECK_LT(i, size_);
#endif
    return qs_[i].template size<T>();
  }

  std::size_t size() const {
    std::size_t c = 0;
    for (int i = 0; i < size_; ++i) {
      c += qs_[i].template size<T>();
    }
    return c;
  }

  const std::list<MemBlock>& getBlocks(int i) const {
    return qs_[i].getBlocks();
  }

  // void splice(int i, ArenaQs* other) { qs_[i].splice(other->qs_[i]); }
  // void splice(int i, MemoryArena3* other_qs) { qs_[i].splice(other_qs[i]); }

  void incrementSize(int i, std::size_t count) {
#ifdef DRT_GLOG_CHECK
    CHECK_LT(i, size_);
#endif
    qs_[i].template incrementSize<T>(count);
  }

  int getNumValidBlocks(int i) const {
#ifdef DRT_GLOG_CHECK
    CHECK_LT(i, size_);
#endif
    return qs_[i].getNumValidBlocks();
  }

  const MemoryArena3& getQ(int i) const { return qs_[i]; }

 private:
  friend class ArenaQ<T>;
  MemoryArena3& get(int i) const { return qs_[i]; }

 private:
  MemoryArena3* qs_;
  int size_;
};

template <typename T>
class ArenaQ {
 public:
  void copy(const T* src) { q_.copy<T>(src); }

  void reset() { q_.reset(); }

  bool empty() const { return q_.empty(); }

  const std::list<MemBlock>& getBlocks() const { return q_.getBlocks(); }

  void splice(ArenaQ* other) { q_.splice(other->q_); }
  void splice(int i, ArenaQs<T>* other_qs) { q_.splice(other_qs->get(i)); }

  void incrementSize(std::size_t count) { q_.template incrementSize<T>(count); }
  std::size_t size() const { return q_.template size<T>(); }

  T* allocate(std::size_t count) {
    return q_.template allocate<T>(count, false);
  }

 private:
  MemoryArena3 q_;
};

}  // namespace spray
