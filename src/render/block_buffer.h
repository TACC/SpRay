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

#include "pbrt/memory3.h"

namespace spray {

class BlockBuffer {
 public:
  BlockBuffer();
  ~BlockBuffer();

 public:
  //! Create 2D array of memory blocks (ndomains * max_nblocks).
  void alloc(int ndomains, int max_nblocks = 8192);
  void reset();
  void reset(int id);
  void parallelResetAll();

 public:
  // NOTE: this function does not preserve existing items.
  void resize(int id, int size);

 public:
  // add and update size
  void push(int id, std::size_t block_size, uint8_t* buf);

  void set(int id, int idx, std::size_t block_size, uint8_t* buf);
  void setLastBlock(int id, std::size_t block_size, uint8_t* buf);

  template <typename T>
  void incrementLastBlockSize(int id, std::size_t count) {
    MemBlock& b = block_[id][size_[id] - 1];
#ifdef SPRAY_GLOG_CHECK
    CHECK_NOTNULL(b.buf);
#endif
    b.size += (count * sizeof(T));
  }

  void append(int id, const MemBlock& block);

  void append(int id, std::size_t block_size, std::size_t block_capacity,
              uint8_t* buf);

  int getNumBlocks(int id) const;

 public:
  template <typename QItemT>
  std::size_t size(int id) const {
#ifdef SPRAY_GLOG_CHECK
    CHECK_LT(id, ndomains_);
#endif
    std::size_t num = 0;
    for (int i = 0; i < size_[id]; ++i) {
      const MemBlock& b = block_[id][i];
      num += ((QItemT*)(b.buf + b.size) - (QItemT*)(b.buf));
    }
    return num;
  }

  template <typename QItemT>
  std::size_t size() const {
    std::size_t num = 0;
    for (int i = 0; i < ndomains_; ++i) {
      num += size<QItemT>(i);
    }
    return num;
  }

  bool empty(int id) const {
    std::size_t num = 0;
    for (int i = 0; i < size_[id]; ++i) {
      const MemBlock& b = block_[id][i];
      if (b.size) return false;
    }
    return true;
  }

  bool empty() const {
    std::size_t num = 0;
    for (int i = 0; i < ndomains_; ++i) {
      if (!empty(i)) return false;
    }
    return true;
  }

 public:
  const MemBlock& getBlock(int id, int idx) const;
  const MemBlock& getLastBlock(int id) const;

 private:
  void cleanup();

 private:
  //! Per-domain memory blocks
  MemBlock** block_;

  //! Per-domain max number of allocated memory blocks.
  int* capacity_;

  //! Per-domain number of valid memory blocks
  int* size_;

  //! Number of domains.
  int ndomains_;
};

}  // namespace spray

