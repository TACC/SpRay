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

#include <omp.h>
#include <numeric>

#include "pbrt/memory3.h"

#include "partition/block_buffer.h"

namespace spray {

BlockBuffer::BlockBuffer()
    : block_(nullptr), capacity_(nullptr), size_(nullptr), ndomains_(0) {}

BlockBuffer::~BlockBuffer() { cleanup(); }

void BlockBuffer::alloc(int ndomains, int max_nblocks) {
  cleanup();

  block_ = AllocAligned<MemBlock*>(ndomains);
  capacity_ = AllocAligned<int>(ndomains);
  size_ = AllocAligned<int>(ndomains);
  ndomains_ = ndomains;

  for (int i = 0; i < ndomains; ++i) {
    block_[i] = AllocAligned<MemBlock>(max_nblocks);
    capacity_[i] = max_nblocks;
    size_[i] = 0;
  }
}

void BlockBuffer::parallelResetAll() {
#pragma omp for collapse(2) schedule(static)
  for (int id = 0; id < ndomains_; ++id) {
    for (int i = 0; i < size_[id]; ++i) {
      MemBlock& blk = block_[id][i];
      blk.size = 0;
      // blk.capacity = 0;
      blk.buf = nullptr;
    }
  }
#pragma omp for schedule(static)
  for (int id = 0; id < ndomains_; ++id) {
    size_[id] = 0;
  }
}

void BlockBuffer::reset() {
  for (int id = 0; id < ndomains_; ++id) {
    for (int i = 0; i < size_[id]; ++i) {
      MemBlock& blk = block_[id][i];
      blk.size = 0;
      // blk.capacity = 0;
      blk.buf = nullptr;
    }
    size_[id] = 0;
  }
}

void BlockBuffer::reset(int id) {
  for (int i = 0; i < size_[id]; ++i) {
    MemBlock& blk = block_[id][i];
    blk.size = 0;
    // blk.capacity = 0;
    blk.buf = nullptr;
  }
  size_[id] = 0;
}

// NOTE: this function does not preserve existing items.
void BlockBuffer::resize(int id, int size) {
  if (capacity_[id] < size) {
    LOG(FATAL) << "capacity not big enough. Performance may degrade. required: "
               << size << " actual: " << capacity_[id];
    FreeAligned(block_[id]);

    block_[id] = AllocAligned<MemBlock>(size);
    capacity_[id] = size;
    size_[id] = size;

  } else {
    size_[id] = size;
    for (int i = 0; i < size; ++i) {
      MemBlock& blk = block_[id][i];
      blk.size = 0;
      blk.buf = nullptr;
    }
  }
}

void BlockBuffer::push(int id, std::size_t block_size, uint8_t* buf) {
  auto size = size_[id];
  CHECK(size < capacity_[id]) << size << " " << capacity_[id] << " " << id;
  MemBlock& b = block_[id][size];
  b.size = block_size;
  // b.capacity = block_capacity;
  b.buf = buf;
  size_[id] = size + 1;
}

void BlockBuffer::set(int id, int idx, std::size_t block_size, uint8_t* buf) {
  CHECK_LT(idx, size_[id]);
  MemBlock& b = block_[id][idx];
  b.size = block_size;
  // b.capacity = block_capacity;
  b.buf = buf;
#ifdef SPRAY_GLOG_CHECK
  CHECK(b.size);
  CHECK_NOTNULL(b.buf);
#endif
}

void BlockBuffer::append(int id, const MemBlock& block) {
  int idx = size_[id];
  CHECK_LT(idx, capacity_[id]);
  block_[id][idx] = block;
  size_[id] = idx + 1;
}

void BlockBuffer::append(int id, std::size_t block_size,
                         std::size_t block_capacity, uint8_t* buf) {
  int idx = size_[id];
  CHECK_LT(idx, capacity_[id]);
  MemBlock& b = block_[id][idx];
  b.size = block_size;
  // b.capacity = block_capacity;
  b.buf = buf;
  size_[id] = idx + 1;
}

void BlockBuffer::setLastBlock(int id, std::size_t block_size, uint8_t* buf) {
#ifdef SPRAY_GLOG_CHECK
  CHECK(size_[id] > 0);
#endif
  MemBlock& b = block_[id][size_[id] - 1];
  b.size = block_size;
#ifdef SPRAY_GLOG_CHECK
  CHECK(b.buf == nullptr);
#endif
  b.buf = buf;
}

int BlockBuffer::getNumBlocks(int id) const { return size_[id]; }

const MemBlock& BlockBuffer::getBlock(int id, int idx) const {
  return block_[id][idx];
}

const MemBlock& BlockBuffer::getLastBlock(int id) const {
  return block_[id][size_[id] - 1];
}

void BlockBuffer::cleanup() {
  for (int i = 0; i < ndomains_; ++i) {
    FreeAligned(block_[i]);
  }
  FreeAligned(block_);
  FreeAligned(capacity_);
  FreeAligned(size_);
}

}  // namespace spray
