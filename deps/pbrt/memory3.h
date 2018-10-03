
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#pragma once

// core/memory.h*
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <list>

#include "glog/logging.h"
#include "pbrt/memory.h"
#include "pbrt/pbrt.h"

namespace spray {

class
#ifdef PBRT_HAVE_ALIGNAS
    alignas(PBRT_L1_CACHE_LINE_SIZE)
#endif
        MemBlock {
 public:
  MemBlock() : size(0), used_capacity(0), total_capacity(0), buf(nullptr) {}

  template <typename T>
  static std::size_t getSize(const MemBlock &b) {
    uint8_t *end = b.buf + b.size;
    return ((T *)end - (T *)b.buf);
  }

  uint8_t *buf;      // never update the address, only updated by MemoryArena3
  std::size_t size;  //!< number of items in bytes
 private:
  friend class MemoryArena3;
  std::size_t used_capacity;   //!< used capacity in bytes
                               /**
                                * Total capacity in bytes.
                                * free capacity = total_capacity - used_capacity.
                                */
  std::size_t total_capacity;  //!< total capacity in bytes
};

class
#ifdef PBRT_HAVE_ALIGNAS
    alignas(PBRT_L1_CACHE_LINE_SIZE)
#endif  // PBRT_HAVE_ALIGNAS
        MemoryArena3 {
 public:
  // typedef std::list<MemBlock>::iterator MemBlockIterator;
  // MemoryArena3 Public Methods
  // MemoryArena3(size_t block_size_in_bytes = 65536)
  MemoryArena3() : block_size_(262144) {}
  // MemoryArena3(size_t block_size_in_bytes = 262144)
  MemoryArena3(size_t block_size_in_bytes) : block_size_(block_size_in_bytes) {}
  ~MemoryArena3() {
    for (auto &block : used_blocks_) FreeAligned((void *)block.buf);
    for (auto &block : available_blocks_) FreeAligned((void *)block.buf);
  }
  void *allocateBlock(std::size_t num_bytes) {
// Round up _nBytes_ to minimum machine alignment
#if __GNUC__ == 4 && __GNUC_MINOR__ < 9
    // gcc bug: max_align_t wasn't in std:: until 4.9.0
    const int align = alignof(::max_align_t);
#elif !defined(PBRT_HAVE_ALIGNOF)
    const int align = 16;
#else
    const int align = alignof(std::max_align_t);
#endif
#ifdef PBRT_HAVE_CONSTEXPR
    static_assert(IsPowerOf2(align), "Minimum alignment not a power of two");
#endif
    std::size_t aligned_nbytes = (num_bytes + align - 1) & ~(align - 1);
#ifdef DRT_GLOG_CHECK
    CHECK_GT(aligned_nbytes, 0);
#endif

    if (used_blocks_.empty()) {
      MemBlock b;
      used_blocks_.push_back(b);
    }

    MemBlock &current_block = used_blocks_.back();

    // current block not big enough
    if ((current_block.used_capacity + aligned_nbytes) >
        current_block.total_capacity) {
      // Get new block of memory for _MemoryArena_
      bool reuse = false;

      // Try to get memory block from _availableBlocks_
      for (auto iter = available_blocks_.begin();
           iter != available_blocks_.end(); ++iter) {
#ifdef DRT_GLOG_CHECK
        CHECK_GT(iter->total_capacity, 0);
        CHECK_NOTNULL(iter->buf);
#endif
        if (iter->total_capacity >= aligned_nbytes) {
          if (current_block.buf) {
            // push available block to used blocks
            // current block is not big enough
            // and delete available block
            iter->size = 0;
            iter->used_capacity = 0;

            used_blocks_.push_back(*iter);
          } else {
            // update current block
            // (i.e. copy available block to current block and delete available
            // block)
            // since current block has not been allocated yet
            current_block.buf = iter->buf;
            current_block.size = 0;
            current_block.used_capacity = 0;
            current_block.total_capacity = iter->total_capacity;
          }
          available_blocks_.erase(iter);
          reuse = true;
          break;
        }
      }

      // allocate current block if it does not exist
      if (!reuse) {
        std::size_t total_capacity = std::max(aligned_nbytes, block_size_);
        uint8_t *buffer = AllocAligned<uint8_t>(total_capacity);
#ifdef DRT_GLOG_CHECK
        CHECK_GT(total_capacity, 0);
        CHECK_NOTNULL(buffer);
#endif
        if (current_block.buf) {
          // push new block, current block is not big enough
          MemBlock b;
          b.buf = buffer;
          b.size = 0;
          b.used_capacity = 0;
          b.total_capacity = total_capacity;

          used_blocks_.push_back(b);
        } else {
          // update current block
          // since current block has not been allocated yet
          current_block.buf = buffer;
          current_block.size = 0;
          current_block.used_capacity = 0;
          current_block.total_capacity = total_capacity;
        }
      }
    }  // else current block big enough

#ifdef DRT_GLOG_CHECK
    CHECK_GT(used_blocks_.size(), 0);
#endif

    MemBlock &b = used_blocks_.back();

    void *ret = b.buf + b.used_capacity;
    b.used_capacity += aligned_nbytes;

#ifdef DRT_GLOG_CHECK
    for (const auto &blk : used_blocks_) {
      CHECK_NOTNULL(blk.buf);
    }
    CHECK_NOTNULL(b.buf);
    CHECK_LE(b.size, b.used_capacity);
    CHECK_LE(b.size, b.total_capacity);
    CHECK_LE(b.used_capacity, b.total_capacity);
#endif
    return ret;
  }

  // increment current block size
  template <typename T>
  void incrementSize(std::size_t num_items) {
#ifdef DRT_GLOG_CHECK
    CHECK_GT(used_blocks_.size(), 0);
#endif

    MemBlock &b = used_blocks_.back();
    b.size += (num_items * sizeof(T));

#ifdef DRT_GLOG_CHECK
    CHECK_NOTNULL(b.buf);
    CHECK_LE(b.size, b.used_capacity);
    CHECK_LE(b.size, b.total_capacity);
    CHECK_LE(b.used_capacity, b.total_capacity);
#endif
  }

  template <typename T>
  T *allocate(size_t n = 1, bool runConstructor = false) {
#ifdef DRT_GLOG_CHECK
    CHECK_GT(n, 0);
    CHECK_GT(n * sizeof(T), 0);
#endif

    T *ret = (T *)allocateBlock(n * sizeof(T));

#ifdef DRT_GLOG_CHECK
    CHECK_NOTNULL(ret);
    CHECK_GT(used_blocks_.back().total_capacity, 0);
#endif
    if (runConstructor)
      for (size_t i = 0; i < n; ++i) new (&ret[i]) T();
    return ret;
  }

  void reset() {
#ifdef DRT_GLOG_CHECK
    for (auto &blk : used_blocks_) {
      CHECK_NOTNULL(blk.buf);
      CHECK_GT(blk.total_capacity, 0);
      CHECK_LE(blk.size, blk.used_capacity);
      CHECK_LE(blk.size, blk.total_capacity);
      CHECK_LE(blk.used_capacity, blk.total_capacity);
    }
#endif
    available_blocks_.splice(available_blocks_.begin(), used_blocks_);
#ifdef DRT_GLOG_CHECK
    CHECK(used_blocks_.empty());
    CHECK(empty());
#endif
  }

  size_t totalAllocated() const {
    size_t total = 0;
    for (const auto &used_blk : used_blocks_) {
#ifdef DRT_GLOG_CHECK
      CHECK_NOTNULL(used_blk.buf);
#endif
      total += used_blk.total_capacity;
    }
    for (const auto &avail_blk : available_blocks_) {
#ifdef DRT_GLOG_CHECK
      CHECK_NOTNULL(avail_blk.buf);
#endif
      total += avail_blk.total_capacity;
    }
    return total;
  }

  bool inResetState() const { return used_blocks_.empty(); }

  //! Returns the total number of items of type T
  template <typename T>
  std::size_t size() const {
    std::size_t total = 0;
    for (const auto &blk : used_blocks_) {
#ifdef DRT_GLOG_CHECK
      CHECK_NOTNULL(blk.buf);
#endif
      total += ((T *)(blk.buf + blk.size) - (T *)blk.buf);
    }
    return total;
  }

  //! Returns the total number of items of type T
  std::size_t sizeInBytes() const {
    std::size_t total = 0;
    for (const auto &blk : used_blocks_) {
      total += blk.size;
    }
    return total;
  }

  bool empty() const {
    for (const auto blk : used_blocks_) {
      if (blk.size) return false;
    }
    return true;
  }

  template <typename T>
  void copy(const T *item) {
    // std::size_t size = sizeof(T);
    // T *dest = (T *)allocateBlock(size);
    T *dest = allocate<T>(1, false);
#ifdef DRT_GLOG_CHECK
    {
      CHECK_NOTNULL(dest);
      CHECK_GT(used_blocks_.back().total_capacity, 0);
    }
#endif
    std::memcpy((void *)dest, (const void *)item, sizeof(T));
    incrementSize<T>(1);
  }

  const std::list<MemBlock> &getBlocks() const { return used_blocks_; }

  int getNumValidBlocks() const {
    int c = 0;
    for (auto &b : used_blocks_) {
      if (b.size) {
#ifdef DRT_GLOG_CHECK
        CHECK_NOTNULL(b.buf);
#endif
        ++c;
      }
    }
    return c;
  }

  void splice(MemoryArena3 &arena) {
#ifdef DRT_GLOG_CHECK
    for (auto &blk : arena.used_blocks_) {
      CHECK_NOTNULL(blk.buf);
      CHECK_GT(blk.total_capacity, 0);
    }
#endif
    used_blocks_.splice(used_blocks_.begin(), arena.used_blocks_);
  }

  // void swap(MemoryArena3 &arena) { used_blocks_.swap(arena.used_blocks); }

 private:
  // MemoryArena3(const MemoryArena3 &) = delete;
  // MemoryArena3 &operator=(const MemoryArena3 &) = delete;

  std::list<MemBlock> used_blocks_;
  std::list<MemBlock> available_blocks_;

  const size_t block_size_;  // block size in bytes.
};

}  // namespace spray

