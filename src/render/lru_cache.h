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

#include <list>
#include <map>
#include <vector>

#include "glog/logging.h"

namespace spray {

struct CacheBlock {
  int block;
  int domain;
};

class LruCache {
 public:
  typedef std::list<CacheBlock>::iterator BlockIter;   // list of domain IDs
  typedef std::map<int, BlockIter>::iterator MapIter;  // ID-to-LRU iterator

 public:
  LruCache();
  // max_aceh_size_ndomains is a don't care
  void init(int num_domains, int cache_size);

  // returns true if hit, false if miss
  bool load(int domid, int* cache_block_id);

  // void setLoaded(int cache_block_id);

  int getCacheSize() const { return capacity_; }
  int getSize() const { return size_; }

  std::size_t getListSize() const { return blocks_.size(); }
  std::size_t getMapSize() const { return id_to_block_.size(); }

  const std::list<CacheBlock>& getBlocks() const { return blocks_; }
  const std::map<int, BlockIter>& getMap() const { return id_to_block_; }

  const std::vector<int>& getStatus() const { return status_; }
  CacheBlock getBlock(int domid) {
    CHECK(id_to_block_.find(domid) != id_to_block_.end());
    BlockIter it = id_to_block_[domid];
    return (*it);
  }

 private:
  void flush();

 private:
  enum Status { HIT = -1, MISS = 0 };

  int size_;
  int capacity_;
  int ndomains_;

  // cache blocks, front (lru) --- back (mru)
  std::list<CacheBlock> blocks_;
  std::map<int, BlockIter> id_to_block_;  ///< domain ID-to-LRU_iterator map
  std::vector<int> status_;  ///< per-domain loaded status (-1 or 0)
};

}  // namespace spray
