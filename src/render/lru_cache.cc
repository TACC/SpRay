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

#include "render/lru_cache.h"

#include "glog/logging.h"

namespace spray {

LruCache::LruCache() : size_(0), capacity_(0) {}

void LruCache::flush() {
  size_ = 0;
  capacity_ = 0;
  ndomains_ = 0;
  blocks_.clear();
  id_to_block_.clear();
  for (std::size_t i = 0; i < status_.size(); ++i) {
    status_[i] = MISS;
  }
}

// max_aceh_size_ndomains is a don't care
void LruCache::init(int num_domains, int cache_size) {
  flush();
  //
  CHECK_NE(cache_size, 0);

  if (cache_size < 0 || cache_size >= num_domains) {
    capacity_ = num_domains;
  } else {
    capacity_ = cache_size;
  }

  CHECK_LE(capacity_, num_domains);
  CHECK_GT(capacity_, 0);

  ndomains_ = num_domains;

  CHECK_GT(ndomains_, 0);

  status_.resize(ndomains_);
  for (std::size_t i = 0; i < status_.size(); ++i) {
    status_[i] = MISS;
  }
}

bool LruCache::load(int domid, int* cache_block_id) {
#ifdef SPRAY_GLOG_CHECK
  CHECK_GT(capacity_, 0);
  CHECK_LT(domid, ndomains_);
  CHECK_GE(domid, 0);
#endif
  bool hit;

  if (status_[domid] == HIT) {  // loaded
    hit = true;

    // get block
#ifdef SPRAY_GLOG_CHECK
    CHECK(id_to_block_.find(domid) != id_to_block_.end()) << domid;
#endif
    BlockIter it = id_to_block_[domid];

    CacheBlock b = *it;
    *cache_block_id = b.block;
#ifdef SPRAY_GLOG_CHECK
    CHECK_EQ(domid, b.domain);
#endif

    // erase block
    blocks_.erase(it);
    // id_to_block_.erase(domid);

    // promote block
    blocks_.push_back(b);
    BlockIter mru_it = blocks_.end();
    --mru_it;

    id_to_block_[domid] = mru_it;

  } else {  // not loaded

    hit = false;
    status_[domid] = HIT;

    if (size_ < capacity_) {  // not full
#ifdef SPRAY_GLOG_CHECK
      CHECK(id_to_block_.find(domid) == id_to_block_.end());
      CHECK_EQ(blocks_.size(), size_);
      CHECK_EQ(id_to_block_.size(), size_);
#endif

      CacheBlock b;
      b.block = size_;
      b.domain = domid;

      blocks_.push_back(b);

      BlockIter mru_it = blocks_.end();
      --mru_it;
      id_to_block_[domid] = mru_it;

      *cache_block_id = b.block;

      ++size_;

    } else {  // full
#ifdef SPRAY_GLOG_CHECK
      CHECK_EQ(blocks_.size(), capacity_);
      CHECK_EQ(id_to_block_.size(), capacity_);
      CHECK_EQ(size_, capacity_);
#endif
      // evict lru block
      CacheBlock old_blk = blocks_.front();

      blocks_.pop_front();
      id_to_block_.erase(old_blk.domain);

      status_[old_blk.domain] = MISS;

      // insert new block
      CacheBlock new_blk;
      new_blk.block = old_blk.block;
      new_blk.domain = domid;

      blocks_.push_back(new_blk);
      BlockIter mru_it = blocks_.end();
      --mru_it;

      id_to_block_[domid] = mru_it;

      *cache_block_id = new_blk.block;
    }
  }

#ifdef SPRAY_GLOG_CHECK
  CHECK_LT(*cache_block_id, capacity_);
  CHECK_GE(*cache_block_id, 0);
  CHECK_LE(size_, capacity_);
  CHECK_EQ(size_, blocks_.size());
  CHECK_EQ(size_, id_to_block_.size());
  CHECK_EQ(status_[domid], HIT);
  CHECK(id_to_block_.find(domid) != id_to_block_.end());

  MapIter map_it = id_to_block_.find(domid);
  CHECK(map_it != id_to_block_.end());
  CacheBlock test_blk = *(map_it->second);
  CHECK_EQ(test_blk.block, *cache_block_id);
  CHECK_EQ(test_blk.domain, domid);
#endif

  return hit;
}

}  // namespace spray
