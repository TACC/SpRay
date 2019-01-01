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

#include "pbrt/memory.h"

#include "insitu/insitu_isector.h"
#include "insitu/insitu_ray.h"
#include "insitu/insitu_vbuf.h"
#include "insitu/insitu_work_stats.h"
#include "render/config.h"
#include "render/qvector.h"
#include "render/scene.h"

namespace spray {

class InsituPartition;

namespace insitu {

class ThreadWorkStats {
 public:
  void registerRadianceRayBlock(int id) { radiance_block_ids_.push(id); }
  void registerCachedRayBlock(bool hasCachedBlock) {
    has_cached_block_ = hasCachedBlock;
  }
  void registerShadowRayBlock(int id) { shadow_block_ids_.push(id); }

  std::queue<int>* getRadianceBlockIds() { return &radiance_block_ids_; }
  std::queue<int>* getShadowBlockIds() { return &shadow_block_ids_; }
  bool hasCachedBlock() const { return has_cached_block_; }

  bool empty() const {
    return (radiance_block_ids_.empty() && shadow_block_ids_.empty());
  }

 private:
  std::queue<int> radiance_block_ids_;
  std::queue<int> shadow_block_ids_;
  bool has_cached_block_;
};

template <typename CacheT, typename ShaderT, typename SceneT>
class TContext {
 public:
  TContext() {
    mem_in_ = &mem_0_;
    mem_out_ = &mem_1_;
  }
  ~TContext() { resetMems(); }

 public:
  void init(const Config& cfg, int ndomains,
            const spray::InsituPartition* partition, SceneT* scene, VBuf* vbuf,
            spray::HdrImage* image);

  void resetMems() {
    mem_in_->Reset();
    mem_out_->Reset();
  }

  Ray* allocMemIn(std::size_t size) {
    Ray* mem = mem_in_->Alloc<Ray>(size, false);
    CHECK_NOTNULL(mem);
    return mem;
  }

  void resetMemIn() { mem_in_->Reset(); }

  void resetAndSwapMems() {
    mem_in_->Reset();
    std::swap(mem_in_, mem_out_);
  }

  spray::MemoryArena* getMemIn() { return mem_in_; }

 private:
  spray::MemoryArena* mem_in_;   // rays being processed
  spray::MemoryArena* mem_out_;  // rays being generated
  spray::MemoryArena mem_0_;
  spray::MemoryArena mem_1_;

  const spray::InsituPartition* partition_;
  SceneT* scene_;  // TODO: use const
  VBuf* vbuf_;
  spray::HdrImage* image_;

  int num_domains_;

 public:
  void isectDomains(Ray* ray) {
#ifdef SPRAY_BACKGROUND_COLOR_BLACK
    isector_.intersect(num_domains_, scene_, ray, &rqs_);
#else
    isector_.intersect(num_domains_, scene_, ray, &rqs_, &bg_retire_q_);
#endif
  }

 public:
  void populateRadWorkStats();
  void populateWorkStats();
  std::queue<int>* getRadianceBlockIds() {
    return work_stats_.getRadianceBlockIds();
  }
  std::queue<int>* getShadowBlockIds() {
    return work_stats_.getShadowBlockIds();
  }
  bool hasCachedBlock() const { return work_stats_.hasCachedBlock(); }

 private:
  ThreadWorkStats work_stats_;  // number of domain blocks to process

 public:
  void setSceneInfo(SceneInfo& sinfo) { sinfo_ = sinfo; }
  void processRays(int id, SceneInfo& sinfo);
  void processRays2() {
    procFrq2();
    procFsq2();
  }
  void updateVisBuf();
  void updateTBuf();
  void updateOBuf();
  void genRays(int id, int ray_depth);
  void isectRecvRad(int id, Ray* ray);
  void occlRecvShad(int id, Ray* ray);
  void isectCachedRq(int ray_depth);
  void procRetireQ();
  void retireBackground();
  void sendRays(bool shadow, int id, Ray* rays);
  void updateTBufWithCached();
  void processCached(int ray_depth);

  void checkQs() const {
    CHECK(rqs_.empty());
    CHECK(sqs_.empty());
    CHECK(rq2_.empty());
    CHECK(sq2_.empty());
    CHECK(frq2_.empty());
    CHECK(fsq2_.empty());
    CHECK(retire_q_.empty());
    CHECK(cached_rq_.empty());
    CHECK(reduced_cached_rq_.empty());
  }

 private:
  void filterSq2(int id);
  void filterRq2(int id);
  void procFsq2();
  void procFrq2();

 private:
  struct IsectInfo {
    Ray* ray;
    spray::RTCRayIntersection* isect;
  };

  struct OcclInfo {
    int samid;
    int light;
  };

  struct IsectCacheItem {
    int domain_id;
    spray::RTCRayIntersection* isect;
    Ray* ray;
  };

  struct OcclCacheItem {
    int domain_id;
    Ray* ray;
  };

  IsectInfo isect_;
  std::queue<IsectInfo> isects_;
  std::queue<IsectInfo> reduced_isects_;

  OcclInfo occl_;
  std::queue<OcclInfo> occls_;

  SceneInfo sinfo_;

  IsectCacheItem icache_item_;
  OcclCacheItem ocache_item_;

 public:
  bool isLocalQsEmpty(int id) const {
    return (rqs_.empty(id) && sqs_.empty(id));
  }

  std::size_t getRqSize(int id) const { return rqs_.size(id); }
  std::size_t getSqSize(int id) const { return sqs_.size(id); }

 public:
  // debug
  void flushRqs() { rqs_.flush(); }

 private:
  OcclInfo occl_info_;
  IsectInfo isect_info_;

  ShaderT shader_;

  Isector<CacheT, SceneT> isector_;
  spray::RTCRayIntersection rtc_isect_;
  RTCRay rtc_ray_;

  spray::QVector<Ray*> rqs_;
  spray::QVector<Ray*> sqs_;

  std::queue<Ray*> rq2_;
  std::queue<Ray*> sq2_;

  std::queue<IsectCacheItem> frq2_;  // filtered
  std::queue<OcclCacheItem> fsq2_;   // filtered

  std::queue<Ray*> retire_q_;     ///< Retire queue for foreground colors.
  std::queue<Ray*> bg_retire_q_;  ///< Retire queue for background colors.

  std::queue<IsectCacheItem> cached_rq_;
  std::queue<IsectCacheItem> reduced_cached_rq_;

  double one_over_num_pixel_samples_;
};

}  // namespace insitu
}  // namespace spray

#define SPRAY_INSITU_TCONTEXT_INL
#include "insitu/insitu_tcontext.inl"
#undef SPRAY_INSITU_TCONTEXT_INL

