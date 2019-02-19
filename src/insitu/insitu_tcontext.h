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

template <typename ShaderT>
class TContext {
  typedef std::queue<Ray*> RayQ;
  typedef spray::QVector<Ray*> Qvector;

 public:
  typedef typename ShaderT::SceneType SceneType;

  TContext() {
    mem_in_ = &mem_0_;
    mem_out_ = &mem_1_;
  }
  ~TContext() { resetMems(); }

  void init(const Config& cfg, int ndomains,
            const spray::InsituPartition* partition, SceneType* scene,
            VBuf* vbuf, spray::HdrImage* image);

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

  void isectDomains(Ray* ray) {
#ifdef SPRAY_GLOG_CHECK
    CHECK_LT(ray->pixid, image_->w * image_->h);
#endif
    isector_.intersect(scene_, ray, &rqs_);
  }

  void processRays(int rank, int ray_depth);
  void pushRadianceRay(int id, Ray* ray) { rqs_.push(id, ray); }
  void pushShadowRay(int id, Ray* ray) { sqs_.push(id, ray); }

  void compositeThreadTbufs(int tid, std::vector<VBuf>* vbufs);
  void compositeThreadObufs(int tid, std::vector<VBuf>* vbufs);

  void populateRadWorkStats();
  void populateWorkStats();
  std::queue<int>* getRadianceBlockIds() {
    return work_stats_.getRadianceBlockIds();
  }
  std::queue<int>* getShadowBlockIds() {
    return work_stats_.getShadowBlockIds();
  }
  bool hasCachedBlock() const { return work_stats_.hasCachedBlock(); }

  void resolveSecondaryRays(const VBuf& vbuf) {
    resolveRadiances(vbuf);
    resolveShadows(vbuf);
  }

  void retireUntouched();
  void retireShadows(const VBuf& vbuf);

  void sendRays(bool shadow, int id, Ray* rays);

  bool isLocalQsEmpty(int id) const {
    return (rqs_.empty(id) && sqs_.empty(id));
  }

  std::size_t getRqSize(int id) const { return rqs_.size(id); }
  std::size_t getSqSize(int id) const { return sqs_.size(id); }

  // debug
  void flushRqs() { rqs_.flush(); }

  void checkQs() const {
    CHECK(rqs_.empty());
    CHECK(sqs_.empty());
    CHECK(rq2_.empty());
    CHECK(sq2_.empty());
    CHECK(frq2_.empty());
    CHECK(fsq2_.empty());
    CHECK(retire_q_.empty());
    CHECK(cached_rq_.empty());
  }

 private:
  void processRadiance(int id, int ray_depth, Ray* ray);
  void processShadow(int id, Ray* ray);

  void filterSq2(int id);
  void filterRq2(int id);

  void resolveRadiances(const VBuf& vbuf);
  void resolveShadows(const VBuf& vbuf);

 private:
  struct IsectInfo {
    int domain_id;
    spray::RTCRayIntersection* isect;
    Ray* ray;
  };

  struct OcclInfo {
    int domain_id;
    Ray* ray;
  };

  spray::MemoryArena* mem_in_;   // rays being processed
  spray::MemoryArena* mem_out_;  // rays being generated
  spray::MemoryArena mem_0_;
  spray::MemoryArena mem_1_;

  const spray::InsituPartition* partition_;
  SceneType* scene_;
  VBuf* vbuf_;
  spray::HdrImage* image_;

  SceneInfo sinfo_;

  ThreadWorkStats work_stats_;  // number of domain blocks to process

  ShaderT shader_;

  Isector<SceneType> isector_;
  spray::RTCRayIntersection rtc_isect_;
  RTCRay rtc_ray_;

  Qvector rqs_;
  Qvector sqs_;

  RayQ rq2_;
  RayQ sq2_;

  RayQ retire_q_;     ///< Retire queue for foreground colors.

  std::queue<IsectInfo> cached_rq_;
  std::queue<IsectInfo> frq2_;
  std::queue<OcclInfo> fsq2_;

  double one_over_num_pixel_samples_;

  int num_domains_;
};

}  // namespace insitu
}  // namespace spray

#define SPRAY_INSITU_TCONTEXT_INL
#include "insitu/insitu_tcontext.inl"
#undef SPRAY_INSITU_TCONTEXT_INL

