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

#include <omp.h>
#include <cstring>
#include <queue>
#include <utility>

// clang-format off
#include <embree2/rtcore_ray.h>
#include <embree2/rtcore_geometry.h>
// clang-format on

#include "glm/glm.hpp"
#include "glog/logging.h"
#include "pbrt/memory.h"

#include "ooc/ooc_domain_stats.h"
#include "ooc/ooc_isector.h"
#include "ooc/ooc_ray.h"
#include "ooc/ooc_vbuf.h"
#include "render/qvector.h"
#include "render/rays.h"
#include "utils/scan.h"
#include "display/image.h"
#include "utils/scan.h"

namespace spray {
namespace ooc {

template <typename CacheT, typename ShaderT, typename SceneT>
class TContext {
 public:
  TContext() {
    mem_in_ = &mem_0_;
    mem_out_ = &mem_1_;

    sqs_in_ = &sqs_0_;
    sqs_out_ = &sqs_1_;

    commit_q_ = &commit_retire_q0_;
    retire_q_ = &commit_retire_q1_;
  }

 public:
  void resize(int ndomains, int num_pixel_samples, const Tile& tile,
              spray::HdrImage* image, int num_bounces);

 private:
  int num_domains_;
  int num_pixel_samples_;
  int num_bounces_;
  double one_over_num_pixel_samples_;

 public:
  void resetMems() {
    mem_in_->Reset();
    mem_out_->Reset();
  }

 private:
  VBuf vbuf_;
  spray::HdrImage* image_;

 public:
  void resetMemIn() { mem_in_->Reset(); }
  void swapMems() { std::swap(mem_in_, mem_out_); }

 public:
  template <typename T>
  T* allocMemIn(std::size_t size) {
    T* mem = mem_in_->Alloc<T>(size, false);
    CHECK_NOTNULL(mem);
    return mem;
  }

  template <typename T>
  T* allocMemOut(std::size_t size) {
    T* mem = mem_out_->Alloc<T>(size, false);
    CHECK_NOTNULL(mem);
    return mem;
  }

  spray::MemoryArena* getMemOut() { return mem_out_; }

 private:
  template <typename T>
  T* allocOut(std::size_t size) {
    T* mem = mem_out_->Alloc<T>(size, false);
    return mem;
  }

 private:
  spray::MemoryArena* mem_in_;   // rays being processed
  spray::MemoryArena* mem_out_;  // rays being generated
  spray::MemoryArena mem_0_;
  spray::MemoryArena mem_1_;

 public:
  void enqRad(SceneT* scene, Ray* ray) {
#ifdef SPRAY_BACKGROUND_COLOR_BLACK
    isector_.intersect(num_domains_, scene, ray, &rqs_, &rstats_);
#else
    isector_.intersect(num_domains_, scene, ray, &rqs_, &rstats_,
                       &bg_retire_q_);
#endif
  }

  spray::RTCRayIntersection& getRTCIsect() { return rtc_isect_; }
  RTCRay& getRTCRay() { return rtc_ray_; }

 private:
  Isector<CacheT, SceneT> isector_;
  spray::RTCRayIntersection rtc_isect_;
  RTCRay rtc_ray_;

 private:
  std::queue<Ray*> sq2_;
  std::queue<Ray*> rq2_;
  std::queue<Ray*> pending_q_;

  std::queue<Ray*> frq_;
  std::queue<Ray*> fsq_in_;
  std::queue<Ray*> fsq_out_;

  std::queue<Ray*> commit_retire_q0_;
  std::queue<Ray*> commit_retire_q1_;

  std::queue<Ray*>* commit_q_;
  std::queue<Ray*>* retire_q_;

  spray::QVector<RayData> rqs_;
  spray::QVector<RayData> sqs_0_;
  spray::QVector<RayData> sqs_1_;

  spray::QVector<RayData>* sqs_in_;
  spray::QVector<RayData>* sqs_out_;

  std::queue<Ray*> bg_retire_q_;  ///< Retire queue for background colors.

 public:
  bool allFilterQsEmpty() const {
    return (frq_.empty() && fsq_in_.empty() && fsq_out_.empty());
  }

  void retire();
  void retireBackground();

  bool sqsInEmpty() const { return sqs_in_->empty(); }
  bool rqsEmpty() const { return rqs_.empty(); }
  bool retireQEmpty() const { return retire_q_->empty(); }
  bool backgroundQEmpty() const { return bg_retire_q_.empty(); }
  bool commitQEmpty() const { return commit_q_->empty(); }
  bool pendingQEmpty() const { return pending_q_.empty(); }

  void swapQs() {
    std::swap(commit_q_, retire_q_);
    std::swap(sqs_in_, sqs_out_);
  }

 public:
  void filterQs(int id) {
#ifdef SPRAY_GLOG_CHECK
    CHECK(frq_.empty());
    CHECK(fsq_in_.empty());
    CHECK(fsq_out_.empty()) << fsq_out_.size();
#endif
    filterRqs(id);
    filterSqs(id, sqs_in_, &fsq_in_);
    filterSqs(id, sqs_out_, &fsq_out_);
  }

 private:
  void filterRqs(int id);
  void filterSqs(int id, QVector<RayData>* sqs, std::queue<Ray*>* fsq);

 public:
  void procFilterQs(int id, SceneT* scene, SceneInfo& sinfo, ShaderT& shader,
                    int ray_depth) {
    procRads(id, scene, sinfo, shader, ray_depth);
    procShads(scene, sinfo, &fsq_in_, retire_q_);
    procShads(scene, sinfo, &fsq_out_, commit_q_);
  }

 private:
  void procRads(int id, SceneT* scene, SceneInfo& sinfo, ShaderT& shader,
                int ray_depth);

  void procShads2(int id, SceneT* scene, SceneInfo& sinfo);

  void procRads2(SceneT* scene, SceneInfo& sinfo);

  void procShads(SceneT* scene, SceneInfo& sinfo, std::queue<Ray*>* qin,
                 std::queue<Ray*>* qout);

 public:
  void procPendingQ(SceneT* scene) {
    while (!pending_q_.empty()) {
      Ray* r = pending_q_.front();
      pending_q_.pop();
      if (vbuf_.correct(*r)) {
        r->depth = 0;
        r->history[0] = SPRAY_FLOAT_INF;
#ifdef SPRAY_BACKGROUND_COLOR_BLACK
        isector_.intersect(num_domains_, scene, r, &rqs_, &rstats_);
#else
        isector_.intersect(num_domains_, scene, r, &rqs_, &rstats_,
                           &bg_retire_q_);
#endif
      }
    }
  }

  void resetVBuf() { vbuf_.reset(); }

 private:
  DomainStats rstats_;

 public:
  void resetRstats() { rstats_.reset(); }
  const DomainStats& getRstats() const { return rstats_; }
};

}  // namespace ooc
}  // namespace spray

#define SPRAY_OOC_TCONTEXT_INL
#include "ooc/ooc_tcontext.inl"
#undef SPRAY_OOC_TCONTEXT_INL

