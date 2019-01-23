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

#if !defined(SPRAY_INSITU_TCONTEXT_INL)
#error An implementation of TContext
#endif

namespace spray {
namespace insitu {

template <typename SceneT, typename ShaderT>
void TContext<SceneT, ShaderT>::init(const Config& cfg, int ndomains,
                                     const spray::InsituPartition* partition,
                                     SceneT* scene, VBuf* vbuf,
                                     spray::HdrImage* image) {
  num_domains_ = ndomains;
  partition_ = partition;
  scene_ = scene;
  vbuf_ = vbuf;
  image_ = image;
  bg_color_ = cfg.bg_color;

  rqs_.resize(ndomains);
  sqs_.resize(ndomains);

  shader_.init(cfg, *scene);

  one_over_num_pixel_samples_ = 1.0 / (double)cfg.pixel_samples;
  isector_.init(ndomains);
}

template <typename SceneT, typename ShaderT>
void TContext<SceneT, ShaderT>::populateRadWorkStats() {
#ifdef SPRAY_GLOG_CHECK
  CHECK(work_stats_.empty());
#endif
  for (int id = 0; id < num_domains_; ++id) {
    if (!rqs_.empty(id)) {
      work_stats_.registerRadianceRayBlock(id);
    }
  }
}

template <typename SceneT, typename ShaderT>
void TContext<SceneT, ShaderT>::populateWorkStats() {
#ifdef SPRAY_GLOG_CHECK
  CHECK(work_stats_.empty());
#endif
  for (int id = 0; id < num_domains_; ++id) {
    if (!rqs_.empty(id)) {
      work_stats_.registerRadianceRayBlock(id);
    }
    if (!sqs_.empty(id)) {
      work_stats_.registerShadowRayBlock(id);
    }
  }
  bool has_cached_block = !cached_rq_.empty();
  work_stats_.registerCachedRayBlock(has_cached_block);
}

// template <typename SceneT, typename ShaderT>
// void TContext<SceneT, ShaderT>::populateRadWorkStats() {
//   work_stats_.reset();
//   for (int id = 0; id < num_domains_; ++id) {
//     int dest = partition_->rank(id);
//     if (!rqs_.empty(id)) {
//       work_stats_.addNumDomains(dest, 1);
//     }
//   }
// }

// template <typename SceneT, typename ShaderT>
// void TContext<SceneT, ShaderT>::populateWorkStats(int rank) {
//   work_stats_.reset();
//
//   int n = 0;
//   n += (!cached_rq_.empty());
//
//   work_stats_.addNumDomains(rank, n);
//
//   for (int id = 0; id < num_domains_; ++id) {
//     int dest = partition_->rank(id);
//     n = 0;
//     n += (!rqs_.empty(id));
//     n += (!sqs_.empty(id));
//     if (n) {
//       work_stats_.addNumDomains(dest, n);
//     }
//   }
// }

template <typename SceneT, typename ShaderT>
void TContext<SceneT, ShaderT>::processRays(int id, SceneInfo& sinfo) {
  sinfo_ = sinfo;

  auto* rq = rqs_.getQ(id);

  while (!rq->empty()) {
    auto* ray = rq->front();
    auto* isect = mem_out_->Alloc<spray::RTCRayIntersection>(1, false);

    bool is_hit = scene_->intersect(sinfo_.rtc_scene, sinfo_.cache_block,
                                    ray->org, ray->dir, isect);
    if (is_hit) {
      isect_.ray = ray;
      isect_.isect = isect;
      isects_.push(isect_);
    }
    rq->pop();
  }

  auto* sq = sqs_.getQ(id);

  while (!sq->empty()) {
    auto* ray = sq->front();

    if (!vbuf_->occluded(ray->samid, ray->light)) {
      bool is_occluded =
          scene_->occluded(sinfo_.rtc_scene, ray->org, ray->dir, &rtc_ray_);

      if (is_occluded) {
        occl_.samid = ray->samid;
        occl_.light = ray->light;
        occls_.push(occl_);
      }
    }
    sq->pop();
  }
}

template <typename SceneT, typename ShaderT>
void TContext<SceneT, ShaderT>::updateVisBuf() {
  while (!isects_.empty()) {
    IsectInfo& info = isects_.front();
    if (vbuf_->updateTBufOutT(info.isect->tfar, info.ray)) {
      reduced_isects_.push(info);
    }
    isects_.pop();
  }

  while (!occls_.empty()) {
    OcclInfo& o = occls_.front();
    vbuf_->setOBuf(o.samid, o.light);
    occls_.pop();
  }
}

template <typename SceneT, typename ShaderT>
void TContext<SceneT, ShaderT>::updateTBuf() {
  while (!isects_.empty()) {
    IsectInfo& info = isects_.front();
    if (vbuf_->updateTBufOutT(info.isect->tfar, info.ray)) {
      reduced_isects_.push(info);
    }
    isects_.pop();
  }
}

template <typename SceneT, typename ShaderT>
void TContext<SceneT, ShaderT>::updateOBuf() {
  while (!occls_.empty()) {
    OcclInfo& o = occls_.front();
    vbuf_->setOBuf(o.samid, o.light);
    occls_.pop();
  }
}

template <typename SceneT, typename ShaderT>
void TContext<SceneT, ShaderT>::genRays(int id, int ray_depth) {
  while (!reduced_isects_.empty()) {
    IsectInfo& info = reduced_isects_.front();

    auto* ray = info.ray;
    auto* isect = info.isect;

    if (vbuf_->equalToTbufOut(ray->samid, isect->tfar)) {
      shader_(id, *ray, *isect, mem_out_, &sq2_, &rq2_, ray_depth);
      filterSq2(id);
      filterRq2(id);
    }
    reduced_isects_.pop();
  }
}

template <typename SceneT, typename ShaderT>
void TContext<SceneT, ShaderT>::filterSq2(int id) {
  ocache_item_.domain_id = id;
  RTCScene rtc_scene = sinfo_.rtc_scene;

  while (!sq2_.empty()) {
    auto* ray = sq2_.front();
    sq2_.pop();

    bool is_occluded =
        scene_->occluded(rtc_scene, ray->org, ray->dir, &rtc_ray_);

    if (is_occluded) {
      ray->occluded = 1;
    }
    ocache_item_.ray = ray;
    fsq2_.push(ocache_item_);
  }
}

template <typename SceneT, typename ShaderT>
void TContext<SceneT, ShaderT>::filterRq2(int id) {
  icache_item_.domain_id = id;

  RTCScene rtc_scene = sinfo_.rtc_scene;
  int cache_block = sinfo_.cache_block;

  while (!rq2_.empty()) {
    auto* ray = rq2_.front();
    rq2_.pop();
    auto* isect = mem_out_->Alloc<spray::RTCRayIntersection>(1, false);
    isect->tfar = SPRAY_FLOAT_INF;
    bool is_hit =
        scene_->intersect(rtc_scene, cache_block, ray->org, ray->dir, isect);
    icache_item_.isect = isect;
    icache_item_.ray = ray;

    frq2_.push(icache_item_);
  }
}

template <typename SceneT, typename ShaderT>
void TContext<SceneT, ShaderT>::isectRecvRad(int id, Ray* ray) {
  //
  auto* isect = mem_out_->Alloc<spray::RTCRayIntersection>(1, false);
  //
  bool is_hit = scene_->intersect(sinfo_.rtc_scene, sinfo_.cache_block,
                                  ray->org, ray->dir, isect);

  if (is_hit) {
    isect_.ray = ray;
    isect_.isect = isect;
    isects_.push(isect_);
  }
}

template <typename SceneT, typename ShaderT>
void TContext<SceneT, ShaderT>::occlRecvShad(int id, Ray* ray) {
  int samid = ray->samid;
  int light = ray->light;

  if (!vbuf_->occluded(samid, light)) {
    bool is_occluded =
        scene_->occluded(sinfo_.rtc_scene, ray->org, ray->dir, &rtc_ray_);

    if (is_occluded) {
      occl_.samid = samid;
      occl_.light = light;
      occls_.push(occl_);
    }
  }
}

template <typename SceneT, typename ShaderT>
void TContext<SceneT, ShaderT>::procFsq2() {
  while (!fsq2_.empty()) {
    auto& item = fsq2_.front();
    auto* ray = item.ray;
    int samid = ray->samid;
    if (vbuf_->correct(samid, ray->t)) {
      if (ray->occluded) {  // set obuf
        vbuf_->setOBuf(samid, ray->light);
      } else {  // to retire q, isect domains exluding its domain
        retire_q_.push(ray);
        isector_.intersect(item.domain_id, scene_, ray, &sqs_);
      }
    }
    fsq2_.pop();
  }
}

template <typename SceneT, typename ShaderT>
void TContext<SceneT, ShaderT>::procFrq2() {
  while (!frq2_.empty()) {
    auto& item = frq2_.front();
    auto* ray = item.ray;
    auto* isect = item.isect;
    if (vbuf_->correct(ray->samid, ray->t)) {
#ifdef SPRAY_BACKGROUND_COLOR
      auto* isect = item.isect;

      if (!std::isinf(isect->tfar)) {  // hit
        cached_rq_.push(item);
        isector_.intersect(item.domain_id, isect->tfar, scene_, ray, &rqs_,
                           &bg_retire_q_);
      } else {
        isector_.intersect(item.domain_id, scene_, ray, &rqs_, &bg_retire_q_);
      }
#else
      auto* isect = item.isect;

      if (!std::isinf(isect->tfar)) {  // hit
        cached_rq_.push(item);
        isector_.intersect(item.domain_id, isect->tfar, scene_, ray, &rqs_);
      } else {
        isector_.intersect(item.domain_id, scene_, ray, &rqs_);
      }
#endif
    }
    frq2_.pop();
  }
}

template <typename SceneT, typename ShaderT>
void TContext<SceneT, ShaderT>::updateTBufWithCached() {
  while (!cached_rq_.empty()) {
    auto& item = cached_rq_.front();
    auto* ray = item.ray;

    // supposed to be correct since it is not speculative
    auto* isect = item.isect;

    if (vbuf_->updateTBufOutT(isect->tfar, ray)) {
      reduced_cached_rq_.push(item);
    }
    cached_rq_.pop();
  }
}

template <typename SceneT, typename ShaderT>
void TContext<SceneT, ShaderT>::processCached(int ray_depth) {
  while (!reduced_cached_rq_.empty()) {
    auto& item = reduced_cached_rq_.front();
    auto* ray = item.ray;
    // supposed to be correct since it is not speculative
    auto* isect = item.isect;

    shader_(item.domain_id, *ray, *isect, mem_out_, &sq2_, &rq2_, ray_depth);

    scene_->load(item.domain_id, &sinfo_);
    filterSq2(item.domain_id);
    filterRq2(item.domain_id);

    reduced_cached_rq_.pop();
  }
}

template <typename SceneT, typename ShaderT>
void TContext<SceneT, ShaderT>::procRetireQ() {
  while (!retire_q_.empty()) {
    auto* ray = retire_q_.front();
    retire_q_.pop();
    if (!vbuf_->occluded(ray->samid, ray->light)) {
      image_->add(ray->pixid, ray->w, one_over_num_pixel_samples_);
    }
  }
}
template <typename SceneT, typename ShaderT>
void TContext<SceneT, ShaderT>::retireBackground() {
  glm::vec3 bgcolor;
  while (!bg_retire_q_.empty()) {
    auto* ray = bg_retire_q_.front();
    bg_retire_q_.pop();
    int oflag = ray->occluded;
    if (vbuf_->tbufOutMiss(ray->samid)) {
      bgcolor = glm::vec3(ray->w[0], ray->w[1], ray->w[2]) *
                spray::computeBackGroundColor(ray->dir, bg_color_);
      image_->add(ray->pixid, &bgcolor[0], one_over_num_pixel_samples_);
    }
  }
}

template <typename SceneT, typename ShaderT>
void TContext<SceneT, ShaderT>::sendRays(bool shadow, int id, Ray* rays) {
  std::queue<Ray*>* q;

  if (shadow) {
    q = sqs_.getQ(id);
  } else {
    q = rqs_.getQ(id);
  }

  std::size_t target = 0;

  while (!q->empty()) {
    auto* ray = q->front();
    memcpy(&rays[target], ray, sizeof(Ray));
    ++target;
    q->pop();
  }
}

}  // namespace insitu
}  // namespace spray
