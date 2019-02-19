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

template <typename ShaderT>
void TContext<ShaderT>::init(const Config& cfg, int ndomains,
                             const spray::InsituPartition* partition,
                             SceneType* scene, VBuf* vbuf,
                             spray::HdrImage* image) {
  num_domains_ = ndomains;
  partition_ = partition;
  scene_ = scene;
  vbuf_ = vbuf;
  image_ = image;

  rqs_.resize(ndomains);
  sqs_.resize(ndomains);

  shader_.init(cfg, scene);

  one_over_num_pixel_samples_ = 1.0 / (double)cfg.pixel_samples;
  isector_.init(ndomains);
}

template <typename ShaderT>
void TContext<ShaderT>::populateRadWorkStats() {
#ifdef SPRAY_GLOG_CHECK
  CHECK(work_stats_.empty());
#endif
  for (int id = 0; id < num_domains_; ++id) {
    if (!rqs_.empty(id)) {
      work_stats_.registerRadianceRayBlock(id);
    }
  }
}

template <typename ShaderT>
void TContext<ShaderT>::populateWorkStats() {
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

template <typename ShaderT>
void TContext<ShaderT>::processRays(int rank, int ray_depth) {
  // cached
  while (!cached_rq_.empty()) {
    auto& info = cached_rq_.front();
    auto* ray = info.ray;

    auto* isect = info.isect;

    if (vbuf_->updateTbufOut(isect->tfar, ray)) {
      shader_(info.domain_id, *ray, *isect, mem_out_, &sq2_, &rq2_, ray_depth);

      scene_->load(info.domain_id, &sinfo_);
      filterSq2(info.domain_id);
      filterRq2(info.domain_id);
    }

    cached_rq_.pop();
  }

  // local
  const auto& ids = partition_->getDomains(rank);

  for (auto id : ids) {
    auto* rq = rqs_.getQ(id);
    auto* sq = sqs_.getQ(id);

    bool rq_empty = rq->empty();
    bool sq_empty = sq->empty();

    if (!(rq_empty && sq_empty)) {
      scene_->load(id, &sinfo_);

      while (!rq->empty()) {
        auto* ray = rq->front();
        rq->pop();
        processRadiance(id, ray_depth, ray);
      }

      while (!sq->empty()) {
        auto* ray = sq->front();
        sq->pop();
        processShadow(id, ray);
      }
    }
  }
}

template <typename ShaderT>
void TContext<ShaderT>::processRadiance(int id, int ray_depth, Ray* ray) {
  bool is_hit = scene_->intersect(sinfo_.rtc_scene, sinfo_.cache_block,
                                  ray->org, ray->dir, &rtc_isect_);
  if (is_hit) {
    if (vbuf_->updateTbufOut(rtc_isect_.tfar, ray)) {
      shader_(id, *ray, rtc_isect_, mem_out_, &sq2_, &rq2_, ray_depth);
      filterSq2(id);
      filterRq2(id);
    }
  }
}

template <typename ShaderT>
void TContext<ShaderT>::processShadow(int id, Ray* ray) {
  if (!vbuf_->occluded(ray->samid, ray->light)) {
    bool is_occluded =
        scene_->occluded(sinfo_.rtc_scene, ray->org, ray->dir, &rtc_ray_);

    if (is_occluded) {
      vbuf_->setObuf(ray->samid, ray->light);
    }
  }
}

template <typename ShaderT>
void TContext<ShaderT>::filterSq2(int id) {
  OcclInfo info;
  info.domain_id = id;
  while (!sq2_.empty()) {
    auto* ray = sq2_.front();
    sq2_.pop();

    bool is_occluded =
        scene_->occluded(sinfo_.rtc_scene, ray->org, ray->dir, &rtc_ray_);

    if (is_occluded) {
      ray->occluded = 1;
    }
    info.ray = ray;
    fsq2_.push(info);
  }
}

template <typename ShaderT>
void TContext<ShaderT>::filterRq2(int id) {
  IsectInfo info;
  info.domain_id = id;

  while (!rq2_.empty()) {
    auto* ray = rq2_.front();
    rq2_.pop();
    auto* isect = mem_out_->Alloc<spray::RTCRayIntersection>(1, false);
    isect->tfar = SPRAY_FLOAT_INF;
    bool is_hit = scene_->intersect(sinfo_.rtc_scene, sinfo_.cache_block,
                                    ray->org, ray->dir, isect);
    info.isect = isect;
    info.ray = ray;

    frq2_.push(info);
  }
}

template <typename ShaderT>
void TContext<ShaderT>::resolveRadiances(const VBuf& vbuf) {
  while (!frq2_.empty()) {
    auto& info = frq2_.front();
    auto* ray = info.ray;
    if (vbuf.correct(ray->samid, ray->t)) {
      auto* isect = info.isect;
      if (!std::isinf(isect->tfar)) {  // hit
        cached_rq_.push(info);
        isector_.intersect(info.domain_id, isect->tfar, scene_, ray, &rqs_);
      } else {
        isector_.intersect(info.domain_id, scene_, ray, &rqs_);
      }
    }
    frq2_.pop();
  }
}

template <typename ShaderT>
void TContext<ShaderT>::resolveShadows(const VBuf& vbuf) {
  while (!fsq2_.empty()) {
    auto& info = fsq2_.front();
    auto* ray = info.ray;
    if (vbuf.correct(ray->samid, ray->t)) {
      if (ray->occluded) {
        vbuf_->setObuf(ray->samid, ray->light);
      } else {  // to retire q, isect domains exluding its domain
        retire_q_.push(ray);
        isector_.intersect(info.domain_id, scene_, ray, &sqs_);
      }
    }
    fsq2_.pop();
  }
}

template <typename ShaderT>
void TContext<ShaderT>::retireUntouched() {
  while (!retire_q_.empty()) {
    auto* ray = retire_q_.front();
    retire_q_.pop();
    image_->add(ray->pixid, ray->w, one_over_num_pixel_samples_);
  }
}

template <typename ShaderT>
void TContext<ShaderT>::retireShadows(const VBuf& vbuf) {
  while (!retire_q_.empty()) {
    auto* ray = retire_q_.front();
    retire_q_.pop();
    if (!vbuf.occluded(ray->samid, ray->light)) {
      image_->add(ray->pixid, ray->w, one_over_num_pixel_samples_);
    }
  }
}

template <typename ShaderT>
void TContext<ShaderT>::sendRays(bool shadow, int id, Ray* rays) {
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

template <typename ShaderT>
void TContext<ShaderT>::compositeThreadTbufs(int tid,
                                             std::vector<VBuf>* vbufs) {
  //
  VBuf* dest_buf = &((*vbufs)[0]);

  std::size_t tbuf_size = dest_buf->getTbufSize();
  std::size_t num_threads = vbufs->size();

  std::size_t step = num_threads << 4;

  for (std::size_t i = 0; i < tbuf_size; i += step) {
    std::size_t base = (tid << 4) + i;

    for (std::size_t j = 0; j < 16; ++j) {
      std::size_t index = base + j;

      if (index < tbuf_size) {
        float min = dest_buf->getTbufOut(index);
        for (std::size_t k = 1; k < num_threads; ++k) {
          float t = (*vbufs)[k].getTbufOut(index);
          if (t < min) min = t;
        }
        dest_buf->setTbufOut(index, min);
      }
    }
  }
}

template <typename ShaderT>
void TContext<ShaderT>::compositeThreadObufs(int tid,
                                             std::vector<VBuf>* vbufs) {
  //
  VBuf* dest_buf = &((*vbufs)[0]);

  std::size_t obuf_size = dest_buf->getObufSize();
  std::size_t num_threads = vbufs->size();

  std::size_t step = num_threads << 4;

  for (std::size_t i = 0; i < obuf_size; i += step) {
    std::size_t base = (tid << 4) + i;

    for (std::size_t j = 0; j < 16; ++j) {
      std::size_t index = base + j;

      if (index < obuf_size) {
        uint32_t reduced_ovalue = dest_buf->getObuf(index);
        for (std::size_t k = 1; k < num_threads; ++k) {
          uint32_t ovalue = (*vbufs)[k].getObuf(index);
          reduced_ovalue |= ovalue;
        }
        dest_buf->setObufByIndex(index, reduced_ovalue);
      }
    }
  }
}

}  // namespace insitu
}  // namespace spray
