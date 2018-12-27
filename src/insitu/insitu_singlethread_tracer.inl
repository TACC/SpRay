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

#if !defined(SPRAY_INSITU_SINGLE_THRAED_TRACER_INL)
#error An implementation of SingleThreadTracer
#endif

namespace spray {
namespace insitu {

template <typename CacheT, typename ShaderT, typename SceneT>
void SingleThreadTracer<CacheT, ShaderT, SceneT>::init(const Config &cfg,
                                                       const Camera &camera,
                                                       SceneT *scene,
                                                       HdrImage *image) {
  int ndomains = static_cast<int>(scene->getNumDomains());
  CHECK_LT(scene->getNumDomains(), std::numeric_limits<int>::max());

  int nranks = mpi::worldSize();
  int rank = mpi::worldRank();

  camera_ = &camera;
  scene_ = scene;
  partition_ = &(scene->getInsituPartition());

  lights_ = scene->getLights();  // copy lights
  image_ = image;

  rank_ = rank;
  num_ranks_ = nranks;
  num_domains_ = ndomains;
  num_pixel_samples_ = cfg.pixel_samples;
  one_over_num_pixel_samples_ = 1.0 / (double)num_pixel_samples_;
  num_bounces_ = cfg.bounces;
  num_threads_ = cfg.nthreads;
  image_w_ = cfg.image_w;
  image_h_ = cfg.image_h;

  CHECK_GT(rank_, -1);
  CHECK_GT(num_ranks_, 0);
  CHECK_GT(num_domains_, 0);
  CHECK_GT(num_pixel_samples_, 0);
  CHECK_GT(num_bounces_, 0);
  CHECK_GT(image_w_, 0);
  CHECK_GT(image_h_, 0);

  tiler_.resize(cfg.image_w, cfg.image_h, cfg.num_tiles, cfg.min_tile_size);

  shader_.init(cfg, scene);

  int total_num_light_samples;
  if (shader_.isAo()) {
    num_lights_ = cfg.ao_samples;
    total_num_light_samples = num_lights_;
  } else {
    CHECK_GT(lights_.size(), 0);
    num_lights_ = lights_.size();
    std::size_t num_lights = lights_.size();

    total_num_light_samples = 0;
    for (std::size_t i = 0; i < lights_.size(); ++i) {
      if (lights_[i]->isAreaLight()) {
        total_num_light_samples += cfg.ao_samples;
      } else {
        ++total_num_light_samples;
      }
    }
  }

  image_tile_.x = 0;
  image_tile_.y = 0;
  image_tile_.w = cfg.image_w;
  image_tile_.h = cfg.image_h;

  mytile_ = RankStriper::make(mpi::size(), mpi::rank(), image_tile_);

  int64_t total_num_samples =
      (int64_t)mytile_.w * mytile_.h * cfg.pixel_samples;
  CHECK_LT(total_num_samples, INT_MAX);

  rqs_.resize(ndomains);
  sqs_.resize(ndomains);
  work_stats_.resize(nranks, cfg.nthreads, ndomains);

  mem_in_ = &mem_0_;
  mem_out_ = &mem_1_;

  vbuf_.resize(image_tile_, cfg.pixel_samples, total_num_light_samples);
}

template <typename CacheT, typename ShaderT, typename SceneT>
void SingleThreadTracer<CacheT, ShaderT, SceneT>::genSingleEyes(
    int image_w, float orgx, float orgy, float orgz, int base_tile_y, Tile tile,
    RayBuf<Ray> *ray_buf) {
  Ray *rays = ray_buf->rays;
  for (int y = tile.y; y < tile.y + tile.h; ++y) {
    for (int x = tile.x; x < tile.x + tile.w; ++x) {
      int y0 = y - tile.y;
      int bufid_offset = y0 * tile.w;
      int pixid_offset = y * image_w;
      int x0 = x - tile.x;
      int bufid = bufid_offset + x0;
      int pixid = pixid_offset + x;
      int samid_offset = (tile.y - base_tile_y) * tile.w;
#ifdef SPRAY_GLOG_CHECK
      CHECK_LT(bufid, ray_buf->num);
#endif
      auto *ray = &rays[bufid];
      //
      ray->org[0] = orgx;
      ray->org[1] = orgy;
      ray->org[2] = orgz;

      ray->pixid = pixid;

      camera_->generateRay((float)x, (float)y, ray->dir);

      ray->samid = bufid + samid_offset;

      ray->w[0] = 1.f;
      ray->w[1] = 1.f;
      ray->w[2] = 1.f;

      ray->t = SPRAY_FLOAT_INF;
    }
  }
}

template <typename CacheT, typename ShaderT, typename SceneT>
void SingleThreadTracer<CacheT, ShaderT, SceneT>::genMultiEyes(
    int image_w, float orgx, float orgy, float orgz, int base_tile_y, Tile tile,
    RayBuf<Ray> *ray_buf) {
  Ray *rays = ray_buf->rays;

  int nsamples = num_pixel_samples_;

  for (int y = tile.y; y < tile.y + tile.h; ++y) {
    for (int x = tile.x; x < tile.x + tile.w; ++x) {
      for (int s = 0; s < nsamples; ++s) {
        int x0 = x - tile.x;
        int y0 = y - tile.y;
        int bufid = nsamples * (y0 * tile.w + x0) + s;
        int pixid = y * image_w + x;
        int samid_offset = (tile.y - base_tile_y) * tile.w * nsamples;
#ifdef SPRAY_GLOG_CHECK
        CHECK_LT(bufid, ray_buf->num);
#endif
        Ray *ray = &rays[bufid];
        //
        ray->org[0] = orgx;
        ray->org[1] = orgy;
        ray->org[2] = orgz;

        ray->pixid = pixid;

        RandomSampler sampler;
        RandomSampler_init(sampler, bufid + samid_offset);

        float fx = (float)(x) + RandomSampler_get1D(sampler);
        float fy = (float)(y) + RandomSampler_get1D(sampler);

        camera_->generateRay(fx, fy, ray->dir);

        ray->samid = bufid + samid_offset;

        ray->w[0] = 1.f;
        ray->w[1] = 1.f;
        ray->w[2] = 1.f;

        ray->t = SPRAY_FLOAT_INF;
      }
    }
  }
}

template <typename CacheT, typename ShaderT, typename SceneT>
void SingleThreadTracer<CacheT, ShaderT, SceneT>::populateRadWorkStats() {
  work_stats_.reset();
  for (int i = 0; i < num_domains_; ++i) {
    int dest = partition_->rank(i);
    if (!rqs_.empty(i)) {
      work_stats_.addNumDomains(dest, 1);
    }
  }
}

template <typename CacheT, typename ShaderT, typename SceneT>
void SingleThreadTracer<CacheT, ShaderT, SceneT>::populateWorkStats() {
  work_stats_.reset();

  int n = 0;
  n += (!cached_rq_.empty());

  work_stats_.addNumDomains(rank_, n);

  for (int i = 0; i < num_domains_; ++i) {
    int dest = partition_->rank(i);
    n = 0;
    n += (!rqs_.empty(i));
    n += (!sqs_.empty(i));
    if (n) {
      work_stats_.addNumDomains(dest, n);
    }
  }
}

template <typename CacheT, typename ShaderT, typename SceneT>
void SingleThreadTracer<CacheT, ShaderT, SceneT>::sendRays() {
  for (int i = 0; i < num_domains_; ++i) {
    int dest = partition_->rank(i);
    if (rank_ != dest) {
      // rads
      auto *rq = rqs_.getQ(i);
      if (!rq->empty()) {
        send(false, i, dest, rq);
      }
      // shads
      auto *sq = sqs_.getQ(i);
      if (!sq->empty()) {
        send(true, i, dest, sq);
      }
    }
  }
}

template <typename CacheT, typename ShaderT, typename SceneT>
void SingleThreadTracer<CacheT, ShaderT, SceneT>::send(bool shadow,
                                                       int domain_id, int dest,
                                                       std::queue<Ray *> *q) {
  MsgHeader hout;
  hout.domain_id = domain_id;
  hout.payload_count = q->size();

  int tag = shadow ? WORK_SEND_SHADS : WORK_SEND_RADS;

  WorkSendMsg<Ray, MsgHeader> *send_work =
      new WorkSendMsg<Ray, MsgHeader>(tag, hout, dest);

  Ray *dest_rays = send_work->getPayload();

  std::size_t target = 0;

  while (!q->empty()) {
    auto *ray = q->front();
    q->pop();

    memcpy(&dest_rays[target], ray, sizeof(Ray));
    ++target;
  }

  comm_.pushSendQ(send_work);
}

template <typename CacheT, typename ShaderT, typename SceneT>
void SingleThreadTracer<CacheT, ShaderT, SceneT>::procLocalQs() {
  const auto &ids = partition_->getDomains(rank_);

  for (auto id : ids) {
    auto *rq = rqs_.getQ(id);
    auto *sq = sqs_.getQ(id);

    bool rq_empty = rq->empty();
    bool sq_empty = sq->empty();

    if (!(rq_empty && sq_empty)) {
      scene_->load(id, &sinfo_);

      while (!rq->empty()) {
        auto *ray = rq->front();
        rq->pop();
        procRad(id, ray);
      }

      while (!sq->empty()) {
        auto *ray = sq->front();
        sq->pop();
        procShad(id, ray);
      }
    }
  }
}

template <typename CacheT, typename ShaderT, typename SceneT>
void SingleThreadTracer<CacheT, ShaderT, SceneT>::procRecvQs() {
  MsgHeader *header;
  Ray *payload;

  while (!recv_rq_.empty()) {
    auto *msg = recv_rq_.front();
    recv_rq_.pop();
    WorkRecvMsg<Ray, MsgHeader>::decode(msg, &header, &payload);
    CHECK_NOTNULL(payload);
    procRads(header->domain_id, payload, header->payload_count);
  }

  while (!recv_sq_.empty()) {
    auto *msg = recv_sq_.front();
    recv_sq_.pop();
    WorkRecvMsg<Ray, MsgHeader>::decode(msg, &header, &payload);
    CHECK_NOTNULL(payload);
    procShads(header->domain_id, payload, header->payload_count);
  }
}

template <typename CacheT, typename ShaderT, typename SceneT>
void SingleThreadTracer<CacheT, ShaderT, SceneT>::procRads(int id, Ray *rays,
                                                           int64_t count) {
  scene_->load(id, &sinfo_);
  for (auto i = 0; i < count; ++i) {
    auto *ray = &rays[i];
    procRad(id, ray);
  }
}

template <typename CacheT, typename ShaderT, typename SceneT>
void SingleThreadTracer<CacheT, ShaderT, SceneT>::procShads(int id, Ray *rays,
                                                            int64_t count) {
  scene_->load(id, &sinfo_);
  for (auto i = 0; i < count; ++i) {
    auto *ray = &rays[i];
    procShad(id, ray);
  }
}

template <typename CacheT, typename ShaderT, typename SceneT>
void SingleThreadTracer<CacheT, ShaderT, SceneT>::procRad(int id, Ray *ray) {
  bool is_hit = scene_->intersect(sinfo_.rtc_scene, sinfo_.cache_block,
                                  ray->org, ray->dir, &rtc_isect_);

  if (is_hit) {
    if (vbuf_.updateTBufOutT(rtc_isect_.tfar, ray)) {
      shader_(id, *ray, rtc_isect_, mem_out_, &sq2_, &rq2_, ray_depth_);
      filterSq2(id);
      filterRq2(id);
    }
  }
}

template <typename CacheT, typename ShaderT, typename SceneT>
void SingleThreadTracer<CacheT, ShaderT, SceneT>::procShad(int id, Ray *ray) {
  if (!vbuf_.occluded(ray->samid, ray->light)) {
    bool is_occluded =
        scene_->occluded(sinfo_.rtc_scene, ray->org, ray->dir, &rtc_ray_);

    if (is_occluded) {
      vbuf_.setOBuf(ray->samid, ray->light);
    }
  }
}

template <typename CacheT, typename ShaderT, typename SceneT>
void SingleThreadTracer<CacheT, ShaderT, SceneT>::filterRq2(int id) {
  IsectInfo info;
  info.domain_id = id;

  while (!rq2_.empty()) {
    auto *ray = rq2_.front();
    rq2_.pop();
    auto *isect = mem_out_->Alloc<spray::RTCRayIntersection>(1, false);
    isect->tfar = SPRAY_FLOAT_INF;
    bool is_hit = scene_->intersect(sinfo_.rtc_scene, sinfo_.cache_block,
                                    ray->org, ray->dir, isect);
    info.isect = isect;
    info.ray = ray;

    frq2_.push(info);
  }
}

template <typename CacheT, typename ShaderT, typename SceneT>
void SingleThreadTracer<CacheT, ShaderT, SceneT>::filterSq2(int id) {
  OcclInfo info;
  info.domain_id = id;
  while (!sq2_.empty()) {
    auto *ray = sq2_.front();
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

template <typename CacheT, typename ShaderT, typename SceneT>
void SingleThreadTracer<CacheT, ShaderT, SceneT>::procFsq2() {
  while (!fsq2_.empty()) {
    auto &info = fsq2_.front();
    auto *ray = info.ray;
    if (vbuf_.correct(ray->samid, ray->t)) {
      if (ray->occluded) {
        vbuf_.setOBuf(ray->samid, ray->light);
      } else {  // to retire q, isect domains exluding its domain
        retire_q_.push(ray);
        isector_.intersect(info.domain_id, num_domains_, scene_, ray, &sqs_);
      }
    }
    fsq2_.pop();
  }
}

template <typename CacheT, typename ShaderT, typename SceneT>
void SingleThreadTracer<CacheT, ShaderT, SceneT>::procFrq2() {
  //
  // DEBUG
  // std::cout << "frq2.size()" << frq2_.size() << "\n";
  //
  while (!frq2_.empty()) {
    auto &info = frq2_.front();
    auto *ray = info.ray;
    if (vbuf_.correct(ray->samid, ray->t)) {
      auto *isect = info.isect;

      if (!std::isinf(isect->tfar)) {  // hit
        cached_rq_.push(info);
        isector_.intersect(info.domain_id, isect->tfar, num_domains_, scene_,
                           ray, &rqs_);
      } else {
        isector_.intersect(info.domain_id, num_domains_, scene_, ray, &rqs_);
      }
    }
    frq2_.pop();
  }
}

template <typename CacheT, typename ShaderT, typename SceneT>
void SingleThreadTracer<CacheT, ShaderT, SceneT>::procCachedRq() {
  while (!cached_rq_.empty()) {
    auto &info = cached_rq_.front();
    auto *ray = info.ray;

    auto *isect = info.isect;

    if (vbuf_.updateTBufOutT(isect->tfar, ray)) {
      shader_(info.domain_id, *ray, *isect, mem_out_, &sq2_, &rq2_, ray_depth_);
      filterSq2(info.domain_id);
      filterRq2(info.domain_id);
    }

    cached_rq_.pop();
  }
}

template <typename CacheT, typename ShaderT, typename SceneT>
void SingleThreadTracer<CacheT, ShaderT, SceneT>::procRetireQ() {
  while (!retire_q_.empty()) {
    auto *ray = retire_q_.front();
    retire_q_.pop();

    if (!vbuf_.occluded(ray->samid, ray->light)) {
      image_->add(ray->pixid, ray->w, one_over_num_pixel_samples_);
    }
  }
}

template <typename CacheT, typename ShaderT, typename SceneT>
void SingleThreadTracer<CacheT, ShaderT, SceneT>::trace() {
  vbuf_.resetTBufOut();
  vbuf_.resetOBuf();

  RayBuf<Ray> shared_eyes;

  shared_eyes.num = (std::size_t)(mytile_.w * mytile_.h) * num_pixel_samples_;
  if (shared_eyes.num) {
    shared_eyes.rays = mem_in_->Alloc<Ray>(shared_eyes.num);
  }

  if (shared_eyes.num) {
    glm::vec3 cam_pos = camera_->getPosition();
    if (num_pixel_samples_ > 1) {  // multi samples
      genMultiEyes(image_w_, cam_pos[0], cam_pos[1], cam_pos[2], image_tile_.y,
                   mytile_, &shared_eyes);

    } else {  // single sample
      genSingleEyes(image_w_, cam_pos[0], cam_pos[1], cam_pos[2], image_tile_.y,
                    mytile_, &shared_eyes);
    }

    isector_.intersect(num_domains_, scene_, shared_eyes, &rqs_);

    populateRadWorkStats();
  }

  ray_depth_ = 0;

  while (1) {
    work_stats_.reduce();

    if (work_stats_.allDone()) {
      procRetireQ();
      comm_.waitForSend();
      break;
    }

    if (num_ranks_ > 1) {
#ifdef SPRAY_GLOG_CHECK
      CHECK(comm_.emptySendQ());
#endif
      sendRays();
      comm_.waitForSend();
      comm_.run(&work_stats_, mem_in_, &recv_rq_, &recv_sq_);
    }

    procCachedRq();
    procLocalQs();
    procRecvQs();

    if (ray_depth_ < num_bounces_ && num_ranks_ > 1) {
      vbuf_.compositeTBuf();
    }

    if (ray_depth_ > 0 && num_ranks_ > 1) {
      vbuf_.compositeOBuf();
    }

    if (ray_depth_ > 0) {
      procRetireQ();
      vbuf_.resetOBuf();
    }

    vbuf_.resetTBufIn();
    vbuf_.swapTBufs();

    procFrq2();
    procFsq2();

    populateWorkStats();

    mem_in_->Reset();
    std::swap(mem_in_, mem_out_);

    ++ray_depth_;
  }
}

}  // namespace insitu
}  // namespace spray

