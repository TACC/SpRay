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

template <typename SceneT, typename ShaderT>
void SingleThreadTracer<SceneT, ShaderT>::init(const Config &cfg,
                                               const Camera &camera,
                                               SceneT *scene, HdrImage *image) {
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

  if (nranks > 0) comm_recv_.set(&recv_rq_, &recv_sq_);

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

  tile_list_.init(cfg.image_w, cfg.image_h, cfg.pixel_samples, nranks, rank,
                  cfg.maximum_num_screen_space_samples_per_rank);

  CHECK(!tile_list_.empty());

  rqs_.resize(ndomains);
  sqs_.resize(ndomains);
  work_stats_.resize(nranks, cfg.nthreads, ndomains);

  mem_in_ = &mem_0_;
  mem_out_ = &mem_1_;

  vbuf_.resize(tile_list_.getLargestBlockingTile(), cfg.pixel_samples,
               total_num_light_samples);

  isector_.init(ndomains);
}

template <typename SceneT, typename ShaderT>
void SingleThreadTracer<SceneT, ShaderT>::populateRadWorkStats() {
  work_stats_.reset();
  for (int i = 0; i < num_domains_; ++i) {
    int dest = partition_->rank(i);
    if (!rqs_.empty(i)) {
      work_stats_.addNumDomains(dest, 1);
    }
  }
}

template <typename SceneT, typename ShaderT>
void SingleThreadTracer<SceneT, ShaderT>::populateWorkStats() {
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

template <typename SceneT, typename ShaderT>
void SingleThreadTracer<SceneT, ShaderT>::sendRays() {
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

template <typename SceneT, typename ShaderT>
void SingleThreadTracer<SceneT, ShaderT>::send(bool shadow, int domain_id,
                                               int dest, std::queue<Ray *> *q) {
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

template <typename SceneT, typename ShaderT>
void SingleThreadTracer<SceneT, ShaderT>::procLocalQs() {
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

template <typename SceneT, typename ShaderT>
void SingleThreadTracer<SceneT, ShaderT>::procRecvQs() {
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

template <typename SceneT, typename ShaderT>
void SingleThreadTracer<SceneT, ShaderT>::procRads(int id, Ray *rays,
                                                   int64_t count) {
  scene_->load(id, &sinfo_);
  for (auto i = 0; i < count; ++i) {
    auto *ray = &rays[i];
    procRad(id, ray);
  }
}

template <typename SceneT, typename ShaderT>
void SingleThreadTracer<SceneT, ShaderT>::procShads(int id, Ray *rays,
                                                    int64_t count) {
  scene_->load(id, &sinfo_);
  for (auto i = 0; i < count; ++i) {
    auto *ray = &rays[i];
    procShad(id, ray);
  }
}

template <typename SceneT, typename ShaderT>
void SingleThreadTracer<SceneT, ShaderT>::procRad(int id, Ray *ray) {
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

template <typename SceneT, typename ShaderT>
void SingleThreadTracer<SceneT, ShaderT>::procShad(int id, Ray *ray) {
  if (!vbuf_.occluded(ray->samid, ray->light)) {
    bool is_occluded =
        scene_->occluded(sinfo_.rtc_scene, ray->org, ray->dir, &rtc_ray_);

    if (is_occluded) {
      vbuf_.setOBuf(ray->samid, ray->light);
    }
  }
}

template <typename SceneT, typename ShaderT>
void SingleThreadTracer<SceneT, ShaderT>::filterRq2(int id) {
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

template <typename SceneT, typename ShaderT>
void SingleThreadTracer<SceneT, ShaderT>::filterSq2(int id) {
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

template <typename SceneT, typename ShaderT>
void SingleThreadTracer<SceneT, ShaderT>::procFsq2() {
  while (!fsq2_.empty()) {
    auto &info = fsq2_.front();
    auto *ray = info.ray;
    if (vbuf_.correct(ray->samid, ray->t)) {
      if (ray->occluded) {
        vbuf_.setOBuf(ray->samid, ray->light);
      } else {  // to retire q, isect domains exluding its domain
        retire_q_.push(ray);
        isector_.intersect(info.domain_id, scene_, ray, &sqs_);
      }
    }
    fsq2_.pop();
  }
}

template <typename SceneT, typename ShaderT>
void SingleThreadTracer<SceneT, ShaderT>::procFrq2() {
  while (!frq2_.empty()) {
    auto &info = frq2_.front();
    auto *ray = info.ray;
    if (vbuf_.correct(ray->samid, ray->t)) {
      auto *isect = info.isect;

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

template <typename SceneT, typename ShaderT>
void SingleThreadTracer<SceneT, ShaderT>::procCachedRq() {
  while (!cached_rq_.empty()) {
    auto &info = cached_rq_.front();
    auto *ray = info.ray;

    auto *isect = info.isect;

    if (vbuf_.updateTBufOutT(isect->tfar, ray)) {
      shader_(info.domain_id, *ray, *isect, mem_out_, &sq2_, &rq2_, ray_depth_);

      scene_->load(info.domain_id, &sinfo_);
      filterSq2(info.domain_id);
      filterRq2(info.domain_id);
    }

    cached_rq_.pop();
  }
}

template <typename SceneT, typename ShaderT>
void SingleThreadTracer<SceneT, ShaderT>::procRetireQ() {
  while (!retire_q_.empty()) {
    auto *ray = retire_q_.front();
    retire_q_.pop();

    if (!vbuf_.occluded(ray->samid, ray->light)) {
      image_->add(ray->pixid, ray->w, one_over_num_pixel_samples_);
    }
  }
}

template <typename SceneT, typename ShaderT>
void SingleThreadTracer<SceneT, ShaderT>::trace() {
  while (!tile_list_.empty()) {
    tile_list_.front(&blocking_tile_, &stripe_);
    tile_list_.pop();

    vbuf_.resetTBufOut();
    vbuf_.resetOBuf();

    shared_eyes_.num =
        (std::size_t)(stripe_.w * stripe_.h) * num_pixel_samples_;
    if (shared_eyes_.num) {
      shared_eyes_.rays = mem_in_->Alloc<Ray>(shared_eyes_.num);
    }

    if (shared_eyes_.num) {
      glm::vec3 cam_pos = camera_->getPosition();
      if (num_pixel_samples_ > 1) {  // multi samples
        spray::insitu::genMultiSampleEyeRays(
            *camera_, image_w_, cam_pos[0], cam_pos[1], cam_pos[2],
            num_pixel_samples_, blocking_tile_, stripe_, &shared_eyes_);

      } else {  // single sample
        spray::insitu::genSingleSampleEyeRays(
            *camera_, image_w_, cam_pos[0], cam_pos[1], cam_pos[2],
            blocking_tile_, stripe_, &shared_eyes_);
      }

      isector_.intersect(scene_, shared_eyes_, &rqs_);

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
        comm_.run(work_stats_, mem_in_, &comm_recv_);
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
  tile_list_.reset();
}

}  // namespace insitu
}  // namespace spray

