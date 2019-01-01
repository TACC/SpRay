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

#if !defined(SPRAY_INSITU_MULTI_THREAD_TRACER_INL)
#error An implementation of MultiThreadTracer
#endif

namespace spray {
namespace insitu {

template <typename CacheT, typename ShaderT, typename SceneT>
void MultiThreadTracer<CacheT, ShaderT, SceneT>::init(const Config &cfg,
                                                      const Camera &camera,
                                                      SceneT *scene,
                                                      HdrImage *image) {
  int ndomains = static_cast<int>(scene->getNumDomains());
  CHECK_LT(scene->getNumDomains(), std::numeric_limits<int>::max());

  int nranks = mpi::worldSize();
  int rank = mpi::worldRank();

  // pointers
  camera_ = &camera;
  scene_ = scene;
  partition_ = &(scene->getInsituPartition());

  lights_ = scene->getLights();  // copy lights
  image_ = image;

  // settings
  rank_ = rank;
  num_ranks_ = nranks;
  num_domains_ = ndomains;
  num_pixel_samples_ = cfg.pixel_samples;
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

  // tiling
  tiler_.resize(cfg.image_w, cfg.image_h, cfg.num_tiles, cfg.min_tile_size);

  // shader
  shader_.init(cfg, scene);

  // process context
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

  work_stats_.resize(nranks, cfg.nthreads, ndomains);

  vbuf_.resize(image_tile_, cfg.pixel_samples, total_num_light_samples);

  tcontexts_.resize(cfg.nthreads);
  for (auto &t : tcontexts_) {
    t.init(cfg, ndomains, partition_, scene_, &vbuf_, image);
  }
  thread_status_.resize(cfg.nthreads);
  scan_.resize(cfg.nthreads);
}

template <typename CacheT, typename ShaderT, typename SceneT>
void MultiThreadTracer<CacheT, ShaderT, SceneT>::genSingleEyes(
    int image_w, float orgx, float orgy, float orgz, int base_tile_y, Tile tile,
    RayBuf<Ray> *ray_buf) {
  Ray *rays = ray_buf->rays;
#pragma omp for collapse(2) schedule(static)
  for (int y = tile.y; y < tile.y + tile.h; ++y) {
    for (int x = tile.x; x < tile.x + tile.w; ++x) {
      int y0 = y - tile.y;
      int bufid_offset = y0 * tile.w;
      int pixid_offset = y * image_w;
      int x0 = x - tile.x;
      int bufid = bufid_offset + x0;
      int pixid = pixid_offset + x;
      int samid_offset = (tile.y - base_tile_y) * tile.w;
//
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
      // ray->samid = bufid;

      ray->w[0] = 1.f;
      ray->w[1] = 1.f;
      ray->w[2] = 1.f;
      // ray->depth = 0;
      // ray->tdom = 0; //unused
      ray->t = SPRAY_FLOAT_INF;

#ifndef SPRAY_BACKGROUND_COLOR_BLACK
      RayUtil::setOccluded(RayUtil::OFLAG_UNDEFINED, ray);
#endif
    }
  }
}

template <typename CacheT, typename ShaderT, typename SceneT>
void MultiThreadTracer<CacheT, ShaderT, SceneT>::genMultiEyes(
    int image_w, float orgx, float orgy, float orgz, int base_tile_y, Tile tile,
    RayBuf<Ray> *ray_buf) {
  Ray *rays = ray_buf->rays;

  int nsamples = num_pixel_samples_;

#pragma omp for collapse(3) schedule(static, 1)
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

        // common origin
        camera_->generateRay(fx, fy, ray->dir);

        ray->samid = bufid + samid_offset;
        // ray->samid = bufid;

        ray->w[0] = 1.f;
        ray->w[1] = 1.f;
        ray->w[2] = 1.f;
        // ray->depth = 0;
        // ray->tdom = 0; //unused
        ray->t = SPRAY_FLOAT_INF;

#ifndef SPRAY_BACKGROUND_COLOR_BLACK
        RayUtil::setOccluded(RayUtil::OFLAG_UNDEFINED, ray);
#endif
      }
    }
  }
}

template <typename CacheT, typename ShaderT, typename SceneT>
void MultiThreadTracer<CacheT, ShaderT, SceneT>::populateRadWorkStats(
    TContextType *tcontext) {
  tcontext->populateRadWorkStats();
#pragma omp barrier
#pragma omp single
  {
    work_stats_.reduceRadianceThreadWorkStats<TContextType>(rank_, partition_,
                                                            tcontexts_);
  }
}

template <typename CacheT, typename ShaderT, typename SceneT>
void MultiThreadTracer<CacheT, ShaderT, SceneT>::populateWorkStats(
    TContextType *tcontext) {
  tcontext->populateWorkStats();
#pragma omp barrier
#pragma omp single
  {
    work_stats_.reduceThreadWorkStats<TContextType>(rank_, partition_,
                                                    tcontexts_);
  }
}

template <typename CacheT, typename ShaderT, typename SceneT>
void MultiThreadTracer<CacheT, ShaderT, SceneT>::sendRays(
    int tid, TContextType *tcontext) {
  for (int id = 0; id < num_domains_; ++id) {
    int dest = partition_->rank(id);
    if (rank_ != dest) {
      auto num_rads = tcontext->getRqSize(id);
      scan_.set(tid, num_rads);
#pragma omp barrier
#pragma omp single
      { scan_.scan(); }

      if (scan_.sum()) {
        send(false, tid, id, dest, num_rads, tcontext);
      }

#pragma omp barrier

      auto num_shads = tcontext->getSqSize(id);
      scan_.set(tid, num_shads);

#pragma omp barrier
#pragma omp single
      { scan_.scan(); }

      if (scan_.sum()) {
        send(true, tid, id, dest, num_shads, tcontext);
      }
    }
#pragma omp barrier
  }
}

template <typename CacheT, typename ShaderT, typename SceneT>
void MultiThreadTracer<CacheT, ShaderT, SceneT>::send(bool shadow, int tid,
                                                      int domain_id, int dest,
                                                      std::size_t num_rays,
                                                      TContextType *tcontext) {
#pragma omp master
  {
    MsgHeader hout;
    hout.domain_id = domain_id;
    hout.payload_count = scan_.sum();
    int tag = shadow ? WORK_SEND_SHADS : WORK_SEND_RADS;
    send_work_ = new WorkSendMsg<Ray, MsgHeader>(tag, hout, dest);
  }

#pragma omp barrier

  Ray *dest_rays = send_work_->getPayload();
  std::size_t target = scan_.get(tid) - num_rays;

  tcontext->sendRays(shadow, domain_id, &dest_rays[target]);

#pragma omp barrier

#pragma omp master
  { comm_.pushSendQ(send_work_); }
}

template <typename CacheT, typename ShaderT, typename SceneT>
void MultiThreadTracer<CacheT, ShaderT, SceneT>::procLocalQs(
    int tid, int ray_depth, TContextType *tcontext) {
  const auto &ids = partition_->getDomains(rank_);

  for (auto id : ids) {
    bool empty = tcontext->isLocalQsEmpty(id);
    if (empty) {
      thread_status_.clear(tid);
    } else {
      thread_status_.set(tid);
    }

#pragma omp barrier

    bool nonempty = thread_status_.isAnySet();
    if (nonempty) {
#pragma omp single
      { scene_->load(id, &sinfo_); }

      tcontext->processRays(id, sinfo_);

#pragma omp barrier

#pragma omp single
      {
        for (auto &t : tcontexts_) {
          t.updateVisBuf();
        }
      }
      tcontext->genRays(id, ray_depth);
    }
#pragma omp barrier
  }
}

template <typename CacheT, typename ShaderT, typename SceneT>
void MultiThreadTracer<CacheT, ShaderT, SceneT>::procRecvQs(
    int ray_depth, TContextType *tcontext) {
  MsgHeader *header;
  Ray *payload;

  while (1) {
#pragma omp single
    {
      if (!recv_rq_.empty()) {
        recv_message_ = recv_rq_.front();
        recv_rq_.pop();
      } else {
        recv_message_ = nullptr;
      }
    }

    if (recv_message_ == nullptr) {
      break;
    }

    WorkRecvMsg<Ray, MsgHeader>::decode(recv_message_, &header, &payload);
    CHECK_NOTNULL(payload);

    procRecvRads(ray_depth, header->domain_id, payload, header->payload_count,
                 tcontext);
#pragma omp barrier
  }

#pragma omp barrier

  while (1) {
#pragma omp single
    {
      if (!recv_sq_.empty()) {
        recv_message_ = recv_sq_.front();
        recv_sq_.pop();
      } else {
        recv_message_ = nullptr;
      }
    }

    if (recv_message_ == nullptr) {
      break;
    }

    WorkRecvMsg<Ray, MsgHeader>::decode(recv_message_, &header, &payload);
    CHECK_NOTNULL(payload);

    procRecvShads(header->domain_id, payload, header->payload_count, tcontext);
#pragma omp barrier
  }
}

template <typename CacheT, typename ShaderT, typename SceneT>
void MultiThreadTracer<CacheT, ShaderT, SceneT>::procCachedRq(
    int ray_depth, TContextType *tcontext) {
#pragma omp single
  {
    for (auto &t : tcontexts_) {
      t.updateTBufWithCached();
    }
  }
  tcontext->processCached(ray_depth);
}

template <typename CacheT, typename ShaderT, typename SceneT>
void MultiThreadTracer<CacheT, ShaderT, SceneT>::procRecvRads(
    int ray_depth, int id, Ray *rays, int64_t count, TContextType *tcontext) {
#pragma omp single
  { scene_->load(id, &sinfo_); }

  tcontext->setSceneInfo(sinfo_);

#pragma omp for schedule(static, 1)
  for (auto i = 0; i < count; ++i) {
    auto *ray = &rays[i];
    tcontext->isectRecvRad(id, ray);
  }

#pragma omp single
  {
    for (auto &t : tcontexts_) {
      t.updateTBuf();
    }
  }

  tcontext->genRays(id, ray_depth);
}

template <typename CacheT, typename ShaderT, typename SceneT>
void MultiThreadTracer<CacheT, ShaderT, SceneT>::procRecvShads(
    int id, Ray *rays, int64_t count, TContextType *tcontext) {
#pragma omp single
  { scene_->load(id, &sinfo_); }

  tcontext->setSceneInfo(sinfo_);

#pragma omp for schedule(static, 1)
  for (auto i = 0; i < count; ++i) {
    auto *ray = &rays[i];
    tcontext->occlRecvShad(id, ray);
  }

#pragma omp single
  {
    for (auto &t : tcontexts_) {
      t.updateOBuf();
    }
  }
}

template <typename CacheT, typename ShaderT, typename SceneT>
void MultiThreadTracer<CacheT, ShaderT, SceneT>::trace() {
  vbuf_.resetTBufOut();
  vbuf_.resetOBuf();

  RayBuf<Ray> shared_eyes;

  for (auto &tc : tcontexts_) {
    tc.resetMems();
  }

  shared_eyes.num = (std::size_t)(mytile_.w * mytile_.h) * num_pixel_samples_;
  if (shared_eyes.num) {
    shared_eyes.rays = tcontexts_[0].allocMemIn(shared_eyes.num);
  }

  int nranks = num_ranks_;
  int nbounces = num_bounces_;
  int done = 0;

#pragma omp parallel firstprivate(shared_eyes, nranks, nbounces)
  {
    int tid = omp_get_thread_num();
    TContextType *tcontext = &tcontexts_[tid];

    // generate eye rays
    if (shared_eyes.num) {
      glm::vec3 cam_pos = camera_->getPosition();
      if (num_pixel_samples_ > 1) {  // multi samples
        genMultiEyes(image_w_, cam_pos[0], cam_pos[1], cam_pos[2],
                     image_tile_.y, mytile_, &shared_eyes);

      } else {  // single sample
        genSingleEyes(image_w_, cam_pos[0], cam_pos[1], cam_pos[2],
                      image_tile_.y, mytile_, &shared_eyes);
      }
#pragma omp barrier

      // isect domains for eyes on shared variables the eyes buffer
#pragma omp for schedule(static, 1)
      for (std::size_t i = 0; i < shared_eyes.num; ++i) {
        Ray *ray = &shared_eyes.rays[i];
        tcontext->isectDomains(ray);
      }

      populateRadWorkStats(tcontext);
    }

    int ray_depth = 0;

    while (1) {
#pragma omp barrier

#pragma omp master
      {
        work_stats_.reduce();

        if (work_stats_.allDone()) {
          done = 1;
          comm_.waitForSend();
        }
      }

#pragma omp barrier

      if (done) {
        tcontext->procRetireQ(num_pixel_samples_);
        break;
      }

#ifdef SPRAY_GLOG_CHECK
      CHECK_LT(ray_depth, nbounces + 1);
#pragma omp barrier
#endif

      // send rays (transfer WorkSendMsg's to the comm q)
      if (nranks > 1) {
#ifdef SPRAY_GLOG_CHECK
        CHECK(comm_.emptySendQ());
#endif
        sendRays(tid, tcontext);

#pragma omp barrier
#pragma omp master
        {
          comm_.waitForSend();
          auto *memin = tcontext->getMemIn();
          comm_.run(&work_stats_, memin, &recv_rq_, &recv_sq_);
        }
#pragma omp barrier
      }

      procCachedRq(ray_depth, tcontext);
#pragma omp barrier

      procLocalQs(tid, ray_depth, tcontext);
#pragma omp barrier

      if (nranks > 1) {
        procRecvQs(ray_depth, tcontext);
#pragma omp barrier
      }

#pragma omp master
      {
        if (ray_depth < nbounces && nranks > 1) {
          vbuf_.compositeTBuf();
        }

        if (ray_depth > 0 && nranks > 1) {
          vbuf_.compositeOBuf();
        }

        if (ray_depth > 0) {
          for (auto &t : tcontexts_) {
            t.procRetireQ(num_pixel_samples_);
          }
          vbuf_.resetOBuf();
        }

#ifndef SPRAY_BACKGROUND_COLOR_BLACK
        for (auto &t : tcontexts_) {
          t.retireBackground();
        }
#endif

        vbuf_.resetTBufIn();
        vbuf_.swapTBufs();
      }
#pragma omp barrier

      // refer to tbuf input for correctness
      tcontext->processRays2();

#pragma omp barrier

      populateWorkStats(tcontext);

      tcontext->resetAndSwapMems();

      ++ray_depth;

#pragma omp barrier
    }  // while (1)
#ifdef SPRAY_GLOG_CHECK
    tcontext->checkQs();
#endif
  }  // # pragma omp parallel
}

template <typename CacheT, typename ShaderT, typename SceneT>
void MultiThreadTracer<CacheT, ShaderT, SceneT>::traceInOmp() {
#pragma omp master
  {
    vbuf_.resetTBufOut();
    vbuf_.resetOBuf();

    for (auto &tc : tcontexts_) {
      tc.resetMems();
    }

    shared_eyes_.num =
        (std::size_t)(mytile_.w * mytile_.h) * num_pixel_samples_;
    if (shared_eyes_.num) {
      shared_eyes_.rays = tcontexts_[0].allocMemIn(shared_eyes_.num);
    } else {
      shared_eyes_.rays = nullptr;
    }
  }
#pragma omp barrier

  const int nranks = num_ranks_;
  const int nbounces = num_bounces_;
  // int done = 0;
#pragma omp single
  done_ = 0;

  int tid = omp_get_thread_num();
  TContextType *tcontext = &tcontexts_[tid];

  // generate eye rays
  if (shared_eyes_.num) {
    glm::vec3 cam_pos = camera_->getPosition();
    if (num_pixel_samples_ > 1) {  // multi samples
      genMultiEyes(image_w_, cam_pos[0], cam_pos[1], cam_pos[2], image_tile_.y,
                   mytile_, &shared_eyes_);

    } else {  // single sample
      genSingleEyes(image_w_, cam_pos[0], cam_pos[1], cam_pos[2], image_tile_.y,
                    mytile_, &shared_eyes_);
    }
#pragma omp barrier

    // isect domains for eyes on shared variables the eyes buffer
#pragma omp for schedule(static, 1)
    for (std::size_t i = 0; i < shared_eyes_.num; ++i) {
      Ray *ray = &shared_eyes_.rays[i];
      tcontext->isectDomains(ray);
    }

    populateRadWorkStats(tcontext);
  }

  int ray_depth = 0;

  while (1) {
#pragma omp barrier

#pragma omp master
    {
      work_stats_.reduce();

      if (work_stats_.allDone()) {
        done_ = 1;
        comm_.waitForSend();
      }
    }

#pragma omp barrier

    if (done_) {
      tcontext->procRetireQ(num_pixel_samples_);
      break;
    }

#ifdef SPRAY_GLOG_CHECK
    CHECK_LT(ray_depth, nbounces + 1);
#pragma omp barrier
#endif

    // send rays (transfer WorkSendMsg's to the comm q)
    if (nranks > 1) {
#ifdef SPRAY_GLOG_CHECK
      CHECK(comm_.emptySendQ());
#endif
      sendRays(tid, tcontext);

#pragma omp barrier
#pragma omp master
      {
        comm_.waitForSend();
        auto *memin = tcontext->getMemIn();
        comm_.run(&work_stats_, memin, &recv_rq_, &recv_sq_);
      }
#pragma omp barrier
    }

    procCachedRq(ray_depth, tcontext);
#pragma omp barrier

    procLocalQs(tid, ray_depth, tcontext);
#pragma omp barrier

    if (nranks > 1) {
      procRecvQs(ray_depth, tcontext);
#pragma omp barrier
    }

#pragma omp master
    {
      if (ray_depth < nbounces && nranks > 1) {
        vbuf_.compositeTBuf();
      }

      if (ray_depth > 0 && nranks > 1) {
        vbuf_.compositeOBuf();
      }

      if (ray_depth > 0) {
        for (auto &t : tcontexts_) {
          t.procRetireQ(num_pixel_samples_);
        }
        vbuf_.resetOBuf();
      }

#ifndef SPRAY_BACKGROUND_COLOR_BLACK
      for (auto &t : tcontexts_) {
        t.retireBackground();
      }
#endif

      vbuf_.resetTBufIn();
      vbuf_.swapTBufs();
    }
#pragma omp barrier

    // refer to tbuf input for correctness
    tcontext->processRays2();

#pragma omp barrier

    populateWorkStats(tcontext);

    tcontext->resetAndSwapMems();

    ++ray_depth;

#pragma omp barrier
  }  // while (1)
}

}  // namespace insitu
}  // namespace spray

