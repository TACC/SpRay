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

#if !defined(SPRAY_BASELINE_INSITU_TRACER_INL)
#error An implementation of InsituTracer
#endif

namespace spray {
namespace baseline {

template <typename SceneT, typename ScheduleT, typename ShaderT>
void InsituTracer<SceneT, ScheduleT, ShaderT>::init(const Config &cfg,
                                                    const Camera &camera,
                                                    SceneT *scene,
                                                    HdrImage *image) {
  initCommon(cfg, camera, scene, image);
  comm_.init(ndomains_, ray_sched_.getSchedule(), &send_rbuf_, &send_sbuf_,
             &recv_rbuf_, &recv_sbuf_);
}

template <typename SceneT, typename ScheduleT, typename ShaderT>
void InsituTracer<SceneT, ScheduleT, ShaderT>::initCommon(const Config &cfg,
                                                          const Camera &camera,
                                                          SceneT *scene,
                                                          HdrImage *image) {
  int ndomains = static_cast<int>(scene->getNumDomains());
  CHECK_LT(scene->getNumDomains(), std::numeric_limits<int>::max());

  int nranks = mpi::worldSize();
  int rank = mpi::worldRank();

  ndomains_ = ndomains;
  nranks_ = nranks;
  rank_ = rank;

  // pointers
  camera_ = &camera;
  scene_ = scene;
  lights_ = scene->getLights();
  image_ = image;

  // config
  bounces_ = cfg.bounces;
  pixel_samples_ = cfg.pixel_samples;

#ifdef SPRAY_GLOG_CHECK
  LOG(INFO) << "image samples: " << pixel_samples_;
#endif

  CHECK(bounces_ > 0) << "invalid bounces<1";

  // queue statistics
  qstats_.init(ndomains);

  // image scheduler
  img_sched_.init(cfg.image_w, cfg.image_h);

  // locks
  domain_locks_.resize(ndomains);
  for (int i = 0; i < ndomains; ++i) {
    omp_init_lock(&domain_locks_[i]);
  }

  // buffers
  send_rbuf_.alloc(mpi::size());
  send_sbuf_.alloc(mpi::size());

  recv_rbuf_.alloc(1);
  recv_sbuf_.alloc(1);

  // ray scheduler
  ray_sched_.init(ndomains, scene->getInsituPartition());

  // shader
  shader_.init(cfg, scene);

  domain_isectors_.resize(cfg.nthreads);
  for (auto &r : domain_isectors_) {
    r.init(ndomains, scene_);
  }
}

template <typename SceneT, typename ScheduleT, typename ShaderT>
void InsituTracer<SceneT, ScheduleT, ShaderT>::trace() {
  Tile tile = img_sched_.schedule();

  DRay *eyeray_buf = nullptr;
  std::size_t num_eyerays = 0;

  if (tile.isValid()) {
    num_eyerays = tile.getArea() * pixel_samples_;
    std::size_t eyeray_buf_size = num_eyerays;
    eyeray_buf = eyeray_arena_.Alloc<DRay>(eyeray_buf_size, false);
    CHECK_NOTNULL(eyeray_buf);
  }

  int num_bounces = bounces_;
  int nsamples = pixel_samples_;

#pragma omp parallel firstprivate(tile, eyeray_buf, num_eyerays, num_bounces, \
                                  nsamples)
  {
    int ndomains = ndomains_;
    int tid = omp_get_thread_num();
    int nthreads = omp_get_num_threads();
#ifdef SPRAY_GLOG_CHECK
    CHECK_GT(nthreads, 0);
    CHECK_LT(tid, nthreads);
#endif

#pragma omp master
    scan_.resize(nthreads);
#pragma omp barrier

    ArenaQs<DRayQItem> qs;
    ArenaQs<DRayQItem> sqs;
    qs.resize(ndomains);
    sqs.resize(ndomains);

    DRayQ temp_q;

    MemoryArena arena;

    CommitBufferB retire_buf;

    QStats stats;
    stats.init(ndomains);

    // DomainIntersector<SceneT> domain_isector(ndomains, scene_);
    DomainIntersector<SceneT> &domain_isector = domain_isectors_[tid];

    if (num_eyerays) {
      genEyeRays(ndomains, nsamples, tile, eyeray_buf);
      isectEyeDomains(ndomains, num_eyerays, eyeray_buf, &qs, &domain_isector);
    }
#pragma omp barrier

    RTCRay rtc_ray;
    RTCRayIntersection isect;

    while (true) {
#pragma omp barrier
      schedule(ndomains, qs, sqs, &stats);
#pragma omp barrier
      int id;
      bool done = ray_sched_.getSchedule(&id);
      if (done) {
        break;
      }

      const std::vector<RayCount> &sched = ray_sched_.getSchedule();

#pragma omp barrier
      spliceQs(tid, id, ndomains, sched, &qs, &send_rbuf_, &recv_rbuf_);
#pragma omp barrier
      spliceQs(tid, id, ndomains, sched, &sqs, &send_sbuf_, &recv_sbuf_);
#pragma omp barrier

#pragma omp master
      {
        if (mpi::size() > 1) {
          {
            comm_.run();
            comm_.waitForSend();
            comm_.resetMsgSendArena();
          }
        }
      }
#pragma omp barrier

      resetSentQs(ndomains, sched, &qs, &sqs);
#pragma omp barrier

      if (id < ndomains) {
        qs.reset(id);
        sqs.reset(id);

        processRays(tid, id, ndomains, num_bounces, &qs, &sqs, &arena,
                    &domain_isector, &rtc_ray, &isect, &retire_buf, &temp_q);
      }
#pragma omp barrier

#pragma omp critical
      { retire_buf.retire(image_, nsamples); }
      retire_buf.reset();
#pragma omp barrier
#pragma omp master
      comm_.resetQitemRecvArena();
#pragma omp barrier
    }

    arena.Reset();
#pragma omp master
    {
      comm_.resetMsgRecvArena();
      eyeray_arena_.Reset();
    }
#pragma omp barrier
  }
}

template <typename SceneT, typename ScheduleT, typename ShaderT>
void InsituTracer<SceneT, ScheduleT, ShaderT>::schedule(
    int ndomains, const ArenaQs<DRayQItem> &qs, const ArenaQs<DRayQItem> &sqs,
    QStats *stats) {
  for (int i = 0; i < ndomains; ++i) {
    int64_t count = qs.size(i) + sqs.size(i);
    stats->set(i, count);
  }
#pragma omp master
  qstats_.reset();
#pragma omp barrier
  aggregateStats(ndomains, *stats, &qstats_);
#pragma omp barrier
#pragma omp master
  {
    ray_sched_.schedule(qstats_);
    int nranks = mpi::size();
    for (int i = 0; i < nranks; ++i) {
      send_rbuf_.reset(i);
      send_sbuf_.reset(i);
    }
    recv_rbuf_.reset(0);
    recv_sbuf_.reset(0);
  }
}

template <typename SceneT, typename ScheduleT, typename ShaderT>
void InsituTracer<SceneT, ScheduleT, ShaderT>::resetSentQs(
    int ndomains, const std::vector<RayCount> &sched, ArenaQs<DRayQItem> *qs,
    ArenaQs<DRayQItem> *sqs) {
  //
  int nranks = mpi::size();
  int myrank = mpi::rank();
  for (int rank = 0; rank < nranks; ++rank) {
    int id = sched[rank].id;
#ifdef SPRAY_GLOG_CHECK
    CHECK_GE(id, 0);
#endif
    if (rank != myrank && id < ndomains) {
      qs->reset(id);
      sqs->reset(id);
#ifdef SPRAY_GLOG_CHECK
      CHECK(qs->empty(id));
      CHECK(sqs->empty(id));
#endif
    }
  }
}

template <typename SceneT, typename ScheduleT, typename ShaderT>
void InsituTracer<SceneT, ScheduleT, ShaderT>::spliceQs(
    int tid, int id, int ndomains, const std::vector<RayCount> &sched,
    ArenaQs<DRayQItem> *qs, BlockBuffer *send_buf, BlockBuffer *recv_buf) {
  //
  int num_ranks = mpi::size();
#ifdef SPRAY_GLOG_CHECK
  int myrank = mpi::rank();
  CHECK_EQ(id, sched[myrank].id) << myrank;
#endif

  if (num_ranks > 1) {
    spliceQsMultiProcesses(tid, id, ndomains, sched, qs, send_buf, recv_buf);
  } else {
    spliceRecvQs(tid, id, ndomains, qs, recv_buf);
  }
}

template <typename SceneT, typename ScheduleT, typename ShaderT>
void InsituTracer<SceneT, ScheduleT, ShaderT>::spliceRecvQs(
    int tid, int id, int ndomains, ArenaQs<DRayQItem> *qs,
    BlockBuffer *recv_buf) {
  // NOTE: scan_ and recv_buf are shared variables
  // qs: thread-local variable

  if (id < ndomains) {
    // count number of blocks
    int nblocks = qs->getNumValidBlocks(id);
    scan_.set(tid, nblocks);
#pragma omp barrier
#pragma omp master
    {
      scan_.scan();
      int sum = scan_.sum();
      recv_buf->resize(0, sum + 1);                    // +1 for incoming rays
      recv_buf->setLastBlock(0, 0 /*size*/, nullptr);  // initialize last block
#ifdef SPRAY_GLOG_CHECK
      for (int k = 0; k < recv_buf->getNumBlocks(0); ++k) {
        const MemBlock &b = recv_buf->getBlock(0, k);
        CHECK(b.buf == nullptr);
        CHECK_EQ(b.size, 0);
      }
#endif
    }
#pragma omp barrier

    // populate recv block buffer
    if (nblocks) {
#ifdef SPRAY_GLOG_CHECK
      CHECK_GE(scan_.get(tid), nblocks);
#endif
      int idx = scan_.get(tid) - nblocks;
      const std::list<MemBlock> &blocks = qs->getBlocks(id);
      for (const MemBlock &b : blocks) {
        if (b.size) {
#ifdef SPRAY_GLOG_CHECK
          CHECK_NOTNULL(b.buf);
#endif
          recv_buf->set(0, idx, b.size, b.buf);
          ++idx;
        }
      }
    }
  }
}

template <typename SceneT, typename ScheduleT, typename ShaderT>
void InsituTracer<SceneT, ScheduleT, ShaderT>::spliceQsMultiProcesses(
    int tid, int id, int ndomains, const std::vector<RayCount> &sched,
    ArenaQs<DRayQItem> *qs, BlockBuffer *send_buf, BlockBuffer *recv_buf) {
  //
  int num_ranks = mpi::size();
  int myrank = mpi::rank();

  // merge recv queues

#ifdef SPRAY_GLOG_CHECK
  CHECK_EQ(id, sched[myrank].id);
  CHECK_GE(id, 0);
  CHECK(recv_buf->empty(0));
#pragma omp barrier
#endif

  // NOTE: scan_, send_buf and recv_buf are shared variables
  // qs: thread-local variable

  // merge recv queues

  spliceRecvQs(tid, id, ndomains, qs, recv_buf);
#pragma omp barrier

  for (int rank = 0; rank < num_ranks; ++rank) {
#pragma omp barrier
    int domain = sched[rank].id;
    if (rank != myrank && domain < ndomains) {
#ifdef SPRAY_GLOG_CHECK
      CHECK_GE(domain, 0);
#endif
      int nblocks = qs->getNumValidBlocks(domain);
      scan_.set(tid, nblocks);
#pragma omp barrier
#pragma omp master
      {
        scan_.scan();
        int sum = scan_.sum();
        send_buf->resize(rank, sum);
#ifdef SPRAY_GLOG_CHECK
        for (int k = 0; k < send_buf->getNumBlocks(rank); ++k) {
          const MemBlock &b = send_buf->getBlock(rank, k);
          CHECK(b.buf == nullptr);
          CHECK_EQ(b.size, 0);
        }
#endif
      }
#pragma omp barrier

      // populate block buffer
      if (nblocks) {
#ifdef SPRAY_GLOG_CHECK
        CHECK_GE(scan_.get(tid), nblocks);
#endif
        int idx = scan_.get(tid) - nblocks;
        const std::list<MemBlock> &blocks = qs->getBlocks(domain);
        for (const MemBlock &b : blocks) {
          if (b.size) {
#ifdef SPRAY_GLOG_CHECK
            CHECK_NOTNULL(b.buf);
#endif
            send_buf->set(rank, idx, b.size, b.buf);
            ++idx;
          }
        }
      }
    }
  }
}

template <typename SceneT, typename ScheduleT, typename ShaderT>
void InsituTracer<SceneT, ScheduleT, ShaderT>::genEyeRays(int ndomains,
                                                          int nsamples,
                                                          Tile tile,
                                                          DRay *ray_buf) {
  const int image_w = image_->w;

  if (nsamples > 1) {
    const int image_h = image_->h;
    const int total_samples = image_w * image_h * nsamples;

    RandomSampler sampler;

#pragma omp for collapse(3) schedule(static)
    for (int y = tile.y; y < tile.y + tile.h; ++y) {
      for (int x = tile.x; x < tile.x + tile.w; ++x) {
        for (int s = 0; s < nsamples; ++s) {
          int x0 = x - tile.x;
          int y0 = y - tile.y;
          int bufid = nsamples * (y0 * tile.w + x0) + s;
          int pixid = y * image_w + x;
          int samid = nsamples * pixid + s;
          // int samid = nsamples * (y * image_w + x) + s;
          DRay &ray = ray_buf[bufid];
          // RandomSampler_init(sampler, g_nframes * total_samples + samid);
          RandomSampler_init(sampler, samid);
          // RandomSampler_init(sampler, pixid, s);
          float fx = static_cast<float>(x) + RandomSampler_get1D(sampler);
          float fy = static_cast<float>(y) + RandomSampler_get1D(sampler);
          camera_->generateRay(fx, fy, ray.dir);
          DRayUtil::makeEyeRay(pixid, samid, &ray);
        }
      }
    }
  } else {  // if (nsamples > 1) {
            // glm::vec3 cam_org = camera_->getPosition();
            // float* org = &cam_org[0];
#pragma omp for collapse(2) schedule(static)
    for (int y = tile.y; y < tile.y + tile.h; ++y) {
      for (int x = tile.x; x < tile.x + tile.w; ++x) {
        int x0 = x - tile.x;
        int y0 = y - tile.y;
        int bufid = y0 * tile.w + x0;
        int pixid = y * image_w + x;
        DRay &ray = ray_buf[bufid];
        camera_->generateRay(static_cast<float>(x), static_cast<float>(y),
                             ray.dir);
        DRayUtil::makeEyeRay(pixid, pixid, &ray);
      }
    }
  }
}

template <typename SceneT, typename ScheduleT, typename ShaderT>
void InsituTracer<SceneT, ScheduleT, ShaderT>::isectEyeDomains(
    int num_domains, std::size_t num_rays, DRay *ray_buf,
    ArenaQs<DRayQItem> *qs, DomainIntersector<SceneT> *domain_isector) {
  DRayQItem data;
  glm::vec3 cam_org = camera_->getPosition();

#pragma omp for schedule(static, 1)
  for (std::size_t i = 0; i < num_rays; ++i) {
    // setup ray
    DRay &ray = ray_buf[i];

    ray.org[0] = cam_org[0];
    ray.org[1] = cam_org[1];
    ray.org[2] = cam_org[2];

    domain_isector->intersectDomains(&ray);

    auto num_hits = domain_isector->getNumHits();

    if (num_hits) {
      ray.domain_pos = 0;
      if (num_hits > 1) {
        ray.next_tdom = domain_isector->getTnear(1);
      } else {
        ray.next_tdom = SPRAY_FLOAT_INF;
      }

      data.ray = &ray;
      qs->copy(domain_isector->getId(0), &data);

    } else {
      ray.domain_pos = INT_MAX;
      ray.next_tdom = SPRAY_FLOAT_INF;
    }
  }
  //   // private to thread
  //   RTCRayExt eray;
  //
  //   DomainList domains;
  //   DomainHit1 hits[SPRAY_RAY_DOMAIN_LIST_SIZE];
  //
  //   DRayQItem data;
  //
  //   glm::vec3 cam_org = camera_->getPosition();

  // #pragma omp for schedule(static, 1)
  //   for (std::size_t i = 0; i < num_rays; ++i) {
  //     // setup ray
  //     DRay &ray = ray_buf[i];
  //
  //     ray.org[0] = cam_org[0];
  //     ray.org[1] = cam_org[1];
  //     ray.org[2] = cam_org[2];
  //
  //     RTCRayUtil::makeEyeRay(ray, &domains, &eray);
  //
  //     // ray-domain intersection tests
  //     scene_->intersectDomains(eray);

  //     // place ray in hit domains
  //     if (domains.count) {
  //       // sort hit domains
  //       RTCRayUtil::sortDomains(domains, hits);
  //
  //       ray.domain_pos = 0;
  //
  //       if (domains.count > 1) {
  // #ifdef SPRAY_GLOG_CHECK
  //         CHECK_LT(hits[1].id, num_domains);
  //         CHECK_GE(hits[1].id, 0);
  // #endif
  //         ray.next_tdom = hits[1].t;
  //       } else {
  //         ray.next_tdom = SPRAY_FLOAT_INF;
  //       }
  //
  //       data.ray = &ray;
  //       qs->copy(hits[0].id, &data);
  //
  //     } else {
  //       ray.domain_pos = INT_MAX;
  //       ray.next_tdom = SPRAY_FLOAT_INF;
  //     }
  //   }
}

template <typename SceneT, typename ScheduleT, typename ShaderT>
void InsituTracer<SceneT, ScheduleT, ShaderT>::processRays(
    int tid, int id, int ndomains, int nbounces, ArenaQs<DRayQItem> *qs,
    ArenaQs<DRayQItem> *sqs, MemoryArena *arena,
    DomainIntersector<SceneT> *domain_isector, RTCRay *rtc_ray,
    RTCRayIntersection *isect, CommitBufferB *retire_buf, DRayQ *temp_q) {
#ifdef SPRAY_GLOG_CHECK
  CHECK(!(recv_rbuf_.empty(0) && recv_sbuf_.empty(0)));
  if (scene_->insitu()) {
    const InsituPartition &partition = scene_->getInsituPartition();
    CHECK_EQ(partition.rank(id), mpi::rank());
  }
#endif
#pragma omp master
  {
#ifdef SPRAY_TIMING
    spray::tStart(spray::TIMER_LOAD);
#endif
    scene_->load(id);
#ifdef SPRAY_TIMING
    spray::tStop(spray::TIMER_LOAD);
#endif
  }
#pragma omp barrier

  // process input shadow rays

  if (!recv_sbuf_.empty(0)) {
    int nblocks = recv_sbuf_.getNumBlocks(0);
    for (int i = 0; i < nblocks; ++i) {
      const MemBlock &shadow_blk = recv_sbuf_.getBlock(0, i);
      if (shadow_blk.size) {
        DRayQItem *buf = (DRayQItem *)shadow_blk.buf;
        std::size_t num = MemBlock::getSize<DRayQItem>(shadow_blk);
#ifdef SPRAY_PROFILE_COUNTERS
#pragma omp single
        { spray::tAgg(spray::COUNTER_RAYS_TESTED, num); }
#endif
#pragma omp for schedule(static, 1)
        for (std::size_t n = 0; n < num; ++n) {
          DRayQItem &data = buf[n];
#ifdef SPRAY_GLOG_CHECK
          CHECK_EQ(DRayUtil::getShadow(*(data.ray)), 0x00000001) << *(data.ray);
#endif
          processShadow(id, ndomains, data.ray, rtc_ray, sqs, arena,
                        domain_isector, retire_buf);
        }
      }
    }
  }
#pragma omp barrier

  // process input radiance rays

  if (!recv_rbuf_.empty(0)) {
    int nblocks = recv_rbuf_.getNumBlocks(0);
    for (int i = 0; i < nblocks; ++i) {
      const MemBlock &blk = recv_rbuf_.getBlock(0, i);
      if (blk.size) {
        DRayQItem *buf = (DRayQItem *)blk.buf;
        std::size_t num = MemBlock::getSize<DRayQItem>(blk);
#ifdef SPRAY_PROFILE_COUNTERS
#pragma omp single
        { spray::tAgg(spray::COUNTER_RAYS_TESTED, num); }
#endif
#pragma omp for schedule(static, 1)
        for (std::size_t n = 0; n < num; ++n) {
          DRayQItem &data = buf[n];
#ifdef SPRAY_GLOG_CHECK
          CHECK_EQ(DRayUtil::getShadow(*(data.ray)), 0x00000000) << *(data.ray);
          CHECK_GE(data.ray->domid, 0) << *(data.ray);
#endif
          processRay2(id, ndomains, data.ray, rtc_ray, isect, qs, sqs, arena,
                      domain_isector, retire_buf, temp_q);
        }
      }
    }
  }
#pragma omp barrier

  while (!temp_q->empty()) {
    DRay *ray = temp_q->front();
#ifdef SPRAY_GLOG_CHECK
    CHECK_EQ(DRayUtil::getShadow(*ray), 0x00000000) << ray;
    CHECK_GE(ray->domid, 0) << *ray;
#endif
    domain_isector->intersect(INT_MAX, ray, nullptr);
#ifdef SPRAY_GLOG_CHECK
    CHECK(qs->empty(id));
#endif
    processRay2(id, ndomains, ray, rtc_ray, isect, qs, sqs, arena,
                domain_isector, retire_buf, temp_q);

    temp_q->pop();
  }
}

template <typename SceneT, typename ScheduleT, typename ShaderT>
void InsituTracer<SceneT, ScheduleT, ShaderT>::processShadow(
    int id, int ndomains, DRay *ray, RTCRay *rtc_ray, ArenaQs<DRayQItem> *sqs,
    MemoryArena *arena, DomainIntersector<SceneT> *domain_isector,
    CommitBufferB *retire_buf) {
  //
  glm::vec3 pos(ray->org[0], ray->org[1], ray->org[2]);
  glm::vec3 wi(ray->dir[0], ray->dir[1], ray->dir[2]);

  // glm::vec3 Lr(ray->w[0], ray->w[1], ray->w[2]);
  // retire_buf->commit(ray->samid, Lr);

  if (!scene_->occluded(pos, wi, rtc_ray)) {  // occluded

    domain_isector->intersect(id, ray, sqs);
#ifdef SPRAY_GLOG_CHECK
    CHECK(sqs->empty(id));
#endif

    // no more domain
    if (!DRayUtil::hasCurrentDomain(ray)) {
      glm::vec3 Lr(ray->w[0], ray->w[1], ray->w[2]);
      retire_buf->commit(ray->pixid, Lr);
    }
  }
}

template <typename SceneT, typename ScheduleT, typename ShaderT>
void InsituTracer<SceneT, ScheduleT, ShaderT>::aggregateStats(
    int ndomains, const QStats &src_stats, QStats *dest_stats) {
  // combine thread-local stats into process-local stats_
  // only counter values of all depths are aggregated
  for (int i = 0; i < ndomains; ++i) {
    omp_set_lock(&domain_locks_[i]);
    dest_stats->add(i, src_stats);
    omp_unset_lock(&domain_locks_[i]);
  }
}

template <typename SceneT, typename ScheduleT, typename ShaderT>
void InsituTracer<SceneT, ScheduleT, ShaderT>::processRay2(
    int id, int ndomains, DRay *ray, RTCRay *rtc_ray, RTCRayIntersection *isect,
    ArenaQs<DRayQItem> *qs, ArenaQs<DRayQItem> *sqs, MemoryArena *arena,
    DomainIntersector<SceneT> *domain_isector, CommitBufferB *retire_buf,
    DRayQ *temp_q) {
  // in case this ray is in this domain for shading only
  if (DRayUtil::getShade(*ray) == 0x00000001) {
    RTCRayUtil::makeIntersection(*ray, isect);

#ifdef SPRAY_GLOG_CHECK
    CHECK_NE(isect->geomID, RTC_INVALID_GEOMETRY_ID);
#endif

    scene_->updateIntersection(isect);

    shader_(id, ndomains, *ray, rtc_ray, isect, qs, sqs, arena, domain_isector,
            retire_buf, temp_q);

  } else if (scene_->intersect(ray->org, ray->dir, isect)) {
#ifdef SPRAY_GLOG_CHECK
    CHECK(!std::isinf(isect->tfar));
#endif
    if (isect->tfar < ray->t) {  // found closer hit point
      // update domid, t, u, v, primid, geomId, etc.
      DRayUtil::updateIntersection(id, *isect, ray);

      // deteremine if this is last domain

      // last domain or not last domain but found closest hit point
      float next_tdom = ray->next_tdom;
      bool has_next_dom = !std::isinf(next_tdom);
      if (!has_next_dom || (has_next_dom && ray->t < next_tdom)) {
        // this domain has true closest
        shader_(id, ndomains, *ray, rtc_ray, isect, qs, sqs, arena,
                domain_isector, retire_buf, temp_q);

      } else {  // has next domain an hit point farther than next domain
                // (overlapped)
        // move to next domain, this domain might not be true closest
        domain_isector->intersect(id, ray, qs);
#ifdef SPRAY_GLOG_CHECK
        CHECK(qs->empty(id));
#endif
      }
    } else {  // did hit in this domain but not closer (other domain has
              // the closest)
#ifdef SPRAY_GLOG_CHECK
      CHECK(!std::isinf(ray->t));
      CHECK_LT(ray->domid, ndomains) << *ray;
      CHECK_NE(id, ray->domid);
#endif
      // move back to the domain where closest is found
      DRayUtil::setShade(ray);
      DRayQItem data;
      data.ray = ray;
      qs->copy(ray->domid, &data);
#ifdef SPRAY_GLOG_CHECK
      CHECK(qs->empty(id));
#endif
    }
  } else {                                // no intersection in this domain
    if (!DRayUtil::hasNextDomain(ray)) {  // last domain

      if (ray->domid < ndomains) {  // has previous hit
#ifdef SPRAY_GLOG_CHECK
        CHECK_GE(ray->domid, 0);
        CHECK_NE(ray->domid, id);
        CHECK(!std::isinf(ray->t)) << *ray;
#endif
        // move to the domain for shading only
        DRayUtil::setShade(ray);
        DRayQItem data;
        data.ray = ray;
        qs->copy(ray->domid, &data);
#ifdef SPRAY_GLOG_CHECK
        CHECK(qs->empty(id));
#endif
      }
    } else {  // not last domain
              // this domain might not be true closest
      domain_isector->intersect(id, ray, qs);
#ifdef SPRAY_GLOG_CHECK
      CHECK(qs->empty(id));
#endif
    }
  }
}

}  // namespace baseline
}  // namespace spray

