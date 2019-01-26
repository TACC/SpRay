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

template <typename SceneT, typename ShaderT>
void MultiThreadTracer<SceneT, ShaderT>::init(const Config &cfg,
                                              const Camera &camera,
                                              SceneT *scene, HdrImage *image) {
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

  if (nranks > 0) comm_recv_.set(&recv_rq_, &recv_sq_);

  // shader
  shader_.init(cfg, *scene);

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

  tile_list_.init(cfg.image_w, cfg.image_h, cfg.pixel_samples, nranks, rank,
                  cfg.maximum_num_screen_space_samples_per_rank);

  CHECK(!tile_list_.empty());

  work_stats_.resize(nranks, cfg.nthreads, ndomains);

  vbuf_.resize(tile_list_.getLargestBlockingTile(), cfg.pixel_samples,
               total_num_light_samples);

  tcontexts_.resize(cfg.nthreads);
  for (auto &t : tcontexts_) {
    t.init(cfg, ndomains, partition_, scene_, &vbuf_, image);
  }
  thread_status_.resize(cfg.nthreads);
  scan_.resize(cfg.nthreads);
}

template <typename SceneT, typename ShaderT>
void MultiThreadTracer<SceneT, ShaderT>::populateRadWorkStats(
    TContextType *tcontext) {
  tcontext->populateRadWorkStats();
#pragma omp barrier
#pragma omp single
  {
    work_stats_.reduceRadianceThreadWorkStats<TContextType>(rank_, partition_,
                                                            tcontexts_);
  }
}

template <typename SceneT, typename ShaderT>
void MultiThreadTracer<SceneT, ShaderT>::populateWorkStats(
    TContextType *tcontext) {
  tcontext->populateWorkStats();
#pragma omp barrier
#pragma omp single
  {
    work_stats_.reduceThreadWorkStats<TContextType>(rank_, partition_,
                                                    tcontexts_);
  }
}

template <typename SceneT, typename ShaderT>
void MultiThreadTracer<SceneT, ShaderT>::sendRays(int tid,
                                                  TContextType *tcontext) {
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

template <typename SceneT, typename ShaderT>
void MultiThreadTracer<SceneT, ShaderT>::send(bool shadow, int tid,
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

template <typename SceneT, typename ShaderT>
void MultiThreadTracer<SceneT, ShaderT>::procLocalQs(int tid, int ray_depth,
                                                     TContextType *tcontext) {
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

template <typename SceneT, typename ShaderT>
void MultiThreadTracer<SceneT, ShaderT>::procRecvQs(int ray_depth,
                                                    TContextType *tcontext) {
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

template <typename SceneT, typename ShaderT>
void MultiThreadTracer<SceneT, ShaderT>::procCachedRq(int ray_depth,
                                                      TContextType *tcontext) {
#pragma omp single
  {
    for (auto &t : tcontexts_) {
      t.updateTBufWithCached();
    }
  }
  tcontext->processCached(ray_depth);
}

template <typename SceneT, typename ShaderT>
void MultiThreadTracer<SceneT, ShaderT>::procRecvRads(int ray_depth, int id,
                                                      Ray *rays, int64_t count,
                                                      TContextType *tcontext) {
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

template <typename SceneT, typename ShaderT>
void MultiThreadTracer<SceneT, ShaderT>::procRecvShads(int id, Ray *rays,
                                                       int64_t count,
                                                       TContextType *tcontext) {
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

template <typename SceneT, typename ShaderT>
void MultiThreadTracer<SceneT, ShaderT>::trace() {
#pragma omp parallel
  {
    int nranks = num_ranks_;
    int nbounces = num_bounces_;
    int image_w = image_w_;
    int num_pixel_samples = num_pixel_samples_;

    while (!tile_list_.empty()) {
#pragma omp barrier
#pragma omp master
      {
        tile_list_.front(&blocking_tile_, &stripe_);
        tile_list_.pop();

        vbuf_.resetTBufOut();
        vbuf_.resetOBuf();

        for (auto &tc : tcontexts_) {
          tc.resetMems();
        }

        shared_eyes_.num =
            (std::size_t)(stripe_.w * stripe_.h) * num_pixel_samples;
        if (shared_eyes_.num) {
          shared_eyes_.rays = tcontexts_[0].allocMemIn(shared_eyes_.num);
        }

        done_ = 0;
      }
#pragma omp barrier

      int tid = omp_get_thread_num();
      TContextType *tcontext = &tcontexts_[tid];

      // generate eye rays
      if (shared_eyes_.num) {
        glm::vec3 cam_pos = camera_->getPosition();
        if (num_pixel_samples > 1) {  // multi samples
          spray::insitu::genMultiSampleEyeRays(
              *camera_, image_w, cam_pos[0], cam_pos[1], cam_pos[2],
              num_pixel_samples, blocking_tile_, stripe_, &shared_eyes_);

        } else {  // single sample
          spray::insitu::genSingleSampleEyeRays(
              *camera_, image_w, cam_pos[0], cam_pos[1], cam_pos[2],
              blocking_tile_, stripe_, &shared_eyes_);
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
#pragma omp single
          {
            for (auto &t : tcontexts_) {
              t.procRetireQ();
#ifdef SPRAY_BACKGROUND_COLOR
              t.retireBackground();
#endif
            }
          }
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
            comm_.run(work_stats_, memin, &comm_recv_);
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
              t.procRetireQ();
            }
            vbuf_.resetOBuf();
          }

#ifdef SPRAY_BACKGROUND_COLOR
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
    }  // while (!tile_list_.empty())
  }    // # pragma omp parallel
  tile_list_.reset();
}

template <typename SceneT, typename ShaderT>
void MultiThreadTracer<SceneT, ShaderT>::traceInOmp() {
  while (!tile_list_.empty()) {
#pragma omp barrier
#pragma omp master
    {
      tile_list_.front(&blocking_tile_, &stripe_);
      tile_list_.pop();

      vbuf_.resetTBufOut();
      vbuf_.resetOBuf();

      for (auto &tc : tcontexts_) {
        tc.resetMems();
      }

      shared_eyes_.num =
          (std::size_t)(stripe_.w * stripe_.h) * num_pixel_samples_;
      if (shared_eyes_.num) {
        shared_eyes_.rays = tcontexts_[0].allocMemIn(shared_eyes_.num);
      } else {
        shared_eyes_.rays = nullptr;
      }
    }
#pragma omp barrier

    const int nranks = num_ranks_;
    const int nbounces = num_bounces_;
#pragma omp single
    done_ = 0;

    int tid = omp_get_thread_num();
    TContextType *tcontext = &tcontexts_[tid];

    // generate eye rays
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
#pragma omp single
        {
          for (auto &t : tcontexts_) {
            t.procRetireQ();
#ifdef SPRAY_BACKGROUND_COLOR
            t.retireBackground();
#endif
          }
        }
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
          comm_.run(work_stats_, memin, &comm_recv_);
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
            t.procRetireQ();
          }
          vbuf_.resetOBuf();
        }

#ifdef SPRAY_BACKGROUND_COLOR
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
  }    // while (!tile_list_.empty())
#pragma omp barrier
#pragma omp master
  tile_list_.reset();
#pragma omp barrier
}

}  // namespace insitu
}  // namespace spray

