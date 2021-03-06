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

#if !defined(SPRAY_OOC_TRACER_INL)
#error An implementation of Tracer
#endif

namespace spray {
namespace ooc {

template <typename ShaderT>
void Tracer<ShaderT>::init(const Config &cfg, const Camera &camera,
                           SceneType *scene, HdrImage *image) {
  int ndomains = static_cast<int>(scene->getNumDomains());
  CHECK_LT(scene->getNumDomains(), std::numeric_limits<int>::max());

  int nranks = mpi::worldSize();
  int rank = mpi::worldRank();

  camera_ = &camera;
  scene_ = scene;

  lights_ = scene->getLights();
  image_ = image;

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

  shader_.init(cfg, scene);

  int num_lights;
  if (shader_.isAo()) {
    num_lights = cfg.ao_samples;
  } else {
    CHECK_GT(lights_.size(), 0);
    num_lights = lights_.size();
  }

  pcontext_.resize(ndomains, cfg.bounces, cfg.nthreads, cfg.pixel_samples,
                   num_lights, image_);

  tile_list_.init(cfg.image_w, cfg.image_h, cfg.pixel_samples, nranks, rank,
                  cfg.maximum_num_screen_space_samples_per_rank);
  CHECK(!tile_list_.empty());

  tcontexts_.resize(cfg.nthreads);
  for (auto &tc : tcontexts_) {
    tc.resize(ndomains, cfg.pixel_samples, tile_list_.getLargestBlockingTile(),
              image_, cfg.bounces);
  }
}

template <typename ShaderT>
void Tracer<ShaderT>::genSingleEyes(int image_w, float orgx, float orgy,
                                    float orgz, spray::Tile tile,
                                    RayBuf<Ray> *ray_buf) {
  Ray *rays = ray_buf->rays;
#pragma omp for collapse(2) schedule(static, 1)
  for (int y = tile.y; y < tile.y + tile.h; ++y) {
    for (int x = tile.x; x < tile.x + tile.w; ++x) {
      int y0 = y - tile.y;
      int bufid_offset = y0 * tile.w;
      int pixid_offset = y * image_w;
      int x0 = x - tile.x;
      int bufid = bufid_offset + x0;
      int pixid = pixid_offset + x;
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

      ray->samid = bufid;

      ray->w[0] = 1.f;
      ray->w[1] = 1.f;
      ray->w[2] = 1.f;

      ray->depth = 0;
      ray->history[0] = SPRAY_FLOAT_INF;
      ray->committed = 0;
    }
  }
}

template <typename ShaderT>
void Tracer<ShaderT>::genMultiEyes(int image_w, float orgx, float orgy,
                                   float orgz, spray::Tile tile,
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
#ifdef SPRAY_GLOG_CHECK
        CHECK_LT(bufid, ray_buf->num);
#endif
        Ray *ray = &rays[bufid];

        ray->org[0] = orgx;
        ray->org[1] = orgy;
        ray->org[2] = orgz;

        ray->pixid = pixid;

        RandomSampler sampler;
        RandomSampler_init(sampler, bufid);

        float fx = (float)(x) + RandomSampler_get1D(sampler);
        float fy = (float)(y) + RandomSampler_get1D(sampler);

        camera_->generateRay(fx, fy, ray->dir);

        ray->samid = bufid;

        ray->w[0] = 1.f;
        ray->w[1] = 1.f;
        ray->w[2] = 1.f;

        ray->depth = 0;

        ray->history[0] = SPRAY_FLOAT_INF;
        ray->committed = 0;
      }
    }
  }
}

template <typename ShaderT>
void Tracer<ShaderT>::isectDomsRads(RayBuf<Ray> buf, TContextType *tc) {
  tc->resetRstats();
#pragma omp for schedule(static, 1)
  for (std::size_t i = 0; i < buf.num; ++i) {
    Ray *ray = &buf.rays[i];
    tc->enqRad(scene_, ray);
  }
}

template <typename ShaderT>
void Tracer<ShaderT>::trace() {
#pragma omp parallel
  {
    while (!tile_list_.empty()) {
#pragma omp barrier
#pragma omp single
      {
        blocking_tile_ = tile_list_.front();
        tile_list_.pop();
      }

      TContextType *tcontext = &tcontexts_[omp_get_thread_num()];
      tcontext->resetMems();

#pragma omp single
      {
        pcontext_.reset();

        shared_eyes_.num = (std::size_t)(blocking_tile_.w * blocking_tile_.h) *
                           num_pixel_samples_;
#ifdef SPRAY_GLOG_CHECK
        CHECK(shared_eyes_.num);
#endif
        if (shared_eyes_.num) {
          shared_eyes_.rays =
              tcontext->template allocMemIn<Ray>(shared_eyes_.num);
        }
      }

      if (shared_eyes_.num) {
        glm::vec3 cam_pos = camera_->getPosition();
        if (num_pixel_samples_ > 1) {
          genMultiEyes(image_w_, cam_pos[0], cam_pos[1], cam_pos[2],
                       blocking_tile_, &shared_eyes_);

        } else {
          genSingleEyes(image_w_, cam_pos[0], cam_pos[1], cam_pos[2],
                        blocking_tile_, &shared_eyes_);
        }

#pragma omp barrier
        isectDomsRads(shared_eyes_, tcontext);
#pragma omp barrier
      }

      pcontext_.isectPrims<SceneType, ShaderT>(scene_, shader_, tcontext);
    }
  }
  tile_list_.reset();
}

template <typename ShaderT>
void Tracer<ShaderT>::traceInOmp() {
  while (!tile_list_.empty()) {
#pragma omp barrier
#pragma omp single
    {
      blocking_tile_ = tile_list_.front();
      tile_list_.pop();
    }

    TContextType *tcontext = &tcontexts_[omp_get_thread_num()];
    tcontext->resetMems();

#pragma omp single
    {
      pcontext_.reset();

      shared_eyes_.num = (std::size_t)(blocking_tile_.w * blocking_tile_.h) *
                         num_pixel_samples_;
#ifdef SPRAY_GLOG_CHECK
      CHECK(shared_eyes_.num);
#endif
      if (shared_eyes_.num) {
        shared_eyes_.rays =
            tcontext->template allocMemIn<Ray>(shared_eyes_.num);
      }
    }

    if (shared_eyes_.num) {
      glm::vec3 cam_pos = camera_->getPosition();
      if (num_pixel_samples_ > 1) {
        genMultiEyes(image_w_, cam_pos[0], cam_pos[1], cam_pos[2],
                     blocking_tile_, &shared_eyes_);

      } else {
        genSingleEyes(image_w_, cam_pos[0], cam_pos[1], cam_pos[2],
                      blocking_tile_, &shared_eyes_);
      }

#pragma omp barrier
      isectDomsRads(shared_eyes_, tcontext);
#pragma omp barrier
    }

    pcontext_.isectPrims<SceneType, ShaderT>(scene_, shader_, tcontext);
  }
#pragma omp single
  tile_list_.reset();
}

}  // namespace ooc
}  // namespace spray

