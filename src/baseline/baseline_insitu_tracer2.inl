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

#if !defined(SPRAY_BASELINE_INSITU_TRACER2_INL)
#error An implementation of InsituTracer2
#endif

namespace spray {
namespace baseline {

template <typename CacheT, typename ShaderT, typename SceneT>
void InsituTracer2<CacheT, ShaderT, SceneT>::init(const Config &cfg,
                                                  const Camera &camera,
                                                  SceneT *scene,
                                                  HdrImage *image) {
  camera_ = &camera;
  scene_ = scene;
  partition_ = &(scene->getInsituPartition());

  lights_ = scene->getLights();  // copy lights
  image_ = image;

  num_ranks_ = mpi::worldSize();
  rank_ = mpi::worldRank();
  num_domains_ = scene->getNumDomains();
  num_pixel_samples_ = cfg.pixel_samples;
  one_over_num_pixel_samples_ = 1.0 / (double)num_pixel_samples_;
  num_bounces_ = cfg.bounces;
  num_threads_ = cfg.nthreads;
  image_w_ = cfg.image_w;
  image_h_ = cfg.image_h;
  bg_color_ = cfg.bg_color;

  CHECK_GT(rank_, -1);
  CHECK_GT(num_ranks_, 0);
  CHECK_GT(num_domains_, 0);
  CHECK_GT(num_pixel_samples_, 0);
  CHECK_GT(num_bounces_, 0);
  CHECK_GT(image_w_, 0);
  CHECK_GT(image_h_, 0);

  // shader_.init(cfg, *scene);

  // int total_num_light_samples;
  // if (shader_.isAo()) {
  //   num_lights_ = cfg.ao_samples;
  //   total_num_light_samples = num_lights_;
  // } else {
  //   CHECK_GT(lights_.size(), 0);
  //   num_lights_ = lights_.size();
  //   std::size_t num_lights = lights_.size();

  //   total_num_light_samples = 0;
  //   for (std::size_t i = 0; i < lights_.size(); ++i) {
  //     if (lights_[i]->isAreaLight()) {
  //       total_num_light_samples += cfg.ao_samples;
  //     } else {
  //       ++total_num_light_samples;
  //     }
  //   }
  // }

  tile_list_.init(cfg.image_w, cfg.image_h, cfg.pixel_samples, num_ranks_,
                  rank_, cfg.maximum_num_screen_space_samples_per_rank);

  CHECK(!tile_list_.empty());

  rqs_.resize(num_domains_);
  sqs_.resize(num_domains_);
  // work_stats_.resize(nranks, cfg.nthreads, num_domains_);

  mem_in_ = &mem_0_;
  mem_out_ = &mem_1_;

  // vbuf_.resize(tile_list_.getLargestBlockingTile(), cfg.pixel_samples,
  //              total_num_light_samples);
}

template <typename CacheT, typename ShaderT, typename SceneT>
void InsituTracer2<CacheT, ShaderT, SceneT>::trace() {
  while (!tile_list_.empty()) {
    tile_list_.front(&blocking_tile_, &stripe_);
    tile_list_.pop();
  }
}

}  // namespace insitu
}  // namespace baseline
