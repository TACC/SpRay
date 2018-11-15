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
#include <string.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <numeric>
#include <queue>
#include <vector>

#include "embree/random_sampler.h"
#include "glog/logging.h"

#include "display/image.h"
#include "materials/reflection.h"
#include "ooc/ooc_pcontext.h"
#include "ooc/ooc_ray.h"
#include "ooc/ooc_tcontext.h"
#include "ooc/ooc_tiler.h"
#include "partition/domain.h"
#include "renderers/spray.h"
#include "scene/camera.h"
#include "scene/light.h"
#include "scene/scene.h"
#include "utils/profiler_util.h"

namespace spray {
namespace ooc {

struct RayBuf {
  RayBuf() { reset(); }
  void reset() {
    num = 0;
    rays = nullptr;
  }
  std::size_t num;
  Ray *rays;
};

template <typename CacheT, typename ShaderT>
class Tracer {
 public:
  void trace();
  void traceInOmpParallel();

 public:
  void init(const Config &cfg, const Camera &camera, Scene<CacheT> *scene,
            HdrImage *image);

 private:
  ShaderT shader_;
  PContext pcontext_;
  std::vector<TContext> tcontexts_;

 private:
  void genSingleEyes(int image_w, float orgx, float orgy, float orgz, Tile tile,
                     RayBuf *ray_buf);
  void genMultiEyes(int image_w, float orgx, float orgy, float orgz, Tile tile,
                    RayBuf *ray_buf);

  void isectDomsRads(RayBuf buf, TContext *tc);
  void isectPrimsRads(TContext *tc);

 private:
  const spray::Camera *camera_;
  std::vector<spray::Light *> lights_;  // copied lights
  Scene<CacheT> *scene_;
  spray::HdrImage *image_;

  Tile image_tile_;
  Tile mytile_;
  RayBuf shared_eyes_;

 private:
  int rank_;
  int num_ranks_;
  int num_domains_;
  int num_pixel_samples_;
  int num_bounces_;
  int num_threads_;
  int image_w_;
  int image_h_;
};

}  // namespace ooc
}  // namespace spray

#define SPRAY_OOC_TRACER_INL
#include "ooc/ooc_tracer.inl"
#undef SPRAY_OOC_TRACER_INL

