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
#include "insitu/insitu_comm.h"
#include "insitu/insitu_isector.h"
#include "insitu/insitu_ray.h"
#include "insitu/insitu_vbuf.h"
#include "insitu/insitu_work.h"
#include "insitu/insitu_work_stats.h"
#include "render/camera.h"
#include "render/data_partition.h"
#include "render/domain.h"
#include "render/light.h"
#include "render/qvector.h"
#include "render/reflection.h"
#include "render/spray.h"
#include "render/tile.h"
#include "utils/profiler_util.h"

namespace spray {
namespace insitu {

template <typename CacheT, typename ShaderT, typename SceneT>
class SingleThreadTracer {
 public:
  void trace();
  void traceInOmp() {
    std::cout << "[warning] tracing in omp parallel region unsupported\n";
  }
  int type() const { return TRACER_TYPE_SPRAY_INSITU_1_THREAD; }

 public:
  void init(const Config &cfg, const Camera &camera, SceneT *scene,
            HdrImage *image);

 private:
  ShaderT shader_;
  spray::TileList tile_list_;
  Comm comm_;
  VBuf vbuf_;

  SceneInfo sinfo_;
  spray::RTCRayIntersection rtc_isect_;
  RTCRay rtc_ray_;

  Tile blocking_tile_, stripe_;
  RayBuf<Ray> shared_eyes_;

 private:
  void sendRays();
  void send(bool shadow, int domain_id, int dest, std::queue<Ray *> *q);
  void procLocalQs();
  void procRecvQs();
  void procRads(int id, Ray *rays, int64_t count);
  void procShads(int id, Ray *rays, int64_t count);

  void procRad(int id, Ray *ray);
  void procShad(int id, Ray *ray);

  void filterRq2(int id);
  void filterSq2(int id);
  void procFrq2();
  void procFsq2();
  void procCachedRq();
  void procRetireQ();
  void retireBackground();

  void populateRadWorkStats();
  void populateWorkStats();

#ifdef SPRAY_GLOG_CHECK
  void checkQs() {
    CHECK(rqs_.empty());
    CHECK(sqs_.empty());
    CHECK(frq2_.empty());
    CHECK(fsq2_.empty());
    CHECK(cached_rq_.empty());
    CHECK(recv_rq_.empty());
    CHECK(recv_sq_.empty());
    CHECK(retire_q_.empty());
    CHECK(bg_retire_q_.empty());
  }
#endif

 private:
  const spray::Camera *camera_;
  const spray::InsituPartition *partition_;
  std::vector<spray::Light *> lights_;
  SceneT *scene_;
  spray::HdrImage *image_;
  Isector<CacheT, SceneT> isector_;

  spray::QVector<Ray *> rqs_;
  spray::QVector<Ray *> sqs_;

  std::queue<msg_word_t *> recv_rq_;
  std::queue<msg_word_t *> recv_sq_;

  std::queue<Ray *> rq2_;
  std::queue<Ray *> sq2_;

  struct IsectInfo {
    int domain_id;
    spray::RTCRayIntersection *isect;
    Ray *ray;
  };

  struct OcclInfo {
    int domain_id;
    Ray *ray;
  };

  std::queue<IsectInfo> cached_rq_;

  std::queue<IsectInfo> frq2_;
  std::queue<OcclInfo> fsq2_;

  std::queue<Ray *> retire_q_;     ///< Retire queue for foreground colors.
  std::queue<Ray *> bg_retire_q_;  ///< Retire queue for background colors.

  WorkStats work_stats_;

  spray::MemoryArena *mem_in_;
  spray::MemoryArena *mem_out_;
  spray::MemoryArena mem_0_;
  spray::MemoryArena mem_1_;

 private:
  int ray_depth_;

 private:
  int rank_;
  int num_ranks_;
  int num_domains_;
  int num_pixel_samples_;
  double one_over_num_pixel_samples_;
  int num_bounces_;
  int num_threads_;
  int num_lights_;
  int image_w_;
  int image_h_;
  glm::vec3 bg_color_;
};

}  // namespace insitu
}  // namespace spray

#define SPRAY_INSITU_SINGLE_THRAED_TRACER_INL
#include "insitu/insitu_singlethread_tracer.inl"
#undef SPRAY_INSITU_SINGLE_THRAED_TRACER_INL

