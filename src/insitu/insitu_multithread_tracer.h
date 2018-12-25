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

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <numeric>
#include <queue>
#include <vector>
#include <string.h>

#include "embree/random_sampler.h"
#include "glog/logging.h"

#include "display/image.h"
#include "insitu/insitu_comm.h"
#include "insitu/insitu_isector.h"
#include "insitu/insitu_ray.h"
#include "insitu/insitu_tcontext.h"
#include "insitu/insitu_tiler.h"
#include "insitu/insitu_vbuf.h"
#include "insitu/insitu_work.h"
#include "insitu/insitu_work_stats.h"
#include "materials/reflection.h"
#include "partition/data_partition.h"
#include "partition/domain.h"
#include "partition/qvector.h"
#include "partition/tile.h"
#include "render/config.h"
#include "render/spray.h"
#include "scene/camera.h"
#include "scene/light.h"
#include "scene/scene.h"
#include "utils/comm.h"
#include "utils/profiler_util.h"
#include "utils/scan.h"

namespace spray {
namespace insitu {

template <typename CacheT, typename ShaderT>
class MultiThreadTracer {
 public:
  void trace();
  void traceInOmp();
  int type() const { return TRACER_TYPE_SPRAY_INSITU_N_THREADS; }

 public:
  void init(const Config &cfg, const Camera &camera, Scene<CacheT> *scene,
            HdrImage *image);

 private:
  typedef TContext<CacheT, ShaderT> TContextType;

  std::vector<TContextType> tcontexts_;

  ShaderT shader_;
  Tiler tiler_;
  Comm comm_;
  VBuf vbuf_;

  SceneInfo sinfo_;
  spray::RTCRayIntersection rtc_isect_;
  RTCRay rtc_ray_;

 private:
  void genSingleEyes(int image_w, float orgx, float orgy, float orgz,
                     int base_tile_y, Tile tile, RayBuf<Ray> *ray_buf);
  void genMultiEyes(int image_w, float orgx, float orgy, float orgz,
                    int base_tile_y, Tile tile, RayBuf<Ray> *ray_buf);

  void sendRays(int tid, TContextType *tcontext);
  void send(bool shadow, int tid, int domain_id, int dest, std::size_t num_rays,
            TContextType *tcontext);
  void procLocalQs(int tid, int ray_depth, TContextType *tcontext);
  void procRecvQs(int ray_depth, TContextType *tcontext);
  void procRecvRads(int ray_depth, int id, Ray *rays, int64_t count,
                    TContextType *tcontext);
  void procRecvShads(int id, Ray *rays, int64_t count, TContextType *tcontext);

  void procCachedRq(int ray_depth, TContextType *tcontext);

  void populateRadWorkStats(TContextType *tcontext);
  void populateWorkStats(TContextType *tcontext);

 private:
  const spray::Camera *camera_;
  const spray::InsituPartition *partition_;
  std::vector<spray::Light *> lights_;  // copied lights
  Scene<CacheT> *scene_;
  spray::HdrImage *image_;

  std::queue<msg_word_t *> recv_rq_;
  std::queue<msg_word_t *> recv_sq_;

  msg_word_t *recv_message_;

  WorkStats work_stats_;  // number of blocks to process

  spray::ThreadStatus thread_status_;
  spray::InclusiveScan<std::size_t> scan_;
  WorkSendMsg<Ray, MsgHeader> *send_work_;

 private:
  Tile mytile_;
  Tile image_tile_;

  RayBuf<Ray> shared_eyes_;
  int done_;

 private:
  int rank_;
  int num_ranks_;
  int num_domains_;
  int num_pixel_samples_;
  int num_bounces_;
  int num_threads_;
  int num_lights_;
  int image_w_;
  int image_h_;
};

}  // namespace insitu
}  // namespace spray

#define SPRAY_INSITU_MULTI_THREAD_TRACER_INL
#include "insitu/insitu_multithread_tracer.inl"
#undef SPRAY_INSITU_MULTI_THREAD_TRACER_INL

