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
#include "pbrt/memory3.h"

#include "baseline/baseline_comm.h"
#include "baseline/baseline_commit_buffer.h"
#include "baseline/baseline_isector.h"
#include "baseline/baseline_qstats.h"
#include "baseline/baseline_ray.h"
#include "baseline/baseline_ray_stats.h"
#include "baseline/baseline_schedulers.h"
#include "baseline/baseline_tilers.h"
#include "display/image.h"
#include "materials/reflection.h"
#include "partition/arena_queue.h"
#include "partition/block_buffer.h"
#include "partition/data_partition.h"
#include "partition/domain.h"
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
namespace baseline {

template <typename CacheT, typename ScheduleT, typename ShaderT>
class InsituTracer {
 public:
  virtual ~InsituTracer() {
    for (std::size_t i = 0; i < domain_locks_.size(); ++i) {
      omp_destroy_lock(&domain_locks_[i]);
    }
  }

  void init(const Config &cfg, const Camera &camera, Scene<CacheT> *scene,
            HdrImage *image);

  void trace();
  void traceInOmp() {
    std::cout << "[warning] tracing in omp parallel region unsupported\n";
  }
  int type() const { return TRACER_TYPE_BASELINE_INSITU; }

  void terminate() { img_sched_.terminate(); }

 protected:
  void initCommon(const Config &cfg, const Camera &camera, Scene<CacheT> *scene,
                  HdrImage *image);

 protected:
  void genEyeRays(int ndomains, int nsamples, Tile tile, DRay *ray_buf);

  void isectEyeDomains(int num_domains, std::size_t num_rays, DRay *ray_buf,
                       ArenaQs<DRayQItem> *qs);

 protected:
  int assignedProcess(unsigned score_index) {
    int rank = score_index % mpi::size();
#ifdef SPRAY_GLOG_CHECK
    CHECK_LT(rank, mpi::size());
#endif
    return rank;
  }

 protected:
  void schedule(int ndomains, const ArenaQs<DRayQItem> &qs,
                const ArenaQs<DRayQItem> &sqs, QStats *stats);

  void aggregateStats(int ndomains, const QStats &src_stats,
                      QStats *dest_stats);

  void spliceQs(int tid, int id, int ndomains,
                const std::vector<RayCount> &sched,
                ArenaQs<DRayQItem> *input_qs, BlockBuffer *send_buf,
                BlockBuffer *recv_buf);

  void spliceQsMultiProcesses(int tid, int id, int ndomains,
                              const std::vector<RayCount> &sched,
                              ArenaQs<DRayQItem> *input_qs,
                              BlockBuffer *send_buf, BlockBuffer *recv_buf);

  void spliceRecvQs(int tid, int id, int ndomains, ArenaQs<DRayQItem> *qs,
                    BlockBuffer *recv_buf);

  void processRays(int tid, int id, int ndomains, int nbounces,
                   ArenaQs<DRayQItem> *qs, ArenaQs<DRayQItem> *sqs,
                   MemoryArena *arena,
                   DomainIntersector<CacheT> *domain_isector, RTCRay *rtc_ray,
                   RTCRayIntersection *isect, CommitBufferB *retire_buf,
                   DRayQ *temp_q);

  void processRay2(int id, int ndomains, DRay *data, RTCRay *rtc_ray,
                   RTCRayIntersection *isect, ArenaQs<DRayQItem> *qs,
                   ArenaQs<DRayQItem> *sqs, MemoryArena *arena,
                   DomainIntersector<CacheT> *domain_isector,
                   CommitBufferB *retire_buf, DRayQ *temp_q);

  void processShadow(int id, int ndomains, DRay *data, RTCRay *rtc_ray,
                     ArenaQs<DRayQItem> *sqs, MemoryArena *arena,
                     DomainIntersector<CacheT> *domain_isector,
                     CommitBufferB *retire_buf);

  void resetSentQs(int ndomains, const std::vector<RayCount> &sched,
                   ArenaQs<DRayQItem> *qs, ArenaQs<DRayQItem> *sqs);

 protected:
  BlockBuffer send_rbuf_;
  BlockBuffer send_sbuf_;

  BlockBuffer recv_rbuf_;
  BlockBuffer recv_sbuf_;

  InclusiveScan<int> scan_;

  Comm<DRayQItem, DRay, DRayOutgoingCopier, DRayIncomingCopier> comm_;

 protected:
  QStats qstats_;

  MemoryArena eyeray_arena_;

  // scheduler
  ImgSchedSingle img_sched_;
  ScheduleT ray_sched_;

  ShaderT shader_;

  // pointers
  std::vector<Light *> lights_;  // copied lights
  const Camera *camera_;
  Scene<CacheT> *scene_;
  HdrImage *image_;

  // parameters
  int pixel_samples_;
  int bounces_;

  int ndomains_;
  int rank_;
  int nranks_;

  std::vector<omp_lock_t> domain_locks_;

#ifdef SPRAY_GLOG_CHECK
  std::vector<std::size_t> qs_counts_;
  std::vector<std::size_t> sqs_counts_;
#endif
#ifdef SPRAY_PROFILE_COUNTERS
 private:
  void profileRaysSpawned(const ArenaQs<DRayQItem> &qs);
#endif
};

}  // namespace baseline
}  // namespace spray

#define SPRAY_BASELINE_INSITU_TRACER_INL
#include "baseline/baseline_insitu_tracer.inl"
#undef SPRAY_BASELINE_INSITU_TRACER_INL

