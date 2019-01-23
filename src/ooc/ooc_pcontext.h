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

#include <cstring>
#include <queue>
#include <utility>

// clang-format off
#include <embree2/rtcore_ray.h>
#include <embree2/rtcore_geometry.h>
// clang-format on
#include <omp.h>

#include "glm/glm.hpp"
#include "glog/logging.h"
#include "pbrt/memory.h"

#include "display/image.h"
#include "ooc/ooc_domain_stats.h"
#include "ooc/ooc_isector.h"
#include "ooc/ooc_ray.h"
#include "ooc/ooc_tcontext.h"
#include "render/qvector.h"
#include "render/tile.h"
#include "render/rays.h"
#include "utils/profiler_util.h"
#include "utils/scan.h"

namespace spray {
namespace ooc {

class PContext {
 public:
  void resize(int ndomains, int max_num_bounces, int num_threads,
              int num_pixel_samples, int num_lights, spray::HdrImage* image);

 private:
  int num_domains_;
  int max_num_bounces_;
  int num_pixel_samples_;
  spray::HdrImage* image_;

 public:
  void reset() { bounce_num_ = 0; }

  void incrementBounceNum() { ++bounce_num_; }
  int getBounceNum() { return bounce_num_; }

 private:
  int bounce_num_;
  std::vector<omp_lock_t> domain_locks_;

 public:
  void resetAndSwap() { ++bounce_num_; }

 private:
  DomainStats rstats_;

 private:
  void mergeStats(const DomainStats& src_stats, DomainStats* dest_stats);

 private:
  std::vector<spray::InclusiveScan<int>> scans_;
  spray::SceneInfo sinfo_;

 public:
  template <typename SceneT, typename ShaderT>
  void isectPrims(SceneT* scene, ShaderT& shader,
                  TContext<SceneT, ShaderT>* tcontext);
};

template <typename SceneT, typename ShaderT>
void PContext::isectPrims(SceneT* scene, ShaderT& shader,
                          TContext<SceneT, ShaderT>* tcontext) {
  int tid = omp_get_thread_num();
  int ray_depth = 0;

  do {
    tcontext->procPendingQ(scene);
    tcontext->resetVBuf();

    while (1) {
      bool sqs_in_empty = tcontext->sqsInEmpty();
      bool rqs_empty = tcontext->rqsEmpty();
      bool retire_q_empty = tcontext->retireQEmpty();

      bool all_empty = sqs_in_empty && rqs_empty && retire_q_empty;
      scans_[0].set(tid, all_empty);
#pragma omp barrier
#pragma omp single
      { scans_[0].scan(); }

      if (scans_[0].sum() == scans_[0].size()) {  // all qs empty
#ifdef SPRAY_GLOG_CHECK
        CHECK(tcontext->commitQEmpty());
#endif
        break;
      }

      bool input_empty = sqs_in_empty && rqs_empty;
      scans_[1].set(tid, input_empty);
#pragma omp barrier
#pragma omp single
      { scans_[1].scan(); }

      if (scans_[1].sum() < scans_[1].size()) {  // input q not empty
#pragma omp single
        { rstats_.reset(); }
        mergeStats(tcontext->getRstats(), &rstats_);
#pragma omp barrier

#pragma omp single
        { rstats_.schedule(); }

        for (int i = 0; i < num_domains_; ++i) {
          //
          int id = rstats_.getDomainId(i);

          tcontext->filterQs(id);
          scans_[2].set(tid, tcontext->allFilterQsEmpty());
#pragma omp barrier
#pragma omp single
          { scans_[2].scan(); }

          if (scans_[2].sum() < scans_[2].size()) {
#pragma omp single
            {
#ifdef SPRAY_TIMING
              spray::tStart(spray::TIMER_LOAD);
#endif
              scene->load(id, &sinfo_);
#ifdef SPRAY_TIMING
              spray::tStop(spray::TIMER_LOAD);
#endif
            }

            tcontext->procFilterQs(id, scene, sinfo_, shader, ray_depth);
          }
#pragma omp barrier
        }
      }

#ifdef SPRAY_GLOG_CHECK
      CHECK(tcontext->sqsInEmpty());
#endif
#pragma omp critical(cs_pcontext_tcontext_retire)
      { tcontext->retire(); }
      tcontext->swapQs();
    }  // while (1)

#pragma omp critical(cs_pcontext_tcontext_retire)
    { tcontext->retire(); }

#ifdef SPRAY_GLOG_CHECK
    tcontext->checkQs();
#endif
    tcontext->resetMemIn();
    tcontext->swapMems();

    ray_depth += SPRAY_HISTORY_SIZE;

    scans_[3].set(tid, tcontext->pendingQEmpty());
#pragma omp barrier
#pragma omp single
    { scans_[3].scan(); }

  } while (scans_[3].sum() <
           scans_[3].size());  // while all pending Qs not empty
}

}  // namespace ooc
}  // namespace spray
