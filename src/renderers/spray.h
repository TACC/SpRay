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
#include <limits>

#include "glog/logging.h"

#include "cmake_config.h"  // auto generated by cmake
#include "utils/profiler.h"

#define SPRAY_ALIGN(...) __attribute__((aligned(__VA_ARGS__)))

#define RANK_THREAD \
  "[rank " << mpi::rank() << "] [thread " << omp_get_thread_num() << "] "

#define SPRAY_INT_MAX std::numeric_limits<int>::max()
#define SPRAY_REAL_MIN std::numeric_limits<float>::min()
#define SPRAY_REAL_MAX std::numeric_limits<float>::max()
#define SPRAY_SIZE_MIN std::numeric_limits<std::size_t>::MIN()
#define SPRAY_SIZE_MAX std::numeric_limits<std::size_t>::max()
#define SPRAY_REAL_INF std::numeric_limits<float>::infinity()
#define SPRAY_FLOAT_INF std::numeric_limits<float>::infinity()
#define SPRAY_FLOAT_MAX std::numeric_limits<float>::max()
#define SPRAY_DOUBLE_INF std::numeric_limits<double>::infinity()

#define SPRAY_RAY_EPSILON 0.001f

#define SPRAY_PI M_PI
#define SPRAY_TWO_PI 6.283185307179586f
#define SPRAY_ONE_OVER_PI 0.3183098861837907f
#define SPRAY_ONE_OVER_TWOPI 0.15915494309189535f
#define SPRAY_INV_PI SPRAY_ONE_OVER_PI
#define SPRAY_INV_TWOPI SPRAY_ONE_OVER_TWOPI
#define SPRAY_1_OVER_255 0.00392156862745098

#ifdef SPRAY_AVX2  // 512
#define SPRAY_RAY_PACKET_SIZE 16
#define SPRAY_RAY_PACKET_ALIGNMENT 64
#define SPRAY_RTC_RAYS RTCRay16
#define SPRAY_RTC_INTERSECT rtcIntersect16
#define SPRAY_RTC_OCCLUDED rtcOccluded16
#elif SPRAY_AVX  // 256
#define SPRAY_RAY_PACKET_SIZE 8
#define SPRAY_RAY_PACKET_ALIGNMENT 32
#define SPRAY_RTC_RAYS RTCRay8
#define SPRAY_RTC_INTERSECT rtcIntersect8
#define SPRAY_RTC_OCCLUDED rtcOccluded8
#elif SPRAY_SSE  // 128
#define SPRAY_RAY_PACKET_SIZE 4
#define SPRAY_RAY_PACKET_ALIGNMENT 16
#define SPRAY_RTC_RAYS RTCRay4
#define SPRAY_RTC_INTERSECT rtcIntersect4
#define SPRAY_RTC_OCCLUDED rtcOccluded4
#else
#define SPRAY_RAY_PACKET_SIZE 8
#define SPRAY_RAY_PACKET_ALIGNMENT 32
#define SPRAY_RTC_RAYS RTCRay8
#define SPRAY_RTC_INTERSECT rtcIntersect8
#define SPRAY_RTC_OCCLUDED rtcOccluded8
#endif

#define SPRAY_ROOT_PROCESS 0
#define SPRAY_HISTORY_SIZE SPRAY_SPECU_HISTORY_SIZE

namespace spray {

enum ViewMode {
  VIEW_MODE_DOMAIN,
  VIEW_MODE_PARTITION,
  VIEW_MODE_FILM,
  VIEW_MODE_GLFW,
  VIEW_MODE_WBVH,
  VIEW_MODE_WBVH_OVERLAY,
  VIEW_MODE_BVH,
  VIEW_MODE_TERMINATE
};

enum TracerType {
  TRACER_TYPE_SPRAY_OOC,
  TRACER_TYPE_SPRAY_INSITU_1_THREAD,
  TRACER_TYPE_SPRAY_INSITU_N_THREADS,
  TRACER_TYPE_BASELINE_IMAGE,
  TRACER_TYPE_BASELINE_INSITU
};

enum CameraCommand {
  CAM_NOP = 0,
  CAM_ZOOM,
  CAM_ROTATE,
  CAM_PAN,
  CAM_RESET,
  CAM_RESET_TO_CFG
};

struct MessageCommand {
  int done;
  int image_w;
  int image_h;
  int view_mode;
  int camera_cmd;
  float rotate_pan_dx;
  float rotate_pan_dy;
  float zoom_offset;
};

struct MpiComm {
  MpiComm() { size = 0; }
  int size;
  int rank;
};

class ThreadStatus {
 public:
  void resize(int num_threads) { status_.resize(num_threads); }
  void set(int tid) {
    CHECK_LT(tid, status_.size());
    status_[tid] = 1;
  }
  void clear(int tid) {
    CHECK_LT(tid, status_.size());
    status_[tid] = 0;
  }

  bool isAnySet() {
    for (std::size_t i = 0; i < status_.size(); ++i) {
      if (status_[i] == 1) {
        return true;
      }
    }
    return false;
  }

 private:
  std::vector<int> status_;
};

extern MpiComm global_mpi_comm;
extern Profiler global_profiler;

extern int global_num_active_frames;
extern int global_num_warmup_frames;
extern unsigned global_num_frames;

}  // namespace spray

