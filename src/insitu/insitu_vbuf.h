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

#include <cstdint>
#include <vector>

#include "glog/logging.h"
#include "pbrt/memory.h"

#include "display/image.h"
#include "insitu/insitu_ray.h"
#include "render/spray.h"
#include "render/tile.h"
#include "utils/profiler_util.h"

namespace spray {
namespace insitu {

struct Ray;

class VBuf {
 public:
  VBuf()
      : tbuf_0_(nullptr),
        tbuf_1_(nullptr),
        tbuf_capacity_(0),
        tbuf_size_(0),
        obuf_(nullptr),
        obuf_capacity_(0),
        obuf_size_(0),
        num_pixel_samples_(0),
        total_num_light_samples_(0) {}

  ~VBuf() {
    FreeAligned(tbuf_0_);
    FreeAligned(tbuf_1_);
    FreeAligned(obuf_);
  }

  void resize(const spray::Tile& tile, int num_pixel_samples,
              int total_num_light_samples);

  std::size_t tbufSize(const Tile& tile, int num_pixel_samples) {
    return (std::size_t)tile.w * tile.h * num_pixel_samples;
  }

  std::size_t obufSize(std::size_t tbuf_size,
                       int total_num_light_samples) const {
    return (tbuf_size * total_num_light_samples + 31) / 32;
  }

  std::size_t getTbufSize() const { return tbuf_size_; }
  std::size_t getObufSize() const { return obuf_size_; }

  void resetTbufIn() {
    for (std::size_t i = 0; i < tbuf_size_; ++i) tbuf_in_[i] = SPRAY_FLOAT_INF;
  }

  void resetTbufOut() {
    for (std::size_t i = 0; i < tbuf_size_; ++i) tbuf_out_[i] = SPRAY_FLOAT_INF;
  }

  void resetObuf() {
    for (std::size_t i = 0; i < obuf_size_; ++i) obuf_[i] = (uint32_t)0;
  }

  void swapTbufs() { std::swap(tbuf_in_, tbuf_out_); }

  void colorTbuf(spray::HdrImage* image);

  bool correct(int samid, float t) const {
#ifdef SPRAY_GLOG_CHECK
    CHECK_LT(samid, tbuf_size_);
    CHECK_LE(tbuf_in_[samid], t);
#endif
    return (t == tbuf_in_[samid]);
  }

  bool tbufOutMiss(int samid) const {
#ifdef SPRAY_GLOG_CHECK
    CHECK_LT(samid, tbuf_size_);
#endif
    return std::isinf(tbuf_out_[samid]);
  }

  void setObuf(int samid, int light);

  bool occluded(int samid, int light) const;

  void compositeTbuf() {
#ifdef SPRAY_TIMING
    spray::tStartMPI(spray::TIMER_SYNC_VBUF);
#endif
    MPI_Allreduce(MPI_IN_PLACE, tbuf_out_, tbuf_size_, MPI_FLOAT, MPI_MIN,
                  MPI_COMM_WORLD);
#ifdef SPRAY_TIMING
    spray::tStop(spray::TIMER_SYNC_VBUF);
#endif
  }

  void compositeObuf() {
#ifdef SPRAY_TIMING
    spray::tStartMPI(spray::TIMER_SYNC_VBUF);
#endif
    MPI_Allreduce(MPI_IN_PLACE, obuf_, obuf_size_, MPI_UINT32_T, MPI_BOR,
                  MPI_COMM_WORLD);
#ifdef SPRAY_TIMING
    spray::tStop(spray::TIMER_SYNC_VBUF);
#endif
  }

  float getTbufOut(int samid) const { return tbuf_out_[samid]; }
  void setTbufOut(int samid, float t) { tbuf_out_[samid] = t; }

  uint32_t getObuf(std::size_t index) const { return obuf_[index]; }
  void setObufByIndex(std::size_t index, uint32_t value) {
#ifdef SPRAY_GLOG_CHECK
    CHECK_LT(index, obuf_size_);
#endif
    obuf_[index] = value;
  }

  bool equalToTbufOut(int samid, float t) { return (tbuf_out_[samid] == t); }

  bool updateTbufOut(float t_new, Ray* ray) {
    auto& t_old = tbuf_out_[ray->samid];
    if (t_new < t_old) {
      t_old = t_new;
      ray->t = t_new;
      return true;
    }
    return false;
  }

 private:
  void reset() {
    for (std::size_t i = 0; i < tbuf_size_; ++i) tbuf_in_[i] = SPRAY_FLOAT_INF;
    for (std::size_t i = 0; i < tbuf_size_; ++i) tbuf_out_[i] = SPRAY_FLOAT_INF;
    for (std::size_t i = 0; i < obuf_size_; ++i) obuf_[i] = (uint32_t)0;
  }

  std::size_t obufIndex(int samid, int light, std::size_t* bit) const;

 private:
  float* tbuf_0_;
  float* tbuf_1_;

  float* tbuf_in_;   // used to commit/retire rays (previous ts)
  float* tbuf_out_;  // used to update t values (current ts)

  std::size_t tbuf_capacity_;
  std::size_t tbuf_size_;

  uint32_t* obuf_;
  std::size_t obuf_capacity_;
  std::size_t obuf_size_;

  spray::Tile tile_;
  int num_pixel_samples_;
  int total_num_light_samples_;
};

}  // namespace insitu
}  // namespace spray
