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

#include "ooc/ooc_vbuf.h"

#include "glog/logging.h"

#include "ooc/ooc_ray.h"

namespace spray {
namespace ooc {

void VBuf::resize(const Tile& tile, int num_pixel_samples) {
  tile_ = tile;
  num_pixel_samples_ = num_pixel_samples;

  std::size_t size = tbufSize(tile, num_pixel_samples);
  CHECK_GT(size, 0);
  tbuf_.resize(size);

  for (auto& t : tbuf_) t = SPRAY_FLOAT_INF;
}

bool VBuf::correct(const Ray& ray) const {
  const auto ray_depth = ray.depth;
  auto offset = (std::size_t)ray.samid * SPRAY_SPECU_HISTORY_SIZE;
  bool correct = true;
  for (int i = 0; i < ray_depth; ++i) {
    if (ray.history[i] != tbuf_[offset + i]) {
      correct = false;
      break;
    }
  }
  return correct;
}

bool VBuf::update(float t, Ray* ray) {
  const auto ray_depth = ray->depth;
  int update_pos = ray_depth;

  const auto offset = (std::size_t)ray->samid * SPRAY_SPECU_HISTORY_SIZE;

  // scan past t values only (history check)
  // so no checking for camera rays
  for (unsigned i = 0; i < ray_depth; ++i) {
    auto local_t = ray->history[i];
    auto global_t = tbuf_[offset + i];

    if (local_t < global_t) {
      update_pos = i;
      break;
    } else if (local_t > global_t) {
      update_pos = -1;
      break;
    }
  }

  // update if needed
  bool updated;

  if (update_pos == -1) {
    updated = false;

  } else if (update_pos < ray_depth) {
    updated = true;

    ray->history[ray_depth] = t;

    // NOTE: ray_depth inclusive (<=, not <)
    for (int i = update_pos; i <= ray_depth; ++i) {
      tbuf_[offset + i] = ray->history[i];
    }

  } else {
    auto idx = offset + ray_depth;

    if (t < tbuf_[idx]) {
      updated = true;

      tbuf_[idx] = t;
      ray->history[ray_depth] = t;

    } else {
      updated = false;
    }
  }

  if (updated) {
    for (unsigned i = ray_depth + 1; i < SPRAY_SPECU_HISTORY_SIZE; ++i) {
      tbuf_[offset + i] = SPRAY_FLOAT_INF;
    }
  }

  return updated;
}

}  // namespace ooc
}  // namespace spray
