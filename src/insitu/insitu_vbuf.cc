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

#include "glog/logging.h"

#include "insitu/insitu_ray.h"
#include "insitu/insitu_vbuf.h"

namespace spray {
namespace insitu {

void VBuf::resize(const spray::Tile& tile, int num_pixel_samples,
                  int total_num_light_samples) {
  tile_ = tile;
  num_pixel_samples_ = num_pixel_samples;
  total_num_light_samples_ = total_num_light_samples;

  // t buffer
  std::size_t size = tbufSize(tile, num_pixel_samples);
  CHECK_GT(size, 0);
  if (size > tbuf_capacity_) {
    FreeAligned(tbuf_0_);
    FreeAligned(tbuf_1_);

    tbuf_0_ = AllocAligned<float>(size);
    tbuf_1_ = AllocAligned<float>(size);
    CHECK_NOTNULL(tbuf_0_);
    CHECK_NOTNULL(tbuf_1_);

    tbuf_capacity_ = size;
    tbuf_size_ = size;
  } else {
    tbuf_size_ = size;
  }

  // o buffer
  size = obufSize(size, total_num_light_samples);

  if (size > obuf_capacity_) {
    FreeAligned(obuf_);

    obuf_ = AllocAligned<uint32_t>(size);
    CHECK_NOTNULL(obuf_);

    obuf_capacity_ = size;
    obuf_size_ = size;
  } else {
    obuf_size_ = size;
  }

  tbuf_in_ = tbuf_0_;
  tbuf_out_ = tbuf_1_;

  reset();
}

std::size_t VBuf::obufIndex(int samid, int light, std::size_t* bit) const {
  //
  std::size_t id = (std::size_t)samid * total_num_light_samples_ + light;

  std::size_t obuf_index = id >> 5;  // id / 32
  *bit = id - (obuf_index << 5);     // id % 32 (bit index)
#ifdef SPRAY_GLOG_CHECK
  CHECK_LT(obuf_index, obuf_size_);
  CHECK_LT(*bit, 32);
  CHECK_GE(*bit, 0);
#endif
  return obuf_index;
}

void VBuf::setOBuf(int samid, int light) {
  std::size_t bit;
  std::size_t i = obufIndex(samid, light, &bit);
  uint32_t o = obuf_[i];
  obuf_[i] = o | (1 << bit);
}

bool VBuf::occluded(int samid, int light) const {
  std::size_t bit;
  std::size_t i = obufIndex(samid, light, &bit);
  uint32_t o = obuf_[i];
  return ((o >> bit) & 1);
}

void VBuf::colorTBuf(spray::HdrImage* image) {
  float color[3];
  color[0] = color[1] = color[2] = 0.5f;
  for (int y = tile_.y; y < tile_.y + tile_.h; ++y) {
    for (int x = tile_.x; x < tile_.x + tile_.w; ++x) {
      int y0 = y - tile_.y;
      int bufid_offset = y0 * tile_.w;
      int pixid_offset = y * image->w;
      int x0 = x - tile_.x;
      int bufid = bufid_offset + x0;
      int pixid = pixid_offset + x;
      if (!std::isinf(tbuf_out_[bufid])) {
        image->add(pixid, color);
      }
    }
  }
}

}  // namespace insitu
}  // namespace spray
