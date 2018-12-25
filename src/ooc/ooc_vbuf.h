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

#include "display/image.h"
#include "ooc/ooc_ray.h"
#include "ooc/ooc_tiler.h"
#include "render/spray.h"

namespace spray {
namespace ooc {

struct Ray;

class VBuf {
 public:
  VBuf() : num_pixel_samples_(0) {}

 private:
  std::vector<float> tbuf_;
  Tile tile_;
  int num_pixel_samples_;

 public:
  void resize(const Tile& tile, int num_pixel_samples);

  std::size_t tbufSize(const Tile& tile, int num_pixel_samples) {
    return (std::size_t)tile.w * tile.h * num_pixel_samples *
           SPRAY_HISTORY_SIZE;
  }

 public:
  void reset() {
    for (auto& t : tbuf_) t = SPRAY_FLOAT_INF;
  }
  bool correct(const Ray& ray) const;
  bool update(float t, Ray* ray);
};

}  // namespace ooc
}  // namespace spray
