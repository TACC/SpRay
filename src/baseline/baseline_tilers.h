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

#include <iostream>
#include <queue>
#include <vector>

#include "glog/logging.h"

#include "partition/tile.h"
#include "render/config.h"

namespace spray {
namespace baseline {

//! Uses entire image plane as a single, large tile
class ImgSchedNop {
 public:
  void terminate() {}
  void initialize(unsigned image_w, unsigned image_h) {
    tile_.x = 0;
    tile_.y = 0;
    tile_.w = image_w;
    tile_.h = image_h;
  }
  Tile schedule() { return tile_; }
  Tile getLargestTile(unsigned image_w, unsigned image_h) { return tile_; }
  Tile getLargestTile() { return tile_; }

 private:
  Tile tile_;
};

//! Assigns horizontal stripes to processes, i.e. one stripe per each process.
class ImgSchedSingle {
 public:
  void terminate() {}
  void initialize(unsigned image_w, unsigned image_h) {
    done_ = false;
    unsigned h = image_h / mpi::size();

    tile_.x = 0;
    tile_.w = image_w;
    tile_.h = 0;

    for (unsigned rank = 0; rank < mpi::size(); ++rank) {
      if (rank == mpi::rank()) {
        tile_.y = rank * h;

        if (rank == mpi::size() - 1) {
          tile_.h = image_h - tile_.y;
        } else {
          tile_.h = h;
        }
#ifdef SPRAY_GLOG_CHECK
        CHECK(image_h >= tile_.y);
#endif
      }
    }
#ifdef SPRAY_GLOG_CHECK
    CHECK(tile_.getArea() > 0);
#endif
  }
  bool isDone() const { return done_; }
  void reset() { done_ = false; }

  Tile schedule() {
    done_ = true;
    return tile_;
  }

  Tile getLargestTile() { return tile_; }

 protected:
  Tile tile_;
  bool done_;
};

}  // namespace baseline
}  // namespace spray
