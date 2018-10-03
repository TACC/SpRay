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

#include "partition/tile.h"

#include <algorithm>

#include "utils/util.h"

namespace spray {

std::ostream& operator<<(std::ostream& os, const Tile& t) {
  os << "tile(" << t.x << "," << t.y << "," << t.w << "," << t.h << ")";
  return os;
}

TileCount::TileCount(int image_w, int image_h, int granularity,
                     int min_tile_size) {
  // int granularity = getNumThreads() * g_mpi_comm.size;

  tile_w = std::max(min_tile_size, image_w / granularity);
  tile_h = std::max(min_tile_size, image_h / granularity);

  int num_tiles_x = (image_w + tile_w - 1) / tile_w;
  int num_tiles_y = (image_h + tile_h - 1) / tile_h;
  num_tiles = num_tiles_x * num_tiles_y;

#ifndef NDEBUG
  printf("[INFO] granularity: %d\n", granularity);
  printf("[INFO] numer of tiles: %d (%d x %d tiles)\n", num_tiles, num_tiles_x,
         num_tiles_y);
#endif
}

}  // namespace spray
