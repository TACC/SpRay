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

#include "insitu/insitu_tiler.h"

#include <algorithm>
#include <vector>

#include "glog/logging.h"

namespace spray {
namespace insitu {

std::ostream& operator<<(std::ostream& os, const Tile& t) {
  os << "tile(" << t.x << "," << t.y << "," << t.w << "," << t.h << ")";
  return os;
}
void Tiler::resize(int image_w, int image_h, int num_tiles_1d,
                   int min_tile_size_1d) {
  //
  int tile_w = std::max(min_tile_size_1d, image_w / num_tiles_1d);
  int tile_h = std::max(min_tile_size_1d, image_h / num_tiles_1d);

  CHECK_GT(tile_w, min_tile_size_1d);
  CHECK_GT(tile_h, min_tile_size_1d);

  int ntiles_w = (image_w + tile_w - 1) / tile_w;
  int ntiles_h = (image_h + tile_h - 1) / tile_h;

  CHECK_GT(ntiles_w, 0);
  CHECK_GT(ntiles_h, 0);

  int ntiles = ntiles_w * ntiles_h;
  tiles_.resize(ntiles);

  int max_area = -1;
  int max_tile_id;

  int i = 0;
  for (int y = 0; y < image_h; y += tile_h) {
    for (int x = 0; x < image_w; x += tile_w) {
      int w = std::min(tile_w, image_w - x);
      int h = std::min(tile_h, image_h - y);

      CHECK_LT(i, ntiles);

      Tile& tile = tiles_[i];
      // tile.id = i;
      tile.x = x;
      tile.y = y;
      tile.w = w;
      tile.h = h;

      CHECK_GT(tile.w, 0);
      CHECK_GT(tile.h, 0);

      int area = w * h;
      if (area > max_area) {
        max_area = area;
        max_tile_id = i;
      }

      ++i;
    }
  }
  CHECK_EQ(ntiles, i);
  CHECK_GT(ntiles, 0);
  CHECK_GT(max_tile_id, -1);

  id_ = 0;
  max_tile_id_ = max_tile_id;

  CHECK_GT(ntiles, 0);

  for (const auto& t : tiles_) {
    CHECK_GT(t.w * t.h, 0);
  }
}

Tile RankStriper::make(int num_ranks, int rank, const Tile& tile_in) {
  Tile tile_out;

  int h = std::max(tile_in.h / num_ranks, 1);
  tile_out.y = tile_in.y + (rank * h);

  int yend = tile_in.y + tile_in.h;

  if (tile_out.y >= yend) {  // invalid
    // tile_out.id = -1;
    tile_out.w = 0;
    tile_out.h = 0;

  } else {  // valid
    tile_out.h = ((tile_out.y + h > yend) || (rank == num_ranks - 1))
                     ? yend - tile_out.y
                     : h;
    tile_out.x = tile_in.x;
    tile_out.w = tile_in.w;
  }

  return tile_out;
}

void RankTiler::resize(const Tile& max_tile, int min_tile_size_1d) {
  // note: rough estimation
  int len = (max_tile.w > max_tile.h) ? max_tile.w : max_tile.h;
  int ntiles_1d = (len / min_tile_size_1d) + 1;

  int ntiles = ntiles_1d * ntiles_1d;

  CHECK_GT(ntiles, 0);

  tiles_.resize(ntiles);
}

void RankTiler::make(int num_ranks, int rank, Tile tile_in, int num_tiles_1d,
                     int min_tile_size_1d) {
  int tile_w = std::max(min_tile_size_1d, tile_in.w / num_tiles_1d);
  int tile_h = std::max(min_tile_size_1d, tile_in.h / num_tiles_1d);

  CHECK_GT(tile_w, min_tile_size_1d);
  CHECK_GT(tile_h, min_tile_size_1d);

  int ntiles_w = (tile_in.w + tile_w - 1) / tile_w;
  int ntiles_h = (tile_in.h + tile_h - 1) / tile_h;

  int ntiles = ntiles_w * ntiles_h;
  CHECK_GT(ntiles, 0);
  CHECK_LE(ntiles, tiles_.size());

  // distribute across cluster
  auto share = std::max(ntiles / num_ranks, 1);
  begin_ = share * rank;
  if (begin_ >= ntiles) {
    begin_ = ntiles;
    end_ = ntiles;
  } else {
    end_ = begin_ + share;
    if (end_ > ntiles) {
      end_ = ntiles;
    }
  }

  CHECK_LE(begin_, end_);

  id_ = begin_;

  int yend = tile_in.y + tile_in.h;
  int xend = tile_in.x + tile_in.w;
  int i = 0;

  for (int y = tile_in.y; y < yend; y += tile_h) {
    for (int x = tile_in.x; x < xend; x += tile_w) {
      //
      int w = std::min(tile_w, xend - x);
      int h = std::min(tile_h, yend - y);

      Tile& tile = tiles_[i];
      // tile.id = i;
      tile.x = x;
      tile.y = y;
      tile.w = w;
      tile.h = h;

      ++i;
    }
  }

  CHECK_EQ(ntiles, i);
  CHECK_GT(ntiles, 0);

  for (const auto& t : tiles_) {
    CHECK_GT(t.w * t.h, 0);
  }

  int area = 0;
  for (int i = begin_; i < end_; ++i) {
    auto& t = tiles_[i];
    area += (t.w * t.h);
  }
  total_area_ = area;
}

}  // namespace insitu
}  // namespace spray
