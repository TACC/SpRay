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

#include "render/tile.h"

#include <algorithm>
#include <cmath>

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

void BlockingTileList::init(int64_t image_w, int64_t image_h,
                            int64_t num_pixel_samples, int64_t num_ranks,
                            int64_t maximum_num_samples_per_rank) {
  //
  int64_t total_num_samples = image_w * image_h * num_pixel_samples;

  int64_t maximum_num_samples_per_cluster =
      maximum_num_samples_per_rank * num_ranks;

  int64_t total_num_tiles =
      (total_num_samples + maximum_num_samples_per_cluster - 1) /
      maximum_num_samples_per_cluster;

  CHECK_GT(total_num_tiles, 0);

  int64_t num_tiles_1d = static_cast<int64_t>(
      std::ceil(std::sqrt(static_cast<double>(total_num_tiles))));

  CHECK_GT(num_tiles_1d, 0);
  CHECK_LE(num_tiles_1d, image_w);
  CHECK_LE(num_tiles_1d, image_h);

  int64_t tile_w = image_w / num_tiles_1d;
  int64_t tile_h = image_h / num_tiles_1d;

  // int tile_w = std::max(min_tile_size_1d, image_w / num_tiles_1d);
  // int tile_h = std::max(min_tile_size_1d, image_h / num_tiles_1d);

  // CHECK_GT(tile_w, min_tile_size_1d);
  // CHECK_GT(tile_h, min_tile_size_1d);

  int64_t num_tiles_w = (image_w + tile_w - 1) / tile_w;
  int64_t num_tiles_h = (image_h + tile_h - 1) / tile_h;

  CHECK_GT(num_tiles_w, 0);
  CHECK_GT(num_tiles_h, 0);

  int64_t estimated_num_tiles = num_tiles_w * num_tiles_h;
  CHECK_GT(estimated_num_tiles, 0);

  tiles_.resize(estimated_num_tiles);

  Tile image;
  image.x = 0;
  image.y = 0;
  image.w = static_cast<int>(image_w);
  image.h = static_cast<int>(image_h);

  Tile tile_size;
  tile_size.x = 0;
  tile_size.y = 0;
  tile_size.w = static_cast<int>(tile_w);
  tile_size.h = static_cast<int>(tile_h);

  int tile_index = 0;
  int max_area = -1;
  int largest_tile_index;

  for (int y = 0; y < image.h; y += tile_size.h) {
    for (int x = 0; x < image.w; x += tile_size.w) {
      int w = std::min(tile_size.w, image.w - x);
      int h = std::min(tile_size.h, image.h - y);

      CHECK_LT(static_cast<int64_t>(tile_index), estimated_num_tiles);

      Tile& tile = tiles_[tile_index];
      tile.x = x;
      tile.y = y;
      tile.w = w;
      tile.h = h;

      CHECK_GE(tile.x, image.x);
      CHECK_LT(tile.x, image.x + image.w);
      CHECK_GE(tile.y, image.y);
      CHECK_LT(tile.y, image.y + image.h);
      CHECK_LE(tile.x + tile.w, image.x + image.w);
      CHECK_LE(tile.y + tile.h, image.y + image.h);

      int area = w * h;

      int64_t num_samples_per_tile =
          static_cast<int64_t>(area) * num_pixel_samples;
      CHECK_LE(num_samples_per_tile, maximum_num_samples_per_cluster);

      if (area > max_area) {
        max_area = area;
        largest_tile_index = tile_index;
      }

      ++tile_index;
    }
  }

  // initialize indices
  tile_index_ = 0;
  largest_tile_index_ = largest_tile_index;

  CHECK_EQ(static_cast<int64_t>(tile_index), estimated_num_tiles);
  CHECK_GT(largest_tile_index, -1);
  CHECK_LT(largest_tile_index, tiles_.size());

#ifdef SPRAY_GLOG_CHECK
  for (auto& t : tiles_) {
    if (t.getArea() > 0) {
      CHECK_GE(t.x, 0) << t;
      CHECK_LT(t.x, image_w) << t;
      CHECK_GE(t.y, 0) << t;
      CHECK_LT(t.y, image_h) << t;
      //
      CHECK_GT(t.x + t.w, 0) << t;
      CHECK_LE(t.x + t.w, image_w) << t;
      CHECK_GT(t.y + t.h, 0) << t;
      CHECK_LE(t.y + t.h, image_h) << t;
    }
  }
#endif
}

void TileList::init(int64_t image_w, int64_t image_h, int64_t num_pixel_samples,
                    int64_t num_ranks, int rank,
                    int64_t maximum_num_samples_per_rank) {
  blocking_tiles_.init(image_w, image_h, num_pixel_samples, num_ranks,
                       maximum_num_samples_per_rank);

  CHECK(!blocking_tiles_.empty());

  // set tiles_

  tiles_.resize(blocking_tiles_.size());

  int i = 0;

  while (!blocking_tiles_.empty()) {
    const Tile& blocking_tile = blocking_tiles_.front();
    blocking_tiles_.pop();

    CHECK_LT(i, blocking_tiles_.size());

    tiles_[i] = makeHorizontalStripe(num_ranks, rank, blocking_tile);
    ++i;
  }

#ifdef SPRAY_GLOG_CHECK
  for (auto& t : tiles_) {
    if (t.getArea() > 0) {
      CHECK_GE(t.x, 0) << t;
      CHECK_LT(t.x, image_w) << t;
      CHECK_GE(t.y, 0) << t;
      CHECK_LT(t.y, image_h) << t;
      //
      CHECK_GT(t.x + t.w, 0) << t;
      CHECK_LE(t.x + t.w, image_w) << t;
      CHECK_GT(t.y + t.h, 0) << t;
      CHECK_LE(t.y + t.h, image_h) << t;
    }
  }
#endif

  CHECK_EQ(i, blocking_tiles_.size());

  blocking_tiles_.reset();

  // set tile index
  tile_index_ = 0;
}

Tile makeHorizontalStripe(int num_ranks, int rank, const Tile& tile_in) {
  Tile tile_out;

  int h = std::max(tile_in.h / num_ranks, 1);
  tile_out.y = tile_in.y + (rank * h);

  int yend = tile_in.y + tile_in.h;

  if (tile_out.y >= yend) {  // invalid
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

Tile makeVerticalStripe(int num_ranks, int rank, const Tile& tile_in) {
  Tile tile_out;

  int w = std::max(tile_in.w / num_ranks, 1);
  tile_out.x = tile_in.x + (rank * w);

  int xend = tile_in.x + tile_in.w;

  if (tile_out.x >= xend) {  // invalid
    tile_out.w = 0;
    tile_out.h = 0;

  } else {  // valid
    tile_out.w = ((tile_out.x + w > xend) || (rank == num_ranks - 1))
                     ? xend - tile_out.x
                     : w;
    tile_out.y = tile_in.y;
    tile_out.h = tile_in.h;
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

void ImageScheduleTileList::init(int64_t image_w, int64_t image_h,
                                 int64_t num_pixel_samples, int64_t num_ranks,
                                 int rank,
                                 int64_t maximum_num_samples_per_rank) {
  Tile image;
  image.x = 0;
  image.y = 0;
  image.w = image_w;
  image.h = image_h;

  Tile vstripe = makeVerticalStripe(num_ranks, rank, image);
  CHECK_GT(vstripe.getArea(), 0);

  // evaluate the size of each blocking tile

  int64_t num_samples = static_cast<int64_t>(vstripe.w) *
                        static_cast<int64_t>(vstripe.h) * num_pixel_samples;

  int64_t estimated_num_tiles =
      (num_samples + maximum_num_samples_per_rank - 1) /
      maximum_num_samples_per_rank;

  CHECK_LT(estimated_num_tiles, INT_MAX);

  int ntiles = static_cast<int>(estimated_num_tiles);
  int tile_h = vstripe.h / ntiles;

  ntiles = (vstripe.h + tile_h - 1) / tile_h;

  CHECK_GT(ntiles, 0);
  tiles_.resize(ntiles);

  // make horizontal stripes
  std::size_t i = 0;

  int max_area = -1;
  largest_tile_index_ = -1;

  for (int y = 0; y < vstripe.h; y += tile_h) {
    int h = std::min(tile_h, vstripe.h - y);

    CHECK_LT(i, tiles_.size());

    Tile& t = tiles_[i];
    t.x = vstripe.x;
    t.y = y;
    t.w = vstripe.w;
    t.h = h;
#ifdef DEBUG_PRINT_TILES
    std::cout << "[" << i << "]" << t << "\n";
#endif
    int area = t.w * t.h;
    CHECK_GT(area, 0);
    if (area > max_area) {
      max_area = area;
      largest_tile_index_ = i;
    }

    ++i;
  }

  CHECK_GE(largest_tile_index_, 0);
  CHECK_GT(getLargestBlockingTile().getArea(), 0) << largest_tile_index_;

#ifdef SPRAY_GLOG_CHECK
  CHECK_EQ(i, tiles_.size());
  for (auto& t : tiles_) {
    CHECK_GT(t.w * t.h, 0);
    CHECK_EQ(t.x, vstripe.x);
    CHECK_GE(t.y, vstripe.y);
    CHECK_LT(t.y, vstripe.y + vstripe.h);
    CHECK_EQ(t.w, vstripe.w);
    CHECK_GT(t.y + t.h, vstripe.y);
    CHECK_LE(t.y + t.h, vstripe.y + vstripe.h);
  }
#endif

  tile_index_ = 0;
}

}  // namespace spray
