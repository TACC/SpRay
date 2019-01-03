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

#include "glog/logging.h"

#include "render/spray.h"
#include "utils/comm.h"
#include "utils/util.h"

#define SPRAY_MIN_TILE_SIZE 8

namespace spray {

struct Tile;

struct TileCount {
  TileCount(int image_w, int image_h,
            int granularity = SPRAY_MIN_TILE_SIZE /* the higher the finer */,
            int min_tile_size = SPRAY_MIN_TILE_SIZE /* 8^2 pixels */);
  int tile_w;
  int tile_h;
  int num_tiles;
};

struct Tile {
  Tile() : x(0), y(0), w(0), h(0) {}
  Tile(int xpos, int ypos, int width, int height)
      : x(xpos), y(ypos), w(width), h(height) {}
  int getArea() const { return w * h; }
  bool isValid() const { return (w != 0); }
  void invalidate() { w = 0; }

  friend std::ostream& operator<<(std::ostream& os, const Tile& t);

  int x, y, w, h;
};

struct TheOwner {
  bool operator()(int tile_id) { return true; }
};

struct SmemTileGrain {
  int operator()() { return util::getNumThreads(); }
};

struct DmemTileGrain {
  int operator()() {
    // int grain = util::getNumThreads() * g_mpi_comm.size;
    int grain = mpi::size();
    if (grain > 32) grain = 32;
    return grain;
  }
};

template <class TileGrainFunctor = SmemTileGrain,
          class IsOwnerFunctor = TheOwner>
std::queue<Tile> makeTiles(int image_w, int image_h,
                           TileGrainFunctor grain_size = SmemTileGrain(),
                           IsOwnerFunctor is_owner = TheOwner()) {
  TileCount tile_count(image_w, image_h, grain_size() /*granularity*/,
                       SPRAY_MIN_TILE_SIZE /*min_tile_size*/);

  std::queue<Tile> tiles;

  // int sum = 0;
  int tile_id = 0;
  Tile tile;
  for (int y = 0; y < image_h; y += tile_count.tile_h) {
    for (int x = 0; x < image_w; x += tile_count.tile_w) {
      if (is_owner(tile_id)) {
        int w = std::min(tile_count.tile_w, image_w - x);
        int h = std::min(tile_count.tile_h, image_h - y);

        tile.x = x;
        tile.y = y;
        tile.w = w;
        tile.h = h;

        tiles.push(tile);
      }
      ++tile_id;
      // sum+=tile.getArea();
    }
  }
  // LOG_VALUE(DBG, "tile area", sum);
  return (std::move(tiles));
}

template <class TileGrainFunctor = SmemTileGrain,
          class IsOwnerFunctor = TheOwner>
void makeTiles(int image_w, int image_h, std::queue<Tile>& tiles,
               int granularity, int min_tile_len) {
  //
  TileCount tile_count(image_w, image_h, granularity, min_tile_len);

  Tile tile;
  for (int y = 0; y < image_h; y += tile_count.tile_h) {
    for (int x = 0; x < image_w; x += tile_count.tile_w) {
      int w = std::min(tile_count.tile_w, image_w - x);
      int h = std::min(tile_count.tile_h, image_h - y);

      tile.x = x;
      tile.y = y;
      tile.w = w;
      tile.h = h;

      tiles.push(tile);

#ifdef SPRAY_GLOG_CHECK
      CHECK(tile.getArea() > 0) << "invalid tile generated";
#endif
    }
  }
}

template <class TileGrainFunctor = SmemTileGrain,
          class IsOwnerFunctor = TheOwner>
std::vector<Tile> makeTileVector(int image_w, int image_h,
                                 TileGrainFunctor grain_size = SmemTileGrain(),
                                 IsOwnerFunctor is_owner = TheOwner()) {
  TileCount tile_count(image_w, image_h, grain_size() /*granularity*/,
                       SPRAY_MIN_TILE_SIZE /*min_tile_size*/);

  std::vector<Tile> tiles;
  tiles.reserve(tile_count.num_tiles);

  // int sum = 0;
  int tile_id = 0;
  Tile tile;
  for (int y = 0; y < image_h; y += tile_count.tile_h) {
    for (int x = 0; x < image_w; x += tile_count.tile_w) {
      if (is_owner(tile_id)) {
        int w = std::min(tile_count.tile_w, image_w - x);
        int h = std::min(tile_count.tile_h, image_h - y);

        tile.x = x;
        tile.y = y;
        tile.w = w;
        tile.h = h;

        tiles.push_back(tile);
      }
      ++tile_id;
      // sum+=tile.getArea();
    }
  }
  // LOG_VALUE(DBG, "tile area", sum);
  return (std::move(tiles));
}

/**
 * Given the number of ranks, create a horizontal stripe of the input tile for
 * the given rank.
 *
 * \param num_ranks Number of ranks in the cluster.
 * \param rank Target rank.
 * \param tile_in Input tile to get a horizontal stripe from.
 * \return Horizontal stripe for the given rank. The stripe can be either valid or
 * invalid.
 */
Tile makeHorizontalStripe(int num_ranks, int rank, const Tile& tile_in);

class BlockingTileList {
 public:
  BlockingTileList() : tile_index_(-1), largest_tile_index_(-1) {}

  void init(int64_t image_w, int64_t image_h, int64_t num_pixel_samples,
            int64_t num_ranks, int64_t maximum_num_samples_per_rank);

  void reset() { tile_index_ = 0; }

  const Tile& front() const {
#ifdef SPRAY_GLOG_CHECK
    CHECK_LT(tile_index_, tiles_.size());
#endif
    return tiles_[tile_index_];
  }

  void pop() {
#ifdef SPRAY_GLOG_CHECK
    CHECK_LT(tile_index_, tiles_.size());
#endif
    ++tile_index_;
  }

  bool empty() const { return (tile_index_ == tiles_.size()); }

  std::size_t size() const { return tiles_.size(); }

  const Tile& getLargestTile() const { return tiles_[largest_tile_index_]; }

 private:
  std::vector<Tile> tiles_;
  int tile_index_;
  int largest_tile_index_;
};

class TileList {
 public:
  void init(int64_t image_w, int64_t image_h, int64_t num_pixel_samples,
            int64_t num_ranks, int rank, int64_t maximum_num_samples_per_rank);

  void reset() {
    tile_index_ = 0;
    blocking_tiles_.reset();
  }

  void front(Tile* blocking_tile, Tile* stripe) const {
#ifdef SPRAY_GLOG_CHECK
    CHECK_LT(tile_index_, tiles_.size());
#endif
    *blocking_tile = blocking_tiles_.front();
    *stripe = tiles_[tile_index_];
  }

  void pop() {
#ifdef SPRAY_GLOG_CHECK
    CHECK_LT(tile_index_, tiles_.size());
#endif
    blocking_tiles_.pop();
    ++tile_index_;
  }

  bool empty() const { return (tile_index_ == tiles_.size()); }

  std::size_t size() const { return tiles_.size(); }

  const Tile& getLargestBlockingTile() const {
    return blocking_tiles_.getLargestTile();
  }

 private:
  BlockingTileList blocking_tiles_;

  std::vector<Tile> tiles_;
  int tile_index_;
};

class RankTiler {
 public:
  RankTiler() : begin_(-1), end_(-1) {}

  void resize(const Tile& max_tile, int min_tile_size_1d);

  void make(int num_ranks, int rank, Tile tile_in, int num_tiles_1d,
            int min_tile_size_1d);

  Tile& front() { return tiles_[id_]; }
  const Tile& front() const { return tiles_[id_]; }

  int size() const { return (end_ - begin_); }
  int totalArea() const { return total_area_; }

  const Tile& getTile(int i) const { return tiles_[i]; }
  int begin() const { return begin_; }
  int end() const { return end_; }

 private:
  std::vector<Tile> tiles_;

  int id_;
  int begin_;
  int end_;
  int total_area_;
};

/**
 * This tile list divides the image plane so that it first creates as many tiles
 * (i.e. vertical stripes) as the number of ranks, and then each rank further
 * divides assigned tiles to create blocking tiles (i.e. horizontal stripes)
 * based on the maximum number of samples allowed.
 */
class ImageScheduleTileList {
 public:
  void init(int64_t image_w, int64_t image_h, int64_t num_pixel_samples,
            int64_t num_ranks, int rank, int64_t maximum_num_samples_per_rank) {
  }

  void reset() { tile_index_ = 0; }

 private:
  std::vector<Tile> tiles_;
  int tile_index_;
};

}  // namespace spray

