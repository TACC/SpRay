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

#include "renderers/spray.h"
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

}  // namespace spray

