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

#include <vector>
#include <iostream>

namespace spray {
namespace ooc {

struct Tile {
  int id;
  int x;
  int y;
  int w;
  int h;
  friend std::ostream& operator<<(std::ostream& os, const Tile& t);
};

class Tiler {
 public:
  Tiler() : id_(-1), max_tile_id_(-1) {}

  void resize(int image_w, int image_h, int num_tiles_1d, int min_tile_size_1d);
  void reset() { id_ = 0; }

  Tile& front() { return tiles_[id_]; }
  const Tile& front() const { return tiles_[id_]; }

  void pop() { ++id_; }

  bool empty() const { return (id_ == tiles_.size()); }

  std::size_t size() const { return tiles_.size(); }

  const Tile& getMaxTile() const { return tiles_[max_tile_id_]; }

 private:
  std::vector<Tile> tiles_;
  int id_;
  int max_tile_id_;
};

struct RankStriper {
  static Tile make(int num_ranks, int rank, const Tile& tile_in);
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

}  // namespace ooc
}  // namespace spray
