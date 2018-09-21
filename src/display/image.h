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

#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>

#include "glm/glm.hpp"
#include "glog/logging.h"
#include "pbrt/memory.h"

#include "renderers/tile.h"
#include "utils/profiler_util.h"

namespace spray {

class TileImage {
 public:
  TileImage() : buf_(nullptr), w_(0), h_(0) {}
  ~TileImage() { FreeAligned(buf_); }

  void resize(int w, int h) {
    FreeAligned(buf_);
    int tile_area = w * h;
    w_ = w;
    h_ = h;
    buf_ = AllocAligned<glm::vec4>(tile_area);
  }
  void add(int pixid, const float rgb[3]) {
    // buf[pixid] += glm::vec4(rgb[0], rgb[1], rgb[2], 0.f);
    buf_[pixid].r += rgb[0];
    buf_[pixid].g += rgb[1];
    buf_[pixid].b += rgb[2];
    // buf[pixid].a = 0.000001f;
  }

 private:
  int w_, h_;
  glm::vec4* buf_;
};

struct HdrImage {
  HdrImage() : w(0), h(0), buf(nullptr) {}
  ~HdrImage() { FreeAligned(buf); }

  void resize(int width, int height) {
    FreeAligned(buf);
    buf = AllocAligned<glm::vec4>(width * height);
    w = width;
    h = height;
  }

  int bytes() const { return w * h * sizeof(glm::vec4); }

  void clear() { memset(buf, 0, w * h * sizeof(glm::vec4)); }
  void clear(std::size_t size) { memset(buf, 0, size * sizeof(glm::vec4)); }

  void update(int pixid, const float rgb[3], float scale) {
#ifdef SPRAY_GLOG_CHECK
    CHECK_LT(pixid, w * h);
#endif
    buf[pixid] =
        glm::vec4(scale * rgb[0], scale * rgb[1], scale * rgb[2], 0.000001f);
  }

  void update(int pixid, const float rgb[3]) {
    buf[pixid] = glm::vec4(rgb[0], rgb[1], rgb[2], 0.000001f);
  }

  void add(int pixid, const float rgb[3], float scale) {
#ifdef SPRAY_GLOG_CHECK
    CHECK_LT(pixid, w * h);
#endif
    buf[pixid].r += (scale * rgb[0]);
    buf[pixid].g += (scale * rgb[1]);
    buf[pixid].b += (scale * rgb[2]);
    // buf[pixid].a = 0.000001f;
  }
  void add(int pixid, const float rgb[3], double scale) {
#ifdef SPRAY_GLOG_CHECK
    CHECK_LT(pixid, w * h);
#endif
    // buf[pixid] +=
    //     glm::vec4(scale * rgb[0], scale * rgb[1], scale * rgb[2], 0.f);
    buf[pixid].r += (scale * rgb[0]);
    buf[pixid].g += (scale * rgb[1]);
    buf[pixid].b += (scale * rgb[2]);
    // buf[pixid].a = 0.000001f;
  }

  void add(int pixid, const float rgb[3]) {
#ifdef SPRAY_GLOG_CHECK
    CHECK_LT(pixid, w * h);
#endif
    // buf[pixid] += glm::vec4(rgb[0], rgb[1], rgb[2], 0.f);
    buf[pixid].r += rgb[0];
    buf[pixid].g += rgb[1];
    buf[pixid].b += rgb[2];
    // buf[pixid].a = 0.000001f;
  }

  void scale(const Tile& tile_in, int num_samples) {
    Tile tile = tile_in;
    int image_w = w;
    int pixid_offset, pixid;
    float scale = 1.f / num_samples;
    for (int y = tile.y; y < tile.y + tile.h; ++y) {
      pixid_offset = y * image_w;
      for (int x = tile.x; x < tile.x + tile.w; ++x) {
        pixid = pixid_offset + x;
        buf[pixid].r *= scale;
        buf[pixid].g *= scale;
        buf[pixid].b *= scale;
      }
    }
  }

  void parallelScale(const Tile& tile_in, int num_samples) {
    Tile tile = tile_in;
#pragma omp for collapse(2) schedule(static, 1)
    for (int y = tile.y; y < tile.y + tile.h; ++y) {
      for (int x = tile.x; x < tile.x + tile.w; ++x) {
        int image_w = w;
        float scale = 1.f / num_samples;
        int pixid_offset = y * image_w;
        int pixid = pixid_offset + x;
        buf[pixid].r *= scale;
        buf[pixid].g *= scale;
        buf[pixid].b *= scale;
      }
    }
  }

  void set(int pixid, const float rgb[3]) {
#ifdef SPRAY_GLOG_CHECK
    CHECK_LT(pixid, w * h);
#endif
    // buf[pixid] += glm::vec4(rgb[0], rgb[1], rgb[2], 0.f);
    buf[pixid].r = rgb[0];
    buf[pixid].g = rgb[1];
    buf[pixid].b = rgb[2];
    // buf[pixid].a = 0.000001f;
  }

  void set(int pixid, const glm::vec4& rgb) { buf[pixid] = rgb; }

  void composite() {
    int count = (w * h) << 2;
#ifdef SPRAY_TIMING
    tStartMPI(TIMER_SYNC_IMAGE);
#endif
    if (mpi::rank() == 0) {
      MPI_Reduce(MPI_IN_PLACE, buf, count, MPI_FLOAT, MPI_SUM, 0,
                 MPI_COMM_WORLD);
    } else {
      MPI_Reduce(buf, nullptr, count, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    }
#ifdef SPRAY_TIMING
    tStop(TIMER_SYNC_IMAGE);
#endif
  }

  void writePpm(const char* filename) {
    std::stringstream header;
    header << "P3" << std::endl;
    header << w << " " << h << std::endl;
    header << "1023" << std::endl;

    std::ofstream file;
    file.open(filename, std::ios::out | std::ios::trunc);
    file << header.str();

    for (int y = h - 1; y > -1; --y) {
      for (int x = 0; x < w; ++x) {
        int i = y * w + x;
        float r = glm::clamp(buf[i].r * 1023.f, 0.0f, 1023.0f);
        float g = glm::clamp(buf[i].g * 1023.f, 0.0f, 1023.0f);
        float b = glm::clamp(buf[i].b * 1023.f, 0.0f, 1023.0f);
        file << (unsigned)r << " " << (unsigned)g << " " << (unsigned)b
             << std::endl;
      }
    }
    file.close();
  }

  static void writePpm(int width, int height, const char* filename,
                       void* rgba) {
    std::stringstream header;
    header << "P3" << std::endl;
    header << width << " " << height << std::endl;
    header << "1023" << std::endl;

    std::ofstream file;
    file.open(filename, std::ios::out | std::ios::trunc);
    file << header.str();

    float* rgba_buf = (float*)rgba;

    // for (unsigned int y = 0; y < height; ++y) {
    for (int y = height - 1; y > -1; --y) {
      int offset = width * y;
      for (int x = 0; x < width; ++x) {
        int i = (offset + x) << 2;
        float r = glm::clamp(rgba_buf[i] * 1023.f, 0.f, 1023.f);
        float g = glm::clamp(rgba_buf[i + 1] * 1023.f, 0.f, 1023.f);
        float b = glm::clamp(rgba_buf[i + 2] * 1023.f, 0.f, 1023.f);
        file << (unsigned)r << " " << (unsigned)g << " " << (unsigned)b
             << std::endl;
      }
    }
    file.close();
  }

  int w;
  int h;
  glm::vec4* buf;
};

}  // namespace spray
