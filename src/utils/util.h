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

#include <getopt.h>
#include <glog/logging.h>
#include <time.h>
#include <cstdint>
#include <cstdlib>
#include <deque>
#include <glm/glm.hpp>
#include <iostream>
#include <queue>
#include <string>
#include <thread>

#include "render/spray.h"

#define DEFAULT_NUM_HARDWARE_THREADS 4

namespace spray {

namespace util {

//! Takes a commandline argument and its length and outputs parsed result.
//! Increments getopt's optind by the amount given by len.
inline void parseTuple(char** argv, int len, float* dout) {
  int index = optind - 1;
  for (int i = 0; i < len; ++i) {
    dout[i] = atof(argv[index + i]);
  }
  optind = index + len;
}

//! Returns a file extension.
inline std::string getFileExtension(const std::string& filename) {
  std::string ext = filename.substr(filename.find_last_of(".") + 1);
  return std::move(ext);
}

inline std::string getFilename(const std::string& s) {
  char sep = '/';

  std::size_t i = s.rfind(sep, s.length());
  if (i != std::string::npos) {
    return (s.substr(i + 1, s.length() - i));
  }
  return ("");
}

inline void unpack(uint32_t color, glm::vec3& rgb) {
  rgb[0] = ((color >> 16) & 0xff) * SPRAY_1_OVER_255;  // red
  rgb[1] = ((color >> 8) & 0xff) * SPRAY_1_OVER_255;   // green
  rgb[2] = (color & 0xff) * SPRAY_1_OVER_255;          // blue
}

inline void unpack(uint32_t color, float rgb[3]) {
  rgb[0] = ((color >> 16) & 0xff) * SPRAY_1_OVER_255;  // red
  rgb[1] = ((color >> 8) & 0xff) * SPRAY_1_OVER_255;   // green
  rgb[2] = (color & 0xff) * SPRAY_1_OVER_255;          // blue
}

inline void unpack(uint32_t color, uint32_t rgb[3]) {
  rgb[0] = (color >> 16) & 0xff;  // red
  rgb[1] = (color >> 8) & 0xff;   // green
  rgb[2] = color & 0xff;          // blue
}

inline void unpack(uint32_t color, glm::vec3* rgb) {
  (*rgb)[0] = ((color >> 16) & 0xff) * SPRAY_1_OVER_255;  // red
  (*rgb)[1] = ((color >> 8) & 0xff) * SPRAY_1_OVER_255;   // green
  (*rgb)[2] = (color & 0xff) * SPRAY_1_OVER_255;          // blue
}

inline uint32_t pack(uint32_t r, uint32_t g, uint32_t b) {
  CHECK_LT(r, 256);
  CHECK_LT(g, 256);
  CHECK_LT(b, 256);
  return ((r << 16) | (g << 8) | b);
}

inline uint32_t normalToColor(float n[3]) {
  glm::vec3 nn = glm::vec3(n[0], n[1], n[2]);
  nn = nn * 255.0f;
  uint32_t r = glm::abs(nn[0]) > 255.f ? 255 : glm::abs(nn[0]);
  uint32_t g = glm::abs(nn[1]) > 255.f ? 255 : glm::abs(nn[1]);
  uint32_t b = glm::abs(nn[2]) > 255.f ? 255 : glm::abs(nn[2]);
  return pack(r, g, b);
}

inline unsigned int getNumThreads() {
  unsigned int concurrency = std::thread::hardware_concurrency();
  return (concurrency == 0) ? DEFAULT_NUM_HARDWARE_THREADS : concurrency;
}

/**
 * A helper taking a string of bytes to evaluate the number of items in type T.
 *
 * \param base A pointer to the string.
 * \param bytes A string size in byte.
 * \return The number of items in type T.
 */
template <class T>
std::size_t getNumOfItems(void* base, std::size_t bytes) {
  uint8_t* end = ((uint8_t*)base) + bytes;
  return ((T*)end - (T*)base);
}

}  // namespace util
}  // namespace spray

