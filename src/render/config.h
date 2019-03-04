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

#include <string>
#include <vector>

#include "glm/glm.hpp"

#include "render/spray.h"

namespace spray {

class Light;

class Config {
  void printUsage(char** argv);

 public:
  Config();

  bool parse(int argc, char** argv);

  // image
  int image_w;
  int image_h;

  // model
  std::string model_descriptor_filename;
  std::string ply_path;

  // camera
  bool has_camera_config;
  glm::vec3 camera_pos;
  glm::vec3 camera_lookat;
  glm::vec3 camera_up;
  float znear;
  float zfar;
  float fov;

  // render
  int nframes;
  std::string output_filename;
  int light_samples;
  int bounces;

  // schedule
  enum Partition { IMAGE, HYBRID, INSITU };
  int partition;
  int num_partitions;  // effective when VIEW_MODE_PARTITION used

  // view mode
  ViewMode view_mode;

  // cache
  int cache_size;

  // ao settings
  int ao_samples;
  int ao_mode;

  // pt settings
  int pixel_samples;

  std::size_t maximum_num_screen_space_samples_per_rank;

  std::string local_disk_path;
  int nthreads;

  bool use_spray_color;

  enum DevMode { DEVMODE_NORMAL, DEVMODE_DEV };
  int dev_mode;

  glm::vec3 bg_color;
};

}  // namespace spray
