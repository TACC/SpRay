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

#include "render/config.h"

#include <getopt.h>

#include "glog/logging.h"

#include "render/light.h"
#include "render/spray.h"
#include "utils/util.h"

namespace spray {

Config::Config() {
  // image
  image_w = 400;
  image_h = 400;

  // camera
  has_camera_config = false;
  znear = 0.000001f;
  zfar = 1000000.0f;
  fov = 90.0f;

  // render
  nframes = -1;
  output_filename = "spray.ppm";
  light_samples = 1;
  bounces = 1;

  // schedule
  partition = IMAGE;
  num_partitions = 1;

  // visualization
  view_mode = VIEW_MODE_GLFW;

  cache_size = -1;

  // ao settings
  ao_samples = 8;
  ao_mode = 0;

  // pt settings
  pixel_samples = 4;

  num_tiles = 1;
  min_tile_size = 64;

  nthreads = 1;

  shading = SPRAY_SHADING_LAMBERT;

  dev_mode = DEVMODE_NORMAL;
}

void Config::printUsage(char** argv) {
  printf("Usage: %s [options] <modelfile>\n", argv[0]);
  printf("Options:\n");
  printf("  (): default setting\n");
  printf("  --output, -o <filename (spray.ppm)>\n");
  printf("     only effective in film mode\n");
  printf("  --local-disk, -l <path to local disk>\n");
  printf("  --ply-path <path to ply files>\n");
  printf("  --mode, -m <film | glfw | domain | partition>\n");
  printf(
      "  --num-partitions <number of partitions>, effective in partition view "
      "mode\n");
  printf("  --partition <image | hybrid | insitu>\n");
  printf("  --cache-size <max. number of domains>\n");
  printf("  --width, -w <image_width>\n");
  printf("  --height, -h <image_height>\n");
  printf("  --frames <number of frames (-1)>\n");
  printf("  --light-samples, <number of light samples (1)>\n");
  printf("  --bounces, <number of bounces (1)>\n");
  printf("  --camera-up <upx upy upz>\n");
  printf("  --camera <posx posy posz lookx looky lookz>\n");
  printf("  --ao-samples <number of samples in AO (8)>\n");
  printf("  --pixel-samples <number of pixel samples (4)>\n");
  printf("  --nthreads <number of threads (1)>\n");
  printf("  --shading <lambert | blinn>\n");
  printf("  --blinn ks_r ks_g ks_b shininess\n");
  printf("  --dev-mode\n");
}

void Config::parse(int argc, char** argv) {
  struct option long_options[] = {
      {"output", required_argument, 0, 'o'},
      {"local-disk", required_argument, 0, 'l'},
      {"mode", required_argument, 0, 'm'},
      {"width", required_argument, 0, 'w'},
      {"height", required_argument, 0, 'h'},
      {"frames", required_argument, 0, 300},
      {"nthreads", required_argument, 0, 302},
      {"camera-up", required_argument, 0, 308},
      {"light-samples", required_argument, 0, 309},
      {"bounces", required_argument, 0, 310},
      {"camera", required_argument, 0, 311},
      {"cache-size", required_argument, 0, 313},
      {"fov", required_argument, 0, 314},
      {"partition", required_argument, 0, 324},
      {"num-partitions", required_argument, 0, 325},
      {"ao-samples", required_argument, 0, 400},
      {"ao-mode", no_argument, 0, 402},
      {"pixel-samples", required_argument, 0, 403},
      {"shading", required_argument, 0, 404},
      {"blinn", required_argument, 0, 405},
      {"num-tiles", required_argument, 0, 406},
      {"min-tile-size", required_argument, 0, 407},
      {"ply-path", required_argument, 0, 408},
      {"dev-mode", no_argument, 0, 1000},
      {0, 0, 0, 0}};

  int option_index = 0;
  int c;
  camera_up[0] = 0.f;
  camera_up[1] = 1.f;
  camera_up[2] = 0.f;

  while ((c = getopt_long(argc, argv, "o:t:m:r:w:h:x", long_options,
                          &option_index)) != -1) {
    switch (c) {
      case 'o':
        output_filename = optarg;
        break;
      case 'l':
        local_disk_path = optarg;
        break;
      case 'm': {
        std::string mode = optarg;
        if (mode == "domain") {
          view_mode = VIEW_MODE_DOMAIN;
        } else if (mode == "partition") {
          view_mode = VIEW_MODE_PARTITION;
        } else if (mode == "film") {
          view_mode = VIEW_MODE_FILM;
        } else if (mode == "glfw") {
          view_mode = VIEW_MODE_GLFW;
        } else {
          LOG(FATAL) << "undefined view mode" << optarg;
        }
      } break;
      case 'w':
        image_w = atoi(optarg);
        break;
      case 'h':
        image_h = atoi(optarg);
        break;

      case 300: {  // --frames
        nframes = atoi(optarg);
      } break;

      case 302: {  // --nthreads
        nthreads = atoi(optarg);
      } break;

      case 308: {  // --camera-upvector
        float data[3];
        util::parseTuple(argv, 3, data);
        camera_up = glm::vec3(data[0], data[1], data[2]);
      } break;

      case 309: {  // --light-samples
        light_samples = atoi(optarg);
      } break;

      case 310: {  // --bounces
        bounces = atoi(optarg);
      } break;

      case 311: {  // --camera
        has_camera_config = true;
        float data[6];
        util::parseTuple(argv, 6, data);
        camera_pos = glm::vec3(data[0], data[1], data[2]);
        camera_lookat = glm::vec3(data[3], data[4], data[5]);
      } break;

      case 313: {  // --cache-size
        cache_size = atoi(optarg);
      } break;
      case 314: {
        fov = atof(optarg);
      } break;
      case 324: {
        std::string cfg_partition = optarg;
        if (cfg_partition == "image") {
          partition = IMAGE;

        } else if (cfg_partition == "hybrid") {
          partition = HYBRID;

        } else if (cfg_partition == "insitu") {
          partition = INSITU;

        } else {
          LOG(FATAL) << "unsupported partitioning type: " << partition;
        }
      } break;

      case 325: { // --num-partitions
        num_partitions = atoi(optarg);
      } break;

      case 400: {  // --ao-samples
        ao_samples = atoi(optarg);
      } break;

      case 402: {  // --ao-mode
        ao_mode = 1;
      } break;

      case 403: {  // --pixel-samples
        pixel_samples = atoi(optarg);
      } break;

      case 404: {  // --shading
        std::string cfg_shading = optarg;
        if (cfg_shading == "blinn") {
          shading = SPRAY_SHADING_BLINN;
        } else {
          shading = SPRAY_SHADING_LAMBERT;
        }
      } break;

      case 405: {  // --blinn
        float data[4];
        util::parseTuple(argv, 4, data);
        ks = glm::vec3(data[0], data[1], data[2]);
        shininess = data[3];
      } break;

      case 406: {  // --num-tiles
        num_tiles = atoi(optarg);
      } break;

      case 407: {  // --min-tile-size
        min_tile_size = atoi(optarg);
      } break;

      case 408: {  // --ply-path
        ply_path = optarg;
      } break;

      case 1000: {  // --dev-mode
        dev_mode = DEVMODE_DEV;
      } break;


      default:
        printUsage(argv);
        LOG(FATAL) << "undefined arg " << c;
        break;
    }
  }  // end of getopt

  CHECK_NE(argc, optind) << "input file not found";

  std::string filename(argv[optind]);
  model_descriptor_filename = filename;

  std::string ext = util::getFileExtension(filename);
  CHECK_EQ(ext, std::string("spray"));
}

}  // namespace spray
