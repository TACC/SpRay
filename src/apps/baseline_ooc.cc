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

#include "glog/logging.h"

#include "baseline/baseline_image_tracer.h"
#include "baseline/baseline_schedulers.h"
#include "baseline/baseline_shader_ao.h"
#include "baseline/baseline_shader_pt.h"
#include "caches/caches.h"
#include "partition/data_partition.h"
#include "renderers/config.h"
#include "renderers/spray.h"
#include "renderers/spray_renderer.h"
#include "utils/comm.h"

int main(int argc, char** argv) {
  int required = MPI_THREAD_FUNNELED;
  int provided;
  MPI_Init_thread(&argc, &argv, required, &provided);

  MPI_Comm_size(MPI_COMM_WORLD, &(spray::global_mpi_comm.size));
  MPI_Comm_rank(MPI_COMM_WORLD, &(spray::global_mpi_comm.rank));

  google::InitGoogleLogging(argv[0]);
#ifdef SPRAY_GLOG_CHECK
  google::InstallFailureSignalHandler();
#endif

  CHECK_EQ(provided, required) << "MPI_THREAD_FUNNELED not available.";

#ifdef SPRAY_GLOG_CHECK
  LOG(INFO) << "rank " << spray::mpi::worldRank()
            << " (world size: " << spray::mpi::worldSize() << ")";
#endif

  spray::Config cfg;
  cfg.parse(argc, argv);

  if (cfg.partition == spray::Config::IMAGE) {
    if (cfg.ao_mode) {
      if (cfg.cache_size < 0) {
        spray::SprayRenderer<
            spray::baseline::ImageTracer<
                spray::InfiniteCache, spray::baseline::LoadAnyOnceImageSched,
                spray::baseline::ShaderAo<spray::InfiniteCache>>,
            spray::InfiniteCache>
            render;
        render.init(cfg);
        render.run();
      } else {
        spray::SprayRenderer<
            spray::baseline::ImageTracer<
                spray::LruCache, spray::baseline::LoadAnyOnceImageSched,
                spray::baseline::ShaderAo<spray::LruCache>>,
            spray::LruCache>
            render;
        render.init(cfg);
        render.run();
      }
    } else {
      if (cfg.cache_size < 0) {
        spray::SprayRenderer<
            spray::baseline::ImageTracer<
                spray::InfiniteCache, spray::baseline::LoadAnyOnceImageSched,
                spray::baseline::ShaderPt<spray::InfiniteCache>>,
            spray::InfiniteCache>
            render;
        render.init(cfg);
        render.run();
      } else {
        spray::SprayRenderer<
            spray::baseline::ImageTracer<
                spray::LruCache, spray::baseline::LoadAnyOnceImageSched,
                spray::baseline::ShaderPt<spray::LruCache>>,
            spray::LruCache>
            render;
        render.init(cfg);
        render.run();
      }
    }
  } else {
    LOG(FATAL) << "unsupported partition " << cfg.partition;
  }

  MPI_Finalize();
  return 0;
}

