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

#if !defined(SPRAY_SPRAY_RENDERER_INL_)
#error An implementation of SprayRenderer
#endif

namespace spray {

template <class TracerT, class CacheT>
void SprayRenderer<TracerT, CacheT>::init(const Config &cfg) {
  // config
  cfg_ = &cfg;

  // scene
  bool insitu_mode = (cfg.partition == spray::Config::INSITU);

  scene_.init(cfg.model_descriptor_filename, cfg.ply_path, cfg.local_disk_path,
              cfg.cache_size, cfg.view_mode, insitu_mode, cfg.num_partitions);

#ifdef SPRAY_GLOG_CHECK
  LOG(INFO) << "scene init done";
#endif

#ifdef SPRAY_TIMING
  global_profiler.init();
#endif

  // build wbvh
  scene_.buildWbvh();

  // frame buffer
  image_.resize(cfg.image_w, cfg.image_h);

  // camera
  camera_.initialize(scene_.getBound(), cfg.image_w, cfg.image_h, cfg.znear,
                     cfg.zfar, cfg.fov);

  if (cfg.has_camera_config) {
    camera_.resetPosition(cfg.camera_pos, cfg.camera_lookat, cfg.camera_up);
  }

  // tracer
  if (!(cfg.view_mode == VIEW_MODE_DOMAIN ||
        cfg.view_mode == VIEW_MODE_PARTITION)) {
    tracer_.init(cfg, camera_, &scene_, &image_);
  }

  // vis
  WbvhObj<WbvhEmbree> wobj;
  wobj.ptr = nullptr;

  Vis<WbvhEmbree>::initialize(wobj);

  // mpi command
  msgcmd_.done = 0;
  msgcmd_.image_w = cfg.image_w;
  msgcmd_.image_h = cfg.image_h;
  msgcmd_.view_mode = cfg.view_mode;
  msgcmd_.camera_cmd = CAM_NOP;

  // glfw
  if (cfg.view_mode != VIEW_MODE_FILM) {
    Glfw<WbvhEmbree, CacheT>::initialize(cfg, mpi::isRootProcess(), cfg.image_w,
                                         cfg.image_h, &camera_, &msgcmd_,
                                         &scene_);
  }

#ifdef SPRAY_GLOG_CHECK
  int ndomains = scene_.getNumDomains();
  LOG_IF(INFO, mpi::rank() == 0) << "number of domains: " << ndomains;
#endif
}

template <class TracerT, class CacheT>
void SprayRenderer<TracerT, CacheT>::run() {
  CHECK(cfg_) << "failed to run. not configured.";

  if (msgcmd_.view_mode == VIEW_MODE_FILM) {
    renderFilm();
  } else if (msgcmd_.view_mode == VIEW_MODE_GLFW) {
    if (mpi::isSingleProcess()) {
      renderGlfwSingleTask();
    } else if (mpi::isRootProcess()) {
      renderGlfwRootTask();
    } else {
      renderGlfwChildTask();
    }
  } else {
    msgcmd_.view_mode = VIEW_MODE_TERMINATE;
    glfwTerminate();
    LOG(FATAL) << "unsupported mode " << msgcmd_.view_mode;
  }
}

template <class TracerT, class CacheT>
void SprayRenderer<TracerT, CacheT>::renderGlfwSingleTask() {
#ifdef SPRAY_TIMING
  tReset();
  tStartMPI(TIMER_TOTAL);
#endif

  int64_t cfg_nframes = cfg_->nframes;
  int image_w = image_.w;
  int image_h = image_.h;
  glm::vec4 *image_buf = image_.buf;

  CHECK_EQ(msgcmd_.view_mode, VIEW_MODE_GLFW);

  int64_t nframes = 0;

  while (nframes < cfg_nframes || (cfg_nframes < 0 && !msgcmd_.done)) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    tracer_.trace();

    glDrawPixels(image_w, image_h, GL_RGBA, GL_FLOAT, image_buf);
    Glfw<WbvhEmbree, CacheT>::swapBuffers();
    glfwPollEvents();
    Glfw<WbvhEmbree, CacheT>::cmdHandler();

    ++nframes;
  }

#ifdef SPRAY_TIMING
  tStop(TIMER_TOTAL);
#endif

  glfwTerminate();

#ifdef SPRAY_TIMING
  tPrint(cfg_nframes);
#endif
}

template <class TracerT, class CacheT>
void SprayRenderer<TracerT, CacheT>::renderGlfwRootTask() {
#ifdef SPRAY_TIMING
  tReset();
  tStartMPI(TIMER_TOTAL);
#endif

  int64_t cfg_nframes = cfg_->nframes;
  int image_w = image_.w;
  int image_h = image_.h;
  glm::vec4 *image_buf = image_.buf;

  CHECK_EQ(msgcmd_.view_mode, VIEW_MODE_GLFW);

  int64_t nframes = 0;

  while (nframes < cfg_nframes || (cfg_nframes < 0 && !msgcmd_.done)) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    image_.clear();
    tracer_.trace();
    image_.composite();

    glDrawPixels(image_w, image_h, GL_RGBA, GL_FLOAT, image_buf);

    Glfw<WbvhEmbree, CacheT>::swapBuffers();
    glfwPollEvents();

    // send command
    MPI_Bcast((void *)&msgcmd_, sizeof(MessageCommand), MPI_UNSIGNED_CHAR, 0,
              MPI_COMM_WORLD);

    // handle command
    Glfw<WbvhEmbree, CacheT>::cmdHandler();

    ++nframes;
  }

#ifdef SPRAY_TIMING
  tStop(TIMER_TOTAL);
#endif

  glfwTerminate();
  msgcmd_.view_mode = VIEW_MODE_TERMINATE;

#ifdef SPRAY_TIMING
  tPrint(cfg_nframes);
#endif
}

template <class TracerT, class CacheT>
void SprayRenderer<TracerT, CacheT>::renderGlfwChildTask() {
#ifdef SPRAY_TIMING
  tReset();
  tStartMPI(TIMER_TOTAL);
#endif

  int64_t cfg_nframes = cfg_->nframes;
  int64_t nframes = 0;

  while (nframes < cfg_nframes || (cfg_nframes < 0 && !msgcmd_.done)) {
    image_.clear();
    tracer_.trace();
    image_.composite();

    // recv command
    MPI_Bcast((void *)&msgcmd_, sizeof(MessageCommand), MPI_UNSIGNED_CHAR, 0,
              MPI_COMM_WORLD);

    // handle command
    Glfw<WbvhEmbree, CacheT>::cmdHandler();

    ++nframes;
  }

#ifdef SPRAY_TIMING
  tStop(TIMER_TOTAL);
  tPrint(cfg_nframes);
#endif

  msgcmd_.view_mode = VIEW_MODE_TERMINATE;
}

template <class TracerT, class CacheT>
void SprayRenderer<TracerT, CacheT>::renderFilm() {
#if defined(SPRAY_TIMING)
  tReset();
  tStartMPI(TIMER_TOTAL);
#endif

  bool cluster = (mpi::size() > 1);
  int64_t cfg_nframes = cfg_->nframes;

  for (int64_t i = 0; i < cfg_nframes; ++i) {
    image_.clear();
    tracer_.trace();

    if (cluster) {
      image_.composite();
    }
  }

#ifdef SPRAY_TIMING
  tStop(TIMER_TOTAL);
#endif

  bool root = (mpi::rank() == 0);

  if (!cluster || root) {
    image_.writePpm(cfg_->output_filename.c_str());
  }

#ifdef SPRAY_TIMING
  tPrint(cfg_nframes);
#endif
}

}  // namespace spray
