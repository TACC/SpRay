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

#include "GLFW/glfw3.h"
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "glm/gtc/type_ptr.hpp"

#include "display/composite.h"
#include "display/image.h"
#include "display/spray_glfw.h"
#include "display/vis.h"
#include "render/aabb.h"
#include "render/camera.h"
#include "render/config.h"
#include "render/domain.h"
#include "render/scene.h"
#include "render/spray.h"
#include "utils/comm.h"
#include "utils/profiler.h"
#include "utils/profiler_util.h"
#include "utils/timer.h"

namespace spray {

template <class TracerT, class SceneT>
class SprayRenderer {
 public:
  SprayRenderer() : cfg_(nullptr) {}

  void init(const Config& cfg);
  void run();

 private:
  void initCamera(const Config& cfg);

  void run_normal();
  void run_dev();
  void renderFilm();
  void renderGlfwSingleTask();
  void renderGlfwRootTask();
  void renderGlfwChildTask();
  void renderGlfwDomainBounds(int view_mode);

  void renderFilmInOmp();
  void renderGlfwInOmp();

 private:
  const Config* cfg_;

  MessageCommand msgcmd_;

  SceneT scene_;
  Camera camera_;
  TracerT tracer_;
  HdrImage image_;
};

}  // namespace spray

#define SPRAY_SPRAY_RENDERER_INL_
#include "render/spray_renderer.inl"
#undef SPRAY_SPRAY_RENDERER_INL_
