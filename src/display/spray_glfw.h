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
#include "glog/logging.h"

#include "display/vis.h"
#include "render/camera.h"
#include "render/config.h"
#include "render/spray.h"

namespace spray {

struct MouseState {
  bool left;
  bool right;
  bool middle;
  int x;
  int y;
};

template <class WbvhT, class SceneT>
class Glfw {
 public:
  static void init(const Config& cfg, bool is_root_process, unsigned image_w,
                   unsigned image_h, Camera* camera, MessageCommand* cmd,
                   SceneT* scene);

  static void setShouldClose() {
    glfwSetWindowShouldClose(glfw_window_, GL_TRUE);
  }
  static bool shouldClose() { return glfwWindowShouldClose(glfw_window_); }
  static void swapBuffers() { glfwSwapBuffers(glfw_window_); }

  static void keyCallback(GLFWwindow* window, int key, int scancode, int action,
                          int mods);

  static void closeCallback(GLFWwindow* window);

  static void scrollCallback(GLFWwindow* window, double xoffset,
                             double yoffset);

  static void mouseButtonCallback(GLFWwindow* window, int button, int action,
                                  int mods);

  static void cursorPosCallback(GLFWwindow* window, double xpos, double ypos);

  static void cmdHandler();

 private:
  static void zoom(float offset) {
    msgcmd_->camera_cmd = CAM_ZOOM;
    msgcmd_->zoom_offset = offset;
  }
  static void rotate(float dx, float dy) {
    msgcmd_->camera_cmd = CAM_ROTATE;
    msgcmd_->rotate_pan_dx = dx;
    msgcmd_->rotate_pan_dy = dy;
  }

  static void pan(float dx, float dy) {
    msgcmd_->camera_cmd = CAM_PAN;
    msgcmd_->rotate_pan_dx = dx;
    msgcmd_->rotate_pan_dy = dy;
  }

  static void reset() { msgcmd_->camera_cmd = CAM_RESET; }
  static void resetToCfg() { msgcmd_->camera_cmd = CAM_RESET_TO_CFG; }

  static MessageCommand* msgcmd_;
  static float rotate_pan_sensitivity_;
  static Camera* camera_;
  static const Config* cfg_;
  static float zoom_sensitivity_;
  static MouseState mouse_state_;
  static GLFWwindow* glfw_window_;
  static SceneT* scene_;
};

}  // namespace spray

#define SPRAY_GLFW_INL
#include "display/spray_glfw.inl"
#undef SPRAY_GLFW_INL

