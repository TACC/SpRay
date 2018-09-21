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

#if !defined(SPRAY_GLFW_INL)
#error An implementation of Glfw
#endif

namespace spray {

template <class WbvhT, class CacheT>
MessageCommand* Glfw<WbvhT, CacheT>::msgcmd_ = nullptr;

template <class WbvhT, class CacheT>
float Glfw<WbvhT, CacheT>::rotate_pan_sensitivity_ = 8.0f;

template <class WbvhT, class CacheT>
const Config* Glfw<WbvhT, CacheT>::cfg_ = nullptr;

template <class WbvhT, class CacheT>
Camera* Glfw<WbvhT, CacheT>::camera_ = nullptr;

template <class WbvhT, class CacheT>
float Glfw<WbvhT, CacheT>::zoom_sensitivity_ = 2.0f;

template <class WbvhT, class CacheT>
MouseState Glfw<WbvhT, CacheT>::mouse_state_;

template <class WbvhT, class CacheT>
GLFWwindow* Glfw<WbvhT, CacheT>::glfw_window_ = nullptr;

template <class WbvhT, class CacheT>
Scene<CacheT>* Glfw<WbvhT, CacheT>::scene_ = nullptr;

template <class WbvhT, class CacheT>
void Glfw<WbvhT, CacheT>::initialize(const Config& cfg, bool is_root_process,
                                     unsigned image_w, unsigned image_h,
                                     Camera* camera, MessageCommand* cmd,
                                     Scene<CacheT>* scene) {
  cfg_ = &cfg;

  camera_ = camera;
  scene_ = scene;

  msgcmd_ = cmd;

  mouse_state_.left = false;
  mouse_state_.right = false;
  mouse_state_.middle = false;

  if (!is_root_process) return;

  CHECK(glfwInit()) << "error failed to initialize glfw";

  glfw_window_ = glfwCreateWindow(image_w, image_h, "spray", NULL, NULL);
  if (!glfw_window_) {
    glfwTerminate();
    LOG(FATAL) << "error failed to create glfw window";
  }

  // int buffer_w, buffer_h;
  // glfwGetFramebufferSize(window_, &buffer_w, &buffer_h);

  glfwMakeContextCurrent(glfw_window_);

  glfwSwapInterval(1);

  glfwSetKeyCallback(glfw_window_, keyCallback);
  glfwSetScrollCallback(glfw_window_, scrollCallback);
  glfwSetMouseButtonCallback(glfw_window_, mouseButtonCallback);
  glfwSetCursorPosCallback(glfw_window_, cursorPosCallback);

  glClearColor(0.0, 0.0, 0.0, 0.0);  // black background
  glShadeModel(GL_FLAT);

  glViewport(0, 0, image_w, image_h);

  glm::mat4 projection =
      glm::perspective(camera_->getVfov(), camera_->getAspectRatio(),
                       camera_->getZnear(), camera_->getZfar());
  glMatrixMode(GL_PROJECTION);
}

template <class WbvhT, class CacheT>
void Glfw<WbvhT, CacheT>::cmdHandler() {
  if (msgcmd_->camera_cmd == CAM_ZOOM) {
    camera_->zoom(msgcmd_->zoom_offset);
  } else if (msgcmd_->camera_cmd == CAM_ROTATE) {
    camera_->rotate(msgcmd_->rotate_pan_dx, msgcmd_->rotate_pan_dy);
  } else if (msgcmd_->camera_cmd == CAM_PAN) {
    camera_->pan(msgcmd_->rotate_pan_dx, msgcmd_->rotate_pan_dy);
  } else if (msgcmd_->camera_cmd == CAM_RESET) {
    camera_->reset();
  } else if (msgcmd_->camera_cmd == CAM_RESET_TO_CFG) {
    camera_->resetPosition(cfg_->camera_pos, cfg_->camera_lookat,
                           cfg_->camera_up);
  }
  msgcmd_->camera_cmd = CAM_NOP;
}

template <class WbvhT, class CacheT>
void Glfw<WbvhT, CacheT>::keyCallback(GLFWwindow* window, int key, int scancode,
                                      int action, int mods) {
  if (action == GLFW_PRESS) {
    switch (key) {
      case GLFW_KEY_ESCAPE:
      case GLFW_KEY_Q:
        msgcmd_->view_mode = VIEW_MODE_TERMINATE;
        msgcmd_->done = 1;
        break;
      case GLFW_KEY_C:
        camera_->reset();
        reset();
        break;
      case GLFW_KEY_SPACE:  // reset camera position according to initial
                            // configuration
        camera_->resetPosition(cfg_->camera_pos, cfg_->camera_lookat,
                               cfg_->camera_up);
        resetToCfg();
        break;
      case GLFW_KEY_EQUAL:  // to ray trace mode
        zoom_sensitivity_ *= 2.f;
        printf("zoom sensitivity: %f\n", zoom_sensitivity_);
        break;

      case GLFW_KEY_MINUS:  // to ray trace mode
        zoom_sensitivity_ /= 2.f;
        printf("zoom sensitivity: %f\n", zoom_sensitivity_);
        break;

      case GLFW_KEY_0:  // to ray trace mode
        rotate_pan_sensitivity_ *= 2.f;
        printf("rotate/pan sensitivity: %f\n", rotate_pan_sensitivity_);
        break;

      case GLFW_KEY_9:  // to ray trace mode
        rotate_pan_sensitivity_ /= 2.f;
        printf("rotate/pan sensitivity: %f\n", rotate_pan_sensitivity_);
        break;

      case GLFW_KEY_G:  // glfw
        msgcmd_->view_mode = VIEW_MODE_GLFW;
        break;

      // case GLFW_KEY_O:  // opengl
      //   msgcmd_->view_mode = VIEW_MODE_OPENGL;
      //   break;

      case GLFW_KEY_W: {  // wbvh
        bool success = Vis<WbvhT>::initializeTraversal();
        std::cout << "wbvh initialization done " << success << "\n";
        if (success) {
          std::cout << "Switching view mode to WBVH" << std::endl;
          msgcmd_->view_mode = VIEW_MODE_WBVH;
        } else {
          msgcmd_->view_mode = VIEW_MODE_TERMINATE;
          LOG(FATAL) << "Failed to switching to WBVH view mode.";
        }
      } break;

      case GLFW_KEY_E: {  // wbvh with overlay
        bool success = Vis<WbvhT>::initializeTraversal();
        if (success) {
          std::cout << "Switching view mode to WBVH_OVERLAY" << std::endl;
          msgcmd_->view_mode = VIEW_MODE_WBVH_OVERLAY;
        } else {
          msgcmd_->view_mode = VIEW_MODE_TERMINATE;
          LOG(FATAL) << "Failed to switching to WBVH_OVERLAY view mode.";
        }
      } break;

      case GLFW_KEY_D:  // dbvh
        msgcmd_->view_mode = VIEW_MODE_BVH;
        break;

      case GLFW_KEY_P: {
        const glm::vec3& p = camera_->getPosition();
        const glm::vec3& c = camera_->getCenter();
        const glm::vec3& u = camera_->getUp();
        const Aabb& bound = camera_->getSceneBound();
        printf("[INFO] Camera (pos,center) --camera %f %f %f %f %f %f\n", p.x,
               p.y, p.z, c.x, c.y, c.z);
        printf("[INFO] Camera (upvector) --up %f %f %f\n", u.x, u.y, u.z);
        printf("[INFO] Camera (upvector cfg) --up %f %f %f\n",
               cfg_->camera_up.x, cfg_->camera_up.y, cfg_->camera_up.z);
        printf("[INFO] Scene bound: min(%f %f %f), max(%f %f %f)\n",
               bound.bounds[0].x, bound.bounds[0].y, bound.bounds[0].z,
               bound.bounds[1].x, bound.bounds[1].y, bound.bounds[1].z);
      } break;

      default:
        break;
    }  // switch (key) {

    int view_mode = msgcmd_->view_mode;

    // DOMAIN mode
    if (view_mode == VIEW_MODE_DOMAIN) {
      int id;
      switch (key) {
        case GLFW_KEY_UP: {
          id = scene_->nextDomain(10);
        } break;

        case GLFW_KEY_DOWN: {
          id = scene_->prevDomain(10);
        } break;

        case GLFW_KEY_LEFT: {
          id = scene_->prevDomain(1);
        } break;

        case GLFW_KEY_RIGHT: {
          id = scene_->nextDomain(1);
        } break;

        default:
          break;
      }
      std::cout << "domain ID : " << id << std::endl;
    } else if (view_mode == VIEW_MODE_PARTITION) {
      int partition_num;
      switch (key) {
        case GLFW_KEY_LEFT: {
          partition_num = scene_->prevPartition();
          std::cout << "partition number : " << partition_num << std::endl;
        } break;

        case GLFW_KEY_RIGHT: {
          partition_num = scene_->nextPartition();
          std::cout << "partition number : " << partition_num << std::endl;
        } break;

        case GLFW_KEY_DOWN: {
          scene_->setPartitionNumber(-1);
          std::cout << "partition number : -1" << std::endl;
        } break;

        case GLFW_KEY_F: {
          scene_->printPartitionSizes();
        } break;

        default:
          break;
      }
    }

    WbvhNode* node = nullptr;
    if (!Vis<WbvhT>::isWbvhStackEmpty()) node = Vis<WbvhT>::wbvhStackTop();

    if (node &&
        (view_mode == VIEW_MODE_WBVH || view_mode == VIEW_MODE_WBVH_OVERLAY)) {
      switch (key) {
        case GLFW_KEY_UP: {
          if (!Vis<WbvhT>::isWbvhRoot(node)) Vis<WbvhT>::wbvhStackPop();
        } break;

        case GLFW_KEY_LEFT: {
          if (node->lchild) Vis<WbvhT>::wbvhStackPush(node->lchild);
        } break;

        case GLFW_KEY_RIGHT: {
          if (node->rchild) Vis<WbvhT>::wbvhStackPush(node->rchild);
        } break;

        default:
          break;
      }
    }
  }  // if (action == GLFW_PRESS) {
}

template <class WbvhT, class CacheT>
void Glfw<WbvhT, CacheT>::scrollCallback(GLFWwindow* window, double xoffset,
                                         double yoffset) {
  float offset = zoom_sensitivity_ * yoffset;
  zoom(offset);
}

template <class WbvhT, class CacheT>
void Glfw<WbvhT, CacheT>::mouseButtonCallback(GLFWwindow* window, int button,
                                              int action, int mods) {
  if (action == GLFW_PRESS) {
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
      mouse_state_.left = true;
    } else if (button == GLFW_MOUSE_BUTTON_RIGHT) {
      mouse_state_.right = true;
    }
  } else if (action == GLFW_RELEASE) {
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
      mouse_state_.left = false;
    } else if (button == GLFW_MOUSE_BUTTON_RIGHT) {
      mouse_state_.right = false;
    }
  }
}

template <class WbvhT, class CacheT>
void Glfw<WbvhT, CacheT>::cursorPosCallback(GLFWwindow* window, double xpos,
                                            double ypos) {
  if (mouse_state_.left && !mouse_state_.middle && !mouse_state_.right) {
    // left click
    float dx = rotate_pan_sensitivity_ * (xpos - mouse_state_.x);
    float dy = rotate_pan_sensitivity_ * (ypos - mouse_state_.y);
    rotate(dx, dy);
  } else if (!mouse_state_.left && !mouse_state_.middle && mouse_state_.right) {
    // right click
    float dx = rotate_pan_sensitivity_ * (xpos - mouse_state_.x);
    float dy = rotate_pan_sensitivity_ * (ypos - mouse_state_.y);
    pan(dx, dy);
  } else if (!mouse_state_.left && !mouse_state_.middle &&
             !mouse_state_.right) {
    // middle click
  }
  mouse_state_.x = xpos;
  mouse_state_.y = ypos;
}

}  // namespace spray

