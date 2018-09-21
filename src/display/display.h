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

namespace spray {

class Renderer;
class Model;
class Aabb;


class Display {
 public:
  Display(unsigned int image_width, unsigned int image_height,
          const Model& model, Renderer* renderer);
  void run();

 private:
  int initGlfw();
  void initGl();

  static void keyCallback(GLFWwindow* window, int key, int scancode, int action,
                          int mods);
  static void scrollCallback(GLFWwindow* window, double xoffset,
                             double yoffset);

  /**
   * This function returns the last state reported for the specified mouse
   * button to the specified window.
   *
   * \param button The mouse button that was pressed or released.
   * \param action One of GLFW_PRESS or GLFW_RELEASE.
   * \param mods Bit field describing which modifier keys were held down.
   */
  static void mouseButtonCallback(GLFWwindow* window, int button, int action,
                                  int mods);

  /**
   * This function sets the cursor position callback of the specified window,
   * which is called when the cursor is moved. The callback is provided with the
   * position, in screen coordinates, relative to the upper-left corner of the
   * client area of the window.
   *
   * \param xpos The new x-coordinate, in screen coordinates, of the cursor.
   * \param ypos The new y-coordinate, in screen coordinates, of the cursor.
   */
  static void cursorPositionCallback(GLFWwindow* window, double xpos,
                                     double ypos);

  void drawAabb(const Aabb& a);

 private:
  static Renderer* renderer_;
  const Model& model_;  ///< scene maintained outside of this display object
  GLFWwindow* window_;
  unsigned int image_w_;
  unsigned int image_h_;
};

}  // namespace spray

