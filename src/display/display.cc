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

#include "display/display.h"

#include <cstdio>

#include "render/renderer.h"
#include "scene/aabb.h"

namespace spray {

Renderer* Display::renderer_ = NULL;

Display::Display(unsigned int image_width, unsigned int image_height,
                 const Model& model, Renderer* renderer)
    : image_w_(image_width), image_h_(image_height), model_(model) {
  Display::renderer_ = renderer;

  if (initGlfw()) {
    printf("error failed to initialize glfw\n");
  }
  printf("done initializing glfw\n");

  initGl();
}

inline int Display::initGlfw() {
  if (!glfwInit()) return -1;

  window_ = glfwCreateWindow(image_w_, image_h_, "Hello World", NULL, NULL);
  if (!window_) {
    glfwTerminate();
    return -1;
  }

  // int buffer_w, buffer_h;
  // glfwGetFramebufferSize(window_, &buffer_w, &buffer_h);

  glfwMakeContextCurrent(window_);
  glfwSwapInterval(1);

  glfwSetKeyCallback(window_, keyCallback);
  glfwSetScrollCallback(window_, scrollCallback);
  glfwSetMouseButtonCallback(window_, mouseButtonCallback);
  glfwSetCursorPosCallback(window_, cursorPositionCallback);

  return 0;
}

inline void Display::initGl() {
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_LIGHTING);
}

void Display::drawAabb(const Aabb& a) {
  glPushMatrix();
  glBegin(GL_LINE_LOOP);
  glVertex3f(a.bounds[0][0], a.bounds[0][1], a.bounds[0][2]);
  glVertex3f(a.bounds[1][0], a.bounds[0][1], a.bounds[0][2]);
  glVertex3f(a.bounds[1][0], a.bounds[1][1], a.bounds[0][2]);
  glVertex3f(a.bounds[0][0], a.bounds[1][1], a.bounds[0][2]);
  glEnd();
  glBegin(GL_LINE_LOOP);
  glVertex3f(a.bounds[0][0], a.bounds[0][1], a.bounds[1][2]);
  glVertex3f(a.bounds[1][0], a.bounds[0][1], a.bounds[1][2]);
  glVertex3f(a.bounds[1][0], a.bounds[1][1], a.bounds[1][2]);
  glVertex3f(a.bounds[0][0], a.bounds[1][1], a.bounds[1][2]);
  glEnd();
  glBegin(GL_LINES);
  glVertex3f(a.bounds[0][0], a.bounds[0][1], a.bounds[0][2]);
  glVertex3f(a.bounds[0][0], a.bounds[0][1], a.bounds[1][2]);
  glVertex3f(a.bounds[1][0], a.bounds[0][1], a.bounds[0][2]);
  glVertex3f(a.bounds[1][0], a.bounds[0][1], a.bounds[1][2]);
  glVertex3f(a.bounds[0][0], a.bounds[1][1], a.bounds[0][2]);
  glVertex3f(a.bounds[0][0], a.bounds[1][1], a.bounds[1][2]);
  glVertex3f(a.bounds[1][0], a.bounds[1][1], a.bounds[0][2]);
  glVertex3f(a.bounds[1][0], a.bounds[1][1], a.bounds[1][2]);
  glEnd();
  glPopMatrix();
}

void Display::run() {
  while (!glfwWindowShouldClose(window_)) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    renderer_->render();
    glfwSwapBuffers(window_);
    glfwPollEvents();
  }
  glfwTerminate();
}

void Display::keyCallback(GLFWwindow* window, int key, int scancode, int action,
                          int mods) {
  if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
    glfwSetWindowShouldClose(window, GL_TRUE);
  } else if (key == GLFW_KEY_R && action == GLFW_PRESS) {
    renderer_->resetCameraPosition();
  }
}

void Display::scrollCallback(GLFWwindow* window, double xoffset,
                             double yoffset) {
  renderer_->scrollCallback(xoffset, yoffset);
}

void Display::mouseButtonCallback(GLFWwindow* window, int button, int action,
                                  int mods) {
  renderer_->mouseButtonCallback(button, action, mods);
}

void Display::cursorPositionCallback(GLFWwindow* window, double xpos,
                                     double ypos) {
  renderer_->cursorPositionCallback(xpos, ypos);
}

}  // namespace spray
