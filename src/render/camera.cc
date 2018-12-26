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

#include "render/camera.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>

#include "glm/glm.hpp"

#include "render/rays.h"

#define EPS_F 0.00001f

namespace spray {

Camera::Camera(const Aabb& scene_aabb, unsigned int image_w,
               unsigned int image_h, float znear, float zfar, float vfov) {
  init(scene_aabb, image_w, image_h, znear, zfar, vfov);
}

void Camera::init(const Aabb& scene_aabb, unsigned int image_w,
                  unsigned int image_h, float znear, float zfar, float vfov) {
  scene_aabb_ = scene_aabb;

  // image size
  image_w_ = image_w;
  image_h_ = image_h;

  // clipping planes
  znear_ = znear;
  zfar_ = zfar;

  ar_ = static_cast<float>(image_w) / image_h;

  vfov_ = glm::radians(vfov);
  hfov_ = ar_ * vfov_;

  tan_half_vfov_ = tanf(vfov_ * 0.5);
  tan_half_hfov_ = tanf(hfov_ * 0.5);

#ifdef DEBUG_CAMERA
  printf("hfov %f vfov %f ar %f tan(hfov/2) %f tan(v_fov/2) %f\n", hfov_, vfov_,
         ar_, tan_half_hfov_, tan_half_vfov_);
#endif

  resetPosition();
}

void Camera::reset() { resetPosition(); }

void Camera::resize(unsigned int w, unsigned int h) {
  image_w_ = w;
  image_h_ = h;
  ar_ = static_cast<float>(w) / h;

  hfov_ = ar_ * vfov_;

  tan_half_hfov_ = tanf(hfov_ * 0.5);
}

inline void Camera::resetPosition() {
  phi_ = glm::radians(90.0);
  theta_ = 0.0;
  center_ = scene_aabb_.getCenter();

  float extent = glm::length(scene_aabb_.getExtent());
  radius_ = extent * 1.5;
  min_radius_ = extent * 0.0001;
  max_radius_ = extent * 10.0;

  pos_ = center_ + (radius_ * up_);

  cam2world_[0] = glm::vec3(1.f, 0.f, 0.f);
  cam2world_[1] = glm::vec3(0.f, 1.f, 0.f);
  cam2world_[2] = glm::vec3(0.f, 0.f, 1.f);
}

void Camera::resetPosition(const glm::vec3& pos, const glm::vec3& center,
                           const glm::vec3& up) {
  center_ = center;

  pos_ = pos;
  up_ = up;

  glm::vec3 dir2cam = pos - center;
  radius_ = glm::length(dir2cam);

  theta_ = glm::atan(dir2cam.y, dir2cam.z);  // y / z

  dir2cam = glm::normalize(dir2cam);

  phi_ = glm::acos(dir2cam.y);

  float sin_phi = sin(phi_);
  if (sin_phi == 0) {
    phi_ += EPS_F;
    sin_phi = sin(phi_);
  }

  // glm::vec3 up(0, sin_phi > 0 ? 1 : -1, 0);
  glm::vec3 up2 = sin_phi > 0 ? up : -up;

  glm::vec3 film_xdir = glm::normalize(cross(up2, dir2cam));
  glm::vec3 film_ydir = glm::normalize(cross(dir2cam, film_xdir));

  cam2world_[0] = film_xdir;
  cam2world_[1] = film_ydir;
  cam2world_[2] = dir2cam;  // normalized
}

void Camera::updatePosition() {
  float sin_phi = sin(phi_);
  if (sin_phi == 0) {
    phi_ += EPS_F;
    sin_phi = sin(phi_);
  }
  glm::vec3 dir2cam(radius_ * sin_phi * sin(theta_), radius_ * cos(phi_),
                    radius_ * sin_phi * cos(theta_));

  pos_ = center_ + dir2cam;

  // glm::vec3 up(0, sin_phi > 0 ? 1 : -1, 0);
  // glm::vec3 up = up_;
  glm::vec3 up = sin_phi > 0 ? up_ : -up_;

#ifdef DEBUG_CAMERA
  printf("up (%f %f %f)\n", up[0], up[1], up[2]);
#endif

  glm::vec3 film_xdir = glm::normalize(cross(up, dir2cam));
  glm::vec3 film_ydir = glm::normalize(cross(dir2cam, film_xdir));

  cam2world_[0] = film_xdir;
  cam2world_[1] = film_ydir;
  cam2world_[2] = glm::normalize(dir2cam);
}

void Camera::generateRay(float x, float y, glm::vec3* org,
                         glm::vec3* dir) const {
  // raster space to ndc space
  // raster space
  //  (0,0)     (W, 0)
  //  (0,H)     (W, H)
  // ndc space
  //  (0,0)     (1, 0)
  //  (0,1)     (1, 1)
  // // add 0.5 to account for the center of a pixel
  // float x_ndc = (x + 0.5) / image_w_;
  // float y_ndc = (y + 0.5) / image_h_;
  float x_ndc = x / image_w_;
  float y_ndc = y / image_h_;

  // ndc space to screen space
  // ndc space
  //  (0,0)     (1, 0)
  //  (0,1)     (1, 1)
  // screen space
  //  (-1, 1)   (1, 1)
  //  (-1,-1)   (1,-1)
  // float x_screen = 2.f * x_ndc - 1.f;
  // // float y_screen = 1.f - (2.f * y_ndc);
  // float y_screen = 2.f * y_ndc - 1.f;
  // // float x_screen = x_ndc - .5f;
  // // float y_screen = y_ndc - .5f;
  float x_screen = 2 * x_ndc - 1;
  float y_screen = 2 * y_ndc - 1;

  // space space to camera space
  // assume the image plane is z=-1 apart from the camera position
  // multiply aspect_ratio to make pixel shapes squares
  float x_camera = x_screen * ar_ * tan_half_vfov_;  // x focal_len = 1 ommitted
  float y_camera = y_screen * tan_half_vfov_;

  glm::vec3 pix_pos(x_camera, y_camera, -1.0);  // camera space

  *org = pos_;
  *dir = glm::normalize(cam2world_ * pix_pos);
}

void Camera::generateRay(float x, float y, float dir[3]) const {
  // raster space to ndc space
  float x_ndc = x / image_w_;
  float y_ndc = y / image_h_;

  // ndc space to screen space
  float x_screen = 2 * x_ndc - 1;
  float y_screen = 2 * y_ndc - 1;

  // space space to camera space
  // assume the image plane is z=-1 apart from the camera position
  // multiply aspect_ratio to make pixel shapes squares
  float x_camera = x_screen * ar_ * tan_half_vfov_;  // x focal_len = 1 ommitted
  float y_camera = y_screen * tan_half_vfov_;

  glm::vec3 pix_pos(x_camera, y_camera, -1.0);  // camera space

  // org[0] = pos_.x;
  // org[1] = pos_.y;
  // org[2] = pos_.z;

  glm::vec3 direction = glm::normalize(cam2world_ * pix_pos);

  dir[0] = direction.x;
  dir[1] = direction.y;
  dir[2] = direction.z;
}

void Camera::zoom(float offset) {
  float new_radius = glm::clamp(radius_ - offset, min_radius_, max_radius_);
  pos_ = (pos_ - center_) * (new_radius / radius_) + center_;
  radius_ = new_radius;
}

void Camera::rotate(float mouse_dx, float mouse_dy) {
  // assume image width/height is 2*pi
  float dtheta = mouse_dx * (M_PI / image_w_);
  float dphi = mouse_dy * (M_PI / image_h_);

  theta_ += dtheta;
  phi_ = glm::clamp(phi_ + dphi, (float)0.0, (float)M_PI);

  updatePosition();
}

void Camera::pan(float mouse_dx, float mouse_dy) {
  float scale = 0.001;
  float dx = -mouse_dx * scale;
  float dy = mouse_dy * scale;
  glm::vec3 disp = (cam2world_[0] * dx) + (cam2world_[1] * dy);
  pos_ += disp;
  center_ += disp;
}

}  // namespace spray

