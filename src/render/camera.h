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

#include <cmath>

#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"

#include "render/aabb.h"

//debug
#include "render/spray.h"
#include "render/shape_buffer.h"

namespace spray {

class Camera {
 public:
  void init(const glm::vec3 campos, const glm::vec3 lookat,
            const glm::vec3 upvec, float vfov, int image_w, int image_h);

  /** Resets the camera position and lookat. */
  void reset(const glm::vec3 campos, const glm::vec3 lookat,
             const glm::vec3 upvec);

  /**
   * Resets the image size.
   *
   * \param w Image width.
   * \param h Image height.
   */
  void resize(unsigned int w, unsigned int h);

  /** Returns the current camera position in world space. */
  const glm::vec3& getPosition() const { return camera_pos_; }

  /** Returns the current lookat position in world space. */
  const glm::vec3& getLookAt() const { return camera_lookat_; }

  /** Returns the current up vector in world space. */
  const glm::vec3& getUpVector() const { return camera_up_; }

  /**
   * Returns the up vector of the camera in world space.
   */
  const glm::vec3& getUp() const { return camera_up_; }

  /** Returns the vertical field of view in radians. */
  float getVfov() const { return camera_vfov_radians_; }

  /** Returns the aspect ratio of the image plane (width / height). */
  float getAspectRatio() const { return image_aspect_ratio_; }

  void generateRay(float x, float y, glm::vec3* org, glm::vec3* dir) const;
  void generateRay(float x, float y, float org[3], float dir[3]) const;
  void generateRay(float x, float y, float dir[3]) const;

  /**
   * Translates the camera along the z-axis upon receiving the mouse scrolling
   * event.
   *
   * \param yoffset The scroll offset along the z-axis in camera space.
   */
  void zoom(float offset);

  /**
   * Rotates the camera with respect to the scene center upon receiving the
   * mouse dragging event while the left mouse button being pushed down.
   * mouse_dx and mouse_dy are internally translated to the angles, theta and
   * phi, respectively, in the spherical coordinate system.
   *
   * \param mouse_dx displacement of mouse cursor in screen space x-drection.
   * \param mouse_dy displacement of mouse cursor in screen space y-drection.
   */
  void rotate(float mouse_dx, float mouse_dy);

  /**
   * Pans the camera upon receiving the mouse dragging event while the right
   * mouse button being pushed down.
   * mouse_dx and mouse_dy are internally translated to the angles, theta and
   * phi, respectively, in the spherical coordinate system.
   *
   * \param mouse_dx displacement of mouse cursor in screen space x-drection.
   * \param mouse_dy displacement of mouse cursor in screen space y-drection.
   */
  void pan(float mouse_dx, float mouse_dy);

  void print(int x, int y) {
    float dir[3];
    generateRay(x, image_h_ - y, dir);

    RTCRay r;
    r.org[0] = camera_pos_[0];
    r.org[1] = camera_pos_[1];
    r.org[2] = camera_pos_[2];
    r.dir[0] = dir[0];
    r.dir[1] = dir[1];
    r.dir[2] = dir[2];
    r.tnear = SPRAY_RAY_EPSILON;
    r.tfar = SPRAY_FLOAT_INF;
    r.instID = RTC_INVALID_GEOMETRY_ID;
    r.geomID = RTC_INVALID_GEOMETRY_ID;
    r.primID = RTC_INVALID_GEOMETRY_ID;
    r.mask = 0xFFFFFFFF;
    r.time = 0.0f;

    glm::vec3 center(-5, 1, 0);
    float radius = 2;
    Sphere sphere(center, radius, nullptr);
    sphere.setGeomId(10);

    ShapeBuffer::intersect(&sphere, r);
    if (r.geomID == 10) {
      std::cout << "(" << x << "," << image_h_ - y << "):" << dir[0] << ","
                << dir[1] << "," << dir[2] << " ,(hit)\n";
    } else {
      std::cout << "(" << x << "," << image_h_ - y << "):" << dir[0] << ","
                << dir[1] << "," << dir[2] << " ,(miss)\n";
    }
  }

 private:
  glm::vec3 getDirection(float u, float v) const;

 private:
  glm::vec3 camera_pos_;       ///< Camera position.
  glm::vec3 camera_lookat_;    ///< Camera lookat position.
  glm::vec3 camera_up_;        ///< Camera up vector.
  float camera_vfov_degrees_;  ///< Vertical field of view in degrees.
  float camera_vfov_radians_;  ///< Vertical field of view in radians.

  glm::vec3 u_vec_;
  glm::vec3 v_vec_;

  /** Position of lower-left corner of the image plane.*/
  glm::vec3 image_lowerleft_pos_;
  glm::vec3 image_w_vec_;  ///< x-axis in the image plane
  glm::vec3 image_h_vec_;  ///< y-axis in the image plane

  float image_w_;             ///< Image width.
  float image_h_;             ///< Image height.
  float image_aspect_ratio_;  ///< Image aspect ratio.
  float image_half_w_;
  float image_half_h_;
};

inline void Camera::init(const glm::vec3 campos, const glm::vec3 lookat,
                         const glm::vec3 upvec, float vfov, int image_w,
                         int image_h) {
  float aspect_ratio =
      static_cast<float>(image_w) / static_cast<float>(image_h);
  float theta = glm::radians(vfov);
  float image_half_h = glm::tan(theta * 0.4f);
  float image_half_w = aspect_ratio * image_half_h;

  glm::vec3 w = campos - lookat;
  glm::vec3 u = glm::cross(upvec, w);
  glm::vec3 v = glm::cross(w, u);

  w = glm::normalize(w);
  u = glm::normalize(u);
  v = glm::normalize(v);

  glm::vec3 image_center = campos - w;
  image_lowerleft_pos_ = image_center - (image_half_w * u) - (image_half_h * v);

  image_w_vec_ = 2.0f * image_half_w * u;
  image_h_vec_ = 2.0f * image_half_h * v;

  // init member variables
  camera_pos_ = campos;
  camera_lookat_ = lookat;
  camera_up_ = upvec;
  camera_vfov_degrees_ = vfov;
  camera_vfov_radians_ = theta;

  image_w_ = static_cast<float>(image_w);
  image_h_ = static_cast<float>(image_h);
  image_aspect_ratio_ = aspect_ratio;
  image_half_w_ = image_half_w;
  image_half_h_ = image_half_h;

  u_vec_ = u;
  v_vec_ = v;
}

inline void Camera::reset(const glm::vec3 campos, const glm::vec3 lookat,
                          const glm::vec3 upvec) {
  init(campos, lookat, upvec, camera_vfov_degrees_, image_w_, image_h_);
}

inline glm::vec3 Camera::getDirection(float u, float v) const {
  glm::vec3 dir = image_lowerleft_pos_ + (u * image_w_vec_) +
                  (v * image_h_vec_) - camera_pos_;
  return glm::normalize(dir);
}

inline void Camera::generateRay(float x, float y, glm::vec3* org,
                                glm::vec3* dir) const {
  *org = camera_pos_;
  float u = x / image_w_;
  float v = y / image_h_;
  *dir = getDirection(u, v);
}

inline void Camera::generateRay(float x, float y, float org[3],
                                float dir[3]) const {
  org[0] = camera_pos_[0];
  org[1] = camera_pos_[1];
  org[2] = camera_pos_[2];

  float u = x / image_w_;
  float v = y / image_h_;
  glm::vec3 direction = getDirection(u, v);

  dir[0] = direction[0];
  dir[1] = direction[1];
  dir[2] = direction[2];
}

inline void Camera::generateRay(float x, float y, float dir[3]) const {
  float u = x / image_w_;
  float v = y / image_h_;
  glm::vec3 direction = getDirection(u, v);
  dir[0] = direction[0];
  dir[1] = direction[1];
  dir[2] = direction[2];
}

inline void Camera::zoom(float offset) {
  glm::vec3 view_dir = camera_lookat_ - camera_pos_;

  glm::vec3 dir = offset * glm::normalize(view_dir);
  glm::vec3 new_campos = camera_pos_ + dir;
  glm::vec3 new_lookat = camera_lookat_ + dir;

  reset(new_campos, new_lookat, camera_up_);
}

inline void Camera::rotate(float mouse_dx, float mouse_dy) {
  float dtheta = mouse_dx * (M_PI / image_w_);
  float dphi = mouse_dy * (M_PI / image_h_);

  glm::vec3 p = camera_pos_ - camera_lookat_;

  glm::mat4 mv = glm::rotate(glm::mat4(1.0f), dtheta, v_vec_);
  glm::mat4 mu = glm::rotate(glm::mat4(1.0f), dphi, u_vec_);
  glm::vec4 p2 = mu * (mv * glm::vec4(p, 1.0f));

  glm::vec3 pos = glm::vec3(p2) + camera_lookat_;

  reset(pos, camera_lookat_, camera_up_);
}

inline void Camera::pan(float mouse_dx, float mouse_dy) {
  float dx = 0.001f * mouse_dx;
  float dy = -0.001f * mouse_dy;

  glm::vec3 offset = (dx * u_vec_) + (dy * v_vec_);
  glm::vec3 pos = camera_pos_ + offset;
  glm::vec3 lookat = camera_lookat_ + offset;

  reset(pos, lookat, camera_up_);
}

}  // namespace spray
