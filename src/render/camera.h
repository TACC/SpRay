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

#include "render/aabb.h"
#include "render/spray.h"

namespace spray {

class Camera {
 public:
  Camera() {}
  /**
   * Creates a simple pinhole camera following the OpenGL coordinate system.
   *
   * \param scene_aabb Scene bounds in world space.
   * \param image_w Image width in the number of pixels.
   * \param image_h Image height in the number of pixels.
   * \param znear Z-coordinate of the near clipping plane in camera space.
   * \param zfar Z-coordinate of the far clipping plane in camera space.
   * \param vfov Vertical field of view in degrees.
   */
  Camera(const Aabb& scene_aabb, unsigned int image_w, unsigned int image_h,
         float znear, float zfar, float vfov);

  void initialize(const Aabb& scene_aabb, unsigned int image_w,
                  unsigned int image_h, float znear, float zfar, float vfov);

  void resetPosition(const glm::vec3& pos, const glm::vec3& center,
                     const glm::vec3& up);

  /**
   * Resets the image size.
   *
   * \param w Image width.
   * \param h Image height.
   */
  void resize(unsigned int w, unsigned int h);

  /** Returns the current position of the camera in world space. */
  const glm::vec3& getPosition() const { return pos_; }

  /** Returns the coordinates of the scene center in world space. */
  const glm::vec3& getCenter() const { return center_; }

  // /** Returns the current up vector of the camera in world space. */
  // const glm::vec3& getUp() const { return up_; }

  /**
   * Returns the up vector of the camera in world space.
   * i.e. the second column of the camera-to-world transform matrix.
   */
  const glm::vec3& getUp() const { return cam2world_[1]; }

  const Aabb& getSceneBound() const { return scene_aabb_; }

  /** Returns the vertical field of view in radians. */
  float getVfov() const { return vfov_; }

  /** Returns the horizontal field of view in radians. */
  float getHfov() const { return hfov_; }

  /** Returns the aspect ratio of the image plane (width / height). */
  float getAspectRatio() const { return ar_; }

  /** Returns the z-coordinate of the near clipping plane in camera space. */
  float getZnear() const { return znear_; }

  /** Returns the z-coordinate of the far clipping plane in camera space. */
  float getZfar() const { return zfar_; }

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

  /**
   * Resets the camera to a default position.
   */
  void reset();

 private:
  /**
   * Updates the position of the camera (#pos_) and the camera-to-world
   * transform matrix (#cam2world_). This function must be called whenever
   * #phi_ or #theta, is updated.
   */
  void updatePosition();

  /**
   * Resets the camera to a default position.
   */
  void resetPosition();

 private:
  Aabb scene_aabb_;  ///< Scene bounds.

  unsigned int image_w_;  ///< Image width.
  unsigned int image_h_;  ///< Image height.
  float znear_;           ///< Near clip plane.
  float zfar_;            ///< Far clip plane.
  float zscreen_;         ///< Distance between camera and screen
  float vfov_;            ///< Vertical field of view in radians.
  float hfov_;            ///< Horizontal field of view in radians.
  float tan_half_vfov_;   ///< Tangent of v_fov_ / 2.
  float tan_half_hfov_;   ///< Tangent of h_fov_ / 2.

  glm::vec3 pos_;     ///< Camera position in world space
  glm::vec3 center_;  ///< Center of scene in world space (viewing target)
  glm::vec3 up_;

  /**
   * 3x3 Camera-to-world space transform matrix. Note that to achieve the full
   * transformation, the current position (#pos_) must be added.
   */
  glm::mat3 cam2world_;

  float ar_;  ///< aspect ratio (image_width / image_height).

  float phi_;         ///< Pitch in radians.
  float theta_;       ///< Yaw in radians.
  float radius_;      ///< Distance between center of scene and camera position.
  float min_radius_;  ///< Minimum radius.
  float max_radius_;  ///< Maximum radius.
};

}  // namespace spray

