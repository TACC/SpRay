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

#include <omp.h>

#include "embree/random_sampler.h"
#include "glm/glm.hpp"
#include "glog/logging.h"

#include "insitu/insitu_ray.h"
#include "render/camera.h"
#include "render/spray.h"
#include "render/tile.h"
#include "utils/math.h"

namespace spray {
namespace insitu {

struct SPRAY_ALIGN(16) Ray {
  float org[3];  ///< Ray origin.
  int pixid;     ///< Pixel ID.
  float dir[3];  ///< Ray direction.
  int samid;     ///< Sample ID of the image plane
  float w[3];    ///< Color weight.
  float t;       ///< Distance to the hit point.
  int light;     ///< Light ID.
  /**
   * -2: background
   * -1: possibly background
   *  0: unoccluded
   *  1: occluded
   */
  int occluded;
#if 0
  friend std::ostream& operator<<(std::ostream& os, const Ray& r) {
    os << "[p " << r.pixid << "][s " << r.samid << "] [l " << r.light << "] [t "
       << r.t << "][occl " << r.occluded << "] [org " << r.org[0] << ","
       << r.org[1] << "," << r.org[2] << "]"
       << "[dir " << r.dir[0] << "," << r.dir[1] << "," << r.dir[2] << "]"
       << "[w " << r.w[0] << "," << r.w[1] << "," << r.w[2] << "]";
    return os;
  }
#endif
};

struct RayData {
  Ray* ray;
  float tdom;
  int dom_depth;
};

struct RayUtil {
  enum ValuesForOccluded {
    OFLAG_UNDEFINED = -3,
    OFLAG_BACKGROUND = -2,
    OFLAG_UNOCCLUDED = 0,
    OFLAG_OCCLUDED = 1
  };

  inline static void setOccluded(int flag, Ray* ray) { ray->occluded = flag; }

  inline static void makeShadow(const Ray& ray, int light, const glm::vec3& pos,
                                const glm::vec3& dir, const glm::vec3& w,
                                float t, Ray* shadow) {
    shadow->org[0] = pos[0];
    shadow->org[1] = pos[1];
    shadow->org[2] = pos[2];

    shadow->pixid = ray.pixid;

    shadow->dir[0] = dir[0];
    shadow->dir[1] = dir[1];
    shadow->dir[2] = dir[2];

    shadow->samid = ray.samid;

    shadow->w[0] = w[0];
    shadow->w[1] = w[1];
    shadow->w[2] = w[2];

    shadow->t = t;

    shadow->light = light;
    shadow->occluded = 0;
  }
  inline static void makeRay(const Ray& rayin, const glm::vec3& pos,
                             const glm::vec3& dir, const glm::vec3& w, float t,
                             Ray* rayout) {
    rayout->org[0] = pos[0];
    rayout->org[1] = pos[1];
    rayout->org[2] = pos[2];

    rayout->pixid = rayin.pixid;

    rayout->dir[0] = dir[0];
    rayout->dir[1] = dir[1];
    rayout->dir[2] = dir[2];

    rayout->samid = rayin.samid;

    rayout->w[0] = w[0];
    rayout->w[1] = w[1];
    rayout->w[2] = w[2];

    rayout->t = t;
  }
};

inline void genSingleSampleEyeRays(const Camera& camera, int image_w,
                                   float orgx, float orgy, float orgz,
                                   Tile blocking_tile, Tile tile,
                                   RayBuf<Ray>* ray_buf) {
  Ray* rays = ray_buf->rays;

#pragma omp for collapse(2) schedule(static, 1)
  for (int y = tile.y; y < tile.y + tile.h; ++y) {
    for (int x = tile.x; x < tile.x + tile.w; ++x) {
      int bufid = tile.w * (y - tile.y) + (x - tile.x);
#ifdef SPRAY_GLOG_CHECK
      CHECK_LT(bufid, ray_buf->num);
#endif
      auto* ray = &rays[bufid];
      //
      ray->org[0] = orgx;
      ray->org[1] = orgy;
      ray->org[2] = orgz;

      ray->pixid = image_w * y + x;

      camera.generateRay((float)x, (float)y, ray->dir);

      ray->samid =
          blocking_tile.w * (y - blocking_tile.y) + (x - blocking_tile.x);

      ray->w[0] = 1.f;
      ray->w[1] = 1.f;
      ray->w[2] = 1.f;

      ray->t = SPRAY_FLOAT_INF;
    }
  }
}

inline void genMultiSampleEyeRays(const Camera& camera, int image_w, float orgx,
                                  float orgy, float orgz, int num_pixel_samples,
                                  spray::Tile blocking_tile, spray::Tile tile,
                                  RayBuf<Ray>* ray_buf) {
  Ray* rays = ray_buf->rays;

#pragma omp for collapse(3) schedule(static, 1)
  for (int y = tile.y; y < tile.y + tile.h; ++y) {
    for (int x = tile.x; x < tile.x + tile.w; ++x) {
      for (int s = 0; s < num_pixel_samples; ++s) {
        int bufid =
            (tile.w * (y - tile.y) + (x - tile.x)) * num_pixel_samples + s;
#ifdef SPRAY_GLOG_CHECK
        CHECK_LT(bufid, ray_buf->num);
#endif
        Ray* ray = &rays[bufid];
        //
        ray->org[0] = orgx;
        ray->org[1] = orgy;
        ray->org[2] = orgz;

        ray->pixid = image_w * y + x;

        RandomSampler sampler;
        RandomSampler_init(sampler, ray->pixid, s);

        float fx = (float)(x) + RandomSampler_get1D(sampler);
        float fy = (float)(y) + RandomSampler_get1D(sampler);

        camera.generateRay(fx, fy, ray->dir);

        ray->samid =
            (blocking_tile.w * (y - blocking_tile.y) + (x - blocking_tile.x)) *
                num_pixel_samples +
            s;

        ray->w[0] = 1.f;
        ray->w[1] = 1.f;
        ray->w[2] = 1.f;

        ray->t = SPRAY_FLOAT_INF;
      }
    }
  }
}

}  // namespace insitu
}  // namespace spray
