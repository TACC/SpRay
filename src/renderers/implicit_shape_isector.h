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

#include "renderers/rays.h"

namespace spray {

class ImplicitShapeIsector {
 public:
  /**
   * An intersector used for processing eye rays.
   *
   * \param ndomains Number of domains.
   * \param scene Target scene to test against.
   * \param ray_buf An allocated ray buffer.
   * \param qs A set of ray queues used to save ray pointers as a result of
   * intersection tests.
   */
  template <typename CacheT>
  void intersect(int ndomains, Scene<CacheT>* scene, RayBuf ray_buf,
                 spray::QVector<Ray*>* qs) {
    Ray* rays = ray_buf.rays;

    for (std::size_t i = 0; i < ray_buf.num; ++i) {
      Ray* ray = &rays[i];

      RTCRayUtil::makeRayForDomainIntersection(ray->org, ray->dir, &domains_,
                                               &eray_);

      scene->intersectDomains(eray_);

      if (domains_.count) {
        RTCRayUtil::sortDomains(domains_, hits_);

        for (int d = 0; d < domains_.count; ++d) {
          int id = hits_[d].id;
#ifdef SPRAY_GLOG_CHECK
          CHECK_LT(id, ndomains);
#endif
          qs->push(id, ray);
        }
      }
    }
  }

 private:
  DomainList domains_;  ///< A fixed-size list of hit domains.
};

}  // namespace spray
