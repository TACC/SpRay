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

#include "baseline/baseline_insitu_tracer.h"

namespace spray {
namespace baseline {

template <typename SceneT, typename ScheduleT, typename ShaderT>
class ImageTracer : public InsituTracer<SceneT, ScheduleT, ShaderT> {
 public:
  typedef InsituTracer<SceneT, ScheduleT, ShaderT> Base;
  void trace();
  int type() const { return TRACER_TYPE_BASELINE_IMAGE; }

 private:
  void schedule(int ndomains, const ArenaQs<DRayQItem> &qs,
                const ArenaQs<DRayQItem> &sqs, QStats *stats);
};

}  // namespace baseline
}  // namespace spray

#define SPRAY_BASELINE_IMAGE_TRACER_INL
#include "baseline/baseline_image_tracer.inl"
#undef SPRAY_BASELINE_IMAGE_TRACER_INL

