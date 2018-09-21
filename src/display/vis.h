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

#include <stack>

#include "accelerators/wbvh_node.h"
#include "display/opengl.h"
#include "renderers/spray.h"

namespace spray {

class WbvhNode;

template <class WbvhT>
struct WbvhObj {
  WbvhObj() : ptr(nullptr) {}
  WbvhT* ptr;
};

template <class WbvhT>
class Vis {
 public:
  static void initialize(const WbvhObj<WbvhT>& wobj);

  static bool initializeTraversal();

  static bool isWbvhStackEmpty();
  static WbvhNode* wbvhStackTop();
  static bool isWbvhRoot(WbvhNode* node);
  static void wbvhStackPop();
  static void wbvhStackPush(WbvhNode* node);

  static void renderOpenGl();
  static void renderWaccel();
  static void renderDaccel();

 private:
  static WbvhObj<WbvhT> wbvh_;
  static std::stack<WbvhNode*> wbvh_stack_;
  static WbvhNode* wbvh_root_;
};

}  // namespace spray

#define SPRAY_VIS_INL_
#include "display/vis.inl"
#undef SPRAY_VIS_INL_

