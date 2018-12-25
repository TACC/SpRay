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

#include "display/opengl.h"
#include "render/aabb.h"

namespace spray {

class WbvhNode {
 public:
  WbvhNode() : lchild(nullptr), rchild(nullptr) {}

  bool isLeaf() const { return lchild == nullptr && rchild == nullptr; }

  void makeInternal(const Aabb &aabb, unsigned begin, unsigned end,
                    WbvhNode *lchild, WbvhNode *rchild) {
    this->aabb = aabb;
    this->begin = begin;
    this->end = end;
    this->lchild = lchild;
    this->rchild = rchild;
  }

  void makeLeaf(const Aabb &aabb, unsigned begin, unsigned end) {
    this->aabb = aabb;
    this->begin = begin;
    this->end = end;
    this->lchild = nullptr;
    this->rchild = nullptr;
  }

  void draw() const {
    glDisable(GL_LIGHTING);
    aabb.draw(glm::vec4(0.0, 1.0, 0.0, 1.0));
    glEnable(GL_LIGHTING);
  }

  WbvhNode *lchild;
  WbvhNode *rchild;
  Aabb aabb;
  unsigned begin;
  unsigned end;
};

}  // namespace spray

