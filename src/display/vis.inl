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

#if !defined(SPRAY_VIS_INL_)
#error An implementation of Vis
#endif

namespace spray {

template <class WbvhT>
WbvhObj<WbvhT> Vis<WbvhT>::wbvh_;

template <class WbvhT>
std::stack<WbvhNode*> Vis<WbvhT>::wbvh_stack_;

template <class WbvhT>
WbvhNode* Vis<WbvhT>::wbvh_root_ = nullptr;

template <class WbvhT>
void Vis<WbvhT>::initialize(const WbvhObj<WbvhT>& wobj) {
  wbvh_ = wobj;
}

template <class WbvhT>
bool Vis<WbvhT>::initializeTraversal() {
  enum { FAIL = 0, SUCCESS = 1 };

  if (!wbvh_.ptr) return (bool)FAIL;

  while (!wbvh_stack_.empty()) {
    wbvh_stack_.pop();
  }
  wbvh_root_ = wbvh_.ptr->getRoot();
  wbvh_stack_.push(wbvh_root_);

  return (bool)SUCCESS;
}

template <class WbvhT>
bool Vis<WbvhT>::isWbvhStackEmpty() {
  return wbvh_stack_.empty();
}

template <class WbvhT>
WbvhNode* Vis<WbvhT>::wbvhStackTop() {
  return wbvh_stack_.top();
}

template <class WbvhT>
bool Vis<WbvhT>::isWbvhRoot(WbvhNode* node) {
  return wbvh_root_ && (node == wbvh_root_);
}

template <class WbvhT>
void Vis<WbvhT>::wbvhStackPop() {
  wbvh_stack_.pop();
}

template <class WbvhT>
void Vis<WbvhT>::wbvhStackPush(WbvhNode* node) {
  wbvh_stack_.push(node);
}

template <class WbvhT>
void Vis<WbvhT>::renderOpenGl() {}

template <class WbvhT>
void Vis<WbvhT>::renderWaccel() {
  glPushAttrib(GL_ENABLE_BIT | GL_LINE_BIT);
  glDisable(GL_LIGHTING);
  glLineWidth(1);
  glEnable(GL_DEPTH_TEST);

  glm::vec4 color_allnodes(.2, .2, .2, 0);
  glm::vec4 color_parent(0, 1, 1, .5);
  glm::vec4 color_lchild(1, 0, 0, .5);
  glm::vec4 color_rchild(0, 0, 1, .5);

  WbvhNode* node = wbvh_stack_.top();

  glPolygonOffset(1.0, 1.0);
  glEnable(GL_POLYGON_OFFSET_FILL);

  WbvhNode* lchild = node->lchild;
  WbvhNode* rchild = node->rchild;

  glDisable(GL_POLYGON_OFFSET_FILL);

  glDepthMask(GL_FALSE);

  std::stack<WbvhNode*> tstack;
  tstack.push(wbvh_root_);

  while (!tstack.empty()) {
    WbvhNode* current = tstack.top();
    tstack.pop();

    current->aabb.draw(color_allnodes);
    if (current->lchild) tstack.push(current->lchild);
    if (current->rchild) tstack.push(current->rchild);
  }
  node->aabb.draw(color_parent);

  glLineWidth(2.f);
  if (lchild) lchild->aabb.draw(color_lchild);
  if (rchild) rchild->aabb.draw(color_rchild);

  glDepthMask(GL_TRUE);
  glPopAttrib();
}

template <class WbvhT>
void Vis<WbvhT>::renderDaccel() {}

}  // namespace spray

