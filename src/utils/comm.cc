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

#include "utils/comm.h"

#include <climits>

#include "pbrt/memory.h"
#include "glog/logging.h"

namespace spray {

msg_word_t* AllocMsg(std::size_t bytes) {
  std::size_t count = (bytes + MSG_WORD_SIZE - 1) / MSG_WORD_SIZE;
  CHECK_LT(count, INT_MAX);
  msg_word_t* msg = AllocAligned<msg_word_t>(count);
  CHECK_NOTNULL(msg);
  // *msg_word_count = (int)count;
  return msg;
}

msg_word_t* AllocMsg(std::size_t bytes, int* msg_word_count) {
  std::size_t count = (bytes + MSG_WORD_SIZE - 1) / MSG_WORD_SIZE;
  CHECK_LT(count, INT_MAX);
  msg_word_t* msg = AllocAligned<msg_word_t>(count);
  CHECK_NOTNULL(msg);
  *msg_word_count = (int)count;
  return msg;
}

int MsgWordCount(std::size_t bytes) {
  std::size_t count = (bytes + MSG_WORD_SIZE - 1) / MSG_WORD_SIZE;
  CHECK_LT(count, INT_MAX);
  return count;
}

}  // namespace spray
