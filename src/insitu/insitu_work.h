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

#include "pbrt/memory.h"

#include "insitu/insitu_comm.h"
#include "insitu/insitu_ray.h"
#include "insitu/insitu_tiler.h"
#include "renderers/spray.h"

namespace spray {
namespace insitu {

enum WorkType {
  WORK_SEND_SHADS,
  WORK_SEND_RADS,
  MSG_TERMINATE,
};

class VBuf;

struct MsgHeader {
  int domain_id;
  int64_t payload_count;
};

class Work {
 public:
  Work(int t) : type(t), msg(nullptr), count(-1), dest(-1) {}

  // for sending (msg is allocated by the derived)
  Work(int t, int dest) : type(t), msg(nullptr), count(-1), dest(dest) {}

  // for receiving (msg is allocated by the communicator)
  Work(int t, void* msg, int count)
      : type(t), msg(msg), count(count), dest(-1) {}

  virtual ~Work() {}

  int type;
  void* msg;
  int count;
  int dest;
};

template <typename PayloadT, typename HeaderT>
struct WorkRecvMsg {
  static void decode(void* msg, HeaderT** header, PayloadT** payload) {
    HeaderT* buf = (HeaderT*)msg;
    *header = buf;
    *payload = buf->payload_count ? (PayloadT*)(buf + 1) : nullptr;
  }
};

template <typename PayloadT, typename HeaderT>
class WorkSendMsg : public Work {
 public:
  WorkSendMsg(int work_type, const HeaderT& header, int dest)
      : Work(work_type, dest) {
    //
    header_ = header;
    allocMsg(header);
  }

  // for receiving
  WorkSendMsg(int work_type, void* msg, int msg_count)
      : Work(work_type, msg, msg_count) {
    decode();
  }

  virtual ~WorkSendMsg() {
    // DEBUG
    // std::cout << RANK_THREAD << "destruct SendRays\n";
    FreeAligned(msg);
  }

 private:
  // output:
  //   header_
  //   payload_
  void decode() {
    HeaderT* buf = (HeaderT*)msg;
    header_ = *buf;
    payload_ = header_.payload_count ? (PayloadT*)(buf + 1) : nullptr;
  }

 public:
  // output:
  //   header.payload_counter
  //   Work::count
  void updateMetaData(int payload_count, int new_block_count) {
    // update decoded header
    header_.payload_count = payload_count;
    // update msg string
    HeaderT* buf = (HeaderT*)msg;
    // CHECK_LT(payload_count, INT_MAX);
    buf->payload_count = payload_count;
    buf->new_block_count = new_block_count;
    Work::count = msgWordCount(payload_count);
  }

  void updatePayloadCount(int payload_count) {
    // update decoded header
    header_.payload_count = payload_count;
    // update msg string
    HeaderT* buf = (HeaderT*)msg;
    buf->payload_count = payload_count;
    Work::count = msgWordCount(payload_count);
  }

  void updateNumNewBlocks(int new_block_count) {
    // update decoded header
    header_.new_block_count = new_block_count;
    // update msg string
    HeaderT* buf = (HeaderT*)msg;
    buf->new_block_count = new_block_count;
  }

 public:
  // message
  //   HeaderT
  //   PayloadT[]
  const HeaderT& getHeader() const { return header_; }
  PayloadT* getPayload() { return payload_; }

 private:
  // decoded data
  HeaderT header_;
  PayloadT* payload_;

 private:
  // output:
  //   Work::msg
  //   Work::count
  //   this->payload_
  void allocMsg(const HeaderT& header) {
    std::size_t bytes =
        sizeof(HeaderT) + (header.payload_count * sizeof(PayloadT));

    Work::msg = AllocMsg(bytes, &(Work::count));

    HeaderT* tmp = (HeaderT*)msg;
    *tmp = header;
    payload_ = (PayloadT*)(tmp + 1);
  }

 private:
  int msgWordCount(std::size_t payload_count) {
    return MsgWordCount(getBytes(payload_count));
  }

 private:
  std::size_t getBytes(std::size_t payload_count) {
    return (sizeof(HeaderT) + (sizeof(PayloadT) * payload_count));
  }
};

}  // namespace insitu
}  // namespace spray
