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

#include <mpi.h>

#include "pbrt/memory.h"
#include "utils/comm.h"

namespace spray {
namespace insitu {

class VBuf;

struct MsgHeader {
  int domain_id;
  int64_t payload_count;
};

class Work {
 public:
  enum Type {
    SEND_SHADOW_RAYS,
    SEND_RADIANCE_RAYS,
    MSG_TERMINATE,
  };

  Work() {}
  Work(int type) : type_(type) {}

  // Work(int type) : type_(type), msg(nullptr), count(-1), dest(-1) {}

  // // for sending (msg is allocated by the derived)
  // Work(int t, int dest) : type(t), msg(nullptr), count(-1), dest(dest) {}

  // // for receiving (msg is allocated by the communicator)
  // Work(int t, void* msg, int count)
  //     : type(t), msg(msg), count(count), dest(-1) {}

  virtual ~Work() {}

  int getType() const { return type_; }

 protected:
  void setType(int type) { type_ = type; }
  int type_;
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
  WorkSendMsg() {}
  virtual ~WorkSendMsg() {}

  void allocate(int work_type, const HeaderT& header, int dest,
                MemoryArena* mem) {
    setType(work_type);
    header_ = header;
    dest_ = dest;
    allocMsg(header, mem);
  }

  const HeaderT& getHeader() const { return header_; }
  PayloadT* getPayload() { return payload_; }

  void isend(MPI_Request* rqst) {
    MPI_Isend(msg_, count_, MPI_WORD_T, dest_, type_, MPI_COMM_WORLD, rqst);
  }

 private:
  // update: msg_, count_, payload_
  void allocMsg(const HeaderT& header, MemoryArena* mem) {
    std::size_t bytes =
        sizeof(HeaderT) + (header.payload_count * sizeof(PayloadT));

    allocMem(bytes, mem);

    HeaderT* tmp = (HeaderT*)msg_;
    *tmp = header;
    payload_ = (PayloadT*)(tmp + 1);
  }

  void allocMem(std::size_t bytes, MemoryArena* mem) {
    std::size_t word_size = sizeof(msg_word_t);
    std::size_t word_count = (bytes + word_size - 1) / word_size;
    CHECK_LT(word_count, INT_MAX);
    msg_ = mem->Alloc<msg_word_t>(word_count, false);
    CHECK_NOTNULL(msg_);
    count_ = static_cast<int>(word_count);
  }

 private:
  HeaderT header_;
  PayloadT* payload_;
  void* msg_;
  int count_;
  int dest_;
};

}  // namespace insitu
}  // namespace spray
