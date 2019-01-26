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
#include <cstdint>
#include <list>
#include <queue>
#include <vector>

#include "glog/logging.h"
#include "pbrt/memory.h"

#include "insitu/insitu_work.h"
#include "insitu/insitu_work_stats.h"
#include "utils/comm.h"
#include "utils/profiler_util.h"

namespace spray {
class MemoryArena;

namespace insitu {

class TileGen;
class InsituPartition;
class WorkStats;

struct Work;
struct SendPixelsWork;

class DefaultReceiver {
  typedef std::queue<msg_word_t*> MessageQ;

 public:
  DefaultReceiver() {}
  DefaultReceiver(MessageQ* radiance_q, MessageQ* shadow_q)
      : rq_(radiance_q), sq_(shadow_q) {}
  void set(MessageQ* radiance_q, MessageQ* shadow_q) {
    rq_ = radiance_q;
    sq_ = shadow_q;
  }

  void operator()(int tag, msg_word_t* msg) {
    if (tag == WORK_SEND_RADS) {
      rq_->push(msg);

    } else if (tag == WORK_SEND_SHADS) {
      sq_->push(msg);

    } else {
      LOG(FATAL) << "unknown mpi tag : " << tag;
    }
  }

 private:
  MessageQ* rq_;
  MessageQ* sq_;
};

template <typename ReceiverT>
class Comm {
 public:
  void init();
  void run(const WorkStats& work_stats, MemoryArena* mem, ReceiverT* receiver);

  bool emptySendQ() const { return send_q_.empty(); }
  void pushSendQ(Work* work) { send_q_.push(work); }

 private:
  void mpiIsendWords(Work* work, void* msg, int count, int dest, int tag);

  void serveRecv(const MPI_Status& status, MemoryArena* mem,
                 ReceiverT* receiver);

 public:
  void waitForSend();

 private:
  void testMpiRqsts();

 private:
  struct MpiRequest {
    MPI_Request req;
    Work* work;
  };

  std::list<MpiRequest> mpi_requests_;
  std::queue<Work*> send_q_;
};

}  // namespace insitu
}  // namespace spray

#define SPRAY_INSITU_COMM_INL
#include "insitu/insitu_comm.inl"
#undef SPRAY_INSITU_COMM_INL

