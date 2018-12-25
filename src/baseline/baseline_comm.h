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
#include <cstring>
#include <queue>
#include <vector>

#include "glog/logging.h"
#include "pbrt/memory.h"
#include "pbrt/memory3.h"

#include "baseline/baseline_ray.h"
#include "baseline/baseline_ray_stats.h"
#include "baseline/baseline_schedulers.h"
#include "partition/arena_queue.h"
#include "partition/block_buffer.h"
#include "render/spray.h"
#include "utils/comm.h"
#include "utils/profiler_util.h"
#include "utils/util.h"

namespace spray {
namespace baseline {

// DRayQItem to send buffer (DRay)
// copy Q items to message (i.e. copy from arena to send buffer)
class DRayOutgoingCopier {
 public:
  void operator()(void* ray_buf, void* qitem_buf, std::size_t n) {
    DRayQItem* qitems = (DRayQItem*)qitem_buf;
    DRay* rays = (DRay*)ray_buf;
#ifdef SPRAY_GLOG_CHECK
    CHECK_NOTNULL(qitems);
    CHECK_NOTNULL(rays);
#endif
    for (std::size_t i = 0; i < n; ++i) {
      std::memcpy(&rays[i], qitems[i].ray, sizeof(DRay));
#ifdef SPRAY_GLOG_CHECK
      DRay* r = &rays[i];
      if (DRayUtil::getShadow(*r) == 0x00000000) {
        CHECK_GE(r->domid, 0) << *r;
      }
#endif
    }
  }
};

// DRay to queue item (DRayQItem)
// copy from recv buffer to memory arena)
class DRayIncomingCopier {
 public:
  void operator()(void* qitem_buf, void* ray_buf, std::size_t n) {
    //
    DRayQItem* qitems = (DRayQItem*)qitem_buf;
    DRay* rays = (DRay*)ray_buf;
#ifdef SPRAY_GLOG_CHECK
    CHECK_NOTNULL(qitems);
    CHECK_NOTNULL(rays);
#endif
    for (std::size_t i = 0; i < n; ++i) {
      qitems[i].ray = &rays[i];
#ifdef SPRAY_GLOG_CHECK
      DRay* r = &rays[i];
      if (DRayUtil::getShadow(*r) == 0x00000000) {
        CHECK_GE(r->domid, 0) << *r;
      }
#endif
    }
  }
};

template <typename T>
struct SPRAY_ALIGN(16) BufferInfo {
  void reset() {
    base = nullptr;
    size = 0;
    capacity = 0;
  }
  T* base;
  std::size_t size;
  std::size_t capacity;
};

template <class QItemT, class MessageT, class OutgoingCopierT,
          class IncomingCopierT>
class Comm {
 public:
  Comm()
      : send_rbuf_(nullptr),
        send_sbuf_(nullptr),
        recv_rbuf_(nullptr),
        recv_sbuf_(nullptr),
        schedule_(nullptr) {}

 public:
  void init(int ndomains, const std::vector<RayCount>& sched,
            BlockBuffer* send_rbuf, BlockBuffer* send_sbuf,
            BlockBuffer* recv_rbuf, BlockBuffer* recv_sbuf);

  void run();

  void waitForSend();
  void resetMsgSendArena() { send_arena_.Reset(); }
  void resetMsgRecvArena() { recv_arena_.Reset(); }
  void resetQitemRecvArena() { qitem_recv_arena_.Reset(); }

  template <typename T>
  T* allocRecvArena(std::size_t size) {
    return recv_arena_.Alloc<T>(size, false);
  }

 private:
  std::size_t recv();
  std::size_t send(int id, int dest);
  void setup();

 private:
  MemoryArena recv_arena_;  // message arena
  MemoryArena send_arena_;  // message arena

  MemoryArena qitem_recv_arena_;

  BlockBuffer* send_rbuf_;
  BlockBuffer* send_sbuf_;

  BlockBuffer* recv_rbuf_;
  BlockBuffer* recv_sbuf_;

  //! Shared ray statistics maintained in the tracer.
  const std::vector<RayCount>* schedule_;

  int ndomains_;  //!< Number of domains.

  BufferInfo<MessageT> outgoing_msgbuf_;
  BufferInfo<MessageT> incoming_msgbuf_;

  BufferInfo<QItemT> incoming_qbuf_;
  BufferInfo<QItemT> incoming_shadow_qbuf_;

  std::queue<MPI_Request> isend_q_;

  OutgoingCopierT senderFunctor_;
  IncomingCopierT recverFunctor_;
};

}  // namespace baseline
}  // namespace spray

#define SPRAY_BASELINE_COMM_INL
#include "baseline/baseline_comm.inl"
#undef SPRAY_BASELINE_COMM_INL

