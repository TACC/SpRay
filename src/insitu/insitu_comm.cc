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

#include "insitu/insitu_comm.h"

#include "glog/logging.h"
#include "pbrt/memory.h"

#include "display/image.h"
#include "insitu/insitu_ray.h"
#include "insitu/insitu_vbuf.h"
#include "insitu/insitu_work.h"
#include "insitu/insitu_work_stats.h"
#include "utils/profiler_util.h"

namespace spray {
namespace insitu {

void Comm::mpiIsendWords(Work* work, void* msg, int count, int dest, int tag) {
  MpiRequest isend_rqst;
  isend_rqst.work = work;

  mpi_requests_.push_back(isend_rqst);
  MpiRequest& r = mpi_requests_.back();

  MPI_Isend(msg, count, MPI_WORD_T, dest, tag, MPI_COMM_WORLD, &r.req);
}

void Comm::serveRecv(const MPI_Status& status, MemoryArena* mem_in,
                     std::queue<msg_word_t*>* recv_rq,
                     std::queue<msg_word_t*>* recv_sq) {
  int tag = status.MPI_TAG;
  int msg_count;
  MPI_Get_count(&status, MPI_WORD_T, &msg_count);

  msg_word_t* msg = mem_in->Alloc<msg_word_t>(msg_count);
  CHECK_NOTNULL(msg);

  MPI_Recv(msg, msg_count, MPI_WORD_T, status.MPI_SOURCE, status.MPI_TAG,
           MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  if (tag == WORK_SEND_RADS) {
    recv_rq->push(msg);

  } else if (tag == WORK_SEND_SHADS) {
    recv_sq->push(msg);

  } else {
    LOG(FATAL) << "unknown mpi tag : " << tag;
  }
}

void Comm::run(WorkStats* work_stats, MemoryArena* mem_in,
               std::queue<msg_word_t*>* recv_rq,
               std::queue<msg_word_t*>* recv_sq) {
  MPI_Status status;
  int flag;

  int num_blocks_recved = 0;
  bool recv_done = work_stats->recvDone(num_blocks_recved);

#ifdef SPRAY_TIMING
  spray::tStart(spray::TIMER_SYNC_RAYS);
#endif
  while (1) {
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);

#ifdef SPRAY_GLOG_CHECK
    CHECK(!(recv_done && flag));
#endif
    if (!recv_done && flag) {
      serveRecv(status, mem_in, recv_rq, recv_sq);
      ++num_blocks_recved;
      recv_done = work_stats->recvDone(num_blocks_recved);
    }

    if (!send_q_.empty()) {
      Work* w = send_q_.front();
      send_q_.pop();

      int t = w->type;
#ifdef SPRAY_GLOG_CHECK
      CHECK((t == WORK_SEND_SHADS) || (t == WORK_SEND_RADS)) << t;
#endif
      mpiIsendWords(w, w->msg, w->count, w->dest, w->type);
    } else if (recv_done) {
      break;
    }
  }
#ifdef SPRAY_TIMING
  spray::tStop(spray::TIMER_SYNC_RAYS);
#endif
}

void Comm::waitForSend() {
  for (auto it = mpi_requests_.begin(); it != mpi_requests_.end(); ++it) {
    MPI_Wait(&it->req, MPI_STATUS_IGNORE);
    delete it->work;
  }
  mpi_requests_.clear();
}

void Comm::testMpiRqsts() {
  int flag;
  for (auto it = mpi_requests_.begin(); it != mpi_requests_.end();) {
    MPI_Test(&it->req, &flag, MPI_STATUS_IGNORE);
    if (flag) {
      delete it->work;
      it = mpi_requests_.erase(it);
    } else {
      ++it;
    }
  }
}

}  // namespace insitu
}  // namespace spray

