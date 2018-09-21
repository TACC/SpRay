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

#include "utils/comm.h"

namespace spray {
class MemoryArena;

namespace insitu {

class TileGen;
class InsituPartition;
class WorkStats;

struct Work;
struct SendPixelsWork;

class Comm {
 public:
  void init();
  void run(WorkStats* work_stats, spray::MemoryArena* mem_in,
           std::queue<msg_word_t*>* recv_rq, std::queue<msg_word_t*>* recv_sq);

  bool emptySendQ() const { return send_q_.empty(); }
  void pushSendQ(Work* work) { send_q_.push(work); }

 private:
  void mpiIsendWords(Work* work, void* msg, int count, int dest, int tag);

  void serveRecv(const MPI_Status& status, MemoryArena* mem_in,
                 std::queue<msg_word_t*>* recv_rq,
                 std::queue<msg_word_t*>* recv_sq);

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

