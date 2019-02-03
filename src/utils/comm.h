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
#include <climits>
#include <queue>

#include "glog/logging.h"
#include "pbrt/memory.h"

#include "render/spray.h"

namespace spray {

struct DomainSchedulerState;

struct Tile;
struct TileSchedule;

typedef uint32_t msg_word_t;
#define MPI_WORD_T MPI_UINT32_T
#define MSG_WORD_SIZE sizeof(msg_word_t)

class TileGen;
class InsituPartition;
class WorkStats;

inline msg_word_t *AllocMsg(std::size_t bytes) {
  std::size_t count = (bytes + MSG_WORD_SIZE - 1) / MSG_WORD_SIZE;
  CHECK_LT(count, INT_MAX);
  msg_word_t *msg = AllocAligned<msg_word_t>(count);
  CHECK_NOTNULL(msg);
  // *msg_word_count = (int)count;
  return msg;
}

inline msg_word_t *AllocMsg(std::size_t bytes, int *msg_word_count) {
  std::size_t count = (bytes + MSG_WORD_SIZE - 1) / MSG_WORD_SIZE;
  CHECK_LT(count, INT_MAX);
  msg_word_t *msg = AllocAligned<msg_word_t>(count);
  CHECK_NOTNULL(msg);
  *msg_word_count = (int)count;
  return msg;
}

inline int MsgWordCount(std::size_t bytes) {
  std::size_t count = (bytes + MSG_WORD_SIZE - 1) / MSG_WORD_SIZE;
  CHECK_LT(count, INT_MAX);
  return count;
}

class Comm {
 public:
  Comm() {}
  virtual ~Comm() {}

  virtual void waitIsends() = 0;
  virtual bool isIsendDone() = 0;

  virtual void isendRays(const DomainSchedulerState &state) = 0;
  virtual unsigned recvRays() = 0;

  virtual void scatterTiles(std::queue<Tile> *tiles, TileSchedule *ts) = 0;

  // helpers

  inline static void bcast(void *data, int count, MPI_Datatype datatype,
                           int root = SPRAY_ROOT_PROCESS,
                           MPI_Comm comm = MPI_COMM_WORLD) {
    MPI_Bcast(data, count, datatype, root, comm);
  }

  inline static void allReduce(const void *sendbuf, void *recvbuf, int count,
                               MPI_Datatype datatype, MPI_Op op = MPI_SUM,
                               MPI_Comm comm = MPI_COMM_WORLD) {
    MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
  }

  inline static void gather(const void *sendbuf, int sendcount,
                            MPI_Datatype sendtype, void *recvbuf, int recvcount,
                            MPI_Datatype recvtype,
                            int root = SPRAY_ROOT_PROCESS,
                            MPI_Comm comm = MPI_COMM_WORLD) {
    MPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root,
               comm);
  }

  inline static void scatter(const void *sendbuf, int sendcount,
                             MPI_Datatype sendtype, void *recvbuf,
                             int recvcount, MPI_Datatype recvtype,
                             int root = SPRAY_ROOT_PROCESS,
                             MPI_Comm comm = MPI_COMM_WORLD) {
    MPI_Scatter(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype,
                root, comm);
  }

  inline static void allgather(const void *sendbuf, int sendcount,
                               MPI_Datatype sendtype, void *recvbuf,
                               int recvcount, MPI_Datatype recvtype,
                               MPI_Comm comm = MPI_COMM_WORLD) {
    MPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype,
                  comm);
  }

  inline static bool isRootProcess() {
    return global_mpi_comm.rank == SPRAY_ROOT_PROCESS;
  }

  inline static bool isRoot(int rank) { return (rank == SPRAY_ROOT_PROCESS); }
  inline static bool isRootProcess(int rank) {
    return (rank == SPRAY_ROOT_PROCESS);
  }

  inline static bool isMyProcess(int i) { return i == global_mpi_comm.rank; }

  inline static bool isMultiNode() { return global_mpi_comm.size > 1; }

  inline static int size() { return global_mpi_comm.size; }
  inline static int rank() { return global_mpi_comm.rank; }
};

}  // namespace spray

namespace spray {
namespace mpi {

enum Rank { kRoot = 0, kChild };

inline bool isMultiProcess() { return global_mpi_comm.size > 1; }
inline bool isSingleProcess() { return global_mpi_comm.size == 1; }
inline bool isRootProcess() { return (global_mpi_comm.rank == kRoot); }
inline int root() { return kRoot; }
inline int size() { return global_mpi_comm.size; }
inline int rank() { return global_mpi_comm.rank; }
inline int worldSize() { return global_mpi_comm.size; }
inline int worldRank() { return global_mpi_comm.rank; }
inline MPI_Comm comm() { return MPI_COMM_WORLD; }

}  // namespace mpi
}  // namespace spray

