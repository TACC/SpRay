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
#include <omp.h>
#include <cstdint>

#include "renderers/spray.h"
#include "utils/profiler.h"

namespace spray {

inline void updateFrameNumber(unsigned i) { global_num_frames = i; }

inline void tStart(int timer) { global_profiler.start(timer); }
inline void tStop(int timer) { global_profiler.stop(timer); }

inline void tStartMPI(int timer) {
  MPI_Barrier(MPI_COMM_WORLD);
  global_profiler.start(timer);
}

inline void tAgg(int counter, uint64_t n) { global_profiler.agg(counter, n); }

inline void tReset() { global_profiler.reset(); }

inline void tPrint(int64_t nframes) { global_profiler.print(nframes); }

}  // namespace spray
