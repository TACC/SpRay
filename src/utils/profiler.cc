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

#include "utils/profiler.h"

#include <mpi.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "glog/logging.h"

#include "render/spray.h"
#include "utils/comm.h"

namespace spray {

std::map<int, std::string> Profiler::timer_names = {
    {TIMER_TOTAL, "total"},         {TIMER_LOAD, "load"},
    {TIMER_SYNC_RAYS, "sync_rays"}, {TIMER_SYNC_SCHED, "sync_sched"},
    {TIMER_SYNC_VBUF, "sync_vbuf"}, {TIMER_SYNC_IMAGE, "sync_image"}};

std::map<int, std::string> Profiler::counter_names = {
    {COUNTER_RAYS_SENT, "rays_sent"},
    {COUNTER_RAYS_SPAWNED, "rays_spawned"},
    {COUNTER_RAYS_TESTED, "rays_tested"}};

void Profiler::aggStats(int rank, const double* timers,
                        std::vector<Stats>* stats) {
  for (int t = 0; t < TIMER_COUNT; ++t) {
    (*stats)[t].aggregate(timers[t]);
  }
}

void Profiler::aggCounterStats(int rank, const uint64_t* counters,
                               std::vector<CounterStats>* counter_stats) {
  for (int t = 0; t < COUNTER_COUNT; ++t) {
    (*counter_stats)[t].aggregate(counters[t]);
  }
}

void Profiler::printStats(const std::vector<Stats>& stats, int64_t nframes) {
  std::cout << "[STATS] " << std::left << std::setw(20) << "items"
            << ": " << std::left << std::setw(12) << "min"
            << " " << std::left << std::setw(12) << "avg"
            << " " << std::left << std::setw(12) << "max"
            << "\n";
  for (std::size_t i = 0; i < stats.size(); ++i) {
    const auto& s = stats[i];
    std::cout << "[STATS] " << std::left << std::setw(20)
              << Profiler::timer_names[i].c_str() << ": " << std::left
              << std::setw(12) << 1000.0 * s.min / nframes << " " << std::left
              << std::setw(12) << 1000.0 * s.avg / nframes << " " << std::left
              << std::setw(12) << 1000.0 * s.max / nframes << " ms/frame\n";
  }
}

void Profiler::printCounterStats(const std::vector<CounterStats>& stats,
                                 int64_t nframes) {
  std::cout << "[STATS] " << std::left << std::setw(20) << "items"
            << ": " << std::left << std::setw(12) << "min"
            << " " << std::left << std::setw(12) << "avg"
            << " " << std::left << std::setw(12) << "max"
            << " " << std::left << std::setw(12) << "sum"
            << "\n";
  for (std::size_t i = 0; i < stats.size(); ++i) {
    const auto& s = stats[i];
    std::cout << "[STATS] " << std::left << std::setw(20)
              << Profiler::counter_names[i].c_str() << ": " << std::left
              << std::setw(12) << s.min / nframes << " " << std::left
              << std::setw(12) << s.avg / nframes << " " << std::left
              << std::setw(12) << s.max / nframes << " " << std::left
              << std::setw(12) << s.sum / nframes << " per frame\n";
  }
}
void Profiler::printTimers(int rank, const double* timers, int64_t nframes) {
  for (int t = 0; t < TIMER_COUNT; ++t) {
    // in seconds
    double total_time = timers[t];
    double frame_time = total_time / nframes;
    std::cout << "[TIME] [PROC " << rank << "] " << std::left << std::setw(20)
              << Profiler::timer_names[t].c_str() << ": " << std::left
              << std::setw(12) << 1000.0 * frame_time << "ms/frame\n";
  }
  double total_time_r = timers[TIMER_TOTAL];
  double frame_time_r = total_time_r / nframes;
  std::cout << "[TIME] [PROC " << rank << "] " << std::left << std::setw(20)
            << "frame rate"
            << ": " << std::left << std::setw(12) << 1.0 / frame_time_r
            << "fps\n";
  std::cout << "[TIME] [PROC " << rank << "] " << std::left << std::setw(20)
            << "frame count"
            << ": " << std::left << std::setw(12) << nframes << "frames\n";
}

void Profiler::printCounters(int rank, const uint64_t* counters,
                             int64_t nframes) {
  for (int t = 0; t < COUNTER_COUNT; ++t) {
    // in seconds
    double total_count = counters[t];
    double perframe_count = total_count / nframes;
    std::cout << "[COUNT] [PROC " << rank << "] " << std::left << std::setw(20)
              << Profiler::counter_names[t].c_str() << ": " << std::left
              << std::setw(12) << total_count << "total\n";
    std::cout << "[COUNT] [PROC " << rank << "] " << std::left << std::setw(20)
              << Profiler::counter_names[t].c_str() << ": " << std::left
              << std::setw(12) << perframe_count << "per frame\n";
  }
}

// this function need not be efficient.
// called once only after rendering all the frames
void Profiler::print(int64_t nframes) {
  int nranks = mpi::size();

  // render timers

  std::vector<double> timers_recvbuf(nranks * TIMER_COUNT);
  std::vector<double> timers_sendbuf(TIMER_COUNT);

  // populate sendbuf with timer values
  for (int i = 0; i < TIMER_COUNT; ++i) {
    timers_sendbuf[i] = timers_[i].getTotalTime();
  }

  MPI_Gather(&timers_sendbuf[0], TIMER_COUNT, MPI_DOUBLE, &timers_recvbuf[0],
             TIMER_COUNT, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // counters
  //
  std::vector<uint64_t> counters_recvbuf(nranks * COUNTER_COUNT);
  std::vector<uint64_t> counters_sendbuf(COUNTER_COUNT);

  // populate sendbuf with counter values
  for (int i = 0; i < COUNTER_COUNT; ++i) {
    counters_sendbuf[i] = counters_[i];
  }

  MPI_Gather(&counters_sendbuf[0], COUNTER_COUNT, MPI_UINT64_T,
             &counters_recvbuf[0], COUNTER_COUNT, MPI_UINT64_T, 0,
             MPI_COMM_WORLD);

  int rank = mpi::rank();
  if (rank > 0) return;

  std::vector<Stats> stats(TIMER_COUNT);
  std::vector<CounterStats> counter_stats(COUNTER_COUNT);
  for (auto& s : stats) s.reset();
  for (auto& cs : counter_stats) cs.reset();

  for (int i = 0; i < nranks; ++i) {
    printf("\n==============[PROC %d]===============\n", i);

    printf("\n");
    aggStats(i, &timers_recvbuf[i * TIMER_COUNT], &stats);
    printTimers(i, &timers_recvbuf[i * TIMER_COUNT], nframes);

    aggCounterStats(i, &counters_recvbuf[i * COUNTER_COUNT], &counter_stats);
    printCounters(i, &counters_recvbuf[i * COUNTER_COUNT], nframes);
  }
  printf("\n");
  printf("\n");
  printf("\n==============[PROC 0-%d]===============\n", nranks - 1);

  for (auto& s : stats) {
    s.average(nranks);
  }
  printStats(stats, nframes);

  for (auto& cs : counter_stats) {
    cs.average(nranks);
  }
  printCounterStats(counter_stats, nframes);
}

}  // namespace spray
