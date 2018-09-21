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
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "utils/timer.h"

namespace spray {

enum TimerNames {
  TIMER_TOTAL = 0,
  TIMER_LOAD,
  TIMER_SYNC_RAYS,
  TIMER_SYNC_SCHED,
  TIMER_SYNC_VBUF,
  TIMER_SYNC_IMAGE,
  TIMER_COUNT
};

enum CounterNames {
  COUNTER_RAYS_SENT = 0,
  COUNTER_RAYS_SPAWNED,
  COUNTER_RAYS_TESTED,
  COUNTER_COUNT
};

class Profiler {
 private:
  struct Stats {
    void reset() {
      min = std::numeric_limits<double>::max();
      max = 0.0;
      avg = 0.0;
    }
    void aggregate(double v) {
      if (v < min) min = v;
      if (v > max) max = v;
      avg += v;
    }
    void average(int n) { avg /= (double)n; }

    double min;
    double max;
    double avg;
  };
  struct CounterStats {
    void reset() {
      min = std::numeric_limits<uint64_t>::max();
      max = (uint64_t)0;
      avg = (uint64_t)0;
      sum = (uint64_t)0;
    }
    void aggregate(uint64_t v) {
      if (v < min) min = v;
      if (v > max) max = v;
      avg += v;
      sum += v;
    }
    void average(int n) {
      // if (avg > std::numeric_limits<uint64_t>::epsilon())
      avg /= (uint64_t)n;
      // else
      //  avg = 0.0;
    }
    uint64_t min;
    uint64_t max;
    uint64_t avg;
    uint64_t sum;
  };

 public:
  void init() {
    // resize timer vectors
    timers_.resize(TIMER_COUNT);
    counters_.resize(COUNTER_COUNT);

    // reset timers and data struct
    reset();
  }

  //! Synchronizes all the measured values across the cluster and prints them on
  //! screen.
  void print(int64_t nframes);

  void start(int timer) { timers_[timer].start(); }
  void stop(int timer) { timers_[timer].stop(); }

  void agg(int counter, uint64_t v) { counters_[counter] += v; }

  static std::map<int, std::string> timer_names;
  static std::map<int, std::string> counter_names;

  void reset() {
    for (auto& timer : timers_) {
      timer.reset();
    }
    for (auto& c : counters_) {
      c = 0;
    }
  }

 private:
  void printTimers(int rank, const double* timers, int64_t nframes);
  void printCounters(int rank, const uint64_t* counters, int64_t nframes);
  void printStats(const std::vector<Stats>& stats, int64_t nframes);
  void printCounterStats(const std::vector<CounterStats>& stats,
                         int64_t nframes);

  void aggStats(int rank, const double* timers, std::vector<Stats>* stats);
  void aggCounterStats(int rank, const uint64_t* counters,
                       std::vector<CounterStats>* stats);

  std::vector<Timer> timers_;
  std::vector<uint64_t> counters_;
};

}  // namespace spray
