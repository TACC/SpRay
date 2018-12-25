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

#include "baseline/baseline_schedulers.h"

#include <algorithm>

#include "glog/logging.h"

#include "baseline/baseline_comm.h"
#include "baseline/baseline_qstats.h"
#include "render/data_partition.h"
#include "utils/profiler_util.h"

namespace spray {
namespace baseline {

LoadAnyOnceSched::LoadAnyOnceSched() : ndomains_(0) {}

void LoadAnyOnceSched::init(int ndomains, const InsituPartition& partition) {
  ndomains_ = ndomains;

  if (mpi::size() > 1) {
    state_.resize(ndomains);
  }

  schedule_.resize(mpi::size());

  ray_counts_.resize(ndomains);
}

void LoadAnyOnceSched::schedule(const QStats& stats) {
  if (mpi::size() > 1) {
    scheduleMulti(stats);
  } else {
    scheduleImage(stats);
  }
#ifdef SPRAY_GLOG_CHECK
  CHECK_EQ(id_, schedule_[mpi::rank()].id);
#endif
}

void LoadAnyOnceSched::scheduleImage(const QStats& stats) {
  // populate local state
  for (int i = 0; i < ndomains_; ++i) {
    ray_counts_[i].id = i;
    ray_counts_[i].count = stats.get(i);
  }

  std::sort(
      ray_counts_.begin(), ray_counts_.end(),
      [](const RayCount& i, const RayCount& j) { return i.count > j.count; });

  id_ = ray_counts_[0].id;
  done_ = (ray_counts_[0].count == 0);

  schedule_[0].id = id_;
  schedule_[0].count = ray_counts_[0].count;
}

void LoadAnyOnceSched::scheduleMulti(const QStats& stats) {
  // populate local state
  for (int i = 0; i < ndomains_; ++i) {
    state_[i] = stats.get(i);
  }

  // aggregate local states

  MPI_Allreduce(MPI_IN_PLACE, &state_[0], ndomains_, MPI_INT64_T, MPI_SUM,
                MPI_COMM_WORLD);

  // aggregated ray counts

  for (int i = 0; i < ndomains_; ++i) {
    ray_counts_[i].id = i;
    ray_counts_[i].count = state_[i];
  }

  std::sort(
      ray_counts_.begin(), ray_counts_.end(),
      [](const RayCount& i, const RayCount& j) { return i.count > j.count; });

  // distribute domains in round robin based on sorted ray counts
  // i.e. rank 0 gets the domain with most rays in it

  unsigned domid = 0;
  int num_ranks = mpi::size();
  for (int rank = 0; rank < num_ranks; ++rank) {
    RayCount& s = schedule_[rank];

    if (domid < ndomains_) {
      const RayCount& c = ray_counts_[domid];
      if (c.count > 0) {
        s.id = c.id;
        s.count = c.count;
      } else {
        s.id = INT_MAX;
        s.count = 0;
      }
      ++domid;
    } else {
      s.id = INT_MAX;
      s.count = 0;
    }
  }

  id_ = schedule_[mpi::rank()].id;
  // root process always gets scheduled first
  // so if its assigned workload is zero, then done is true
  done_ = (ray_counts_[0].count == 0);

#ifdef SPRAY_GLOG_CHECK
  if (done_) {
    int64_t count = 0;
    for (std::size_t i = 0; i < schedule_.size(); ++i) {
      count += schedule_[i].count;
    }
    CHECK_EQ(count, 0);
  }
#endif
}

void LoadAnyOnceInsituSched::schedule(const QStats& stats) {
  if (mpi::size() > 1) {
    scheduleMulti(stats);
  } else {
    scheduleImage(stats);
  }
#ifdef SPRAY_GLOG_CHECK
  CHECK_EQ(id_, schedule_[mpi::rank()].id);
#endif
}

void LoadAnyOnceInsituSched::scheduleMulti(const QStats& stats) {
  // populate local state
  for (int i = 0; i < ndomains_; ++i) {
    state_[i] = stats.get(i);
  }

// aggregate local states
#ifdef SPRAY_TIMING
  spray::tStartMPI(spray::TIMER_SYNC_SCHED);
#endif
  MPI_Allreduce(MPI_IN_PLACE, &state_[0], ndomains_, MPI_INT64_T, MPI_SUM,
                MPI_COMM_WORLD);
#ifdef SPRAY_TIMING
  spray::tStop(spray::TIMER_SYNC_SCHED);
#endif

  // aggregated ray counts

  int num_ranks = mpi::size();

  int64_t count = 0;
  for (auto s : state_) {
    count += s;
  }
  done_ = (count == 0);

  if (done_) {
    id_ = INT_MAX;
    for (auto& sched : schedule_) {
      sched.id = INT_MAX;
      sched.count = 0;
    }
  } else {
    for (int rank = 0; rank < num_ranks; ++rank) {
      const std::list<int>& domains = partition_->getDomains(rank);

      int64_t max = 0;
      int max_id = INT_MAX;

      for (int id : domains) {
        int64_t ray_count = state_[id];

        if (ray_count > max) {
          max = ray_count;
          max_id = id;
        }
      }

      schedule_[rank].id = max_id;
      schedule_[rank].count = max;
    }

    id_ = schedule_[mpi::rank()].id;
  }

#ifdef SPRAY_GLOG_CHECK
  if (done_) {
    int64_t count = 0;
    for (std::size_t i = 0; i < schedule_.size(); ++i) {
      count += schedule_[i].count;
    }
    CHECK_EQ(count, 0);

    int done_count = 0;
    int64_t ray_count = 0;

    for (auto& s : schedule_) {
      done_count += (s.id >= ndomains_);
      ray_count += s.count;
    }
    if (done_) {
      CHECK_EQ(done_count, num_ranks);
      CHECK_EQ(ray_count, 0);
    } else {
      CHECK_LT(done_count, num_ranks);
      CHECK_GT(ray_count, 0);
    }
  }
#endif
}

}  // namespace baseline
}  // namespace spray

