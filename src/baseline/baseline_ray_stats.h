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

#include <algorithm>
#include <cstdint>
#include <string>

#include "glog/logging.h"
#include "pbrt/memory.h"

#include "cmake_config.h"

#include "render/spray.h"
#include "utils/comm.h"

namespace spray {

class InsituPartition;

namespace baseline {

struct SPRAY_ALIGN(16) Score {
  int id;
  int64_t score;
};

class RayStats {
 public:
  RayStats()
      : buf_(nullptr),
        scores_(nullptr),
        counts_(nullptr),
        schedule_(nullptr),
        stats_(nullptr),
        scores_and_counts_(nullptr),
        counts_and_schedule_(nullptr) {}

  virtual ~RayStats() { FreeAligned(buf_); }

 public:
  /**
   * Resizes internal data structures.
   * \param [in] ndomains Number of domains.
   * \param [in] max_domain_depth Maximum domain depth.
   * \param [in] allocate_stats_only Flag to allocate buffer for statistics
   * only.
   */
  void resize(int ndomains, int max_domain_depth, bool allocate_stats_only) {
    //
    ndomains_ = ndomains;
    int stats_size = ndomains * SPRAY_RAY_DOMAIN_LIST_SIZE;

    if (allocate_stats_only) {
      buf_ = AllocAligned<uint32_t>(stats_size);
      CHECK_NOTNULL(buf_);
      scores_ = nullptr;
      counts_ = nullptr;
      schedule_ = nullptr;
      stats_ = buf_;
      scores_and_counts_ = nullptr;
      counts_and_schedule_ = nullptr;
    } else {
      int buf_size = 3 * ndomains + stats_size;
      buf_ = AllocAligned<uint32_t>(buf_size);
      CHECK_NOTNULL(buf_);
      scores_ = buf_;
      counts_ = &buf_[ndomains];
      schedule_ = &buf_[2 * ndomains];
      stats_ = &buf_[3 * ndomains];
      scores_and_counts_ = scores_;
      counts_and_schedule_ = counts_;
      //
      sorted_scores_.resize(ndomains);
      id_to_rank_map_.resize(ndomains);
    }
  }

  /**
   * Resets statistics of all domains to zero.
   */
  void reset() {
    int size = SPRAY_RAY_DOMAIN_LIST_SIZE * ndomains_;
    for (int i = 0; i < size; ++i) {
      stats_[i] = 0;
    }
  }

  /**
   * Resets statistics of a domain to zero.
   * \param [in] id Domain ID.
   */
  void reset(int id) {
    for (int d = 0; d < SPRAY_RAY_DOMAIN_LIST_SIZE; ++d) {
      stats_[statsIndex(id, d)] = 0;
    }
  }

  /**
   * Add other statistics to this object for a given domain.
   * \param [in] id Domain ID.
   * \param [in] stats Ray statistics to add.
   */
  void addStats(int id, const RayStats& stats) {
#ifdef SPRAY_GLOG_CHECK
    CHECK_NOTNULL(stats_);
    CHECK_NOTNULL(stats.stats_);
#endif
    int offset = id * SPRAY_RAY_DOMAIN_LIST_SIZE;

    for (int d = 0; d < SPRAY_RAY_DOMAIN_LIST_SIZE; ++d) {
      int i = offset + d;
      stats_[i] += stats.stats_[i];
    }
  }

  /** Schedule rays in different domains based on statistics. */
  void schedule() {
    evaluateTraversalOrderAcrossCluster();
    bcastRayCountsAndSchedule();
  }

 protected:
  /**
   * Populates the queue sizes array with ray counts for all domains.
   * - The queue size for each domain extracted from statistics
   */
  void evaluateQSizesUsingStats() {
#ifdef SPRAY_GLOG_CHECK
    CHECK_NOTNULL(counts_);
#endif
    for (int id = 0; id < ndomains_; ++id) {
      counts_[id] = getRayCountFromStats(id);
    }
  }

 public:
  /** Increments the ray count of the domain at a given depth. */
  void increment(int id, int depth) {
#ifdef SPRAY_GLOG_CHECK
    CHECK_NOTNULL(stats_);
#endif
    ++stats_[statsIndex(id, depth)];
  }

  /** Decrements the ray count of the domain at a given depth. */
  void decrement(int id, int depth) {
#ifdef SPRAY_GLOG_CHECK
    CHECK_NOTNULL(stats_);
#endif
    int i = statsIndex(id, depth);
    auto count = stats_[i];
    CHECK_GT(count, 0);
    stats_[i] = count - 1;
  }

 public:
  /** Returns the number of domains. */
  int getNumDomains() const { return ndomains_; }

  /**
   * Returns the domain ID of a given traversal index into the traversal order.
   */
  int getDomainId(int schedule_index) const {
#ifdef SPRAY_GLOG_CHECK
    CHECK_NOTNULL(schedule_);
    CHECK_LT(schedule_index, ndomains_);
#endif
    return schedule_[schedule_index];
  }

  /**
   * Returns whether all queues are empty or not by scanning RayStats::counts_
   * array.
   */
  bool allQsEmpty() const {
#ifdef SPRAY_GLOG_CHECK
    CHECK_NOTNULL(counts_);
#endif
    for (int id = 0; id < ndomains_; ++id) {
      if (counts_[id] != 0) return false;
    }
    return true;
  }

  /** Returns true if ray count for a given domain is zero. */
  bool empty(int id) const {
#ifdef SPRAY_GLOG_CHECK
    CHECK_NOTNULL(counts_);
#endif
    return (counts_[id] == 0);
  }

 public:
  int rank(int id) const {
#ifdef SPRAY_GLOG_CHECK
    CHECK(id_to_rank_map_.size());
    CHECK_LT(id, ndomains_);
#endif
    return id_to_rank_map_[id];
  }

 public:
  /** Return the queue size (number of rays) of a given domain. */
  uint32_t getRayCount(int id) const {
#ifdef SPRAY_GLOG_CHECK
    CHECK_NOTNULL(counts_);
#endif
    return counts_[id];
  }

  /**
   * Extract ray count from the array of statistics.
   * \param [in] id Domain ID.
   */
  uint32_t getRayCountFromStats(int id) const {
#ifdef SPRAY_GLOG_CHECK
    CHECK_NOTNULL(stats_);
#endif
    uint32_t count = 0;
    int offset = statsOffset(id);
    for (int d = 0; d < SPRAY_RAY_DOMAIN_LIST_SIZE; ++d) {
      int idx = offset + d;
#ifdef SPRAY_GLOG_CHECK
      CHECK_LT(idx, ndomains_ * SPRAY_RAY_DOMAIN_LIST_SIZE);
#endif
      count += stats_[idx];
    }
    return count;
  }

 protected:
  void evaluateTraversalOrderAcrossCluster() {
    // #ifdef DEBUG_RAY_STATS
    //   fileoutQSizes("before_evalscore");
    //   fileoutStats("before_evalscore");
    // #endif

    // evaluate traversal score of all domains
    evalScores();

    // add up all ray counts at all depth levels
    // to get aggregated ray count of each domain
    updateRayCounts();

    // add up all scores across the cluster
    // sort them in descending order
    reduceScoresAndSort();

    // update traversal buffer by reading sorted scores in order
    updateSchedule();
  }

  /**
   * Broadcast cluster-level queue sizes and traversal domain IDs.
   *  First n elements are queue sizes.
   *  Next n elements are domain IDs.
   *  n is number of domains.
   */
  void bcastRayCountsAndSchedule() {
  //
#ifdef SPRAY_GLOG_CHECK
    CHECK_NOTNULL(counts_);
    CHECK_NOTNULL(schedule_);
    CHECK_NOTNULL(counts_and_schedule_);
#endif
    int count = ndomains_ << 1;
    MPI_Bcast(counts_and_schedule_, count, MPI_UINT32_T, 0, MPI_COMM_WORLD);

    int destination = 0;
    int num_ranks = mpi::size();
    for (int i = 0; i < ndomains_; ++i) {
      int id = schedule_[i];
#ifdef SPRAY_GLOG_CHECK
      CHECK_LT(id, ndomains_);
      CHECK_GE(id, 0);
#endif
      if (counts_[id]) {
        id_to_rank_map_[id] = destination;
        ++destination;
        if (destination == num_ranks) destination = 0;
      } else {
        // TODO: change to INT_MAX
        id_to_rank_map_[id] = -1;
      }
    }
  }

 protected:
  /**
   * Aggregate unsorted scores and queue sizes from the cluster.
   * Then sort the aggregated scores in descending order.
   */
  void reduceScoresAndSort() {
#ifdef SPRAY_GLOG_CHECK
    CHECK_NOTNULL(scores_and_counts_);
#endif
    if (mpi::isMultiProcess()) {
      int count = ndomains_ << 1;
      if (mpi::rank() == 0) {
        MPI_Reduce(MPI_IN_PLACE, scores_and_counts_, count, MPI_UINT32_T,
                   MPI_SUM, 0, MPI_COMM_WORLD);
      } else {
        MPI_Reduce(scores_and_counts_, nullptr, count, MPI_UINT32_T, MPI_SUM, 0,
                   MPI_COMM_WORLD);
      }
    }
    // sort
    sortScoresInDescendingOrder();
  }

  /**
   * Copy Domain IDs from RayStats::sorted_scores_ to
   * RayStats::schedule_. Update with invalid domain ID if the domain has
   * no ray.
   */
  void updateSchedule() {
#ifdef SPRAY_GLOG_CHECK
    CHECK_NOTNULL(schedule_);
#endif
    for (int i = 0; i < ndomains_; ++i) {
      int id = sorted_scores_[i].id;
      schedule_[i] = id;
#ifdef SPRAY_GLOG_CHECK
      CHECK_LT(id, ndomains_);
      CHECK_GE(id, 0);
#endif
    }
  }

  /** Extract queue sizes from collected statistics for all domains. */
  void updateRayCounts() {
#ifdef SPRAY_GLOG_CHECK
    CHECK_NOTNULL(counts_);
#endif
    for (int id = 0; id < ndomains_; ++id) {
      counts_[id] = getRayCountFromStats(id);
    }
  }

 protected:
  /**
   * Evaluates the base index of RayStats::stats_ for a given domain.
   * \param [in] id Domain ID.
   */
  int statsOffset(int id) const {
    int offset = id * SPRAY_RAY_DOMAIN_LIST_SIZE;
#ifdef SPRAY_GLOG_CHECK
    CHECK_LT(offset, ndomains_ * SPRAY_RAY_DOMAIN_LIST_SIZE);
#endif
    return offset;
  }

  /**
   * Evaluates index of RayStats::stats_.
   * \param [in] id Domain ID.
   * \param [in] depth Domain depth.
   */
  int statsIndex(int id, int depth) const {
    int index = id * SPRAY_RAY_DOMAIN_LIST_SIZE + depth;
#ifdef SPRAY_GLOG_CHECK
    CHECK_LT(index, ndomains_ * SPRAY_RAY_DOMAIN_LIST_SIZE);
#endif
    return index;
  }

 protected:
  uint32_t getStats(int id, int depth) const {
    int idx = statsIndex(id, depth);
#ifdef SPRAY_GLOG_CHECK
    CHECK_GE(id, 0);
    CHECK_LT(id, ndomains_);
    CHECK_LT(depth, SPRAY_RAY_DOMAIN_LIST_SIZE);
    CHECK_NOTNULL(stats_);
    CHECK_LT(idx, ndomains_ * SPRAY_RAY_DOMAIN_LIST_SIZE);
#endif
    return stats_[idx];
  }

 protected:
  /** Evaluate scores based on statistics and save them in RayStats::scores_. */
  void evalScores() {
#ifdef SPRAY_GLOG_CHECK
    CHECK_NOTNULL(scores_);
#endif
    for (int id = 0; id < ndomains_; ++id) {
      uint32_t score = 0;
      for (int depth = 0; depth < SPRAY_RAY_DOMAIN_LIST_SIZE; ++depth) {
        int w = SPRAY_RAY_DOMAIN_LIST_SIZE - depth;
        uint32_t c = getStats(id, depth);
        score += (c * w);
      }
      scores_[id] = score;
    }
  }

 protected:
  void sortScoresInDescendingOrder() {
#ifdef SPRAY_GLOG_CHECK
    CHECK_NOTNULL(scores_);
    CHECK(sorted_scores_.size());
#endif
    // populate scores with IDs
    for (int i = 0; i < ndomains_; ++i) {
      sorted_scores_[i].id = i;
      sorted_scores_[i].score = scores_[i];
    }

    // WARNING: seg faults if begin()/end() used inside openmp parallel region
    // e.g. the following seg faults
    // std::sort(unsorted_scores_.begin(), unsorted_scores_.end(),
    //           [](const Score & a, const Score & b)
    //               ->bool { return a.score > b.score; });

    // sort scores in descending order
    Score* s = &sorted_scores_[0];
    std::sort(s, s + ndomains_, [](const Score& a, const Score& b) -> bool {
      return a.score > b.score;
    });
  }

  // public:
  //  void printScores() const;
  //  void printIds() const;

 protected:
  /** Number of domains in the scene. */
  int ndomains_;

  /**
   * Dynamically allocated multi-purpose buffer.
   * - N: number of domains.
   * - 1st N: array of traversal scores.
   * - 2nd N: array of ray counts.
   * - 3rd N: array of schedlue.
   * - Rest: array of statistics.
   */
  uint32_t* buf_;

  /**
   * Array of per-domain ray statistics.
   * \sa #resize
   * - Pointer to RayStats:buf_[0] or RayStats:buf_[3N].
   * - N: number of domains.
   * - Size: N * SPRAY_RAY_DOMAIN_LIST_SIZE
   * - 1D array as 2D (1D: domain depth, 2D: domain).
   * - Ray count at each doman depth.
   */
  uint32_t* stats_;

  /**
   * Array of traversal scores.
   * - Pointer to RayStats:buf_[2N] (N: number of domains)
   * - Size: N
   */
  uint32_t* scores_;

  /**
   * Array of ray counts.
   * - Pointer to RayStats:buf_[N] (N: number of domains)
   * - Size: N
   */
  uint32_t* counts_;

  /**
   * Array of traversal order.
   * - Pointer to RayStats:buf_[2N] (N: number of domains)
   * - Size: N
   */
  uint32_t* schedule_;

  /**
   * Array of ray counts and scores.
   * - Pointer to RayStats:buf_[0]
   * - Size: 2N (N: number of domains)
   */
  uint32_t* scores_and_counts_;

  /**
   * Array of ray counts and scores.
   * - Pointer to RayStats:buf_[N] (N: number of domains)
   * - Size: 2N
   */
  uint32_t* counts_and_schedule_;

  /**
   * Array of sorted scores.
   * - Size: N (number of domains)
   */
  std::vector<Score> sorted_scores_;

  /**
   * Array of domain to rank mapping.
   * - Size: number of domains.
   * - Non-empty domains distributed in round robin based on the traversal
   * order.
   */
  std::vector<int> id_to_rank_map_;
};

}  // namespace baseline
}  // namespace spray

