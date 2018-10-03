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

#include <cstdint>
#include <vector>

#include "glog/logging.h"

#include "renderers/spray.h"

namespace spray {

class InsituPartition;

namespace baseline {

class QStats;

struct SPRAY_ALIGN(16) RayCount {
  int id;
  int64_t count;
};

class LoadAnyOnceSched {
 public:
  LoadAnyOnceSched();

  void init(int ndomains, const InsituPartition& partition);

  // returns true if done
  void schedule(const QStats& stats);
  bool getSchedule(int* id) const {
    *id = id_;
    return done_;
  }

  const std::vector<RayCount>& getSchedule() const { return schedule_; }

 protected:
  void scheduleImage(const QStats& stats);
  void scheduleMulti(const QStats& stats);

 protected:
  bool done_;
  int id_;

  // all processes' states (root proc only)
  std::vector<int64_t> state_;

  // per-process schedules (rank-to-domain ID mapping)
  std::vector<RayCount> schedule_;

  std::vector<RayCount> ray_counts_;

  int ndomains_;
};

class LoadAnyOnceInsituSched : public LoadAnyOnceSched {
 public:
  typedef LoadAnyOnceSched Base;
  void init(int ndomains, const InsituPartition& partition) {
    Base::init(ndomains, partition);
    partition_ = &partition;
  }

  void schedule(const QStats& stats);

 private:
  void scheduleMulti(const QStats& stats);

  const InsituPartition* partition_;
};

class LoadAnyOnceImageSched : public LoadAnyOnceSched {
 public:
  typedef LoadAnyOnceSched Base;
  void init(int ndomains, const InsituPartition& partition) {
    Base::init(ndomains, partition);
  }

  void schedule(const QStats& stats) { scheduleImage(stats); }
};

}  // namespace baseline
}  // namespace spray

