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

#include <chrono>

namespace spray {

class Timer {
 public:
  Timer() { reset(); }

  void reset() {
    total_time_ = 0.0;
    tstart_ = std::chrono::high_resolution_clock::now();
    tend_ = tstart_;
  }

  void start() { tstart_ = std::chrono::high_resolution_clock::now(); }

  void stop() {
    tend_ = std::chrono::high_resolution_clock::now();
    total_time_ += getElapsed();
  }

  /**
   * Get elapsed time in seconds.
   */
  double getElapsed() {
    std::chrono::duration<double> elapsed =
        std::chrono::duration_cast<std::chrono::duration<double>>(tend_ -
                                                                  tstart_);
    return elapsed.count();
  }

  double getTotalTime() const { return total_time_; }

 private:
  std::chrono::high_resolution_clock::time_point tstart_;
  std::chrono::high_resolution_clock::time_point tend_;
  double total_time_;
};

}  // namespace spray

