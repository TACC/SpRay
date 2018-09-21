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

#include "ooc/ooc_pcontext.h"

#include "utils/scan.h"

namespace spray {
namespace ooc {

void PContext::resize(int ndomains, int max_num_bounces, int num_threads,
                      const Tile& tile, int num_pixel_samples, int num_lights,
                      spray::HdrImage* image) {
  num_domains_ = ndomains;
  max_num_bounces_ = max_num_bounces;
  num_pixel_samples_ = num_pixel_samples;
  image_ = image;

  tile_.x = tile.x;
  tile_.y = tile.y;
  tile_.w = tile.w;
  tile_.h = tile.h;

  rstats_.resize(ndomains, false /*stats_only*/);

  scans_.resize(4);
  for (auto& s : scans_) {
    s.resize(num_threads);
  }

  domain_locks_.resize(ndomains);
  for (int i = 0; i < ndomains; ++i) {
    omp_init_lock(&domain_locks_[i]);
  }
}

void PContext::mergeStats(const DomainStats& src_stats,
                          DomainStats* dest_stats) {
  //
  for (int id = 0; id < num_domains_; ++id) {
    omp_set_lock(&domain_locks_[id]);
    dest_stats->addStats(id, src_stats);
    omp_unset_lock(&domain_locks_[id]);
  }
}

}  // namespace ooc
}  // namespace spray
