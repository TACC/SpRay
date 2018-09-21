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

#include "ooc/ooc_tcontext.h"

#include "display/image.h"
#include "utils/scan.h"

namespace spray {
namespace ooc {

void TContext::resize(int ndomains, int num_pixel_samples, const Tile& tile,
                      spray::HdrImage* image, int num_bounces) {
  // tid_ = tid;
  num_domains_ = ndomains;
  num_pixel_samples_ = num_pixel_samples;
  num_bounces_ = num_bounces;

  vbuf_.resize(tile, num_pixel_samples);
  image_ = image;

  rqs_.resize(ndomains);
  sqs_in_->resize(ndomains);
  sqs_out_->resize(ndomains);

  rstats_.resize(ndomains, true /*stats_only*/);
}

void TContext::retire() {
  double scale = 1.0 / (double)num_pixel_samples_;
  while (!retire_q_->empty()) {
    Ray* r = retire_q_->front();
    retire_q_->pop();
    if (!r->occluded) {
      if (vbuf_.correct(*r)) {
        image_->add(r->pixid, r->w, scale);
      }
    }
  }
}

void TContext::filterRqs(int id) {
  auto* rq = rqs_.getQ(id);

  while (!rq->empty()) {
    RayData& data = rq->front();
    Ray* r = data.ray;

    rstats_.decrement(id, data.dom_depth);

    if (data.tdom <= r->history[r->depth] && vbuf_.correct(*r)) {
      frq_.push(r);
    }
    rq->pop();
  }
}

void TContext::filterSqs(int id, QVector<RayData>* sqs, std::queue<Ray*>* fsq) {
  auto* sq = sqs->getQ(id);

  while (!sq->empty()) {
    RayData& data = sq->front();
    Ray* r = data.ray;

    rstats_.decrement(id, data.dom_depth);

    if (!r->occluded && vbuf_.correct(*r)) {
      fsq->push(r);
    }
    sq->pop();
  }
}

}  // namespace ooc
}  // namespace spray
