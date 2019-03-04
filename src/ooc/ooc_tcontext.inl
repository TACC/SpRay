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

#if !defined(SPRAY_OOC_TCONTEXT_INL)
#error An implementation of ooc::TContext
#endif

namespace spray {
namespace ooc {

template <typename SceneT, typename ShaderT>
void TContext<SceneT, ShaderT>::procRads(int id, SceneT* scene,
                                         SceneInfo& sinfo, ShaderT& shader,
                                         int ray_depth) {
  //
  while (!frq_.empty()) {
    Ray* r = frq_.front();
    frq_.pop();

    bool is_hit = scene->intersect(sinfo.rtc_scene, sinfo.cache_block, r->org,
                                   r->dir, &rtc_isect_);

    if (is_hit) {
      if (vbuf_.update(rtc_isect_.tfar, r)) {
        shader(id, *r, rtc_isect_, mem_out_, &sq2_, &rq2_, &pending_q_,
               ray_depth);
        procShads2(id, scene, sinfo);
        procRads2(scene, sinfo);
      }
    }
  }
}

template <typename SceneT, typename ShaderT>
void TContext<SceneT, ShaderT>::procShads(SceneT* scene, SceneInfo& sinfo,
                                          std::queue<Ray*>* qin,
                                          std::queue<Ray*>* qout) {
  while (!qin->empty()) {
    Ray* r = qin->front();
    qin->pop();
    bool is_occluded =
        scene->occluded(sinfo.rtc_scene, r->org, r->dir, &rtc_ray_);

    if (is_occluded) {
      r->occluded = 1;
    }

    if (!is_occluded && !r->committed) {
      r->committed = 1;
      qout->push(r);
    }
  }
}

template <typename SceneT, typename ShaderT>
void TContext<SceneT, ShaderT>::procShads2(int id, SceneT* scene,
                                           SceneInfo& sinfo) {
  while (!sq2_.empty()) {
    Ray* r = sq2_.front();
    sq2_.pop();
    bool is_occluded =
        scene->occluded(sinfo.rtc_scene, r->org, r->dir, &rtc_ray_);
    if (!is_occluded) {
      if (!isector_.intersect(id, scene, r, sqs_out_,
                              &rstats_)) {  // unoccluded
#ifdef SPRAY_GLOG_CHECK
        CHECK_EQ(r->occluded, 0);
#endif
        r->committed = 1;
        commit_q_->push(r);
      }
    }
  }
}

template <typename SceneT, typename ShaderT>
void TContext<SceneT, ShaderT>::procRads2(SceneT* scene, SceneInfo& sinfo) {
  while (!rq2_.empty()) {
    Ray* r = rq2_.front();
    rq2_.pop();
#ifdef SPRAY_BACKGROUND_COLOR
    isector_.intersect(num_domains_, scene, r, &rqs_, &rstats_, bg_commit_q_);
#else
    isector_.intersect(num_domains_, scene, r, &rqs_, &rstats_);
#endif
  }
}

template <typename SceneT, typename ShaderT>
void TContext<SceneT, ShaderT>::resize(int ndomains, int num_pixel_samples,
                                       const Tile& tile, spray::HdrImage* image,
                                       int num_bounces,
                                       const glm::vec3& bg_color) {
  // tid_ = tid;
  num_domains_ = ndomains;
  num_pixel_samples_ = num_pixel_samples;
  num_bounces_ = num_bounces;
  one_over_num_pixel_samples_ = 1.0 / (double)num_pixel_samples_;

  vbuf_.resize(tile, num_pixel_samples);
  image_ = image;

  rqs_.resize(ndomains);
  sqs_in_->resize(ndomains);
  sqs_out_->resize(ndomains);

  rstats_.resize(ndomains, true /*stats_only*/);
  bg_color_ = bg_color;
  isector_.init(ndomains);
}

template <typename SceneT, typename ShaderT>
void TContext<SceneT, ShaderT>::retire() {
  while (!retire_q_->empty()) {
    Ray* r = retire_q_->front();
    retire_q_->pop();
    if (!r->occluded) {
      if (vbuf_.correct(*r)) {
        image_->add(r->pixid, r->w, one_over_num_pixel_samples_);
      }
    }
  }
#ifdef SPRAY_BACKGROUND_COLOR
  retireBackground();
#endif
}

template <typename SceneT, typename ShaderT>
void TContext<SceneT, ShaderT>::retireBackground() {
  glm::vec3 bgcolor;
  while (!bg_retire_q_->empty()) {
    Ray* ray = bg_retire_q_->front();
    bg_retire_q_->pop();
    if (vbuf_.correctAndMiss(*ray)) {
      bgcolor = glm::vec3(ray->w[0], ray->w[1], ray->w[2]) *
                spray::computeBackGroundColor(ray->dir, bg_color_);
      image_->add(ray->pixid, &bgcolor[0], one_over_num_pixel_samples_);
    }
  }
}

template <typename SceneT, typename ShaderT>
void TContext<SceneT, ShaderT>::filterRqs(int id) {
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

template <typename SceneT, typename ShaderT>
void TContext<SceneT, ShaderT>::filterSqs(int id, QVector<RayData>* sqs,
                                          std::queue<Ray*>* fsq) {
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
