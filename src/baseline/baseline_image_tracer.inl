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

#if !defined(SPRAY_BASELINE_IMAGE_TRACER_INL)
#error An implementation of ImageTracer
#endif

namespace spray {
namespace baseline {

template <typename CacheT, typename ScheduleT, typename ShaderT>
void ImageTracer<CacheT, ScheduleT, ShaderT>::trace() {
  Base::image_->clear();

  Tile tile = Base::img_sched_.schedule();

  DRay *eyeray_buf = nullptr;
  std::size_t num_eyerays = 0;

  if (tile.isValid()) {
    num_eyerays = tile.getArea() * Base::pixel_samples_;
    std::size_t eyeray_buf_size = num_eyerays;
    eyeray_buf =
        Base::eyeray_arena_.template Alloc<DRay>(eyeray_buf_size, false);
    CHECK_NOTNULL(eyeray_buf);
  }

  int num_bounces = Base::bounces_;
  int nsamples = Base::pixel_samples_;

#pragma omp parallel firstprivate(tile, eyeray_buf, num_eyerays, num_bounces, \
                                  nsamples)
  {
    int ndomains = Base::ndomains_;
    int tid = omp_get_thread_num();
    int nthreads = omp_get_num_threads();
#ifdef SPRAY_GLOG_CHECK
    CHECK_GT(nthreads, 0);
    CHECK_LT(tid, nthreads);
#endif

#pragma omp master
    Base::scan_.resize(nthreads);
#pragma omp barrier

    ArenaQs<DRayQItem> qs;
    ArenaQs<DRayQItem> sqs;
    qs.resize(ndomains);
    sqs.resize(ndomains);

    DRayQ temp_q;

    MemoryArena arena;

    CommitBufferB retire_buf;

    QStats stats;
    stats.init(ndomains);

    DomainIntersector<CacheT> domain_isector(ndomains, Base::scene_);

    if (num_eyerays) {
      Base::genEyeRays(ndomains, nsamples, tile, eyeray_buf);
      Base::isectEyeDomains(ndomains, num_eyerays, eyeray_buf, &qs);
    }
#pragma omp barrier

    RTCRay rtc_ray;
    RTCRayIntersection isect;

    while (true) {
#pragma omp barrier
      schedule(ndomains, qs, sqs, &stats);
#pragma omp barrier
      int id;
      bool done = Base::ray_sched_.getSchedule(&id);
      if (done) break;

      const std::vector<RayCount> &sched = Base::ray_sched_.getSchedule();

#pragma omp barrier
      Base::spliceRecvQs(tid, id, ndomains, &qs, &(Base::recv_rbuf_));
#pragma omp barrier
      Base::spliceRecvQs(tid, id, ndomains, &sqs, &(Base::recv_sbuf_));
#pragma omp barrier

      if (id < ndomains) {
        qs.reset(id);
        sqs.reset(id);

        Base::processRays(tid, id, ndomains, num_bounces, &qs, &sqs, &arena,
                          &domain_isector, &rtc_ray, &isect, &retire_buf,
                          &temp_q);
      }
#pragma omp barrier

#pragma omp critical
      { retire_buf.retire(Base::image_, nsamples); }
      retire_buf.reset();
#pragma omp barrier
#pragma omp master
      Base::comm_.resetQitemRecvArena();
#pragma omp barrier
    }  // while (true)

    arena.Reset();
#pragma omp master
    {
      Base::comm_.resetMsgRecvArena();
      Base::eyeray_arena_.Reset();
    }
#pragma omp barrier
  }  // end of omp parallel
}

template <typename CacheT, typename ScheduleT, typename ShaderT>
void ImageTracer<CacheT, ScheduleT, ShaderT>::schedule(
    int ndomains, const ArenaQs<DRayQItem> &qs, const ArenaQs<DRayQItem> &sqs,
    QStats *stats) {
  for (int i = 0; i < ndomains; ++i) {
    int64_t count = qs.size(i) + sqs.size(i);
    stats->set(i, count);
  }
#pragma omp master
  Base::qstats_.reset();
#pragma omp barrier
  Base::aggregateStats(ndomains, *stats, &(Base::qstats_));
#pragma omp barrier

#pragma omp master
  {
    Base::ray_sched_.schedule(Base::qstats_);
    int nranks = mpi::size();
    Base::recv_rbuf_.reset(0);
    Base::recv_sbuf_.reset(0);
  }
}

}  // namespace baseline
}  // namespace spray
