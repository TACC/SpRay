#if !defined(DRT_TRACER_BASELINE_INSITU_SYNC_INL)
#error An implementation of TracerBaselineInsitu
#endif
// DEPRECATED !!!!!!!!!!!!!!!!!!!!!!!!!!
#error deprecated

namespace drt {

//////////////////////////////////////////////////////////////////////////

// insitu


//////////////////////////////////////////////////////////////////////////

template <typename CacheT, typename ScheduleT>
void TracerBaseline<CacheT, ScheduleT, InsituPartition>::schedule(
    int ndomains, const ArenaQs<DRayQItem> &qs, const ArenaQs<DRayQItem> &sqs,
    QStats *stats) {
  // place barrier outside function
  // #pragma omp barrier

  for (int i = 0; i < ndomains; ++i) {
    int64_t count = qs.size(i) + sqs.size(i);
    // int64_t count = qs[i].size<DRayQItem>() + sqs[i].size<DRayQItem>();
    stats->set(i, count);
  }
#pragma omp master
  Base::qstats_.reset();
#pragma omp barrier
  Base::aggregateStats(ndomains, *stats, &qstats_);
#pragma omp barrier

#pragma omp master
  {
    Base::ray_sched_.schedule(qstats_);

    // reset send/recv queues
    for (int i = 0; i < ndomains; ++i) {
      Base::send_qs_.reset(i);
      Base::send_sqs_.reset(i);
    }
    Base::recv_q_.reset();
    Base::recv_sq_.reset();
  }
  // place barrier outside function
  // #pragma omp barrier
}

//////////////////////////////////////////////////////////////////////////

}  // namespace drt
