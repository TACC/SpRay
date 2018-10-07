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

#if !defined(SPRAY_BASELINE_COMM_INL)
#error An implementation of Comm
#endif

namespace spray {
namespace baseline {

// bit 0: shadow (1) or radiance(0)
// bit 31-1: domain id
struct MpiTagUtil {
  static int encode(bool shadow, int domid) {
    int tag = (int)shadow;
    tag = tag | (domid << 1);
    return tag;
  }
  static bool isShadow(int tag) { return tag & 1; }
  static int domainId(int tag) { return tag >> 1; }
};

template <class QItemT, class MessageT, class OutgoingCopierT,
          class IncomingCopierT>
void Comm<QItemT, MessageT, OutgoingCopierT, IncomingCopierT>::init(
    int ndomains, const std::vector<RayCount>& sched, BlockBuffer* send_rbuf,
    BlockBuffer* send_sbuf, BlockBuffer* recv_rbuf, BlockBuffer* recv_sbuf) {
  ndomains_ = ndomains;

  schedule_ = &sched;

  send_rbuf_ = send_rbuf;
  send_sbuf_ = send_sbuf;

  recv_rbuf_ = recv_rbuf;
  recv_sbuf_ = recv_sbuf;
}

template <class QItemT, class MessageT, class OutgoingCopierT,
          class IncomingCopierT>
void Comm<QItemT, MessageT, OutgoingCopierT, IncomingCopierT>::run() {
#ifdef SPRAY_GLOG_CHECK
  CHECK_EQ(isend_q_.size(), 0);
#endif
  // evaluate number of items to recv
  setup();

  std::size_t num_msgs_recved = 0, num_msgs_sent = 0;
  bool all_recved = (num_msgs_recved == incoming_msgbuf_.capacity);
  bool all_sent = (num_msgs_sent == outgoing_msgbuf_.capacity);

#ifdef SPRAY_GLOG_CHECK
  LOG(INFO) << RANK_THREAD << "num_msgs_recved: " << num_msgs_recved
            << " num_msgs2recv: " << incoming_msgbuf_.capacity
            << " num_msgs_sent: " << num_msgs_sent << " num_msgs2send "
            << outgoing_msgbuf_.capacity;
#endif

  int num_ranks = mpi::size();
  int my_rank = mpi::rank();
  int rank = 0;

  const std::vector<RayCount>& sched = *schedule_;

#ifdef SPRAY_TIMING
  spray::tStart(spray::TIMER_SYNC_RAYS);
#endif
  while (!(all_sent && all_recved)) {
    // recv
    if (!all_recved) {
      std::size_t num_msgs = recv();
      if (num_msgs) {
        num_msgs_recved += num_msgs;
        all_recved = (num_msgs_recved == incoming_msgbuf_.capacity);

#ifdef SPRAY_GLOG_CHECK
        CHECK_LE(num_msgs_recved, incoming_msgbuf_.capacity);
        LOG(INFO) << RANK_THREAD << "num_msgs_recved: " << num_msgs_recved
                  << " num_msgs2recv: " << incoming_msgbuf_.capacity;
#endif
      }
    }

    // send (rank by rank)
    if (!all_sent) {
      if (rank == my_rank) ++rank;
#ifdef SPRAY_GLOG_CHECK
      CHECK_LT(rank, num_ranks);
      CHECK_GE(rank, 0);
#endif
      int id = sched[rank].id;
#ifdef SPRAY_GLOG_CHECK
      CHECK_GE(id, 0);
      CHECK_LT(rank, num_ranks);
#endif
      if (id < ndomains_) {
        std::size_t num_msgs = send(id, rank);
        if (num_msgs) {
          num_msgs_sent += num_msgs;
          all_sent = (num_msgs_sent == outgoing_msgbuf_.capacity);

#ifdef SPRAY_GLOG_CHECK
          CHECK_LE(num_msgs_sent, outgoing_msgbuf_.capacity);
          LOG(INFO) << RANK_THREAD << "num_msgs_sent: " << num_msgs_sent
                    << " num_msgs2send: " << outgoing_msgbuf_.capacity;
#endif
        }
      }
      ++rank;
    }
  }
#ifdef SPRAY_TIMING
  spray::tStop(spray::TIMER_SYNC_RAYS);
#endif
}

template <class QItemT, class MessageT, class OutgoingCopierT,
          class IncomingCopierT>
std::size_t Comm<QItemT, MessageT, OutgoingCopierT, IncomingCopierT>::recv() {
  int flag = 0;
  MPI_Status status;

  MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);

  std::size_t num_msgs = 0;

  if (flag) {
    int count;
    // count number of rays to recv
    MPI_Get_count(&status, MPI_UINT8_T, &count);
    std::size_t msg_bytes = sizeof(uint8_t) * count;

    MessageT* recvbuf = incoming_msgbuf_.base + incoming_msgbuf_.size;
    num_msgs = util::getNumOfItems<MessageT>(recvbuf, msg_bytes);

#ifdef SPRAY_GLOG_CHECK
    CHECK_GT(count, 0);
    CHECK_EQ(count % sizeof(MessageT), 0);
    CHECK_LE(incoming_msgbuf_.size + num_msgs, incoming_msgbuf_.capacity)
        << "incoming_msgbuf_.size " << incoming_msgbuf_.size << " num_msgs "
        << num_msgs;
#endif

    int tag = status.MPI_TAG;

    MPI_Recv((void*)recvbuf, count, MPI_UINT8_T, status.MPI_SOURCE, tag,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    incoming_msgbuf_.size += num_msgs;

    // decode tag (domain id)
    int id = MpiTagUtil::domainId(tag);
    bool is_shadow = MpiTagUtil::isShadow(tag);

#ifdef SPRAY_GLOG_CHECK
    CHECK_EQ(id, (*schedule_)[mpi::rank()].id);
    CHECK_LE(incoming_msgbuf_.size, incoming_msgbuf_.capacity);
    CHECK_LT(id, ndomains_);
#endif

    BufferInfo<QItemT>* incoming_buffer;
    BlockBuffer* recv_blocks;
    if (is_shadow) {
      incoming_buffer = &incoming_shadow_qbuf_;
      recv_blocks = recv_sbuf_;
    } else {
      incoming_buffer = &incoming_qbuf_;
      recv_blocks = recv_rbuf_;
    }

    // populate queue buffer with queue items
    // copy message_buf to blk.buf
    BufferInfo<QItemT>& buf_info = *incoming_buffer;
#ifdef SPRAY_GLOG_CHECK
    CHECK_NOTNULL(buf_info.base);
    CHECK_LE(buf_info.size, buf_info.capacity);
#endif
    recverFunctor_(buf_info.base + buf_info.size /*qitem_buf*/,
                   recvbuf /*ray_buf*/, num_msgs);

    buf_info.size += num_msgs;
    recv_blocks->incrementLastBlockSize<QItemT>(0, num_msgs);

#ifdef SPRAY_GLOG_CHECK
    CHECK_LE(buf_info.size, buf_info.capacity);
#endif
  }
  return num_msgs;
}

template <class QItemT, class MessageT, class OutgoingCopierT,
          class IncomingCopierT>
std::size_t Comm<QItemT, MessageT, OutgoingCopierT, IncomingCopierT>::send(
    int id, int dest) {
#ifdef SPRAY_GLOG_CHECK
  CHECK_NE(dest, mpi::rank());
  CHECK_LE(id, ndomains_);
#endif
  auto num_items =
      send_rbuf_->size<QItemT>(dest) + send_sbuf_->size<QItemT>(dest);
  if (num_items == 0) return 0;

#ifdef SPRAY_PROFILE_COUNTERS
  spray::tAgg(spray::COUNTER_RAYS_SENT, num_items);
#endif

#ifdef SPRAY_GLOG_CHECK
  CHECK_LE(outgoing_msgbuf_.size, outgoing_msgbuf_.capacity);
#endif
  MPI_Request isend_rqst;

  std::size_t total_num_msgs = 0;

  MessageT* sendbuf = outgoing_msgbuf_.base + outgoing_msgbuf_.size;
  MessageT* sendbuf_tmp = sendbuf;

  for (int i = 0; i < 2; ++i) {
    MessageT* mpi_sendbuf = sendbuf_tmp;
    const BlockBuffer* blocks;
    if (i == 0) {
      blocks = send_rbuf_;
    } else {
      blocks = send_sbuf_;
    }
#ifdef SPRAY_GLOG_CHECK
    int block_id = 0;
#endif
    std::size_t total_num_items = 0;
    int nblocks = blocks->getNumBlocks(dest);

    for (int j = 0; j < nblocks; ++j) {
      const MemBlock& blk = blocks->getBlock(dest, j);
      if (blk.size) {
#ifdef SPRAY_GLOG_CHECK
        CHECK_NOTNULL(blk.buf);
#endif
        // get the current block
        std::size_t num_items =
            (QItemT*)(blk.buf + blk.size) - (QItemT*)blk.buf;

#ifdef SPRAY_GLOG_CHECK
        if (i == 0) {
          LOG(INFO) << RANK_THREAD << "sending radiance block [" << block_id
                    << "] num_rays: " << num_items << " size: " << blk.size
                    << " blk.buf" << (std::size_t)blk.buf << " dest " << dest;
        } else {
          LOG(INFO) << RANK_THREAD << "sending shadow block [" << block_id
                    << "] num_rays(msgs): " << num_items
                    << " size(bytes): " << blk.size << " blk.buf"
                    << (std::size_t)blk.buf << " dest " << dest;
        }
        CHECK_GT(blk.size, 0);
        ++block_id;
#endif

        // populate send buffer with block items
        senderFunctor_(sendbuf_tmp /*ray_buf*/, blk.buf /*qitem_buf*/,
                       num_items);

        sendbuf_tmp += num_items;
        total_num_items += num_items;
      }
    }

    if (total_num_items) {
      total_num_msgs += total_num_items;

      // save the request object
      isend_q_.push(isend_rqst);
      MPI_Request* isend_rqst_ptr = &(isend_q_.back());

      bool is_shadow = (i == 1);
      int tag = MpiTagUtil::encode(is_shadow, id);

      // nonblocking send
      MPI_Isend(mpi_sendbuf, total_num_items * sizeof(MessageT), MPI_UINT8_T,
                dest, tag, MPI_COMM_WORLD, isend_rqst_ptr);
#ifdef SPRAY_GLOG_CHECK
      LOG(INFO) << RANK_THREAD << "[sent : " << total_num_items
                << " rays(msgs) to rank " << dest << "]";
#endif
    }
  }  // end of for (int i = 0; i < 2; ++i) {

  outgoing_msgbuf_.size += total_num_msgs;
#ifdef SPRAY_GLOG_CHECK
  CHECK_LE(outgoing_msgbuf_.size, outgoing_msgbuf_.capacity);
#endif

  return total_num_msgs;
}

template <class QItemT, class MessageT, class OutgoingCopierT,
          class IncomingCopierT>
void Comm<QItemT, MessageT, OutgoingCopierT, IncomingCopierT>::setup() {
#ifdef SPRAY_GLOG_CHECK
  std::size_t check_bytes = 0;
#endif

  // evaluate num of items to send/recv, allocate/set memory block for ray
  // queuing
  std::size_t num_items2recv = 0;
  int myrank = mpi::rank();
  int num_ranks = mpi::size();

  const std::vector<RayCount>& sched = *schedule_;

  // setup recv
  int myid = sched[myrank].id;

  if (myid >= ndomains_) {  // unmapped
#ifdef SPRAY_GLOG_CHECK
    CHECK_EQ(sched[myrank].count, 0);
#endif
    incoming_qbuf_.reset();
    incoming_shadow_qbuf_.reset();

  } else {  // mapped
#ifdef SPRAY_GLOG_CHECK
    CHECK_GE(myid, 0);
    CHECK_GT(sched[myrank].count, 0);
#endif
    auto cluster_items = sched[myrank].count;
    auto process_items =
        recv_rbuf_->size<QItemT>(0) + recv_sbuf_->size<QItemT>(0);

    // number of items to receive
    auto num_incoming_items = cluster_items - process_items;
    num_items2recv += num_incoming_items;

#ifdef SPRAY_GLOG_CHECK
    LOG(INFO) << RANK_THREAD << "[recv] mapped rank " << myrank << " cluster "
              << cluster_items << " already owned " << process_items
              << " recv_rbuf " << recv_rbuf_->size<QItemT>(0) << " recv_sbuf "
              << recv_sbuf_->size<QItemT>(0) << " num items to recv "
              << num_items2recv;
    CHECK_GE(cluster_items, process_items) << "domain " << myid;
    CHECK_EQ(send_rbuf_->size<QItemT>(myrank), 0);
    CHECK_EQ(send_sbuf_->size<QItemT>(myrank), 0);
#endif
    // TODO: exact resizing

    if (num_incoming_items) {  // other nodes having rays this domain
      // allocate and set memory block for ray queuing
      BufferInfo<QItemT>& buf_info = incoming_qbuf_;
      buf_info.base =
          qitem_recv_arena_.Alloc<QItemT>(num_incoming_items, false);
      buf_info.size = 0;
      buf_info.capacity = num_incoming_items;

      CHECK_NOTNULL(buf_info.base);

#ifdef SPRAY_GLOG_CHECK
      CHECK_GT(recv_rbuf_->getNumBlocks(0), 0);
      CHECK_EQ(recv_rbuf_->getLastBlock(0).size, 0);
      CHECK(recv_rbuf_->getLastBlock(0).buf == nullptr);
#endif
      // set size upon receving (since we don't know how many to receive)
      // set last block
      recv_rbuf_->setLastBlock(0 /*id*/, 0 /*size*/, (uint8_t*)buf_info.base);

      BufferInfo<QItemT>& shadow_buf_info = incoming_shadow_qbuf_;
      shadow_buf_info.base =
          qitem_recv_arena_.Alloc<QItemT>(num_incoming_items, false);
      shadow_buf_info.size = 0;
      shadow_buf_info.capacity = num_incoming_items;

      CHECK_NOTNULL(shadow_buf_info.base);

#ifdef SPRAY_GLOG_CHECK
      CHECK_GT(recv_sbuf_->getNumBlocks(0), 0);
#endif
      // set size upon receving (since we don't know how many to receive)
      // set last block
      recv_sbuf_->setLastBlock(0 /*id*/, 0 /*size*/,
                               (uint8_t*)shadow_buf_info.base);

    } else {  // this node owning all rays in this domain
      incoming_qbuf_.reset();
      incoming_shadow_qbuf_.reset();

#ifdef SPRAY_GLOG_CHECK
      CHECK_GT(recv_rbuf_->getNumBlocks(0), 0);
      CHECK_EQ(recv_rbuf_->getLastBlock(0).size, 0);
      CHECK(recv_rbuf_->getLastBlock(0).buf == nullptr);

      CHECK_GT(recv_sbuf_->getNumBlocks(0), 0);
      CHECK_EQ(recv_sbuf_->getLastBlock(0).size, 0);
      CHECK(recv_sbuf_->getLastBlock(0).buf == nullptr);
#endif
      recv_rbuf_->setLastBlock(0, 0, nullptr);
      recv_sbuf_->setLastBlock(0, 0, nullptr);
    }
  }

  // setup send

  std::size_t num_items2send = 0;

  for (int rank = 0; rank < num_ranks; ++rank) {
    int id = sched[rank].id;
#ifdef SPRAY_GLOG_CHECK
    CHECK_GE(id, 0);
    if (id >= ndomains_) {
      CHECK_EQ(sched[rank].count, 0);
    }
#endif
    if (rank != myrank && id < ndomains_) {  // not mine and mapped
#ifdef SPRAY_GLOG_CHECK
      LOG(INFO) << RANK_THREAD << " id " << id << " [send] mapped rank " << rank
                << " send_rbuf " << send_rbuf_->size<QItemT>(rank)
                << " send_sbuf " << send_sbuf_->size<QItemT>(rank);
#endif
      num_items2send +=
          (send_rbuf_->size<QItemT>(rank) + send_sbuf_->size<QItemT>(rank));
    }
  }

  if (num_items2recv) {
    // allocate incoming message buffer
    incoming_msgbuf_.base = recv_arena_.Alloc<MessageT>(num_items2recv, false);
    incoming_msgbuf_.size = 0;
    incoming_msgbuf_.capacity = num_items2recv;

    CHECK_NOTNULL(incoming_msgbuf_.base);
  } else {
    incoming_msgbuf_.reset();
  }

  if (num_items2send) {
    // allocate message buffer
    outgoing_msgbuf_.base = send_arena_.Alloc<MessageT>(num_items2send, false);
    outgoing_msgbuf_.size = 0;
    outgoing_msgbuf_.capacity = num_items2send;

    CHECK_NOTNULL(outgoing_msgbuf_.base);

#ifdef SPRAY_GLOG_CHECK
    std::size_t bytes2send = num_items2send * sizeof(MessageT);
    CHECK_EQ(check_bytes, bytes2send);
#endif
  } else {
    outgoing_msgbuf_.reset();
  }
}

template <class QItemT, class MessageT, class OutgoingCopierT,
          class IncomingCopierT>
void Comm<QItemT, MessageT, OutgoingCopierT, IncomingCopierT>::waitForSend() {
  while (!isend_q_.empty()) {
    MPI_Request* r = &(isend_q_.front());
    MPI_Wait(r, MPI_STATUS_IGNORE);
    isend_q_.pop();
  }
}

}  // namespace baseline
}  // namespace spray

