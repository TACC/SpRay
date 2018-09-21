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

#include "display/composite.h"

#include <cstddef>

#include "renderers/spray.h"
#include "utils/comm.h"
#include "utils/math.h"

namespace spray {

void ImgCompGatherByte::initialize(void* buf, MPI_Comm comm) {
  buf_ = (unsigned char*)buf;
  mpi_comm_ = comm;
}

void* ImgCompGatherByte::run(int image_w, int image_h) {
  int image_area = image_w * image_h;
  int send_count = image_area << 2;  // rgba

  unsigned char* recvbuf = nullptr;

  bool is_worker_root = (mpi::rank() == 0);

  if (is_worker_root) {
    recvbuf = static_cast<unsigned char*>(
        malloc(send_count * global_mpi_comm.size * sizeof(unsigned char)));
  }

  MPI_Gather(buf_, send_count, MPI_UNSIGNED_CHAR, recvbuf, send_count,
             MPI_UNSIGNED_CHAR, 0, mpi_comm_);

  if (is_worker_root) {
    for (int i = 0; i < image_area; ++i) {
      glm::ivec3 acc(0);
      int framebuf_pos = i << 2;
      for (int rank = 0; rank < global_mpi_comm.size; ++rank) {
        int recvbuf_pos = rank * send_count + framebuf_pos;
        acc[0] += recvbuf[recvbuf_pos];
        acc[1] += recvbuf[recvbuf_pos + 1];
        acc[2] += recvbuf[recvbuf_pos + 2];
      }
      buf_[framebuf_pos] =
          static_cast<unsigned char>(spray::clamp(acc[0], 0x0, 0xff));
      buf_[framebuf_pos + 1] =
          static_cast<unsigned char>(spray::clamp(acc[1], 0x0, 0xff));
      buf_[framebuf_pos + 2] =
          static_cast<unsigned char>(spray::clamp(acc[2], 0x0, 0xff));
      buf_[framebuf_pos + 3] = 0xff;
    }
    free(recvbuf);
  }

  return (void*)buf_;
}

void ImgCompGatherFloat::initialize(void* buf, MPI_Comm comm) {
  buf_ = (float*)buf;
  // size_ = size;
  mpi_comm_ = comm;
}

void* ImgCompGatherFloat::run(int image_w, int image_h) {
  int image_area = image_w * image_h;
  int send_count = image_area << 2;  // rgba

  float* recvbuf = nullptr;

  bool is_worker_root = (mpi::rank() == 0);

  if (is_worker_root) {
    recvbuf =
        static_cast<float*>(malloc(send_count * Comm::size() * sizeof(float)));
  }

  MPI_Gather(buf_, send_count, MPI_FLOAT, recvbuf, send_count, MPI_FLOAT, 0,
             mpi_comm_);

  if (is_worker_root) {
    for (int i = 0; i < image_area; ++i) {
      glm::vec3 acc(0.f);
      int framebuf_pos = i << 2;
      for (int rank = 0; rank < global_mpi_comm.size; ++rank) {
        int recvbuf_pos = rank * send_count + framebuf_pos;
        acc[0] += recvbuf[recvbuf_pos];
        acc[1] += recvbuf[recvbuf_pos + 1];
        acc[2] += recvbuf[recvbuf_pos + 2];
      }
      buf_[framebuf_pos] = static_cast<float>(spray::clamp(acc[0], 0.f, 1.f));
      buf_[framebuf_pos + 1] =
          static_cast<float>(spray::clamp(acc[1], 0.f, 1.f));
      buf_[framebuf_pos + 2] =
          static_cast<float>(spray::clamp(acc[2], 0.f, 1.f));
      buf_[framebuf_pos + 3] = 1.f;
    }
    free(recvbuf);
  }

  return (void*)buf_;
}

void ImgCompIceTFloat::initialize(void* buf, MPI_Comm comm,
                                  const IceTConfig& cfg) {
  mpi_comm_ = comm;

  buf_ = (float*)buf;

  icet_comm_ = icetCreateMPICommunicator(mpi_comm_);
  icet_context_ = icetCreateContext(icet_comm_);

#ifndef NDEBUG
  int rank, num_proc;
  icetGetIntegerv(ICET_RANK, &rank);
  icetGetIntegerv(ICET_NUM_PROCESSES, &num_proc);

  printf("[ICET INFO] rank %d, %d processes\n", rank, num_proc);
#endif

  icetCompositeMode(cfg.mode);
  icetStrategy(cfg.strategy);
  icetSetColorFormat(ICET_IMAGE_COLOR_RGBA_FLOAT);
  icetSetDepthFormat(cfg.depth);
}

void* ImgCompIceTFloat::run(int image_w, int image_h) {
  // make this always black, do not use this feature
  // if we do, black pixels (INCLUDING shadow regions) all become background
  // color
  const IceTFloat background_color[] = {0.0, 0.0, 0.0, 0.0};

  // icetGLDrawFrame();

  // layout of tiled display
  icetResetTiles();
  icetAddTile(0, 0, image_w, image_h, 0 /*root process*/);

  IceTImage image = icetCompositeImage(
      buf_, nullptr /* depth_buffer */, nullptr /* valid_pixels_viewport */,
      nullptr /* projection_matrix */, nullptr /* modelview_matrix */,
      background_color);

  IceTFloat* icet_buf = nullptr;

  bool is_worker_root = (mpi::rank() == 0);

  if (is_worker_root) {
    icet_buf = icetImageGetColorf(image);
  }
  return (void*)icet_buf;
}

void ImgCompIceTByte::initialize(void* buf, MPI_Comm comm,
                                 const IceTConfig& cfg) {
  mpi_comm_ = comm;

  buf_ = (unsigned char*)buf;

  icet_comm_ = icetCreateMPICommunicator(mpi_comm_);
  icet_context_ = icetCreateContext(icet_comm_);
  icetCompositeMode(cfg.mode);
  icetStrategy(cfg.strategy);
  icetSetColorFormat(ICET_IMAGE_COLOR_RGBA_UBYTE);
  icetSetDepthFormat(cfg.depth);
}

void* ImgCompIceTByte::run(int image_w, int image_h) {
  // make this always black, do not use this feature
  // if we do, black pixels (INCLUDING shadow regions) all become background
  // color
  const IceTFloat background_color[] = {0.0, 0.0, 0.0, 0.0};

  // icetGLDrawFrame();
  bool is_worker_root = (mpi::rank() == 0);

  // layout of tiled display
  icetResetTiles();
  icetAddTile(0, 0, image_w, image_h, 0 /* root process */);

  IceTImage image = icetCompositeImage(
      buf_, nullptr /* depth_buffer */, nullptr /* valid_pixels_viewport */,
      nullptr /* projection_matrix */, nullptr /* modelview_matrix */,
      background_color);

  IceTUByte* icet_buf = nullptr;

  if (is_worker_root) {
    icet_buf = icetImageGetColorub(image);
  }
  return (void*)icet_buf;
}

}  // namespace spray

