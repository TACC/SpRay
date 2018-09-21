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

#include <mpi.h>
#include <cstddef>

#include "icet/src/include/IceT.h"
#include "icet/src/include/IceTMPI.h"

namespace spray {

class ImgCompGatherByte {
 public:
  void initialize(void* buf, MPI_Comm comm = MPI_COMM_WORLD);
  void* run(int image_w, int image_h);

 private:
  unsigned char* buf_;  // expects rgba, externally managed
  MPI_Comm mpi_comm_;
};

class ImgCompGatherFloat {
 public:
  void initialize(void* buf, MPI_Comm comm = MPI_COMM_WORLD);
  void* run(int image_w, int image_h);

 private:
  float* buf_;  // expects rgba, externally managed
  MPI_Comm mpi_comm_;
};

struct IceTConfig {
  IceTEnum mode = ICET_COMPOSITE_MODE_BLEND;
  // IceTEnum strategy = ICET_STRATEGY_DIRECT;
  // IceTEnum strategy = ICET_STRATEGY_SEQUENTIAL; // best for single tile
  // IceTEnum strategy = ICET_STRATEGY_SPLIT;
  IceTEnum strategy = ICET_STRATEGY_REDUCE;  // all-around good performer
  // IceTEnum strategy = ICET_STRATEGY_VTREE;
  IceTEnum depth = ICET_IMAGE_DEPTH_NONE;
};

class ImgCompIceTFloat {
 public:
  void initialize(void* buf, MPI_Comm comm = MPI_COMM_WORLD,
                  const IceTConfig& cfg = IceTConfig());

  void* run(int image_w, int image_h);

 private:
  IceTCommunicator icet_comm_;
  IceTContext icet_context_;
  float* buf_;  // expects rgba, externally managed
  MPI_Comm mpi_comm_;
};

class ImgCompIceTByte {
 public:
  void initialize(void* buf, MPI_Comm comm = MPI_COMM_WORLD,
                  const IceTConfig& cfg = IceTConfig());

  void* run(int image_w, int image_h);

 private:
  IceTCommunicator icet_comm_;
  IceTContext icet_context_;
  unsigned char* buf_;  // expects rgba, externally managed
  MPI_Comm mpi_comm_;
};

class ImgCompBypass {
 public:
  void initialize(void* buf, MPI_Comm comm = MPI_COMM_WORLD,
                  const IceTConfig& cfg = IceTConfig()) {
    buf_ = buf;
  }

  void* run(int image_w, int image_h) { return buf_; }

 private:
  void* buf_;  // expects rgba, externally managed
  MPI_Comm mpi_comm_;
};

}  // namespace spray

