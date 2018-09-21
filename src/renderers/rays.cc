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

#include "renderers/rays.h"

#include <iostream>

namespace spray {

#ifdef SPRAY_GLOG_CHECK
std::ostream& operator<<(std::ostream& os, const DRay& r) {
  os << "[org " << r.org[0] << " " << r.org[1] << " " << r.org[2] << "] [pixid "
     << r.pixid << "][dir " << r.dir[0] << " " << r.dir[1] << " " << r.dir[2]
     << "] [samid " << r.samid << "][depth " << r.depth << "][w " << r.w[0]
     << " " << r.w[1] << " " << r.w[2] << "] [t" << r.t << "] [u" << r.u
     << "] [v" << r.v << "][primID " << r.primID << "][flag" << r.flag
     << "] [domid " << r.domid << "][domain_pos " << r.domain_pos
     << "][next_tdom " << r.next_tdom << "]";
  return os;
}
#endif

}  // namespace spray
