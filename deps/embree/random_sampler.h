// ======================================================================== //
// Copyright 2009-2017 Intel Corporation                                    //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

// updates for spray:
// - changed namespace name 
// - commented out RandomSampler_get3D function
// - commented out inclusion of vec.h

#pragma once

// #include "../math/vec.h"
#include <cstdint>
#include <glm/glm.hpp>

namespace spray {

struct RandomSampler {
  uint32_t s;
};

__forceinline uint32_t MurmurHash3_mix(uint32_t hash, uint32_t k) {
  const uint32_t c1 = 0xcc9e2d51;
  const uint32_t c2 = 0x1b873593;
  const uint32_t r1 = 15;
  const uint32_t r2 = 13;
  const uint32_t m = 5;
  const uint32_t n = 0xe6546b64;

  k *= c1;
  k = (k << r1) | (k >> (32 - r1));
  k *= c2;

  hash ^= k;
  hash = ((hash << r2) | (hash >> (32 - r2))) * m + n;

  return hash;
}

__forceinline uint32_t MurmurHash3_finalize(uint32_t hash) {
  hash ^= hash >> 16;
  hash *= 0x85ebca6b;
  hash ^= hash >> 13;
  hash *= 0xc2b2ae35;
  hash ^= hash >> 16;

  return hash;
}

__forceinline uint32_t LCG_next(uint32_t value) {
  const uint32_t m = 1664525;
  const uint32_t n = 1013904223;

  return value * m + n;
}

__forceinline void RandomSampler_init(RandomSampler& self, int id) {
  uint32_t hash = 0;
  hash = MurmurHash3_mix(hash, id);
  hash = MurmurHash3_finalize(hash);

  self.s = hash;
}

__forceinline void RandomSampler_init(RandomSampler& self, int pixelId,
                                      int sampleId) {
  uint32_t hash = 0;
  hash = MurmurHash3_mix(hash, pixelId);
  hash = MurmurHash3_mix(hash, sampleId);
  hash = MurmurHash3_finalize(hash);

  self.s = hash;
}

__forceinline void RandomSampler_init(RandomSampler& self, int x, int y,
                                      int sampleId) {
  RandomSampler_init(self, x | (y << 16), sampleId);
}

__forceinline int RandomSampler_getInt(RandomSampler& self) {
  self.s = LCG_next(self.s);
  return self.s >> 1;
}

__forceinline uint32_t RandomSampler_getUInt(RandomSampler& self) {
  self.s = LCG_next(self.s);
  return self.s;
}

__forceinline float RandomSampler_getFloat(RandomSampler& self) {
  return (float)RandomSampler_getInt(self) * 4.656612873077392578125e-10f;
}

__forceinline float RandomSampler_get1D(RandomSampler& self) {
  return RandomSampler_getFloat(self);
}

__forceinline glm::vec2 RandomSampler_get2D(RandomSampler& self) {
  const float u = RandomSampler_get1D(self);
  const float v = RandomSampler_get1D(self);
  return glm::vec2(u, v);
}

// __forceinline glm::vec3 RandomSampler_get3D(RandomSampler& self) {
//   /*
//   const float u = RandomSampler_get1D(self);
//   const float v = RandomSampler_get1D(self);
//   const float w = RandomSampler_get1D(self);
//   return glm::vec3(u,v,w);
//   */
//   const int u = RandomSampler_getUInt(self);
//   const int v = RandomSampler_getUInt(self);
//   const int w = RandomSampler_getUInt(self);
//   return glm::vec3(srl(Vec3ia(u, v, w), 1)) * 4.656612873077392578125e-10f;
// }

}  // namespace spray
