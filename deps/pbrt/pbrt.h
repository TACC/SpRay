
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

// #if defined(_MSC_VER)
// #define NOMINMAX
// #pragma once
// #endif
//
// #ifndef PBRT_CORE_PBRT_H
// #define PBRT_CORE_PBRT_H

#pragma once

// // core/pbrt.h*
// // Global Include Files
// #include <algorithm>
// #include <cinttypes>
// #include <cmath>
// #include <iostream>
// #include <limits>
// #include <memory>
// #include <string>
// #include <type_traits>
// #include <vector>
// #include "pbrt/error.h"
// #ifdef PBRT_HAVE_MALLOC_H
// #include <malloc.h>  // for _alloca, memalign
// #endif
// #ifdef PBRT_HAVE_ALLOCA_H
// #include <alloca.h>
// #endif
// #include <assert.h>
// #include <glog/logging.h>
// #include <string.h>

// // Platform-specific definitions
// #if defined(_WIN32) || defined(_WIN64)
// #define PBRT_IS_WINDOWS
// #endif

// #if defined(_MSC_VER)
// #define PBRT_IS_MSVC
// #if _MSC_VER == 1800
// #define snprintf _snprintf
// #endif
// #endif

#ifndef PBRT_L1_CACHE_LINE_SIZE
#define PBRT_L1_CACHE_LINE_SIZE 64
#endif

// #include <stdint.h>
// #if defined(PBRT_IS_MSVC)
// #include <float.h>
// #include <intrin.h>
// #pragma warning(disable : 4305)  // double constant assigned to float
// #pragma warning(disable : 4244)  // int -> float conversion
// #pragma warning(disable : 4843)  // double -> float conversion
// #endif

// Global Macros
#define ALLOCA(TYPE, COUNT) (TYPE *)alloca((COUNT) * sizeof(TYPE))

namespace spray {

template <typename T>
inline PBRT_CONSTEXPR bool IsPowerOf2(T v) {
  return v && !(v & (v - 1));
}

}  // namespace spray

// #endif  // PBRT_CORE_PBRT_H
