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

#include <cstdint>
#include <fstream>
#include <iostream>
#include <list>
#include <string>
#include <vector>

#include "glog/logging.h"

#include "renderers/spray.h"

namespace spray {

class PlyLoader {
 public:
  struct Header {
    std::size_t num_vertices;  // out
    std::size_t num_faces;     // out
    bool has_color;            // out
  };

  struct Data {
    std::size_t vertices_capacity;  // in
    std::size_t faces_capacity;     // in
    std::size_t colors_capacity;    // in

    float *vertices;           // in/out
    std::size_t num_vertices;  // out

    uint32_t *faces;        // in/out
    std::size_t num_faces;  // out

    uint32_t *colors;  // rgb, in/out
  };

 public:
  enum DataFormat { kUNDEFINED_FORMAT, kASCII, kBIG_ENDIAN, kLITTLE_ENDIAN };
  enum ElementName { kUNDEFINED_ELEMENT_NAME, kVERTEX, kFACE };
  enum PropertyName {
    kUNDEFINED_PROPERTY,
    kX,
    kY,
    kZ,
    kRED,
    kGREEN,
    kBLUE,
    kVERTEX_INDICES
  };

  enum DataType {
    kUNDEFINED_DATA_TYPE,
    kCHAR,    // 1 byte(s)
    kUCHAR,   // 1
    kSHORT,   // 2
    kUSHORT,  // 2
    kINT,     // 4
    kUINT,    // 4
    kFLOAT,   // 4
    kDOUBLE   // 8
  };

 public:
  struct Element {
    int name;                  // vertex or face
    std::size_t num_elements;  // num of vertices or number of faces

    // vertex specific
    bool has_color;
    int vertex_dtype;
    int color_dtype;
    int property;  // temporary dummy variable to check for correct ordering

    // face specific
    int list_size_dtype;
    int index_dtype;
  };

 public:
  PlyLoader();

  // for test purposes
  static void quickHeaderRead(const std::string &filename, Header *header);

  void load(const std::string &filename, Data *d);

 private:
  void parseHeader();

 public:
  const std::vector<Element> &getElements() const { return elements_; }
  int getFormat() const { return format_; }
  bool isFileOpen() const { return file_.is_open(); }

 private:
  void parseElement(std::size_t elem_idx, std::istringstream &ss);
  void parseProperty(std::size_t elem_idx, std::istringstream &ss);

  int getDataType(const std::string &s) const;
  int getFormat(const std::string &s) const;

  void parseVertices(const Element &e, Data *d);
  void parseFaces(const Element &e, Data *d);

 private:
  int getDataSize(int type);

  void getLittleEndian(void *value, std::size_t bytes) {
    file_.read((char *)value, bytes);
  }

  void getBigEndian(void *value, std::size_t bytes) {
    for (std::size_t i = 0; i < bytes; ++i) {
      char *v = (char *)value;
      file_.read(((char *)(v + bytes - i - 1)), 1);
    }
  }

 private:
  std::ifstream file_;

  int format_;
  std::size_t num_vertices_;
  std::size_t num_faces_;

  std::vector<Element> elements_;
};

}  // namespace spray

