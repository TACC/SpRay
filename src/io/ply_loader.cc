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

#include <cstddef>
#include <cstring>
#include <iostream>
#include <vector>

#include "glog/logging.h"
#include "pbrt/memory.h"

#include "partition/aabb.h"
#include "io/ply_loader.h"

#define DEBUG_PLY_LOADER
#undef DEBUG_PLY_LOADER

#define MAX_NUM_ELEMENTS 2

namespace spray {

PlyLoader::PlyLoader() { elements_.resize(MAX_NUM_ELEMENTS); }

void PlyLoader::quickHeaderRead(const std::string &filename, Header *header) {
  // open file
  std::ifstream file;
  file.open(filename.c_str(), std::fstream::in | std::fstream::binary);
  CHECK(file.is_open());

  // init
  std::size_t num_vertices = 0, num_faces = 0;
  std::string line, word;

  // verify
  std::getline(file, line);
  CHECK(line == "ply") << "unknown file type";

  // search for
  //   "element vertex N"
  //   "element face N"
  int element_name = kUNDEFINED_ELEMENT_NAME;
  int color_components = 0;

  while (!file.eof()) {
    // read a line
    std::getline(file, line);
    std::istringstream ss(line);

    // read the first word
    ss >> word;

    if (word == "element") {
      // read the second word
      ss >> word;

      // read the third word
      if (word == "vertex") {
        ss >> num_vertices;
        element_name = kVERTEX;
      } else if (word == "face") {
        ss >> num_faces;
        element_name = kFACE;
      } else {
        element_name = kUNDEFINED_ELEMENT_NAME;
      }
    } else if (word == "property") {
      if (element_name == kVERTEX) {
        std::string data_type, property_name;
        ss >> data_type >> property_name;
        if (property_name == "red" || property_name == "green" ||
            property_name == "blue") {
          ++color_components;
        }
      }
    } else if (word == "end_header") {
      break;
    }
  }

  CHECK_GT(num_vertices, 0);
  CHECK_GT(num_faces, 0);

  header->num_vertices = num_vertices;
  header->num_faces = num_faces;

  if (color_components == 3) {
    header->has_color = true;
  } else {
    header->has_color = false;
  }

  // close file
  file.close();
}

void PlyLoader::parseHeader() {
  // assume file is open

  // verify
  std::string line;
  std::getline(file_, line);
  CHECK(line == "ply") << "unknown file type";

  // initialize
  format_ = kUNDEFINED_FORMAT;  // ascii/little/big format
  num_vertices_ = 0;
  num_faces_ = 0;

  std::string word;
  std::size_t next_elem_idx = 0, current_elem_idx = 0;
  bool end_header = false;

  // parse header
  while (!file_.eof()) {
    std::getline(file_, line);
    std::istringstream ss(line);
    // std::cout << ss.str() << std::endl;

    // read the first word
    ss >> word;

    if (word == "format") {
      CHECK_EQ(end_header, false);
      // read the second word
      ss >> word;

      // update format
      format_ = getFormat(word);
      // TODO: yet to test
      CHECK_NE(format_, kBIG_ENDIAN);

    } else if (word == "element") {
      CHECK_EQ(end_header, false);
      CHECK_NE(format_, kUNDEFINED_FORMAT);
      CHECK_LT(next_elem_idx, elements_.size());

      current_elem_idx = next_elem_idx;

      parseElement(current_elem_idx, ss);
      ++next_elem_idx;

    } else if (word == "property") {
      CHECK_EQ(end_header, false);
      CHECK_NE(format_, kUNDEFINED_FORMAT);
      CHECK_LE(current_elem_idx, elements_.size());

      parseProperty(current_elem_idx, ss);

    } else if (word == "end_header") {
      end_header = true;
      break;
    }
  }

  CHECK_EQ(end_header, true);
  CHECK_NE(format_, kUNDEFINED_FORMAT);
  CHECK_GT(num_vertices_, 0);
  CHECK_GT(num_faces_, 0);
}

void PlyLoader::load(const std::string &filename, Data *d) {
  //
  file_.open(filename.c_str(), std::fstream::in | std::fstream::binary);
  CHECK(file_.is_open()) << filename;

  // update header_
  parseHeader();

  // load elements
  d->num_vertices = num_vertices_;
  d->num_faces = num_faces_;
  CHECK_NOTNULL(d->vertices);
  CHECK_NOTNULL(d->faces);

  for (auto &e : elements_) {
    if (e.name == kVERTEX) {
      parseVertices(e, d);
    } else if (e.name == kFACE) {
      parseFaces(e, d);
    } else {
      LOG(FATAL) << "unknown element name " << e.name;
    }
  }

  file_.close();
}

void PlyLoader::parseVertices(const Element &e, Data *d) {
  // check memory allocations
  CHECK_LE(e.num_elements * 3, d->vertices_capacity);
  if (e.has_color && d->colors) {
    CHECK_LE(e.num_elements, d->colors_capacity);
  }

  float x, y, z;
  uint32_t r, g, b, rgb;
  std::string line;

  if (format_ == kLITTLE_ENDIAN) {
    for (std::size_t n = 0, idx = 0; n < e.num_elements; ++n, idx += 3) {
      x = y = z = 0.0f;

      file_.read((char *)&x, 4);  // assume float
      file_.read((char *)&y, 4);
      file_.read((char *)&z, 4);

      d->vertices[idx] = x;
      d->vertices[idx + 1] = y;
      d->vertices[idx + 2] = z;

#ifdef DEBUG_PLY_LOADER
      LOG_FIRST_N(INFO, 1) << " has_color " << has_color;
      LOG_FIRST_N(INFO, 20) << x << " " << y << " " << z;
#endif
      // printf("%.9g %.9g %.9g\n", x,y,z);

      // user may have assigned nullptr to d->colors
      if (e.has_color) {
        r = g = b = 0;
        file_.read((char *)&r, 1);  // assume unsigned char
        file_.read((char *)&g, 1);
        file_.read((char *)&b, 1);

        CHECK_LT(r, 256);
        CHECK_LT(g, 256);
        CHECK_LT(b, 256);

        if (d->colors) {
          rgb = (r << 16) | (g << 8) | b;
          d->colors[n] = rgb;
        }
      }
    }
  } else if (format_ == kASCII) {
    for (std::size_t n = 0, idx = 0; n < e.num_elements; ++n, idx += 3) {
      std::getline(file_, line);
      std::istringstream ss(line);
      ss >> x >> y >> z;

      d->vertices[idx] = x;
      d->vertices[idx + 1] = y;
      d->vertices[idx + 2] = z;

      // user may have assigned nullptr to d->colors
      if (e.has_color) {
        ss >> r >> g >> b;

        if (d->colors) {
          rgb = (r << 16) | (g << 8) | b;
          d->colors[n] = rgb;
        }
      }
    }
  } else {
    LOG(FATAL) << "unknown format " << format_;
  }
}

void PlyLoader::parseFaces(const Element &e, Data *d) {
  // check memory allocations
  CHECK_LE(e.num_elements * 3, d->faces_capacity);

  int num_indices;
  uint32_t a, b, c;
  std::string line;

  if (format_ == kLITTLE_ENDIAN) {
    for (std::size_t n = 0, idx = 0; n < e.num_elements; ++n, idx += 3) {
      num_indices = 0;
      file_.read((char *)&num_indices, getDataSize(e.list_size_dtype));
      CHECK_EQ(num_indices, 3) << getDataSize(e.list_size_dtype);
#ifdef DEBUG_PLY_LOADER
      LOG_FIRST_N(INFO, 20) << "num_indices: " << num_indices;
#endif
      int bytes = getDataSize(e.index_dtype);
      a = b = c = 0;
      file_.read((char *)&a, bytes);
      file_.read((char *)&b, bytes);
      file_.read((char *)&c, bytes);

      d->faces[idx] = a;
      d->faces[idx + 1] = b;
      d->faces[idx + 2] = c;
#ifdef DEBUG_PLY_LOADER
      LOG_FIRST_N(INFO, 60)
          << "index[" << n << "]" << a << " " << b << " " << c;
#endif
    }
  } else if (format_ == kASCII) {
    for (std::size_t n = 0, idx = 0; n < e.num_elements; ++n, idx += 3) {
      std::getline(file_, line);
      std::istringstream ss(line);

      ss >> num_indices;
      CHECK_EQ(num_indices, 3);

      ss >> a >> b >> c;

      d->faces[idx] = a;
      d->faces[idx + 1] = b;
      d->faces[idx + 2] = c;
    }
  } else {
    LOG(FATAL) << "unknown format " << format_;
  }
}

void PlyLoader::parseElement(std::size_t elem_idx, std::istringstream &ss) {
  // the first word of the line has already read from ss.
  // ss is pointing to the second word at this point.
  // e.g. element [face] 6

  Element &e = elements_[elem_idx];
  std::string name;

  ss >> name >> e.num_elements;

  if (name == "vertex") {
    e.name = kVERTEX;
    e.has_color = false;
    e.vertex_dtype = kUNDEFINED_DATA_TYPE;
    e.color_dtype = kUNDEFINED_DATA_TYPE;
    e.property = kUNDEFINED_PROPERTY;

    num_vertices_ = e.num_elements;

  } else if (name == "face") {
    e.name = kFACE;
    e.list_size_dtype = kUNDEFINED_DATA_TYPE;
    e.index_dtype = kUNDEFINED_DATA_TYPE;

    num_faces_ = e.num_elements;

  } else {
    LOG(FATAL) << "unknown element: " << name;
  }
}

void PlyLoader::parseProperty(std::size_t elem_idx, std::istringstream &ss) {
  // the first word of the line has already read from ss.
  // ss is pointing to the second word at this point.
  // e.g. property [float] x,

  // current element
  Element &e = elements_[elem_idx];

  if (e.name == kVERTEX) {
    // property float x
    // property float y
    // property float z
    // property uchar red
    // property uchar green
    // property uchar blue

    // first word after the property identifier
    std::string dtype, name;
    ss >> dtype >> name;

    if (name == "x" || name == "y" || name == "z") {
      if (e.vertex_dtype == kUNDEFINED_DATA_TYPE) {
        e.vertex_dtype = getDataType(dtype);
        CHECK_EQ(e.vertex_dtype, kFLOAT);
      } else {  // check consistent data types for all xyz coordinates
        CHECK_EQ(e.vertex_dtype, getDataType(dtype));
      }
      // make sure the ordering (x->y->z)
      if (name == "x") {
        e.property = kX;
      } else if (name == "y") {
        CHECK_EQ(e.property, kX);
        e.property = kY;
      } else {
        CHECK_EQ(e.property, kY);
        e.property = kZ;
      }
    } else if (name == "red" || name == "green" || name == "blue") {
      e.has_color = true;
      if (e.color_dtype == kUNDEFINED_DATA_TYPE) {
        e.color_dtype = getDataType(dtype);
        CHECK_EQ(e.color_dtype, kUCHAR);
      } else {  // check consistent data types for all rgb color components
        CHECK_EQ(e.color_dtype, getDataType(dtype));
      }
      // make sure the ordering (r->g->b)
      if (name == "red") {
        e.property = kRED;
      } else if (name == "green") {
        CHECK_EQ(e.property, kRED);
        e.property = kGREEN;
      } else {
        CHECK_EQ(e.property, kGREEN);
        e.property = kBLUE;
      }
    } else {
      LOG(FATAL) << "unknown property " << name;
    }
  } else if (e.name == kFACE) {
    // property list uchar int vertex_index

    // list
    std::string word;
    ss >> word;
    CHECK_EQ(word, "list");

    // list size type
    ss >> word;
    e.list_size_dtype = getDataType(word);

    // index type
    ss >> word;
    e.index_dtype = getDataType(word);

    // property name
    ss >> word;
    CHECK(word == "vertex_index" || word == "vertex_indices")
        << "unknown property " << word;
  } else {
    LOG(FATAL) << "unknown element name " << e.name;
  }
}

int PlyLoader::getFormat(const std::string &s) const {
  int format;
  if (s == "ascii") {
    format = kASCII;
  } else if (s == "binary_little_endian") {
    format = kLITTLE_ENDIAN;
  } else if (s == "binary_big_endian") {
    format = kBIG_ENDIAN;
  }
  return format;
}

int PlyLoader::getDataType(const std::string &s) const {
  int dtype;
  // if (s == "list") {
  //   dtype = kLIST;
  if (s == "char") {
    dtype = kCHAR;
  } else if (s == "uchar") {
    dtype = kUCHAR;
  } else if (s == "short") {
    dtype = kSHORT;
  } else if (s == "ushort") {
    dtype = kUSHORT;
  } else if (s == "int") {
    dtype = kINT;
  } else if (s == "uint") {
    dtype = kUINT;
  } else if (s == "float") {
    dtype = kFLOAT;
  } else if (s == "double") {
    dtype = kDOUBLE;
  } else {
    LOG(FATAL) << "unknown data type " << s;
  }
  return dtype;
}

int PlyLoader::getDataSize(int type) {
  int bytes;
  if (type == kCHAR) {
    bytes = 1;  // 1 byte(s)
  } else if (type == kUCHAR) {
    bytes = 1;  // 1
  } else if (type == kSHORT) {
    bytes = 2;  // 2
  } else if (type == kUSHORT) {
    bytes = 2;  // 2
  } else if (type == kINT) {
    bytes = 4;  // 4
  } else if (type == kUINT) {
    bytes = 4;  // 4
  } else if (type == kFLOAT) {
    bytes = 4;  // 4
  } else if (type == kDOUBLE) {
    bytes = 8;  // 8
  } else {
    LOG(FATAL) << "unknown data type " << type;
  }
  return bytes;
}

}  // namespace spray

