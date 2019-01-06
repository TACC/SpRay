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

#include "io/obj_loader.h"

#include <cstring>
#include <fstream>
#include <iostream>
#include <string>

#include "glog/logging.h"

namespace spray {

void ObjLoader::load(const std::string &filename, Data *data) {
  data_ = data;

  std::ifstream infile(filename);
  CHECK(infile.is_open()) << "unable to open input file " << filename;

  char delim[] = " \t\r\n\v\f/";

  std::queue<std::string> tokens;
  std::string line, original_line;
  while (infile.good()) {
    getline(infile, line);
    original_line = line;
#ifdef SPRAY_PRINT_LINES
    std::cout << line << "\n";
#endif
    char *token = std::strtok(&line[0], delim);
    CHECK(tokens.empty());

    while (token != NULL) {
#ifdef SPRAY_PRINT_TOKENS
      std::cout << token << " token len:" << std::string(token).size() << "\n";
#endif
      tokens.push(token);
      token = std::strtok(NULL, delim);
    }

    if (!tokens.empty()) parseTokens(original_line, &tokens);
  }

  infile.close();
  postProcessing();

  CHECK(vertices_.empty());
  CHECK(normals_.empty());
  CHECK(vertex_indices_.empty());
  CHECK(texture_indices_.empty());
  CHECK(normal_indices_.empty());
}

void ObjLoader::parseTokens(const std::string &line, StringQ *tokens) {
  std::string &token = tokens->front();
  tokens->pop();

  if (token == "v") {
    parseVertices(tokens);
  } else if (token == "vn") {
    parseNormals(tokens);
  } else if (token == "f") {
    parseFaces(line, tokens);
  }
  flush(tokens);
}

void ObjLoader::flush(StringQ *tokens) {
  while (!tokens->empty()) {
    tokens->pop();
  }
}

void ObjLoader::parseVertices(StringQ *tokens) {
  CHECK_EQ(tokens->size(), 3);
  while (!tokens->empty()) {
    std::string &token = tokens->front();
    tokens->pop();
    vertices_.push(atof(token.c_str()));
  }
}

void ObjLoader::parseNormals(StringQ *tokens) {
  CHECK_EQ(tokens->size(), 3);
  while (!tokens->empty()) {
    std::string &token = tokens->front();
    tokens->pop();
    normals_.push(atof(token.c_str()));
  }
}

void ObjLoader::parseFaces(const std::string &line, StringQ *tokens) {
  std::size_t num_tokens = tokens->size();
  CHECK(num_tokens == 3 || num_tokens == 6 || num_tokens == 9) << num_tokens;
  if (num_tokens == 3) {
    while (!tokens->empty()) {
      std::string &token = tokens->front();
      tokens->pop();
      vertex_indices_.push(atoi(token.c_str()));
    }
  } else if (num_tokens == 6) {
    std::size_t found = line.find("//");
    bool has_normals = (found != std::string::npos);
    int i = 0;
    while (!tokens->empty()) {
      std::string &token = tokens->front();
      tokens->pop();
      // std::cout << "token: " << token << "\n";
      if (i == 0) {
        vertex_indices_.push(atoi(token.c_str()) - 1);
      } else {
        if (has_normals) {
          normal_indices_.push(atoi(token.c_str()) - 1);
        } else {
          texture_indices_.push(atoi(token.c_str()) - 1);
        }
      }
      ++i;
      if (i == 2) i = 0;
    }
    CHECK_EQ(i, 0);
  } else if (num_tokens == 9) {
    int i = 0;
    while (!tokens->empty()) {
      std::string &token = tokens->front();
      tokens->pop();
      if (i == 0) {
        vertex_indices_.push(atoi(token.c_str()));
      } else if (i == 1) {
        texture_indices_.push(atoi(token.c_str()));
      } else {
        normal_indices_.push(atoi(token.c_str()));
      }
      if (i == 3) i = 0;
    }
    CHECK_EQ(i, 0);
  }
}

void ObjLoader::postProcessing() {
  // vertices
  std::vector<float> *vout = data_->vertices;
  vout->reserve(vertices_.size());

  while (!vertices_.empty()) {
    vout->push_back(vertices_.front());
    vertices_.pop();
  }

  // normals
  vout = data_->normals;
  vout->reserve(normals_.size());

  while (!normals_.empty()) {
    vout->push_back(normals_.front());
    normals_.pop();
  }

  // vertex indices
  std::vector<int> *iout = data_->vertex_indices;
  iout->reserve(vertex_indices_.size());

  while (!vertex_indices_.empty()) {
    iout->push_back(vertex_indices_.front());
    vertex_indices_.pop();
  }

  // texture indices
  iout = data_->texture_indices;
  iout->reserve(texture_indices_.size());

  while (!texture_indices_.empty()) {
    iout->push_back(texture_indices_.front());
    texture_indices_.pop();
  }

  // normal indices
  iout = data_->normal_indices;
  iout->reserve(normal_indices_.size());

  while (!normal_indices_.empty()) {
    iout->push_back(normal_indices_.front());
    normal_indices_.pop();
  }
}

}  // namespace spray
