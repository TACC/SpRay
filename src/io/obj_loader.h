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

#include <cstdlib>
#include <queue>
#include <string>

namespace spray {

class ObjLoader {
 public:
  struct Data {
    std::vector<float> *vertices;       // out
    std::vector<float> *normals;        // out
    std::vector<int> *vertex_indices;   // out
    std::vector<int> *texture_indices;  // out (optional)
    std::vector<int> *normal_indices;   // out
  };

  void load(const std::string &filename, Data *data);

 private:
  typedef std::queue<std::string> StringQ;

  void parseObjTokens(const std::string &line, StringQ *tokens);
  void parseVertices(StringQ *tokens);
  void parseNormals(StringQ *tokens);
  void parseFaces(const std::string &line, StringQ *tokens);
  void parseMaterialFile(StringQ *tokens);
  void loadMaterialFile(const std::string &filename);
  void postProcessing();

  void parseMaterialTokens(StringQ *tokens);
  void parseNewMtl(StringQ *tokens);

  template <typename T>
  void flush(std::queue<T> *q) {
    while (!q->empty()) q->pop();
  }

 private:
  std::queue<float> vertices_;
  std::queue<float> normals_;
  std::queue<int> vertex_indices_;
  std::queue<int> texture_indices_;
  std::queue<int> normal_indices_;
  Data *data_;

  // std::map<std::string, Material> materials_;
};

}  // namespace spray

