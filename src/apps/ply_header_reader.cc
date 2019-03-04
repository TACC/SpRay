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

#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>

#include "glm/glm.hpp"
#include "glog/logging.h"

#include "io/ply_loader.h"
#include "render/aabb.h"

void printUsage(char** argv) {
  printf("Usage: %s <infile> <outfile>\n", argv[0]);
}

int main(int argc, char** argv) {
  google::InitGoogleLogging(argv[0]);

  if (argc < 2) {
    printUsage(argv);
    std::cout << "[error] invalid commandline\n";
    return 0;
  }

  std::string infile(argv[1]);

  std::string outfile;
  if (argc > 2) {
    outfile = std::string(argv[2]);
  } else {
    outfile = std::string("./ply_info.txt");
  }

  std::vector<float> vertices;
  std::vector<uint32_t> faces;

  std::cout << "[info] loading " << infile << "\n";
  // std::cout << "[info] reading header\n";

  spray::PlyLoader::Header h;
  spray::PlyLoader::quickHeaderRead(infile, &h);

  CHECK_GT(h.num_vertices, 0);
  CHECK_GT(h.num_faces, 0);

  vertices.resize(h.num_vertices * 3);
  faces.resize(h.num_faces * 3);

  spray::PlyLoader::Data data;
  data.vertices_capacity = vertices.size();  // in
  data.faces_capacity = faces.size();        // in
  data.colors_capacity = 0;                  // in
  data.vertices = &vertices[0];              // in/out
  // std::size_t num_vertices;     // out
  data.faces = &faces[0];  // in/out
  // std::size_t num_faces;  // out
  data.colors = nullptr;  // rgb, in/out

  // std::cout << "[info] loading..\n";
  spray::PlyLoader loader;
  loader.load(infile, &data);

  CHECK_EQ(data.num_vertices, h.num_vertices);
  CHECK_EQ(data.num_faces, h.num_faces);

  // merge vertices
  spray::Aabb aabb;
  for (std::size_t i = 0; i < vertices.size(); i += 3) {
    aabb.merge(glm::vec3(vertices[i], vertices[i + 1], vertices[i + 2]));
  }

  std::cout << "[info] dumping ply info " << outfile << "\n";

  std::ofstream fout;
  fout.open(outfile);

  fout << "file " << infile << "\n";
  fout << "vertex " << h.num_vertices << "\n";
  fout << "face " << h.num_faces << "\n";
  fout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
       << "bounds " << aabb.bounds[0].x << " " << aabb.bounds[0].y << " "
       << aabb.bounds[0].z << " " << aabb.bounds[1].x << " " << aabb.bounds[1].y
       << " " << aabb.bounds[1].z << "\n";

  fout.close();

  return 0;
}
