#include <fstream>
#include <cstdio>

#include "glog/logging.h"
#include "gtest/gtest.h"

#include "io/obj_loader.h"

namespace spray {

class CubeTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
    std::ofstream file("cube.obj");

    // clang-format off
    std::string obj_string=
      "# cube.obj\n"
      "#\n"
      "\n"
      "g cube\n"
      "\n"
      "v  0.0  0.0  0.0\n"
      "v  0.0  0.0  1.0\n"
      "v  0.0  1.0  0.0\n"
      "v  0.0  1.0  1.0\n"
      "v  1.0  0.0  0.0\n"
      "v  1.0  0.0  1.0\n"
      "v  1.0  1.0  0.0\n"
      "v  1.0  1.0  1.0\n"
      "\n"
      "vn  0.0  0.0  1.0\n"
      "vn  0.0  0.0 -1.0\n"
      "vn  0.0  1.0  0.0\n"
      "vn  0.0 -1.0  0.0\n"
      "vn  1.0  0.0  0.0\n"
      "vn -1.0  0.0  0.0\n"
      "\n"
      "f  1//2  7//2  5//2\n"
      "f  1//2  3//2  7//2\n"
      "f  1//6  4//6  3//6\n"
      "f  1//6  2//6  4//6\n"
      "f  3//3  8//3  7//3\n"
      "f  3//3  4//3  8//3\n"
      "f  5//5  7//5  8//5\n"
      "f  5//5  8//5  6//5\n"
      "f  1//4  5//4  6//4\n"
      "f  1//4  6//4  2//4\n"
      "f  2//1  6//1  8//1\n"
      "f  2//1  8//1  4//1\n";

    // clang-format on

    file << obj_string;
    file.close();

    data_.vertices = &vertices_;
    data_.normals = &normals_;
    data_.vertex_indices = &vertex_indices_;
    data_.texture_indices = &texture_indices_;
    data_.normal_indices = &normal_indices_;

    loader_.load("cube.obj", &data_);
  }

  virtual void TearDown() { std::remove("cube.obj"); }

  ObjLoader loader_;
  ObjLoader::Data data_;

  std::vector<float> vertices_;
  std::vector<float> normals_;
  std::vector<int> vertex_indices_;
  std::vector<int> texture_indices_;
  std::vector<int> normal_indices_;
};

TEST_F(CubeTest, Sizes) {
  ASSERT_EQ(vertices_.size(), 8 * 3);
  ASSERT_EQ(normals_.size(), 6 * 3);
  ASSERT_EQ(vertex_indices_.size(), 12 * 3);
  ASSERT_EQ(texture_indices_.size(), 0);
  ASSERT_EQ(normal_indices_.size(), 12 * 3);
}

TEST_F(CubeTest, VertexValues) {
  // clang-format off
  std::vector<float> v
    ={0.0,  0.0,  0.0,
      0.0,  0.0,  1.0,
      0.0,  1.0,  0.0,
      0.0,  1.0,  1.0,
      1.0,  0.0,  0.0,
      1.0,  0.0,  1.0,
      1.0,  1.0,  0.0,
      1.0,  1.0,  1.0};
  // clang-format on

  ASSERT_EQ(vertices_.size(), v.size());

  for (std::size_t i = 0; i < vertices_.size(); ++i) {
    ASSERT_EQ(vertices_[i], v[i]);
  }
}

TEST_F(CubeTest, NormalValues) {
  // clang-format off
  std::vector<float> v
     ={0.0,  0.0,  1.0,
       0.0,  0.0, -1.0,
       0.0,  1.0,  0.0,
       0.0, -1.0,  0.0,
       1.0,  0.0,  0.0,
      -1.0,  0.0,  0.0};
  // clang-format on

  ASSERT_EQ(normals_.size(), v.size());

  for (std::size_t i = 0; i < normals_.size(); ++i) {
    ASSERT_EQ(normals_[i], v[i]);
  }
}

TEST_F(CubeTest, FaceVertexIndices) {
  // clang-format off
  std::vector<float> v
    ={1,  7,  5,
      1,  3,  7,
      1,  4,  3,
      1,  2,  4,
      3,  8,  7,
      3,  4,  8,
      5,  7,  8,
      5,  8,  6,
      1,  5,  6,
      1,  6,  2,
      2,  6,  8,
      2,  8,  4};
  // clang-format on

  ASSERT_EQ(vertex_indices_.size(), v.size());

  for (std::size_t i = 0; i < vertex_indices_.size(); ++i) {
    ASSERT_EQ(vertex_indices_[i], v[i]);
  }
}

TEST_F(CubeTest, FaceNormalIndices) {
  // clang-format off
  std::vector<float> v
    ={2,  2,  2,
      2,  2,  2,
      6,  6,  6,
      6,  6,  6,
      3,  3,  3,
      3,  3,  3,
      5,  5,  5,
      5,  5,  5,
      4,  4,  4,
      4,  4,  4,
      1,  1,  1,
      1,  1,  1};
  // clang-format on

  ASSERT_EQ(normal_indices_.size(), v.size());

  for (std::size_t i = 0; i < normal_indices_.size(); ++i) {
    ASSERT_EQ(normal_indices_[i], v[i]);
  }
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  google::InitGoogleLogging(argv[0]);
  google::InstallFailureSignalHandler();
  return RUN_ALL_TESTS();
}

}  // namespace spray

