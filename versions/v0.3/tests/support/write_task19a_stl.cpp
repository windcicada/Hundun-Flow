// SPDX-License-Identifier: Apache-2.0

#include "tests/support/stage3_stl_fixture.hpp"

#include "hundun/rt_error.hpp"

#include <filesystem>
#include <string>
#include <vector>

namespace {

hundun::runtime::Real3 midpoint(hundun::runtime::Real3 left,
                                hundun::runtime::Real3 right) {
  return {0.5 * (left.x + right.x), 0.5 * (left.y + right.y),
          0.5 * (left.z + right.z)};
}

std::vector<hundun::test::StlFixtureTriangle> refined_internal_cube() {
  auto coarse = hundun::test::outward_cube();
  for (auto& triangle : coarse) {
    for (auto& vertex : triangle.vertices) {
      vertex.x = 0.30 + 0.40 * vertex.x;
      vertex.y = 0.30 + 0.40 * vertex.y;
      vertex.z = 0.30 + 0.40 * vertex.z;
    }
  }
  std::vector<hundun::test::StlFixtureTriangle> result;
  result.reserve(4U * coarse.size());
  for (const auto& triangle : coarse) {
    const auto ab = midpoint(triangle.vertices[0], triangle.vertices[1]);
    const auto bc = midpoint(triangle.vertices[1], triangle.vertices[2]);
    const auto ca = midpoint(triangle.vertices[2], triangle.vertices[0]);
    result.push_back({triangle.file_normal, {triangle.vertices[0], ab, ca}});
    result.push_back({triangle.file_normal, {ab, triangle.vertices[1], bc}});
    result.push_back({triangle.file_normal, {ca, bc, triangle.vertices[2]}});
    result.push_back({triangle.file_normal, {ab, bc, ca}});
  }
  return result;
}

} // namespace

int main(int argc, char** argv) {
  if (argc != 2) {
    throw hundun::runtime::Error("usage: write_task19a_stl <path>");
  }
  const std::filesystem::path path(argv[1]);
  hundun::test::write_text(
      path, hundun::test::ascii_stl(refined_internal_cube(), "body"));
  return 0;
}
