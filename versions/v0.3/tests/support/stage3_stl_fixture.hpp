// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/rt_types.hpp"

#include <array>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace hundun::test {

struct StlFixtureTriangle final {
  runtime::Real3 file_normal{};
  std::array<runtime::Real3, 3> vertices{};
};

class Stage3TemporaryDirectory final {
public:
  explicit Stage3TemporaryDirectory(std::string name);
  ~Stage3TemporaryDirectory() noexcept;
  Stage3TemporaryDirectory(const Stage3TemporaryDirectory &) = delete;
  Stage3TemporaryDirectory &
  operator=(const Stage3TemporaryDirectory &) = delete;

  const std::filesystem::path &path() const noexcept { return path_; }

private:
  std::filesystem::path path_;
};

std::vector<StlFixtureTriangle> outward_tetrahedron();
std::vector<StlFixtureTriangle> outward_cube();
std::vector<StlFixtureTriangle>
projected_octahedral_sphere(runtime::Real3 center, double radius,
                            unsigned refinement_levels);
std::vector<StlFixtureTriangle>
translated(std::vector<StlFixtureTriangle> triangles, runtime::Real3 offset);

std::string ascii_stl(const std::vector<StlFixtureTriangle> &triangles,
                      std::string name = "hundun");
std::vector<std::uint8_t>
binary_stl(const std::vector<StlFixtureTriangle> &triangles,
           bool header_starts_with_solid = false);

void write_text(const std::filesystem::path &path, const std::string &text);
void write_bytes(const std::filesystem::path &path,
                 const std::vector<std::uint8_t> &bytes);

} // namespace hundun::test
