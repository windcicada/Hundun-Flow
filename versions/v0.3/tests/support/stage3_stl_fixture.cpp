// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "tests/support/stage3_stl_fixture.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

namespace hundun::test {
namespace {

runtime::Real3 subtract(runtime::Real3 a, runtime::Real3 b) noexcept {
  return {a.x - b.x, a.y - b.y, a.z - b.z};
}

runtime::Real3 cross(runtime::Real3 a, runtime::Real3 b) noexcept {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

runtime::Real3 unit_normal(const StlFixtureTriangle &triangle) noexcept {
  const runtime::Real3 normal =
      cross(subtract(triangle.vertices[1], triangle.vertices[0]),
            subtract(triangle.vertices[2], triangle.vertices[0]));
  const double length = std::hypot(normal.x, normal.y, normal.z);
  return {normal.x / length, normal.y / length, normal.z / length};
}

void append_u16(std::vector<std::uint8_t> &bytes, std::uint16_t value) {
  bytes.push_back(static_cast<std::uint8_t>(value & 0xffU));
  bytes.push_back(static_cast<std::uint8_t>((value >> 8U) & 0xffU));
}

void append_u32(std::vector<std::uint8_t> &bytes, std::uint32_t value) {
  for (unsigned shift = 0U; shift < 32U; shift += 8U) {
    bytes.push_back(static_cast<std::uint8_t>((value >> shift) & 0xffU));
  }
}

void append_f32(std::vector<std::uint8_t> &bytes, double source) {
  const float value = static_cast<float>(source);
  std::uint32_t bits = 0U;
  static_assert(sizeof(bits) == sizeof(value));
  std::memcpy(&bits, &value, sizeof(bits));
  append_u32(bytes, bits);
}

StlFixtureTriangle face(runtime::Real3 a, runtime::Real3 b, runtime::Real3 c) {
  StlFixtureTriangle triangle;
  triangle.vertices = {a, b, c};
  triangle.file_normal = unit_normal(triangle);
  return triangle;
}

} // namespace

Stage3TemporaryDirectory::Stage3TemporaryDirectory(std::string name) {
  static std::atomic<std::uint64_t> next{0U};
  const auto timestamp = static_cast<std::uint64_t>(
      std::chrono::steady_clock::now().time_since_epoch().count());
  path_ =
      std::filesystem::temp_directory_path() /
      ("hundun-stage3-" + std::move(name) + "-" + std::to_string(timestamp) +
       "-" + std::to_string(next.fetch_add(1U)));
  std::filesystem::create_directories(path_);
}

Stage3TemporaryDirectory::~Stage3TemporaryDirectory() noexcept {
  std::error_code error;
  std::filesystem::remove_all(path_, error);
}

std::vector<StlFixtureTriangle> outward_tetrahedron() {
  const runtime::Real3 v0{0.0, 0.0, 0.0};
  const runtime::Real3 v1{1.0, 0.0, 0.0};
  const runtime::Real3 v2{0.0, 1.0, 0.0};
  const runtime::Real3 v3{0.0, 0.0, 1.0};
  return {face(v0, v2, v1), face(v0, v1, v3), face(v0, v3, v2),
          face(v1, v2, v3)};
}

std::vector<StlFixtureTriangle> outward_cube() {
  const std::array<runtime::Real3, 8> v{
      runtime::Real3{0.0, 0.0, 0.0}, runtime::Real3{1.0, 0.0, 0.0},
      runtime::Real3{1.0, 1.0, 0.0}, runtime::Real3{0.0, 1.0, 0.0},
      runtime::Real3{0.0, 0.0, 1.0}, runtime::Real3{1.0, 0.0, 1.0},
      runtime::Real3{1.0, 1.0, 1.0}, runtime::Real3{0.0, 1.0, 1.0}};
  return {
      face(v[0], v[2], v[1]), face(v[0], v[3], v[2]), face(v[4], v[5], v[6]),
      face(v[4], v[6], v[7]), face(v[0], v[1], v[5]), face(v[0], v[5], v[4]),
      face(v[3], v[7], v[6]), face(v[3], v[6], v[2]), face(v[0], v[4], v[7]),
      face(v[0], v[7], v[3]), face(v[1], v[2], v[6]), face(v[1], v[6], v[5])};
}

std::vector<StlFixtureTriangle>
projected_octahedral_sphere(runtime::Real3 center, double radius,
                            unsigned refinement_levels) {
  if (!(radius > 0.0) || !std::isfinite(radius))
    throw std::invalid_argument(
        "sphere fixture radius must be finite and positive");
  const runtime::Real3 px{center.x + radius, center.y, center.z};
  const runtime::Real3 nx{center.x - radius, center.y, center.z};
  const runtime::Real3 py{center.x, center.y + radius, center.z};
  const runtime::Real3 ny{center.x, center.y - radius, center.z};
  const runtime::Real3 pz{center.x, center.y, center.z + radius};
  const runtime::Real3 nz{center.x, center.y, center.z - radius};
  std::vector<StlFixtureTriangle> triangles{
      face(pz, px, py), face(pz, py, nx), face(pz, nx, ny),
      face(pz, ny, px), face(nz, py, px), face(nz, nx, py),
      face(nz, ny, nx), face(nz, px, ny)};
  const auto project = [center, radius](runtime::Real3 first,
                                        runtime::Real3 second) {
    const runtime::Real3 midpoint{(first.x + second.x) * 0.5,
                                  (first.y + second.y) * 0.5,
                                  (first.z + second.z) * 0.5};
    const auto radial = subtract(midpoint, center);
    const double length = std::hypot(radial.x, radial.y, radial.z);
    if (!(length > 0.0) || !std::isfinite(length))
      throw std::runtime_error("sphere fixture projection is degenerate");
    return runtime::Real3{center.x + radius * radial.x / length,
                          center.y + radius * radial.y / length,
                          center.z + radius * radial.z / length};
  };
  for (unsigned level = 0U; level < refinement_levels; ++level) {
    std::vector<StlFixtureTriangle> refined;
    refined.reserve(4U * triangles.size());
    for (const auto &triangle : triangles) {
      const auto ab = project(triangle.vertices[0], triangle.vertices[1]);
      const auto bc = project(triangle.vertices[1], triangle.vertices[2]);
      const auto ca = project(triangle.vertices[2], triangle.vertices[0]);
      refined.push_back(face(triangle.vertices[0], ab, ca));
      refined.push_back(face(ab, triangle.vertices[1], bc));
      refined.push_back(face(ca, bc, triangle.vertices[2]));
      refined.push_back(face(ab, bc, ca));
    }
    triangles = std::move(refined);
  }
  return triangles;
}

std::vector<StlFixtureTriangle>
translated(std::vector<StlFixtureTriangle> triangles, runtime::Real3 offset) {
  for (StlFixtureTriangle &triangle : triangles) {
    for (runtime::Real3 &vertex : triangle.vertices) {
      vertex.x += offset.x;
      vertex.y += offset.y;
      vertex.z += offset.z;
    }
  }
  return triangles;
}

std::string ascii_stl(const std::vector<StlFixtureTriangle> &triangles,
                      std::string name) {
  std::ostringstream stream;
  stream.imbue(std::locale::classic());
  stream << std::setprecision(std::numeric_limits<double>::max_digits10)
         << std::scientific;
  stream << "solid " << name << '\n';
  for (const StlFixtureTriangle &triangle : triangles) {
    stream << "  facet normal " << triangle.file_normal.x << ' '
           << triangle.file_normal.y << ' ' << triangle.file_normal.z << '\n'
           << "    outer loop\n";
    for (const runtime::Real3 vertex : triangle.vertices) {
      stream << "      vertex " << vertex.x << ' ' << vertex.y << ' '
             << vertex.z << '\n';
    }
    stream << "    endloop\n"
           << "  endfacet\n";
  }
  stream << "endsolid " << name << '\n';
  return stream.str();
}

std::vector<std::uint8_t>
binary_stl(const std::vector<StlFixtureTriangle> &triangles,
           bool header_starts_with_solid) {
  if (triangles.size() > std::numeric_limits<std::uint32_t>::max()) {
    throw std::runtime_error("fixture triangle count exceeds uint32");
  }
  std::vector<std::uint8_t> bytes(80U, 0U);
  const std::string header =
      header_starts_with_solid ? "solid binary hundun" : "HUNDUN binary STL";
  std::copy(header.begin(), header.end(), bytes.begin());
  append_u32(bytes, static_cast<std::uint32_t>(triangles.size()));
  for (const StlFixtureTriangle &triangle : triangles) {
    append_f32(bytes, triangle.file_normal.x);
    append_f32(bytes, triangle.file_normal.y);
    append_f32(bytes, triangle.file_normal.z);
    for (const runtime::Real3 vertex : triangle.vertices) {
      append_f32(bytes, vertex.x);
      append_f32(bytes, vertex.y);
      append_f32(bytes, vertex.z);
    }
    append_u16(bytes, 0U);
  }
  return bytes;
}

void write_text(const std::filesystem::path &path, const std::string &text) {
  std::ofstream stream(path, std::ios::binary | std::ios::trunc);
  stream.write(text.data(), static_cast<std::streamsize>(text.size()));
  if (!stream) {
    throw std::runtime_error("unable to write ASCII STL fixture");
  }
}

void write_bytes(const std::filesystem::path &path,
                 const std::vector<std::uint8_t> &bytes) {
  std::ofstream stream(path, std::ios::binary | std::ios::trunc);
  if (!bytes.empty()) {
    stream.write(reinterpret_cast<const char *>(bytes.data()),
                 static_cast<std::streamsize>(bytes.size()));
  }
  if (!stream) {
    throw std::runtime_error("unable to write binary STL fixture");
  }
}

} // namespace hundun::test
