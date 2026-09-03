// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/ib_surface.hpp"

#include "hundun/rt_collective_status.hpp"
#include "hundun/rt_error.hpp"

#include "ib_stl_reader_detail.hpp"
#include "ib_surface_bvh_detail.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <ios>
#include <limits>
#include <locale>
#include <map>
#include <memory>
#include <queue>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace hundun::immersed::detail {
namespace {

using Real3 = runtime::Real3;

constexpr std::uint64_t fnv_offset = UINT64_C(14695981039346656037);
constexpr std::uint64_t fnv_prime = UINT64_C(1099511628211);
constexpr std::array<std::uint8_t, 8> payload_magic{'H', 'D', 'S', 'U',
                                                    'R', 'F', '3', '\0'};
constexpr std::uint32_t payload_version = 1U;

struct RawTriangle final {
  std::array<Real3, 3> vertices{};
  Real3 file_normal{};
};

struct EdgeOccurrence final {
  std::uint64_t from{};
  std::uint64_t to{};
  std::uint64_t triangle{};
};

[[noreturn]] void fail(std::string_view message) {
  throw runtime::Error(std::string{"immersed surface: "} +
                       std::string{message});
}

std::uint64_t double_bits(double value) noexcept {
  std::uint64_t bits = 0;
  static_assert(sizeof(bits) == sizeof(value));
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

double double_from_bits(std::uint64_t bits) noexcept {
  double value = 0.0;
  std::memcpy(&value, &bits, sizeof(value));
  return value;
}

void hash_byte(std::uint64_t &hash, std::uint8_t value) noexcept {
  hash ^= value;
  hash *= fnv_prime;
}

void hash_u64(std::uint64_t &hash, std::uint64_t value) noexcept {
  for (unsigned shift = 0; shift < 64U; shift += 8U) {
    hash_byte(hash, static_cast<std::uint8_t>((value >> shift) & 0xffU));
  }
}

std::uint32_t read_u32(const std::vector<std::uint8_t> &bytes,
                       std::size_t offset) {
  if (offset > bytes.size() || bytes.size() - offset < 4U) {
    fail("truncated binary STL");
  }
  std::uint32_t value = 0;
  for (unsigned index = 0; index < 4U; ++index) {
    value |= static_cast<std::uint32_t>(bytes[offset + index]) << (8U * index);
  }
  return value;
}

float read_f32(const std::vector<std::uint8_t> &bytes, std::size_t offset) {
  const std::uint32_t bits = read_u32(bytes, offset);
  float value = 0.0F;
  static_assert(sizeof(value) == sizeof(bits));
  std::memcpy(&value, &bits, sizeof(value));
  return value;
}

std::vector<std::string> split_words(std::string_view line) {
  std::istringstream stream{std::string{line}};
  stream.imbue(std::locale::classic());
  std::vector<std::string> words;
  for (std::string word; stream >> word;) {
    words.push_back(std::move(word));
  }
  return words;
}

bool decimal_syntax(std::string_view token) noexcept {
  std::size_t index = 0;
  if (index < token.size() && (token[index] == '+' || token[index] == '-')) {
    ++index;
  }
  bool leading_digits = false;
  while (index < token.size() && token[index] >= '0' && token[index] <= '9') {
    leading_digits = true;
    ++index;
  }
  bool trailing_digits = false;
  if (index < token.size() && token[index] == '.') {
    ++index;
    while (index < token.size() && token[index] >= '0' && token[index] <= '9') {
      trailing_digits = true;
      ++index;
    }
  }
  if (!leading_digits && !trailing_digits) {
    return false;
  }
  if (index < token.size() && (token[index] == 'e' || token[index] == 'E')) {
    ++index;
    if (index < token.size() && (token[index] == '+' || token[index] == '-')) {
      ++index;
    }
    bool exponent_digits = false;
    while (index < token.size() && token[index] >= '0' && token[index] <= '9') {
      exponent_digits = true;
      ++index;
    }
    if (!exponent_digits) {
      return false;
    }
  }
  return index == token.size();
}

double parse_decimal(std::string_view token) {
  if (!decimal_syntax(token)) {
    fail("ASCII STL contains an invalid numeric token");
  }
  std::istringstream stream{std::string{token}};
  stream.imbue(std::locale::classic());
  double value = 0.0;
  stream >> value;
  if (!stream || !stream.eof() || !std::isfinite(value)) {
    fail("ASCII STL contains a non-finite or invalid number");
  }
  return value;
}

std::vector<std::string>
nonempty_lines(const std::vector<std::uint8_t> &bytes) {
  if (std::find(bytes.begin(), bytes.end(), std::uint8_t{0}) != bytes.end()) {
    fail("ASCII STL contains a NUL byte");
  }
  const std::string text(bytes.begin(), bytes.end());
  std::istringstream stream{text};
  stream.imbue(std::locale::classic());
  std::vector<std::string> lines;
  for (std::string line; std::getline(stream, line);) {
    if (!line.empty() && line.back() == '\r') {
      line.pop_back();
    }
    const auto first = line.find_first_not_of(" \t");
    if (first == std::string::npos) {
      continue;
    }
    const auto last = line.find_last_not_of(" \t");
    lines.push_back(line.substr(first, last - first + 1U));
  }
  return lines;
}

std::vector<RawTriangle>
parse_ascii_stl(const std::vector<std::uint8_t> &bytes) {
  const auto lines = nonempty_lines(bytes);
  if (lines.size() < 2U) {
    fail("ASCII STL is incomplete");
  }
  const auto first_words = split_words(lines.front());
  if (first_words.empty() || first_words.front() != "solid") {
    fail("ASCII STL must begin with solid");
  }

  std::vector<RawTriangle> triangles;
  std::size_t line = 1U;
  bool ended = false;
  while (line < lines.size()) {
    const auto words = split_words(lines[line]);
    if (!words.empty() && words.front() == "endsolid") {
      ++line;
      ended = true;
      break;
    }
    if (words.size() != 5U || words[0] != "facet" || words[1] != "normal") {
      fail("ASCII STL expected facet normal");
    }
    RawTriangle triangle;
    triangle.file_normal = {parse_decimal(words[2]), parse_decimal(words[3]),
                            parse_decimal(words[4])};
    ++line;
    if (line >= lines.size() ||
        split_words(lines[line]) != std::vector<std::string>{"outer", "loop"}) {
      fail("ASCII STL expected outer loop");
    }
    ++line;
    for (std::size_t vertex = 0; vertex < 3U; ++vertex, ++line) {
      if (line >= lines.size()) {
        fail("ASCII STL is truncated in a facet");
      }
      const auto vertex_words = split_words(lines[line]);
      if (vertex_words.size() != 4U || vertex_words[0] != "vertex") {
        fail("ASCII STL expected vertex");
      }
      triangle.vertices[vertex] = {parse_decimal(vertex_words[1]),
                                   parse_decimal(vertex_words[2]),
                                   parse_decimal(vertex_words[3])};
    }
    if (line >= lines.size() ||
        split_words(lines[line]) != std::vector<std::string>{"endloop"}) {
      fail("ASCII STL expected endloop");
    }
    ++line;
    if (line >= lines.size() ||
        split_words(lines[line]) != std::vector<std::string>{"endfacet"}) {
      fail("ASCII STL expected endfacet");
    }
    ++line;
    triangles.push_back(triangle);
  }
  if (!ended) {
    fail("ASCII STL is missing endsolid");
  }
  if (line != lines.size()) {
    fail("ASCII STL contains trailing tokens");
  }
  if (triangles.empty()) {
    fail("STL contains no triangles");
  }
  return triangles;
}

std::vector<RawTriangle>
parse_binary_stl(const std::vector<std::uint8_t> &bytes,
                 std::uint32_t triangle_count) {
  constexpr std::size_t header_bytes = 84U;
  constexpr std::size_t record_bytes = 50U;
  if (triangle_count >
      static_cast<std::uint32_t>(std::numeric_limits<std::int32_t>::max())) {
    fail("binary STL triangle count exceeds 2^31-1");
  }
  const std::uint64_t expected_u64 =
      UINT64_C(84) + static_cast<std::uint64_t>(triangle_count) * UINT64_C(50);
  if constexpr (sizeof(std::size_t) < sizeof(std::uint64_t)) {
    if (expected_u64 > std::numeric_limits<std::size_t>::max()) {
      fail("binary STL size overflows the host size type");
    }
  }
  if (expected_u64 > std::numeric_limits<std::streamsize>::max()) {
    fail("binary STL size overflows the host size type");
  }
  const std::size_t expected = static_cast<std::size_t>(expected_u64);
  if (bytes.size() != expected) {
    fail(bytes.size() < expected ? "binary STL is truncated"
                                 : "binary STL contains trailing bytes");
  }
  if (triangle_count == 0U) {
    fail("STL contains no triangles");
  }

  std::vector<RawTriangle> triangles;
  triangles.reserve(triangle_count);
  for (std::size_t index = 0; index < triangle_count; ++index) {
    const std::size_t base = header_bytes + index * record_bytes;
    RawTriangle triangle;
    triangle.file_normal = {static_cast<double>(read_f32(bytes, base)),
                            static_cast<double>(read_f32(bytes, base + 4U)),
                            static_cast<double>(read_f32(bytes, base + 8U))};
    for (std::size_t vertex = 0; vertex < 3U; ++vertex) {
      const std::size_t position = base + 12U + vertex * 12U;
      triangle.vertices[vertex] = {
          static_cast<double>(read_f32(bytes, position)),
          static_cast<double>(read_f32(bytes, position + 4U)),
          static_cast<double>(read_f32(bytes, position + 8U))};
    }
    triangles.push_back(triangle);
  }
  return triangles;
}

bool appears_binary(const std::vector<std::uint8_t> &bytes) noexcept {
  if (bytes.size() < 5U ||
      std::string_view{reinterpret_cast<const char *>(bytes.data()), 5U} !=
          "solid") {
    return true;
  }
  return std::any_of(bytes.begin(), bytes.end(), [](std::uint8_t value) {
    return value == 0U || (value < 0x09U || (value > 0x0dU && value < 0x20U));
  });
}

std::vector<RawTriangle> parse_stl(const std::vector<std::uint8_t> &bytes) {
  if (bytes.empty()) {
    fail("STL file is empty");
  }
  if (bytes.size() >= 84U) {
    const std::uint32_t count = read_u32(bytes, 80U);
    const std::uint64_t expected =
        UINT64_C(84) + static_cast<std::uint64_t>(count) * UINT64_C(50);
    if ((expected == bytes.size()) || appears_binary(bytes)) {
      return parse_binary_stl(bytes, count);
    }
  }
  return parse_ascii_stl(bytes);
}

std::array<std::int64_t, 3> bucket_key(Real3 value, Real3 minimum,
                                       double tolerance) {
  std::array<std::int64_t, 3> key{};
  const std::array<double, 3> shifted{(value.x - minimum.x) / tolerance,
                                      (value.y - minimum.y) / tolerance,
                                      (value.z - minimum.z) / tolerance};
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    if (!std::isfinite(shifted[axis]) ||
        shifted[axis] >
            static_cast<double>(std::numeric_limits<std::int64_t>::max() - 2) ||
        shifted[axis] <
            static_cast<double>(std::numeric_limits<std::int64_t>::min() + 2)) {
      fail("vertex welding bucket index is not representable");
    }
    key[axis] = static_cast<std::int64_t>(std::floor(shifted[axis]));
  }
  return key;
}

double coordinate(runtime::Real3 value, std::size_t axis) noexcept {
  if (axis == 0U) {
    return value.x;
  }
  if (axis == 1U) {
    return value.y;
  }
  return value.z;
}

struct Point2 final {
  double x{};
  double y{};
};

struct Point2Long final {
  long double x{};
  long double y{};
};

Point2 project(Real3 point, std::size_t dropped_axis) noexcept {
  if (dropped_axis == 0U) {
    return {point.y, point.z};
  }
  if (dropped_axis == 1U) {
    return {point.x, point.z};
  }
  return {point.x, point.y};
}

double orient2d(Point2 a, Point2 b, Point2 c) noexcept {
  return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

bool point_in_triangle_2d(Point2 point, const std::array<Point2, 3> &triangle,
                          double tolerance) noexcept {
  const double o0 = orient2d(triangle[0], triangle[1], point);
  const double o1 = orient2d(triangle[1], triangle[2], point);
  const double o2 = orient2d(triangle[2], triangle[0], point);
  const bool nonnegative =
      o0 >= -tolerance && o1 >= -tolerance && o2 >= -tolerance;
  const bool nonpositive =
      o0 <= tolerance && o1 <= tolerance && o2 <= tolerance;
  return nonnegative || nonpositive;
}

bool segments_intersect_2d(Point2 a, Point2 b, Point2 c, Point2 d,
                           double tolerance) noexcept {
  const double o1 = orient2d(a, b, c);
  const double o2 = orient2d(a, b, d);
  const double o3 = orient2d(c, d, a);
  const double o4 = orient2d(c, d, b);
  const auto opposite = [tolerance](double lhs, double rhs) {
    return (lhs > tolerance && rhs < -tolerance) ||
           (lhs < -tolerance && rhs > tolerance);
  };
  if (opposite(o1, o2) && opposite(o3, o4)) {
    return true;
  }
  const auto within = [tolerance](double value, double first, double second) {
    return value >= std::min(first, second) - tolerance &&
           value <= std::max(first, second) + tolerance;
  };
  const auto on_segment = [&](Point2 p, Point2 q, Point2 r,
                              double orientation) {
    return std::abs(orientation) <= tolerance && within(r.x, p.x, q.x) &&
           within(r.y, p.y, q.y);
  };
  return on_segment(a, b, c, o1) || on_segment(a, b, d, o2) ||
         on_segment(c, d, a, o3) || on_segment(c, d, b, o4);
}

bool point_on_segment_2d(Point2 point, Point2 first, Point2 second,
                         double tolerance) noexcept {
  return std::abs(orient2d(first, second, point)) <= tolerance &&
         point.x >= std::min(first.x, second.x) - tolerance &&
         point.x <= std::max(first.x, second.x) + tolerance &&
         point.y >= std::min(first.y, second.y) - tolerance &&
         point.y <= std::max(first.y, second.y) + tolerance;
}

struct SharedVertices final {
  std::array<std::uint64_t, 3> ids{};
  std::size_t count{};

  bool contains(std::uint64_t id) const noexcept {
    return std::find(ids.begin(),
                     ids.begin() + static_cast<std::ptrdiff_t>(count),
                     id) != ids.begin() + static_cast<std::ptrdiff_t>(count);
  }
};

SharedVertices shared_vertices(const IndexedTriangle &a,
                               const IndexedTriangle &b) noexcept {
  SharedVertices shared;
  for (const std::uint64_t a_vertex : a.vertices) {
    if (std::find(b.vertices.begin(), b.vertices.end(), a_vertex) !=
        b.vertices.end()) {
      shared.ids[shared.count++] = a_vertex;
    }
  }
  return shared;
}

double triangle_edge_scale(const SurfaceTriangle &triangle) noexcept {
  double scale = 0.0;
  for (std::size_t edge = 0U; edge < 3U; ++edge) {
    scale = std::max(scale, norm(subtract(triangle.vertices_m[(edge + 1U) % 3U],
                                          triangle.vertices_m[edge])));
  }
  return scale;
}

std::size_t projection_axis(Real3 normal) noexcept {
  std::size_t dropped = 0U;
  if (std::abs(normal.y) > std::abs(normal.x)) {
    dropped = 1U;
  }
  if (std::abs(coordinate(normal, 2U)) >
      std::abs(coordinate(normal, dropped))) {
    dropped = 2U;
  }
  return dropped;
}

bool coplanar_segment_overlaps_triangle(Real3 first, Real3 second,
                                        const SurfaceTriangle &triangle,
                                        double coincidence) noexcept {
  const std::size_t dropped =
      projection_axis(triangle.geometric_outward_normal);
  const Point2 projected_first = project(first, dropped);
  const Point2 projected_second = project(second, dropped);
  std::array<Point2, 3> projected_triangle{};
  for (std::size_t index = 0U; index < 3U; ++index) {
    projected_triangle[index] = project(triangle.vertices_m[index], dropped);
  }
  const double tolerance =
      coincidence *
      std::max(triangle_edge_scale(triangle), norm(subtract(second, first)));
  if (point_in_triangle_2d(projected_first, projected_triangle, tolerance) ||
      point_in_triangle_2d(projected_second, projected_triangle, tolerance)) {
    return true;
  }
  for (std::size_t edge = 0U; edge < 3U; ++edge) {
    if (segments_intersect_2d(
            projected_first, projected_second, projected_triangle[edge],
            projected_triangle[(edge + 1U) % 3U], tolerance)) {
      return true;
    }
  }
  return false;
}

bool same_point_bits(Real3 first, Real3 second) noexcept {
  return double_bits(first.x) == double_bits(second.x) &&
         double_bits(first.y) == double_bits(second.y) &&
         double_bits(first.z) == double_bits(second.z);
}

bool segment_provably_leaves_triangle_at_vertex(
    Real3 shared_vertex, Real3 other_endpoint,
    const SurfaceTriangle &triangle) noexcept {
  std::size_t shared_position = triangle.vertices_m.size();
  for (std::size_t position = 0U; position < triangle.vertices_m.size();
       ++position) {
    if (same_point_bits(shared_vertex, triangle.vertices_m[position])) {
      shared_position = position;
      break;
    }
  }
  if (shared_position == triangle.vertices_m.size()) {
    return false;
  }

  const std::size_t dropped =
      projection_axis(triangle.geometric_outward_normal);
  std::array<Point2, 3> projected_triangle{};
  for (std::size_t position = 0U; position < projected_triangle.size();
       ++position) {
    projected_triangle[position] =
        project(triangle.vertices_m[position], dropped);
  }
  const Point2 projected_shared = project(shared_vertex, dropped);
  const Point2 projected_other = project(other_endpoint, dropped);
  const Point2Long direction{
      static_cast<long double>(projected_other.x) -
          static_cast<long double>(projected_shared.x),
      static_cast<long double>(projected_other.y) -
          static_cast<long double>(projected_shared.y)};

  double dominant_component =
      coordinate(triangle.geometric_outward_normal, dropped);
  if (dropped == 1U) {
    dominant_component = -dominant_component;
  }
  if (dominant_component == 0.0) {
    return false;
  }
  const long double orientation_sign =
      dominant_component > 0.0 ? 1.0L : -1.0L;
  constexpr long double error_multiplier = 64.0L;
  constexpr long double epsilon =
      static_cast<long double>(std::numeric_limits<double>::epsilon());

  for (std::size_t edge = 0U; edge < projected_triangle.size(); ++edge) {
    const std::size_t next = (edge + 1U) % projected_triangle.size();
    if (edge != shared_position && next != shared_position) {
      continue;
    }
    const Point2Long projected_edge{
        static_cast<long double>(projected_triangle[next].x) -
            static_cast<long double>(projected_triangle[edge].x),
        static_cast<long double>(projected_triangle[next].y) -
            static_cast<long double>(projected_triangle[edge].y)};
    const long double first_product = projected_edge.x * direction.y;
    const long double second_product = projected_edge.y * direction.x;
    const long double derivative =
        orientation_sign * (first_product - second_product);
    const long double error_bound =
        error_multiplier * epsilon *
        (std::abs(first_product) + std::abs(second_product));
    if (derivative < -error_bound) {
      return true;
    }
  }
  return false;
}

bool coplanar_contact_is_forbidden(
    const SurfaceTriangle &target, Real3 first, Real3 second,
    bool first_is_canonical_shared_vertex,
    bool second_is_canonical_shared_vertex) noexcept {
  if (first_is_canonical_shared_vertex ==
      second_is_canonical_shared_vertex) {
    return true;
  }
  if (first_is_canonical_shared_vertex) {
    return !segment_provably_leaves_triangle_at_vertex(first, second, target);
  }
  return !segment_provably_leaves_triangle_at_vertex(second, first, target);
}

bool coplanar_triangles_have_forbidden_overlap(const SurfaceTriangle &a,
                                               const IndexedTriangle &a_indexed,
                                               const SurfaceTriangle &b,
                                               const IndexedTriangle &b_indexed,
                                               double coincidence) noexcept {
  const std::size_t dropped = projection_axis(a.geometric_outward_normal);
  std::array<Point2, 3> projected_a{};
  std::array<Point2, 3> projected_b{};
  for (std::size_t index = 0; index < 3U; ++index) {
    projected_a[index] = project(a.vertices_m[index], dropped);
    projected_b[index] = project(b.vertices_m[index], dropped);
  }
  const double tolerance =
      coincidence * std::max(triangle_edge_scale(a), triangle_edge_scale(b));
  const SharedVertices shared = shared_vertices(a_indexed, b_indexed);
  for (std::size_t vertex = 0U; vertex < 3U; ++vertex) {
    if (!shared.contains(a_indexed.vertices[vertex]) &&
        point_in_triangle_2d(projected_a[vertex], projected_b, tolerance)) {
      return true;
    }
    if (!shared.contains(b_indexed.vertices[vertex]) &&
        point_in_triangle_2d(projected_b[vertex], projected_a, tolerance)) {
      return true;
    }
  }
  for (std::size_t edge_a = 0; edge_a < 3U; ++edge_a) {
    for (std::size_t edge_b = 0; edge_b < 3U; ++edge_b) {
      const std::size_t a_next = (edge_a + 1U) % 3U;
      const std::size_t b_next = (edge_b + 1U) % 3U;
      if (!segments_intersect_2d(projected_a[edge_a], projected_a[a_next],
                                 projected_b[edge_b], projected_b[b_next],
                                 tolerance)) {
        continue;
      }
      const bool same_shared_edge =
          shared.count == 2U && shared.contains(a_indexed.vertices[edge_a]) &&
          shared.contains(a_indexed.vertices[a_next]) &&
          shared.contains(b_indexed.vertices[edge_b]) &&
          shared.contains(b_indexed.vertices[b_next]);
      if (same_shared_edge) {
        continue;
      }
      if ((!shared.contains(a_indexed.vertices[edge_a]) &&
           point_on_segment_2d(projected_a[edge_a], projected_b[edge_b],
                               projected_b[b_next], tolerance)) ||
          (!shared.contains(a_indexed.vertices[a_next]) &&
           point_on_segment_2d(projected_a[a_next], projected_b[edge_b],
                               projected_b[b_next], tolerance)) ||
          (!shared.contains(b_indexed.vertices[edge_b]) &&
           point_on_segment_2d(projected_b[edge_b], projected_a[edge_a],
                               projected_a[a_next], tolerance)) ||
          (!shared.contains(b_indexed.vertices[b_next]) &&
           point_on_segment_2d(projected_b[b_next], projected_a[edge_a],
                               projected_a[a_next], tolerance))) {
        return true;
      }
      const double a_first = orient2d(projected_a[edge_a], projected_a[a_next],
                                      projected_b[edge_b]);
      const double a_second = orient2d(projected_a[edge_a], projected_a[a_next],
                                       projected_b[b_next]);
      const double b_first = orient2d(projected_b[edge_b], projected_b[b_next],
                                      projected_a[edge_a]);
      const double b_second = orient2d(projected_b[edge_b], projected_b[b_next],
                                       projected_a[a_next]);
      if (((a_first > tolerance && a_second < -tolerance) ||
           (a_first < -tolerance && a_second > tolerance)) &&
          ((b_first > tolerance && b_second < -tolerance) ||
           (b_first < -tolerance && b_second > tolerance))) {
        return true;
      }
    }
  }
  return false;
}

bool segment_hits_triangle(Real3 first, Real3 second,
                           const SurfaceTriangle &triangle, double coincidence,
                           bool &coplanar, Real3 &intersection) noexcept {
  const Real3 direction = subtract(second, first);
  const Real3 edge1 = subtract(triangle.vertices_m[1], triangle.vertices_m[0]);
  const Real3 edge2 = subtract(triangle.vertices_m[2], triangle.vertices_m[0]);
  const double first_distance =
      std::abs(dot(triangle.geometric_outward_normal,
                   subtract(first, triangle.vertices_m[0])));
  const double second_distance =
      std::abs(dot(triangle.geometric_outward_normal,
                   subtract(second, triangle.vertices_m[0])));
  if (first_distance <= coincidence && second_distance <= coincidence) {
    coplanar = coplanar_segment_overlaps_triangle(first, second, triangle,
                                                  coincidence);
    return false;
  }
  const Real3 p = cross(direction, edge2);
  const double determinant = dot(edge1, p);
  const double determinant_scale = norm(direction) * norm(edge1) * norm(edge2);
  const double determinant_tolerance =
      64.0 * std::numeric_limits<double>::epsilon() * determinant_scale;
  if (std::abs(determinant) <= determinant_tolerance) {
    coplanar = first_distance <= coincidence &&
               second_distance <= coincidence &&
               coplanar_segment_overlaps_triangle(first, second, triangle,
                                                  coincidence);
    return false;
  }
  coplanar = false;
  const double inverse = 1.0 / determinant;
  const Real3 offset = subtract(first, triangle.vertices_m[0]);
  const double u = dot(offset, p) * inverse;
  const Real3 q = cross(offset, edge1);
  const double v = dot(direction, q) * inverse;
  const double fraction = dot(edge2, q) * inverse;
  const double barycentric_tolerance =
      coincidence /
      std::max({norm(edge1), norm(edge2), norm(direction), coincidence});
  const bool hit = u >= -barycentric_tolerance && v >= -barycentric_tolerance &&
                   u + v <= 1.0 + barycentric_tolerance &&
                   fraction >= -barycentric_tolerance &&
                   fraction <= 1.0 + barycentric_tolerance;
  if (hit) {
    intersection =
        add(first, multiply(direction, std::clamp(fraction, 0.0, 1.0)));
  }
  return hit;
}

bool boxes_overlap(const SurfaceTriangle &a, const SurfaceTriangle &b,
                   double tolerance) noexcept {
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    double a_min = coordinate(a.vertices_m[0], axis);
    double a_max = a_min;
    double b_min = coordinate(b.vertices_m[0], axis);
    double b_max = b_min;
    for (std::size_t vertex = 1U; vertex < 3U; ++vertex) {
      a_min = std::min(a_min, coordinate(a.vertices_m[vertex], axis));
      a_max = std::max(a_max, coordinate(a.vertices_m[vertex], axis));
      b_min = std::min(b_min, coordinate(b.vertices_m[vertex], axis));
      b_max = std::max(b_max, coordinate(b.vertices_m[vertex], axis));
    }
    if (a_max < b_min - tolerance || b_max < a_min - tolerance) {
      return false;
    }
  }
  return true;
}

bool triangles_have_forbidden_intersection(const SurfaceTriangle &a,
                                           const IndexedTriangle &a_indexed,
                                           const SurfaceTriangle &b,
                                           const IndexedTriangle &b_indexed,
                                           double coincidence) noexcept {
  if (!boxes_overlap(a, b, coincidence)) {
    return false;
  }
  const double plane_distance = std::abs(dot(
      a.geometric_outward_normal, subtract(b.vertices_m[0], a.vertices_m[0])));
  const double normal_alignment =
      std::abs(dot(a.geometric_outward_normal, b.geometric_outward_normal));
  if (plane_distance <= coincidence &&
      normal_alignment >=
          1.0 - 128.0 * std::numeric_limits<double>::epsilon()) {
    return coplanar_triangles_have_forbidden_overlap(a, a_indexed, b, b_indexed,
                                                     coincidence);
  }
  const SharedVertices shared = shared_vertices(a_indexed, b_indexed);
  const auto allowed_shared_contact = [&](Real3 point) {
    for (std::size_t index = 0U; index < shared.count; ++index) {
      const std::uint64_t vertex = shared.ids[index];
      for (std::size_t position = 0U; position < 3U; ++position) {
        if (a_indexed.vertices[position] == vertex &&
            norm(subtract(point, a.vertices_m[position])) <= coincidence) {
          return true;
        }
      }
    }
    return false;
  };
  const auto edge_is_shared = [&](const IndexedTriangle &triangle,
                                  std::size_t edge) {
    return shared.count == 2U && shared.contains(triangle.vertices[edge]) &&
           shared.contains(triangle.vertices[(edge + 1U) % 3U]);
  };
  for (std::size_t edge = 0; edge < 3U; ++edge) {
    bool coplanar = false;
    Real3 intersection{};
    if (segment_hits_triangle(a.vertices_m[edge],
                              a.vertices_m[(edge + 1U) % 3U], b, coincidence,
                              coplanar, intersection)) {
      if (!allowed_shared_contact(intersection)) {
        return true;
      }
    } else if (coplanar && !edge_is_shared(a_indexed, edge) &&
               coplanar_contact_is_forbidden(
                   b, a.vertices_m[edge],
                   a.vertices_m[(edge + 1U) % 3U],
                   shared.contains(a_indexed.vertices[edge]),
                   shared.contains(a_indexed.vertices[(edge + 1U) % 3U]))) {
      return true;
    }
    coplanar = false;
    if (segment_hits_triangle(b.vertices_m[edge],
                              b.vertices_m[(edge + 1U) % 3U], a, coincidence,
                              coplanar, intersection)) {
      if (!allowed_shared_contact(intersection)) {
        return true;
      }
    } else if (coplanar && !edge_is_shared(b_indexed, edge) &&
               coplanar_contact_is_forbidden(
                   a, b.vertices_m[edge],
                   b.vertices_m[(edge + 1U) % 3U],
                   shared.contains(b_indexed.vertices[edge]),
                   shared.contains(b_indexed.vertices[(edge + 1U) % 3U]))) {
      return true;
    }
  }
  return false;
}

void canonicalize_vertex_ids(SurfaceStorage &surface) {
  const auto absent = std::numeric_limits<std::uint64_t>::max();
  std::vector<std::uint64_t> remap(surface.vertices.size(), absent);
  std::vector<Real3> canonical_vertices;
  canonical_vertices.reserve(surface.vertices.size());
  for (IndexedTriangle &triangle : surface.indexed_triangles) {
    for (std::uint64_t &vertex : triangle.vertices) {
      if (remap[vertex] == absent) {
        remap[vertex] = static_cast<std::uint64_t>(canonical_vertices.size());
        canonical_vertices.push_back(surface.vertices[vertex]);
      }
      vertex = remap[vertex];
    }
  }
  surface.vertices = std::move(canonical_vertices);
}

void update_normalized_surface_scale(SurfaceStorage &surface) {
  Real3 minimum = surface.vertices.front();
  Real3 maximum = surface.vertices.front();
  for (const Real3 vertex : surface.vertices) {
    minimum.x = std::min(minimum.x, vertex.x);
    minimum.y = std::min(minimum.y, vertex.y);
    minimum.z = std::min(minimum.z, vertex.z);
    maximum.x = std::max(maximum.x, vertex.x);
    maximum.y = std::max(maximum.y, vertex.y);
    maximum.z = std::max(maximum.z, vertex.z);
  }
  const double reference_length = norm(subtract(maximum, minimum));
  const double epsilon = std::numeric_limits<double>::epsilon();
  const double weld_tolerance = 128.0 * epsilon * reference_length;
  const double minimum_area =
      1024.0 * epsilon * reference_length * reference_length;
  const double coincidence = 512.0 * epsilon * reference_length;
  if (!std::isfinite(reference_length) || reference_length <= 0.0 ||
      !std::isfinite(weld_tolerance) || weld_tolerance <= 0.0 ||
      !std::isfinite(minimum_area) || minimum_area <= 0.0 ||
      !std::isfinite(coincidence) || coincidence <= 0.0) {
    fail("normalized surface scale is not positive and finite");
  }
  surface.bounding_box_min = minimum;
  surface.bounding_box_max = maximum;
  surface.reference_length_m = reference_length;
  surface.weld_tolerance_m = weld_tolerance;
  surface.minimum_triangle_area_m2 = minimum_area;
  surface.intersection_coincidence_m = coincidence;
}

double oriented_volume(const SurfaceStorage &surface) noexcept {
  const Real3 reference = surface.bounding_box_min;
  double volume = 0.0;
  for (const IndexedTriangle &triangle : surface.indexed_triangles) {
    const Real3 a = subtract(surface.vertices[triangle.vertices[0]], reference);
    const Real3 b = subtract(surface.vertices[triangle.vertices[1]], reference);
    const Real3 c = subtract(surface.vertices[triangle.vertices[2]], reference);
    volume += dot(a, cross(b, c)) / 6.0;
  }
  return volume;
}

std::uint64_t surface_fingerprint(const SurfaceStorage &surface) noexcept {
  std::uint64_t hash = fnv_offset;
  hash_u64(hash, UINT64_C(0x48554e44554e5333));
  hash_u64(hash, surface.vertices.size());
  hash_u64(hash, surface.triangles.size());
  hash_u64(hash, double_bits(surface.bounding_box_min.x));
  hash_u64(hash, double_bits(surface.bounding_box_min.y));
  hash_u64(hash, double_bits(surface.bounding_box_min.z));
  hash_u64(hash, double_bits(surface.bounding_box_max.x));
  hash_u64(hash, double_bits(surface.bounding_box_max.y));
  hash_u64(hash, double_bits(surface.bounding_box_max.z));
  hash_u64(hash, double_bits(surface.reference_length_m));
  hash_u64(hash, double_bits(surface.weld_tolerance_m));
  hash_u64(hash, double_bits(surface.minimum_triangle_area_m2));
  hash_u64(hash, double_bits(surface.intersection_coincidence_m));
  for (const Real3 vertex : surface.vertices) {
    hash_u64(hash, double_bits(vertex.x));
    hash_u64(hash, double_bits(vertex.y));
    hash_u64(hash, double_bits(vertex.z));
  }
  for (const IndexedTriangle &triangle : surface.indexed_triangles) {
    hash_u64(hash, triangle.id);
    for (const std::uint64_t vertex : triangle.vertices) {
      hash_u64(hash, vertex);
    }
  }
  hash_u64(hash, double_bits(surface.closed_volume_m3));
  return hash;
}

std::shared_ptr<const SurfaceStorage>
normalize_triangles(std::vector<RawTriangle> raw, double length_scale_to_m) {
  if (!std::isfinite(length_scale_to_m) || length_scale_to_m <= 0.0) {
    fail("length scale must be positive and finite");
  }
  if (raw.empty()) {
    fail("STL contains no triangles");
  }
  if (raw.size() >
      static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max())) {
    fail("STL triangle count exceeds 2^31-1");
  }

  Real3 minimum{std::numeric_limits<double>::infinity(),
                std::numeric_limits<double>::infinity(),
                std::numeric_limits<double>::infinity()};
  Real3 maximum{-std::numeric_limits<double>::infinity(),
                -std::numeric_limits<double>::infinity(),
                -std::numeric_limits<double>::infinity()};
  for (RawTriangle &triangle : raw) {
    if (!finite(triangle.file_normal)) {
      fail("STL contains a non-finite file normal");
    }
    for (Real3 &vertex : triangle.vertices) {
      if (!finite(vertex)) {
        fail("STL contains a non-finite coordinate");
      }
      vertex = multiply(vertex, length_scale_to_m);
      if (!finite(vertex)) {
        fail("scaled STL coordinate is non-finite");
      }
      minimum.x = std::min(minimum.x, vertex.x);
      minimum.y = std::min(minimum.y, vertex.y);
      minimum.z = std::min(minimum.z, vertex.z);
      maximum.x = std::max(maximum.x, vertex.x);
      maximum.y = std::max(maximum.y, vertex.y);
      maximum.z = std::max(maximum.z, vertex.z);
    }
  }
  const double reference_length = norm(subtract(maximum, minimum));
  if (!std::isfinite(reference_length) || reference_length <= 0.0) {
    fail("surface bounding-box diagonal must be positive and finite");
  }
  const double epsilon = std::numeric_limits<double>::epsilon();
  const double weld_tolerance = 128.0 * epsilon * reference_length;
  const double minimum_area =
      1024.0 * epsilon * reference_length * reference_length;
  const double coincidence = 512.0 * epsilon * reference_length;
  if (!std::isfinite(weld_tolerance) || weld_tolerance <= 0.0 ||
      !std::isfinite(minimum_area) || minimum_area <= 0.0 ||
      !std::isfinite(coincidence) || coincidence <= 0.0) {
    fail("derived surface tolerances are not positive and finite");
  }

  auto surface = std::make_shared<SurfaceStorage>();
  surface->bounding_box_min = minimum;
  surface->bounding_box_max = maximum;
  surface->reference_length_m = reference_length;
  surface->weld_tolerance_m = weld_tolerance;
  surface->minimum_triangle_area_m2 = minimum_area;
  surface->intersection_coincidence_m = coincidence;

  using Bucket = std::array<std::int64_t, 3>;
  std::map<Bucket, std::vector<std::uint64_t>> buckets;
  const auto weld_vertex = [&](Real3 vertex) {
    const Bucket key = bucket_key(vertex, minimum, weld_tolerance);
    std::uint64_t match = std::numeric_limits<std::uint64_t>::max();
    for (std::int64_t dx = -1; dx <= 1; ++dx) {
      for (std::int64_t dy = -1; dy <= 1; ++dy) {
        for (std::int64_t dz = -1; dz <= 1; ++dz) {
          const Bucket candidate_key{key[0] + dx, key[1] + dy, key[2] + dz};
          const auto found = buckets.find(candidate_key);
          if (found == buckets.end()) {
            continue;
          }
          for (const std::uint64_t candidate : found->second) {
            if (norm(subtract(vertex, surface->vertices[candidate])) <=
                    weld_tolerance &&
                candidate < match) {
              match = candidate;
            }
          }
        }
      }
    }
    if (match != std::numeric_limits<std::uint64_t>::max()) {
      return match;
    }
    if (surface->vertices.size() >= std::numeric_limits<std::uint64_t>::max()) {
      fail("surface vertex count exceeds uint64 range");
    }
    const auto id = static_cast<std::uint64_t>(surface->vertices.size());
    surface->vertices.push_back(vertex);
    buckets[key].push_back(id);
    return id;
  };

  surface->indexed_triangles.reserve(raw.size());
  std::map<std::array<std::uint64_t, 3>, TriangleId> triangle_keys;
  for (std::size_t index = 0; index < raw.size(); ++index) {
    IndexedTriangle triangle;
    triangle.id = static_cast<TriangleId>(index);
    for (std::size_t vertex = 0; vertex < 3U; ++vertex) {
      triangle.vertices[vertex] = weld_vertex(raw[index].vertices[vertex]);
    }
    if (triangle.vertices[0] == triangle.vertices[1] ||
        triangle.vertices[1] == triangle.vertices[2] ||
        triangle.vertices[2] == triangle.vertices[0]) {
      fail("triangle is degenerate after deterministic welding");
    }
    std::array<std::uint64_t, 3> key = triangle.vertices;
    std::sort(key.begin(), key.end());
    if (!triangle_keys.emplace(key, triangle.id).second) {
      fail("surface contains a duplicate triangle");
    }
    const Real3 edge1 = subtract(surface->vertices[triangle.vertices[1]],
                                 surface->vertices[triangle.vertices[0]]);
    const Real3 edge2 = subtract(surface->vertices[triangle.vertices[2]],
                                 surface->vertices[triangle.vertices[0]]);
    const double area = 0.5 * norm(cross(edge1, edge2));
    if (!std::isfinite(area) || area <= minimum_area) {
      fail("triangle area is below the derived minimum");
    }
    surface->indexed_triangles.push_back(triangle);
  }

  std::map<std::array<std::uint64_t, 2>, std::vector<EdgeOccurrence>> edges;
  for (const IndexedTriangle &triangle : surface->indexed_triangles) {
    for (std::size_t edge = 0; edge < 3U; ++edge) {
      const std::uint64_t from = triangle.vertices[edge];
      const std::uint64_t to = triangle.vertices[(edge + 1U) % 3U];
      edges[{std::min(from, to), std::max(from, to)}].push_back(
          {from, to, triangle.id});
    }
  }

  std::vector<std::vector<TriangleId>> adjacency(raw.size());
  for (const auto &entry : edges) {
    const auto &occurrences = entry.second;
    if (occurrences.size() != 2U) {
      fail(occurrences.size() < 2U ? "surface contains an open edge"
                                   : "surface contains a non-manifold edge");
    }
    if (occurrences[0].from != occurrences[1].to ||
        occurrences[0].to != occurrences[1].from) {
      fail("adjacent triangle orientation is inconsistent");
    }
    adjacency[occurrences[0].triangle].push_back(occurrences[1].triangle);
    adjacency[occurrences[1].triangle].push_back(occurrences[0].triangle);
  }

  std::vector<bool> visited(raw.size(), false);
  std::queue<TriangleId> pending;
  pending.push(0U);
  visited[0] = true;
  std::size_t visited_count = 0U;
  while (!pending.empty()) {
    const TriangleId current = pending.front();
    pending.pop();
    ++visited_count;
    for (const TriangleId next : adjacency[current]) {
      if (!visited[next]) {
        visited[next] = true;
        pending.push(next);
      }
    }
  }
  if (visited_count != raw.size()) {
    fail("surface contains multiple connected components");
  }

  double signed_volume = oriented_volume(*surface);
  const double minimum_volume = minimum_area * reference_length;
  if (!std::isfinite(signed_volume) ||
      std::abs(signed_volume) <= minimum_volume) {
    fail("surface has zero enclosed volume");
  }
  if (signed_volume < 0.0) {
    for (IndexedTriangle &triangle : surface->indexed_triangles) {
      std::swap(triangle.vertices[1], triangle.vertices[2]);
    }
  }
  canonicalize_vertex_ids(*surface);
  update_normalized_surface_scale(*surface);
  signed_volume = oriented_volume(*surface);
  const double normalized_minimum_volume =
      surface->minimum_triangle_area_m2 * surface->reference_length_m;
  if (!std::isfinite(signed_volume) ||
      signed_volume <= normalized_minimum_volume) {
    fail("surface has zero enclosed volume");
  }
  surface->closed_volume_m3 = signed_volume;

  surface->triangles.reserve(raw.size());
  for (const IndexedTriangle &indexed : surface->indexed_triangles) {
    SurfaceTriangle triangle;
    triangle.id = indexed.id;
    for (std::size_t vertex = 0; vertex < 3U; ++vertex) {
      triangle.vertices_m[vertex] = surface->vertices[indexed.vertices[vertex]];
    }
    const Real3 normal =
        cross(subtract(triangle.vertices_m[1], triangle.vertices_m[0]),
              subtract(triangle.vertices_m[2], triangle.vertices_m[0]));
    const double twice_area = norm(normal);
    triangle.area_m2 = 0.5 * twice_area;
    triangle.geometric_outward_normal = multiply(normal, 1.0 / twice_area);
    surface->triangles.push_back(triangle);
  }

  const auto intersection_query = build_surface_query(surface);
  for (std::size_t first = 0; first < surface->triangles.size(); ++first) {
    const SurfaceTriangle &triangle = surface->triangles[first];
    Real3 triangle_minimum = triangle.vertices_m[0];
    Real3 triangle_maximum = triangle.vertices_m[0];
    for (std::size_t vertex = 1U; vertex < 3U; ++vertex) {
      triangle_minimum.x =
          std::min(triangle_minimum.x, triangle.vertices_m[vertex].x);
      triangle_minimum.y =
          std::min(triangle_minimum.y, triangle.vertices_m[vertex].y);
      triangle_minimum.z =
          std::min(triangle_minimum.z, triangle.vertices_m[vertex].z);
      triangle_maximum.x =
          std::max(triangle_maximum.x, triangle.vertices_m[vertex].x);
      triangle_maximum.y =
          std::max(triangle_maximum.y, triangle.vertices_m[vertex].y);
      triangle_maximum.z =
          std::max(triangle_maximum.z, triangle.vertices_m[vertex].z);
    }
    triangle_minimum = {triangle_minimum.x - coincidence,
                        triangle_minimum.y - coincidence,
                        triangle_minimum.z - coincidence};
    triangle_maximum = {triangle_maximum.x + coincidence,
                        triangle_maximum.y + coincidence,
                        triangle_maximum.z + coincidence};
    const auto candidates = bounded_candidates(
        *intersection_query, triangle_minimum, triangle_maximum);
    for (const TriangleId second : candidates) {
      if (second <= first) {
        continue;
      }
      if (triangles_have_forbidden_intersection(
              surface->triangles[first], surface->indexed_triangles[first],
              surface->triangles[second], surface->indexed_triangles[second],
              coincidence)) {
        fail("surface contains a self-intersection");
      }
    }
  }

  surface->fingerprint = surface_fingerprint(*surface);
  return surface;
}

void append_u32(std::vector<std::uint8_t> &bytes, std::uint32_t value) {
  for (unsigned shift = 0; shift < 32U; shift += 8U) {
    bytes.push_back(static_cast<std::uint8_t>((value >> shift) & 0xffU));
  }
}

void append_u64(std::vector<std::uint8_t> &bytes, std::uint64_t value) {
  for (unsigned shift = 0; shift < 64U; shift += 8U) {
    bytes.push_back(static_cast<std::uint8_t>((value >> shift) & 0xffU));
  }
}

std::uint64_t consume_u64(const std::vector<std::uint8_t> &bytes,
                          std::size_t &offset) {
  if (offset > bytes.size() || bytes.size() - offset < 8U) {
    fail("normalized surface payload is truncated");
  }
  std::uint64_t value = 0;
  for (unsigned index = 0; index < 8U; ++index) {
    value |= static_cast<std::uint64_t>(bytes[offset + index]) << (8U * index);
  }
  offset += 8U;
  return value;
}

std::uint32_t consume_u32(const std::vector<std::uint8_t> &bytes,
                          std::size_t &offset) {
  if (offset > bytes.size() || bytes.size() - offset < 4U) {
    fail("normalized surface payload is truncated");
  }
  const std::uint32_t value = read_u32(bytes, offset);
  offset += 4U;
  return value;
}

std::vector<std::uint8_t> read_file(const std::filesystem::path &path) {
  std::ifstream stream(path, std::ios::binary | std::ios::ate);
  if (!stream) {
    fail("unable to open STL file");
  }
  const std::streampos end = stream.tellg();
  if (end < 0) {
    fail("unable to determine STL file size");
  }
  const auto unsigned_end = static_cast<std::uintmax_t>(end);
  if (unsigned_end > std::numeric_limits<std::size_t>::max() ||
      unsigned_end > static_cast<std::uintmax_t>(
                         std::numeric_limits<std::streamsize>::max())) {
    fail("STL file exceeds local size range");
  }
  std::vector<std::uint8_t> bytes(static_cast<std::size_t>(unsigned_end));
  stream.seekg(0, std::ios::beg);
  if (!bytes.empty()) {
    stream.read(reinterpret_cast<char *>(bytes.data()),
                static_cast<std::streamsize>(bytes.size()));
  }
  if (!stream || stream.peek() != std::char_traits<char>::eof()) {
    fail("unable to read exact STL file bytes");
  }
  return bytes;
}

std::string mpi_error_message(int code, std::string_view operation) {
  std::array<char, MPI_MAX_ERROR_STRING> buffer{};
  int length = 0;
  if (MPI_Error_string(code, buffer.data(), &length) != MPI_SUCCESS ||
      length < 0 || static_cast<std::size_t>(length) > buffer.size()) {
    return std::string{operation} + " failed";
  }
  return std::string{operation} + " failed: " +
         std::string{buffer.data(), static_cast<std::size_t>(length)};
}

void check_mpi(int code, std::string_view operation) {
  if (code != MPI_SUCCESS) {
    throw runtime::Error(mpi_error_message(code, operation));
  }
}

void require_collective(const runtime::MpiContext &mpi, bool local_ok,
                        std::string_view message) {
  const runtime::CollectiveStatus status =
      runtime::collective_status(mpi, local_ok, message);
  if (!status.ok) {
    throw runtime::Error(status.message + " (lowest failing rank " +
                         std::to_string(status.failing_rank) + ")");
  }
}

void broadcast_bytes(std::vector<std::uint8_t> &bytes,
                     const runtime::MpiContext &mpi, int root,
                     std::size_t maximum_chunk_bytes) {
  std::uint64_t size = 0U;
  if (mpi.rank() == root) {
    size = bytes.size();
  }
  check_mpi(MPI_Bcast(&size, 1, MPI_UINT64_T, root, mpi.comm()),
            "MPI_Bcast surface byte length");
  const bool size_ok = size <= std::numeric_limits<std::size_t>::max();
  require_collective(mpi, size_ok,
                     "immersed surface: payload exceeds local size range");

  bool allocation_ok = true;
  if (mpi.rank() != root) {
    try {
      bytes.resize(static_cast<std::size_t>(size));
    } catch (...) {
      allocation_ok = false;
    }
  }
  require_collective(mpi, allocation_ok,
                     "immersed surface: unable to allocate payload");
  std::size_t offset = 0U;
  while (offset < bytes.size()) {
    const std::size_t chunk =
        std::min(maximum_chunk_bytes, bytes.size() - offset);
    check_mpi(MPI_Bcast(bytes.data() + offset, static_cast<int>(chunk),
                        MPI_BYTE, root, mpi.comm()),
              "MPI_Bcast surface bytes");
    offset += chunk;
  }
}

} // namespace

#if defined(HUNDUN_IMMERSED_ENABLE_TEST_ACCESS)
bool coplanar_contact_is_forbidden_for_test(
    const SurfaceTriangle &target, runtime::Real3 first, runtime::Real3 second,
    bool first_is_canonical_shared_vertex,
    bool second_is_canonical_shared_vertex) noexcept {
  return coplanar_contact_is_forbidden(
      target, first, second, first_is_canonical_shared_vertex,
      second_is_canonical_shared_vertex);
}
#endif

runtime::Real3 add(runtime::Real3 a, runtime::Real3 b) noexcept {
  return {a.x + b.x, a.y + b.y, a.z + b.z};
}

runtime::Real3 subtract(runtime::Real3 a, runtime::Real3 b) noexcept {
  return {a.x - b.x, a.y - b.y, a.z - b.z};
}

runtime::Real3 multiply(runtime::Real3 value, double factor) noexcept {
  return {value.x * factor, value.y * factor, value.z * factor};
}

double dot(runtime::Real3 a, runtime::Real3 b) noexcept {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

runtime::Real3 cross(runtime::Real3 a, runtime::Real3 b) noexcept {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

double squared_norm(runtime::Real3 value) noexcept { return dot(value, value); }

double norm(runtime::Real3 value) noexcept {
  return std::hypot(value.x, value.y, value.z);
}

bool finite(runtime::Real3 value) noexcept {
  return std::isfinite(value.x) && std::isfinite(value.y) &&
         std::isfinite(value.z);
}

std::shared_ptr<const SurfaceStorage>
read_and_normalize_stl(const std::vector<std::uint8_t> &bytes,
                       double length_scale_to_m) {
  return normalize_triangles(parse_stl(bytes), length_scale_to_m);
}

std::vector<std::uint8_t>
encode_normalized_surface(const SurfaceStorage &surface) {
  constexpr std::size_t fixed_bytes = 124U;
  constexpr std::size_t vertex_bytes = 24U;
  constexpr std::size_t triangle_bytes = 32U;
  std::vector<std::uint8_t> bytes;
  if (surface.vertices.size() >
      (std::numeric_limits<std::size_t>::max() - fixed_bytes) / vertex_bytes) {
    fail("normalized surface payload size overflows");
  }
  std::size_t payload_bytes =
      fixed_bytes + surface.vertices.size() * vertex_bytes;
  if (surface.indexed_triangles.size() >
      (std::numeric_limits<std::size_t>::max() - payload_bytes) /
          triangle_bytes) {
    fail("normalized surface payload size overflows");
  }
  payload_bytes += surface.indexed_triangles.size() * triangle_bytes;
  bytes.reserve(payload_bytes);
  bytes.insert(bytes.end(), payload_magic.begin(), payload_magic.end());
  append_u32(bytes, payload_version);
  append_u64(bytes, surface.vertices.size());
  append_u64(bytes, surface.indexed_triangles.size());
  for (const Real3 value :
       {surface.bounding_box_min, surface.bounding_box_max}) {
    append_u64(bytes, double_bits(value.x));
    append_u64(bytes, double_bits(value.y));
    append_u64(bytes, double_bits(value.z));
  }
  for (const double value :
       {surface.reference_length_m, surface.weld_tolerance_m,
        surface.minimum_triangle_area_m2, surface.intersection_coincidence_m,
        surface.closed_volume_m3}) {
    append_u64(bytes, double_bits(value));
  }
  for (const Real3 vertex : surface.vertices) {
    append_u64(bytes, double_bits(vertex.x));
    append_u64(bytes, double_bits(vertex.y));
    append_u64(bytes, double_bits(vertex.z));
  }
  for (const IndexedTriangle &triangle : surface.indexed_triangles) {
    append_u64(bytes, triangle.id);
    for (const std::uint64_t vertex : triangle.vertices) {
      append_u64(bytes, vertex);
    }
  }
  append_u64(bytes, surface.fingerprint);
  return bytes;
}

std::shared_ptr<const SurfaceStorage>
decode_normalized_surface(const std::vector<std::uint8_t> &bytes) {
  constexpr std::size_t fixed_bytes = 124U;
  constexpr std::size_t vertex_bytes = 24U;
  constexpr std::size_t triangle_bytes = 32U;
  if (bytes.size() < fixed_bytes ||
      !std::equal(payload_magic.begin(), payload_magic.end(), bytes.begin())) {
    fail("normalized surface payload has invalid magic");
  }
  std::size_t offset = payload_magic.size();
  if (consume_u32(bytes, offset) != payload_version) {
    fail("normalized surface payload has unsupported version");
  }
  const std::uint64_t vertex_count = consume_u64(bytes, offset);
  const std::uint64_t triangle_count = consume_u64(bytes, offset);
  if (vertex_count == 0U || triangle_count == 0U ||
      vertex_count > std::numeric_limits<std::size_t>::max() ||
      triangle_count > static_cast<std::uint64_t>(
                           std::numeric_limits<std::int32_t>::max()) ||
      triangle_count > std::numeric_limits<std::size_t>::max()) {
    fail("normalized surface counts are not representable");
  }
  const std::size_t vertices = static_cast<std::size_t>(vertex_count);
  const std::size_t count = static_cast<std::size_t>(triangle_count);
  if (vertices >
      (std::numeric_limits<std::size_t>::max() - fixed_bytes) / vertex_bytes) {
    fail("normalized surface payload size overflows");
  }
  std::size_t expected = fixed_bytes + vertices * vertex_bytes;
  if (count >
      (std::numeric_limits<std::size_t>::max() - expected) / triangle_bytes) {
    fail("normalized surface payload size overflows");
  }
  expected += count * triangle_bytes;
  if (bytes.size() != expected) {
    fail("normalized surface payload has invalid exact size");
  }

  auto surface = std::make_shared<SurfaceStorage>();
  const auto consume_real = [&] {
    return double_from_bits(consume_u64(bytes, offset));
  };
  surface->bounding_box_min = {consume_real(), consume_real(), consume_real()};
  surface->bounding_box_max = {consume_real(), consume_real(), consume_real()};
  surface->reference_length_m = consume_real();
  surface->weld_tolerance_m = consume_real();
  surface->minimum_triangle_area_m2 = consume_real();
  surface->intersection_coincidence_m = consume_real();
  surface->closed_volume_m3 = consume_real();
  if (!finite(surface->bounding_box_min) ||
      !finite(surface->bounding_box_max) ||
      surface->bounding_box_min.x > surface->bounding_box_max.x ||
      surface->bounding_box_min.y > surface->bounding_box_max.y ||
      surface->bounding_box_min.z > surface->bounding_box_max.z ||
      !std::isfinite(surface->reference_length_m) ||
      surface->reference_length_m <= 0.0 ||
      !std::isfinite(surface->weld_tolerance_m) ||
      surface->weld_tolerance_m <= 0.0 ||
      !std::isfinite(surface->minimum_triangle_area_m2) ||
      surface->minimum_triangle_area_m2 <= 0.0 ||
      !std::isfinite(surface->intersection_coincidence_m) ||
      surface->intersection_coincidence_m <= 0.0 ||
      !std::isfinite(surface->closed_volume_m3) ||
      surface->closed_volume_m3 <= 0.0) {
    fail("normalized surface payload has invalid metadata");
  }

  surface->vertices.resize(vertices);
  for (Real3 &vertex : surface->vertices) {
    vertex = {consume_real(), consume_real(), consume_real()};
    if (!finite(vertex)) {
      fail("normalized surface payload has a non-finite vertex");
    }
  }

  surface->indexed_triangles.resize(count);
  surface->triangles.reserve(count);
  for (std::size_t index = 0U; index < count; ++index) {
    IndexedTriangle &indexed = surface->indexed_triangles[index];
    indexed.id = consume_u64(bytes, offset);
    for (std::uint64_t &vertex : indexed.vertices) {
      vertex = consume_u64(bytes, offset);
    }
    if (indexed.id != index ||
        std::any_of(indexed.vertices.begin(), indexed.vertices.end(),
                    [vertex_count](std::uint64_t vertex) {
                      return vertex >= vertex_count;
                    }) ||
        indexed.vertices[0] == indexed.vertices[1] ||
        indexed.vertices[1] == indexed.vertices[2] ||
        indexed.vertices[2] == indexed.vertices[0]) {
      fail("normalized surface payload has invalid triangle indices");
    }
    SurfaceTriangle triangle;
    triangle.id = indexed.id;
    for (std::size_t vertex = 0U; vertex < 3U; ++vertex) {
      triangle.vertices_m[vertex] = surface->vertices[indexed.vertices[vertex]];
    }
    const Real3 normal =
        cross(subtract(triangle.vertices_m[1], triangle.vertices_m[0]),
              subtract(triangle.vertices_m[2], triangle.vertices_m[0]));
    const double twice_area = norm(normal);
    triangle.area_m2 = 0.5 * twice_area;
    if (!std::isfinite(triangle.area_m2) ||
        triangle.area_m2 <= surface->minimum_triangle_area_m2) {
      fail("normalized surface payload has an invalid triangle");
    }
    triangle.geometric_outward_normal = multiply(normal, 1.0 / twice_area);
    surface->triangles.push_back(triangle);
  }

  const std::uint64_t expected_fingerprint = consume_u64(bytes, offset);
  const double signed_volume = oriented_volume(*surface);
  if (!std::isfinite(signed_volume) || signed_volume <= 0.0 ||
      double_bits(signed_volume) != double_bits(surface->closed_volume_m3)) {
    fail("normalized surface payload has invalid oriented volume");
  }
  surface->fingerprint = surface_fingerprint(*surface);
  if (surface->fingerprint != expected_fingerprint) {
    fail("normalized surface payload fingerprint mismatch");
  }
  return surface;
}

} // namespace hundun::immersed::detail

namespace hundun::immersed {

ImmersedSurface
ImmersedSurface::load_collective(const std::filesystem::path &path,
                                 double length_scale_to_m,
                                 const runtime::MpiContext &mpi, int root) {
  return load_collective_impl(
      path, length_scale_to_m, mpi, root,
      static_cast<std::size_t>(std::numeric_limits<int>::max()));
}

ImmersedSurface ImmersedSurface::load_collective_impl(
    const std::filesystem::path &path, double length_scale_to_m,
    const runtime::MpiContext &mpi, int root, std::size_t maximum_chunk_bytes) {
  if (mpi.comm() == MPI_COMM_NULL) {
    throw runtime::Error("immersed surface: MPI context is invalid");
  }

  int minimum_root = root;
  int maximum_root = root;
  detail::check_mpi(
      MPI_Allreduce(&root, &minimum_root, 1, MPI_INT, MPI_MIN, mpi.comm()),
      "MPI_Allreduce surface root minimum");
  detail::check_mpi(
      MPI_Allreduce(&root, &maximum_root, 1, MPI_INT, MPI_MAX, mpi.comm()),
      "MPI_Allreduce surface root maximum");
  const bool root_ok =
      minimum_root == maximum_root && root >= 0 && root < mpi.size();
  detail::require_collective(
      mpi, root_ok,
      "immersed surface: collective root is invalid or inconsistent");
  root = minimum_root;

  std::uint64_t minimum_chunk = maximum_chunk_bytes;
  std::uint64_t maximum_chunk = maximum_chunk_bytes;
  detail::check_mpi(MPI_Allreduce(MPI_IN_PLACE, &minimum_chunk, 1, MPI_UINT64_T,
                                  MPI_MIN, mpi.comm()),
                    "MPI_Allreduce surface chunk minimum");
  detail::check_mpi(MPI_Allreduce(MPI_IN_PLACE, &maximum_chunk, 1, MPI_UINT64_T,
                                  MPI_MAX, mpi.comm()),
                    "MPI_Allreduce surface chunk maximum");
  const bool chunk_ok =
      minimum_chunk == maximum_chunk && maximum_chunk_bytes > 0U &&
      maximum_chunk_bytes <=
          static_cast<std::size_t>(std::numeric_limits<int>::max());
  detail::require_collective(
      mpi, chunk_ok,
      "immersed surface: collective chunk limit is invalid or inconsistent");

  std::uint64_t root_scale_bits = 0U;
  if (mpi.rank() == root) {
    root_scale_bits = detail::double_bits(length_scale_to_m);
  }
  detail::check_mpi(
      MPI_Bcast(&root_scale_bits, 1, MPI_UINT64_T, root, mpi.comm()),
      "MPI_Bcast surface length scale");
  detail::require_collective(
      mpi, detail::double_bits(length_scale_to_m) == root_scale_bits,
      "immersed surface: length scale differs between ranks");

  std::vector<std::uint8_t> root_path;
  const std::string local_path = path.lexically_normal().generic_string();
  if (mpi.rank() == root) {
    root_path.assign(local_path.begin(), local_path.end());
  }
  detail::broadcast_bytes(root_path, mpi, root, maximum_chunk_bytes);
  detail::require_collective(
      mpi, std::string(root_path.begin(), root_path.end()) == local_path,
      "immersed surface: STL path differs between ranks");

  std::vector<std::uint8_t> payload;
  bool root_ok_parse = true;
  std::string root_error;
  if (mpi.rank() == root) {
    try {
      const auto surface = detail::read_and_normalize_stl(
          detail::read_file(path), length_scale_to_m);
      payload = detail::encode_normalized_surface(*surface);
    } catch (const std::exception &error) {
      root_ok_parse = false;
      root_error = error.what();
    } catch (...) {
      root_ok_parse = false;
      root_error = "immersed surface: unknown STL load failure";
    }
  }
  detail::require_collective(mpi, mpi.rank() != root || root_ok_parse,
                             root_error);
  detail::broadcast_bytes(payload, mpi, root, maximum_chunk_bytes);

  std::shared_ptr<const detail::SurfaceStorage> storage;
  bool decode_ok = true;
  std::string decode_error;
  try {
    storage = detail::decode_normalized_surface(payload);
  } catch (const std::exception &error) {
    decode_ok = false;
    decode_error = error.what();
  } catch (...) {
    decode_ok = false;
    decode_error = "immersed surface: unknown normalized decode failure";
  }
  detail::require_collective(mpi, decode_ok, decode_error);

  std::uint64_t minimum_fingerprint = storage->fingerprint;
  std::uint64_t maximum_fingerprint = storage->fingerprint;
  detail::check_mpi(MPI_Allreduce(MPI_IN_PLACE, &minimum_fingerprint, 1,
                                  MPI_UINT64_T, MPI_MIN, mpi.comm()),
                    "MPI_Allreduce surface fingerprint minimum");
  detail::check_mpi(MPI_Allreduce(MPI_IN_PLACE, &maximum_fingerprint, 1,
                                  MPI_UINT64_T, MPI_MAX, mpi.comm()),
                    "MPI_Allreduce surface fingerprint maximum");
  detail::require_collective(
      mpi, minimum_fingerprint == maximum_fingerprint,
      "immersed surface: normalized fingerprints differ between ranks");
  return ImmersedSurface(std::move(storage));
}

std::size_t ImmersedSurface::vertex_count() const noexcept {
  return storage_->vertices.size();
}

std::size_t ImmersedSurface::triangle_count() const noexcept {
  return storage_->triangles.size();
}

const SurfaceTriangle &ImmersedSurface::triangle(TriangleId id) const {
  if (id >= storage_->triangles.size()) {
    throw runtime::Error("immersed surface: triangle ID is out of range");
  }
  return storage_->triangles[static_cast<std::size_t>(id)];
}

runtime::Real3 ImmersedSurface::bounding_box_min_m() const noexcept {
  return storage_->bounding_box_min;
}

runtime::Real3 ImmersedSurface::bounding_box_max_m() const noexcept {
  return storage_->bounding_box_max;
}

double ImmersedSurface::closed_volume_m3() const noexcept {
  return storage_->closed_volume_m3;
}

std::uint64_t ImmersedSurface::fingerprint() const noexcept {
  return storage_->fingerprint;
}

} // namespace hundun::immersed
