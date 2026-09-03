// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_ibm.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <new>
#include <numeric>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace hundun::v04 {
namespace {

constexpr std::uint64_t kFnvOffset = UINT64_C(14695981039346656037);
constexpr std::uint64_t kFnvPrime = UINT64_C(1099511628211);
constexpr std::uint32_t kSurfaceDetailInput = 1U;
constexpr std::uint32_t kSurfaceDetailTopology = 2U;
constexpr std::uint32_t kSurfaceDetailWinding = 3U;
constexpr std::uint32_t kSurfaceDetailVolume = 4U;
constexpr std::uint32_t kSurfaceDetailQuery = 5U;
constexpr std::uint32_t kSurfaceDetailAllocation = 6U;
constexpr std::uint32_t kBvhLeafTriangles = 4U;
constexpr double kQueryRoundoff = 128.0;

Status invalid_surface(std::uint32_t detail) noexcept {
  return {StatusCode::invalid_plan, detail};
}

Status allocation_failure() noexcept {
  return {StatusCode::allocation_failure, kSurfaceDetailAllocation};
}

bool finite(Real3 value) noexcept {
  return std::isfinite(value.x) && std::isfinite(value.y) &&
         std::isfinite(value.z);
}

double canonical_zero(double value) noexcept {
  return value == 0.0 ? 0.0 : value;
}

Real3 canonical_zero(Real3 value) noexcept {
  return {canonical_zero(value.x), canonical_zero(value.y),
          canonical_zero(value.z)};
}

Real3 add(Real3 lhs, Real3 rhs) noexcept {
  return {lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z};
}

Real3 subtract(Real3 lhs, Real3 rhs) noexcept {
  return {lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z};
}

Real3 multiply(double scalar, Real3 value) noexcept {
  return {scalar * value.x, scalar * value.y, scalar * value.z};
}

double dot(Real3 lhs, Real3 rhs) noexcept {
  return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

Real3 cross(Real3 lhs, Real3 rhs) noexcept {
  return {lhs.y * rhs.z - lhs.z * rhs.y,
          lhs.z * rhs.x - lhs.x * rhs.z,
          lhs.x * rhs.y - lhs.y * rhs.x};
}

double norm_squared(Real3 value) noexcept { return dot(value, value); }

Real3 minimum(Real3 lhs, Real3 rhs) noexcept {
  return {std::min(lhs.x, rhs.x), std::min(lhs.y, rhs.y),
          std::min(lhs.z, rhs.z)};
}

Real3 maximum(Real3 lhs, Real3 rhs) noexcept {
  return {std::max(lhs.x, rhs.x), std::max(lhs.y, rhs.y),
          std::max(lhs.z, rhs.z)};
}

double component(Real3 value, int axis) noexcept {
  return axis == 0 ? value.x : (axis == 1 ? value.y : value.z);
}

bool coordinate_less(Real3 lhs, Real3 rhs) noexcept {
  if (lhs.x != rhs.x) {
    return lhs.x < rhs.x;
  }
  if (lhs.y != rhs.y) {
    return lhs.y < rhs.y;
  }
  return lhs.z < rhs.z;
}

bool coordinate_equal(Real3 lhs, Real3 rhs) noexcept {
  return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
}

bool vertices_less(const std::array<Real3, 3U>& lhs,
                   const std::array<Real3, 3U>& rhs) noexcept {
  for (std::size_t vertex = 0U; vertex < 3U; ++vertex) {
    if (coordinate_less(lhs[vertex], rhs[vertex])) {
      return true;
    }
    if (coordinate_less(rhs[vertex], lhs[vertex])) {
      return false;
    }
  }
  return false;
}

bool vertices_equal(const std::array<Real3, 3U>& lhs,
                    const std::array<Real3, 3U>& rhs) noexcept {
  return coordinate_equal(lhs[0], rhs[0]) &&
         coordinate_equal(lhs[1], rhs[1]) &&
         coordinate_equal(lhs[2], rhs[2]);
}

std::array<Real3, 3U> rotate_to_minimum(
    std::array<Real3, 3U> vertices) noexcept {
  std::size_t minimum_index = 0U;
  if (coordinate_less(vertices[1], vertices[minimum_index])) {
    minimum_index = 1U;
  }
  if (coordinate_less(vertices[2], vertices[minimum_index])) {
    minimum_index = 2U;
  }
  return {vertices[minimum_index], vertices[(minimum_index + 1U) % 3U],
          vertices[(minimum_index + 2U) % 3U]};
}

class Hash64 {
 public:
  void byte(std::uint8_t value) noexcept {
    value_ ^= value;
    value_ *= kFnvPrime;
  }

  template <class Integer>
  void integer(Integer value) noexcept {
    using Unsigned = std::make_unsigned_t<Integer>;
    const Unsigned bits = static_cast<Unsigned>(value);
    for (std::size_t shift = 0U; shift < sizeof(bits) * 8U; shift += 8U) {
      byte(static_cast<std::uint8_t>((bits >> shift) & Unsigned{0xffU}));
    }
  }

  void real(double value) noexcept {
    value = canonical_zero(value);
    std::uint64_t bits{};
    static_assert(sizeof(bits) == sizeof(value));
    std::memcpy(&bits, &value, sizeof(bits));
    integer(bits);
  }

  void real3(Real3 value) noexcept {
    real(value.x);
    real(value.y);
    real(value.z);
  }

  PlanFingerprint finish() const noexcept {
    return value_ == 0U ? PlanFingerprint{1U} : value_;
  }

 private:
  std::uint64_t value_{kFnvOffset};
};

struct WorkingTriangle {
  std::array<Real3, 3U> oriented{};
  std::array<Real3, 3U> unordered_key{};
};

struct EdgeRecord {
  Real3 lower{};
  Real3 upper{};
  std::uint32_t triangle{};
  std::uint8_t local_edge{};
  bool lower_to_upper{};
};

bool edge_less(const EdgeRecord& lhs, const EdgeRecord& rhs) noexcept {
  if (coordinate_less(lhs.lower, rhs.lower)) {
    return true;
  }
  if (coordinate_less(rhs.lower, lhs.lower)) {
    return false;
  }
  if (coordinate_less(lhs.upper, rhs.upper)) {
    return true;
  }
  if (coordinate_less(rhs.upper, lhs.upper)) {
    return false;
  }
  if (lhs.triangle != rhs.triangle) {
    return lhs.triangle < rhs.triangle;
  }
  return lhs.local_edge < rhs.local_edge;
}

bool same_edge(const EdgeRecord& lhs, const EdgeRecord& rhs) noexcept {
  return coordinate_equal(lhs.lower, rhs.lower) &&
         coordinate_equal(lhs.upper, rhs.upper);
}

struct VertexRecord {
  Real3 point{};
  std::uint32_t triangle{};
  std::uint8_t local_vertex{};
};

bool vertex_less(const VertexRecord& lhs, const VertexRecord& rhs) noexcept {
  if (coordinate_less(lhs.point, rhs.point)) {
    return true;
  }
  if (coordinate_less(rhs.point, lhs.point)) {
    return false;
  }
  if (lhs.triangle != rhs.triangle) {
    return lhs.triangle < rhs.triangle;
  }
  return lhs.local_vertex < rhs.local_vertex;
}

long double translated_six_volume(
    const std::array<Real3, 3U>& triangle, Real3 reference) noexcept {
  const long double ax =
      static_cast<long double>(triangle[0].x) - reference.x;
  const long double ay =
      static_cast<long double>(triangle[0].y) - reference.y;
  const long double az =
      static_cast<long double>(triangle[0].z) - reference.z;
  const long double bx =
      static_cast<long double>(triangle[1].x) - reference.x;
  const long double by =
      static_cast<long double>(triangle[1].y) - reference.y;
  const long double bz =
      static_cast<long double>(triangle[1].z) - reference.z;
  const long double cx =
      static_cast<long double>(triangle[2].x) - reference.x;
  const long double cy =
      static_cast<long double>(triangle[2].y) - reference.y;
  const long double cz =
      static_cast<long double>(triangle[2].z) - reference.z;
  return ax * (by * cz - bz * cy) + ay * (bz * cx - bx * cz) +
         az * (bx * cy - by * cx);
}

struct ComponentGeometry {
  Real3 lower{};
  Real3 upper{};
  Real3 representative{};
  double scale_squared{};
};

enum class PointComponentRelation : std::uint8_t {
  outside,
  inside,
  boundary,
  invalid
};

Real3 closest_triangle_point(
    Real3 point, const std::array<Real3, 3U>& vertices) noexcept;

PointComponentRelation classify_component_point(
    Real3 point, const std::vector<WorkingTriangle>& triangles,
    const std::vector<std::uint32_t>& component_members,
    std::size_t begin, std::size_t end, double scale_squared) noexcept {
  if (!finite(point) || begin >= end || end > component_members.size() ||
      !(scale_squared > 0.0) || !std::isfinite(scale_squared)) {
    return PointComponentRelation::invalid;
  }
  const double distance_tolerance =
      4096.0 * std::numeric_limits<double>::epsilon() *
      std::sqrt(scale_squared);
  const double distance_tolerance_squared =
      distance_tolerance * distance_tolerance;
  long double solid_angle = 0.0L;
  for (std::size_t member = begin; member < end; ++member) {
    const std::uint32_t triangle_id = component_members[member];
    if (triangle_id >= triangles.size()) {
      return PointComponentRelation::invalid;
    }
    const auto& vertices = triangles[triangle_id].oriented;
    const Real3 closest = closest_triangle_point(point, vertices);
    const double distance_squared = norm_squared(subtract(point, closest));
    if (!std::isfinite(distance_squared)) {
      return PointComponentRelation::invalid;
    }
    if (distance_squared <= distance_tolerance_squared) {
      return PointComponentRelation::boundary;
    }

    const long double ax =
        static_cast<long double>(vertices[0].x) - point.x;
    const long double ay =
        static_cast<long double>(vertices[0].y) - point.y;
    const long double az =
        static_cast<long double>(vertices[0].z) - point.z;
    const long double bx =
        static_cast<long double>(vertices[1].x) - point.x;
    const long double by =
        static_cast<long double>(vertices[1].y) - point.y;
    const long double bz =
        static_cast<long double>(vertices[1].z) - point.z;
    const long double cx =
        static_cast<long double>(vertices[2].x) - point.x;
    const long double cy =
        static_cast<long double>(vertices[2].y) - point.y;
    const long double cz =
        static_cast<long double>(vertices[2].z) - point.z;
    const long double a_norm = std::sqrt(ax * ax + ay * ay + az * az);
    const long double b_norm = std::sqrt(bx * bx + by * by + bz * bz);
    const long double c_norm = std::sqrt(cx * cx + cy * cy + cz * cz);
    if (!(a_norm > 0.0L) || !(b_norm > 0.0L) || !(c_norm > 0.0L) ||
        !std::isfinite(a_norm) || !std::isfinite(b_norm) ||
        !std::isfinite(c_norm)) {
      return PointComponentRelation::boundary;
    }
    const long double numerator =
        ax * (by * cz - bz * cy) + ay * (bz * cx - bx * cz) +
        az * (bx * cy - by * cx);
    const long double denominator =
        a_norm * b_norm * c_norm +
        (ax * bx + ay * by + az * bz) * c_norm +
        (bx * cx + by * cy + bz * cz) * a_norm +
        (cx * ax + cy * ay + cz * az) * b_norm;
    if (!std::isfinite(numerator) || !std::isfinite(denominator)) {
      return PointComponentRelation::invalid;
    }
    solid_angle += 2.0L * std::atan2(numerator, denominator);
  }
  if (!std::isfinite(solid_angle)) {
    return PointComponentRelation::invalid;
  }
  const long double pi = std::acos(-1.0L);
  const long double magnitude = std::abs(solid_angle);
  if (magnitude > 3.0L * pi) {
    return PointComponentRelation::inside;
  }
  if (magnitude < pi) {
    return PointComponentRelation::outside;
  }
  return PointComponentRelation::invalid;
}

double closest_box_distance_squared(Real3 point, Real3 lower,
                                    Real3 upper) noexcept {
  double distance = 0.0;
  for (int axis = 0; axis < 3; ++axis) {
    const double coordinate = component(point, axis);
    const double minimum_value = component(lower, axis);
    const double maximum_value = component(upper, axis);
    const double delta = coordinate < minimum_value
                             ? minimum_value - coordinate
                             : (coordinate > maximum_value
                                    ? coordinate - maximum_value
                                    : 0.0);
    distance += delta * delta;
  }
  return distance;
}

Real3 closest_triangle_point(Real3 point,
                             const std::array<Real3, 3U>& vertices) noexcept {
  const Real3 a = vertices[0];
  const Real3 b = vertices[1];
  const Real3 c = vertices[2];
  const Real3 ab = subtract(b, a);
  const Real3 ac = subtract(c, a);
  const Real3 ap = subtract(point, a);
  const double d1 = dot(ab, ap);
  const double d2 = dot(ac, ap);
  if (d1 <= 0.0 && d2 <= 0.0) {
    return a;
  }

  const Real3 bp = subtract(point, b);
  const double d3 = dot(ab, bp);
  const double d4 = dot(ac, bp);
  if (d3 >= 0.0 && d4 <= d3) {
    return b;
  }

  const double vc = d1 * d4 - d3 * d2;
  if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
    const double fraction = d1 / (d1 - d3);
    return add(a, multiply(fraction, ab));
  }

  const Real3 cp = subtract(point, c);
  const double d5 = dot(ab, cp);
  const double d6 = dot(ac, cp);
  if (d6 >= 0.0 && d5 <= d6) {
    return c;
  }

  const double vb = d5 * d2 - d1 * d6;
  if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
    const double fraction = d2 / (d2 - d6);
    return add(a, multiply(fraction, ac));
  }

  const double va = d3 * d6 - d5 * d4;
  if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
    const Real3 bc = subtract(c, b);
    const double fraction = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    return add(b, multiply(fraction, bc));
  }

  const double denominator = 1.0 / (va + vb + vc);
  const double v = vb * denominator;
  const double w = vc * denominator;
  return add(a, add(multiply(v, ab), multiply(w, ac)));
}

double comparison_tolerance(double lhs, double rhs,
                            double scene_scale_squared) noexcept {
  return kQueryRoundoff * std::numeric_limits<double>::epsilon() *
         std::max({scene_scale_squared, std::abs(lhs), std::abs(rhs),
                   std::numeric_limits<double>::min()});
}

bool segment_box_interval(Real3 begin, Real3 direction, Real3 lower,
                          Real3 upper, double& entry) noexcept {
  double minimum_fraction = 0.0;
  double maximum_fraction = 1.0;
  for (int axis = 0; axis < 3; ++axis) {
    const double origin = component(begin, axis);
    const double delta = component(direction, axis);
    const double minimum_value = component(lower, axis);
    const double maximum_value = component(upper, axis);
    if (delta == 0.0) {
      if (origin < minimum_value || origin > maximum_value) {
        return false;
      }
      continue;
    }
    double first = (minimum_value - origin) / delta;
    double second = (maximum_value - origin) / delta;
    if (first > second) {
      std::swap(first, second);
    }
    minimum_fraction = std::max(minimum_fraction, first);
    maximum_fraction = std::min(maximum_fraction, second);
    if (minimum_fraction > maximum_fraction) {
      return false;
    }
  }
  entry = minimum_fraction;
  return maximum_fraction >= 0.0 && minimum_fraction <= 1.0;
}

bool segment_triangle_fraction(Real3 begin, Real3 direction,
                               const SurfaceTriangle& triangle,
                               double& fraction) noexcept {
  const Real3 edge1 =
      subtract(triangle.vertices[1], triangle.vertices[0]);
  const Real3 edge2 =
      subtract(triangle.vertices[2], triangle.vertices[0]);
  const Real3 p = cross(direction, edge2);
  const double determinant = dot(edge1, p);
  const double determinant_scale =
      std::sqrt(norm_squared(direction) * norm_squared(edge1) *
                norm_squared(edge2));
  const double determinant_tolerance =
      kQueryRoundoff * std::numeric_limits<double>::epsilon() *
      determinant_scale;
  if (!(std::abs(determinant) > determinant_tolerance)) {
    return false;
  }
  const double inverse = 1.0 / determinant;
  const Real3 from_vertex = subtract(begin, triangle.vertices[0]);
  const double u = dot(from_vertex, p) * inverse;
  const Real3 q = cross(from_vertex, edge1);
  const double v = dot(direction, q) * inverse;
  const double candidate = dot(edge2, q) * inverse;
  const double barycentric_tolerance =
      kQueryRoundoff * std::numeric_limits<double>::epsilon();
  if (u < -barycentric_tolerance || v < -barycentric_tolerance ||
      u + v > 1.0 + barycentric_tolerance ||
      candidate < -barycentric_tolerance ||
      candidate > 1.0 + barycentric_tolerance) {
    return false;
  }
  fraction = std::clamp(candidate, 0.0, 1.0);
  return true;
}

}  // namespace

Status ImmersedSurfaceCompiler::compile(const StlScanPlan& scan,
                                        ImmersedSurfacePlan& out) noexcept {
  try {
    const TriangleSoA& source = scan.triangles();
    const std::size_t triangle_count = source.size();
    const auto valid_size = [triangle_count](Span<const double> values) {
      return values.size == triangle_count &&
             (triangle_count == 0U || values.data != nullptr);
    };
    if (scan.fingerprint() == 0U || source.fingerprint() == 0U ||
        triangle_count == 0U ||
        triangle_count > static_cast<std::size_t>(UINT32_MAX) ||
        !valid_size(source.ax()) || !valid_size(source.ay()) ||
        !valid_size(source.az()) || !valid_size(source.e1x()) ||
        !valid_size(source.bx()) || !valid_size(source.by()) ||
        !valid_size(source.bz()) || !valid_size(source.cx()) ||
        !valid_size(source.cy()) || !valid_size(source.cz()) ||
        !valid_size(source.e1y()) || !valid_size(source.e1z()) ||
        !valid_size(source.e2x()) || !valid_size(source.e2y()) ||
        !valid_size(source.e2z()) || !valid_size(source.nx()) ||
        !valid_size(source.ny()) || !valid_size(source.nz()) ||
        !valid_size(source.min_x()) || !valid_size(source.min_y()) ||
        !valid_size(source.min_z()) || !valid_size(source.max_x()) ||
        !valid_size(source.max_y()) || !valid_size(source.max_z())) {
      return invalid_surface(kSurfaceDetailInput);
    }

    std::vector<WorkingTriangle> triangles;
    triangles.reserve(triangle_count);
    for (std::size_t index = 0U; index < triangle_count; ++index) {
      const Real3 a = canonical_zero(
          {source.ax().data[index], source.ay().data[index],
           source.az().data[index]});
      const Real3 b = canonical_zero(
          {source.bx().data[index], source.by().data[index],
           source.bz().data[index]});
      const Real3 c = canonical_zero(
          {source.cx().data[index], source.cy().data[index],
           source.cz().data[index]});
      const Real3 stored_normal{source.nx().data[index],
                                source.ny().data[index],
                                source.nz().data[index]};
      const Real3 raw_normal = cross(subtract(b, a), subtract(c, a));
      const double normal_squared = norm_squared(raw_normal);
      if (!finite(a) || !finite(b) || !finite(c) ||
          !finite(stored_normal) || !(normal_squared > 0.0) ||
          !std::isfinite(normal_squared) ||
          !(dot(raw_normal, stored_normal) > 0.0)) {
        return invalid_surface(kSurfaceDetailInput);
      }
      WorkingTriangle triangle;
      triangle.oriented = rotate_to_minimum({a, b, c});
      triangle.unordered_key = {a, b, c};
      std::sort(triangle.unordered_key.begin(),
                triangle.unordered_key.end(), coordinate_less);
      triangles.push_back(triangle);
    }
    std::sort(triangles.begin(), triangles.end(),
              [](const WorkingTriangle& lhs,
                 const WorkingTriangle& rhs) noexcept {
                return vertices_less(lhs.unordered_key, rhs.unordered_key);
              });
    for (std::size_t index = 1U; index < triangles.size(); ++index) {
      if (vertices_equal(triangles[index - 1U].unordered_key,
                         triangles[index].unordered_key)) {
        return invalid_surface(kSurfaceDetailTopology);
      }
    }

    if (triangle_count > std::numeric_limits<std::size_t>::max() / 3U) {
      return invalid_surface(kSurfaceDetailInput);
    }
    std::vector<EdgeRecord> edges;
    edges.reserve(triangle_count * 3U);
    std::vector<VertexRecord> vertices;
    vertices.reserve(triangle_count * 3U);
    for (std::size_t triangle_id = 0U; triangle_id < triangle_count;
         ++triangle_id) {
      const auto& triangle = triangles[triangle_id].oriented;
      for (std::size_t local = 0U; local < 3U; ++local) {
        vertices.push_back({triangle[local],
                            static_cast<std::uint32_t>(triangle_id),
                            static_cast<std::uint8_t>(local)});
        const Real3 first = triangle[local];
        const Real3 second = triangle[(local + 1U) % 3U];
        const bool forward = coordinate_less(first, second);
        edges.push_back({forward ? first : second,
                         forward ? second : first,
                         static_cast<std::uint32_t>(triangle_id),
                         static_cast<std::uint8_t>(local), forward});
      }
    }
    std::sort(edges.begin(), edges.end(), edge_less);
    std::vector<std::array<std::uint32_t, 3U>> neighbours(
        triangle_count,
        {kInvalidIbmIndex, kInvalidIbmIndex, kInvalidIbmIndex});
    for (std::size_t begin = 0U; begin < edges.size();) {
      std::size_t end = begin + 1U;
      while (end < edges.size() && same_edge(edges[begin], edges[end])) {
        ++end;
      }
      if (end - begin != 2U ||
          edges[begin].lower_to_upper == edges[begin + 1U].lower_to_upper ||
          edges[begin].triangle == edges[begin + 1U].triangle) {
        return invalid_surface(end - begin == 2U ? kSurfaceDetailWinding
                                                 : kSurfaceDetailTopology);
      }
      const EdgeRecord& first = edges[begin];
      const EdgeRecord& second = edges[begin + 1U];
      neighbours[first.triangle][first.local_edge] = second.triangle;
      neighbours[second.triangle][second.local_edge] = first.triangle;
      begin = end;
    }

    std::sort(vertices.begin(), vertices.end(), vertex_less);
    std::vector<std::uint64_t> vertex_visit(triangle_count, 0U);
    std::vector<std::uint32_t> stack;
    stack.reserve(triangle_count);
    std::uint64_t visit_token = 0U;
    for (std::size_t begin = 0U; begin < vertices.size();) {
      std::size_t end = begin + 1U;
      while (end < vertices.size() &&
             coordinate_equal(vertices[begin].point, vertices[end].point)) {
        ++end;
      }
      if (end - begin < 3U || visit_token == UINT64_MAX) {
        return invalid_surface(kSurfaceDetailTopology);
      }
      ++visit_token;
      stack.clear();
      stack.push_back(vertices[begin].triangle);
      vertex_visit[vertices[begin].triangle] = visit_token;
      std::size_t visited_count = 0U;
      while (!stack.empty()) {
        const std::uint32_t triangle_id = stack.back();
        stack.pop_back();
        ++visited_count;
        const auto& triangle = triangles[triangle_id].oriented;
        std::size_t local_vertex = 3U;
        for (std::size_t local = 0U; local < 3U; ++local) {
          if (coordinate_equal(triangle[local], vertices[begin].point)) {
            local_vertex = local;
            break;
          }
        }
        if (local_vertex == 3U) {
          return invalid_surface(kSurfaceDetailTopology);
        }
        const std::array<std::uint32_t, 2U> around_vertex{
            neighbours[triangle_id][local_vertex],
            neighbours[triangle_id][(local_vertex + 2U) % 3U]};
        if (around_vertex[0] == kInvalidIbmIndex ||
            around_vertex[1] == kInvalidIbmIndex ||
            around_vertex[0] == around_vertex[1]) {
          return invalid_surface(kSurfaceDetailTopology);
        }
        for (const std::uint32_t adjacent : around_vertex) {
          bool contains_vertex = false;
          for (const Real3 point : triangles[adjacent].oriented) {
            contains_vertex |= coordinate_equal(point, vertices[begin].point);
          }
          if (!contains_vertex) {
            return invalid_surface(kSurfaceDetailTopology);
          }
          if (vertex_visit[adjacent] != visit_token) {
            vertex_visit[adjacent] = visit_token;
            stack.push_back(adjacent);
          }
        }
      }
      if (visited_count != end - begin) {
        return invalid_surface(kSurfaceDetailTopology);
      }
      begin = end;
    }

    std::vector<std::uint32_t> component_ids(triangle_count,
                                             kInvalidIbmIndex);
    std::vector<std::size_t> component_offsets{0U};
    std::vector<std::uint32_t> component_members;
    component_members.reserve(triangle_count);
    for (std::size_t start = 0U; start < triangle_count; ++start) {
      if (component_ids[start] != kInvalidIbmIndex) {
        continue;
      }
      if (component_offsets.size() - 1U >=
          static_cast<std::size_t>(UINT32_MAX)) {
        return invalid_surface(kSurfaceDetailTopology);
      }
      const std::uint32_t component_id =
          static_cast<std::uint32_t>(component_offsets.size() - 1U);
      stack.clear();
      stack.push_back(static_cast<std::uint32_t>(start));
      component_ids[start] = component_id;
      while (!stack.empty()) {
        const std::uint32_t triangle_id = stack.back();
        stack.pop_back();
        component_members.push_back(triangle_id);
        for (const std::uint32_t adjacent : neighbours[triangle_id]) {
          if (adjacent == kInvalidIbmIndex) {
            return invalid_surface(kSurfaceDetailTopology);
          }
          if (component_ids[adjacent] == kInvalidIbmIndex) {
            component_ids[adjacent] = component_id;
            stack.push_back(adjacent);
          } else if (component_ids[adjacent] != component_id) {
            return invalid_surface(kSurfaceDetailTopology);
          }
        }
      }
      component_offsets.push_back(component_members.size());
    }

    Real3 global_lower = triangles.front().unordered_key.front();
    Real3 global_upper = global_lower;
    for (const WorkingTriangle& triangle : triangles) {
      for (const Real3 point : triangle.oriented) {
        global_lower = minimum(global_lower, point);
        global_upper = maximum(global_upper, point);
      }
    }
    const Real3 global_extent = subtract(global_upper, global_lower);
    if (!(global_extent.x > 0.0) || !(global_extent.y > 0.0) ||
        !(global_extent.z > 0.0) || !finite(global_extent)) {
      return invalid_surface(kSurfaceDetailVolume);
    }
    long double total_volume = 0.0L;
    std::vector<ComponentGeometry> component_geometry;
    component_geometry.reserve(component_offsets.size() - 1U);
    for (std::size_t component_id = 0U;
         component_id + 1U < component_offsets.size(); ++component_id) {
      const std::size_t begin = component_offsets[component_id];
      const std::size_t end = component_offsets[component_id + 1U];
      Real3 reference = triangles[component_members[begin]].oriented[0];
      Real3 component_lower = reference;
      Real3 component_upper = reference;
      for (std::size_t member = begin; member < end; ++member) {
        for (const Real3 point :
             triangles[component_members[member]].oriented) {
          component_lower = minimum(component_lower, point);
          component_upper = maximum(component_upper, point);
          if (coordinate_less(point, reference)) {
            reference = point;
          }
        }
      }
      const Real3 component_extent =
          subtract(component_upper, component_lower);
      if (!(component_extent.x > 0.0) || !(component_extent.y > 0.0) ||
          !(component_extent.z > 0.0) || !finite(component_extent)) {
        return invalid_surface(kSurfaceDetailVolume);
      }
      const long double component_box_volume =
          static_cast<long double>(component_extent.x) * component_extent.y *
          component_extent.z;
      const long double volume_tolerance =
          4096.0L * std::numeric_limits<double>::epsilon() *
          component_box_volume;
      component_geometry.push_back(
          {component_lower, component_upper,
           triangles[component_members[begin]].oriented[0],
           norm_squared(component_extent)});
      long double six_volume = 0.0L;
      for (std::size_t member = begin; member < end; ++member) {
        six_volume += translated_six_volume(
            triangles[component_members[member]].oriented, reference);
      }
      if (!(std::abs(six_volume) > 6.0L * volume_tolerance) ||
          !std::isfinite(six_volume)) {
        return invalid_surface(kSurfaceDetailVolume);
      }
      if (six_volume < 0.0L) {
        for (std::size_t member = begin; member < end; ++member) {
          std::swap(triangles[component_members[member]].oriented[1],
                    triangles[component_members[member]].oriented[2]);
        }
        six_volume = -six_volume;
      }
      total_volume += six_volume / 6.0L;
    }
    if (!(total_volume > 0.0L) || !std::isfinite(total_volume) ||
        total_volume > std::numeric_limits<double>::max()) {
      return invalid_surface(kSurfaceDetailVolume);
    }

    const auto point_can_be_inside_box = [](Real3 point,
                                            const ComponentGeometry& target) {
      const double tolerance =
          4096.0 * std::numeric_limits<double>::epsilon() *
          std::sqrt(target.scale_squared);
      return point.x >= target.lower.x - tolerance &&
             point.x <= target.upper.x + tolerance &&
             point.y >= target.lower.y - tolerance &&
             point.y <= target.upper.y + tolerance &&
             point.z >= target.lower.z - tolerance &&
             point.z <= target.upper.z + tolerance;
    };
    for (std::size_t first = 0U; first < component_geometry.size(); ++first) {
      for (std::size_t second = first + 1U;
           second < component_geometry.size(); ++second) {
        const auto classify_if_candidate =
            [&](std::size_t point_component,
                std::size_t target_component) noexcept {
              const Real3 point =
                  component_geometry[point_component].representative;
              const ComponentGeometry& target =
                  component_geometry[target_component];
              if (!point_can_be_inside_box(point, target)) {
                return PointComponentRelation::outside;
              }
              return classify_component_point(
                  point, triangles, component_members,
                  component_offsets[target_component],
                  component_offsets[target_component + 1U],
                  target.scale_squared);
            };
        const PointComponentRelation first_in_second =
            classify_if_candidate(first, second);
        const PointComponentRelation second_in_first =
            classify_if_candidate(second, first);
        if (first_in_second != PointComponentRelation::outside ||
            second_in_first != PointComponentRelation::outside) {
          return invalid_surface(kSurfaceDetailTopology);
        }
      }
    }

    ImmersedSurfacePlan candidate;
    candidate.triangles_.reserve(triangle_count);
    for (std::size_t triangle_id = 0U; triangle_id < triangle_count;
         ++triangle_id) {
      const auto& vertices_out = triangles[triangle_id].oriented;
      const Real3 raw_normal =
          cross(subtract(vertices_out[1], vertices_out[0]),
                subtract(vertices_out[2], vertices_out[0]));
      const double magnitude = std::sqrt(norm_squared(raw_normal));
      if (!(magnitude > 0.0) || !std::isfinite(magnitude)) {
        return invalid_surface(kSurfaceDetailInput);
      }
      SurfaceTriangle triangle;
      triangle.id = static_cast<SurfaceTriangleId>(triangle_id);
      triangle.vertices = vertices_out;
      triangle.geometric_outward_normal = canonical_zero(
          multiply(1.0 / magnitude, raw_normal));
      triangle.centroid = canonical_zero(multiply(
          1.0 / 3.0,
          add(vertices_out[0], add(vertices_out[1], vertices_out[2]))));
      triangle.area = 0.5 * magnitude;
      candidate.triangles_.push_back(triangle);
    }
    for (auto& adjacent : neighbours) {
      std::sort(adjacent.begin(), adjacent.end());
    }
    candidate.triangle_neighbours_ = std::move(neighbours);
    candidate.bounding_box_min_ = global_lower;
    candidate.bounding_box_max_ = global_upper;
    candidate.closed_volume_ = static_cast<double>(total_volume);
    candidate.source_triangle_fingerprint_ = source.fingerprint();

    candidate.bvh_triangle_ids_.resize(triangle_count);
    std::iota(candidate.bvh_triangle_ids_.begin(),
              candidate.bvh_triangle_ids_.end(), 0U);
    if (triangle_count >
        (std::numeric_limits<std::size_t>::max() - 1U) / 2U) {
      return invalid_surface(kSurfaceDetailInput);
    }
    candidate.bvh_nodes_.reserve(triangle_count * 2U - 1U);
    const auto triangle_bounds = [&candidate](std::uint32_t triangle_id,
                                               Real3& lower,
                                               Real3& upper) noexcept {
      const auto& vertices = candidate.triangles_[triangle_id].vertices;
      lower = minimum(vertices[0], minimum(vertices[1], vertices[2]));
      upper = maximum(vertices[0], maximum(vertices[1], vertices[2]));
    };
    const auto build_bvh = [&candidate, &triangle_bounds](
                               auto&& self, std::uint32_t begin,
                               std::uint32_t end) -> std::uint32_t {
      if (candidate.bvh_nodes_.size() >=
          static_cast<std::size_t>(UINT32_MAX)) {
        throw std::length_error("surface BVH node limit");
      }
      const std::uint32_t node_id =
          static_cast<std::uint32_t>(candidate.bvh_nodes_.size());
      candidate.bvh_nodes_.push_back({});
      Real3 lower{};
      Real3 upper{};
      triangle_bounds(candidate.bvh_triangle_ids_[begin], lower, upper);
      for (std::uint32_t index = begin + 1U; index < end; ++index) {
        Real3 triangle_lower{};
        Real3 triangle_upper{};
        triangle_bounds(candidate.bvh_triangle_ids_[index], triangle_lower,
                        triangle_upper);
        lower = minimum(lower, triangle_lower);
        upper = maximum(upper, triangle_upper);
      }
      candidate.bvh_nodes_[node_id].lower = lower;
      candidate.bvh_nodes_[node_id].upper = upper;
      candidate.bvh_nodes_[node_id].begin = begin;
      const std::uint32_t count = end - begin;
      if (count <= kBvhLeafTriangles) {
        candidate.bvh_nodes_[node_id].count = count;
        return node_id;
      }
      const Real3 extent = subtract(upper, lower);
      int axis = 0;
      if (extent.y > extent.x) {
        axis = 1;
      }
      if (component(extent, 2) > component(extent, axis)) {
        axis = 2;
      }
      std::sort(candidate.bvh_triangle_ids_.begin() + begin,
                candidate.bvh_triangle_ids_.begin() + end,
                [&candidate, axis](std::uint32_t lhs,
                                   std::uint32_t rhs) noexcept {
                  const double first = component(
                      candidate.triangles_[lhs].centroid, axis);
                  const double second = component(
                      candidate.triangles_[rhs].centroid, axis);
                  return first != second ? first < second : lhs < rhs;
                });
      const std::uint32_t middle = begin + count / 2U;
      const std::uint32_t left = self(self, begin, middle);
      const std::uint32_t right = self(self, middle, end);
      candidate.bvh_nodes_[node_id].left = left;
      candidate.bvh_nodes_[node_id].right = right;
      return node_id;
    };
    static_cast<void>(build_bvh(build_bvh, 0U,
                                static_cast<std::uint32_t>(triangle_count)));

    Hash64 hash;
    hash.integer(static_cast<std::uint64_t>(candidate.triangles_.size()));
    hash.real3(candidate.bounding_box_min_);
    hash.real3(candidate.bounding_box_max_);
    hash.real(candidate.closed_volume_);
    for (const SurfaceTriangle& triangle : candidate.triangles_) {
      hash.integer(triangle.id);
      for (const Real3 point : triangle.vertices) {
        hash.real3(point);
      }
      hash.real3(triangle.geometric_outward_normal);
      hash.real3(triangle.centroid);
      hash.real(triangle.area);
    }
    for (const auto& adjacent : candidate.triangle_neighbours_) {
      hash.integer(adjacent[0U]);
      hash.integer(adjacent[1U]);
      hash.integer(adjacent[2U]);
    }
    candidate.fingerprint_ = hash.finish();
    out = std::move(candidate);
    return {};
  } catch (const std::bad_alloc&) {
    return allocation_failure();
  } catch (const std::length_error&) {
    return invalid_surface(kSurfaceDetailInput);
  } catch (...) {
    return invalid_surface(kSurfaceDetailInput);
  }
}

Status ImmersedSurfacePlan::closest_point(
    Real3 point, ClosestSurfacePoint& out) const noexcept {
  if (!finite(point) || fingerprint_ == 0U || triangles_.empty() ||
      bvh_triangle_ids_.size() != triangles_.size() || bvh_nodes_.empty()) {
    return invalid_surface(kSurfaceDetailQuery);
  }
  const Real3 extent = subtract(bounding_box_max_, bounding_box_min_);
  const double scene_scale_squared = norm_squared(extent);
  if (!(scene_scale_squared > 0.0) || !std::isfinite(scene_scale_squared)) {
    return invalid_surface(kSurfaceDetailQuery);
  }

  ClosestSurfacePoint candidate;
  candidate.squared_distance = std::numeric_limits<double>::infinity();
  const auto visit = [this, point, scene_scale_squared, &candidate](
                         auto&& self, std::uint32_t node_id) noexcept -> void {
    if (node_id >= bvh_nodes_.size()) {
      return;
    }
    const BvhNode& node = bvh_nodes_[node_id];
    const double bound =
        closest_box_distance_squared(point, node.lower, node.upper);
    if (candidate.triangle != kInvalidSurfaceTriangle &&
        bound > candidate.squared_distance +
                    comparison_tolerance(bound, candidate.squared_distance,
                                         scene_scale_squared)) {
      return;
    }
    if (node.count != 0U) {
      if (node.begin > bvh_triangle_ids_.size() ||
          node.count > bvh_triangle_ids_.size() - node.begin) {
        return;
      }
      for (std::uint32_t offset = 0U; offset < node.count; ++offset) {
        const std::uint32_t triangle_id =
            bvh_triangle_ids_[node.begin + offset];
        if (triangle_id >= triangles_.size()) {
          continue;
        }
        const SurfaceTriangle& triangle = triangles_[triangle_id];
        const Real3 closest = closest_triangle_point(point, triangle.vertices);
        const double distance = norm_squared(subtract(point, closest));
        const double tolerance = comparison_tolerance(
            distance, candidate.squared_distance, scene_scale_squared);
        if (candidate.triangle == kInvalidSurfaceTriangle ||
            distance < candidate.squared_distance - tolerance ||
            (std::abs(distance - candidate.squared_distance) <= tolerance &&
             triangle.id < candidate.triangle)) {
          candidate = {triangle.id, canonical_zero(closest),
                       triangle.geometric_outward_normal, distance};
        }
      }
      return;
    }
    if (node.left >= bvh_nodes_.size() || node.right >= bvh_nodes_.size()) {
      return;
    }
    const double left_bound = closest_box_distance_squared(
        point, bvh_nodes_[node.left].lower, bvh_nodes_[node.left].upper);
    const double right_bound = closest_box_distance_squared(
        point, bvh_nodes_[node.right].lower, bvh_nodes_[node.right].upper);
    if (right_bound < left_bound) {
      self(self, node.right);
      self(self, node.left);
    } else {
      self(self, node.left);
      self(self, node.right);
    }
  };
  visit(visit, 0U);
  if (candidate.triangle == kInvalidSurfaceTriangle ||
      !std::isfinite(candidate.squared_distance)) {
    return invalid_surface(kSurfaceDetailQuery);
  }
  out = candidate;
  return {};
}

Status ImmersedSurfacePlan::first_segment_intersection(
    Real3 begin, Real3 end, SurfaceSegmentIntersection& out) const noexcept {
  const Real3 direction = subtract(end, begin);
  const double direction_squared = norm_squared(direction);
  if (!finite(begin) || !finite(end) || !(direction_squared > 0.0) ||
      !std::isfinite(direction_squared) || fingerprint_ == 0U ||
      triangles_.empty() || bvh_triangle_ids_.size() != triangles_.size() ||
      bvh_nodes_.empty()) {
    return invalid_surface(kSurfaceDetailQuery);
  }

  SurfaceSegmentIntersection candidate;
  double best_fraction = std::numeric_limits<double>::infinity();
  const double fraction_tolerance =
      kQueryRoundoff * std::numeric_limits<double>::epsilon();
  const auto visit = [this, begin, direction, fraction_tolerance, &candidate,
                      &best_fraction](auto&& self,
                                      std::uint32_t node_id) noexcept -> void {
    if (node_id >= bvh_nodes_.size()) {
      return;
    }
    const BvhNode& node = bvh_nodes_[node_id];
    double entry = 0.0;
    if (!segment_box_interval(begin, direction, node.lower, node.upper,
                              entry) ||
        entry > best_fraction + fraction_tolerance) {
      return;
    }
    if (node.count != 0U) {
      if (node.begin > bvh_triangle_ids_.size() ||
          node.count > bvh_triangle_ids_.size() - node.begin) {
        return;
      }
      for (std::uint32_t offset = 0U; offset < node.count; ++offset) {
        const std::uint32_t triangle_id =
            bvh_triangle_ids_[node.begin + offset];
        if (triangle_id >= triangles_.size()) {
          continue;
        }
        const SurfaceTriangle& triangle = triangles_[triangle_id];
        double fraction = 0.0;
        if (!segment_triangle_fraction(begin, direction, triangle,
                                       fraction)) {
          continue;
        }
        if (candidate.triangle == kInvalidSurfaceTriangle ||
            fraction < best_fraction - fraction_tolerance ||
            (std::abs(fraction - best_fraction) <= fraction_tolerance &&
             triangle.id < candidate.triangle)) {
          best_fraction = fraction;
          candidate = {triangle.id,
                       canonical_zero(add(begin,
                                          multiply(fraction, direction))),
                       fraction};
        }
      }
      return;
    }
    if (node.left >= bvh_nodes_.size() || node.right >= bvh_nodes_.size()) {
      return;
    }
    double left_entry = 0.0;
    double right_entry = 0.0;
    const bool left_hit = segment_box_interval(
        begin, direction, bvh_nodes_[node.left].lower,
        bvh_nodes_[node.left].upper, left_entry);
    const bool right_hit = segment_box_interval(
        begin, direction, bvh_nodes_[node.right].lower,
        bvh_nodes_[node.right].upper, right_entry);
    if (left_hit && right_hit) {
      if (right_entry < left_entry) {
        self(self, node.right);
        self(self, node.left);
      } else {
        self(self, node.left);
        self(self, node.right);
      }
    } else if (left_hit) {
      self(self, node.left);
    } else if (right_hit) {
      self(self, node.right);
    }
  };
  visit(visit, 0U);
  out = candidate;
  return {};
}

}  // namespace hundun::v04
