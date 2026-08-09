// SPDX-License-Identifier: Apache-2.0

#include "hundun/ib_surface_query.hpp"

#include "hundun/cfg_resolved_case_v3.hpp"
#include "hundun/rt_error.hpp"

#include "ib_surface_bvh_detail.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <functional>
#include <limits>
#include <memory>
#include <queue>
#include <string>
#include <utility>
#include <vector>

namespace hundun::immersed::detail {
namespace {

constexpr std::uint64_t fnv_offset = UINT64_C(14695981039346656037);
constexpr std::uint64_t fnv_prime = UINT64_C(1099511628211);

[[noreturn]] void fail(const char *message) {
  throw runtime::Error(std::string{"surface query: "} + message);
}

std::uint64_t double_bits(double value) noexcept {
  std::uint64_t bits = 0U;
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

void hash_u64(std::uint64_t &hash, std::uint64_t value) noexcept {
  for (unsigned shift = 0U; shift < 64U; shift += 8U) {
    hash ^= static_cast<std::uint8_t>((value >> shift) & 0xffU);
    hash *= fnv_prime;
  }
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

AxisAlignedBox triangle_box(const SurfaceTriangle &triangle) noexcept {
  AxisAlignedBox box{triangle.vertices_m[0], triangle.vertices_m[0]};
  for (std::size_t vertex = 1U; vertex < 3U; ++vertex) {
    box.minimum.x = std::min(box.minimum.x, triangle.vertices_m[vertex].x);
    box.minimum.y = std::min(box.minimum.y, triangle.vertices_m[vertex].y);
    box.minimum.z = std::min(box.minimum.z, triangle.vertices_m[vertex].z);
    box.maximum.x = std::max(box.maximum.x, triangle.vertices_m[vertex].x);
    box.maximum.y = std::max(box.maximum.y, triangle.vertices_m[vertex].y);
    box.maximum.z = std::max(box.maximum.z, triangle.vertices_m[vertex].z);
  }
  return box;
}

AxisAlignedBox merge(AxisAlignedBox a, AxisAlignedBox b) noexcept {
  return {
      {std::min(a.minimum.x, b.minimum.x), std::min(a.minimum.y, b.minimum.y),
       std::min(a.minimum.z, b.minimum.z)},
      {std::max(a.maximum.x, b.maximum.x), std::max(a.maximum.y, b.maximum.y),
       std::max(a.maximum.z, b.maximum.z)}};
}

runtime::Real3 centroid(const SurfaceTriangle &triangle) noexcept {
  const runtime::Real3 &origin = triangle.vertices_m[0];
  return add(origin, multiply(add(subtract(triangle.vertices_m[1], origin),
                                  subtract(triangle.vertices_m[2], origin)),
                              1.0 / 3.0));
}

bool boxes_overlap(const AxisAlignedBox &first,
                   const AxisAlignedBox &second) noexcept {
  return first.maximum.x >= second.minimum.x &&
         second.maximum.x >= first.minimum.x &&
         first.maximum.y >= second.minimum.y &&
         second.maximum.y >= first.minimum.y &&
         first.maximum.z >= second.minimum.z &&
         second.maximum.z >= first.minimum.z;
}

double box_distance_squared(const AxisAlignedBox &box,
                            runtime::Real3 point) noexcept {
  double result = 0.0;
  for (std::size_t axis = 0U; axis < 3U; ++axis) {
    const double value = coordinate(point, axis);
    const double minimum = coordinate(box.minimum, axis);
    const double maximum = coordinate(box.maximum, axis);
    const double difference = value < minimum
                                  ? minimum - value
                                  : (value > maximum ? value - maximum : 0.0);
    result += difference * difference;
  }
  return result;
}

runtime::Real3 closest_on_triangle(runtime::Real3 point,
                                   const SurfaceTriangle &triangle) noexcept {
  const runtime::Real3 &a = triangle.vertices_m[0];
  const runtime::Real3 &b = triangle.vertices_m[1];
  const runtime::Real3 &c = triangle.vertices_m[2];
  const runtime::Real3 ab = subtract(b, a);
  const runtime::Real3 ac = subtract(c, a);
  const runtime::Real3 ap = subtract(point, a);
  const double d1 = dot(ab, ap);
  const double d2 = dot(ac, ap);
  if (d1 <= 0.0 && d2 <= 0.0) {
    return a;
  }

  const runtime::Real3 bp = subtract(point, b);
  const double d3 = dot(ab, bp);
  const double d4 = dot(ac, bp);
  if (d3 >= 0.0 && d4 <= d3) {
    return b;
  }

  const double vc = d1 * d4 - d3 * d2;
  if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
    const double fraction = d1 / (d1 - d3);
    return add(a, multiply(ab, fraction));
  }

  const runtime::Real3 cp = subtract(point, c);
  const double d5 = dot(ab, cp);
  const double d6 = dot(ac, cp);
  if (d6 >= 0.0 && d5 <= d6) {
    return c;
  }

  const double vb = d5 * d2 - d1 * d6;
  if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
    const double fraction = d2 / (d2 - d6);
    return add(a, multiply(ac, fraction));
  }

  const double va = d3 * d6 - d5 * d4;
  if (va <= 0.0 && d4 - d3 >= 0.0 && d5 - d6 >= 0.0) {
    const runtime::Real3 bc = subtract(c, b);
    const double fraction = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    return add(b, multiply(bc, fraction));
  }

  const double denominator = 1.0 / (va + vb + vc);
  const double v = vb * denominator;
  const double w = vc * denominator;
  return add(a, add(multiply(ab, v), multiply(ac, w)));
}

enum class HitDisposition { no_hit, hit, coplanar };

struct TriangleHit final {
  TriangleId triangle{};
  double fraction{};
  runtime::Real3 point{};
};

HitDisposition segment_triangle_hit(runtime::Real3 first, runtime::Real3 second,
                                    const SurfaceTriangle &triangle,
                                    double coincidence,
                                    TriangleHit &result) noexcept {
  const runtime::Real3 direction = subtract(second, first);
  const runtime::Real3 edge1 =
      subtract(triangle.vertices_m[1], triangle.vertices_m[0]);
  const runtime::Real3 edge2 =
      subtract(triangle.vertices_m[2], triangle.vertices_m[0]);
  const runtime::Real3 p = cross(direction, edge2);
  const double determinant = dot(edge1, p);
  const double determinant_scale = norm(direction) * norm(edge1) * norm(edge2);
  const double determinant_tolerance =
      64.0 * std::numeric_limits<double>::epsilon() * determinant_scale;
  if (std::abs(determinant) <= determinant_tolerance) {
    const double first_distance =
        std::abs(dot(triangle.geometric_outward_normal,
                     subtract(first, triangle.vertices_m[0])));
    const double second_distance =
        std::abs(dot(triangle.geometric_outward_normal,
                     subtract(second, triangle.vertices_m[0])));
    return first_distance <= coincidence && second_distance <= coincidence
               ? HitDisposition::coplanar
               : HitDisposition::no_hit;
  }

  const double inverse = 1.0 / determinant;
  const runtime::Real3 offset = subtract(first, triangle.vertices_m[0]);
  const double u = dot(offset, p) * inverse;
  const runtime::Real3 q = cross(offset, edge1);
  const double v = dot(direction, q) * inverse;
  const double fraction = dot(edge2, q) * inverse;
  const double fraction_tolerance =
      coincidence / std::max(norm(direction), coincidence);
  const double barycentric_tolerance =
      coincidence /
      std::max({norm(edge1), norm(edge2), norm(direction), coincidence});
  if (u < -barycentric_tolerance || v < -barycentric_tolerance ||
      u + v > 1.0 + barycentric_tolerance || fraction < -fraction_tolerance ||
      fraction > 1.0 + fraction_tolerance) {
    return HitDisposition::no_hit;
  }
  result.triangle = triangle.id;
  result.fraction = std::clamp(fraction, 0.0, 1.0);
  result.point = add(first, multiply(direction, result.fraction));
  return HitDisposition::hit;
}

std::uint64_t query_fingerprint(const SurfaceQueryStorage &query) noexcept {
  std::uint64_t hash = fnv_offset;
  hash_u64(hash, UINT64_C(0x48554e44554e5133));
  hash_u64(hash, query.surface->fingerprint);
  hash_u64(hash, query.triangle_order.size());
  for (const TriangleId id : query.triangle_order) {
    hash_u64(hash, id);
  }
  hash_u64(hash, query.nodes.size());
  for (const BvhNode &node : query.nodes) {
    hash_u64(hash, double_bits(node.bounds.minimum.x));
    hash_u64(hash, double_bits(node.bounds.minimum.y));
    hash_u64(hash, double_bits(node.bounds.minimum.z));
    hash_u64(hash, double_bits(node.bounds.maximum.x));
    hash_u64(hash, double_bits(node.bounds.maximum.y));
    hash_u64(hash, double_bits(node.bounds.maximum.z));
    hash_u64(hash, node.left);
    hash_u64(hash, node.right);
    hash_u64(hash, node.begin);
    hash_u64(hash, node.end);
    hash_u64(hash, node.leaf ? 1U : 0U);
  }
  return hash;
}

} // namespace

namespace {

std::shared_ptr<const SurfaceQueryStorage>
build_surface_query_impl(std::shared_ptr<const SurfaceStorage> surface,
                         std::vector<TriangleId> initial_order) {
  if (!surface || surface->triangles.empty()) {
    fail("surface storage is empty");
  }
  if (initial_order.size() != surface->triangles.size()) {
    fail("initial BVH order has the wrong size");
  }
  auto sorted_order = initial_order;
  std::sort(sorted_order.begin(), sorted_order.end());
  for (std::size_t index = 0U; index < sorted_order.size(); ++index) {
    if (sorted_order[index] != index) {
      fail("initial BVH order is not a triangle-ID permutation");
    }
  }
  auto query = std::make_shared<SurfaceQueryStorage>();
  query->surface = std::move(surface);
  query->triangle_bounds.reserve(query->surface->triangles.size());
  query->triangle_order = std::move(initial_order);
  for (const SurfaceTriangle &triangle : query->surface->triangles) {
    query->triangle_bounds.push_back(triangle_box(triangle));
  }

  constexpr std::size_t leaf_size = 4U;
  const std::function<std::uint32_t(std::size_t, std::size_t)> build =
      [&](std::size_t begin, std::size_t end) -> std::uint32_t {
    if (query->nodes.size() >= std::numeric_limits<std::uint32_t>::max()) {
      fail("BVH node count exceeds uint32 range");
    }
    AxisAlignedBox bounds =
        query->triangle_bounds[query->triangle_order[begin]];
    for (std::size_t index = begin + 1U; index < end; ++index) {
      bounds =
          merge(bounds, query->triangle_bounds[query->triangle_order[index]]);
    }
    const auto node_index = static_cast<std::uint32_t>(query->nodes.size());
    query->nodes.push_back({bounds, 0U, 0U, begin, end, false});
    if (end - begin <= leaf_size) {
      query->nodes[node_index].leaf = true;
      return node_index;
    }

    const runtime::Real3 extent = subtract(bounds.maximum, bounds.minimum);
    std::size_t axis = 0U;
    if (extent.y > extent.x) {
      axis = 1U;
    }
    if (extent.z > coordinate(extent, axis)) {
      axis = 2U;
    }
    std::stable_sort(
        query->triangle_order.begin() + static_cast<std::ptrdiff_t>(begin),
        query->triangle_order.begin() + static_cast<std::ptrdiff_t>(end),
        [&](TriangleId first, TriangleId second) {
          const double first_key =
              coordinate(centroid(query->surface->triangles[first]), axis);
          const double second_key =
              coordinate(centroid(query->surface->triangles[second]), axis);
          return first_key < second_key ||
                 (first_key == second_key && first < second);
        });
    const std::size_t middle = begin + (end - begin) / 2U;
    const std::uint32_t left = build(begin, middle);
    const std::uint32_t right = build(middle, end);
    query->nodes[node_index].left = left;
    query->nodes[node_index].right = right;
    return node_index;
  };
  (void)build(0U, query->triangle_order.size());
  query->fingerprint = query_fingerprint(*query);
  return query;
}

} // namespace

std::shared_ptr<const SurfaceQueryStorage>
build_surface_query(std::shared_ptr<const SurfaceStorage> surface) {
  std::vector<TriangleId> initial_order;
  initial_order.reserve(surface ? surface->triangles.size() : 0U);
  if (surface) {
    for (const SurfaceTriangle &triangle : surface->triangles) {
      initial_order.push_back(triangle.id);
    }
  }
  return build_surface_query_impl(std::move(surface), std::move(initial_order));
}

#if defined(HUNDUN_IMMERSED_ENABLE_TEST_ACCESS)
std::shared_ptr<const SurfaceQueryStorage>
build_surface_query_with_initial_order(
    std::shared_ptr<const SurfaceStorage> surface,
    std::vector<TriangleId> initial_order) {
  return build_surface_query_impl(std::move(surface), std::move(initial_order));
}
#endif

std::vector<TriangleId> bounded_candidates(const SurfaceQueryStorage &query,
                                           runtime::Real3 minimum_m,
                                           runtime::Real3 maximum_m) {
  if (!finite(minimum_m) || !finite(maximum_m) || minimum_m.x > maximum_m.x ||
      minimum_m.y > maximum_m.y || minimum_m.z > maximum_m.z) {
    fail("candidate bounds are invalid");
  }
  const AxisAlignedBox requested{minimum_m, maximum_m};
  std::vector<TriangleId> result;
  std::vector<std::uint32_t> pending{0U};
  while (!pending.empty()) {
    const std::uint32_t node_index = pending.back();
    pending.pop_back();
    const BvhNode &node = query.nodes[node_index];
    if (!boxes_overlap(node.bounds, requested)) {
      continue;
    }
    if (node.leaf) {
      for (std::size_t index = static_cast<std::size_t>(node.begin);
           index < static_cast<std::size_t>(node.end); ++index) {
        const TriangleId triangle = query.triangle_order[index];
        if (boxes_overlap(query.triangle_bounds[triangle], requested)) {
          result.push_back(triangle);
        }
      }
    } else {
      pending.push_back(node.right);
      pending.push_back(node.left);
    }
  }
  std::sort(result.begin(), result.end());
  result.erase(std::unique(result.begin(), result.end()), result.end());
  return result;
}

bool require_consistent_parity(const std::array<bool, 3> &votes) {
  if (votes[0] != votes[1] || votes[0] != votes[2]) {
    fail("fixed ray parity disagrees");
  }
  return votes[0];
}

} // namespace hundun::immersed::detail

namespace hundun::immersed {

namespace {

void increment(std::atomic<std::uint64_t> &counter) {
  auto value = counter.load(std::memory_order_relaxed);
  for (;;) {
    if (value == std::numeric_limits<std::uint64_t>::max())
      throw runtime::Error("surface query performance counter would overflow");
    if (counter.compare_exchange_weak(value, value + 1U,
                                      std::memory_order_relaxed,
                                      std::memory_order_relaxed))
      return;
  }
}

} // namespace

SurfaceQuery SurfaceQuery::create(const ImmersedSurface &surface) {
  return SurfaceQuery(detail::build_surface_query(surface.storage_));
}

ClosestPointResult SurfaceQuery::closest_point(runtime::Real3 point_m) const {
  if (!detail::finite(point_m)) {
    throw runtime::Error("surface query: closest-point input is non-finite");
  }
  increment(storage_->closest_calls);

  struct PendingNode final {
    double distance_squared{};
    std::uint32_t node{};
  };
  const auto farther = [](const PendingNode &first, const PendingNode &second) {
    return first.distance_squared > second.distance_squared ||
           (first.distance_squared == second.distance_squared &&
            first.node > second.node);
  };
  std::priority_queue<PendingNode, std::vector<PendingNode>, decltype(farther)>
      pending{farther};
  pending.push(
      {detail::box_distance_squared(storage_->nodes.front().bounds, point_m),
       0U});

  ClosestPointResult result;
  result.squared_distance_m2 = std::numeric_limits<double>::infinity();
  result.triangle = std::numeric_limits<TriangleId>::max();
  while (!pending.empty()) {
    const PendingNode current = pending.top();
    pending.pop();
    if (current.distance_squared > result.squared_distance_m2) {
      continue;
    }
    const detail::BvhNode &node = storage_->nodes[current.node];
    if (node.leaf) {
      for (std::size_t index = static_cast<std::size_t>(node.begin);
           index < static_cast<std::size_t>(node.end); ++index) {
        const TriangleId id = storage_->triangle_order[index];
        const SurfaceTriangle &triangle = storage_->surface->triangles[id];
        const runtime::Real3 closest =
            detail::closest_on_triangle(point_m, triangle);
        const double distance_squared =
            detail::squared_norm(detail::subtract(point_m, closest));
        if (distance_squared < result.squared_distance_m2 ||
            (distance_squared == result.squared_distance_m2 &&
             id < result.triangle)) {
          result = {id, closest, triangle.geometric_outward_normal,
                    distance_squared};
        }
      }
      continue;
    }
    const double left_distance = detail::box_distance_squared(
        storage_->nodes[node.left].bounds, point_m);
    const double right_distance = detail::box_distance_squared(
        storage_->nodes[node.right].bounds, point_m);
    if (left_distance <= result.squared_distance_m2) {
      pending.push({left_distance, node.left});
    }
    if (right_distance <= result.squared_distance_m2) {
      pending.push({right_distance, node.right});
    }
  }
  return result;
}

std::vector<SegmentIntersection>
SurfaceQuery::segment_intersections(runtime::Real3 a_m,
                                    runtime::Real3 b_m) const {
  if (!detail::finite(a_m) || !detail::finite(b_m)) {
    throw runtime::Error("surface query: segment input is non-finite");
  }
  increment(storage_->segment_calls);
  const runtime::Real3 direction = detail::subtract(b_m, a_m);
  const double length = detail::norm(direction);
  if (!std::isfinite(length) ||
      length <= storage_->surface->intersection_coincidence_m) {
    throw runtime::Error("surface query: segment length is zero or ambiguous");
  }
  const runtime::Real3 minimum{std::min(a_m.x, b_m.x), std::min(a_m.y, b_m.y),
                               std::min(a_m.z, b_m.z)};
  const runtime::Real3 maximum{std::max(a_m.x, b_m.x), std::max(a_m.y, b_m.y),
                               std::max(a_m.z, b_m.z)};
  const auto candidates =
      detail::bounded_candidates(*storage_, minimum, maximum);
  std::vector<detail::TriangleHit> raw_hits;
  for (const TriangleId id : candidates) {
    detail::TriangleHit hit;
    const auto disposition = detail::segment_triangle_hit(
        a_m, b_m, storage_->surface->triangles[id],
        storage_->surface->intersection_coincidence_m, hit);
    if (disposition == detail::HitDisposition::coplanar) {
      throw runtime::Error(
          "surface query: segment is coplanar with a surface triangle");
    }
    if (disposition == detail::HitDisposition::hit) {
      raw_hits.push_back(hit);
    }
  }
  std::sort(
      raw_hits.begin(), raw_hits.end(),
      [](const detail::TriangleHit &first, const detail::TriangleHit &second) {
        return first.fraction < second.fraction ||
               (first.fraction == second.fraction &&
                first.triangle < second.triangle);
      });

  const double fraction_tolerance =
      storage_->surface->intersection_coincidence_m / length;
  std::vector<SegmentIntersection> result;
  for (std::size_t begin = 0U; begin < raw_hits.size();) {
    std::size_t end = begin + 1U;
    detail::TriangleHit owned = raw_hits[begin];
    while (end < raw_hits.size() &&
           std::abs(raw_hits[end].fraction - raw_hits[begin].fraction) <=
               fraction_tolerance &&
           detail::norm(
               detail::subtract(raw_hits[end].point, raw_hits[begin].point)) <=
               storage_->surface->intersection_coincidence_m) {
      if (raw_hits[end].triangle < owned.triangle) {
        owned = raw_hits[end];
      }
      ++end;
    }
    result.push_back({owned.triangle, owned.point, owned.fraction});
    begin = end;
  }
  return result;
}

CellRegion SurfaceQuery::classify(runtime::Real3 point_m,
                                  config::ImmersedFluidSide fluid_side) const {
  if (fluid_side != config::ImmersedFluidSide::outside &&
      fluid_side != config::ImmersedFluidSide::inside) {
    throw runtime::Error("surface query: fluid side is invalid");
  }
  if (!detail::finite(point_m)) {
    throw runtime::Error("surface query: classification input is non-finite");
  }
  const auto &minimum = storage_->surface->bounding_box_min;
  const auto &bounding_maximum = storage_->surface->bounding_box_max;
  const double coincidence = storage_->surface->intersection_coincidence_m;
  const bool separated_from_bounding_box =
      point_m.x < minimum.x - coincidence ||
      point_m.x > bounding_maximum.x + coincidence ||
      point_m.y < minimum.y - coincidence ||
      point_m.y > bounding_maximum.y + coincidence ||
      point_m.z < minimum.z - coincidence ||
      point_m.z > bounding_maximum.z + coincidence;
  if (separated_from_bounding_box) {
    return fluid_side == config::ImmersedFluidSide::outside ? CellRegion::fluid
                                                            : CellRegion::solid;
  }
  const ClosestPointResult closest = closest_point(point_m);
  if (closest.squared_distance_m2 <= coincidence * coincidence) {
    throw runtime::Error(
        "surface query: point lies on the classification boundary");
  }

  const std::array<runtime::Real3, 3> unnormalized{
      runtime::Real3{1.0, std::sqrt(2.0), std::sqrt(3.0)},
      runtime::Real3{std::sqrt(5.0), 1.0, std::sqrt(7.0)},
      runtime::Real3{std::sqrt(11.0), std::sqrt(13.0), 1.0}};
  std::array<bool, 3> inside{};
  for (std::size_t ray = 0U; ray < unnormalized.size(); ++ray) {
    const runtime::Real3 direction = detail::multiply(
        unnormalized[ray], 1.0 / detail::norm(unnormalized[ray]));
    double exit_fraction = 0.0;
    for (std::size_t axis = 0U; axis < 3U; ++axis) {
      const double direction_value = detail::coordinate(direction, axis);
      const double maximum =
          detail::coordinate(storage_->surface->bounding_box_max, axis);
      const double point_value = detail::coordinate(point_m, axis);
      exit_fraction =
          std::max(exit_fraction, (maximum - point_value) / direction_value);
    }
    const double ray_length = std::max(0.0, exit_fraction) +
                              2.0 * storage_->surface->reference_length_m;
    const runtime::Real3 endpoint =
        detail::add(point_m, detail::multiply(direction, ray_length));
    inside[ray] = segment_intersections(point_m, endpoint).size() % 2U == 1U;
  }
  const bool geometrically_inside = detail::require_consistent_parity(inside);

  const bool fluid = fluid_side == config::ImmersedFluidSide::outside
                         ? !geometrically_inside
                         : geometrically_inside;
  return fluid ? CellRegion::fluid : CellRegion::solid;
}

std::uint64_t SurfaceQuery::fingerprint() const noexcept {
  return storage_->fingerprint;
}

diagnostics::Stage3PerformanceCounters
SurfaceQuery::performance_counters() const noexcept {
  diagnostics::Stage3PerformanceCounters result;
  result.init_query_closest_calls =
      storage_->closest_calls.load(std::memory_order_relaxed);
  result.init_query_segment_calls =
      storage_->segment_calls.load(std::memory_order_relaxed);
  return result;
}

} // namespace hundun::immersed
