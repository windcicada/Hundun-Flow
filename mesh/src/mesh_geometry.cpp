// SPDX-License-Identifier: Apache-2.0

#include "hundun/mesh/mesh_geometry.hpp"

#include "geometry_arithmetic.hpp"
#include "hundun/runtime/error.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <new>
#include <stdexcept>
#include <utility>
#include <variant>
#include <vector>

namespace hundun::mesh {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kRadiansToDegrees = 180.0 / kPi;
constexpr double kClosureFactor = 256.0;

using Mapping = std::variant<UniformBoxMapping, AnalyticWarpedBoxMapping>;

runtime::Real3 add(runtime::Real3 lhs, runtime::Real3 rhs) noexcept {
  return {lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z};
}

runtime::Real3 subtract(runtime::Real3 lhs, runtime::Real3 rhs) noexcept {
  return {lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z};
}

runtime::Real3 multiply(runtime::Real3 value, double scalar) noexcept {
  return {value.x * scalar, value.y * scalar, value.z * scalar};
}

runtime::Real3 negate(runtime::Real3 value) noexcept {
  return {-value.x, -value.y, -value.z};
}

double dot(runtime::Real3 lhs, runtime::Real3 rhs) noexcept {
  return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

runtime::Real3 cross(runtime::Real3 lhs, runtime::Real3 rhs) noexcept {
  return {lhs.y * rhs.z - lhs.z * rhs.y, lhs.z * rhs.x - lhs.x * rhs.z,
          lhs.x * rhs.y - lhs.y * rhs.x};
}

double norm(runtime::Real3 value) noexcept {
  const double scale =
      std::max({std::abs(value.x), std::abs(value.y), std::abs(value.z)});
  if (scale == 0.0 || !std::isfinite(scale)) {
    return scale;
  }
  const double x = value.x / scale;
  const double y = value.y / scale;
  const double z = value.z / scale;
  return scale * std::sqrt(x * x + y * y + z * z);
}

bool finite(runtime::Real3 value) noexcept {
  return std::isfinite(value.x) && std::isfinite(value.y) &&
         std::isfinite(value.z);
}

double range_safe_determinant(runtime::Real3 first,
                              runtime::Real3 second,
                              runtime::Real3 third,
                              double factor = 1.0) noexcept {
  if (!finite(first) || !finite(second) || !finite(third) ||
      !std::isfinite(factor)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (factor == 0.0) {
    return 0.0;
  }

  const double x_scale =
      std::max({std::abs(first.x), std::abs(second.x), std::abs(third.x)});
  const double y_scale =
      std::max({std::abs(first.y), std::abs(second.y), std::abs(third.y)});
  const double z_scale =
      std::max({std::abs(first.z), std::abs(second.z), std::abs(third.z)});
  if (x_scale == 0.0 || y_scale == 0.0 || z_scale == 0.0) {
    return 0.0;
  }

  const runtime::Real3 normalized_first{first.x / x_scale,
                                         first.y / y_scale,
                                         first.z / z_scale};
  const runtime::Real3 normalized_second{second.x / x_scale,
                                          second.y / y_scale,
                                          second.z / z_scale};
  const runtime::Real3 normalized_third{third.x / x_scale,
                                         third.y / y_scale,
                                         third.z / z_scale};
  // Row normalization bounds every cross-product operand.  Restore the three
  // physical row scales by exponent so their intermediate product is never
  // formed in binary64.
  const double normalized_determinant =
      dot(normalized_first, cross(normalized_second, normalized_third));
  if (!std::isfinite(normalized_determinant) ||
      normalized_determinant == 0.0) {
    return normalized_determinant;
  }

  int x_exponent = 0;
  int y_exponent = 0;
  int z_exponent = 0;
  int determinant_exponent = 0;
  int factor_exponent = 0;
  const double mantissa =
      std::frexp(x_scale, &x_exponent) *
      std::frexp(y_scale, &y_exponent) *
      std::frexp(z_scale, &z_exponent) *
      std::frexp(normalized_determinant, &determinant_exponent) *
      std::frexp(factor, &factor_exponent);
  return std::scalbn(mantissa, x_exponent + y_exponent + z_exponent +
                                   determinant_exponent + factor_exponent);
}

bool positive(runtime::Real3 value) noexcept {
  return value.x > 0.0 && value.y > 0.0 && value.z > 0.0;
}

bool same(runtime::Int3 lhs, runtime::Int3 rhs) noexcept {
  return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
}

bool same(runtime::Box3 lhs, runtime::Box3 rhs) noexcept {
  return same(lhs.begin, rhs.begin) && same(lhs.end, rhs.end);
}

bool same(LogicalFace lhs, LogicalFace rhs) noexcept {
  return lhs.axis == rhs.axis && same(lhs.coordinate, rhs.coordinate);
}

void validate_box_parameters(runtime::Real3 origin, runtime::Real3 length) {
  if (!finite(origin)) {
    throw runtime::Error("box mapping origin components must be finite");
  }
  if (!finite(length) || !positive(length)) {
    throw runtime::Error("box mapping lengths must be finite and positive");
  }
  const runtime::Real3 upper = add(origin, length);
  if (!finite(upper)) {
    throw runtime::Error("box mapping upper endpoints must be finite");
  }
}

void validate_logical(runtime::Real3 logical) {
  if (!finite(logical) || logical.x < 0.0 || logical.x > 1.0 ||
      logical.y < 0.0 || logical.y > 1.0 || logical.z < 0.0 ||
      logical.z > 1.0) {
    throw runtime::Error("logical coordinate is outside [0, 1]^3");
  }
}

runtime::Real3 mapping_origin(const Mapping& mapping) {
  return std::visit([](const auto& value) { return value.origin_m(); },
                    mapping);
}

runtime::Real3 mapping_length(const Mapping& mapping) {
  return std::visit([](const auto& value) { return value.length_m(); },
                    mapping);
}

runtime::Real3 map_point(const Mapping& mapping, runtime::Real3 logical) {
  return std::visit([&](const auto& value) { return value.map(logical); },
                    mapping);
}

MappingJacobian map_jacobian(const Mapping& mapping, runtime::Real3 logical) {
  return std::visit([&](const auto& value) { return value.jacobian(logical); },
                    mapping);
}

bool is_affine_mapping(const Mapping& mapping) noexcept {
  if (std::holds_alternative<UniformBoxMapping>(mapping)) {
    return true;
  }
  const auto amplitude =
      std::get<AnalyticWarpedBoxMapping>(mapping).amplitude();
  return amplitude.x == 0.0 && amplitude.y == 0.0 && amplitude.z == 0.0;
}

runtime::Real3 logical_vertex(runtime::Int3 vertex,
                              runtime::Int3 extent) noexcept {
  return {static_cast<double>(vertex.x) / static_cast<double>(extent.x),
          static_cast<double>(vertex.y) / static_cast<double>(extent.y),
          static_cast<double>(vertex.z) / static_cast<double>(extent.z)};
}

std::array<runtime::Int3, 4> face_vertex_coordinates(LogicalFace face) {
  const auto c = face.coordinate;
  switch (face.axis) {
    case FaceAxis::x:
      return {runtime::Int3{c.x, c.y, c.z}, runtime::Int3{c.x, c.y + 1, c.z},
              runtime::Int3{c.x, c.y + 1, c.z + 1},
              runtime::Int3{c.x, c.y, c.z + 1}};
    case FaceAxis::y:
      return {runtime::Int3{c.x, c.y, c.z}, runtime::Int3{c.x, c.y, c.z + 1},
              runtime::Int3{c.x + 1, c.y, c.z + 1},
              runtime::Int3{c.x + 1, c.y, c.z}};
    case FaceAxis::z:
      return {runtime::Int3{c.x, c.y, c.z}, runtime::Int3{c.x + 1, c.y, c.z},
              runtime::Int3{c.x + 1, c.y + 1, c.z},
              runtime::Int3{c.x, c.y + 1, c.z}};
    default:
      throw runtime::Error("invalid face axis in mesh geometry");
  }
}

std::array<LogicalFace, 6> cell_faces(runtime::Int3 cell) {
  return {LogicalFace{FaceAxis::x, {cell.x, cell.y, cell.z}},
          LogicalFace{FaceAxis::x, {cell.x + 1, cell.y, cell.z}},
          LogicalFace{FaceAxis::y, {cell.x, cell.y, cell.z}},
          LogicalFace{FaceAxis::y, {cell.x, cell.y + 1, cell.z}},
          LogicalFace{FaceAxis::z, {cell.x, cell.y, cell.z}},
          LogicalFace{FaceAxis::z, {cell.x, cell.y, cell.z + 1}}};
}

bool is_logical_minimum(LogicalFace face) {
  switch (face.axis) {
    case FaceAxis::x:
      return face.coordinate.x == 0;
    case FaceAxis::y:
      return face.coordinate.y == 0;
    case FaceAxis::z:
      return face.coordinate.z == 0;
    default:
      throw runtime::Error("invalid face axis in mesh geometry");
  }
}

struct CellIdentity {
  GlobalCellId global_id{};
  runtime::Int3 global_cell{};
  EntityOwnership ownership{EntityOwnership::ghost};
};

struct FaceIdentity {
  GlobalFaceId global_id{};
  LogicalFace logical{};
  EntityOwnership ownership{EntityOwnership::ghost};
  LocalCellId owner{};
  std::optional<LocalCellId> neighbour{};
  std::optional<std::uint32_t> patch{};
  std::optional<GlobalFaceId> periodic_pair{};
};

struct CellMetrics {
  runtime::Real3 center{};
  double volume{};
  double minimum_jacobian{};
  runtime::Real3 closure{};
};

struct FaceMetrics {
  runtime::Real3 center{};
  runtime::Real3 owner_area{};
  double area{};
  double skewness{};
  double non_orthogonality_degrees{};
};

struct PolyhedralCell {
  runtime::Real3 center{};
  double volume{};
};

std::array<runtime::Real3, 4> mapped_face_vertices(const Mapping& mapping,
                                                   runtime::Int3 extent,
                                                   LogicalFace face) {
  const auto coordinates = face_vertex_coordinates(face);
  std::array<runtime::Real3, 4> vertices{};
  for (std::size_t index = 0; index < vertices.size(); ++index) {
    vertices[index] =
        map_point(mapping, logical_vertex(coordinates[index], extent));
    if (!finite(vertices[index])) {
      throw runtime::Error("mapped face vertex is non-finite");
    }
  }
  return vertices;
}

runtime::Real3 face_center(
    const std::array<runtime::Real3, 4>& vertices, bool affine) noexcept {
  const auto anchor = vertices[0];
  if (affine) {
    // For a parallelogram this is the same arithmetic mean, expressed without
    // tangential cancellation across extreme physical aspect ratios.
    return add(anchor, multiply(subtract(vertices[2], anchor), 0.5));
  }
  auto center = anchor;
  for (std::size_t index = 1; index < vertices.size(); ++index) {
    center = add(
        center, multiply(subtract(vertices[index], anchor), 0.25));
  }
  return center;
}

runtime::Real3 positive_axis_area(
    const std::array<runtime::Real3, 4>& vertices) noexcept {
  const auto first = multiply(cross(subtract(vertices[1], vertices[0]),
                                    subtract(vertices[2], vertices[0])),
                              0.5);
  const auto second = multiply(cross(subtract(vertices[2], vertices[0]),
                                     subtract(vertices[3], vertices[0])),
                               0.5);
  return add(first, second);
}

std::array<runtime::Real3, 8> mapped_cell_vertices(const Mapping& mapping,
                                                   runtime::Int3 extent,
                                                   runtime::Int3 cell) {
  const std::array<runtime::Int3, 8> coordinates{
      runtime::Int3{cell.x, cell.y, cell.z},
      runtime::Int3{cell.x + 1, cell.y, cell.z},
      runtime::Int3{cell.x + 1, cell.y + 1, cell.z},
      runtime::Int3{cell.x, cell.y + 1, cell.z},
      runtime::Int3{cell.x, cell.y, cell.z + 1},
      runtime::Int3{cell.x + 1, cell.y, cell.z + 1},
      runtime::Int3{cell.x + 1, cell.y + 1, cell.z + 1},
      runtime::Int3{cell.x, cell.y + 1, cell.z + 1}};
  std::array<runtime::Real3, 8> result{};
  for (std::size_t index = 0; index < result.size(); ++index) {
    result[index] =
        map_point(mapping, logical_vertex(coordinates[index], extent));
  }
  return result;
}

PolyhedralCell calculate_warped_cell(const Mapping& mapping,
                                     runtime::Int3 extent, runtime::Int3 cell) {
  const auto vertices = mapped_cell_vertices(mapping, extent, cell);
  const auto anchor = vertices[0];
  auto reference = anchor;
  const bool affine = is_affine_mapping(mapping);
  for (std::size_t index = 0; index < vertices.size(); ++index) {
    const auto vertex = vertices[index];
    if (!finite(vertex)) {
      throw runtime::Error("mapped cell vertex is non-finite");
    }
    if (!affine && index != 0) {
      reference = add(
          reference, multiply(subtract(vertex, anchor), 0.125));
    }
  }
  if (affine) {
    // The opposite-vertex midpoint equals the eight-vertex mean for this
    // affine hexahedron and stays aligned with the affine face-centre path.
    reference = add(anchor, multiply(subtract(vertices[6], anchor), 0.5));
  }
  if (!finite(reference)) {
    throw runtime::Error("mesh cell reference point is non-finite");
  }

  struct TetrahedronContribution {
    double signed_volume{};
    runtime::Real3 relative_centroid{};
  };
  std::array<TetrahedronContribution, 12> contributions{};
  std::size_t contribution_count = 0;
  double volume = 0.0;
  const auto accumulate_triangle = [&](runtime::Real3 a, runtime::Real3 b,
                                       runtime::Real3 c) {
    const auto relative_a = subtract(a, reference);
    const auto relative_b = subtract(b, reference);
    const auto relative_c = subtract(c, reference);
    const double signed_volume = range_safe_determinant(
        relative_a, relative_b, relative_c, 1.0 / 6.0);
    if (!std::isfinite(signed_volume)) {
      throw runtime::Error("cell tetrahedron volume is non-finite");
    }
    runtime::Real3 relative_centroid{};
    for (const auto relative :
         std::array<runtime::Real3, 3>{relative_a, relative_b, relative_c}) {
      relative_centroid =
          add(relative_centroid, multiply(relative, 0.25));
    }
    if (!finite(relative_a) || !finite(relative_b) || !finite(relative_c) ||
        !finite(relative_centroid) ||
        contribution_count >= contributions.size()) {
      throw runtime::Error("cell tetrahedron centroid is non-finite");
    }
    contributions[contribution_count++] =
        {signed_volume, relative_centroid};
    volume += signed_volume;
    if (!std::isfinite(volume)) {
      throw runtime::Error("cell tetrahedron volume sum is non-finite");
    }
  };

  const auto accumulate_quad = [&](LogicalFace face, bool reverse) {
    const auto quad = mapped_face_vertices(mapping, extent, face);
    if (reverse) {
      accumulate_triangle(quad[0], quad[3], quad[2]);
      accumulate_triangle(quad[0], quad[2], quad[1]);
    } else {
      accumulate_triangle(quad[0], quad[1], quad[2]);
      accumulate_triangle(quad[0], quad[2], quad[3]);
    }
  };

  const auto faces = cell_faces(cell);
  accumulate_quad(faces[0], true);
  accumulate_quad(faces[1], false);
  accumulate_quad(faces[2], true);
  accumulate_quad(faces[3], false);
  accumulate_quad(faces[4], true);
  accumulate_quad(faces[5], false);

  if (contribution_count != contributions.size() || !std::isfinite(volume) ||
      volume <= 0.0) {
    throw runtime::Error("mesh cell volume must be finite and positive");
  }
  if (affine) {
    return {reference, volume};
  }
  runtime::Real3 relative_center{};
  for (const auto contribution : contributions) {
    const double weight = contribution.signed_volume / volume;
    const auto term = multiply(contribution.relative_centroid, weight);
    if (!std::isfinite(weight) || !finite(term)) {
      throw runtime::Error("mesh cell relative centroid is non-finite");
    }
    relative_center = add(relative_center, term);
    if (!finite(relative_center)) {
      throw runtime::Error("mesh cell relative centroid sum is non-finite");
    }
  }
  const runtime::Real3 center = add(reference, relative_center);
  if (!finite(center)) {
    throw runtime::Error("mesh cell centre must be finite");
  }
  return {center, volume};
}

double minimum_sampled_jacobian(const Mapping& mapping, runtime::Int3 extent,
                                runtime::Int3 cell) {
  double minimum = std::numeric_limits<double>::infinity();
  const auto sample = [&](double x, double y, double z) {
    const runtime::Real3 logical{
        (static_cast<double>(cell.x) + x) / static_cast<double>(extent.x),
        (static_cast<double>(cell.y) + y) / static_cast<double>(extent.y),
        (static_cast<double>(cell.z) + z) / static_cast<double>(extent.z)};
    const double determinant = map_jacobian(mapping, logical).determinant_m3();
    if (!std::isfinite(determinant) || determinant <= 0.0) {
      throw runtime::Error(
          "mesh mapping Jacobian must be finite and strictly positive");
    }
    minimum = std::min(minimum, determinant);
  };

  for (int z = 0; z <= 1; ++z) {
    for (int y = 0; y <= 1; ++y) {
      for (int x = 0; x <= 1; ++x) {
        sample(static_cast<double>(x), static_cast<double>(y),
               static_cast<double>(z));
      }
    }
  }
  sample(0.5, 0.5, 0.5);
  const double offset = 1.0 / (2.0 * std::sqrt(3.0));
  for (int z = -1; z <= 1; z += 2) {
    for (int y = -1; y <= 1; y += 2) {
      for (int x = -1; x <= 1; x += 2) {
        sample(0.5 + static_cast<double>(x) * offset,
               0.5 + static_cast<double>(y) * offset,
               0.5 + static_cast<double>(z) * offset);
      }
    }
  }
  return minimum;
}

runtime::Real3 periodic_displacement(const FaceIdentity& face,
                                     runtime::Real3 length) {
  runtime::Real3 shift{};
  if (!face.periodic_pair.has_value()) {
    return shift;
  }
  const bool minimum = is_logical_minimum(face.logical);
  switch (face.logical.axis) {
    case FaceAxis::x:
      shift.x = minimum ? -length.x : length.x;
      break;
    case FaceAxis::y:
      shift.y = minimum ? -length.y : length.y;
      break;
    case FaceAxis::z:
      shift.z = minimum ? -length.z : length.z;
      break;
    default:
      throw runtime::Error("invalid face axis in mesh geometry");
  }
  return shift;
}

template <class Function>
auto translate_allocation_failures(Function&& function) {
  try {
    return function();
  } catch (const runtime::Error&) {
    throw;
  } catch (const std::bad_alloc&) {
    throw runtime::Error("mesh geometry allocation failed");
  } catch (const std::length_error&) {
    throw runtime::Error("mesh geometry allocation size is unsupported");
  }
}

}  // namespace

double MappingJacobian::determinant_m3() const noexcept {
  return range_safe_determinant(d_xi_m, d_eta_m, d_zeta_m);
}

UniformBoxMapping::UniformBoxMapping(runtime::Real3 origin_m,
                                     runtime::Real3 length_m)
    : origin_m_(origin_m), length_m_(length_m) {
  validate_box_parameters(origin_m_, length_m_);
}

runtime::Real3 UniformBoxMapping::origin_m() const noexcept {
  return origin_m_;
}

runtime::Real3 UniformBoxMapping::length_m() const noexcept {
  return length_m_;
}

runtime::Real3 UniformBoxMapping::map(runtime::Real3 logical) const {
  validate_logical(logical);
  const runtime::Real3 mapped{origin_m_.x + length_m_.x * logical.x,
                              origin_m_.y + length_m_.y * logical.y,
                              origin_m_.z + length_m_.z * logical.z};
  if (!finite(mapped)) {
    throw runtime::Error("mapped coordinate is non-finite");
  }
  return mapped;
}

MappingJacobian UniformBoxMapping::jacobian(runtime::Real3 logical) const {
  validate_logical(logical);
  return {{length_m_.x, 0.0, 0.0},
          {0.0, length_m_.y, 0.0},
          {0.0, 0.0, length_m_.z}};
}

AnalyticWarpedBoxMapping::AnalyticWarpedBoxMapping(runtime::Real3 origin_m,
                                                   runtime::Real3 length_m,
                                                   runtime::Real3 amplitude)
    : origin_m_(origin_m), length_m_(length_m), amplitude_(amplitude) {
  validate_box_parameters(origin_m_, length_m_);
  if (!finite(amplitude_) || std::abs(amplitude_.x) > 0.02 ||
      std::abs(amplitude_.y) > 0.02 || std::abs(amplitude_.z) > 0.02) {
    throw runtime::Error(
        "analytic warped-box amplitude must be finite and "
        "within [-0.02, 0.02]");
  }
}

runtime::Real3 AnalyticWarpedBoxMapping::origin_m() const noexcept {
  return origin_m_;
}

runtime::Real3 AnalyticWarpedBoxMapping::length_m() const noexcept {
  return length_m_;
}

runtime::Real3 AnalyticWarpedBoxMapping::amplitude() const noexcept {
  return amplitude_;
}

runtime::Real3 AnalyticWarpedBoxMapping::map(runtime::Real3 logical) const {
  validate_logical(logical);
  const double sin_x = std::sin(kPi * logical.x);
  const double sin_y = std::sin(kPi * logical.y);
  const double sin_z = std::sin(kPi * logical.z);
  const runtime::Real3 mapped{
      origin_m_.x +
          length_m_.x *
              (logical.x +
               amplitude_.x * std::sin(2.0 * kPi * logical.x) * sin_y * sin_z),
      origin_m_.y + length_m_.y * (logical.y +
                                   amplitude_.y * sin_x *
                                       std::sin(2.0 * kPi * logical.y) * sin_z),
      origin_m_.z +
          length_m_.z * (logical.z + amplitude_.z * sin_x * sin_y *
                                         std::sin(2.0 * kPi * logical.z))};
  if (!finite(mapped)) {
    throw runtime::Error("mapped coordinate is non-finite");
  }
  return mapped;
}

MappingJacobian AnalyticWarpedBoxMapping::jacobian(
    runtime::Real3 logical) const {
  validate_logical(logical);
  const double sin_x = std::sin(kPi * logical.x);
  const double sin_y = std::sin(kPi * logical.y);
  const double sin_z = std::sin(kPi * logical.z);
  const double cos_x = std::cos(kPi * logical.x);
  const double cos_y = std::cos(kPi * logical.y);
  const double cos_z = std::cos(kPi * logical.z);
  const double sin_2x = std::sin(2.0 * kPi * logical.x);
  const double sin_2y = std::sin(2.0 * kPi * logical.y);
  const double sin_2z = std::sin(2.0 * kPi * logical.z);
  const double cos_2x = std::cos(2.0 * kPi * logical.x);
  const double cos_2y = std::cos(2.0 * kPi * logical.y);
  const double cos_2z = std::cos(2.0 * kPi * logical.z);

  MappingJacobian result{};
  result.d_xi_m = {
      length_m_.x * (1.0 + 2.0 * kPi * amplitude_.x * cos_2x * sin_y * sin_z),
      length_m_.y *
          (kPi * amplitude_.y * cos_x * sin_2y * sin_z),
      length_m_.z *
          (kPi * amplitude_.z * cos_x * sin_y * sin_2z)};
  result.d_eta_m = {
      length_m_.x *
          (kPi * amplitude_.x * sin_2x * cos_y * sin_z),
      length_m_.y * (1.0 + 2.0 * kPi * amplitude_.y * sin_x * cos_2y * sin_z),
      length_m_.z *
          (kPi * amplitude_.z * sin_x * cos_y * sin_2z)};
  result.d_zeta_m = {
      length_m_.x *
          (kPi * amplitude_.x * sin_2x * sin_y * cos_z),
      length_m_.y *
          (kPi * amplitude_.y * sin_x * sin_2y * cos_z),
      length_m_.z * (1.0 + 2.0 * kPi * amplitude_.z * sin_x * sin_y * cos_2z)};
  if (!finite(result.d_xi_m) || !finite(result.d_eta_m) ||
      !finite(result.d_zeta_m)) {
    throw runtime::Error("mapping Jacobian is non-finite");
  }
  return result;
}

struct MeshGeometry::Impl {
  MappingKind kind{MappingKind::uniform_box};
  runtime::Int3 extent{};
  runtime::Box3 owned_box{};
  runtime::Real3 origin{};
  runtime::Real3 length{};
  std::optional<runtime::Real3> spacing{};
  std::size_t owned_cell_count{};
  std::vector<CellIdentity> cell_identity{};
  std::vector<FaceIdentity> face_identity{};
  std::vector<CellMetrics> cells{};
  std::vector<FaceMetrics> faces{};
  Mapping mapping;

  Impl(const MeshTopology& topology, Mapping supplied_mapping,
       MappingKind mapping_kind)
      : kind(mapping_kind),
        extent(topology.global_extent()),
        owned_box(topology.owned_global_box()),
        origin(mapping_origin(supplied_mapping)),
        length(mapping_length(supplied_mapping)),
        owned_cell_count(topology.owned_cell_count()),
        mapping(std::move(supplied_mapping)) {
    if (extent.x <= 0 || extent.y <= 0 || extent.z <= 0) {
      throw runtime::Error("mesh geometry requires a positive global extent");
    }
    cell_identity.reserve(topology.local_cell_count());
    cells.reserve(topology.local_cell_count());
    face_identity.reserve(topology.local_face_count());
    faces.reserve(topology.local_face_count());

    if (kind == MappingKind::uniform_box) {
      const runtime::Real3 uniform_spacing{
          length.x / static_cast<double>(extent.x),
          length.y / static_cast<double>(extent.y),
          length.z / static_cast<double>(extent.z)};
      if (!finite(uniform_spacing) || !positive(uniform_spacing)) {
        throw runtime::Error("mesh spacing must be finite and positive");
      }
      spacing = uniform_spacing;
    }

    const double uniform_volume =
        spacing.has_value() ? detail::range_safe_product(*spacing) : 0.0;
    if (spacing.has_value() &&
        (!std::isfinite(uniform_volume) || uniform_volume <= 0.0)) {
      throw runtime::Error("mesh cell volume must be finite and positive");
    }

    for (LocalCellId local = 0; local < topology.local_cell_count(); ++local) {
      const auto global_cell = topology.global_cell(local);
      cell_identity.push_back({topology.global_cell_id(local), global_cell,
                               topology.cell_ownership(local)});
      CellMetrics metrics{};
      if (spacing.has_value()) {
        metrics.center = {
            origin.x + (static_cast<double>(global_cell.x) + 0.5) * spacing->x,
            origin.y + (static_cast<double>(global_cell.y) + 0.5) * spacing->y,
            origin.z + (static_cast<double>(global_cell.z) + 0.5) * spacing->z};
        metrics.volume = uniform_volume;
      } else {
        const auto polyhedron =
            calculate_warped_cell(mapping, extent, global_cell);
        metrics.center = polyhedron.center;
        metrics.volume = polyhedron.volume;
      }
      if (!finite(metrics.center) || !std::isfinite(metrics.volume) ||
          metrics.volume <= 0.0) {
        throw runtime::Error("mesh cell geometry is non-finite or degenerate");
      }
      metrics.minimum_jacobian =
          minimum_sampled_jacobian(mapping, extent, global_cell);
      cells.push_back(metrics);
    }

    for (LocalFaceId local = 0; local < topology.local_face_count(); ++local) {
      FaceIdentity identity{
          topology.global_face_id(local), topology.logical_face(local),
          topology.face_ownership(local), topology.owner(local),
          topology.neighbour(local),      topology.patch_id(local),
          topology.periodic_pair(local)};
      const auto vertices =
          mapped_face_vertices(mapping, extent, identity.logical);
      FaceMetrics metrics{};
      metrics.center = face_center(vertices, is_affine_mapping(mapping));
      const auto positive_area = positive_axis_area(vertices);
      metrics.owner_area = is_logical_minimum(identity.logical)
                               ? negate(positive_area)
                               : positive_area;
      metrics.area = norm(metrics.owner_area);
      if (!finite(metrics.center)) {
        throw runtime::Error("mesh face centre is non-finite");
      }
      if (!finite(metrics.owner_area)) {
        throw runtime::Error("mesh face area vector is non-finite");
      }
      if (!std::isfinite(metrics.area)) {
        throw runtime::Error("mesh face area magnitude is non-finite");
      }
      if (metrics.area <= 0.0) {
        throw runtime::Error("mesh face area is degenerate");
      }
      face_identity.push_back(identity);
      faces.push_back(metrics);
    }

    for (LocalCellId local = 0; local < owned_cell_count; ++local) {
      runtime::Real3 closure{};
      std::array<double, 6> area_magnitudes{};
      double area_scale = 0.0;
      std::size_t area_index = 0;
      for (const auto logical : cell_faces(cell_identity[local].global_cell)) {
        const auto face =
            topology.find_local_face(topology.global_face_id(logical));
        if (!face.has_value()) {
          throw runtime::Error(
              "owned cell face is missing from local topology");
        }
        runtime::Real3 outward{};
        if (face_identity[*face].owner == local) {
          outward = faces[*face].owner_area;
        } else if (face_identity[*face].neighbour.has_value() &&
                   *face_identity[*face].neighbour == local) {
          outward = negate(faces[*face].owner_area);
        } else {
          throw runtime::Error("owned cell face connectivity is inconsistent");
        }
        closure = add(closure, outward);
        const double area = faces[*face].area;
        if (!std::isfinite(area) || area <= 0.0 ||
            area_index >= area_magnitudes.size()) {
          throw runtime::Error("mesh cell face area is invalid");
        }
        area_magnitudes[area_index++] = area;
        area_scale = std::max(area_scale, area);
      }
      if (!finite(closure) || area_index != area_magnitudes.size() ||
          !std::isfinite(area_scale) || area_scale <= 0.0) {
        throw runtime::Error("mesh cell face-area closure check failed");
      }

      double scaled_area_sum = 0.0;
      for (const double area : area_magnitudes) {
        const double normalized_area = area / area_scale;
        // A positive area may underflow relative to an extreme finite scale;
        // zero then has no representable contribution to the scaled sum.
        if (!std::isfinite(normalized_area) || normalized_area < 0.0) {
          throw runtime::Error("mesh cell scaled face area is invalid");
        }
        scaled_area_sum += normalized_area;
        if (!std::isfinite(scaled_area_sum) || scaled_area_sum < 0.0) {
          throw runtime::Error("mesh cell scaled face-area sum is invalid");
        }
      }
      if (scaled_area_sum <= 0.0) {
        throw runtime::Error("mesh cell scaled face-area sum is invalid");
      }
      const runtime::Real3 scaled_closure{
          closure.x / area_scale, closure.y / area_scale,
          closure.z / area_scale};
      const double scaled_closure_norm = norm(scaled_closure);
      if (!finite(scaled_closure) || !std::isfinite(scaled_closure_norm) ||
          scaled_closure_norm > kClosureFactor *
                                    std::numeric_limits<double>::epsilon() *
                                    scaled_area_sum) {
        throw runtime::Error("mesh cell face-area closure check failed");
      }
      cells[local].closure = closure;
    }

    for (LocalFaceId local = 0; local < faces.size(); ++local) {
      auto& metrics = faces[local];
      const auto& identity = face_identity[local];
      const auto owner_center = cells[identity.owner].center;
      runtime::Real3 displacement{};
      const bool physical_boundary = !identity.neighbour.has_value();
      if (physical_boundary) {
        displacement = subtract(metrics.center, owner_center);
      } else {
        displacement =
            add(subtract(cells[*identity.neighbour].center, owner_center),
                periodic_displacement(identity, length));
      }
      const double displacement_norm = norm(displacement);
      const double projection = dot(metrics.owner_area, displacement);
      if (!finite(displacement) || !std::isfinite(displacement_norm) ||
          displacement_norm <= 0.0 || !std::isfinite(projection) ||
          projection <= 0.0) {
        throw runtime::Error("mesh face centre displacement is invalid");
      }
      const runtime::Real3 unit_area{
          metrics.owner_area.x / metrics.area,
          metrics.owner_area.y / metrics.area,
          metrics.owner_area.z / metrics.area};
      const runtime::Real3 unit_displacement{
          displacement.x / displacement_norm,
          displacement.y / displacement_norm,
          displacement.z / displacement_norm};
      if (!finite(unit_area) || !finite(unit_displacement)) {
        throw runtime::Error("mesh face normalized direction is invalid");
      }
      const double normalized_projection = dot(unit_area, unit_displacement);
      if (!std::isfinite(normalized_projection)) {
        throw runtime::Error("mesh face normalized projection is invalid");
      }
      const double cosine = std::clamp(normalized_projection, -1.0, 1.0);
      metrics.non_orthogonality_degrees =
          spacing.has_value() ? 0.0 : std::acos(cosine) * kRadiansToDegrees;
      if (physical_boundary || spacing.has_value()) {
        metrics.skewness = 0.0;
      } else {
        const double interpolation =
            dot(metrics.owner_area, subtract(metrics.center, owner_center)) /
            projection;
        const auto intersection =
            add(owner_center, multiply(displacement, interpolation));
        metrics.skewness =
            norm(subtract(intersection, metrics.center)) / displacement_norm;
      }
      if (!std::isfinite(metrics.non_orthogonality_degrees) ||
          metrics.non_orthogonality_degrees < 0.0 ||
          metrics.non_orthogonality_degrees >= 90.0 ||
          !std::isfinite(metrics.skewness) || metrics.skewness < 0.0) {
        throw runtime::Error("mesh face quality metric is invalid");
      }
    }
  }
};

std::unique_ptr<MeshGeometry::Impl> MeshGeometry::create_impl(
    const MeshTopology& topology, UniformBoxMapping mapping) {
  return translate_allocation_failures([&] {
    return std::make_unique<Impl>(topology, Mapping{std::move(mapping)},
                                  MappingKind::uniform_box);
  });
}

std::unique_ptr<MeshGeometry::Impl> MeshGeometry::create_impl(
    const MeshTopology& topology, AnalyticWarpedBoxMapping mapping) {
  return translate_allocation_failures([&] {
    return std::make_unique<Impl>(topology, Mapping{std::move(mapping)},
                                  MappingKind::analytic_warped_box);
  });
}

MeshGeometry::MeshGeometry(const MeshTopology& topology,
                           UniformBoxMapping mapping)
    : impl_(create_impl(topology, std::move(mapping))) {}

MeshGeometry::MeshGeometry(const MeshTopology& topology,
                           AnalyticWarpedBoxMapping mapping)
    : impl_(create_impl(topology, std::move(mapping))) {}

MeshGeometry::~MeshGeometry() = default;
MeshGeometry::MeshGeometry(MeshGeometry&&) noexcept = default;
MeshGeometry& MeshGeometry::operator=(MeshGeometry&&) noexcept = default;

MappingKind MeshGeometry::mapping_kind() const noexcept { return impl_->kind; }

runtime::Int3 MeshGeometry::global_extent() const noexcept {
  return impl_->extent;
}

runtime::Box3 MeshGeometry::owned_global_box() const noexcept {
  return impl_->owned_box;
}

runtime::Real3 MeshGeometry::origin_m() const noexcept { return impl_->origin; }

runtime::Real3 MeshGeometry::length_m() const noexcept { return impl_->length; }

std::optional<runtime::Real3> MeshGeometry::uniform_spacing_m() const noexcept {
  return impl_->spacing;
}

bool MeshGeometry::compatible(const MeshTopology& topology) const {
  if (!same(impl_->extent, topology.global_extent()) ||
      !same(impl_->owned_box, topology.owned_global_box()) ||
      impl_->owned_cell_count != topology.owned_cell_count() ||
      impl_->cell_identity.size() != topology.local_cell_count() ||
      impl_->face_identity.size() != topology.local_face_count()) {
    return false;
  }
  for (LocalCellId local = 0; local < impl_->cell_identity.size(); ++local) {
    const auto& identity = impl_->cell_identity[local];
    if (identity.global_id != topology.global_cell_id(local) ||
        !same(identity.global_cell, topology.global_cell(local)) ||
        identity.ownership != topology.cell_ownership(local)) {
      return false;
    }
  }
  for (LocalFaceId local = 0; local < impl_->face_identity.size(); ++local) {
    const auto& identity = impl_->face_identity[local];
    if (identity.global_id != topology.global_face_id(local) ||
        !same(identity.logical, topology.logical_face(local)) ||
        identity.ownership != topology.face_ownership(local) ||
        identity.owner != topology.owner(local) ||
        identity.neighbour != topology.neighbour(local) ||
        identity.patch != topology.patch_id(local) ||
        identity.periodic_pair != topology.periodic_pair(local)) {
      return false;
    }
  }
  return true;
}

void MeshGeometry::require_compatible(const MeshTopology& topology) const {
  if (!compatible(topology)) {
    throw runtime::Error("mesh geometry is incompatible with mesh topology");
  }
}

runtime::Real3 MeshGeometry::vertex_position_m(
    runtime::Int3 global_vertex) const {
  if (global_vertex.x < 0 || global_vertex.x > impl_->extent.x ||
      global_vertex.y < 0 || global_vertex.y > impl_->extent.y ||
      global_vertex.z < 0 || global_vertex.z > impl_->extent.z) {
    throw runtime::Error("global vertex is outside the mesh lattice");
  }
  return map_point(impl_->mapping,
                   logical_vertex(global_vertex, impl_->extent));
}

template <class Value>
const Value& require_local(const std::vector<Value>& values, std::size_t local,
                           const char* subject) {
  if (local >= values.size()) {
    throw runtime::Error(subject);
  }
  return values[local];
}

runtime::Real3 MeshGeometry::cell_center_m(LocalCellId local) const {
  return require_local(impl_->cells, local,
                       "local cell index is outside mesh geometry")
      .center;
}

double MeshGeometry::cell_volume_m3(LocalCellId local) const {
  return require_local(impl_->cells, local,
                       "local cell index is outside mesh geometry")
      .volume;
}

double MeshGeometry::minimum_jacobian_determinant_m3(LocalCellId local) const {
  return require_local(impl_->cells, local,
                       "local cell index is outside mesh geometry")
      .minimum_jacobian;
}

runtime::Real3 MeshGeometry::face_center_m(LocalFaceId local) const {
  return require_local(impl_->faces, local,
                       "local face index is outside mesh geometry")
      .center;
}

runtime::Real3 MeshGeometry::face_area_vector_m2(LocalFaceId local,
                                                 FaceSide side) const {
  const auto& metrics = require_local(
      impl_->faces, local, "local face index is outside mesh geometry");
  switch (side) {
    case FaceSide::owner:
      return metrics.owner_area;
    case FaceSide::neighbour:
      if (!impl_->face_identity[local].neighbour.has_value()) {
        throw runtime::Error("physical boundary face has no neighbour side");
      }
      return negate(metrics.owner_area);
    default:
      throw runtime::Error("invalid face side");
  }
}

double MeshGeometry::face_area_m2(LocalFaceId local) const {
  return require_local(impl_->faces, local,
                       "local face index is outside mesh geometry")
      .area;
}

double MeshGeometry::face_skewness(LocalFaceId local) const {
  return require_local(impl_->faces, local,
                       "local face index is outside mesh geometry")
      .skewness;
}

double MeshGeometry::face_non_orthogonality_degrees(LocalFaceId local) const {
  return require_local(impl_->faces, local,
                       "local face index is outside mesh geometry")
      .non_orthogonality_degrees;
}

runtime::Real3 MeshGeometry::cell_closure_m2(LocalCellId local) const {
  if (local >= impl_->owned_cell_count) {
    throw runtime::Error("cell closure is available for owned cells only");
  }
  return impl_->cells[local].closure;
}

}  // namespace hundun::mesh
