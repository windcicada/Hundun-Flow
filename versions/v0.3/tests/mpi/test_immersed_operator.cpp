// SPDX-License-Identifier: Apache-2.0

#include "hundun/fvm_immersed_operator.hpp"

#include "src/fvm_immersed_boundary_authority_detail.hpp"
#include "tests/support/fvm_immersed_operator_test_access.hpp"
#include "hundun/bc_basic_boundary.hpp"
#include "hundun/cfg_resolved_case.hpp"
#include "hundun/cfg_resolved_case_v3.hpp"
#include "hundun/fvm_immersed_reconstruction.hpp"
#include "hundun/ib_ghost_stencil_plan.hpp"
#include "hundun/ib_local_flow_pattern.hpp"
#include "hundun/mesh_geometry.hpp"
#include "hundun/rt_exchange_plan.hpp"
#include "hundun/rt_field_descriptor.hpp"
#include "hundun/rt_field_registry.hpp"
#include "hundun/rt_field_storage.hpp"
#include "hundun/rt_halo_exchange.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "hundun/rt_structured_decomposition.hpp"
#include "src/ib_deterministic_qr_detail.hpp"
#include "src/ib_quadratic_reconstruction_detail.hpp"
#include "tests/support/stage3_stl_fixture.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <string>
#include <system_error>
#include <tuple>
#include <type_traits>
#include <vector>

namespace {

using namespace hundun;

constexpr runtime::Int3 kExtent{12, 12, 12};
constexpr runtime::PhaseId kPhase = 81U;
constexpr runtime::ActorId kActor = 23U;

std::uint64_t bits(double value) {
  std::uint64_t result{};
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

runtime::Real3 add(runtime::Real3 lhs, runtime::Real3 rhs) {
  return {lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z};
}

runtime::Real3 subtract(runtime::Real3 lhs, runtime::Real3 rhs) {
  return {lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z};
}

runtime::Real3 multiply(runtime::Real3 value, double scale) {
  return {value.x * scale, value.y * scale, value.z * scale};
}

double dot(runtime::Real3 lhs, runtime::Real3 rhs) {
  return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

double max_norm(runtime::Real3 value) {
  return std::max({std::abs(value.x), std::abs(value.y), std::abs(value.z)});
}

void prove_periodic_checkerboard_operator_selection() {
  constexpr int extent = 4;
  double checkerboard_norm_squared = 0.0;
  double constant_difference_squared = 0.0;
  double compact_checkerboard_energy = 0.0;
  double compact_constant_action_squared = 0.0;
  for (int k = 0; k < extent; ++k)
    for (int j = 0; j < extent; ++j)
      for (int i = 0; i < extent; ++i) {
        const auto parity = [](int x, int y, int z) {
          return (x + y + z) % 2 == 0 ? 1.0 : -1.0;
        };
        const auto wrap = [=](int value) {
          return (value + extent) % extent;
        };
        const double center = parity(i, j, k);
        constant_difference_squared += (center - 1.0) * (center - 1.0);
        const std::array<double, 3> centered_gradient{
            0.5 * (parity(wrap(i + 1), j, k) -
                   parity(wrap(i - 1), j, k)),
            0.5 * (parity(i, wrap(j + 1), k) -
                   parity(i, wrap(j - 1), k)),
            0.5 * (parity(i, j, wrap(k + 1)) -
                   parity(i, j, wrap(k - 1)))};
        for (const double value : centered_gradient)
          checkerboard_norm_squared += value * value;
        const double compact_action =
            6.0 * center - parity(wrap(i + 1), j, k) -
            parity(wrap(i - 1), j, k) - parity(i, wrap(j + 1), k) -
            parity(i, wrap(j - 1), k) - parity(i, j, wrap(k + 1)) -
            parity(i, j, wrap(k - 1));
        compact_checkerboard_energy += center * compact_action;
        const double compact_constant_action = 6.0 - 1.0 - 1.0 - 1.0 -
                                               1.0 - 1.0 - 1.0;
        compact_constant_action_squared +=
            compact_constant_action * compact_constant_action;
      }
  HUNDUN_CHECK(constant_difference_squared > 0.0);
  HUNDUN_CHECK(bits(checkerboard_norm_squared) == bits(0.0));
  HUNDUN_CHECK(compact_checkerboard_energy > 0.0);
  HUNDUN_CHECK(bits(compact_constant_action_squared) == bits(0.0));
}

  using SnapshotDonorKey =
    std::pair<mesh::GlobalCellId, std::uint32_t>;
using SnapshotWallKey =
    std::pair<immersed::ImmersedLinkId, std::uint32_t>;

template <class Term, class KeyFunction>
bool snapshot_terms_to_map(const std::vector<Term> &terms,
                           KeyFunction key_function,
                           std::map<decltype(key_function(terms.front())),
                                    double> &result) {
  result.clear();
  for (const auto &term : terms) {
    if (!std::isfinite(term.coefficient) ||
        !result.emplace(key_function(term), term.coefficient).second)
      return false;
  }
  return true;
}

template <class Key>
std::map<Key, double> subtract_snapshot_maps(const std::map<Key, double> &lhs,
                                             const std::map<Key, double> &rhs) {
  auto result = lhs;
  for (const auto &[key, value] : rhs)
    result[key] -= value;
  for (auto item = result.begin(); item != result.end();) {
    if (item->second == 0.0)
      item = result.erase(item);
    else
      ++item;
  }
  return result;
}

template <class Key>
bool snapshot_maps_near(const std::map<Key, double> &left,
                        const std::map<Key, double> &right) {
  if (left.size() != right.size())
    return false;
  auto lhs = left.begin();
  auto rhs = right.begin();
  for (; lhs != left.end(); ++lhs, ++rhs) {
    if (lhs->first != rhs->first)
      return false;
    const double tolerance =
        64.0 * std::numeric_limits<double>::epsilon() *
        std::max({1.0, std::abs(lhs->second), std::abs(rhs->second)});
    if (std::abs(lhs->second - rhs->second) > tolerance)
      return false;
  }
  return true;
}

template <class Key>
bool snapshot_maps_bitwise_equal(const std::map<Key, double> &left,
                                 const std::map<Key, double> &right) {
  if (left.size() != right.size())
    return false;
  auto lhs = left.begin();
  auto rhs = right.begin();
  for (; lhs != left.end(); ++lhs, ++rhs)
    if (lhs->first != rhs->first || bits(lhs->second) != bits(rhs->second))
      return false;
  return true;
}

 bool interface_pressure_snapshot_is_consistent(
    const finite_volume::test::InterfacePressureRowSnapshot &row) {
  if (row.authority_fingerprint == 0U || row.a22_donor_terms.empty() ||
      row.background_donor_terms.empty() ||
      row.difference_donor_terms.empty() || row.a22_wall_terms.empty() ||
      row.difference_wall_terms.empty())
    return false;
  const auto donor_key = [](const auto &term) {
    return SnapshotDonorKey{term.pressure_cell, term.output_component};
  };
  const auto wall_key = [](const auto &term) {
    return SnapshotWallKey{term.link, term.output_component};
  };
  std::map<SnapshotDonorKey, double> background;
  std::map<SnapshotDonorKey, double> a22;
  std::map<SnapshotDonorKey, double> legacy;
  std::map<SnapshotDonorKey, double> difference;
  std::map<SnapshotWallKey, double> background_wall;
  std::map<SnapshotWallKey, double> a22_wall;
  std::map<SnapshotWallKey, double> legacy_wall;
  std::map<SnapshotWallKey, double> difference_wall;
  if (!snapshot_terms_to_map(row.background_donor_terms, donor_key,
                             background) ||
      !snapshot_terms_to_map(row.a22_donor_terms, donor_key, a22) ||
      !snapshot_terms_to_map(row.legacy_unconstrained_lfp_donor_terms,
                             donor_key, legacy) ||
      !snapshot_terms_to_map(row.difference_donor_terms, donor_key,
                             difference) ||
      !snapshot_terms_to_map(row.background_wall_terms, wall_key,
                             background_wall) ||
      !snapshot_terms_to_map(row.a22_wall_terms, wall_key, a22_wall) ||
      !snapshot_terms_to_map(row.legacy_unconstrained_lfp_wall_terms,
                             wall_key, legacy_wall) ||
      !snapshot_terms_to_map(row.difference_wall_terms, wall_key,
                             difference_wall))
    return false;
  if (!snapshot_maps_near(subtract_snapshot_maps(a22, background),
                          difference) ||
      !snapshot_maps_near(
          subtract_snapshot_maps(a22_wall, background_wall),
          difference_wall))
    return false;
  for (std::uint32_t component = 0U; component < 3U; ++component) {
    double sum = 0.0;
    double scale = 0.0;
    for (const auto &[key, coefficient] : difference) {
      if (key.second != component)
        continue;
      sum += coefficient;
      scale += std::abs(coefficient);
    }
    const double tolerance =
        64.0 * std::numeric_limits<double>::epsilon() *
        std::max(1.0, scale);
    if (std::abs(sum) > tolerance)
      return false;
  }
  return true;
}

void prove_interface_pressure_difference_snapshot(
    const finite_volume::ImmersedOperatorAdapter &adapter) {
  const auto product_rows =
      finite_volume::test::ImmersedOperatorTestAccess::rows(adapter);
  const auto rows = finite_volume::test::ImmersedOperatorTestAccess::
      interface_pressure_rows(adapter);
  const auto expected_count = static_cast<std::size_t>(std::count_if(
      product_rows.begin(), product_rows.end(),
      [](const auto &row) { return !row.links.empty(); }));
  HUNDUN_CHECK(!rows.empty());
  HUNDUN_CHECK(rows.size() == expected_count);
  bool distinguished_from_legacy = false;
  for (const auto &row : rows) {
    HUNDUN_CHECK(interface_pressure_snapshot_is_consistent(row));
    HUNDUN_CHECK(row.background_wall_terms.empty());
    HUNDUN_CHECK(row.legacy_unconstrained_lfp_wall_terms.empty());

    const auto donor_key = [](const auto &term) {
      return SnapshotDonorKey{term.pressure_cell, term.output_component};
    };
    std::map<SnapshotDonorKey, double> a22;
    std::map<SnapshotDonorKey, double> background;
    std::map<SnapshotDonorKey, double> legacy;
    std::map<SnapshotDonorKey, double> difference;
    HUNDUN_CHECK(snapshot_terms_to_map(row.a22_donor_terms, donor_key, a22));
    HUNDUN_CHECK(snapshot_terms_to_map(row.background_donor_terms, donor_key,
                                       background));
    HUNDUN_CHECK(snapshot_terms_to_map(
        row.legacy_unconstrained_lfp_donor_terms, donor_key, legacy));
    HUNDUN_CHECK(snapshot_terms_to_map(row.difference_donor_terms, donor_key,
                                       difference));
    distinguished_from_legacy =
        distinguished_from_legacy ||
        !snapshot_maps_near(subtract_snapshot_maps(a22, legacy), difference);

    auto coefficient_mutation = row;
    coefficient_mutation.difference_donor_terms.front().coefficient +=
        0.125 * std::max(
                    1.0, std::abs(coefficient_mutation.difference_donor_terms
                                      .front()
                                      .coefficient));
    HUNDUN_CHECK(
        !interface_pressure_snapshot_is_consistent(coefficient_mutation));

    auto nested_size_mutation = row;
    nested_size_mutation.difference_donor_terms.pop_back();
    HUNDUN_CHECK(
        !interface_pressure_snapshot_is_consistent(nested_size_mutation));

    auto background_mutation = row;
    background_mutation.background_donor_terms.front().coefficient +=
        0.25 * std::max(
                   1.0, std::abs(background_mutation.background_donor_terms
                                     .front()
                                     .coefficient));
    HUNDUN_CHECK(!interface_pressure_snapshot_is_consistent(
        background_mutation));

    auto donor_identity_mutation = row;
    donor_identity_mutation.difference_donor_terms.front().pressure_cell ^=
        1U;
    HUNDUN_CHECK(!interface_pressure_snapshot_is_consistent(
        donor_identity_mutation));

    auto wall_mutation = row;
    wall_mutation.difference_wall_terms.front().coefficient +=
        0.125 * std::max(
                    1.0, std::abs(wall_mutation.difference_wall_terms.front()
                                      .coefficient));
    HUNDUN_CHECK(!interface_pressure_snapshot_is_consistent(wall_mutation));
  }
  HUNDUN_CHECK(distinguished_from_legacy);
}

 runtime::Real3 cross(runtime::Real3 lhs, runtime::Real3 rhs) {
  return {lhs.y * rhs.z - lhs.z * rhs.y, lhs.z * rhs.x - lhs.x * rhs.z,
          lhs.x * rhs.y - lhs.y * rhs.x};
}

double quadratic_average(runtime::Int3 cell,
                         const mesh::MeshGeometry &geometry) {
  constexpr std::array<runtime::Int3, 8> offsets{
      runtime::Int3{0, 0, 0}, runtime::Int3{1, 0, 0}, runtime::Int3{1, 1, 0},
      runtime::Int3{0, 1, 0}, runtime::Int3{0, 0, 1}, runtime::Int3{1, 0, 1},
      runtime::Int3{1, 1, 1}, runtime::Int3{0, 1, 1}};
  std::array<runtime::Real3, 8> vertices{};
  runtime::Real3 reference{};
  for (std::size_t index = 0U; index < vertices.size(); ++index) {
    vertices[index] = geometry.vertex_position_m({cell.x + offsets[index].x,
                                                  cell.y + offsets[index].y,
                                                  cell.z + offsets[index].z});
    reference = add(reference, multiply(vertices[index], 0.125));
  }
  constexpr std::array<std::array<std::size_t, 3>, 12> triangles{{
      {{0, 4, 7}},
      {{0, 7, 3}},
      {{1, 2, 6}},
      {{1, 6, 5}},
      {{0, 1, 5}},
      {{0, 5, 4}},
      {{3, 7, 6}},
      {{3, 6, 2}},
      {{0, 3, 2}},
      {{0, 2, 1}},
      {{4, 5, 6}},
      {{4, 6, 7}},
  }};
  double integral = 0.0;
  double volume = 0.0;
  for (const auto triangle : triangles) {
    const std::array<runtime::Real3, 4> tetra{reference, vertices[triangle[0]],
                                              vertices[triangle[1]],
                                              vertices[triangle[2]]};
    const double tetra_volume =
        dot(subtract(tetra[1], tetra[0]),
            cross(subtract(tetra[2], tetra[0]), subtract(tetra[3], tetra[0]))) /
        6.0;
    HUNDUN_CHECK(tetra_volume > 0.0);
    runtime::Real3 mean{};
    for (const auto point : tetra)
      mean = add(mean, multiply(point, 0.25));
    const auto product_mean = [&](int lhs, int rhs) {
      double lhs_sum = 0.0;
      double rhs_sum = 0.0;
      double diagonal = 0.0;
      for (const auto point : tetra) {
        const std::array<double, 3> values{point.x, point.y, point.z};
        lhs_sum += values[static_cast<std::size_t>(lhs)];
        rhs_sum += values[static_cast<std::size_t>(rhs)];
        diagonal += values[static_cast<std::size_t>(lhs)] *
                    values[static_cast<std::size_t>(rhs)];
      }
      return (lhs_sum * rhs_sum + diagonal) / 20.0;
    };
    double value = 2.0 + 0.4 * mean.x - 0.7 * mean.y + 0.3 * mean.z +
                   0.2 * product_mean(0, 0) - 0.1 * product_mean(0, 1) +
                   0.15 * product_mean(1, 1) + 0.12 * product_mean(2, 2);
    integral += tetra_volume * value;
    volume += tetra_volume;
  }
  return integral / volume;
}

double quadratic_value(runtime::Real3 point) {
  return 2.0 + 0.4 * point.x - 0.7 * point.y + 0.3 * point.z +
         0.2 * point.x * point.x - 0.1 * point.x * point.y +
         0.15 * point.y * point.y + 0.12 * point.z * point.z;
}

runtime::Real3 quadratic_gradient(runtime::Real3 point) {
  return {0.4 + 0.4 * point.x - 0.1 * point.y,
          -0.7 - 0.1 * point.x + 0.3 * point.y, 0.3 + 0.24 * point.z};
}

std::array<runtime::Int3, 4> logical_face_vertices(mesh::LogicalFace face) {
  const auto coordinate = face.coordinate;
  switch (face.axis) {
  case mesh::FaceAxis::x:
    return {runtime::Int3{coordinate.x, coordinate.y, coordinate.z},
            runtime::Int3{coordinate.x, coordinate.y + 1, coordinate.z},
            runtime::Int3{coordinate.x, coordinate.y + 1, coordinate.z + 1},
            runtime::Int3{coordinate.x, coordinate.y, coordinate.z + 1}};
  case mesh::FaceAxis::y:
    return {runtime::Int3{coordinate.x, coordinate.y, coordinate.z},
            runtime::Int3{coordinate.x, coordinate.y, coordinate.z + 1},
            runtime::Int3{coordinate.x + 1, coordinate.y, coordinate.z + 1},
            runtime::Int3{coordinate.x + 1, coordinate.y, coordinate.z}};
  case mesh::FaceAxis::z:
    return {runtime::Int3{coordinate.x, coordinate.y, coordinate.z},
            runtime::Int3{coordinate.x + 1, coordinate.y, coordinate.z},
            runtime::Int3{coordinate.x + 1, coordinate.y + 1, coordinate.z},
            runtime::Int3{coordinate.x, coordinate.y + 1, coordinate.z}};
  }
  throw runtime::Error("invalid immersed-operator face axis");
}

double quadratic_triangle_average(runtime::Real3 a, runtime::Real3 b,
                                  runtime::Real3 c) {
  constexpr std::array<std::array<double, 3>, 3> barycentric{{
      {{2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0}},
      {{1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0}},
      {{1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0}},
  }};
  double result = 0.0;
  for (const auto &weight : barycentric)
    result += quadratic_value(
        add(add(multiply(a, weight[0]), multiply(b, weight[1])),
            multiply(c, weight[2])));
  return result / 3.0;
}

double uniform_face_quadratic_average(const mesh::MeshTopology &topology,
                                      const mesh::MeshGeometry &geometry,
                                      mesh::LocalFaceId face) {
  const auto logical = logical_face_vertices(topology.logical_face(face));
  std::array<runtime::Real3, 4> vertices{};
  for (std::size_t index = 0U; index < vertices.size(); ++index)
    vertices[index] = geometry.vertex_position_m(logical[index]);
  constexpr std::array<std::array<std::size_t, 3>, 2> triangles{{
      {{0U, 1U, 2U}},
      {{0U, 2U, 3U}},
  }};
  double integral = 0.0;
  double area = 0.0;
  for (const auto &triangle : triangles) {
    const auto &a = vertices[triangle[0]];
    const auto &b = vertices[triangle[1]];
    const auto &c = vertices[triangle[2]];
    const double triangle_area =
        0.5 * std::sqrt(dot(cross(subtract(b, a), subtract(c, a)),
                            cross(subtract(b, a), subtract(c, a))));
    HUNDUN_CHECK(triangle_area > 0.0 && std::isfinite(triangle_area));
    integral += triangle_area * quadratic_triangle_average(a, b, c);
    area += triangle_area;
  }
  return integral / area;
}

using Matrix7 = std::array<std::array<double, 7>, 7>;

std::array<double, 7> multiply(const Matrix7 &matrix,
                               const std::array<double, 7> &values) {
  std::array<double, 7> result{};
  for (std::size_t row = 0U; row < result.size(); ++row)
    for (std::size_t column = 0U; column < values.size(); ++column)
      result[row] += matrix[row][column] * values[column];
  return result;
}

runtime::Real3 normalized(runtime::Real3 value) {
  const double magnitude = std::sqrt(dot(value, value));
  HUNDUN_CHECK(magnitude > 0.0 && std::isfinite(magnitude));
  return multiply(value, 1.0 / magnitude);
}

std::array<double, 3> independent_rotation_factors(runtime::Real3 normal) {
  constexpr double pi = 3.141592653589793238462643383279502884;
  const std::array<double, 3> absolute{std::abs(normal.x), std::abs(normal.y),
                                       std::abs(normal.z)};
  std::size_t axis = 0U;
  if (absolute[1] < absolute[axis])
    axis = 1U;
  if (absolute[2] < absolute[axis])
    axis = 2U;
  const std::array<runtime::Real3, 3> axes{runtime::Real3{1.0, 0.0, 0.0},
                                           runtime::Real3{0.0, 1.0, 0.0},
                                           runtime::Real3{0.0, 0.0, 1.0}};
  const auto tangent_one = normalized(cross(axes[axis], normal));
  const auto tangent_two = cross(normal, tangent_one);
  const double pitch = std::asin(std::clamp(-tangent_one.z, -1.0, 1.0));
  const double roll = std::atan2(tangent_two.z, normal.z);
  const double yaw = std::atan2(tangent_one.y, tangent_one.x);
  const auto factor = [=](double angle) {
    double folded = std::abs(std::remainder(angle, pi));
    if (folded > 0.5 * pi)
      folded = pi - folded;
    return std::clamp(1.0 - 2.0 * folded / pi, 0.0, 1.0);
  };
  return {factor(roll), factor(pitch), factor(yaw)};
}

immersed::LocalCoefficientRow
independent_paper_map(immersed::LocalCoefficientRow row, double k0,
                      runtime::Real3 solid_to_fluid_normal) {
  const auto normal = normalized(solid_to_fluid_normal);
  if (normal.x < 0.0)
    std::swap(row.neighbour[2], row.neighbour[3]);
  if (normal.y < 0.0)
    std::swap(row.neighbour[0], row.neighbour[1]);
  if (normal.z < 0.0)
    std::swap(row.neighbour[4], row.neighbour[5]);
  const auto factors = independent_rotation_factors(
      {std::abs(normal.x), std::abs(normal.y), std::abs(normal.z)});
  const double k1 = factors[0];
  const double k2 = factors[1];
  const double k3 = factors[2];
  const std::array<double, 7> input{
      row.neighbour[0], row.neighbour[1], row.neighbour[2], row.neighbour[3],
      row.neighbour[4], row.neighbour[5], row.diagonal};
  Matrix7 equation_16{};
  equation_16[0][0] = 1.0;
  equation_16[1][0] = 1.0;
  equation_16[2][2] = k1;
  equation_16[2][5] = 1.0 - k1;
  equation_16[3][3] = k1;
  equation_16[3][4] = 1.0 - k1;
  equation_16[4][4] = k1;
  equation_16[4][2] = 1.0 - k1;
  equation_16[5][5] = k1;
  equation_16[5][3] = 1.0 - k1;
  equation_16[6][6] = 1.0;
  Matrix7 equation_15{};
  equation_15[0][0] = k2;
  equation_15[0][5] = 1.0 - k2;
  equation_15[1][1] = k2;
  equation_15[1][4] = 1.0 - k2;
  equation_15[2][2] = 1.0;
  equation_15[3][3] = 1.0;
  equation_15[4][4] = k2;
  equation_15[4][0] = 1.0 - k2;
  equation_15[5][5] = k2;
  equation_15[5][1] = 1.0 - k2;
  equation_15[6][6] = 1.0;
  Matrix7 equation_14{};
  equation_14[0][0] = k3;
  equation_14[0][2] = 1.0 - k3;
  equation_14[1][1] = k3;
  equation_14[1][3] = 1.0 - k3;
  equation_14[2][2] = k3;
  equation_14[2][1] = 1.0 - k3;
  equation_14[3][3] = k3;
  equation_14[3][0] = 1.0 - k3;
  equation_14[4][4] = 1.0;
  equation_14[5][5] = 1.0;
  equation_14[6][6] = 1.0;
  const auto rotated = multiply(
      equation_14, multiply(equation_15, multiply(equation_16, input)));
  Matrix7 equation_13{};
  for (std::size_t neighbour = 0U; neighbour < 6U; ++neighbour) {
    equation_13[neighbour][neighbour] = k0;
    equation_13[6][neighbour] = 1.0 - k0;
  }
  equation_13[6][6] = 1.0;
  const auto transformed = multiply(equation_13, rotated);
  immersed::LocalCoefficientRow result{};
  for (std::size_t index = 0U; index < 6U; ++index)
    result.neighbour[index] = transformed[index];
  result.diagonal = transformed[6];
  result.source = row.source;
  if (normal.x < 0.0)
    std::swap(result.neighbour[2], result.neighbour[3]);
  if (normal.y < 0.0)
    std::swap(result.neighbour[0], result.neighbour[1]);
  if (normal.z < 0.0)
    std::swap(result.neighbour[4], result.neighbour[5]);
  return result;
}

runtime::Int3 logical_cell(mesh::GlobalCellId id, runtime::Int3 extent) {
  const auto nx = static_cast<std::uint64_t>(extent.x);
  const auto ny = static_cast<std::uint64_t>(extent.y);
  const auto plane = nx * ny;
  return {static_cast<int>(id % nx), static_cast<int>((id / nx) % ny),
          static_cast<int>(id / plane)};
}

runtime::Int3 neighbour_offset(std::size_t occurrence) {
  constexpr std::array<runtime::Int3, 6> offsets{
      runtime::Int3{-1, 0, 0}, runtime::Int3{1, 0, 0},  runtime::Int3{0, -1, 0},
      runtime::Int3{0, 1, 0},  runtime::Int3{0, 0, -1}, runtime::Int3{0, 0, 1}};
  HUNDUN_CHECK(occurrence < offsets.size());
  return offsets[occurrence];
}

runtime::Int3 process_grid(int ranks) {
  if (ranks == 1)
    return {1, 1, 1};
  if (ranks == 2)
    return {2, 1, 1};
  return {2, 2, 1};
}

config::FlowCaseConfig case_config(int ranks, runtime::Int3 grid, bool warped) {
  config::FlowCaseConfig config{};
  config.schema_version = 2;
  config.case_name = "task8-immersed-operator";
  config.simulation_type = config::SimulationType::variable_density_flow;
  config.density_model = config::DensityModel::constant;
  config.resources.expected_ranks = ranks;
  config.resources.process_grid = grid;
  config.mesh.cells = kExtent;
  config.mesh.origin_m = {0.0, 0.0, 0.0};
  config.mesh.length_m = {1.0, 1.0, 1.0};
  config.mesh.mapping = warped ? config::MeshMapping::analytic_warped_box
                               : config::MeshMapping::uniform_box;
  if (warped)
    config.mesh.warp_amplitude = runtime::Real3{0.02, -0.015, 0.01};
  config.time.mode = config::TimeMode::fixed;
  config.time.steps = 1;
  config.time.initial_dt_s = 0.01;
  config.time.min_dt_s = 0.01;
  config.time.max_dt_s = 0.01;
  config.physics.rho_ref_kg_per_m3 = 1.0;
  config.physics.dynamic_viscosity_pa_s = 0.01;
  config.physics.inlet_consistency_rtol = 1.0e-12;
  config.scalars.push_back(config::FlowScalarConfig{"alpha", 0.02});
  constexpr std::array<config::PatchName, 6> names{
      config::PatchName::x_min, config::PatchName::x_max,
      config::PatchName::y_min, config::PatchName::y_max,
      config::PatchName::z_min, config::PatchName::z_max};
  for (std::size_t patch = 0U; patch < names.size(); ++patch) {
    config.boundaries[patch].patch = names[patch];
    config.boundaries[patch].type = config::BoundaryType::no_slip_wall;
  }
  return config;
}

test::StlFixtureTriangle sphere_triangle(runtime::Real3 first,
                                         runtime::Real3 second,
                                         runtime::Real3 third,
                                         runtime::Real3 center) {
  auto normal = cross(subtract(second, first), subtract(third, first));
  const auto centroid =
      multiply(add(add(first, second), third), 1.0 / 3.0);
  if (dot(normal, subtract(centroid, center)) < 0.0) {
    std::swap(second, third);
    normal = cross(subtract(second, first), subtract(third, first));
  }
  const double length = std::sqrt(dot(normal, normal));
  HUNDUN_CHECK(length > 0.0);
  return {multiply(normal, 1.0 / length), {first, second, third}};
}

runtime::Real3 project_to_sphere(runtime::Real3 point,
                                 runtime::Real3 center, double radius) {
  const auto radial = subtract(point, center);
  const double length = std::sqrt(dot(radial, radial));
  HUNDUN_CHECK(length > 0.0);
  return add(center, multiply(radial, radius / length));
}

std::vector<test::StlFixtureTriangle> wall_authority_sphere() {
  constexpr runtime::Real3 center{1.25, 1.25, 1.25};
  constexpr double radius = 0.45;
  const runtime::Real3 px{center.x + radius, center.y, center.z};
  const runtime::Real3 nx{center.x - radius, center.y, center.z};
  const runtime::Real3 py{center.x, center.y + radius, center.z};
  const runtime::Real3 ny{center.x, center.y - radius, center.z};
  const runtime::Real3 pz{center.x, center.y, center.z + radius};
  const runtime::Real3 nz{center.x, center.y, center.z - radius};
  std::vector<test::StlFixtureTriangle> triangles{
      sphere_triangle(pz, px, py, center), sphere_triangle(pz, py, nx, center),
      sphere_triangle(pz, nx, ny, center), sphere_triangle(pz, ny, px, center),
      sphere_triangle(nz, py, px, center), sphere_triangle(nz, nx, py, center),
      sphere_triangle(nz, ny, nx, center), sphere_triangle(nz, px, ny, center)};
  for (int level = 0; level < 2; ++level) {
    std::vector<test::StlFixtureTriangle> refined;
    refined.reserve(4U * triangles.size());
    for (const auto &triangle : triangles) {
      const auto ab = project_to_sphere(
          multiply(add(triangle.vertices[0], triangle.vertices[1]), 0.5),
          center, radius);
      const auto bc = project_to_sphere(
          multiply(add(triangle.vertices[1], triangle.vertices[2]), 0.5),
          center, radius);
      const auto ca = project_to_sphere(
          multiply(add(triangle.vertices[2], triangle.vertices[0]), 0.5),
          center, radius);
      refined.push_back(
          sphere_triangle(triangle.vertices[0], ab, ca, center));
      refined.push_back(
          sphere_triangle(ab, triangle.vertices[1], bc, center));
      refined.push_back(
          sphere_triangle(ca, bc, triangle.vertices[2], center));
      refined.push_back(sphere_triangle(ab, bc, ca, center));
    }
    triangles = std::move(refined);
  }
  return triangles;
}

runtime::Real3 midpoint(runtime::Real3 first, runtime::Real3 second) {
  return {(first.x + second.x) * 0.5, (first.y + second.y) * 0.5,
          (first.z + second.z) * 0.5};
}

std::vector<test::StlFixtureTriangle>
partition_conforming_cube(runtime::Real3 translation) {
  auto refined = test::translated(test::outward_cube(), translation);
  for (unsigned level = 0U; level < 3U; ++level) {
    std::vector<test::StlFixtureTriangle> next;
    next.reserve(4U * refined.size());
    for (const auto &triangle : refined) {
      const auto ab = midpoint(triangle.vertices[0], triangle.vertices[1]);
      const auto bc = midpoint(triangle.vertices[1], triangle.vertices[2]);
      const auto ca = midpoint(triangle.vertices[2], triangle.vertices[0]);
      next.push_back({triangle.file_normal, {triangle.vertices[0], ab, ca}});
      next.push_back({triangle.file_normal, {ab, triangle.vertices[1], bc}});
      next.push_back({triangle.file_normal, {ca, bc, triangle.vertices[2]}});
      next.push_back({triangle.file_normal, {ab, bc, ca}});
    }
    refined = std::move(next);
  }
  return refined;
}

class FixtureFile final {
public:
  explicit FixtureFile(const runtime::MpiContext &mpi, bool permuted = false,
                       bool wall_authority_fixture = false)
      : mpi_(&mpi) {
    std::string text;
    if (mpi.rank() == 0) {
      path_ =
          std::filesystem::temp_directory_path() /
          ("hundun-task8-" +
           std::to_string(
               std::chrono::steady_clock::now().time_since_epoch().count()) +
           ".stl");
      auto triangles = wall_authority_fixture
                           ? wall_authority_sphere()
                           : partition_conforming_cube({0.75, 0.75, 0.75});
      if (permuted) {
        std::reverse(triangles.begin(), triangles.end());
      }
      if (wall_authority_fixture) {
        test::write_text(path_, test::ascii_stl(
                                    triangles, "task11-wall-authority-sphere"));
        text = path_.string();
      } else {
        constexpr double cosine_z = 0.90630778703664996;
        constexpr double sine_z = 0.42261826174069944;
        constexpr double cosine_y = 0.96592582628906831;
        constexpr double sine_y = 0.25881904510252074;
        for (auto &triangle : triangles) {
          for (auto &vertex : triangle.vertices) {
            const runtime::Real3 centered{vertex.x - 1.25, vertex.y - 1.25,
                                          vertex.z - 1.25};
            const runtime::Real3 z_rotated{
                cosine_z * centered.x - sine_z * centered.y,
                sine_z * centered.x + cosine_z * centered.y, centered.z};
            const runtime::Real3 yz_rotated{
                cosine_y * z_rotated.x + sine_y * z_rotated.z, z_rotated.y,
                -sine_y * z_rotated.x + cosine_y * z_rotated.z};
            vertex = add(yz_rotated, {1.25, 1.25, 1.25});
          }
          const auto first =
              subtract(triangle.vertices[1], triangle.vertices[0]);
          const auto second =
              subtract(triangle.vertices[2], triangle.vertices[0]);
          auto normal = cross(first, second);
          const double magnitude = std::sqrt(dot(normal, normal));
          HUNDUN_CHECK(magnitude > 0.0 && std::isfinite(magnitude));
          triangle.file_normal = multiply(normal, 1.0 / magnitude);
        }
        test::write_text(
            path_, test::ascii_stl(triangles, "task8-rotated-cube"));
        text = path_.string();
      }
    }
    std::uint64_t size = text.size();
    HUNDUN_CHECK(MPI_Bcast(&size, 1, MPI_UINT64_T, 0, mpi.comm()) ==
                 MPI_SUCCESS);
    text.resize(static_cast<std::size_t>(size));
    HUNDUN_CHECK(MPI_Bcast(text.data(), static_cast<int>(text.size()), MPI_BYTE,
                           0, mpi.comm()) == MPI_SUCCESS);
    path_ = text;
    HUNDUN_CHECK(MPI_Barrier(mpi.comm()) == MPI_SUCCESS);
  }
  ~FixtureFile() {
    MPI_Barrier(mpi_->comm());
    if (mpi_->rank() == 0) {
      std::error_code error;
      std::filesystem::remove(path_, error);
    }
  }
  const std::filesystem::path &path() const noexcept { return path_; }

private:
  const runtime::MpiContext *mpi_;
  std::filesystem::path path_;
};

runtime::FieldDescriptor cell_descriptor(const char *name,
                                         std::uint32_t components, int ghost) {
  return {name,
          "1",
          "task8",
          runtime::FunctionSpace::cell_average,
          runtime::ScalarType::float64,
          components,
          ghost,
          false,
          runtime::RestartPolicy::transient,
          runtime::OutputPolicy::never};
}

runtime::FieldDescriptor face_descriptor(const char *name,
                                         std::uint32_t components) {
  return {name,
          "1",
          "task8",
          runtime::FunctionSpace::face_value,
          runtime::ScalarType::float64,
          components,
          0,
          false,
          runtime::RestartPolicy::transient,
          runtime::OutputPolicy::never};
}

struct Fields final {
  runtime::FieldId velocity{};
  runtime::FieldId pressure{};
  runtime::FieldId scalar{};
  runtime::FieldId velocity_gradient{};
  runtime::FieldId scalar_gradient{};
  runtime::FieldId face_velocity{};
  runtime::FieldId face_scalar{};
  runtime::FieldId viscosity{};
  runtime::FieldId gamma{};
  runtime::FieldId mass_flux{};
  runtime::FieldId momentum_residual{};
  runtime::FieldId scalar_residual{};
};

Fields declare_fields(runtime::FieldRegistry &registry, int ghost) {
  Fields result{};
  result.velocity =
      registry.declare_field(cell_descriptor("velocity", 3U, ghost));
  result.pressure =
      registry.declare_field(cell_descriptor("pressure", 1U, ghost));
  result.scalar = registry.declare_field(cell_descriptor("scalar", 1U, ghost));
  result.velocity_gradient =
      registry.declare_field(cell_descriptor("grad_u", 9U, ghost));
  result.scalar_gradient =
      registry.declare_field(cell_descriptor("grad_q", 3U, ghost));
  result.face_velocity = registry.declare_field(face_descriptor("face_u", 3U));
  result.face_scalar = registry.declare_field(face_descriptor("face_q", 1U));
  result.viscosity = registry.declare_field(face_descriptor("mu", 1U));
  result.gamma = registry.declare_field(face_descriptor("gamma", 1U));
  result.mass_flux = finite_volume::declare_face_mass_flux(registry);
  result.momentum_residual =
      registry.declare_field(cell_descriptor("momentum_r", 3U, 0));
  result.scalar_residual =
      registry.declare_field(cell_descriptor("scalar_r", 1U, 0));
  return result;
}

finite_volume::ReconstructionFieldBinding
binding(runtime::FieldStorage &storage, const runtime::FieldAccessPlan &access,
        runtime::FieldId field) {
  return {storage, access, kPhase, kActor, field};
}

std::vector<std::uint64_t>
owned_active_bits(const mesh::MeshTopology &topology,
                  const immersed::ImmersedDomain &domain,
                  const runtime::FieldView<const double> &field,
                  std::uint32_t components, runtime::Box3 box) {
  std::vector<std::uint64_t> result;
  for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count();
       ++cell) {
    if (domain.region(cell) != immersed::CellRegion::fluid)
      continue;
    const auto global = topology.global_cell(cell);
    for (std::uint32_t component = 0U; component < components; ++component)
      result.push_back(
          bits(field(global.x - box.begin.x, global.y - box.begin.y,
                     global.z - box.begin.z, static_cast<int>(component))));
  }
  return result;
}

void fill_zero(runtime::FieldView<double> view) {
  const auto extent = view.interior_extent();
  for (int k = 0; k < extent.z; ++k)
    for (int j = 0; j < extent.y; ++j)
      for (int i = 0; i < extent.x; ++i)
        for (std::uint32_t component = 0U; component < view.components();
             ++component)
          view(i, j, k, static_cast<int>(component)) = 0.0;
}

bool report_bits_equal(const finite_volume::ImmersedOperatorReport &lhs,
                       const finite_volume::ImmersedOperatorReport &rhs) {
  if (lhs.active_row_count != rhs.active_row_count ||
      lhs.replacement_group_count != rhs.replacement_group_count ||
      lhs.simultaneous_substitution_count !=
          rhs.simultaneous_substitution_count)
    return false;
  const auto parts_equal = [](const finite_volume::ImmersedResidualParts &a,
                              const finite_volume::ImmersedResidualParts &b) {
    for (std::size_t component = 0U; component < 3U; ++component)
      if (bits(a.pressure[component]) != bits(b.pressure[component]) ||
          bits(a.viscous[component]) != bits(b.viscous[component]) ||
          bits(a.convective[component]) != bits(b.convective[component]))
        return false;
    return true;
  };
  return parts_equal(lhs.budget_reaction_N, rhs.budget_reaction_N);
}

bool operator_rows_bits_equal(
    const std::vector<finite_volume::test::ImmersedOperatorRowSnapshot> &left,
    const std::vector<finite_volume::test::ImmersedOperatorRowSnapshot>
        &right) {
  if (left.size() != right.size())
    return false;
  for (std::size_t row = 0U; row < left.size(); ++row) {
    const auto &a = left[row];
    const auto &b = right[row];
    if (a.active_cell != b.active_cell ||
        a.row_replacement_fingerprint != b.row_replacement_fingerprint ||
        a.replacement_group_count != b.replacement_group_count ||
        a.links.size() != b.links.size() ||
        a.covered_physical_terms.size() != b.covered_physical_terms.size())
      return false;
    for (std::size_t link = 0U; link < a.links.size(); ++link) {
      const auto &x = a.links[link];
      const auto &y = b.links[link];
      if (x.id != y.id || x.occurrence != y.occurrence ||
          bits(x.normal_scale) != bits(y.normal_scale) ||
          bits(x.solid_to_fluid_normal.x) !=
              bits(y.solid_to_fluid_normal.x) ||
          bits(x.solid_to_fluid_normal.y) !=
              bits(y.solid_to_fluid_normal.y) ||
          bits(x.solid_to_fluid_normal.z) !=
              bits(y.solid_to_fluid_normal.z) ||
          x.replacement_fingerprint != y.replacement_fingerprint ||
          bits(x.wall_intercept_m.x) != bits(y.wall_intercept_m.x) ||
          bits(x.wall_intercept_m.y) != bits(y.wall_intercept_m.y) ||
          bits(x.wall_intercept_m.z) != bits(y.wall_intercept_m.z) ||
          bits(x.area_from_fluid_m2.x) != bits(y.area_from_fluid_m2.x) ||
          bits(x.area_from_fluid_m2.y) != bits(y.area_from_fluid_m2.y) ||
          bits(x.area_from_fluid_m2.z) != bits(y.area_from_fluid_m2.z) ||
          bits(x.signed_wall_measure_m2) !=
              bits(y.signed_wall_measure_m2) ||
          bits(x.pressure_quadrature_m.x) !=
              bits(y.pressure_quadrature_m.x) ||
          bits(x.pressure_quadrature_m.y) !=
              bits(y.pressure_quadrature_m.y) ||
          bits(x.pressure_quadrature_m.z) !=
              bits(y.pressure_quadrature_m.z))
        return false;
    }
    for (std::size_t term = 0U; term < a.covered_physical_terms.size();
         ++term) {
      const auto &x = a.covered_physical_terms[term];
      const auto &y = b.covered_physical_terms[term];
      if (x.stable_term_id != y.stable_term_id || x.link != y.link ||
          x.kind != y.kind ||
          x.algebraic_occurrence != y.algebraic_occurrence ||
          x.output_component != y.output_component ||
          bits(x.coefficient) != bits(y.coefficient) ||
          x.evaluation_group_id != y.evaluation_group_id ||
          x.source_term_ids != y.source_term_ids)
        return false;
    }
  }
  return true;
}

void test_report_comparison_oracle() {
  finite_volume::ImmersedOperatorReport reference{};
  reference.active_row_count = 4U;
  reference.replacement_group_count = 7U;
  reference.simultaneous_substitution_count = 2U;
  reference.budget_reaction_N.pressure = {1.0, -2.0, 3.0};
  reference.budget_reaction_N.viscous = {-4.0, 5.0, -6.0};
  reference.budget_reaction_N.convective = {7.0, -8.0, 9.0};
  HUNDUN_CHECK(report_bits_equal(reference, reference));
  auto ordinary_mutation = reference;
  ordinary_mutation.active_row_count += 1U;
  HUNDUN_CHECK(!report_bits_equal(reference, ordinary_mutation));
  auto nested_mutation = reference;
  nested_mutation.budget_reaction_N.viscous[1] =
      std::nextafter(nested_mutation.budget_reaction_N.viscous[1],
                     std::numeric_limits<double>::infinity());
  HUNDUN_CHECK(!report_bits_equal(reference, nested_mutation));
}

bool physical_term_coverage_is_complete(
    const finite_volume::test::ImmersedOperatorRowSnapshot &row) {
  using finite_volume::test::ImmersedPhysicalTermKind;
  if (row.links.empty())
    return row.covered_physical_terms.empty();
  if (!std::is_sorted(row.covered_physical_terms.begin(),
                      row.covered_physical_terms.end(),
                      [](const auto &left, const auto &right) {
                        return left.stable_term_id < right.stable_term_id;
                      }))
    return false;
  for (std::size_t index = 1U; index < row.covered_physical_terms.size();
       ++index)
    if (row.covered_physical_terms[index - 1U].stable_term_id ==
        row.covered_physical_terms[index].stable_term_id)
      return false;
  std::uint64_t joint_group = 0U;
  for (const auto &link : row.links) {
    std::array<std::size_t, 4> kind_counts{};
    std::uint64_t link_group = 0U;
    for (const auto &term : row.covered_physical_terms) {
      if (term.link != link.id)
        continue;
      if (term.algebraic_occurrence != link.occurrence ||
          term.output_component >= 3U || !std::isfinite(term.coefficient) ||
          term.evaluation_group_id == 0U ||
          !std::is_sorted(term.source_term_ids.begin(),
                          term.source_term_ids.end()) ||
          std::adjacent_find(term.source_term_ids.begin(),
                             term.source_term_ids.end()) !=
              term.source_term_ids.end() ||
          (term.coefficient != 0.0 && term.source_term_ids.empty()))
        return false;
      if (link_group == 0U)
        link_group = term.evaluation_group_id;
      else if (link_group != term.evaluation_group_id)
        return false;
      const auto slot =
          term.kind == ImmersedPhysicalTermKind::convective_direct   ? 0U
          : term.kind == ImmersedPhysicalTermKind::pressure_direct  ? 1U
          : term.kind == ImmersedPhysicalTermKind::viscous_orthogonal
              ? 2U
          : term.kind ==
                    ImmersedPhysicalTermKind::viscous_deferred_gradient
              ? 3U
              : kind_counts.size();
      if (slot >= kind_counts.size())
        return false;
      ++kind_counts[slot];
    }
    if (link_group == 0U ||
        kind_counts != std::array<std::size_t, 4>{3U, 3U, 3U, 3U})
      return false;
    if (joint_group == 0U)
      joint_group = link_group;
    else if (joint_group != link_group)
      return false;
  }
  std::vector<std::pair<std::uint64_t, immersed::ImmersedLinkId>> groups;
  for (const auto &term : row.covered_physical_terms)
    groups.push_back({term.evaluation_group_id, term.link});
  std::sort(groups.begin(), groups.end());
  groups.erase(std::unique(groups.begin(), groups.end()), groups.end());
  if (groups.size() != row.links.size() || row.replacement_group_count != 1U)
    return false;
  for (std::size_t index = 1U; index < groups.size(); ++index)
    if (groups[index - 1U].first != groups[index].first ||
        groups[index - 1U].second == groups[index].second)
      return false;
  return true;
}

void test_physical_background_term_oracle() {
  using Kind = finite_volume::test::ImmersedPhysicalTermKind;
  finite_volume::test::ImmersedOperatorRowSnapshot reference{};
  reference.active_cell = 17U;
  reference.links.resize(2U);
  reference.links[0].id = 3U;
  reference.links[0].occurrence = 2U;
  reference.links[1].id = 9U;
  reference.links[1].occurrence = 4U;
  reference.replacement_group_count = 1U;
  std::uint64_t stable_id = 1U;
  for (std::size_t link_index = 0U; link_index < reference.links.size();
       ++link_index) {
    const auto &link = reference.links[link_index];
    const std::uint64_t group = 71U;
    for (std::uint32_t component = 0U; component < 3U; ++component)
      for (const auto kind :
           {Kind::convective_direct, Kind::pressure_direct,
            Kind::viscous_orthogonal, Kind::viscous_deferred_gradient}) {
        const double coefficient =
            kind == Kind::convective_direct
                ? 0.0
                : 0.125 * static_cast<double>(stable_id + 1U);
        std::vector<std::uint64_t> sources;
        if (coefficient != 0.0)
          sources.push_back(1000U + stable_id);
        reference.covered_physical_terms.push_back(
            {stable_id++, link.id, kind,
             static_cast<std::uint32_t>(link.occurrence), component,
             coefficient, group, std::move(sources)});
      }
  }
  std::sort(reference.covered_physical_terms.begin(),
            reference.covered_physical_terms.end(),
            [](const auto &left, const auto &right) {
              return left.stable_term_id < right.stable_term_id;
            });
  HUNDUN_CHECK(physical_term_coverage_is_complete(reference));
  auto wrong_component = reference;
  wrong_component.covered_physical_terms.front().output_component = 3U;
  HUNDUN_CHECK(!physical_term_coverage_is_complete(wrong_component));
  auto wrong_kind = reference;
  wrong_kind.covered_physical_terms.front().kind = Kind::neighbour;
  HUNDUN_CHECK(!physical_term_coverage_is_complete(wrong_kind));
  auto wrong_group = reference;
  wrong_group.covered_physical_terms.front().evaluation_group_id = 999U;
  HUNDUN_CHECK(!physical_term_coverage_is_complete(wrong_group));
  auto missing_source = reference;
  const auto nonzero = std::find_if(
      missing_source.covered_physical_terms.begin(),
      missing_source.covered_physical_terms.end(),
      [](const auto &term) { return term.coefficient != 0.0; });
  HUNDUN_CHECK(nonzero != missing_source.covered_physical_terms.end());
  nonzero->source_term_ids.clear();
  HUNDUN_CHECK(!physical_term_coverage_is_complete(missing_source));
}

void prove_physical_term_coverage_oracle_is_mutation_sensitive(
    const finite_volume::test::ImmersedOperatorRowSnapshot &reference) {
  HUNDUN_CHECK(physical_term_coverage_is_complete(reference));
  auto omitted = reference;
  omitted.covered_physical_terms.pop_back();
  HUNDUN_CHECK(!physical_term_coverage_is_complete(omitted));
  auto duplicated = reference;
  duplicated.covered_physical_terms.push_back(
      duplicated.covered_physical_terms.front());
  std::sort(duplicated.covered_physical_terms.begin(),
            duplicated.covered_physical_terms.end(),
            [](const auto &left, const auto &right) {
              return left.stable_term_id < right.stable_term_id;
            });
  HUNDUN_CHECK(!physical_term_coverage_is_complete(duplicated));
  auto wrong_link = reference;
  wrong_link.covered_physical_terms.front().link =
      std::numeric_limits<immersed::ImmersedLinkId>::max();
  HUNDUN_CHECK(!physical_term_coverage_is_complete(wrong_link));
}

bool replacement_term_partition_is_complete(
    const finite_volume::test::BoundaryRowEvaluationSnapshot &evaluation,
    const finite_volume::test::ImmersedOperatorRowSnapshot &row) {
  using Kind = finite_volume::test::BoundaryReplacementTermKind;
  if (row.links.empty())
    return evaluation.replacement_terms.empty();
  if (evaluation.replacement_terms.empty() ||
      evaluation.evaluated_group_count == 0U ||
      !std::is_sorted(
          evaluation.replacement_terms.begin(),
          evaluation.replacement_terms.end(), [](const auto &left,
                                                  const auto &right) {
            return left.stable_term_id < right.stable_term_id;
          }))
    return false;
  std::vector<std::uint64_t> groups;
  std::vector<std::pair<std::uint64_t, immersed::ImmersedLinkId>>
      group_links;
  for (std::size_t index = 0U; index < evaluation.replacement_terms.size();
       ++index) {
    const auto &term = evaluation.replacement_terms[index];
    if (term.stable_term_id == 0U || term.evaluation_group_id == 0U ||
        term.component >= 3U || !std::isfinite(term.value) ||
        (index > 0U &&
         evaluation.replacement_terms[index - 1U].stable_term_id ==
             term.stable_term_id))
      return false;
    const auto found_link = std::lower_bound(
        row.links.begin(), row.links.end(), term.link,
        [](const auto &candidate, immersed::ImmersedLinkId target) {
          return candidate.id < target;
        });
    if (found_link == row.links.end() || found_link->id != term.link)
      return false;
    switch (term.kind) {
    case Kind::pressure_face:
    case Kind::viscous_wall:
      if (term.occurrence != 7U)
        return false;
      break;
    case Kind::pressure_diagonal_defect:
    case Kind::viscous_diagonal_defect:
      if (term.occurrence != 6U)
        return false;
      break;
    case Kind::pressure_neighbour_defect:
    case Kind::viscous_neighbour_defect:
      if (term.occurrence >= 6U)
        return false;
      break;
    }
    groups.push_back(term.evaluation_group_id);
    group_links.push_back({term.evaluation_group_id, term.link});
  }
  std::sort(groups.begin(), groups.end());
  groups.erase(std::unique(groups.begin(), groups.end()), groups.end());
  std::sort(group_links.begin(), group_links.end());
  group_links.erase(std::unique(group_links.begin(), group_links.end()),
                    group_links.end());
  if (groups.size() != evaluation.evaluated_group_count ||
      groups.size() != 1U ||
      group_links.size() != row.links.size())
    return false;
  for (std::size_t index = 1U; index < group_links.size(); ++index)
    if (group_links[index - 1U].first != group_links[index].first ||
        group_links[index - 1U].second == group_links[index].second)
      return false;
  for (const auto &link : row.links) {
    for (std::uint32_t component = 0U; component < 3U; ++component) {
      std::array<std::size_t, 6> counts{};
      std::array<std::size_t, 6> pressure_neighbours{};
      std::array<std::size_t, 6> viscous_neighbours{};
      for (const auto &term : evaluation.replacement_terms) {
        if (term.link != link.id || term.component != component)
          continue;
        switch (term.kind) {
        case Kind::pressure_face:
          ++counts[0];
          break;
        case Kind::pressure_diagonal_defect:
          ++counts[1];
          break;
        case Kind::pressure_neighbour_defect:
          ++pressure_neighbours[term.occurrence];
          break;
        case Kind::viscous_wall:
          ++counts[3];
          break;
        case Kind::viscous_diagonal_defect:
          ++counts[4];
          break;
        case Kind::viscous_neighbour_defect:
          ++viscous_neighbours[term.occurrence];
          break;
        }
      }
      counts[2] = std::accumulate(pressure_neighbours.begin(),
                                  pressure_neighbours.end(), std::size_t{0});
      counts[5] = std::accumulate(viscous_neighbours.begin(),
                                  viscous_neighbours.end(), std::size_t{0});
      if (counts != std::array<std::size_t, 6>{1U, 1U, 6U, 1U, 1U, 6U} ||
          std::any_of(pressure_neighbours.begin(), pressure_neighbours.end(),
                      [](std::size_t value) { return value != 1U; }) ||
          std::any_of(viscous_neighbours.begin(), viscous_neighbours.end(),
                      [](std::size_t value) { return value != 1U; }))
        return false;
    }
  }
  return true;
}

void test_replacement_term_partition_oracle() {
  using Kind = finite_volume::test::BoundaryReplacementTermKind;
  finite_volume::test::ImmersedOperatorRowSnapshot row{};
  row.active_cell = 17U;
  row.links.resize(2U);
  row.links[0].id = 3U;
  row.links[0].occurrence = 2U;
  row.links[1].id = 9U;
  row.links[1].occurrence = 4U;
  finite_volume::test::BoundaryRowEvaluationSnapshot reference{};
  reference.active_cell = row.active_cell;
  reference.evaluated_group_count = 1U;
  std::uint64_t stable_id = 1U;
  const auto add = [&](immersed::ImmersedLinkId link, Kind kind,
                       std::uint32_t occurrence, std::uint32_t component,
                       std::uint64_t group) {
    reference.replacement_terms.push_back(
        {stable_id++, link, kind, occurrence, component,
         0.125 * static_cast<double>(stable_id), group});
  };
  for (std::size_t link_index = 0U; link_index < row.links.size(); ++link_index) {
    const auto &link = row.links[link_index];
    const auto group = 73U;
    for (std::uint32_t component = 0U; component < 3U; ++component) {
      add(link.id, Kind::pressure_face, 7U, component, group);
      add(link.id, Kind::pressure_diagonal_defect, 6U, component, group);
      for (std::uint32_t occurrence = 0U; occurrence < 6U; ++occurrence)
        add(link.id, Kind::pressure_neighbour_defect, occurrence, component,
            group);
      add(link.id, Kind::viscous_wall, 7U, component, group);
      add(link.id, Kind::viscous_diagonal_defect, 6U, component, group);
      for (std::uint32_t occurrence = 0U; occurrence < 6U; ++occurrence)
        add(link.id, Kind::viscous_neighbour_defect, occurrence, component,
            group);
    }
  }
  std::sort(reference.replacement_terms.begin(),
            reference.replacement_terms.end(), [](const auto &left,
                                                   const auto &right) {
              return left.stable_term_id < right.stable_term_id;
            });
  HUNDUN_CHECK(replacement_term_partition_is_complete(reference, row));
  auto omitted = reference;
  omitted.replacement_terms.pop_back();
  HUNDUN_CHECK(!replacement_term_partition_is_complete(omitted, row));
  auto duplicated = reference;
  duplicated.replacement_terms.push_back(duplicated.replacement_terms.back());
  HUNDUN_CHECK(!replacement_term_partition_is_complete(duplicated, row));
  auto wrong_link = reference;
  wrong_link.replacement_terms.front().link = 99U;
  HUNDUN_CHECK(!replacement_term_partition_is_complete(wrong_link, row));
  auto wrong_occurrence = reference;
  wrong_occurrence.replacement_terms.front().occurrence = 6U;
  HUNDUN_CHECK(
      !replacement_term_partition_is_complete(wrong_occurrence, row));
  auto missing_group = reference;
  missing_group.replacement_terms.front().evaluation_group_id = 0U;
  HUNDUN_CHECK(!replacement_term_partition_is_complete(missing_group, row));
}

bool affine_source_trace_is_complete(
    const finite_volume::test::BoundaryRowEvaluationSnapshot &evaluation) {
  std::vector<std::uint64_t> declared;
  std::vector<std::uint64_t> required;
  for (const auto &term : evaluation.replacement_terms) {
    declared.push_back(term.stable_term_id);
    if (term.value != 0.0)
      required.push_back(term.stable_term_id);
  }
  std::sort(declared.begin(), declared.end());
  std::sort(required.begin(), required.end());
  if (std::adjacent_find(declared.begin(), declared.end()) != declared.end())
    return false;
  std::vector<std::uint64_t> observed;
  const auto append_sources = [&](const auto &term) {
    if (!std::is_sorted(term.source_term_ids.begin(),
                        term.source_term_ids.end()) ||
        std::adjacent_find(term.source_term_ids.begin(),
                           term.source_term_ids.end()) !=
            term.source_term_ids.end())
      return false;
    observed.insert(observed.end(), term.source_term_ids.begin(),
                    term.source_term_ids.end());
    return true;
  };
  for (const auto &term : evaluation.affine_donor_terms)
    if (!append_sources(term))
      return false;
  for (const auto &term : evaluation.affine_wall_gradient_terms)
    if (!append_sources(term))
      return false;
  std::sort(observed.begin(), observed.end());
  observed.erase(std::unique(observed.begin(), observed.end()), observed.end());
  return std::includes(declared.begin(), declared.end(), observed.begin(),
                       observed.end()) &&
         std::includes(observed.begin(), observed.end(), required.begin(),
                       required.end());
}

bool background_affine_source_trace_is_complete(
    const finite_volume::test::BoundaryRowEvaluationSnapshot &evaluation,
    const finite_volume::test::ImmersedOperatorRowSnapshot &row) {
  std::vector<std::uint64_t> declared;
  std::vector<std::uint64_t> required;
  for (const auto &term : row.covered_physical_terms) {
    declared.push_back(term.stable_term_id);
    if (term.coefficient != 0.0)
      required.push_back(term.stable_term_id);
    const std::vector<std::uint64_t> expected_sources =
        term.coefficient == 0.0
            ? std::vector<std::uint64_t>{}
            : std::vector<std::uint64_t>{term.stable_term_id};
    if (term.source_term_ids != expected_sources)
      return false;
  }
  std::sort(declared.begin(), declared.end());
  std::sort(required.begin(), required.end());
  if (std::adjacent_find(declared.begin(), declared.end()) != declared.end())
    return false;
  std::vector<std::uint64_t> observed;
  const auto append_sources = [&](const auto &term) {
    if (!std::is_sorted(term.source_term_ids.begin(),
                        term.source_term_ids.end()) ||
        std::adjacent_find(term.source_term_ids.begin(),
                           term.source_term_ids.end()) !=
            term.source_term_ids.end())
      return false;
    observed.insert(observed.end(), term.source_term_ids.begin(),
                    term.source_term_ids.end());
    return true;
  };
  for (const auto &term : evaluation.background_affine_donor_terms)
    if (!append_sources(term))
      return false;
  for (const auto &term : evaluation.background_affine_wall_gradient_terms)
    if (!append_sources(term))
      return false;
  std::sort(observed.begin(), observed.end());
  observed.erase(std::unique(observed.begin(), observed.end()), observed.end());
  return std::includes(declared.begin(), declared.end(), observed.begin(),
                       observed.end()) &&
         std::includes(observed.begin(), observed.end(), required.begin(),
                       required.end());
}

void test_row_plan_snapshot(
    const finite_volume::ImmersedOperatorAdapter &adapter,
    const immersed::LocalFlowPatternTransform &transform) {
  const auto rows =
      finite_volume::test::ImmersedOperatorTestAccess::rows(adapter);
  HUNDUN_CHECK(std::is_sorted(rows.begin(), rows.end(),
                              [](const auto &lhs, const auto &rhs) {
                                return lhs.active_cell < rhs.active_cell;
                              }));
  std::size_t multi_link_rows = 0U;
  std::size_t transformed_quadrature_points = 0U;
  const finite_volume::test::ImmersedOperatorRowSnapshot *coverage_sample =
      nullptr;
  for (const auto &row : rows) {
    HUNDUN_CHECK(row.links.empty() == (row.row_replacement_fingerprint == 0U));
    if (!row.links.empty()) {
      HUNDUN_CHECK(row.row_replacement_fingerprint != 0U);
      HUNDUN_CHECK(row.replacement_group_count > 0U);
    }
    HUNDUN_CHECK(std::is_sorted(
        row.links.begin(), row.links.end(),
        [](const auto &lhs, const auto &rhs) { return lhs.id < rhs.id; }));
    std::array<bool, 6> occurrence_seen{};
    for (const auto &snapshot : row.links) {
      HUNDUN_CHECK(snapshot.occurrence < occurrence_seen.size());
      HUNDUN_CHECK(!occurrence_seen[snapshot.occurrence]);
      occurrence_seen[snapshot.occurrence] = true;
      HUNDUN_CHECK(snapshot.normal_scale > 0.0 && snapshot.normal_scale < 1.0 &&
                   std::isfinite(snapshot.normal_scale));
      HUNDUN_CHECK(std::isfinite(snapshot.solid_to_fluid_normal.x) &&
                   std::isfinite(snapshot.solid_to_fluid_normal.y) &&
                   std::isfinite(snapshot.solid_to_fluid_normal.z));
      HUNDUN_CHECK(std::isfinite(snapshot.pressure_quadrature_m.x));
      HUNDUN_CHECK(std::isfinite(snapshot.pressure_quadrature_m.y));
      HUNDUN_CHECK(std::isfinite(snapshot.pressure_quadrature_m.z));
      const double independent_signed_measure =
          -dot(snapshot.area_from_fluid_m2, snapshot.solid_to_fluid_normal);
      HUNDUN_CHECK(snapshot.signed_wall_measure_m2 > 0.0);
      HUNDUN_CHECK(bits(snapshot.signed_wall_measure_m2) ==
                   bits(independent_signed_measure));
      const auto displacement =
          subtract(snapshot.pressure_quadrature_m, snapshot.wall_intercept_m);
      if (dot(displacement, displacement) >
          256.0 * std::numeric_limits<double>::epsilon())
        ++transformed_quadrature_points;
      HUNDUN_CHECK(snapshot.replacement_fingerprint ==
                   row.row_replacement_fingerprint);
    }
    if (row.links.size() > 1U)
      ++multi_link_rows;
    if (coverage_sample == nullptr && !row.links.empty())
      coverage_sample = &row;
    HUNDUN_CHECK(physical_term_coverage_is_complete(row));
  }
  HUNDUN_CHECK(multi_link_rows > 0U);
  HUNDUN_CHECK(transformed_quadrature_points > 0U);
  HUNDUN_CHECK(coverage_sample != nullptr);
  prove_physical_term_coverage_oracle_is_mutation_sensitive(*coverage_sample);
  static_cast<void>(transform);
}

void prove_wall_pressure_authority_matches_operator_row(
    const finite_volume::ImmersedOperatorAdapter &adapter,
    const immersed::WallQuadraturePlan &wall_plan,
    const immersed::GhostStencilPlan &ghost_plan,
    const runtime::MpiContext &mpi) {
  struct AuthorityProbe final {
    std::uint64_t link{};
    runtime::Real3 point_m{};
    double value{};
    runtime::Real3 gradient{};
    double condition_estimate{};
    std::uint64_t donor_fingerprint{};
    std::uint64_t donor_count{};
  };
  static_assert(std::is_trivially_copyable_v<AuthorityProbe>);

  const auto donor_value = [](mesh::GlobalCellId donor) {
    const double gid = static_cast<double>(donor);
    return std::sin(0.013 * (gid + 1.0)) +
           0.375 * std::cos(0.029 * (gid + 7.0)) +
           0.1 * std::sin(0.0007 * (gid + 3.0) * (gid + 3.0));
  };
  const auto weighted_value = [&](const auto &weights) {
    double result = 0.0;
    for (const auto &donor : weights)
      result += donor.weight * donor_value(donor.global_cell);
    return result;
  };
  const auto evaluate = [&](const immersed::QuadraticReconstruction &authority,
                            runtime::Real3 point) {
    std::array<double, 4> result{};
    result[0] = weighted_value(
        immersed::detail::QuadraticReconstructionWeights::value_weights(
            authority, point));
    constexpr std::array<runtime::Real3, 3> directions{
        runtime::Real3{1.0, 0.0, 0.0}, runtime::Real3{0.0, 1.0, 0.0},
        runtime::Real3{0.0, 0.0, 1.0}};
    for (std::size_t direction = 0U; direction < directions.size(); ++direction)
      result[direction + 1U] = weighted_value(
          immersed::detail::QuadraticReconstructionWeights::
              directional_gradient_weights(authority, point,
                                           directions[direction]));
    return result;
  };
  const auto donor_fingerprint = [](const auto &authority) {
    std::uint64_t result = UINT64_C(1469598103934665603);
    for (const auto donor :
         immersed::detail::QuadraticReconstructionWeights::donor_global_ids(
             authority)) {
      result ^= donor;
      result *= UINT64_C(1099511628211);
    }
    return result;
  };

  std::vector<AuthorityProbe> local_probes;
  local_probes.reserve(wall_plan.local_points().size());
  for (const auto &point : wall_plan.local_points()) {
    const auto authority =
        immersed::detail::boundary_authority_reconstruction(
            point.reconstruction);
    const auto sample = evaluate(authority, point.position_m);
    local_probes.push_back(
        {immersed::detail::boundary_authority_link(point.reconstruction),
         point.position_m,
         sample[0],
         {sample[1], sample[2], sample[3]},
         authority.quality().condition_estimate,
         donor_fingerprint(authority),
         static_cast<std::uint64_t>(
             immersed::detail::QuadraticReconstructionWeights::
                 donor_global_ids(authority)
                     .size())});
  }
  HUNDUN_CHECK(local_probes.size() <=
               static_cast<std::size_t>(std::numeric_limits<int>::max()) /
                   sizeof(AuthorityProbe));
  const int local_bytes =
      static_cast<int>(local_probes.size() * sizeof(AuthorityProbe));
  std::vector<int> byte_counts(static_cast<std::size_t>(mpi.size()));
  HUNDUN_CHECK(MPI_Allgather(&local_bytes, 1, MPI_INT, byte_counts.data(), 1,
                             MPI_INT, mpi.comm()) == MPI_SUCCESS);
  std::vector<int> displacements(byte_counts.size());
  int total_bytes = 0;
  for (std::size_t rank = 0U; rank < byte_counts.size(); ++rank) {
    HUNDUN_CHECK(byte_counts[rank] >= 0 &&
                 byte_counts[rank] % static_cast<int>(sizeof(AuthorityProbe)) ==
                     0);
    HUNDUN_CHECK(total_bytes <=
                 std::numeric_limits<int>::max() - byte_counts[rank]);
    displacements[rank] = total_bytes;
    total_bytes += byte_counts[rank];
  }
  HUNDUN_CHECK(total_bytes % static_cast<int>(sizeof(AuthorityProbe)) == 0);
  std::vector<AuthorityProbe> probes(
      static_cast<std::size_t>(total_bytes) / sizeof(AuthorityProbe));
  HUNDUN_CHECK(MPI_Allgatherv(local_probes.data(), local_bytes, MPI_BYTE,
                              probes.data(), byte_counts.data(),
                              displacements.data(), MPI_BYTE, mpi.comm()) ==
               MPI_SUCCESS);

  std::vector<immersed::ImmersedLinkId> multi_link_row_links;
  for (const auto &row :
       finite_volume::test::ImmersedOperatorTestAccess::rows(adapter)) {
    if (row.links.size() <= 1U)
      continue;
    for (const auto &link : row.links)
      multi_link_row_links.push_back(link.id);
  }
  std::sort(multi_link_row_links.begin(), multi_link_row_links.end());

  std::uint64_t local_comparisons = 0U;
  int local_legacy_mutation_detected = 0;
  for (const auto &probe : probes) {
    if (!std::binary_search(multi_link_row_links.begin(),
                            multi_link_row_links.end(), probe.link))
      continue;

    const auto &row_authority =
        finite_volume::detail::ImmersedBoundaryAuthorityAccess::
            row_reconstruction(adapter, probe.link);
    const auto row_sample = evaluate(row_authority, probe.point_m);
    const double row_value = row_sample[0];
    const runtime::Real3 row_gradient{row_sample[1], row_sample[2],
                                      row_sample[3]};
    HUNDUN_CHECK(donor_fingerprint(row_authority) ==
                 probe.donor_fingerprint);
    HUNDUN_CHECK(
        immersed::detail::QuadraticReconstructionWeights::donor_global_ids(
            row_authority)
            .size() == probe.donor_count);
    const double scale =
        std::max({1.0, std::abs(row_value), std::abs(probe.value),
                  max_norm(row_gradient), max_norm(probe.gradient)});
    const double tolerance =
        8192.0 * std::numeric_limits<double>::epsilon() * scale *
        (1.0 + std::max(row_authority.quality().condition_estimate,
                        probe.condition_estimate));
    HUNDUN_CHECK(std::abs(row_value - probe.value) <= tolerance);
    HUNDUN_CHECK(max_norm(subtract(row_gradient, probe.gradient)) <=
                 tolerance);

    const auto &legacy = ghost_plan.reconstruction(probe.link);
    const auto legacy_sample = evaluate(legacy, probe.point_m);
    const double legacy_value = legacy_sample[0];
    const runtime::Real3 legacy_gradient{legacy_sample[1], legacy_sample[2],
                                         legacy_sample[3]};
    if (std::abs(legacy_value - row_value) > tolerance ||
        max_norm(subtract(legacy_gradient, row_gradient)) > tolerance)
      local_legacy_mutation_detected = 1;
    ++local_comparisons;
  }

  std::uint64_t global_comparisons = local_comparisons;
  int global_legacy_mutation_detected = local_legacy_mutation_detected;
  HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, &global_comparisons, 1,
                             MPI_UINT64_T, MPI_SUM, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE,
                             &global_legacy_mutation_detected, 1, MPI_INT,
                             MPI_MAX, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(global_comparisons > 0U);
  HUNDUN_CHECK(global_legacy_mutation_detected == 1);
}

bool boundary_evaluation_is_complete(
    const finite_volume::test::BoundaryRowEvaluationSnapshot &evaluation,
    const finite_volume::test::ImmersedOperatorRowSnapshot &row) {
  using InputKind = finite_volume::test::BoundaryAffineInputKind;
  if (evaluation.active_cell != row.active_cell ||
      evaluation.row_fingerprint != row.row_replacement_fingerprint ||
      evaluation.evaluated_group_count == 0U ||
      evaluation.simultaneous_substitution_count != 1U ||
      evaluation.affine_plan_fingerprint == 0U ||
      evaluation.affine_donor_terms.empty() ||
      evaluation.canonical_affine_row_evaluation_count != 1U ||
      evaluation.link_local_runtime_evaluation_count != 0U ||
      evaluation.immutable_input_snapshot_count != 1U ||
      evaluation.background_functional_evaluation_count != 1U ||
      evaluation.background_removal_count != 1U ||
      evaluation.final_row_write_count != 1U ||
      !replacement_term_partition_is_complete(evaluation, row) ||
      !background_affine_source_trace_is_complete(evaluation, row) ||
      !affine_source_trace_is_complete(evaluation))
    return false;
  for (std::size_t index = 0U; index < evaluation.affine_donor_terms.size();
       ++index) {
    const auto &term = evaluation.affine_donor_terms[index];
    if (term.stable_term_id == 0U || term.output_component >= 3U ||
        term.input_component >=
            (term.input_kind == InputKind::pressure ? 1U : 3U) ||
        !std::isfinite(term.coefficient) ||
        term.contributing_link_count == 0U ||
        term.contributing_link_count > row.links.size())
      return false;
    if (index > 0U) {
      const auto &previous = evaluation.affine_donor_terms[index - 1U];
      if (std::tie(previous.input_kind, previous.donor_global_cell,
                   previous.input_component, previous.output_component) >=
          std::tie(term.input_kind, term.donor_global_cell,
                   term.input_component, term.output_component))
        return false;
    }
  }
  for (std::size_t index = 0U;
       index < evaluation.affine_wall_gradient_terms.size(); ++index) {
    const auto &term = evaluation.affine_wall_gradient_terms[index];
    const auto found_link = std::lower_bound(
        row.links.begin(), row.links.end(), term.link,
        [](const auto &candidate, immersed::ImmersedLinkId target) {
          return candidate.id < target;
        });
    if (term.stable_term_id == 0U || term.output_component >= 3U ||
        !std::isfinite(term.coefficient) ||
        found_link == row.links.end() || found_link->id != term.link)
      return false;
    if (index > 0U) {
      const auto &previous =
          evaluation.affine_wall_gradient_terms[index - 1U];
      if (std::tie(previous.link, previous.output_component) >=
          std::tie(term.link, term.output_component))
        return false;
    }
  }
  for (std::size_t component = 0U; component < 3U; ++component) {
    if (bits(evaluation.background_contribution.convective[component]) !=
            bits(evaluation.removed_background_contribution
                     .convective[component]) ||
        bits(evaluation.background_contribution.pressure[component]) !=
            bits(evaluation.removed_background_contribution
                     .pressure[component]) ||
        bits(evaluation.background_contribution.viscous[component]) !=
            bits(evaluation.removed_background_contribution
                     .viscous[component]))
      return false;
    const double background =
        evaluation.background_contribution.convective[component] +
        evaluation.background_contribution.pressure[component] +
        evaluation.background_contribution.viscous[component];
    const double removed =
        evaluation.removed_background_contribution.convective[component] +
        evaluation.removed_background_contribution.pressure[component] +
        evaluation.removed_background_contribution.viscous[component];
    const double expected =
        (evaluation.residual_before_wall[component] + background) - removed +
        evaluation.wall_contribution.pressure[component] +
        evaluation.wall_contribution.viscous[component];
    if (bits(evaluation.wall_contribution.convective[component]) != bits(0.0) ||
        bits(evaluation.residual_after_wall[component]) != bits(expected))
      return false;
  }
  return true;
}

bool affine_terms_match_quadratic_oracle(
    const std::vector<finite_volume::test::BoundaryAffineDonorTermSnapshot>
        &donor_terms,
    const std::vector<
        finite_volume::test::BoundaryAffineWallGradientTermSnapshot>
        &wall_gradient_terms,
    const finite_volume::test::BoundaryResidualPartsSnapshot &contribution,
    const mesh::MeshGeometry &geometry,
    const std::vector<finite_volume::detail::ImmersedWallNormalGradient>
        *wall_normal_gradients = nullptr) {
  using InputKind = finite_volume::test::BoundaryAffineInputKind;
  constexpr std::array<double, 3> velocity_scale{0.3, -0.2, 0.1};
  std::array<double, 3> pressure{};
  std::array<double, 3> viscous{};
  for (const auto &term : donor_terms) {
    const double average =
        quadratic_average(logical_cell(term.donor_global_cell, kExtent),
                          geometry);
    auto &target = term.input_kind == InputKind::pressure ? pressure : viscous;
    const double value =
        term.input_kind == InputKind::pressure
            ? average
            : velocity_scale[term.input_component] * average;
    target[term.output_component] += term.coefficient * value;
  }
  for (const auto &term : wall_gradient_terms) {
    if (wall_normal_gradients == nullptr)
      return false;
    const auto condition = std::lower_bound(
        wall_normal_gradients->begin(), wall_normal_gradients->end(), term.link,
        [](const auto &candidate, immersed::ImmersedLinkId target) {
          return candidate.link < target;
        });
    if (condition == wall_normal_gradients->end() ||
        condition->link != term.link)
      return false;
    pressure[term.output_component] += term.coefficient * condition->value;
  }
  for (std::size_t component = 0U; component < 3U; ++component) {
    const double pressure_scale =
        std::max({1.0, std::abs(pressure[component]),
                  std::abs(contribution.pressure[component])});
    const double viscous_scale =
        std::max({1.0, std::abs(viscous[component]),
                  std::abs(contribution.viscous[component])});
    if (bits(contribution.convective[component]) != bits(0.0) ||
        std::abs(pressure[component] - contribution.pressure[component]) >
            64.0 * std::numeric_limits<double>::epsilon() * pressure_scale ||
        std::abs(viscous[component] - contribution.viscous[component]) >
            64.0 * std::numeric_limits<double>::epsilon() * viscous_scale)
      return false;
  }
  return true;
}

bool affine_evaluation_matches_quadratic_oracle(
    const finite_volume::test::BoundaryRowEvaluationSnapshot &evaluation,
    const mesh::MeshGeometry &geometry,
    const std::vector<finite_volume::detail::ImmersedWallNormalGradient>
        *wall_normal_gradients = nullptr) {
  return affine_terms_match_quadratic_oracle(
      evaluation.affine_donor_terms, evaluation.affine_wall_gradient_terms,
      evaluation.wall_contribution, geometry, wall_normal_gradients);
}

bool background_affine_evaluation_matches_quadratic_oracle(
    const finite_volume::test::BoundaryRowEvaluationSnapshot &evaluation,
    const mesh::MeshGeometry &geometry,
    const std::vector<finite_volume::detail::ImmersedWallNormalGradient>
        *wall_normal_gradients = nullptr) {
  return affine_terms_match_quadratic_oracle(
      evaluation.background_affine_donor_terms,
      evaluation.background_affine_wall_gradient_terms,
      evaluation.background_contribution, geometry, wall_normal_gradients);
}

void prove_boundary_evaluation_oracle_is_mutation_sensitive(
    const finite_volume::test::BoundaryRowEvaluationSnapshot &reference,
    const finite_volume::test::ImmersedOperatorRowSnapshot &row,
    const mesh::MeshGeometry &geometry) {
  HUNDUN_CHECK(boundary_evaluation_is_complete(reference, row));
  HUNDUN_CHECK(affine_evaluation_matches_quadratic_oracle(reference, geometry));
  HUNDUN_CHECK(background_affine_evaluation_matches_quadratic_oracle(
      reference, geometry));
  auto omitted_affine = reference;
  omitted_affine.affine_donor_terms.pop_back();
  HUNDUN_CHECK(
      !affine_evaluation_matches_quadratic_oracle(omitted_affine, geometry));
  HUNDUN_CHECK(!reference.background_affine_donor_terms.empty());
  auto background_coefficient_mutation = reference;
  background_coefficient_mutation.background_affine_donor_terms.front()
      .coefficient += 0.125 *
                      (1.0 + std::abs(background_coefficient_mutation
                                          .background_affine_donor_terms
                                          .front()
                                          .coefficient));
  HUNDUN_CHECK(!background_affine_evaluation_matches_quadratic_oracle(
      background_coefficient_mutation, geometry));
  auto omitted_background_affine = reference;
  omitted_background_affine.background_affine_donor_terms.pop_back();
  HUNDUN_CHECK(!background_affine_evaluation_matches_quadratic_oracle(
      omitted_background_affine, geometry));
  auto duplicate_affine = reference;
  duplicate_affine.affine_donor_terms.push_back(
      duplicate_affine.affine_donor_terms.front());
  HUNDUN_CHECK(!boundary_evaluation_is_complete(duplicate_affine, row));
  auto coefficient_mutation = reference;
  coefficient_mutation.affine_donor_terms.front().coefficient +=
      0.125 *
      (1.0 + std::abs(
                 coefficient_mutation.affine_donor_terms.front().coefficient));
  HUNDUN_CHECK(!affine_evaluation_matches_quadratic_oracle(
      coefficient_mutation, geometry));
  auto donor_mutation = reference;
  donor_mutation.affine_donor_terms.front().donor_global_cell ^= 1U;
  HUNDUN_CHECK(
      !affine_evaluation_matches_quadratic_oracle(donor_mutation, geometry));
  auto input_mutation = reference;
  input_mutation.affine_donor_terms.front().input_kind =
      input_mutation.affine_donor_terms.front().input_kind ==
              finite_volume::test::BoundaryAffineInputKind::pressure
          ? finite_volume::test::BoundaryAffineInputKind::velocity
          : finite_volume::test::BoundaryAffineInputKind::pressure;
  input_mutation.affine_donor_terms.front().input_component = 0U;
  HUNDUN_CHECK(
      !affine_evaluation_matches_quadratic_oracle(input_mutation, geometry));
  auto output_mutation = reference;
  output_mutation.affine_donor_terms.front().output_component =
      (output_mutation.affine_donor_terms.front().output_component + 1U) % 3U;
  HUNDUN_CHECK(
      !affine_evaluation_matches_quadratic_oracle(output_mutation, geometry));
  auto missing_contributor = reference;
  missing_contributor.affine_donor_terms.front().contributing_link_count = 0U;
  HUNDUN_CHECK(!boundary_evaluation_is_complete(missing_contributor, row));
  auto missing_source = reference;
  const auto required_term = std::find_if(
      missing_source.replacement_terms.begin(),
      missing_source.replacement_terms.end(),
      [](const auto &term) { return term.value != 0.0; });
  HUNDUN_CHECK(required_term != missing_source.replacement_terms.end());
  for (auto &term : missing_source.affine_donor_terms)
    term.source_term_ids.erase(
        std::remove(term.source_term_ids.begin(), term.source_term_ids.end(),
                    required_term->stable_term_id),
        term.source_term_ids.end());
  for (auto &term : missing_source.affine_wall_gradient_terms)
    term.source_term_ids.erase(
        std::remove(term.source_term_ids.begin(), term.source_term_ids.end(),
                    required_term->stable_term_id),
        term.source_term_ids.end());
  HUNDUN_CHECK(!affine_source_trace_is_complete(missing_source));
  auto unknown_source = reference;
  unknown_source.affine_donor_terms.front().source_term_ids.push_back(
      std::numeric_limits<std::uint64_t>::max());
  std::sort(unknown_source.affine_donor_terms.front().source_term_ids.begin(),
            unknown_source.affine_donor_terms.front().source_term_ids.end());
  HUNDUN_CHECK(!affine_source_trace_is_complete(unknown_source));
  auto changed_removal = reference;
  changed_removal.removed_background_contribution.pressure[0] =
      std::nextafter(changed_removal.removed_background_contribution
                         .pressure[0],
                     std::numeric_limits<double>::infinity());
  HUNDUN_CHECK(!boundary_evaluation_is_complete(changed_removal, row));
  auto subsequent_snapshot = reference;
  subsequent_snapshot.immutable_input_snapshot_count = 2U;
  HUNDUN_CHECK(!boundary_evaluation_is_complete(subsequent_snapshot, row));
  auto duplicate_write = reference;
  duplicate_write.final_row_write_count = 2U;
  HUNDUN_CHECK(!boundary_evaluation_is_complete(duplicate_write, row));
}

std::vector<finite_volume::test::BoundaryRowEvaluationSnapshot>
test_boundary_row_evaluations(
    const finite_volume::ImmersedOperatorAdapter &adapter,
    const mesh::MeshGeometry &geometry) {
  const auto rows =
      finite_volume::test::ImmersedOperatorTestAccess::rows(adapter);
  const auto evaluations = finite_volume::test::ImmersedOperatorTestAccess::
      last_boundary_row_evaluations(adapter);
  std::size_t expected = 0U;
  bool observed_multi_link = false;
  bool mutation_proved = false;
  for (const auto &row : rows) {
    if (row.links.empty())
      continue;
    HUNDUN_CHECK(expected < evaluations.size());
    if (!boundary_evaluation_is_complete(evaluations[expected], row)) {
      std::cerr << "row_closure_diag cell=" << row.active_cell
                << " links=" << row.links.size()
                << " replacement_terms="
                << evaluations[expected].replacement_terms.size()
                << " affine_terms="
                << evaluations[expected].affine_donor_terms.size()
                << " wall_terms="
                << evaluations[expected].affine_wall_gradient_terms.size()
                << " groups=" << evaluations[expected].evaluated_group_count
                << " partition="
                << replacement_term_partition_is_complete(
                       evaluations[expected], row)
                << " trace="
                << affine_source_trace_is_complete(evaluations[expected])
                << '\n';
    }
    HUNDUN_CHECK(boundary_evaluation_is_complete(evaluations[expected], row));
    HUNDUN_CHECK(affine_evaluation_matches_quadratic_oracle(
        evaluations[expected], geometry));
    HUNDUN_CHECK(background_affine_evaluation_matches_quadratic_oracle(
        evaluations[expected], geometry));
    observed_multi_link = observed_multi_link || row.links.size() > 1U;
    if (!mutation_proved && row.links.size() > 1U) {
      prove_boundary_evaluation_oracle_is_mutation_sensitive(
          evaluations[expected], row, geometry);
      mutation_proved = true;
    }
    ++expected;
  }
  HUNDUN_CHECK(expected == evaluations.size());
  HUNDUN_CHECK(observed_multi_link);
  HUNDUN_CHECK(mutation_proved);
  return evaluations;
}

bool boundary_evaluations_bits_equal(
    const std::vector<finite_volume::test::BoundaryRowEvaluationSnapshot> &left,
    const std::vector<finite_volume::test::BoundaryRowEvaluationSnapshot>
        &right) {
  if (left.size() != right.size())
    return false;
  for (std::size_t row = 0U; row < left.size(); ++row) {
    if (left[row].active_cell != right[row].active_cell ||
        left[row].row_fingerprint != right[row].row_fingerprint ||
        left[row].affine_plan_fingerprint !=
            right[row].affine_plan_fingerprint ||
        left[row].replacement_terms.size() !=
            right[row].replacement_terms.size() ||
        left[row].background_affine_donor_terms.size() !=
            right[row].background_affine_donor_terms.size() ||
        left[row].background_affine_wall_gradient_terms.size() !=
            right[row].background_affine_wall_gradient_terms.size() ||
        left[row].affine_donor_terms.size() !=
            right[row].affine_donor_terms.size() ||
        left[row].affine_wall_gradient_terms.size() !=
            right[row].affine_wall_gradient_terms.size())
      return false;
    for (std::size_t term = 0U; term < left[row].replacement_terms.size();
         ++term) {
      const auto &a = left[row].replacement_terms[term];
      const auto &b = right[row].replacement_terms[term];
      if (a.stable_term_id != b.stable_term_id || a.link != b.link ||
          a.kind != b.kind || a.occurrence != b.occurrence ||
          a.component != b.component || bits(a.value) != bits(b.value) ||
          a.evaluation_group_id != b.evaluation_group_id)
        return false;
    }
    const auto donor_terms_equal = [](const auto &a, const auto &b) {
      return a.stable_term_id == b.stable_term_id &&
             a.input_kind == b.input_kind &&
             a.donor_global_cell == b.donor_global_cell &&
             a.input_component == b.input_component &&
             a.output_component == b.output_component &&
             bits(a.coefficient) == bits(b.coefficient) &&
             a.contributing_link_count == b.contributing_link_count &&
             a.source_term_ids == b.source_term_ids;
    };
    for (std::size_t term = 0U;
         term < left[row].background_affine_donor_terms.size(); ++term)
      if (!donor_terms_equal(left[row].background_affine_donor_terms[term],
                             right[row].background_affine_donor_terms[term]))
        return false;
    for (std::size_t term = 0U; term < left[row].affine_donor_terms.size();
         ++term) {
      const auto &a = left[row].affine_donor_terms[term];
      const auto &b = right[row].affine_donor_terms[term];
      if (!donor_terms_equal(a, b))
        return false;
    }
    const auto wall_terms_equal = [](const auto &a, const auto &b) {
      return a.stable_term_id == b.stable_term_id && a.link == b.link &&
             a.output_component == b.output_component &&
             bits(a.coefficient) == bits(b.coefficient) &&
             a.source_term_ids == b.source_term_ids;
    };
    for (std::size_t term = 0U;
         term < left[row].background_affine_wall_gradient_terms.size();
         ++term)
      if (!wall_terms_equal(
              left[row].background_affine_wall_gradient_terms[term],
              right[row].background_affine_wall_gradient_terms[term]))
        return false;
    for (std::size_t term = 0U;
         term < left[row].affine_wall_gradient_terms.size(); ++term) {
      const auto &a = left[row].affine_wall_gradient_terms[term];
      const auto &b = right[row].affine_wall_gradient_terms[term];
      if (!wall_terms_equal(a, b))
        return false;
    }
    if (left[row].evaluated_group_count != right[row].evaluated_group_count ||
        left[row].simultaneous_substitution_count !=
            right[row].simultaneous_substitution_count ||
        left[row].canonical_affine_row_evaluation_count !=
            right[row].canonical_affine_row_evaluation_count ||
        left[row].link_local_runtime_evaluation_count !=
            right[row].link_local_runtime_evaluation_count ||
        left[row].immutable_input_snapshot_count !=
            right[row].immutable_input_snapshot_count ||
        left[row].background_functional_evaluation_count !=
            right[row].background_functional_evaluation_count ||
        left[row].background_removal_count !=
            right[row].background_removal_count ||
        left[row].final_row_write_count != right[row].final_row_write_count)
      return false;
    for (std::size_t component = 0U; component < 3U; ++component) {
      if (bits(left[row].residual_before_wall[component]) !=
              bits(right[row].residual_before_wall[component]) ||
          bits(left[row].wall_contribution.convective[component]) !=
              bits(right[row].wall_contribution.convective[component]) ||
          bits(left[row].wall_contribution.pressure[component]) !=
              bits(right[row].wall_contribution.pressure[component]) ||
          bits(left[row].wall_contribution.viscous[component]) !=
              bits(right[row].wall_contribution.viscous[component]) ||
          bits(left[row].background_contribution.convective[component]) !=
              bits(right[row].background_contribution.convective[component]) ||
          bits(left[row].background_contribution.pressure[component]) !=
              bits(right[row].background_contribution.pressure[component]) ||
          bits(left[row].background_contribution.viscous[component]) !=
              bits(right[row].background_contribution.viscous[component]) ||
          bits(left[row].removed_background_contribution
                   .convective[component]) !=
              bits(right[row].removed_background_contribution
                       .convective[component]) ||
          bits(left[row].removed_background_contribution.pressure[component]) !=
              bits(right[row].removed_background_contribution
                       .pressure[component]) ||
          bits(left[row].removed_background_contribution.viscous[component]) !=
              bits(right[row].removed_background_contribution
                       .viscous[component]) ||
          bits(left[row].residual_after_wall[component]) !=
              bits(right[row].residual_after_wall[component]))
        return false;
    }
  }
  return true;
}

finite_volume::ImmersedResidualParts independent_wall_reaction(
    const finite_volume::ImmersedOperatorAdapter &adapter) {
  finite_volume::ImmersedResidualParts expected{};
  constexpr std::array<double, 3> velocity_scale{0.3, -0.2, 0.1};
  constexpr double mu = 0.01;
  for (const auto &row :
       finite_volume::test::ImmersedOperatorTestAccess::rows(adapter)) {
    for (const auto &link : row.links) {
      const double pressure = quadratic_value(link.pressure_quadrature_m);
      const auto gradient = quadratic_gradient(link.wall_intercept_m);
      std::array<std::array<double, 3>, 3> velocity_gradient{};
      for (std::size_t component = 0U; component < 3U; ++component) {
        velocity_gradient[component] = {velocity_scale[component] * gradient.x,
                                        velocity_scale[component] * gradient.y,
                                        velocity_scale[component] * gradient.z};
      }
      const double divergence = velocity_gradient[0][0] +
                                velocity_gradient[1][1] +
                                velocity_gradient[2][2];
      const std::array<double, 3> area{link.area_from_fluid_m2.x,
                                       link.area_from_fluid_m2.y,
                                       link.area_from_fluid_m2.z};
      for (std::size_t component = 0U; component < 3U; ++component) {
        expected.pressure[component] -= pressure * area[component];
        double traction = 0.0;
        for (std::size_t direction = 0U; direction < 3U; ++direction) {
          double stress = mu * (velocity_gradient[component][direction] +
                                velocity_gradient[direction][component]);
          if (component == direction)
            stress -= mu * (2.0 / 3.0) * divergence;
          traction += stress * area[direction];
        }
        expected.viscous[component] += traction;
      }
    }
  }
  return expected;
}

struct IndependentPressureHybridRow final {
  std::array<double, 3> reaction{};
  std::vector<std::array<double, 3>> reaction_by_row;
  std::size_t direction_mixing_link_count{};
  double legacy_scalar_mutation_max_difference{};
};

std::array<long double, 10>
uniform_cell_quadratic_basis_average(runtime::Int3 cell,
                                     const mesh::MeshGeometry &geometry) {
  const auto lo = geometry.vertex_position_m(cell);
  const auto hi = geometry.vertex_position_m(
      {cell.x + 1, cell.y + 1, cell.z + 1});
  const auto average = [](long double a, long double b) {
    return (a + b) / 2.0L;
  };
  const auto square_average = [](long double a, long double b) {
    return (a * a + a * b + b * b) / 3.0L;
  };
  const long double x = average(lo.x, hi.x);
  const long double y = average(lo.y, hi.y);
  const long double z = average(lo.z, hi.z);
  return {1.0L, x, y, z, square_average(lo.x, hi.x), x * y, x * z,
          square_average(lo.y, hi.y), y * z,
          square_average(lo.z, hi.z)};
}

double global_quadratic_mode_value(runtime::Real3 point, std::size_t mode) {
  const std::array<double, 10> basis{1.0,
                                     point.x,
                                     point.y,
                                     point.z,
                                     point.x * point.x,
                                     point.x * point.y,
                                     point.x * point.z,
                                     point.y * point.y,
                                     point.y * point.z,
                                     point.z * point.z};
  HUNDUN_CHECK(mode < basis.size());
  return basis[mode];
}

double uniform_face_quadratic_mode_average(
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    mesh::LocalFaceId face, std::size_t mode) {
  const auto logical = logical_face_vertices(topology.logical_face(face));
  std::array<runtime::Real3, 4> vertices{};
  for (std::size_t index = 0U; index < vertices.size(); ++index)
    vertices[index] = geometry.vertex_position_m(logical[index]);
  constexpr std::array<std::array<std::size_t, 3>, 2> triangles{{
      {{0U, 1U, 2U}},
      {{0U, 2U, 3U}},
  }};
  constexpr std::array<std::array<double, 3>, 3> barycentric{{
      {{2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0}},
      {{1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0}},
      {{1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0}},
  }};
  double integral = 0.0;
  double area = 0.0;
  for (const auto &triangle : triangles) {
    const auto &a = vertices[triangle[0]];
    const auto &b = vertices[triangle[1]];
    const auto &c = vertices[triangle[2]];
    const double triangle_area =
        0.5 * std::sqrt(dot(cross(subtract(b, a), subtract(c, a)),
                            cross(subtract(b, a), subtract(c, a))));
    double average = 0.0;
    for (const auto &weight : barycentric)
      average += global_quadratic_mode_value(
          add(add(multiply(a, weight[0]), multiply(b, weight[1])),
              multiply(c, weight[2])),
          mode);
    integral += triangle_area * average / 3.0;
    area += triangle_area;
  }
  return integral / area;
}

std::array<long double, 10>
global_quadratic_normal_derivative_basis(runtime::Real3 point,
                                         runtime::Real3 normal) {
  return {0.0L,
          normal.x,
          normal.y,
          normal.z,
          2.0L * point.x * normal.x,
          point.y * normal.x + point.x * normal.y,
          point.z * normal.x + point.x * normal.z,
          2.0L * point.y * normal.y,
          point.z * normal.y + point.y * normal.z,
          2.0L * point.z * normal.z};
}

std::vector<long double>
solve_dense_independent(std::vector<long double> matrix,
                        std::vector<long double> rhs) {
  const std::size_t size = rhs.size();
  HUNDUN_CHECK(matrix.size() == size * size);
  for (std::size_t pivot = 0U; pivot < size; ++pivot) {
    std::size_t selected = pivot;
    for (std::size_t row = pivot + 1U; row < size; ++row)
      if (std::abs(matrix[row * size + pivot]) >
          std::abs(matrix[selected * size + pivot]))
        selected = row;
    HUNDUN_CHECK(std::abs(matrix[selected * size + pivot]) > 1.0e-24L);
    if (selected != pivot) {
      for (std::size_t column = 0U; column < size; ++column)
        std::swap(matrix[pivot * size + column],
                  matrix[selected * size + column]);
      std::swap(rhs[pivot], rhs[selected]);
    }
    const long double diagonal = matrix[pivot * size + pivot];
    for (std::size_t row = pivot + 1U; row < size; ++row) {
      const long double factor = matrix[row * size + pivot] / diagonal;
      matrix[row * size + pivot] = 0.0L;
      for (std::size_t column = pivot + 1U; column < size; ++column)
        matrix[row * size + column] -=
            factor * matrix[pivot * size + column];
      rhs[row] -= factor * rhs[pivot];
    }
  }
  std::vector<long double> result(size, 0.0L);
  for (std::size_t reverse = size; reverse-- > 0U;) {
    long double value = rhs[reverse];
    for (std::size_t column = reverse + 1U; column < size; ++column)
      value -= matrix[reverse * size + column] * result[column];
    result[reverse] = value / matrix[reverse * size + reverse];
  }
  return result;
}

std::array<std::uint64_t, 10> direct_row_constrained_quadratic_bits(
    const finite_volume::test::ImmersedOperatorRowSnapshot &row,
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    const immersed::GhostStencilPlan &ghost_plan, bool reverse_inputs) {
  std::vector<mesh::GlobalCellId> donors;
  std::vector<finite_volume::test::ImmersedWallLinkSnapshot> links =
      row.links;
  if (reverse_inputs)
    std::reverse(links.begin(), links.end());
  for (const auto &link : links) {
    const auto &ids = immersed::detail::QuadraticReconstructionWeights::
        donor_global_ids(ghost_plan.reconstruction(link.id));
    donors.insert(donors.end(), ids.begin(), ids.end());
  }
  std::sort(donors.begin(), donors.end());
  donors.erase(std::unique(donors.begin(), donors.end()), donors.end());
  std::sort(links.begin(), links.end(),
            [](const auto &left, const auto &right) {
              return left.id < right.id;
            });
  HUNDUN_CHECK(donors.size() >= 10U && !links.empty() && links.size() < 10U);

  std::vector<std::array<long double, 10>> moments;
  moments.reserve(donors.size());
  for (const auto donor : donors)
    moments.push_back(uniform_cell_quadratic_basis_average(
        logical_cell(donor, topology.global_extent()), geometry));
  std::vector<std::array<long double, 10>> constraints;
  constraints.reserve(links.size());
  for (const auto &link : links)
    constraints.push_back(global_quadratic_normal_derivative_basis(
        link.wall_intercept_m, link.solid_to_fluid_normal));

  constexpr std::size_t basis_size = 10U;
  const std::size_t system_size = basis_size + constraints.size();
  std::vector<long double> kkt(system_size * system_size, 0.0L);
  for (std::size_t lhs = 0U; lhs < basis_size; ++lhs)
    for (std::size_t rhs = 0U; rhs < basis_size; ++rhs)
      for (const auto &moment : moments)
        kkt[lhs * system_size + rhs] += moment[lhs] * moment[rhs];
  for (std::size_t constraint = 0U; constraint < constraints.size();
       ++constraint)
    for (std::size_t basis = 0U; basis < basis_size; ++basis) {
      kkt[basis * system_size + basis_size + constraint] =
          constraints[constraint][basis];
      kkt[(basis_size + constraint) * system_size + basis] =
          constraints[constraint][basis];
    }

  std::array<std::uint64_t, basis_size> fingerprints{};
  for (std::size_t mode = 0U; mode < basis_size; ++mode) {
    std::vector<long double> rhs(system_size, 0.0L);
    for (std::size_t basis = 0U; basis < basis_size; ++basis)
      for (const auto &moment : moments)
        rhs[basis] += moment[basis] * moment[mode];
    for (std::size_t constraint = 0U; constraint < constraints.size();
         ++constraint)
      rhs[basis_size + constraint] = constraints[constraint][mode];
    const auto solution = solve_dense_independent(kkt, rhs);
    for (std::size_t coefficient = 0U; coefficient < basis_size;
         ++coefficient)
      HUNDUN_CHECK_NEAR(
          static_cast<double>(solution[coefficient]),
          coefficient == mode ? 1.0 : 0.0,
          4096.0 * std::numeric_limits<double>::epsilon());
    std::uint64_t fingerprint = 1469598103934665603ULL;
    for (std::size_t coefficient = 0U; coefficient < basis_size;
         ++coefficient) {
      const auto value = static_cast<double>(solution[coefficient]);
      const auto word = bits(value);
      fingerprint ^= word;
      fingerprint *= 1099511628211ULL;
    }
    fingerprints[mode] = fingerprint;
  }
  return fingerprints;
}

void prove_direct_row_constrained_quadratic_oracle(
    const finite_volume::ImmersedOperatorAdapter &adapter,
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    const immersed::GhostStencilPlan &ghost_plan) {
  bool observed_single_link = false;
  bool observed_multi_link = false;
  for (const auto &row :
       finite_volume::test::ImmersedOperatorTestAccess::rows(adapter)) {
    if (row.links.empty())
      continue;
    const auto canonical = direct_row_constrained_quadratic_bits(
        row, topology, geometry, ghost_plan, false);
    const auto reversed = direct_row_constrained_quadratic_bits(
        row, topology, geometry, ghost_plan, true);
    HUNDUN_CHECK(canonical == reversed);
    observed_single_link = observed_single_link || row.links.size() == 1U;
    observed_multi_link = observed_multi_link || row.links.size() > 1U;
  }
  HUNDUN_CHECK(observed_single_link);
  HUNDUN_CHECK(observed_multi_link);
}

std::array<double, 3> independent_complete_pressure_boundary_replacement(
    const finite_volume::test::ImmersedOperatorRowSnapshot &row,
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    const immersed::ImmersedDomain &domain) {
  const auto local = topology.find_local_cell(row.active_cell);
  HUNDUN_CHECK(local.has_value());
  std::array<double, 3> exact_complete_row{};
  std::array<double, 3> active_shared_row{};
  std::size_t incident_face_count = 0U;
  for (mesh::LocalFaceId face = 0U; face < topology.local_face_count(); ++face) {
    const auto owner = topology.owner(face);
    const auto neighbour = topology.neighbour(face);
    const bool row_is_owner = owner == *local;
    const bool row_is_neighbour = neighbour.has_value() && *neighbour == *local;
    if (!row_is_owner && !row_is_neighbour)
      continue;
    ++incident_face_count;
    const auto side = row_is_owner ? mesh::FaceSide::owner
                                   : mesh::FaceSide::neighbour;
    const auto area = geometry.face_area_vector_m2(face, side);
    const double exact_face =
        uniform_face_quadratic_average(topology, geometry, face);
    const std::array<double, 3> oriented_area{area.x, area.y, area.z};
    for (std::size_t component = 0U; component < 3U; ++component)
      exact_complete_row[component] +=
          exact_face * oriented_area[component];

    if (!neighbour.has_value())
      continue;
    const auto other = row_is_owner ? *neighbour : owner;
    if (domain.region(other) != immersed::CellRegion::fluid)
      continue;
    const auto owner_center = geometry.cell_center_m(owner);
    const auto neighbour_center = geometry.cell_center_m(*neighbour);
    const double distance = std::sqrt(dot(subtract(neighbour_center, owner_center),
                                          subtract(neighbour_center, owner_center)));
    const double owner_to_face = std::sqrt(
        dot(subtract(geometry.face_center_m(face), owner_center),
            subtract(geometry.face_center_m(face), owner_center)));
    HUNDUN_CHECK(distance > 0.0 && owner_to_face > 0.0 &&
                 owner_to_face < distance);
    const double owner_pressure =
        quadratic_average(topology.global_cell(owner), geometry);
    const double neighbour_pressure =
        quadratic_average(topology.global_cell(*neighbour), geometry);
    const double shared_face_pressure =
        ((distance - owner_to_face) / distance) * owner_pressure +
        (owner_to_face / distance) * neighbour_pressure;
    for (std::size_t component = 0U; component < 3U; ++component)
      active_shared_row[component] +=
          shared_face_pressure * oriented_area[component];
  }
  HUNDUN_CHECK(incident_face_count == 6U);
  std::array<double, 3> result{};
  for (std::size_t component = 0U; component < 3U; ++component)
    result[component] =
        exact_complete_row[component] - active_shared_row[component];
  return result;
}

std::array<double, 3> independent_complete_pressure_boundary_replacement_mode(
    const finite_volume::test::ImmersedOperatorRowSnapshot &row,
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    const immersed::ImmersedDomain &domain, std::size_t mode) {
  const auto local = topology.find_local_cell(row.active_cell);
  HUNDUN_CHECK(local.has_value());
  std::array<double, 3> result{};
  std::size_t incident_face_count = 0U;
  for (mesh::LocalFaceId face = 0U; face < topology.local_face_count(); ++face) {
    const auto owner = topology.owner(face);
    const auto neighbour = topology.neighbour(face);
    const bool row_is_owner = owner == *local;
    const bool row_is_neighbour = neighbour.has_value() && *neighbour == *local;
    if (!row_is_owner && !row_is_neighbour)
      continue;
    ++incident_face_count;
    const auto side = row_is_owner ? mesh::FaceSide::owner
                                   : mesh::FaceSide::neighbour;
    const auto area = geometry.face_area_vector_m2(face, side);
    const double exact_face = uniform_face_quadratic_mode_average(
        topology, geometry, face, mode);
    const double components[3]{area.x, area.y, area.z};
    for (std::size_t component = 0U; component < 3U; ++component)
      result[component] += exact_face * components[component];
    if (!neighbour.has_value())
      continue;
    const auto other = row_is_owner ? *neighbour : owner;
    if (domain.region(other) != immersed::CellRegion::fluid)
      continue;
    const auto owner_center = geometry.cell_center_m(owner);
    const auto neighbour_center = geometry.cell_center_m(*neighbour);
    const auto center_displacement = subtract(neighbour_center, owner_center);
    const auto face_displacement =
        subtract(geometry.face_center_m(face), owner_center);
    const double distance = std::sqrt(dot(center_displacement,
                                          center_displacement));
    const double owner_to_face =
        std::sqrt(dot(face_displacement, face_displacement));
    const auto owner_basis = uniform_cell_quadratic_basis_average(
        topology.global_cell(owner), geometry);
    const auto neighbour_basis = uniform_cell_quadratic_basis_average(
        topology.global_cell(*neighbour), geometry);
    const double shared =
        ((distance - owner_to_face) / distance) *
            static_cast<double>(owner_basis[mode]) +
        (owner_to_face / distance) *
            static_cast<double>(neighbour_basis[mode]);
    for (std::size_t component = 0U; component < 3U; ++component)
      result[component] -= shared * components[component];
  }
  HUNDUN_CHECK(incident_face_count == 6U);
  return result;
}

IndependentPressureHybridRow independent_pressure_hybrid_row(
    const finite_volume::ImmersedOperatorAdapter &adapter,
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    const immersed::GhostStencilPlan &ghost_plan,
    const runtime::FieldView<const double> &pressure,
    const std::vector<finite_volume::detail::ImmersedWallNormalGradient>
        *wall_normal_gradients = nullptr) {
  IndependentPressureHybridRow expected{};
  const auto extent = topology.global_extent();
  for (const auto &row :
       finite_volume::test::ImmersedOperatorTestAccess::rows(adapter)) {
    if (row.links.empty())
      continue;
    std::array<double, 3> row_reaction{};
    const auto fluid_logical = logical_cell(row.active_cell, extent);
    const auto fluid_local = topology.find_local_cell(row.active_cell);
    HUNDUN_CHECK(fluid_local.has_value());
    const auto fluid_center = geometry.cell_center_m(*fluid_local);
    std::array<runtime::Real3, 6> neighbour_centers{};
    for (std::size_t occurrence = 0U; occurrence < 6U; ++occurrence) {
      const auto offset = neighbour_offset(occurrence);
      const runtime::Int3 logical{fluid_logical.x + offset.x,
                                  fluid_logical.y + offset.y,
                                  fluid_logical.z + offset.z};
      const auto global = topology.global_cell_id(logical);
      const auto local = topology.find_local_cell(global);
      HUNDUN_CHECK(local.has_value());
      neighbour_centers[occurrence] = geometry.cell_center_m(*local);
    }
    for (const auto &link : row.links) {
      const auto solid_center = neighbour_centers[link.occurrence];
      const double distance =
          std::sqrt(dot(subtract(solid_center, fluid_center),
                        subtract(solid_center, fluid_center)));
      const double owner_to_face =
          std::sqrt(dot(subtract(link.pressure_quadrature_m, fluid_center),
                        subtract(link.pressure_quadrature_m, fluid_center)));
      HUNDUN_CHECK(distance > 0.0 && owner_to_face > 0.0 &&
                   owner_to_face < distance);
      immersed::LocalCoefficientRow background{};
      background.neighbour[link.occurrence] = owner_to_face / distance;
      background.diagonal = 1.0 - background.neighbour[link.occurrence];
      const auto transformed = independent_paper_map(
          background, link.normal_scale, link.solid_to_fluid_normal);
      const auto &reconstruction = ghost_plan.reconstruction(link.id);
      std::size_t nonzero_neighbours = 0U;
      const bool constrained = wall_normal_gradients != nullptr;
      const finite_volume::detail::ImmersedWallNormalGradient *condition =
          nullptr;
      if (constrained) {
        const auto found = std::lower_bound(
            wall_normal_gradients->begin(), wall_normal_gradients->end(),
            link.id, [](const auto &candidate, std::uint64_t id) {
              return candidate.link < id;
            });
        HUNDUN_CHECK(found != wall_normal_gradients->end() &&
                     found->link == link.id);
        condition = &*found;
      }
      const double wall_value =
          constrained
              ? immersed::detail::value_with_origin_normal_gradient(
                    reconstruction, link.wall_intercept_m, pressure, 0U,
                    condition->value)
              : reconstruction.value(link.wall_intercept_m, pressure, 0U);
      const auto wall_gradient =
          constrained
              ? immersed::detail::gradient_with_origin_normal_gradient(
                    reconstruction, link.wall_intercept_m, pressure, 0U,
                    condition->value)
              : reconstruction.gradient(link.wall_intercept_m, pressure, 0U);
      const auto remainder = [&](runtime::Real3 point_m) {
        const double value =
            constrained
                ? immersed::detail::value_with_origin_normal_gradient(
                      reconstruction, point_m, pressure, 0U, condition->value)
                : reconstruction.value(point_m, pressure, 0U);
        return value - wall_value -
               dot(wall_gradient, subtract(point_m, link.wall_intercept_m));
      };
      const auto evaluate_remainder_row =
          [&](const immersed::LocalCoefficientRow &coefficients) {
            double value = coefficients.diagonal * remainder(fluid_center) +
                           coefficients.source;
            for (std::size_t occurrence = 0U; occurrence < 6U; ++occurrence)
              value += coefficients.neighbour[occurrence] *
                       remainder(neighbour_centers[occurrence]);
            return value;
          };
      const double transformed_pressure =
          (constrained ? immersed::detail::value_with_origin_normal_gradient(
                             reconstruction, link.pressure_quadrature_m,
                             pressure, 0U, condition->value)
                       : reconstruction.value(link.pressure_quadrature_m,
                                              pressure, 0U)) +
          evaluate_remainder_row(transformed) -
          evaluate_remainder_row(background);
      for (std::size_t occurrence = 0U; occurrence < 6U; ++occurrence) {
        if (transformed.neighbour[occurrence] != 0.0)
          ++nonzero_neighbours;
      }
      if (nonzero_neighbours > 1U)
        ++expected.direction_mixing_link_count;
      const std::array<double, 3> area{link.area_from_fluid_m2.x,
                                       link.area_from_fluid_m2.y,
                                       link.area_from_fluid_m2.z};
      const double legacy_pressure =
          quadratic_value(link.pressure_quadrature_m);
      for (std::size_t component = 0U; component < 3U; ++component) {
        const double contribution = -transformed_pressure * area[component];
        expected.reaction[component] += contribution;
        row_reaction[component] += contribution;
        expected.legacy_scalar_mutation_max_difference = std::max(
            expected.legacy_scalar_mutation_max_difference,
            std::abs((legacy_pressure - transformed_pressure) *
                     area[component]));
      }
    }
    expected.reaction_by_row.push_back(row_reaction);
  }
  return expected;
}

struct IndependentPressureProjectedDefectRow final {
  std::array<double, 3> reaction{};
  std::size_t link_count{};
  std::size_t multi_link_row_count{};
};

enum class ProjectedDefectMutation {
  none,
  reversed_projection_sign,
  half_signed_measure,
  rotated_normal,
  face_wall_swap,
  omitted_defect,
  whole_vector_projection,
};

IndependentPressureProjectedDefectRow
independent_pressure_projected_defect_row(
    const finite_volume::ImmersedOperatorAdapter &adapter,
    const immersed::GhostStencilPlan &ghost_plan,
    const runtime::FieldView<const double> &pressure,
    const std::vector<finite_volume::detail::ImmersedWallNormalGradient>
        *wall_normal_gradients = nullptr,
    ProjectedDefectMutation mutation = ProjectedDefectMutation::none) {
  IndependentPressureProjectedDefectRow expected{};
  for (const auto &row :
       finite_volume::test::ImmersedOperatorTestAccess::rows(adapter)) {
    if (row.links.size() > 1U)
      ++expected.multi_link_row_count;
    for (const auto &link : row.links) {
      const auto &reconstruction = ghost_plan.reconstruction(link.id);
      const finite_volume::detail::ImmersedWallNormalGradient *condition =
          nullptr;
      if (wall_normal_gradients != nullptr) {
        const auto found = std::lower_bound(
            wall_normal_gradients->begin(), wall_normal_gradients->end(),
            link.id, [](const auto &candidate, std::uint64_t id) {
              return candidate.link < id;
            });
        HUNDUN_CHECK(found != wall_normal_gradients->end() &&
                     found->link == link.id);
        condition = &*found;
      }
      const auto value_at = [&](runtime::Real3 point_m) {
        return condition == nullptr
                   ? reconstruction.value(point_m, pressure, 0U)
                   : immersed::detail::value_with_origin_normal_gradient(
                         reconstruction, point_m, pressure, 0U,
                         condition->value);
      };
      const double face_pressure = value_at(link.pressure_quadrature_m);
      const double wall_pressure = value_at(link.wall_intercept_m);
      auto projected_normal = link.solid_to_fluid_normal;
      if (mutation == ProjectedDefectMutation::rotated_normal)
        projected_normal = {projected_normal.y, projected_normal.z,
                            projected_normal.x};
      double projected_scale = -link.signed_wall_measure_m2;
      if (mutation == ProjectedDefectMutation::reversed_projection_sign)
        projected_scale = -projected_scale;
      if (mutation == ProjectedDefectMutation::half_signed_measure)
        projected_scale *= 0.5;
      const auto projected_area =
          multiply(projected_normal, projected_scale);
      const std::array<double, 3> background_area{
          link.area_from_fluid_m2.x, link.area_from_fluid_m2.y,
          link.area_from_fluid_m2.z};
      const std::array<double, 3> projected{
          projected_area.x, projected_area.y, projected_area.z};
      for (std::size_t component = 0U; component < 3U; ++component) {
        double wall_row = face_pressure * background_area[component] +
                          (wall_pressure - face_pressure) *
                              projected[component];
        if (mutation == ProjectedDefectMutation::face_wall_swap)
          wall_row = wall_pressure * background_area[component] +
                     (face_pressure - wall_pressure) * projected[component];
        if (mutation == ProjectedDefectMutation::omitted_defect)
          wall_row = face_pressure * background_area[component];
        if (mutation == ProjectedDefectMutation::whole_vector_projection)
          wall_row = wall_pressure * projected[component];
        expected.reaction[component] -= wall_row;
      }
      ++expected.link_count;
    }
  }
  return expected;
}

std::array<double, 3> independent_viscous_traction(
    const immersed::QuadraticReconstruction &reconstruction,
    runtime::Real3 point_m, const runtime::FieldView<const double> &velocity,
    double mu, runtime::Real3 area_from_fluid_m2) {
  std::array<runtime::Real3, 3> gradient{};
  for (std::size_t component = 0U; component < 3U; ++component)
    gradient[component] = immersed::detail::gradient_with_origin_constraint(
        reconstruction, point_m, velocity, component, 0.0);
  const auto component = [](runtime::Real3 value, std::size_t direction) {
    return direction == 0U ? value.x : direction == 1U ? value.y : value.z;
  };
  const std::array<double, 3> area{area_from_fluid_m2.x, area_from_fluid_m2.y,
                                   area_from_fluid_m2.z};
  const double divergence = gradient[0].x + gradient[1].y + gradient[2].z;
  std::array<double, 3> traction{};
  for (std::size_t row = 0U; row < 3U; ++row)
    for (std::size_t direction = 0U; direction < 3U; ++direction) {
      double stress = mu * (component(gradient[row], direction) +
                            component(gradient[direction], row));
      if (row == direction)
        stress -= mu * (2.0 / 3.0) * divergence;
      traction[row] += stress * area[direction];
    }
  return traction;
}

std::array<double, 3> independent_viscous_hybrid_row(
    const finite_volume::ImmersedOperatorAdapter &adapter,
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    const immersed::GhostStencilPlan &ghost_plan,
    const runtime::FieldView<const double> &velocity, double mu) {
  std::array<double, 3> expected{};
  const auto extent = topology.global_extent();
  for (const auto &row :
       finite_volume::test::ImmersedOperatorTestAccess::rows(adapter)) {
    if (row.links.empty())
      continue;
    const auto fluid_logical = logical_cell(row.active_cell, extent);
    const auto fluid_local = topology.find_local_cell(row.active_cell);
    HUNDUN_CHECK(fluid_local.has_value());
    const auto fluid_center = geometry.cell_center_m(*fluid_local);
    std::array<runtime::Real3, 6> neighbour_centers{};
    for (std::size_t occurrence = 0U; occurrence < 6U; ++occurrence) {
      const auto offset = neighbour_offset(occurrence);
      const runtime::Int3 logical{fluid_logical.x + offset.x,
                                  fluid_logical.y + offset.y,
                                  fluid_logical.z + offset.z};
      const auto local =
          topology.find_local_cell(topology.global_cell_id(logical));
      HUNDUN_CHECK(local.has_value());
      neighbour_centers[occurrence] = geometry.cell_center_m(*local);
    }
    for (const auto &link : row.links) {
      const auto solid_center = neighbour_centers[link.occurrence];
      const double distance =
          std::sqrt(dot(subtract(solid_center, fluid_center),
                        subtract(solid_center, fluid_center)));
      const double owner_to_face =
          std::sqrt(dot(subtract(link.pressure_quadrature_m, fluid_center),
                        subtract(link.pressure_quadrature_m, fluid_center)));
      immersed::LocalCoefficientRow background{};
      background.neighbour[link.occurrence] = owner_to_face / distance;
      background.diagonal = 1.0 - background.neighbour[link.occurrence];
      const auto transformed = independent_paper_map(
          background, link.normal_scale, link.solid_to_fluid_normal);
      const auto &reconstruction = ghost_plan.reconstruction(link.id);
      const auto wall_traction =
          independent_viscous_traction(reconstruction, link.wall_intercept_m,
                                       velocity, mu, link.area_from_fluid_m2);
      const auto diagonal_traction = independent_viscous_traction(
          reconstruction, fluid_center, velocity, mu, link.area_from_fluid_m2);
      std::array<double, 3> transformed_traction = wall_traction;
      for (std::size_t component = 0U; component < 3U; ++component)
        transformed_traction[component] +=
            (transformed.diagonal - background.diagonal) *
            diagonal_traction[component];
      for (std::size_t occurrence = 0U; occurrence < 6U; ++occurrence) {
        const auto neighbour_traction = independent_viscous_traction(
            reconstruction, neighbour_centers[occurrence], velocity, mu,
            link.area_from_fluid_m2);
        for (std::size_t component = 0U; component < 3U; ++component)
          transformed_traction[component] +=
              (transformed.neighbour[occurrence] -
               background.neighbour[occurrence]) *
              neighbour_traction[component];
      }
      for (std::size_t component = 0U; component < 3U; ++component)
        expected[component] += transformed_traction[component];
    }
  }
  return expected;
}

std::array<double, 3> independent_shared_row_viscous_row(
    const finite_volume::ImmersedOperatorAdapter &adapter,
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    const immersed::GhostStencilPlan &ghost_plan,
    const runtime::FieldView<const double> &velocity, double mu) {
  std::array<double, 3> expected{};
  const auto extent = topology.global_extent();
  for (const auto &row :
       finite_volume::test::ImmersedOperatorTestAccess::rows(adapter)) {
    if (row.links.empty())
      continue;
    const auto fluid_logical = logical_cell(row.active_cell, extent);
    const auto fluid_local = topology.find_local_cell(row.active_cell);
    HUNDUN_CHECK(fluid_local.has_value());
    const auto fluid_center = geometry.cell_center_m(*fluid_local);
    std::array<runtime::Real3, 6> neighbour_centers{};
    for (std::size_t occurrence = 0U; occurrence < 6U; ++occurrence) {
      const auto offset = neighbour_offset(occurrence);
      const runtime::Int3 logical{fluid_logical.x + offset.x,
                                  fluid_logical.y + offset.y,
                                  fluid_logical.z + offset.z};
      const auto local =
          topology.find_local_cell(topology.global_cell_id(logical));
      HUNDUN_CHECK(local.has_value());
      neighbour_centers[occurrence] = geometry.cell_center_m(*local);
    }
    std::vector<mesh::GlobalCellId> donor_ids;
    for (const auto &link : row.links) {
      const auto &ids =
          immersed::detail::QuadraticReconstructionWeights::donor_global_ids(
              ghost_plan.reconstruction(link.id));
      donor_ids.insert(donor_ids.end(), ids.begin(), ids.end());
    }
    std::sort(donor_ids.begin(), donor_ids.end());
    donor_ids.erase(std::unique(donor_ids.begin(), donor_ids.end()),
                    donor_ids.end());
    std::vector<runtime::Int3> donor_cells;
    donor_cells.reserve(donor_ids.size());
    for (const auto donor : donor_ids)
      donor_cells.push_back(logical_cell(donor, extent));
    immersed::detail::BoundaryAuthorityCoverageScope coverage;
    const auto shared_reconstruction = immersed::QuadraticReconstruction::create(
        fluid_center, runtime::Real3{1.0, 0.0, 0.0},
        runtime::Real3{0.0, 1.0, 0.0}, runtime::Real3{0.0, 0.0, 1.0},
        std::cbrt(geometry.cell_volume_m3(*fluid_local)), fluid_logical,
        donor_cells, topology, geometry);
    for (const auto &link : row.links) {
      const auto solid_center = neighbour_centers[link.occurrence];
      const double distance =
          std::sqrt(dot(subtract(solid_center, fluid_center),
                        subtract(solid_center, fluid_center)));
      const double owner_to_face =
          std::sqrt(dot(subtract(link.pressure_quadrature_m, fluid_center),
                        subtract(link.pressure_quadrature_m, fluid_center)));
      immersed::LocalCoefficientRow background{};
      background.neighbour[link.occurrence] = owner_to_face / distance;
      background.diagonal = 1.0 - background.neighbour[link.occurrence];
      const auto transformed = independent_paper_map(
          background, link.normal_scale, link.solid_to_fluid_normal);
      const auto wall_traction = independent_viscous_traction(
          shared_reconstruction, link.wall_intercept_m, velocity, mu,
          link.area_from_fluid_m2);
      const auto diagonal_traction = independent_viscous_traction(
          shared_reconstruction, fluid_center, velocity, mu,
          link.area_from_fluid_m2);
      std::array<double, 3> transformed_traction = wall_traction;
      for (std::size_t component = 0U; component < 3U; ++component)
        transformed_traction[component] +=
            (transformed.diagonal - background.diagonal) *
            diagonal_traction[component];
      for (std::size_t occurrence = 0U; occurrence < 6U; ++occurrence) {
        const auto neighbour_traction = independent_viscous_traction(
            shared_reconstruction, neighbour_centers[occurrence], velocity, mu,
            link.area_from_fluid_m2);
        for (std::size_t component = 0U; component < 3U; ++component)
          transformed_traction[component] +=
              (transformed.neighbour[occurrence] -
               background.neighbour[occurrence]) *
              neighbour_traction[component];
      }
      for (std::size_t component = 0U; component < 3U; ++component)
        expected[component] += transformed_traction[component];
    }
  }
  return expected;
}

void run_case(const runtime::MpiContext &mpi, bool warped) {
  const auto grid = process_grid(mpi.size());
  auto decomposition = runtime::StructuredDecomposition::create(
      mpi, kExtent, {false, false, false}, runtime::DecompositionOptions{grid});
  mesh::MeshTopology topology(decomposition);
  mesh::MeshGeometry geometry =
      warped
          ? mesh::MeshGeometry(
                topology, mesh::AnalyticWarpedBoxMapping{{0.0, 0.0, 0.0},
                                                         {1.0, 1.0, 1.0},
                                                         {0.02, -0.015, 0.01}})
          : mesh::MeshGeometry(topology, mesh::UniformBoxMapping{
                                             {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}});
  auto boundaries = boundary::BoundaryRegistry::create(
      case_config(mpi.size(), grid, warped), topology);
  FixtureFile fixture(mpi);
  const auto surface =
      immersed::ImmersedSurface::load_collective(fixture.path(), 0.4, mpi, 0);
  const auto query = immersed::SurfaceQuery::create(surface);
  const auto domain = immersed::ImmersedDomain::create(
      surface, query, config::ImmersedFluidSide::outside, topology, geometry,
      boundaries, mpi);
  const auto ghost_plan = immersed::GhostStencilPlan::create(
      surface, query, domain, topology, geometry, decomposition, mpi);
  const int ghost = static_cast<int>(ghost_plan.maximum_halo_reach());
  auto halo = runtime::HaloExchange::create(
      decomposition, runtime::ExchangePlan::create(
                         decomposition, decomposition.local_extent(), ghost));
  auto reconstruction = finite_volume::ImmersedReconstruction::create(
      topology, geometry, boundaries, domain, ghost_plan, decomposition, mpi,
      halo);
  const immersed::LocalFlowPatternTransform transform;
  auto adapter = finite_volume::ImmersedOperatorAdapter::create(
      topology, geometry, domain, ghost_plan, transform, reconstruction);

  FixtureFile wall_authority_fixture(mpi, false, true);
  const auto wall_authority_surface = immersed::ImmersedSurface::load_collective(
      wall_authority_fixture.path(), 0.4, mpi, 0);
  const auto wall_authority_query =
      immersed::SurfaceQuery::create(wall_authority_surface);
  const auto wall_authority_domain = immersed::ImmersedDomain::create(
      wall_authority_surface, wall_authority_query,
      config::ImmersedFluidSide::outside, topology, geometry, boundaries, mpi);
  const auto wall_authority_ghost_plan = immersed::GhostStencilPlan::create(
      wall_authority_surface, wall_authority_query, wall_authority_domain,
      topology, geometry, decomposition, mpi);
  HUNDUN_CHECK(wall_authority_ghost_plan.maximum_halo_reach() <=
               ghost_plan.maximum_halo_reach());
  auto wall_authority_halo = runtime::HaloExchange::create(
      decomposition, runtime::ExchangePlan::create(
                         decomposition, decomposition.local_extent(), ghost));
  auto wall_authority_reconstruction =
      finite_volume::ImmersedReconstruction::create(
          topology, geometry, boundaries, wall_authority_domain,
          wall_authority_ghost_plan, decomposition, mpi, wall_authority_halo);
  auto wall_authority_adapter = finite_volume::ImmersedOperatorAdapter::create(
      topology, geometry, wall_authority_domain, wall_authority_ghost_plan,
      transform, wall_authority_reconstruction);
  const auto wall_authority_plan = immersed::WallQuadraturePlan::create(
      wall_authority_surface, wall_authority_query, wall_authority_domain,
      topology, geometry, mpi);
  prove_interface_pressure_difference_snapshot(adapter);
  test_row_plan_snapshot(adapter, transform);
  if (!warped)
    prove_direct_row_constrained_quadratic_oracle(
        adapter, topology, geometry, ghost_plan);
  FixtureFile permuted_fixture(mpi, true);
  const auto permuted_surface = immersed::ImmersedSurface::load_collective(
      permuted_fixture.path(), 0.4, mpi, 0);
  const auto permuted_query =
      immersed::SurfaceQuery::create(permuted_surface);
  const auto permuted_domain = immersed::ImmersedDomain::create(
      permuted_surface, permuted_query, config::ImmersedFluidSide::outside,
      topology, geometry, boundaries, mpi);
  const auto permuted_ghost_plan = immersed::GhostStencilPlan::create(
      permuted_surface, permuted_query, permuted_domain, topology, geometry,
      decomposition, mpi);
  HUNDUN_CHECK(permuted_ghost_plan.maximum_halo_reach() ==
               ghost_plan.maximum_halo_reach());
  auto permuted_halo = runtime::HaloExchange::create(
      decomposition, runtime::ExchangePlan::create(
                         decomposition, decomposition.local_extent(), ghost));
  auto permuted_reconstruction = finite_volume::ImmersedReconstruction::create(
      topology, geometry, boundaries, permuted_domain, permuted_ghost_plan,
      decomposition, mpi, permuted_halo);
  auto permuted_adapter = finite_volume::ImmersedOperatorAdapter::create(
      topology, geometry, permuted_domain, permuted_ghost_plan, transform,
      permuted_reconstruction);
  const auto canonical_rows =
      finite_volume::test::ImmersedOperatorTestAccess::rows(adapter);
  const auto permuted_rows =
      finite_volume::test::ImmersedOperatorTestAccess::rows(permuted_adapter);
  if (!operator_rows_bits_equal(canonical_rows, permuted_rows)) {
    std::cerr << "row_permutation_diag canonical_rows="
              << canonical_rows.size()
              << " permuted_rows=" << permuted_rows.size() << '\n';
    const auto count = std::min(canonical_rows.size(), permuted_rows.size());
    for (std::size_t row = 0U; row < count; ++row) {
      if (canonical_rows[row].active_cell != permuted_rows[row].active_cell ||
          canonical_rows[row].row_replacement_fingerprint !=
              permuted_rows[row].row_replacement_fingerprint ||
          canonical_rows[row].links.size() != permuted_rows[row].links.size() ||
          canonical_rows[row].covered_physical_terms.size() !=
              permuted_rows[row].covered_physical_terms.size()) {
        std::cerr << "row_permutation_first row=" << row
                  << " canonical_cell=" << canonical_rows[row].active_cell
                  << " permuted_cell=" << permuted_rows[row].active_cell
                  << " canonical_fingerprint="
                  << canonical_rows[row].row_replacement_fingerprint
                  << " permuted_fingerprint="
                  << permuted_rows[row].row_replacement_fingerprint
                  << " canonical_links=" << canonical_rows[row].links.size()
                  << " permuted_links=" << permuted_rows[row].links.size()
                  << '\n';
        const auto link_count = std::min(canonical_rows[row].links.size(),
                                         permuted_rows[row].links.size());
        for (std::size_t link = 0U; link < link_count; ++link) {
          const auto &a = canonical_rows[row].links[link];
          const auto &b = permuted_rows[row].links[link];
          std::cerr << "row_permutation_link slot=" << link
                    << " ids=" << a.id << ',' << b.id
                    << " occurrences=" << a.occurrence << ',' << b.occurrence
                    << " normal_bits=" << bits(a.solid_to_fluid_normal.x)
                    << ',' << bits(b.solid_to_fluid_normal.x) << ':'
                    << bits(a.solid_to_fluid_normal.y) << ','
                    << bits(b.solid_to_fluid_normal.y) << ':'
                    << bits(a.solid_to_fluid_normal.z) << ','
                    << bits(b.solid_to_fluid_normal.z)
                    << " intercept_bits=" << bits(a.wall_intercept_m.x) << ','
                    << bits(b.wall_intercept_m.x) << ':'
                    << bits(a.wall_intercept_m.y) << ','
                    << bits(b.wall_intercept_m.y) << ':'
                    << bits(a.wall_intercept_m.z) << ','
                    << bits(b.wall_intercept_m.z) << '\n';
        }
        break;
      }
    }
  }
  HUNDUN_CHECK(operator_rows_bits_equal(canonical_rows, permuted_rows));
  std::uint64_t planned_link_count = 0U;
  for (const auto &row :
       finite_volume::test::ImmersedOperatorTestAccess::rows(adapter))
    planned_link_count += static_cast<std::uint64_t>(row.links.size());
  std::uint64_t published_link_count =
      static_cast<std::uint64_t>(domain.links().size());
  HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, &planned_link_count, 1, MPI_UINT64_T,
                             MPI_SUM, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, &published_link_count, 1,
                             MPI_UINT64_T, MPI_SUM, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(planned_link_count == published_link_count);

  runtime::FieldRegistry registry;
  const auto ids = declare_fields(registry, ghost);
  registry.freeze();
  runtime::FieldStorage storage(
      registry, runtime::FieldLayoutSet{decomposition.local_extent(),
                                        topology.local_face_count()});
  runtime::FieldAccessPlan access(registry);
  const std::array<runtime::FieldId, 12> fields{
      ids.velocity,          ids.pressure,          ids.scalar,
      ids.velocity_gradient, ids.scalar_gradient,   ids.face_velocity,
      ids.face_scalar,       ids.viscosity,         ids.gamma,
      ids.mass_flux,         ids.momentum_residual, ids.scalar_residual};
  for (const auto field : fields)
    access.declare_access(kPhase, kActor, field,
                          runtime::AccessMode::read_write);
  access.freeze();

  auto velocity = storage.view<double>(ids.velocity);
  auto pressure = storage.view<double>(ids.pressure);
  auto scalar = storage.view<double>(ids.scalar);
  const auto box = decomposition.owned_box();
  const auto local = decomposition.local_extent();
  for (int k = -ghost; k < local.z + ghost; ++k)
    for (int j = -ghost; j < local.y + ghost; ++j)
      for (int i = -ghost; i < local.x + ghost; ++i) {
        const runtime::Int3 global{box.begin.x + i, box.begin.y + j,
                                   box.begin.z + k};
        double value = 0.0;
        if (global.x >= 0 && global.y >= 0 && global.z >= 0 &&
            global.x < kExtent.x && global.y < kExtent.y &&
            global.z < kExtent.z)
          value = quadratic_average(global, geometry);
        pressure(i, j, k, 0) = value;
        scalar(i, j, k, 0) = 0.75 * value;
        velocity(i, j, k, 0) = 0.3 * value;
        velocity(i, j, k, 1) = -0.2 * value;
        velocity(i, j, k, 2) = 0.1 * value;
      }
  prove_wall_pressure_authority_matches_operator_row(
      wall_authority_adapter, wall_authority_plan, wall_authority_ghost_plan,
      mpi);
  for (mesh::LocalCellId cell = 0U; cell < topology.local_cell_count();
       ++cell) {
    if (domain.region(cell) != immersed::CellRegion::solid)
      continue;
    const auto global = topology.global_cell(cell);
    const int i = global.x - box.begin.x;
    const int j = global.y - box.begin.y;
    const int k = global.z - box.begin.z;
    if (i < -ghost || i >= local.x + ghost || j < -ghost ||
        j >= local.y + ghost || k < -ghost || k >= local.z + ghost)
      continue;
    pressure(i, j, k, 0) = 11111.0;
    scalar(i, j, k, 0) = -22222.0;
    velocity(i, j, k, 0) = 33333.0;
    velocity(i, j, k, 1) = -44444.0;
    velocity(i, j, k, 2) = 55555.0;
  }
  auto mass_flux_writer =
      storage.acquire_face_write<double>(access, kPhase, kActor, ids.mass_flux);
  auto viscosity =
      storage.acquire_face_write<double>(access, kPhase, kActor, ids.viscosity);
  auto gamma =
      storage.acquire_face_write<double>(access, kPhase, kActor, ids.gamma);
  const auto initialize_mass_flux = [&] {
    for (mesh::LocalFaceId face = 0U; face < topology.local_face_count();
         ++face) {
      const auto neighbour = topology.neighbour(face);
      const bool active_pair =
          neighbour.has_value() &&
          domain.region(topology.owner(face)) == immersed::CellRegion::fluid &&
          domain.region(*neighbour) == immersed::CellRegion::fluid;
      mass_flux_writer(face, 0) =
          active_pair ? (topology.global_face_id(face) % 2U == 0U ? 0.2 : -0.15)
                      : 0.0;
    }
  };
  initialize_mass_flux();
  for (mesh::LocalFaceId face = 0U; face < topology.local_face_count();
       ++face) {
    viscosity(face, 0) = 0.01;
    gamma(face, 0) = 0.02;
  }
  auto mass_flux = finite_volume::FaceMassFlux::acquire(
      registry, storage, access, kPhase, kActor, ids.mass_flux, topology);

  reconstruction.compute_gradient(
      finite_volume::GradientScheme::weighted_least_squares,
      finite_volume::FiniteVolumeQuantity::velocity(),
      binding(storage, access, ids.velocity),
      binding(storage, access, ids.velocity_gradient));
  halo.begin(storage, ids.velocity_gradient);
  halo.wait(storage, ids.velocity_gradient);
  reconstruction.compute_gradient(
      finite_volume::GradientScheme::weighted_least_squares,
      finite_volume::FiniteVolumeQuantity::scalar(0U),
      binding(storage, access, ids.scalar),
      binding(storage, access, ids.scalar_gradient));
  halo.begin(storage, ids.scalar_gradient);
  halo.wait(storage, ids.scalar_gradient);
  reconstruction.reconstruct_momentum_faces(
      mass_flux, binding(storage, access, ids.velocity),
      binding(storage, access, ids.face_velocity));
  reconstruction.reconstruct_transport_faces(
      finite_volume::FiniteVolumeQuantity::scalar(0U), mass_flux,
      binding(storage, access, ids.scalar),
      binding(storage, access, ids.face_scalar));

  auto momentum_residual = storage.view<double>(ids.momentum_residual);
  auto scalar_residual = storage.view<double>(ids.scalar_residual);
  fill_zero(momentum_residual);
  fill_zero(scalar_residual);
  const auto const_velocity =
      static_cast<const runtime::FieldStorage &>(storage).view<double>(
          ids.velocity);
  const auto const_pressure =
      static_cast<const runtime::FieldStorage &>(storage).view<double>(
          ids.pressure);
  const auto const_scalar =
      static_cast<const runtime::FieldStorage &>(storage).view<double>(
          ids.scalar);
  const auto velocity_gradient =
      static_cast<const runtime::FieldStorage &>(storage).view<double>(
          ids.velocity_gradient);
  const auto scalar_gradient =
      static_cast<const runtime::FieldStorage &>(storage).view<double>(
          ids.scalar_gradient);
  const auto face_velocity =
      static_cast<const runtime::FieldStorage &>(storage)
          .acquire_face_read<double>(access, kPhase, kActor, ids.face_velocity);
  const auto face_scalar =
      static_cast<const runtime::FieldStorage &>(storage)
          .acquire_face_read<double>(access, kPhase, kActor, ids.face_scalar);
  const auto const_viscosity =
      static_cast<const runtime::FieldStorage &>(storage)
          .acquire_face_read<double>(access, kPhase, kActor, ids.viscosity);
  const auto const_gamma =
      static_cast<const runtime::FieldStorage &>(storage)
          .acquire_face_read<double>(access, kPhase, kActor, ids.gamma);

  adapter.accumulate_momentum(mass_flux, face_velocity, const_velocity,
                              const_pressure, velocity_gradient,
                              const_viscosity, momentum_residual);
  std::uint64_t local_immersed_rows = 0U;
  for (const auto &row :
       finite_volume::test::ImmersedOperatorTestAccess::rows(adapter))
    if (!row.links.empty())
      ++local_immersed_rows;
  HUNDUN_CHECK(finite_volume::test::ImmersedOperatorTestAccess::
                   last_wall_functional_evaluation_count(adapter) ==
               local_immersed_rows);
  HUNDUN_CHECK(finite_volume::test::ImmersedOperatorTestAccess::
                   last_boundary_row_evaluations(adapter)
                       .size() == local_immersed_rows);
  const auto first_row_evaluations =
      test_boundary_row_evaluations(adapter, geometry);
  const auto first_momentum_before_permutation = owned_active_bits(
      topology, domain,
      static_cast<const runtime::FieldStorage &>(storage).view<double>(
          ids.momentum_residual),
      3U, box);
  const auto first_report_before_permutation = adapter.report();
  fill_zero(momentum_residual);
  permuted_adapter.accumulate_momentum(
      mass_flux, face_velocity, const_velocity, const_pressure,
      velocity_gradient, const_viscosity, momentum_residual);
  HUNDUN_CHECK(first_momentum_before_permutation ==
               owned_active_bits(
                   topology, permuted_domain,
                   static_cast<const runtime::FieldStorage &>(storage)
                       .view<double>(ids.momentum_residual),
                   3U, box));
  HUNDUN_CHECK(report_bits_equal(first_report_before_permutation,
                                 permuted_adapter.report()));
  HUNDUN_CHECK(boundary_evaluations_bits_equal(
      first_row_evaluations,
      finite_volume::test::ImmersedOperatorTestAccess::
          last_boundary_row_evaluations(permuted_adapter)));
  fill_zero(momentum_residual);
  adapter.accumulate_momentum(mass_flux, face_velocity, const_velocity,
                              const_pressure, velocity_gradient,
                              const_viscosity, momentum_residual);
  HUNDUN_CHECK(boundary_evaluations_bits_equal(
      first_row_evaluations,
      finite_volume::test::ImmersedOperatorTestAccess::
          last_boundary_row_evaluations(adapter)));
  adapter.accumulate_transport(finite_volume::FiniteVolumeQuantity::scalar(0U),
                               mass_flux, face_scalar, const_scalar,
                               scalar_gradient, const_gamma, scalar_residual);
  const auto report = adapter.report();
  HUNDUN_CHECK(report.active_row_count ==
               domain.active_cells().owned_active_count());
  HUNDUN_CHECK(report.replacement_group_count > 0U);
  std::uint64_t global_simultaneous = report.simultaneous_substitution_count;
  HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, &global_simultaneous, 1,
                             MPI_UINT64_T, MPI_SUM, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(global_simultaneous > 0U);
  for (const double value : report.budget_reaction_N.convective)
    HUNDUN_CHECK(bits(value) == bits(0.0));
  std::array<double, 3> actual_pressure_wall_contribution{};
  std::array<double, 3> actual_viscous_wall_contribution{};
  for (const auto &row : first_row_evaluations)
    for (std::size_t component = 0U; component < 3U; ++component) {
      actual_pressure_wall_contribution[component] +=
          row.wall_contribution.pressure[component];
      actual_viscous_wall_contribution[component] +=
          row.wall_contribution.viscous[component];
    }
  for (std::size_t component = 0U; component < 3U; ++component) {
    HUNDUN_CHECK(bits(report.budget_reaction_N.pressure[component]) ==
                 bits(-actual_pressure_wall_contribution[component]));
    HUNDUN_CHECK(bits(report.budget_reaction_N.viscous[component]) ==
                 bits(-actual_viscous_wall_contribution[component]));
  }
  const auto expected_reaction = independent_wall_reaction(adapter);
  for (std::size_t component = 0U; component < 3U; ++component) {
    HUNDUN_CHECK(std::isfinite(expected_reaction.pressure[component]));
    HUNDUN_CHECK(std::isfinite(expected_reaction.viscous[component]));
  }
  const auto expected_hybrid_row = independent_pressure_hybrid_row(
      adapter, topology, geometry, ghost_plan, const_pressure);
  const auto expected_projected_defect =
      independent_pressure_projected_defect_row(adapter, ghost_plan,
                                                const_pressure);
  const auto expected_viscous_hybrid_row = independent_viscous_hybrid_row(
      adapter, topology, geometry, ghost_plan, const_velocity, 0.01);
  const auto expected_shared_viscous_row = independent_shared_row_viscous_row(
      adapter, topology, geometry, ghost_plan, const_velocity, 0.01);
  double shared_vs_product_viscous = 0.0;
  double shared_vs_perlink_viscous = 0.0;
  for (std::size_t component = 0U; component < 3U; ++component) {
    shared_vs_product_viscous =
        std::max(shared_vs_product_viscous,
                 std::abs(report.budget_reaction_N.viscous[component] -
                          expected_shared_viscous_row[component]));
    shared_vs_perlink_viscous =
        std::max(shared_vs_perlink_viscous,
                 std::abs(expected_viscous_hybrid_row[component] -
                          expected_shared_viscous_row[component]));
  }
  std::cerr << "SHAREDVISC rank=" << mpi.rank()
            << " product=" << report.budget_reaction_N.viscous[0] << ','
            << report.budget_reaction_N.viscous[1] << ','
            << report.budget_reaction_N.viscous[2]
            << " perlink=" << expected_viscous_hybrid_row[0] << ','
            << expected_viscous_hybrid_row[1] << ','
            << expected_viscous_hybrid_row[2]
            << " shared=" << expected_shared_viscous_row[0] << ','
            << expected_shared_viscous_row[1] << ','
            << expected_shared_viscous_row[2]
            << " shared_vs_product=" << shared_vs_product_viscous
            << " shared_vs_perlink=" << shared_vs_perlink_viscous << '\n';
  HUNDUN_CHECK(expected_hybrid_row.direction_mixing_link_count > 0U);
  HUNDUN_CHECK(expected_projected_defect.link_count > 0U);
  HUNDUN_CHECK(expected_projected_defect.multi_link_row_count > 0U);
  if (mpi.size() == 1) {
    constexpr std::array<ProjectedDefectMutation, 6> mutations{
        ProjectedDefectMutation::reversed_projection_sign,
        ProjectedDefectMutation::half_signed_measure,
        ProjectedDefectMutation::rotated_normal,
        ProjectedDefectMutation::face_wall_swap,
        ProjectedDefectMutation::omitted_defect,
        ProjectedDefectMutation::whole_vector_projection};
    for (const auto mutation : mutations) {
      const auto mutated = independent_pressure_projected_defect_row(
          adapter, ghost_plan, const_pressure, nullptr, mutation);
      double difference = 0.0;
      for (std::size_t component = 0U; component < 3U; ++component)
        difference =
            std::max(difference,
                     std::abs(mutated.reaction[component] -
                              expected_projected_defect.reaction[component]));
      HUNDUN_CHECK(difference > 1.0e-8);
    }
  }
  double retired_hybrid_difference = 0.0;
  for (std::size_t component = 0U; component < 3U; ++component)
    {
      retired_hybrid_difference =
          std::max(retired_hybrid_difference,
                   std::abs(expected_projected_defect.reaction[component] -
                            expected_hybrid_row.reaction[component]));
    }
  HUNDUN_CHECK(expected_hybrid_row.legacy_scalar_mutation_max_difference >
               1.0e-8);
  HUNDUN_CHECK(retired_hybrid_difference > 1.0e-8);
  double viscous_reaction_max = 0.0;
  for (std::size_t component = 0U; component < 3U; ++component) {
    HUNDUN_CHECK_NEAR(report.budget_reaction_N.pressure[component],
                      expected_hybrid_row.reaction[component], 1.0e-10);
    HUNDUN_CHECK_NEAR(report.budget_reaction_N.viscous[component],
                      expected_viscous_hybrid_row[component], 1.0e-10);
    viscous_reaction_max =
        std::max(viscous_reaction_max,
                 std::abs(report.budget_reaction_N.viscous[component]));
  }
  HUNDUN_CHECK(viscous_reaction_max > 1.0e-8);

  std::vector<finite_volume::detail::ImmersedWallNormalGradient>
      analytic_wall_gradients;
  for (const auto &row :
       finite_volume::test::ImmersedOperatorTestAccess::rows(adapter)) {
    for (const auto &link : row.links) {
      analytic_wall_gradients.push_back(
          {link.id, dot(quadratic_gradient(link.wall_intercept_m),
                        link.solid_to_fluid_normal)});
    }
  }
  std::sort(analytic_wall_gradients.begin(), analytic_wall_gradients.end(),
            [](const auto &left, const auto &right) {
              return left.link < right.link;
            });
  HUNDUN_CHECK(std::adjacent_find(analytic_wall_gradients.begin(),
                                  analytic_wall_gradients.end(),
                                  [](const auto &left, const auto &right) {
                                    return left.link == right.link;
                                  }) == analytic_wall_gradients.end());
  fill_zero(momentum_residual);
  finite_volume::detail::accumulate_momentum_with_wall_normal_constraints(
      adapter, mass_flux, face_velocity, const_velocity, const_pressure,
      analytic_wall_gradients, velocity_gradient, const_viscosity,
      momentum_residual);
  const auto constrained_report = adapter.report();
  const auto expected_constrained_hybrid_row = independent_pressure_hybrid_row(
      adapter, topology, geometry, ghost_plan, const_pressure,
      &analytic_wall_gradients);
  std::vector<finite_volume::test::BoundaryRowEvaluationSnapshot>
      constrained_evaluations;
  std::array<double, 3> expected_complete_reaction{};
  double retired_hybrid_row_max_difference_after_a22 = 0.0;
  {
    const auto rows =
        finite_volume::test::ImmersedOperatorTestAccess::rows(adapter);
    const auto evaluations = finite_volume::test::ImmersedOperatorTestAccess::
        last_boundary_row_evaluations(adapter);
    std::size_t evaluation = 0U;
    for (const auto &row : rows) {
      if (row.links.empty())
        continue;
      HUNDUN_CHECK(evaluation < evaluations.size());
      HUNDUN_CHECK(boundary_evaluation_is_complete(evaluations[evaluation],
                                                   row));
      HUNDUN_CHECK(affine_evaluation_matches_quadratic_oracle(
          evaluations[evaluation], geometry, &analytic_wall_gradients));
      HUNDUN_CHECK(background_affine_evaluation_matches_quadratic_oracle(
          evaluations[evaluation], geometry, &analytic_wall_gradients));
      if (!warped) {
        const auto expected_complete_replacement =
            independent_complete_pressure_boundary_replacement(
                row, topology, geometry, domain);
        HUNDUN_CHECK(
            evaluation < expected_constrained_hybrid_row.reaction_by_row.size());
        for (std::size_t component = 0U; component < 3U; ++component)
          {
            HUNDUN_CHECK_NEAR(
                evaluations[evaluation].wall_contribution.pressure[component],
                expected_complete_replacement[component], 1.0e-12);
            expected_complete_reaction[component] -=
                expected_complete_replacement[component];
            retired_hybrid_row_max_difference_after_a22 = std::max(
                retired_hybrid_row_max_difference_after_a22,
                std::abs(-expected_complete_replacement[component] -
                         expected_constrained_hybrid_row
                             .reaction_by_row[evaluation][component]));
          }
        for (std::size_t mode = 0U; mode < 10U; ++mode) {
          const auto expected_mode =
              independent_complete_pressure_boundary_replacement_mode(
                  row, topology, geometry, domain, mode);
          std::array<double, 3> observed_mode{};
          std::array<double, 3> absolute_sum{};
          for (const auto &term : evaluations[evaluation].affine_donor_terms) {
            if (term.input_kind !=
                finite_volume::test::BoundaryAffineInputKind::pressure)
              continue;
            const auto average = uniform_cell_quadratic_basis_average(
                logical_cell(term.donor_global_cell, topology.global_extent()),
                geometry);
            const double contribution =
                term.coefficient * static_cast<double>(average[mode]);
            observed_mode[term.output_component] += contribution;
            absolute_sum[term.output_component] += std::abs(contribution);
          }
          for (const auto &term :
               evaluations[evaluation].affine_wall_gradient_terms) {
            const auto link = std::lower_bound(
                row.links.begin(), row.links.end(), term.link,
                [](const auto &candidate, std::uint64_t id) {
                  return candidate.id < id;
                });
            HUNDUN_CHECK(link != row.links.end() && link->id == term.link);
            const auto derivative = global_quadratic_normal_derivative_basis(
                link->wall_intercept_m, link->solid_to_fluid_normal);
            const double contribution =
                term.coefficient * static_cast<double>(derivative[mode]);
            observed_mode[term.output_component] += contribution;
            absolute_sum[term.output_component] += std::abs(contribution);
          }
          for (std::size_t component = 0U; component < 3U; ++component) {
            const double tolerance =
                64.0 * std::numeric_limits<double>::epsilon() *
                std::max({1.0, absolute_sum[component],
                          std::abs(expected_mode[component])});
            HUNDUN_CHECK_NEAR(observed_mode[component],
                              expected_mode[component], tolerance);
          }
        }
      }
      ++evaluation;
    }
    HUNDUN_CHECK(evaluation == evaluations.size());
    constrained_evaluations = evaluations;
  }
  if (!warped)
    for (std::size_t component = 0U; component < 3U; ++component)
      HUNDUN_CHECK_NEAR(
          constrained_report.budget_reaction_N.pressure[component],
          expected_complete_reaction[component], 1.0e-12);
  double retired_hybrid_difference_after_a22 = 0.0;
  for (std::size_t component = 0U; component < 3U; ++component)
    retired_hybrid_difference_after_a22 = std::max(
        retired_hybrid_difference_after_a22,
        std::abs(constrained_report.budget_reaction_N.pressure[component] -
                 expected_constrained_hybrid_row.reaction[component]));
  HUNDUN_CHECK((warped ? retired_hybrid_difference_after_a22
                       : retired_hybrid_row_max_difference_after_a22) >
               1.0e-8);

  auto perturbed_wall_gradients = analytic_wall_gradients;
  for (auto &condition : perturbed_wall_gradients)
    condition.value +=
        0.25 * std::sin(0.031 * static_cast<double>(condition.link + 1U));
  fill_zero(momentum_residual);
  finite_volume::detail::accumulate_momentum_with_wall_normal_constraints(
      adapter, mass_flux, face_velocity, const_velocity, const_pressure,
      perturbed_wall_gradients, velocity_gradient, const_viscosity,
      momentum_residual);
  const auto perturbed_report = adapter.report();
  const auto perturbed_evaluations =
      finite_volume::test::ImmersedOperatorTestAccess::
          last_boundary_row_evaluations(adapter);
  HUNDUN_CHECK(perturbed_evaluations.size() == constrained_evaluations.size());
  double affine_response = 0.0;
  double expected_affine_response = 0.0;
  std::array<double, 3> expected_reaction_increment{};
  for (std::size_t row = 0U; row < constrained_evaluations.size(); ++row) {
    for (std::size_t component = 0U; component < 3U; ++component) {
      double expected_row_increment = 0.0;
      for (const auto &term :
           constrained_evaluations[row].affine_wall_gradient_terms) {
        if (term.output_component != component)
          continue;
        const auto baseline = std::lower_bound(
            analytic_wall_gradients.begin(), analytic_wall_gradients.end(),
            term.link, [](const auto &candidate, std::uint64_t link) {
              return candidate.link < link;
            });
        const auto perturbed = std::lower_bound(
            perturbed_wall_gradients.begin(), perturbed_wall_gradients.end(),
            term.link, [](const auto &candidate, std::uint64_t link) {
              return candidate.link < link;
            });
        HUNDUN_CHECK(baseline != analytic_wall_gradients.end() &&
                     perturbed != perturbed_wall_gradients.end() &&
                     baseline->link == term.link && perturbed->link == term.link);
        expected_row_increment +=
            term.coefficient * (perturbed->value - baseline->value);
      }
      const double observed_row_increment =
          perturbed_evaluations[row].wall_contribution.pressure[component] -
          constrained_evaluations[row].wall_contribution.pressure[component];
      HUNDUN_CHECK_NEAR(observed_row_increment, expected_row_increment,
                        1.0e-12);
      expected_reaction_increment[component] -= expected_row_increment;
    }
  }
  for (std::size_t component = 0U; component < 3U; ++component) {
    HUNDUN_CHECK_NEAR(
        perturbed_report.budget_reaction_N.pressure[component] -
            constrained_report.budget_reaction_N.pressure[component],
        expected_reaction_increment[component], 1.0e-12);
    affine_response = std::max(
        affine_response,
        std::abs(perturbed_report.budget_reaction_N.pressure[component] -
                 constrained_report.budget_reaction_N.pressure[component]));
    expected_affine_response = std::max(
        expected_affine_response,
        std::abs(expected_reaction_increment[component]));
  }
  mpi.allreduce_fp64_in_place(&affine_response, 1U,
                              runtime::Fp64ReductionOperation::maximum);
  mpi.allreduce_fp64_in_place(&expected_affine_response, 1U,
                              runtime::Fp64ReductionOperation::maximum);
  HUNDUN_CHECK(expected_affine_response > 1.0e-8);
  HUNDUN_CHECK_NEAR(affine_response, expected_affine_response, 1.0e-10);
  int failure_rank = analytic_wall_gradients.empty() ? mpi.size() : mpi.rank();
  HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, &failure_rank, 1, MPI_INT, MPI_MIN,
                             mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(failure_rank >= 0 && failure_rank < mpi.size());
  const auto expect_collective_wall_condition_rejection =
      [&](std::vector<finite_volume::detail::ImmersedWallNormalGradient>
              conditions) {
        fill_zero(momentum_residual);
        bool rejected = false;
        std::string message;
        try {
          finite_volume::detail::
              accumulate_momentum_with_wall_normal_constraints(
                  adapter, mass_flux, face_velocity, const_velocity,
                  const_pressure, conditions, velocity_gradient,
                  const_viscosity, momentum_residual);
        } catch (const runtime::Error &error) {
          rejected = true;
          message = error.what();
        }
        HUNDUN_CHECK(rejected);
        HUNDUN_CHECK(message.find("lowest failing rank " +
                                  std::to_string(failure_rank)) !=
                     std::string::npos);
      };
  {
    auto missing = analytic_wall_gradients;
    if (mpi.rank() == failure_rank)
      missing.pop_back();
    expect_collective_wall_condition_rejection(std::move(missing));
  }
  {
    auto duplicate = analytic_wall_gradients;
    if (mpi.rank() == failure_rank)
      duplicate.push_back(duplicate.front());
    expect_collective_wall_condition_rejection(std::move(duplicate));
  }
  {
    auto nonfinite = analytic_wall_gradients;
    if (mpi.rank() == failure_rank)
      nonfinite.front().value = std::numeric_limits<double>::quiet_NaN();
    expect_collective_wall_condition_rejection(std::move(nonfinite));
  }
  fill_zero(momentum_residual);
  adapter.accumulate_momentum(mass_flux, face_velocity, const_velocity,
                              const_pressure, velocity_gradient,
                              const_viscosity, momentum_residual);

  std::array<double, 3> momentum_sum{};
  double scalar_sum = 0.0;
  for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count();
       ++cell) {
    const auto global = topology.global_cell(cell);
    const int i = global.x - box.begin.x;
    const int j = global.y - box.begin.y;
    const int k = global.z - box.begin.z;
    if (domain.region(cell) == immersed::CellRegion::solid) {
      HUNDUN_CHECK(bits(momentum_residual(i, j, k, 0)) == bits(0.0));
      HUNDUN_CHECK(bits(momentum_residual(i, j, k, 1)) == bits(0.0));
      HUNDUN_CHECK(bits(momentum_residual(i, j, k, 2)) == bits(0.0));
      HUNDUN_CHECK(bits(scalar_residual(i, j, k, 0)) == bits(0.0));
      continue;
    }
    for (int component = 0; component < 3; ++component)
      momentum_sum[static_cast<std::size_t>(component)] +=
          momentum_residual(i, j, k, component);
    scalar_sum += scalar_residual(i, j, k, 0);
  }
  HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, momentum_sum.data(), 3, MPI_DOUBLE,
                             MPI_SUM, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, &scalar_sum, 1, MPI_DOUBLE, MPI_SUM,
                             mpi.comm()) == MPI_SUCCESS);
  std::array<double, 3> reaction{};
  for (std::size_t component = 0U; component < 3U; ++component)
    reaction[component] = report.budget_reaction_N.pressure[component] +
                          report.budget_reaction_N.viscous[component];
  HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, reaction.data(), 3, MPI_DOUBLE,
                             MPI_SUM, mpi.comm()) == MPI_SUCCESS);
  const double tolerance = 1.0e-10;
  for (std::size_t component = 0U; component < 3U; ++component)
    HUNDUN_CHECK_NEAR(momentum_sum[component] + reaction[component], 0.0,
                      tolerance);
  HUNDUN_CHECK_NEAR(scalar_sum, 0.0, tolerance);

  const auto first_momentum = owned_active_bits(
      topology, domain,
      static_cast<const runtime::FieldStorage &>(storage).view<double>(
          ids.momentum_residual),
      3U, box);
  const auto first_scalar = owned_active_bits(
      topology, domain,
      static_cast<const runtime::FieldStorage &>(storage).view<double>(
          ids.scalar_residual),
      1U, box);
  const auto first_report = adapter.report();
  for (mesh::LocalCellId cell = 0U; cell < topology.local_cell_count();
       ++cell) {
    if (domain.region(cell) != immersed::CellRegion::solid)
      continue;
    const auto global = topology.global_cell(cell);
    const int i = global.x - box.begin.x;
    const int j = global.y - box.begin.y;
    const int k = global.z - box.begin.z;
    if (i < -ghost || i >= local.x + ghost || j < -ghost ||
        j >= local.y + ghost || k < -ghost || k >= local.z + ghost)
      continue;
    pressure(i, j, k, 0) = -77777.0;
    scalar(i, j, k, 0) = 88888.0;
    velocity(i, j, k, 0) = -99991.0;
    velocity(i, j, k, 1) = 99992.0;
    velocity(i, j, k, 2) = -99993.0;
  }
  fill_zero(momentum_residual);
  fill_zero(scalar_residual);
  adapter.accumulate_momentum(mass_flux, face_velocity, const_velocity,
                              const_pressure, velocity_gradient,
                              const_viscosity, momentum_residual);
  adapter.accumulate_transport(finite_volume::FiniteVolumeQuantity::scalar(0U),
                               mass_flux, face_scalar, const_scalar,
                               scalar_gradient, const_gamma, scalar_residual);
  HUNDUN_CHECK(
      first_momentum ==
      owned_active_bits(
          topology, domain,
          static_cast<const runtime::FieldStorage &>(storage).view<double>(
              ids.momentum_residual),
          3U, box));
  HUNDUN_CHECK(
      first_scalar ==
      owned_active_bits(
          topology, domain,
          static_cast<const runtime::FieldStorage &>(storage).view<double>(
              ids.scalar_residual),
          1U, box));
  const auto repeated_report = adapter.report();
  HUNDUN_CHECK(report_bits_equal(first_report, repeated_report));
  HUNDUN_CHECK(boundary_evaluations_bits_equal(
      first_row_evaluations, finite_volume::test::ImmersedOperatorTestAccess::
                                 last_boundary_row_evaluations(adapter)));

  const auto before_failure = owned_active_bits(
      topology, domain,
      static_cast<const runtime::FieldStorage &>(storage).view<double>(
          ids.momentum_residual),
      3U, box);
  const auto scalar_before_wall_flux_failure = owned_active_bits(
      topology, domain,
      static_cast<const runtime::FieldStorage &>(storage).view<double>(
          ids.scalar_residual),
      1U, box);
  std::string invalid_quantity_failure;
  try {
    adapter.accumulate_transport(
        {static_cast<finite_volume::FiniteVolumeQuantityKind>(255U), 0U},
        mass_flux, face_scalar, const_scalar, scalar_gradient, const_gamma,
        scalar_residual);
  } catch (const runtime::Error &error) {
    invalid_quantity_failure = error.what();
  }
  HUNDUN_CHECK(invalid_quantity_failure.find("quantity") != std::string::npos);
  HUNDUN_CHECK(
      scalar_before_wall_flux_failure ==
      owned_active_bits(
          topology, domain,
          static_cast<const runtime::FieldStorage &>(storage).view<double>(
              ids.scalar_residual),
          1U, box));
  std::string invalid_index_failure;
  try {
    adapter.accumulate_transport(
        {finite_volume::FiniteVolumeQuantityKind::density, 1U}, mass_flux,
        face_scalar, const_scalar, scalar_gradient, const_gamma,
        scalar_residual);
  } catch (const runtime::Error &error) {
    invalid_index_failure = error.what();
  }
  HUNDUN_CHECK(invalid_index_failure.find("quantity") != std::string::npos);
  HUNDUN_CHECK(
      scalar_before_wall_flux_failure ==
      owned_active_bits(
          topology, domain,
          static_cast<const runtime::FieldStorage &>(storage).view<double>(
              ids.scalar_residual),
          1U, box));

  for (mesh::LocalFaceId face = 0U; face < topology.local_face_count(); ++face)
    mass_flux_writer(face, 0) = -0.0;
  std::string wall_flux_failure;
  try {
    adapter.accumulate_transport(
        finite_volume::FiniteVolumeQuantity::scalar(0U), mass_flux, face_scalar,
        const_scalar, scalar_gradient, const_gamma, scalar_residual);
  } catch (const runtime::Error &error) {
    wall_flux_failure = error.what();
  }
  initialize_mass_flux();
  HUNDUN_CHECK(wall_flux_failure.find("positive zero") != std::string::npos);
  HUNDUN_CHECK(
      scalar_before_wall_flux_failure ==
      owned_active_bits(
          topology, domain,
          static_cast<const runtime::FieldStorage &>(storage).view<double>(
              ids.scalar_residual),
          1U, box));

  const int injected_failure_rank = mpi.size() > 1 ? 1 : 0;
  if (mpi.rank() == injected_failure_rank)
    for (mesh::LocalFaceId face = 0U; face < topology.local_face_count();
         ++face)
      viscosity(face, 0) = std::numeric_limits<double>::quiet_NaN();
  std::string failure_message;
  try {
    adapter.accumulate_momentum(mass_flux, face_velocity, const_velocity,
                                const_pressure, velocity_gradient,
                                const_viscosity, momentum_residual);
  } catch (const runtime::Error &error) {
    failure_message = error.what();
  }
  if (mpi.rank() == injected_failure_rank)
    for (mesh::LocalFaceId face = 0U; face < topology.local_face_count();
         ++face)
      viscosity(face, 0) = 0.01;
  int every_rank_failed = failure_message.empty() ? 0 : 1;
  HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, &every_rank_failed, 1, MPI_INT,
                             MPI_MIN, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(every_rank_failed == 1);
  HUNDUN_CHECK(failure_message.find("viscosity") != std::string::npos);
  HUNDUN_CHECK(failure_message.find("lowest failing rank " +
                                    std::to_string(injected_failure_rank)) !=
               std::string::npos);
  HUNDUN_CHECK(
      before_failure ==
      owned_active_bits(
          topology, domain,
          static_cast<const runtime::FieldStorage &>(storage).view<double>(
              ids.momentum_residual),
          3U, box));
  const auto after_failure_report = adapter.report();
  HUNDUN_CHECK(report_bits_equal(repeated_report, after_failure_report));
}

} // namespace

int main(int argc, char **argv) {
  runtime::MpiEnvironment environment(argc, argv);
  auto mpi = runtime::MpiContext::duplicate(MPI_COMM_WORLD);
  return test::run([&] {
    prove_periodic_checkerboard_operator_selection();
    test_report_comparison_oracle();
    test_physical_background_term_oracle();
    test_replacement_term_partition_oracle();
    const std::string mode = argc > 1 ? argv[1] : "uniform";
    if (mode == "uniform")
      run_case(mpi, false);
    else if (mode == "warped")
      run_case(mpi, true);
    else
      throw runtime::Error("unknown immersed operator test mode");
  });
}
