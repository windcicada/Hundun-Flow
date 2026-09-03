// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/ib_quadratic_reconstruction.hpp"
#include "hundun/mesh_geometry.hpp"
#include "hundun/mesh_topology.hpp"
#include "hundun/rt_field_descriptor.hpp"
#include "hundun/rt_field_registry.hpp"
#include "hundun/rt_field_storage.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "hundun/rt_structured_decomposition.hpp"
#include "src/ib_deterministic_qr_detail.hpp"
#include "src/ib_quadratic_reconstruction_detail.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using hundun::immersed::QuadraticReconstruction;
using hundun::immersed::detail::geometry_validation_is_scoped;
using hundun::immersed::detail::QuadraticReconstructionWeights;
using hundun::immersed::detail::ValidatedGeometryScope;
using hundun::mesh::AnalyticWarpedBoxMapping;
using hundun::mesh::MeshGeometry;
using hundun::mesh::MeshTopology;
using hundun::mesh::UniformBoxMapping;
using hundun::runtime::DecompositionOptions;
using hundun::runtime::FieldDescriptor;
using hundun::runtime::FieldRegistry;
using hundun::runtime::FieldStorage;
using hundun::runtime::FunctionSpace;
using hundun::runtime::Int3;
using hundun::runtime::MpiContext;
using hundun::runtime::MpiEnvironment;
using hundun::runtime::OutputPolicy;
using hundun::runtime::Real3;
using hundun::runtime::RestartPolicy;
using hundun::runtime::ScalarType;
using hundun::runtime::StructuredDecomposition;

constexpr Int3 kExtent{12, 10, 8};
constexpr Real3 kOrigin{-1.0, -0.5, -0.25};
constexpr Real3 kLength{2.0, 1.5, 1.0};
constexpr Real3 kWall{0.2, 0.325, 0.3125};
constexpr double kScale = 0.15;

std::vector<Int3> donors() {
  return {{7, 5, 4}, {8, 5, 4}, {9, 4, 3}, {7, 6, 4}, {7, 4, 4},
          {7, 5, 5}, {7, 5, 3}, {8, 6, 4}, {8, 4, 4}, {8, 5, 5},
          {8, 5, 3}, {9, 6, 5}, {9, 4, 5}, {9, 6, 3}};
}

std::vector<Int3> thirty_two_donors() {
  std::vector<Int3> result;
  result.reserve(32U);
  for (int x = 7; x <= 10; ++x) {
    for (int y = 4; y <= 6; ++y) {
      for (int z = 3; z <= 5; ++z) {
        if (result.size() == 32U) {
          break;
        }
        result.push_back({x, y, z});
      }
    }
  }
  HUNDUN_CHECK(result.size() == 32U);
  return result;
}

Real3 add(Real3 lhs, Real3 rhs) {
  return {lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z};
}

Real3 subtract(Real3 lhs, Real3 rhs) {
  return {lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z};
}

Real3 multiply(Real3 value, double factor) {
  return {value.x * factor, value.y * factor, value.z * factor};
}

double dot(Real3 lhs, Real3 rhs) {
  return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

Real3 cross(Real3 lhs, Real3 rhs) {
  return {lhs.y * rhs.z - lhs.z * rhs.y, lhs.z * rhs.x - lhs.x * rhs.z,
          lhs.x * rhs.y - lhs.y * rhs.x};
}

std::array<double, 10> basis(Real3 point) {
  const auto offset = subtract(point, kWall);
  const double n = offset.x / kScale;
  const double t1 = offset.y / kScale;
  const double t2 = offset.z / kScale;
  return {1.0, n, t1, t2, n * n, n * t1, n * t2, t1 * t1, t1 * t2, t2 * t2};
}

constexpr std::array<double, 10> kCoefficients{1.0,  2.0, -3.0, 0.5,  0.7,
                                               -0.2, 0.4, 0.3,  -0.6, 0.9};

double polynomial(Real3 point) {
  const auto values = basis(point);
  double result = 0.0;
  for (std::size_t index = 0; index < values.size(); ++index) {
    result += kCoefficients[index] * values[index];
  }
  return result;
}

double independent_cell_average(Int3 cell, const MeshGeometry &geometry) {
  const std::array<Int3, 8> offsets{Int3{0, 0, 0}, Int3{1, 0, 0}, Int3{1, 1, 0},
                                    Int3{0, 1, 0}, Int3{0, 0, 1}, Int3{1, 0, 1},
                                    Int3{1, 1, 1}, Int3{0, 1, 1}};
  std::array<Real3, 8> vertices{};
  Real3 reference{};
  for (std::size_t index = 0; index < vertices.size(); ++index) {
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
    const std::array<Real3, 4> tetra{reference, vertices[triangle[0]],
                                     vertices[triangle[1]],
                                     vertices[triangle[2]]};
    const double tetra_volume =
        dot(subtract(tetra[1], tetra[0]),
            cross(subtract(tetra[2], tetra[0]), subtract(tetra[3], tetra[0]))) /
        6.0;
    HUNDUN_CHECK(tetra_volume > 0.0);
    std::array<double, 10> average{};
    average[0] = 1.0;
    std::array<std::array<double, 3>, 4> coordinates{};
    for (std::size_t vertex = 0; vertex < tetra.size(); ++vertex) {
      const auto values = basis(tetra[vertex]);
      coordinates[vertex] = {values[1], values[2], values[3]};
      average[1] += coordinates[vertex][0] / 4.0;
      average[2] += coordinates[vertex][1] / 4.0;
      average[3] += coordinates[vertex][2] / 4.0;
    }
    const auto product_average = [&](std::size_t lhs, std::size_t rhs) {
      double sum_lhs = 0.0;
      double sum_rhs = 0.0;
      double diagonal = 0.0;
      for (const auto coordinate : coordinates) {
        sum_lhs += coordinate[lhs];
        sum_rhs += coordinate[rhs];
        diagonal += coordinate[lhs] * coordinate[rhs];
      }
      return (sum_lhs * sum_rhs + diagonal) / 20.0;
    };
    average[4] = product_average(0, 0);
    average[5] = product_average(0, 1);
    average[6] = product_average(0, 2);
    average[7] = product_average(1, 1);
    average[8] = product_average(1, 2);
    average[9] = product_average(2, 2);
    double tetra_average = 0.0;
    for (std::size_t term = 0; term < average.size(); ++term) {
      tetra_average += kCoefficients[term] * average[term];
    }
    integral += tetra_volume * tetra_average;
    volume += tetra_volume;
  }
  HUNDUN_CHECK(volume > 0.0);
  return integral / volume;
}

template <class Function> std::string expect_error(Function &&function) {
  try {
    function();
  } catch (const hundun::runtime::Error &error) {
    return error.what();
  }
  throw std::runtime_error("expected hundun::runtime::Error");
}

bool contains(Int3 begin, Int3 end, Int3 point) {
  return point.x >= begin.x && point.x < end.x && point.y >= begin.y &&
         point.y < end.y && point.z >= begin.z && point.z < end.z;
}

void run_success_case(const MpiContext &context, bool warped,
                      std::optional<Int3> process_grid) {
  DecompositionOptions options;
  options.process_grid = process_grid;
  auto decomposition = StructuredDecomposition::create(
      context, kExtent, {false, false, false}, options);
  MeshTopology topology(decomposition);
  std::optional<MeshGeometry> geometry;
  if (warped) {
    geometry.emplace(
        topology,
        AnalyticWarpedBoxMapping(kOrigin, kLength, Real3{0.02, -0.015, 0.01}));
  } else {
    geometry.emplace(topology, UniformBoxMapping(kOrigin, kLength));
  }

  const auto reconstruction = QuadraticReconstruction::create(
      kWall, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, kScale,
      {7, 5, 4}, donors(), topology, *geometry);
  HUNDUN_CHECK(reconstruction.quality().rank == 10U);
  HUNDUN_CHECK(reconstruction.quality().halo_reach <= 4U);
  auto reordered_donors = donors();
  std::reverse(reordered_donors.begin(), reordered_donors.end());
  const auto reordered = QuadraticReconstruction::create(
      kWall, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, kScale,
      {7, 5, 4}, reordered_donors, topology, *geometry);
  HUNDUN_CHECK(reordered.quality().rank == reconstruction.quality().rank);
  HUNDUN_CHECK(reordered.quality().halo_reach ==
               reconstruction.quality().halo_reach);
  HUNDUN_CHECK(reordered.quality().pivot_fingerprint ==
               reconstruction.quality().pivot_fingerprint);
  HUNDUN_CHECK(std::memcmp(&reordered.quality().condition_estimate,
                           &reconstruction.quality().condition_estimate,
                           sizeof(double)) == 0);
  const auto maximum = QuadraticReconstruction::create(
      kWall, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, kScale,
      {7, 5, 4}, thirty_two_donors(), topology, *geometry);
  HUNDUN_CHECK(maximum.quality().rank == 10U);
  auto reach_four_donors = donors();
  reach_four_donors.back() = {11, 6, 3};
  const auto reach_four = QuadraticReconstruction::create(
      kWall, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, kScale,
      {7, 5, 4}, reach_four_donors, topology, *geometry);
  HUNDUN_CHECK(reach_four.quality().halo_reach == 4U);
  const auto &global_ids =
      QuadraticReconstructionWeights::donor_global_ids(reconstruction);
  HUNDUN_CHECK(global_ids.size() == donors().size());
  HUNDUN_CHECK(std::is_sorted(global_ids.begin(), global_ids.end()));
  HUNDUN_CHECK(global_ids ==
               QuadraticReconstructionWeights::donor_global_ids(reordered));
  unsigned long long fingerprint = static_cast<unsigned long long>(
      reconstruction.quality().pivot_fingerprint);
  unsigned long long fingerprint_min = 0U;
  unsigned long long fingerprint_max = 0U;
  HUNDUN_CHECK(MPI_Allreduce(&fingerprint, &fingerprint_min, 1,
                             MPI_UNSIGNED_LONG_LONG, MPI_MIN,
                             context.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(MPI_Allreduce(&fingerprint, &fingerprint_max, 1,
                             MPI_UNSIGNED_LONG_LONG, MPI_MAX,
                             context.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(fingerprint_min == fingerprint_max);
  unsigned long long condition_bits = 0U;
  static_assert(sizeof(condition_bits) == sizeof(double));
  std::memcpy(&condition_bits, &reconstruction.quality().condition_estimate,
              sizeof(condition_bits));
  unsigned long long condition_bits_min = 0U;
  unsigned long long condition_bits_max = 0U;
  HUNDUN_CHECK(MPI_Allreduce(&condition_bits, &condition_bits_min, 1,
                             MPI_UNSIGNED_LONG_LONG, MPI_MIN,
                             context.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(MPI_Allreduce(&condition_bits, &condition_bits_max, 1,
                             MPI_UNSIGNED_LONG_LONG, MPI_MAX,
                             context.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(condition_bits_min == condition_bits_max);

  FieldRegistry registry;
  const auto id = registry.declare_field(FieldDescriptor{
      "q", "1", "task4", FunctionSpace::cell_average, ScalarType::float64, 1U,
      4, false, RestartPolicy::transient, OutputPolicy::never});
  registry.freeze();
  FieldStorage storage(registry, decomposition.local_extent());
  auto view = storage.view<double>(id);
  const auto box = decomposition.owned_box();
  int local_point_sample_mutation = 0;
  for (int k = -4; k < decomposition.local_extent().z + 4; ++k) {
    for (int j = -4; j < decomposition.local_extent().y + 4; ++j) {
      for (int i = -4; i < decomposition.local_extent().x + 4; ++i) {
        const Int3 global{box.begin.x + i, box.begin.y + j, box.begin.z + k};
        if (global.x < 0 || global.y < 0 || global.z < 0 ||
            global.x >= kExtent.x || global.y >= kExtent.y ||
            global.z >= kExtent.z) {
          view(i, j, k, 0) = 0.0;
          continue;
        }
        const double average = independent_cell_average(global, *geometry);
        view(i, j, k, 0) = average;
        if (const auto local =
                topology.find_local_cell(topology.global_cell_id(global));
            local.has_value() &&
            std::abs(average - polynomial(geometry->cell_center_m(*local))) >
                1.0e-8) {
          local_point_sample_mutation = 1;
        }
      }
    }
  }
  int global_point_sample_mutation = 0;
  HUNDUN_CHECK(MPI_Allreduce(&local_point_sample_mutation,
                             &global_point_sample_mutation, 1, MPI_INT, MPI_MAX,
                             context.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(global_point_sample_mutation == 1);
  if (contains(box.begin, box.end, {7, 5, 4})) {
    const auto read =
        static_cast<const FieldStorage &>(storage).view<double>(id);
    const Real3 point{0.31, 0.27, 0.36};
    const double expected = polynomial(point);
    const double value = reconstruction.value(point, read, 0U);
    const double reordered_value = reordered.value(point, read, 0U);
    const double tolerance = 4096.0 * std::numeric_limits<double>::epsilon() *
                             std::max(1.0, std::abs(expected));
    HUNDUN_CHECK_NEAR(value, expected, tolerance);
    HUNDUN_CHECK(std::memcmp(&value, &reordered_value, sizeof(double)) == 0);
    const auto value_weights =
        QuadraticReconstructionWeights::value_weights(reconstruction, point);
    const auto reordered_weights =
        QuadraticReconstructionWeights::value_weights(reordered, point);
    HUNDUN_CHECK(value_weights.size() == global_ids.size());
    HUNDUN_CHECK(reordered_weights.size() == value_weights.size());
    double weighted_value = 0.0;
    bool has_remote_donor = false;
    for (std::size_t donor = 0; donor < value_weights.size(); ++donor) {
      HUNDUN_CHECK(value_weights[donor].global_cell == global_ids[donor]);
      HUNDUN_CHECK(value_weights[donor].global_cell ==
                   reordered_weights[donor].global_cell);
      HUNDUN_CHECK(std::memcmp(&value_weights[donor].weight,
                               &reordered_weights[donor].weight,
                               sizeof(double)) == 0);
      const auto local_cell =
          topology.find_local_cell(value_weights[donor].global_cell);
      HUNDUN_CHECK(local_cell.has_value());
      has_remote_donor =
          has_remote_donor ||
          !contains(box.begin, box.end, topology.global_cell(*local_cell));
      weighted_value += value_weights[donor].weight *
                        independent_cell_average(
                            topology.global_cell(*local_cell), *geometry);
    }
    if (context.size() == 4) {
      HUNDUN_CHECK(has_remote_donor);
    }
    HUNDUN_CHECK_NEAR(weighted_value, expected, tolerance);
    const auto values = basis(point);
    const double n = values[1];
    const double t1 = values[2];
    const double t2 = values[3];
    const Real3 expected_gradient{
        (2.0 + 1.4 * n - 0.2 * t1 + 0.4 * t2) / kScale,
        (-3.0 - 0.2 * n + 0.6 * t1 - 0.6 * t2) / kScale,
        (0.5 + 0.4 * n - 0.6 * t1 + 1.8 * t2) / kScale};
    const Real3 gradient = reconstruction.gradient(point, read, 0U);
    HUNDUN_CHECK_NEAR(gradient.x, expected_gradient.x,
                      4.0 * tolerance / kScale);
    HUNDUN_CHECK_NEAR(gradient.y, expected_gradient.y,
                      4.0 * tolerance / kScale);
    HUNDUN_CHECK_NEAR(gradient.z, expected_gradient.z,
                      4.0 * tolerance / kScale);
    const double constrained_value =
        hundun::immersed::detail::value_with_origin_constraint(
            reconstruction, point, read, 0U, kCoefficients[0]);
    const Real3 constrained_gradient =
        hundun::immersed::detail::gradient_with_origin_constraint(
            reconstruction, point, read, 0U, kCoefficients[0]);
    HUNDUN_CHECK_NEAR(constrained_value, expected, tolerance);
    HUNDUN_CHECK_NEAR(constrained_gradient.x, expected_gradient.x,
                      4.0 * tolerance / kScale);
    HUNDUN_CHECK_NEAR(constrained_gradient.y, expected_gradient.y,
                      4.0 * tolerance / kScale);
    HUNDUN_CHECK_NEAR(constrained_gradient.z, expected_gradient.z,
                      4.0 * tolerance / kScale);
    const double exact_wall_normal_gradient = kCoefficients[1] / kScale;
    const double normal_gradient_constrained_value =
        hundun::immersed::detail::value_with_origin_normal_gradient(
            reconstruction, point, read, 0U, exact_wall_normal_gradient);
    const Real3 normal_gradient_constrained_gradient =
        hundun::immersed::detail::gradient_with_origin_normal_gradient(
            reconstruction, point, read, 0U, exact_wall_normal_gradient);
    HUNDUN_CHECK_NEAR(normal_gradient_constrained_value, expected, tolerance);
    HUNDUN_CHECK_NEAR(normal_gradient_constrained_gradient.x,
                      expected_gradient.x, 4.0 * tolerance / kScale);
    HUNDUN_CHECK_NEAR(normal_gradient_constrained_gradient.y,
                      expected_gradient.y, 4.0 * tolerance / kScale);
    HUNDUN_CHECK_NEAR(normal_gradient_constrained_gradient.z,
                      expected_gradient.z, 4.0 * tolerance / kScale);
    const auto zero_normal_gradient =
        hundun::immersed::detail::gradient_with_origin_normal_gradient(
            reconstruction, kWall, read, 0U, 0.0);
    HUNDUN_CHECK_NEAR(zero_normal_gradient.x, 0.0, 4.0 * tolerance / kScale);
    HUNDUN_CHECK(std::abs(zero_normal_gradient.x - kCoefficients[1] / kScale) >
                 1.0e-6);
    HUNDUN_CHECK_NEAR(hundun::immersed::detail::value_with_origin_constraint(
                          reconstruction, kWall, read, 0U, 0.0),
                      0.0, tolerance);
    const auto zero_wall_gradient =
        hundun::immersed::detail::gradient_with_origin_constraint(
            reconstruction, kWall, read, 0U, 0.0);
    const auto unconstrained_wall_gradient =
        reconstruction.gradient(kWall, read, 0U);
    HUNDUN_CHECK(
        std::abs(zero_wall_gradient.x - unconstrained_wall_gradient.x) +
            std::abs(zero_wall_gradient.y - unconstrained_wall_gradient.y) +
            std::abs(zero_wall_gradient.z - unconstrained_wall_gradient.z) >
        1.0e-6);
    for (const auto direction : std::array<Real3, 3>{
             {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}}) {
      const auto weights =
          QuadraticReconstructionWeights::directional_gradient_weights(
              reconstruction, point, direction);
      double observed = 0.0;
      for (const auto weight : weights) {
        const auto local_cell = topology.find_local_cell(weight.global_cell);
        HUNDUN_CHECK(local_cell.has_value());
        observed += weight.weight *
                    independent_cell_average(topology.global_cell(*local_cell),
                                             *geometry);
      }
      const double expected_directional = direction.x * expected_gradient.x +
                                          direction.y * expected_gradient.y +
                                          direction.z * expected_gradient.z;
      HUNDUN_CHECK_NEAR(observed, expected_directional,
                        4.0 * tolerance / kScale);
    }

    HUNDUN_CHECK_NEAR(reconstruction.value(kWall, read, 0U), kCoefficients[0],
                      tolerance);
    const Real3 wall_gradient = reconstruction.gradient(kWall, read, 0U);
    HUNDUN_CHECK_NEAR(wall_gradient.x, kCoefficients[1] / kScale,
                      4.0 * tolerance / kScale);
    HUNDUN_CHECK_NEAR(wall_gradient.y, kCoefficients[2] / kScale,
                      4.0 * tolerance / kScale);
    HUNDUN_CHECK_NEAR(wall_gradient.z, kCoefficients[3] / kScale,
                      4.0 * tolerance / kScale);

    const auto ghost_local =
        topology.find_local_cell(topology.global_cell_id({6, 5, 4}));
    HUNDUN_CHECK(ghost_local.has_value());
    const Real3 ghost_center = geometry->cell_center_m(*ghost_local);
    const double ghost_expected = polynomial(ghost_center);
    HUNDUN_CHECK_NEAR(reconstruction.value(ghost_center, read, 0U),
                      ghost_expected,
                      tolerance * std::max(1.0, std::abs(ghost_expected)));

    HUNDUN_CHECK(expect_error([&] {
                   reconstruction.value(point, read, 1U);
                 }).find("component") != std::string::npos);

    FieldRegistry narrow_registry;
    const auto narrow_id = narrow_registry.declare_field(
        FieldDescriptor{"narrow", "1", "task4", FunctionSpace::cell_average,
                        ScalarType::float64, 1U, 1, false,
                        RestartPolicy::transient, OutputPolicy::never});
    narrow_registry.freeze();
    FieldStorage narrow_storage(narrow_registry, decomposition.local_extent());
    auto narrow_read = static_cast<const FieldStorage &>(narrow_storage)
                           .view<double>(narrow_id);
    if (decomposition.process_grid().x == 4) {
      HUNDUN_CHECK(expect_error([&] {
                     maximum.value(point, narrow_read, 0U);
                   }).find("spatial index") != std::string::npos);
    }

    storage.begin_rebuild();
    HUNDUN_CHECK(expect_error([&] {
                   reconstruction.value(point, read, 0U);
                 }).find("stale") != std::string::npos);
  }
  context.barrier();
}

void run_failures(const MpiContext &context, std::optional<Int3> process_grid) {
  DecompositionOptions options;
  options.process_grid = process_grid;
  auto decomposition = StructuredDecomposition::create(
      context, kExtent, {false, false, false}, options);
  MeshTopology topology(decomposition);
  MeshGeometry geometry(topology, UniformBoxMapping(kOrigin, kLength));
  const auto make = [&](const std::vector<Int3> &selected, Int3 anchor) {
    return QuadraticReconstruction::create(
        kWall, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, kScale,
        anchor, selected, topology, geometry);
  };
  auto too_few = donors();
  too_few.pop_back();
  HUNDUN_CHECK(expect_error([&] {
                 make(too_few, {7, 5, 4});
               }).find("fourteen") != std::string::npos);
  auto too_many = thirty_two_donors();
  too_many.push_back({11, 5, 4});
  HUNDUN_CHECK(expect_error([&] {
                 make(too_many, {7, 5, 4});
               }).find("thirty-two") != std::string::npos);
  auto duplicate = donors();
  duplicate.back() = duplicate.front();
  HUNDUN_CHECK(expect_error([&] {
                 make(duplicate, {7, 5, 4});
               }).find("unique") != std::string::npos);
  auto reach = donors();
  reach.back() = {11, 6, 3};
  HUNDUN_CHECK(expect_error([&] {
                 make(reach, {6, 5, 4});
               }).find("reach") != std::string::npos);
  auto nonpositive = donors();
  nonpositive.front() = {6, 5, 4};
  HUNDUN_CHECK(expect_error([&] {
                 make(nonpositive, {7, 5, 4});
               }).find("positive-normal") != std::string::npos);
  std::vector<Int3> two_bands;
  for (int x : std::array<int, 2>{7, 8}) {
    for (int y : std::array<int, 3>{4, 5, 6}) {
      for (int z : std::array<int, 3>{3, 4, 5}) {
        two_bands.push_back({x, y, z});
      }
    }
  }
  two_bands.resize(14U);
  HUNDUN_CHECK(expect_error([&] {
                 make(two_bands, {7, 5, 4});
               }).find("three positive-normal") != std::string::npos);
  const std::vector<Int3> no_quadrant{
      {7, 5, 4}, {7, 5, 5}, {7, 6, 4}, {7, 6, 5}, {8, 5, 4},
      {8, 5, 5}, {8, 6, 4}, {8, 6, 5}, {9, 5, 4}, {9, 5, 5},
      {9, 6, 4}, {9, 6, 5}, {7, 7, 4}, {8, 7, 5}};
  HUNDUN_CHECK(expect_error([&] {
                 make(no_quadrant, {7, 5, 4});
               }).find("quadrant") != std::string::npos);
  HUNDUN_CHECK(expect_error([&] {
                 QuadraticReconstruction::create(
                     kWall, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0},
                     0.0, {7, 5, 4}, donors(), topology, geometry);
               }).find("positive scale") != std::string::npos);
  HUNDUN_CHECK(expect_error([&] {
                 QuadraticReconstruction::create(
                     kWall, {1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 0.0, 1.0},
                     kScale, {7, 5, 4}, donors(), topology, geometry);
               }).find("orthonormal") != std::string::npos);
  auto incompatible_decomposition = StructuredDecomposition::create(
      context, {16, 10, 8}, {false, false, false}, options);
  MeshTopology incompatible_topology(incompatible_decomposition);
  MeshGeometry incompatible_geometry(incompatible_topology,
                                     UniformBoxMapping(kOrigin, kLength));
  HUNDUN_CHECK(!geometry_validation_is_scoped(topology, geometry));
  {
    ValidatedGeometryScope validated(topology, geometry);
    HUNDUN_CHECK(geometry_validation_is_scoped(topology, geometry));
    HUNDUN_CHECK(
        !geometry_validation_is_scoped(topology, incompatible_geometry));
    HUNDUN_CHECK(expect_error([&] {
                   ValidatedGeometryScope nested(topology, geometry);
                 }).find("nested") != std::string::npos);
    HUNDUN_CHECK(geometry_validation_is_scoped(topology, geometry));
  }
  HUNDUN_CHECK(!geometry_validation_is_scoped(topology, geometry));
  HUNDUN_CHECK(expect_error([&] {
                 QuadraticReconstruction::create(
                     kWall, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0},
                     kScale, {7, 5, 4}, donors(), topology,
                     incompatible_geometry);
               }).find("incompatible") != std::string::npos);
  context.barrier();
}

void run_tests(const MpiContext &context, const std::string &mode,
               const std::string &selector) {
  std::optional<Int3> grid;
  if (selector == "2x2x1") {
    grid = Int3{2, 2, 1};
  } else if (selector == "4x1x1") {
    grid = Int3{4, 1, 1};
  }
  if (mode == "failures") {
    run_failures(context, grid);
  } else {
    run_success_case(context, false, grid);
    run_success_case(context, true, grid);
  }
}

[[noreturn]] void abort_active_mpi_test(const MpiContext &context) noexcept {
  (void)MPI_Abort(context.comm(), EXIT_FAILURE);
  std::abort();
}

} // namespace

int main(int argc, char **argv) {
  int result = EXIT_FAILURE;
  {
    MpiEnvironment environment(argc, argv);
    MpiContext context = MpiContext::duplicate(MPI_COMM_WORLD);
    const std::string mode = argc > 1 ? argv[1] : "success";
    const std::string selector = argc > 2 ? argv[2] : "default";
    result = hundun::test::run([&] { run_tests(context, mode, selector); });
    if (result != EXIT_SUCCESS) {
      abort_active_mpi_test(context);
    }
  }
  return result;
}
