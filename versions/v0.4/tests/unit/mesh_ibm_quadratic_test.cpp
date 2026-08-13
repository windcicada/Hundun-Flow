// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_ibm.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <string_view>
#include <utility>
#include <vector>

namespace {

using namespace hundun::v04;

constexpr double kInvSqrt2 = 0.707106781186547524400844362104849039;
constexpr double kInvSqrt3 = 0.577350269189625764509148780501957456;
constexpr double kInvSqrt6 = 0.408248290463863016366214012450981899;

bool expect(bool condition, std::string_view description) {
  if (!condition) {
    std::cerr << "FAIL: " << description << '\n';
  }
  return condition;
}

bool close(double actual, double expected, double tolerance = 2.0e-11) {
  return std::isfinite(actual) &&
         std::abs(actual - expected) <=
             tolerance * std::max({1.0, std::abs(actual), std::abs(expected)});
}

double dot(Real3 left, Real3 right) {
  return left.x * right.x + left.y * right.y + left.z * right.z;
}

Real3 add(Real3 left, Real3 right) {
  return {left.x + right.x, left.y + right.y, left.z + right.z};
}

Real3 multiply(double scale, Real3 value) {
  return {scale * value.x, scale * value.y, scale * value.z};
}

struct Polynomial {
  // [1,n,t1,t2,n^2,n*t1,n*t2,t1^2,t1*t2,t2^2]
  std::array<double, 10U> coefficient{
      2.0, -0.4, 0.7, -1.2, 0.13,
      -0.21, 0.09, 0.31, -0.17, 0.24};
};

std::array<double, 10U> point_basis(double n, double t1, double t2) {
  return {1.0, n, t1, t2, n * n, n * t1, n * t2,
          t1 * t1, t1 * t2, t2 * t2};
}

double polynomial_value(const Polynomial& polynomial, double n, double t1,
                        double t2) {
  const auto basis = point_basis(n, t1, t2);
  double value = 0.0;
  for (std::size_t column = 0U; column < basis.size(); ++column) {
    value += polynomial.coefficient[column] * basis[column];
  }
  return value;
}

double polynomial_directional_derivative(const Polynomial& polynomial,
                                         double n, double t1, double t2,
                                         double dn, double dt1,
                                         double dt2) {
  const std::array<double, 10U> derivative{
      0.0,
      dn,
      dt1,
      dt2,
      2.0 * n * dn,
      n * dt1 + t1 * dn,
      n * dt2 + t2 * dn,
      2.0 * t1 * dt1,
      t1 * dt2 + t2 * dt1,
      2.0 * t2 * dt2};
  double value = 0.0;
  for (std::size_t column = 0U; column < derivative.size(); ++column) {
    value += polynomial.coefficient[column] * derivative[column];
  }
  return value;
}

std::array<double, 10U> cell_average_basis(const QuadraticFrame& frame,
                                           const QuadraticDonorCell& donor) {
  // Independent moment oracle: split the axis-aligned hexahedron into the
  // twelve equal-volume tetrahedra formed by its centre and two triangles on
  // each face.  Exact simplex first/second moments are used, rather than the
  // Cartesian covariance identity used by the product implementation.
  std::array<Real3, 8U> corner{};
  for (std::size_t index = 0U; index < corner.size(); ++index) {
    corner[index] = {
        donor.centre.x + ((index & 1U) == 0U ? -0.5 : 0.5) * donor.widths.x,
        donor.centre.y + ((index & 2U) == 0U ? -0.5 : 0.5) * donor.widths.y,
        donor.centre.z + ((index & 4U) == 0U ? -0.5 : 0.5) * donor.widths.z};
  }
  constexpr std::array<std::array<std::size_t, 4U>, 6U> faces{{
      {{0U, 4U, 6U, 2U}}, {{1U, 3U, 7U, 5U}},
      {{0U, 1U, 5U, 4U}}, {{2U, 6U, 7U, 3U}},
      {{0U, 2U, 3U, 1U}}, {{4U, 5U, 7U, 6U}},
  }};
  const auto local = [&](Real3 point) {
    const Real3 offset{point.x - frame.origin.x, point.y - frame.origin.y,
                       point.z - frame.origin.z};
    return std::array<double, 3U>{dot(offset, frame.normal) / frame.scale,
                                  dot(offset, frame.tangent1) / frame.scale,
                                  dot(offset, frame.tangent2) / frame.scale};
  };
  const auto tetra_average_product = [](const std::array<double, 4U>& left,
                                        const std::array<double, 4U>& right) {
    double diagonal = 0.0;
    double sum_left = 0.0;
    double sum_right = 0.0;
    for (std::size_t vertex = 0U; vertex < 4U; ++vertex) {
      diagonal += left[vertex] * right[vertex];
      sum_left += left[vertex];
      sum_right += right[vertex];
    }
    return (diagonal + sum_left * sum_right) / 20.0;
  };
  std::array<double, 10U> average{};
  for (const auto& face : faces) {
    for (std::size_t half = 0U; half < 2U; ++half) {
      const std::array<Real3, 4U> vertices{
          donor.centre, corner[face[0U]], corner[face[half + 1U]],
          corner[face[half + 2U]]};
      std::array<double, 4U> n{};
      std::array<double, 4U> t1{};
      std::array<double, 4U> t2{};
      for (std::size_t vertex = 0U; vertex < vertices.size(); ++vertex) {
        const auto value = local(vertices[vertex]);
        n[vertex] = value[0U];
        t1[vertex] = value[1U];
        t2[vertex] = value[2U];
      }
      average[0U] += 1.0 / 12.0;
      average[1U] +=
          (n[0U] + n[1U] + n[2U] + n[3U]) / 48.0;
      average[2U] +=
          (t1[0U] + t1[1U] + t1[2U] + t1[3U]) / 48.0;
      average[3U] +=
          (t2[0U] + t2[1U] + t2[2U] + t2[3U]) / 48.0;
      average[4U] += tetra_average_product(n, n) / 12.0;
      average[5U] += tetra_average_product(n, t1) / 12.0;
      average[6U] += tetra_average_product(n, t2) / 12.0;
      average[7U] += tetra_average_product(t1, t1) / 12.0;
      average[8U] += tetra_average_product(t1, t2) / 12.0;
      average[9U] += tetra_average_product(t2, t2) / 12.0;
    }
  }
  return average;
}

double donor_average(const Polynomial& polynomial,
                     const QuadraticFrame& frame,
                     const QuadraticDonorCell& donor) {
  const auto basis = cell_average_basis(frame, donor);
  double value = 0.0;
  for (std::size_t column = 0U; column < basis.size(); ++column) {
    value += polynomial.coefficient[column] * basis[column];
  }
  return value;
}

QuadraticFrame rotated_frame() {
  QuadraticFrame frame;
  frame.origin = {1.25, -0.75, 2.5};
  frame.normal = {kInvSqrt3, kInvSqrt3, kInvSqrt3};
  frame.tangent1 = {kInvSqrt2, -kInvSqrt2, 0.0};
  frame.tangent2 = {kInvSqrt6, kInvSqrt6, -2.0 * kInvSqrt6};
  frame.scale = 0.75;
  frame.anchor_global_cell = {10, 10, 10};
  return frame;
}

std::vector<QuadraticDonorCell> complete_donors(
    const QuadraticFrame& frame) {
  constexpr std::array<std::array<double, 3U>, 20U> local{{
      {{0.45, 0.55, 0.65}}, {{0.55, -0.75, 0.85}},
      {{0.65, -0.95, -0.55}}, {{0.75, 0.85, -0.75}},
      {{1.15, 1.25, 0.45}}, {{1.25, -0.45, 1.35}},
      {{1.35, -1.15, -0.85}}, {{1.45, 0.65, -1.25}},
      {{2.05, 0.35, 1.55}}, {{2.15, -1.45, 0.35}},
      {{2.25, -0.65, -1.45}}, {{2.35, 1.55, -0.45}},
      {{3.05, 1.15, 1.25}}, {{3.15, -0.35, 0.65}},
      {{3.25, -1.55, -0.35}}, {{3.35, 0.45, -1.55}},
      {{0.95, 1.65, 0.25}}, {{1.75, -1.65, 0.55}},
      {{2.75, -0.25, -1.65}}, {{3.55, 1.35, -0.95}},
  }};
  constexpr std::array<Int3, 20U> offset{{
      {1, 1, 1}, {1, -1, 1}, {1, -1, -1}, {1, 1, -1},
      {2, 2, 1}, {2, -1, 2}, {2, -2, -1}, {2, 1, -2},
      {3, 1, 2}, {3, -2, 1}, {3, -1, -2}, {3, 2, -1},
      {4, 1, 1}, {4, -1, 1}, {4, -2, -1}, {4, 1, -2},
      {1, 2, 0}, {2, -2, 1}, {3, 0, -2}, {4, 2, -1},
  }};
  std::vector<QuadraticDonorCell> donors;
  donors.reserve(local.size());
  for (std::size_t index = 0U; index < local.size(); ++index) {
    const double n = local[index][0U];
    const double t1 = local[index][1U];
    const double t2 = local[index][2U];
    Real3 centre = frame.origin;
    centre = add(centre, multiply(frame.scale * n, frame.normal));
    centre = add(centre, multiply(frame.scale * t1, frame.tangent1));
    centre = add(centre, multiply(frame.scale * t2, frame.tangent2));
    QuadraticDonorCell donor;
    donor.global_cell = 1000U + index * 17U;
    donor.global_index = {frame.anchor_global_cell.x + offset[index].x,
                          frame.anchor_global_cell.y + offset[index].y,
                          frame.anchor_global_cell.z + offset[index].z};
    donor.local_index = {static_cast<std::int32_t>(index), 0, 0};
    donor.centre = centre;
    donor.widths = {0.18 + 0.007 * static_cast<double>(index % 3U),
                    0.23 + 0.011 * static_cast<double>(index % 5U),
                    0.31 + 0.013 * static_cast<double>(index % 7U)};
    donors.push_back(donor);
  }
  return donors;
}

struct OwnedField {
  std::vector<double> storage;
  ConstFieldView view{};
};

OwnedField donor_field(const Polynomial& polynomial,
                       const QuadraticFrame& frame,
                       const std::vector<QuadraticDonorCell>& donors,
                       std::uint8_t components = 2U) {
  OwnedField field;
  field.storage.assign(donors.size() * components, -91.0);
  field.view.base = field.storage.data();
  field.view.interior = {static_cast<std::int32_t>(donors.size()), 1, 1};
  field.view.ghosts = {0, 0, 0};
  field.view.components = components;
  field.view.stride_y = donors.size();
  field.view.stride_z = donors.size();
  field.view.component_stride = donors.size();
  field.view.field = 7U;
  field.view.revision = 1U;
  field.view.storage_identity = 41U;
  field.view.revision_domain = 43U;
  for (std::size_t index = 0U; index < donors.size(); ++index) {
    field.storage[index] = donor_average(polynomial, frame, donors[index]);
    field.storage[donors.size() + index] =
        3.0 * field.storage[index] - 0.25;
  }
  return field;
}

std::vector<QuadraticFunctionalRequest> functionals(
    const QuadraticFrame& frame) {
  const Real3 value_point =
      add(add(add(frame.origin, multiply(frame.scale * -0.6, frame.normal)),
              multiply(frame.scale * 0.35, frame.tangent1)),
          multiply(frame.scale * -0.25, frame.tangent2));
  const Real3 derivative_point =
      add(add(add(frame.origin, multiply(frame.scale * -0.3, frame.normal)),
              multiply(frame.scale * -0.2, frame.tangent1)),
          multiply(frame.scale * 0.4, frame.tangent2));
  const Real3 derivative_direction =
      add(add(multiply(0.5, frame.normal), multiply(-0.25, frame.tangent1)),
          multiply(0.75, frame.tangent2));
  return {
      {QuadraticFunctionalKind::value, QuadraticConstraint::none,
       value_point, {}},
      {QuadraticFunctionalKind::value, QuadraticConstraint::origin_value,
       value_point, {}},
      {QuadraticFunctionalKind::directional_derivative,
       QuadraticConstraint::origin_normal_gradient, derivative_point,
       derivative_direction},
  };
}

bool compile_one(const QuadraticFrame& frame,
                 const std::vector<QuadraticDonorCell>& donors,
                 const std::vector<QuadraticFunctionalRequest>& rows,
                 QuadraticStencilPlan& plan,
                 QuadraticStencilLimits limits = {}) {
  const QuadraticStencilRequest request{
      frame, {donors.data(), donors.size()}, {rows.data(), rows.size()}};
  return static_cast<bool>(QuadraticStencilCompiler::compile(
      {&request, 1U}, limits, plan));
}

bool test_cell_average_quadratic_and_affine_constraints() {
  const QuadraticFrame frame = rotated_frame();
  const auto donors = complete_donors(frame);
  const auto rows = functionals(frame);
  QuadraticStencilPlan plan;
  bool passed = expect(compile_one(frame, donors, rows, plan),
                       "valid rotated quadratic stencil compiles");
  if (!passed) {
    return false;
  }
  passed &= expect(plan.groups().size == 1U && plan.rows().size == 3U &&
                       plan.donor_global_cells().size == donors.size() &&
                       plan.weights().size == donors.size() * rows.size(),
                   "one compact group owns all donor weights and rows");
  const auto quality = plan.groups().data[0U].quality;
  passed &= expect(quality.donor_count == donors.size() &&
                       quality.normal_band_count >= 3U &&
                       quality.quadrant_mask == 0x0fU &&
                       quality.reach == 4U && quality.rank == 10U &&
                       quality.condition_estimate > 1.0 &&
                       quality.condition_estimate <= 1.0e8 &&
                       quality.functional_l1 > 0.0 &&
                       quality.pivot_fingerprint != 0U &&
                       plan.fingerprint() != 0U,
                   "published quality records every hard stencil gate");

  const Polynomial polynomial;
  const OwnedField field = donor_field(polynomial, frame, donors);
  const auto local_value_coordinates = std::array<double, 3U>{-0.6, 0.35,
                                                               -0.25};
  const double expected_value =
      polynomial_value(polynomial, local_value_coordinates[0U],
                       local_value_coordinates[1U],
                       local_value_coordinates[2U]);
  const double wall_value = polynomial.coefficient[0U];
  const double wall_gradient = polynomial.coefficient[1U] / frame.scale;

  double actual = -777.0;
  passed &= expect(static_cast<bool>(evaluate_quadratic_row(
                       plan, 0U, field.view, 0U, 99.0, -88.0, actual)) &&
                       close(actual, expected_value),
                   "unconstrained row exactly reconstructs a cell-average "
                   "quadratic");
  actual = -777.0;
  passed &= expect(static_cast<bool>(evaluate_quadratic_row(
                       plan, 1U, field.view, 0U, wall_value, -88.0,
                       actual)) &&
                       close(actual, expected_value),
                   "Dirichlet affine row uses the exact origin value");

  const std::array<double, 3U> derivative_coordinates{-0.3, -0.2, 0.4};
  const std::array<double, 3U> local_direction{
      0.5 / frame.scale, -0.25 / frame.scale, 0.75 / frame.scale};
  const double expected_derivative = polynomial_directional_derivative(
      polynomial, derivative_coordinates[0U], derivative_coordinates[1U],
      derivative_coordinates[2U], local_direction[0U], local_direction[1U],
      local_direction[2U]);
  actual = -777.0;
  passed &= expect(static_cast<bool>(evaluate_quadratic_row(
                       plan, 2U, field.view, 0U, -99.0, wall_gradient,
                       actual)) &&
                       close(actual, expected_derivative),
                   "Neumann affine row uses the physical normal gradient");

  actual = -777.0;
  const double expected_component_one = 3.0 * expected_value - 0.25;
  passed &= expect(static_cast<bool>(evaluate_quadratic_row(
                       plan, 0U, field.view, 1U, 0.0, 0.0, actual)) &&
                       close(actual, expected_component_one),
                   "row evaluation respects component-major field layout");
  return passed;
}

bool same_bytes(Span<const double> left, Span<const double> right) {
  return left.size == right.size &&
         (left.size == 0U ||
          std::memcmp(left.data, right.data, left.size * sizeof(double)) == 0);
}

template <class T>
bool same_bytes(Span<const T> left, Span<const T> right) {
  return left.size == right.size &&
         (left.size == 0U ||
          std::memcmp(left.data, right.data, left.size * sizeof(T)) == 0);
}

bool test_donor_permutation_is_bitwise_deterministic() {
  const QuadraticFrame frame = rotated_frame();
  const auto donors = complete_donors(frame);
  auto reversed = donors;
  std::reverse(reversed.begin(), reversed.end());
  const auto rows = functionals(frame);
  QuadraticStencilPlan forward_plan;
  QuadraticStencilPlan reverse_plan;
  bool passed = expect(compile_one(frame, donors, rows, forward_plan) &&
                           compile_one(frame, reversed, rows, reverse_plan),
                       "both donor orders compile");
  if (!passed) {
    return false;
  }
  passed &= expect(forward_plan.fingerprint() == reverse_plan.fingerprint(),
                   "canonical donor order produces one plan identity");
  passed &= expect(same_bytes(forward_plan.donor_global_cells(),
                              reverse_plan.donor_global_cells()) &&
                       same_bytes(forward_plan.donor_local_indices(),
                                  reverse_plan.donor_local_indices()) &&
                       same_bytes(forward_plan.weights(),
                                  reverse_plan.weights()),
                   "canonical donor order produces bitwise-identical payload");
  passed &= expect(
      forward_plan.groups().data[0U].quality.pivot_fingerprint ==
              reverse_plan.groups().data[0U].quality.pivot_fingerprint &&
          forward_plan.groups().data[0U].fingerprint ==
              reverse_plan.groups().data[0U].fingerprint,
      "QR pivot and group fingerprints ignore caller permutation");
  return passed;
}

std::vector<QuadraticDonorCell> scored_donors(
    const QuadraticFrame& frame) {
  std::uint64_t state = 1U;
  auto uniform = [&](double lower, double upper) {
    state = state * UINT64_C(6364136223846793005) +
            UINT64_C(1442695040888963407);
    const double unit = static_cast<double>(state >> 11U) /
                        static_cast<double>(UINT64_C(1) << 53U);
    return lower + (upper - lower) * unit;
  };
  std::vector<QuadraticDonorCell> donors;
  donors.reserve(32U);
  for (std::size_t index = 0U; index < 32U; ++index) {
    const double n = uniform(0.15, 4.5);
    const double t1 = uniform(-4.0, 4.0);
    const double t2 = uniform(-4.0, 4.0);
    Real3 centre = frame.origin;
    centre = add(centre, multiply(frame.scale * n, frame.normal));
    centre = add(centre, multiply(frame.scale * t1, frame.tangent1));
    centre = add(centre, multiply(frame.scale * t2, frame.tangent2));
    const auto offset = [](std::size_t value) {
      return static_cast<std::int32_t>(value);
    };
    donors.push_back(
        {1000U + index,
         {frame.anchor_global_cell.x + offset(index % 9U) - 4,
          frame.anchor_global_cell.y + offset((index / 9U) % 4U) - 2,
          frame.anchor_global_cell.z + offset(index / 27U) - 1},
         {static_cast<std::int32_t>(index), 0, 0},
         centre,
         {uniform(0.05, 0.8), uniform(0.05, 0.8), uniform(0.05, 0.8)}});
  }
  return donors;
}

bool test_subset_search_rejects_bad_prefixes_and_scores_all_legal_counts() {
  const QuadraticFrame frame = rotated_frame();
  const auto rows = functionals(frame);
  const auto donors = scored_donors(frame);

  QuadraticStencilPlan selected;
  bool passed = expect(compile_one(frame, donors, rows, selected),
                       "deterministic prefix family compiles");
  if (!selected.fingerprint()) {
    return false;
  }
  const auto& selected_quality = selected.groups().data[0U].quality;
  passed &= expect(selected_quality.donor_count > 14U &&
                       selected_quality.donor_count < donors.size(),
                   "subset search publishes neither the minimum nor "
                   "an unconditional all-donor stencil");
  const auto selected_ids = selected.donor_global_cells();
  passed &= expect(std::is_sorted(selected_ids.data,
                                  selected_ids.data + selected_ids.size) &&
                       std::adjacent_find(selected_ids.data,
                                          selected_ids.data +
                                              selected_ids.size) ==
                           selected_ids.data + selected_ids.size,
                   "published subset donor IDs are strictly increasing");

  QuadraticStencilLimits full_only;
  full_only.minimum_donors = 32U;
  QuadraticStencilPlan full;
  passed &= expect(compile_one(frame, donors, rows, full, full_only),
                   "all-donor comparison candidate compiles");
  if (!full.fingerprint()) {
    return false;
  }
  passed &= expect(selected_quality.functional_l1 <
                       full.groups().data[0U].quality.functional_l1,
                   "maximum functional L1 score selects a better legal "
                   "subset than the 32-donor candidate");

  auto reversed = donors;
  std::reverse(reversed.begin(), reversed.end());
  QuadraticStencilPlan reversed_plan;
  passed &= expect(compile_one(frame, reversed, rows, reversed_plan) &&
                       reversed_plan.fingerprint() == selected.fingerprint() &&
                       same_bytes(reversed_plan.donor_global_cells(),
                                  selected.donor_global_cells()) &&
                       same_bytes(reversed_plan.weights(), selected.weights()),
                   "subset family, score and payload are independent of "
                   "input donor permutation");
  return passed;
}

bool test_geometry_mutation_changes_plan_identity() {
  const QuadraticFrame frame = rotated_frame();
  const auto rows = functionals(frame);
  const auto donors = complete_donors(frame);
  QuadraticStencilPlan reference;
  bool passed = expect(compile_one(frame, donors, rows, reference),
                       "geometry identity reference compiles");
  auto changed = donors;
  changed.back().widths.x *= 1.03125;
  QuadraticStencilPlan mutated;
  passed &= expect(
      compile_one(frame, changed, rows, mutated) &&
          mutated.fingerprint() != reference.fingerprint() &&
          !same_bytes(mutated.weights(), reference.weights()),
      "cell-width mutation changes weights and public plan identity");
  return passed;
}

bool rejected_without_replacing(
    const QuadraticFrame& frame,
    const std::vector<QuadraticDonorCell>& invalid_donors,
    const std::vector<QuadraticFunctionalRequest>& rows,
    QuadraticStencilLimits limits, QuadraticStencilPlan& retained,
    PlanFingerprint expected_fingerprint, std::string_view description) {
  const QuadraticStencilRequest request{
      frame, {invalid_donors.data(), invalid_donors.size()},
      {rows.data(), rows.size()}};
  const Status status = QuadraticStencilCompiler::compile(
      {&request, 1U}, limits, retained);
  return expect(!status && retained.fingerprint() == expected_fingerprint,
                description);
}

bool test_hard_quality_gates_and_atomic_publication() {
  const QuadraticFrame frame = rotated_frame();
  const auto donors = complete_donors(frame);
  const auto rows = functionals(frame);
  QuadraticStencilPlan retained;
  bool passed = expect(compile_one(frame, donors, rows, retained),
                       "baseline stencil compiles before rejection tests");
  if (!passed) {
    return false;
  }
  const PlanFingerprint fingerprint = retained.fingerprint();

  auto too_few = donors;
  too_few.resize(13U);
  passed &= rejected_without_replacing(
      frame, too_few, rows, {}, retained, fingerprint,
      "fewer than fourteen donors is rejected atomically");

  auto duplicate = donors;
  duplicate.back().global_cell = duplicate.front().global_cell;
  passed &= rejected_without_replacing(
      frame, duplicate, rows, {}, retained, fingerprint,
      "duplicate global donor is rejected atomically");

  auto duplicate_global_index = donors;
  duplicate_global_index.back().global_index =
      duplicate_global_index.front().global_index;
  passed &= rejected_without_replacing(
      frame, duplicate_global_index, rows, {}, retained, fingerprint,
      "duplicate physical donor index is rejected atomically");

  auto duplicate_local_index = donors;
  duplicate_local_index.back().local_index =
      duplicate_local_index.front().local_index;
  passed &= rejected_without_replacing(
      frame, duplicate_local_index, rows, {}, retained, fingerprint,
      "duplicate field-view donor index is rejected atomically");

  auto nonpositive = donors;
  nonpositive.front().centre =
      add(frame.origin, multiply(-0.1 * frame.scale, frame.normal));
  passed &= rejected_without_replacing(
      frame, nonpositive, rows, {}, retained, fingerprint,
      "non-positive normal donor is rejected atomically");

  auto two_bands = donors;
  for (std::size_t index = 0U; index < two_bands.size(); ++index) {
    const Real3 original_offset{two_bands[index].centre.x - frame.origin.x,
                                two_bands[index].centre.y - frame.origin.y,
                                two_bands[index].centre.z - frame.origin.z};
    const double original_n = dot(original_offset, frame.normal) / frame.scale;
    const double replacement_n = index % 2U == 0U ? 0.5 : 1.5;
    two_bands[index].centre = add(
        two_bands[index].centre,
        multiply((replacement_n - original_n) * frame.scale, frame.normal));
  }
  passed &= rejected_without_replacing(
      frame, two_bands, rows, {}, retained, fingerprint,
      "fewer than three normal bands is rejected atomically");

  auto missing_quadrants = donors;
  for (auto& donor : missing_quadrants) {
    const Real3 offset{donor.centre.x - frame.origin.x,
                       donor.centre.y - frame.origin.y,
                       donor.centre.z - frame.origin.z};
    const double t1 = dot(offset, frame.tangent1) / frame.scale;
    if (t1 < 0.0) {
      donor.centre = add(donor.centre,
                         multiply(-2.0 * t1 * frame.scale, frame.tangent1));
    }
  }
  passed &= rejected_without_replacing(
      frame, missing_quadrants, rows, {}, retained, fingerprint,
      "missing tangential quadrants is rejected atomically");

  auto excessive_reach = donors;
  excessive_reach.front().global_index.x =
      frame.anchor_global_cell.x + 5;
  passed &= rejected_without_replacing(
      frame, excessive_reach, rows, {}, retained, fingerprint,
      "donor reach beyond four cells is rejected atomically");

  auto rank_deficient = donors;
  for (std::size_t index = 1U; index < rank_deficient.size(); ++index) {
    rank_deficient[index].centre = rank_deficient[index % 4U].centre;
    rank_deficient[index].widths = rank_deficient[index % 4U].widths;
  }
  passed &= rejected_without_replacing(
      frame, rank_deficient, rows, {}, retained, fingerprint,
      "rank-deficient quadratic basis is rejected atomically");

  auto ill_conditioned = donors;
  for (auto& donor : ill_conditioned) {
    const Real3 offset{donor.centre.x - frame.origin.x,
                       donor.centre.y - frame.origin.y,
                       donor.centre.z - frame.origin.z};
    const double t1 = dot(offset, frame.tangent1) / frame.scale;
    const double t2 = dot(offset, frame.tangent2) / frame.scale;
    donor.centre = add(
        donor.centre,
        add(multiply((-t1 + t1 * 1.0e-10) * frame.scale, frame.tangent1),
            multiply((-t2 + t2 * 1.0e-10) * frame.scale, frame.tangent2)));
    donor.widths = {1.0e-12, 1.0e-12, 1.0e-12};
  }
  passed &= rejected_without_replacing(
      frame, ill_conditioned, rows, {}, retained, fingerprint,
      "condition estimate above 1e8 is rejected atomically");

  QuadraticStencilLimits weak_limits;
  weak_limits.minimum_donors = 10U;
  passed &= rejected_without_replacing(
      frame, donors, rows, weak_limits, retained, fingerprint,
      "limits cannot weaken the fourteen-donor hard floor");
  weak_limits = {};
  weak_limits.maximum_reach = 5U;
  passed &= rejected_without_replacing(
      frame, donors, rows, weak_limits, retained, fingerprint,
      "limits cannot weaken the four-cell hard reach");
  weak_limits = {};
  weak_limits.condition_limit = 1.0e9;
  passed &= rejected_without_replacing(
      frame, donors, rows, weak_limits, retained, fingerprint,
      "limits cannot weaken the 1e8 hard condition ceiling");

  QuadraticFrame left_handed = frame;
  left_handed.tangent2 = multiply(-1.0, left_handed.tangent2);
  passed &= rejected_without_replacing(
      left_handed, donors, rows, {}, retained, fingerprint,
      "left-handed orthonormal frame is rejected atomically");

  auto zero_local_scale = donors;
  zero_local_scale.front().widths.x = 0.0;
  passed &= rejected_without_replacing(
      frame, zero_local_scale, rows, {}, retained, fingerprint,
      "non-positive donor width is rejected atomically");
  return passed;
}

bool test_evaluation_rejects_invalid_views_atomically() {
  const QuadraticFrame frame = rotated_frame();
  const auto donors = complete_donors(frame);
  const auto rows = functionals(frame);
  QuadraticStencilPlan plan;
  bool passed = expect(compile_one(frame, donors, rows, plan),
                       "evaluation fixture compiles");
  if (!passed) {
    return false;
  }
  const Polynomial polynomial;
  OwnedField field = donor_field(polynomial, frame, donors);
  double output = 12345.0;
  passed &= expect(!evaluate_quadratic_row(plan, 99U, field.view, 0U, 0.0,
                                           0.0, output) &&
                       output == 12345.0,
                   "invalid row preserves caller output");
  passed &= expect(!evaluate_quadratic_row(plan, 0U, field.view, 2U, 0.0,
                                           0.0, output) &&
                       output == 12345.0,
                   "invalid component preserves caller output");
  field.storage[0U] = std::numeric_limits<double>::quiet_NaN();
  passed &= expect(!evaluate_quadratic_row(plan, 0U, field.view, 0U, 0.0,
                                           0.0, output) &&
                       output == 12345.0,
                   "non-finite donor value preserves caller output");
  return passed;
}

bool test_multiple_groups_keep_offsets_independent() {
  const QuadraticFrame first_frame = rotated_frame();
  QuadraticFrame second_frame = first_frame;
  second_frame.origin = add(first_frame.origin, {3.0, -1.0, 0.5});
  second_frame.anchor_global_cell = {30, 20, 15};
  const auto first_donors = complete_donors(first_frame);
  auto second_donors = complete_donors(second_frame);
  for (std::size_t index = 0U; index < second_donors.size(); ++index) {
    second_donors[index].global_cell += 10000U;
    second_donors[index].local_index = {
        static_cast<std::int32_t>(first_donors.size() + index), 0, 0};
  }
  const auto first_rows = functionals(first_frame);
  const auto second_rows = functionals(second_frame);
  const std::array<QuadraticStencilRequest, 2U> requests{{
      {first_frame, {first_donors.data(), first_donors.size()},
       {first_rows.data(), first_rows.size()}},
      {second_frame, {second_donors.data(), second_donors.size()},
       {second_rows.data(), second_rows.size()}},
  }};
  QuadraticStencilPlan plan;
  bool passed = expect(static_cast<bool>(QuadraticStencilCompiler::compile(
                           {requests.data(), requests.size()}, {}, plan)),
                       "two quadratic groups compile in one immutable plan");
  if (!passed) {
    return false;
  }
  passed &= expect(plan.groups().size == 2U && plan.rows().size == 6U &&
                       plan.groups().data[1U].donor_begin ==
                           first_donors.size() &&
                       plan.groups().data[1U].row_begin == first_rows.size() &&
                       plan.rows().data[3U].group == 1U,
                   "second group records nonzero donor and row offsets");

  const Polynomial polynomial;
  OwnedField field;
  field.storage.assign((first_donors.size() + second_donors.size()) * 2U,
                       0.0);
  for (std::size_t index = 0U; index < first_donors.size(); ++index) {
    field.storage[index] =
        donor_average(polynomial, first_frame, first_donors[index]);
  }
  for (std::size_t index = 0U; index < second_donors.size(); ++index) {
    field.storage[first_donors.size() + index] =
        donor_average(polynomial, second_frame, second_donors[index]);
  }
  field.view.base = field.storage.data();
  field.view.interior = {
      static_cast<std::int32_t>(first_donors.size() + second_donors.size()),
      1, 1};
  field.view.ghosts = {0, 0, 0};
  field.view.components = 2U;
  field.view.stride_y = first_donors.size() + second_donors.size();
  field.view.stride_z = first_donors.size() + second_donors.size();
  field.view.component_stride = first_donors.size() + second_donors.size();
  field.view.field = 7U;
  field.view.revision = 1U;
  field.view.storage_identity = 41U;
  field.view.revision_domain = 43U;
  double actual = -1.0;
  const double expected = polynomial_value(polynomial, -0.6, 0.35, -0.25);
  passed &= expect(static_cast<bool>(evaluate_quadratic_row(
                       plan, 3U, field.view, 0U, 0.0, 0.0, actual)) &&
                       close(actual, expected),
                   "second group row reads its own donor slice");
  return passed;
}

}  // namespace

int main() {
  bool passed = true;
  passed &= test_cell_average_quadratic_and_affine_constraints();
  passed &= test_donor_permutation_is_bitwise_deterministic();
  passed &=
      test_subset_search_rejects_bad_prefixes_and_scores_all_legal_counts();
  passed &= test_geometry_mutation_changes_plan_identity();
  passed &= test_hard_quality_gates_and_atomic_publication();
  passed &= test_evaluation_rejects_invalid_views_atomically();
  passed &= test_multiple_groups_keep_offsets_independent();
  if (!passed) {
    return 1;
  }
  std::cout << "v0.4 immersed quadratic reconstruction tests passed\n";
  return 0;
}
