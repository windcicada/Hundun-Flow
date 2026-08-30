// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_linear.hpp"

#include "solver_krylov_test_detail.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <string_view>
#include <vector>

namespace {
using namespace hundun::v04;

constexpr Int3 kOddCells{5, 3, 7};
constexpr std::uint8_t kGhosts = 1U;
constexpr std::size_t kMaximumCount =
    detail::kMaximumFgmresBasisUpdateCountForTest;

struct OwnedField {
  std::vector<double> storage;
  FieldView view{};
};

OwnedField make_field(Int3 cells, std::uint8_t ghosts,
                      std::size_t row_padding,
                      std::size_t plane_padding,
                      std::uint64_t storage_identity) {
  OwnedField result;
  const std::size_t padded_x =
      static_cast<std::size_t>(cells.x) + 2U * ghosts + row_padding;
  const std::size_t padded_y =
      static_cast<std::size_t>(cells.y) + 2U * ghosts;
  const std::size_t padded_z = static_cast<std::size_t>(cells.z) + 2U * ghosts;
  const std::size_t stride_z = padded_x * padded_y + plane_padding;
  const std::size_t component_stride = stride_z * padded_z;
  result.storage.assign(component_stride, 0.0);
  result.view.base = result.storage.data() + ghosts +
                     static_cast<std::size_t>(ghosts) *
                         (padded_x + stride_z);
  result.view.interior = cells;
  result.view.ghosts = {ghosts, ghosts, ghosts};
  result.view.components = 1U;
  result.view.stride_y = padded_x;
  result.view.stride_z = stride_z;
  result.view.component_stride = component_stride;
  result.view.field = 1U;
  result.view.revision = 1U;
  result.view.storage_identity = storage_identity;
  result.view.revision_domain = 1U;
  return result;
}

OwnedField make_field(std::uint64_t storage_identity) {
  return make_field(kOddCells, kGhosts, 0U, 0U, storage_identity);
}

void fill_field(FieldView field, double seed) {
  for (std::int32_t z = 0; z < field.interior.z; ++z) {
    for (std::int32_t y = 0; y < field.interior.y; ++y) {
      for (std::int32_t x = 0; x < field.interior.x; ++x) {
        field.unchecked({x, y, z}, 0U) =
            seed + 0.03125 * static_cast<double>(x + 3 * y + 5 * z);
      }
    }
  }
}

Status scalar_reference_update(FieldView destination,
                               const ConstFieldView* sources,
                               const double* scales,
                               std::size_t count) noexcept {
  // This is the legacy sequence: every add_scaled call scans the complete
  // field before the next basis row starts.  Keeping row outside the cell
  // loops is essential to the finite bitwise oracle and to the documented
  // partial-workspace behavior on a nonfinite update.
  for (std::size_t row = 0U; row < count; ++row) {
    for (std::int32_t z = 0; z < destination.interior.z; ++z) {
      for (std::int32_t y = 0; y < destination.interior.y; ++y) {
        for (std::int32_t x = 0; x < destination.interior.x; ++x) {
          const double value =
              destination.unchecked({x, y, z}, 0U) +
              scales[row] * sources[row].unchecked({x, y, z}, 0U);
          if (!std::isfinite(value)) {
            return {StatusCode::numerical_failure, 605U};
          }
          destination.unchecked({x, y, z}, 0U) = value;
        }
      }
    }
  }
  return {};
}

bool same_storage(const OwnedField& left, const OwnedField& right) noexcept {
  return left.storage.size() == right.storage.size() &&
         std::memcmp(left.storage.data(), right.storage.data(),
                     left.storage.size() * sizeof(double)) == 0;
}

bool expect(bool condition, std::string_view description) {
  if (!condition) {
    std::cerr << "FAIL: " << description << '\n';
  }
  return condition;
}

double scale_for(std::size_t row) noexcept {
  switch (row % 6U) {
    case 0U:
      return 0.125 * static_cast<double>(row + 1U);
    case 1U:
      return -0.0625 * static_cast<double>(row + 1U);
    case 2U:
      return 0.0;
    case 3U:
      return -0.0;
    case 4U:
      return 1.0;
    default:
      return -1.0;
  }
}

bool test_finite_bitwise_equivalence() {
  std::vector<OwnedField> basis;
  basis.reserve(kMaximumCount);
  std::array<ConstFieldView, kMaximumCount> sources{};
  std::array<double, kMaximumCount> scales{};
  bool has_positive = false;
  bool has_negative = false;
  bool has_positive_zero = false;
  bool has_negative_zero = false;
  for (std::size_t row = 0U; row < kMaximumCount; ++row) {
    basis.push_back(make_field(1000U + row));
    fill_field(basis.back().view, 0.125 + 0.001 * static_cast<double>(row));
    sources[row] = as_const(basis.back().view);
    switch (row % 6U) {
      case 0U:
        scales[row] = 0.125 * static_cast<double>(row + 1U);
        has_positive = true;
        break;
      case 1U:
        scales[row] = -0.0625 * static_cast<double>(row + 1U);
        has_negative = true;
        break;
      case 2U:
        scales[row] = 0.0;
        has_positive_zero = true;
        break;
      case 3U:
        scales[row] = -0.0;
        has_negative_zero = true;
        break;
      case 4U:
        scales[row] = 1.0;
        has_positive = true;
        break;
      default:
        scales[row] = -1.0;
        has_negative = true;
        break;
    }
  }

  bool passed = expect(kMaximumCount >= 12U,
                       "validated FGMRES update bound covers count 12") &&
                expect(has_positive && has_negative && has_positive_zero &&
                           has_negative_zero,
                       "finite oracle exercises positive, negative, and zero scales");
  for (std::size_t count = 1U; count <= kMaximumCount; ++count) {
    OwnedField fused = make_field(2000U + count);
    OwnedField scalar = make_field(3000U + count);
    fill_field(fused.view, -0.375);
    fill_field(scalar.view, -0.375);
    const Status fused_status = detail::fused_krylov_basis_update_for_test(
        fused.view, sources.data(), scales.data(), count);
    const Status scalar_status = scalar_reference_update(
        scalar.view, sources.data(), scales.data(), count);
    passed &= expect(static_cast<bool>(fused_status) &&
                         fused_status.code == scalar_status.code &&
                         fused_status.detail == scalar_status.detail,
                     "fused and scalar statuses match for finite count");
    passed &= expect(same_storage(fused, scalar),
                     "fused basis update is bitwise equal to scalar reference");
  }
  return passed;
}

bool test_strip_boundaries_and_distinct_layouts() {
  constexpr std::array<std::int32_t, 5U> widths{1, 63, 64, 65, 129};
  bool passed = true;
  for (const std::int32_t width : widths) {
    const Int3 cells{width, 3, 2};
    std::vector<OwnedField> basis;
    basis.reserve(kMaximumCount);
    std::array<ConstFieldView, kMaximumCount> sources{};
    std::array<double, kMaximumCount> scales{};
    for (std::size_t row = 0U; row < kMaximumCount; ++row) {
      basis.push_back(make_field(
          cells, static_cast<std::uint8_t>(1U + row % 2U), row % 5U,
          row % 7U, 6000U + 1000U * static_cast<std::uint64_t>(width) + row));
      fill_field(basis.back().view,
                 0.0625 + 0.00025 * static_cast<double>(row));
      sources[row] = as_const(basis.back().view);
      scales[row] = scale_for(row);
    }

    for (std::size_t count = 1U; count <= kMaximumCount; ++count) {
      OwnedField striped = make_field(cells, 2U, 3U, 5U,
                                      700000U + count + width);
      OwnedField scalar = make_field(cells, 2U, 3U, 5U,
                                     800000U + count + width);
      fill_field(striped.view, -0.4375);
      fill_field(scalar.view, -0.4375);
      const Status striped_status = detail::fused_krylov_basis_update_for_test(
          striped.view, sources.data(), scales.data(), count);
      const Status scalar_status = scalar_reference_update(
          scalar.view, sources.data(), scales.data(), count);
      passed &= expect(static_cast<bool>(striped_status) &&
                           striped_status.code == scalar_status.code &&
                           striped_status.detail == scalar_status.detail,
                       "strip-boundary status equals the scalar oracle");
      passed &= expect(
          same_storage(striped, scalar),
          "all counts and strip-boundary layouts are bitwise scalar-equal");
    }
  }
  return passed;
}

bool row_range_matches_snapshot(const OwnedField& field,
                                const std::vector<double>& before,
                                std::int32_t begin_x,
                                std::int32_t end_x) noexcept {
  for (std::int32_t x = begin_x; x < end_x; ++x) {
    const double* const value = field.view.base + x;
    const std::size_t offset =
        static_cast<std::size_t>(value - field.storage.data());
    if (std::memcmp(value, &before[offset], sizeof(double)) != 0) {
      return false;
    }
  }
  return true;
}

bool test_nonfinite_active_strip_is_not_published() {
  struct FailureCase {
    std::int32_t x{};
    std::size_t basis{};
    double value{};
  };
  const std::array<FailureCase, 4U> cases{
      FailureCase{0, 1U, std::numeric_limits<double>::infinity()},
      FailureCase{63, 2U, -std::numeric_limits<double>::infinity()},
      FailureCase{64, 1U, std::numeric_limits<double>::quiet_NaN()},
      FailureCase{128, 3U, std::numeric_limits<double>::infinity()}};
  const Int3 cells{129, 2, 2};
  bool passed = true;
  for (std::size_t case_index = 0U; case_index < cases.size(); ++case_index) {
    std::array<OwnedField, 4U> basis{};
    std::array<ConstFieldView, 4U> sources{};
    std::array<double, 4U> scales{1.0, -0.5, 0.25, 2.0};
    for (std::size_t row = 0U; row < basis.size(); ++row) {
      basis[row] = make_field(cells, static_cast<std::uint8_t>(1U + row % 2U),
                              row + 1U, row + 2U,
                              900000U + 100U * case_index + row);
      fill_field(basis[row].view, 0.5 + 0.125 * static_cast<double>(row));
      sources[row] = as_const(basis[row].view);
    }
    basis[cases[case_index].basis].view.unchecked(
        {cases[case_index].x, 0, 0}, 0U) = cases[case_index].value;
    sources[cases[case_index].basis] =
        as_const(basis[cases[case_index].basis].view);

    OwnedField striped =
        make_field(cells, 2U, 3U, 5U, 910000U + case_index);
    fill_field(striped.view, -0.75);
    const std::vector<double> before = striped.storage;
    const Status status = detail::fused_krylov_basis_update_for_test(
        striped.view, sources.data(), scales.data(), sources.size());
    const std::int32_t strip_begin = (cases[case_index].x / 64) * 64;
    const std::int32_t strip_end = std::min(strip_begin + 64, cells.x);
    passed &= expect(status.code == StatusCode::numerical_failure &&
                         status.detail == 605U,
                     "nonfinite strip returns the frozen status and detail");
    passed &= expect(
        row_range_matches_snapshot(striped, before, strip_begin, strip_end),
        "the active x strip is not published after a nonfinite basis row");
  }
  return passed;
}

bool test_nonfinite_failure_is_not_published() {
  std::vector<OwnedField> basis;
  basis.reserve(kMaximumCount);
  std::array<ConstFieldView, kMaximumCount> sources{};
  std::array<double, kMaximumCount> scales{};
  for (std::size_t row = 0U; row < kMaximumCount; ++row) {
    basis.push_back(make_field(4000U + row));
    fill_field(basis.back().view, 0.5);
    sources[row] = as_const(basis.back().view);
    scales[row] = row % 2U == 0U ? 1.0 : -1.0;
  }
  basis[1U].view.unchecked({0, 0, 0}, 0U) =
      std::numeric_limits<double>::infinity();
  sources[1U] = as_const(basis[1U].view);

  OwnedField fused = make_field(5001U);
  OwnedField scalar = make_field(5002U);
  fill_field(fused.view, -0.75);
  fill_field(scalar.view, -0.75);
  const std::vector<double> fused_before = fused.storage;
  const std::vector<double> scalar_before = scalar.storage;
  const Status fused_status = detail::fused_krylov_basis_update_for_test(
      fused.view, sources.data(), scales.data(), 12U);
  const Status scalar_status = scalar_reference_update(
      scalar.view, sources.data(), scales.data(), 12U);
  bool passed = expect(fused_status.code == StatusCode::numerical_failure &&
                           fused_status.code == scalar_status.code &&
                           fused_status.detail == scalar_status.detail,
                       "nonfinite mutation preserves scalar failure code/detail") &&
                expect(fused.storage == fused_before,
                       "fused nonfinite failure does not publish the partial cell") &&
                expect(scalar.storage != scalar_before,
                       "legacy nonfinite failure records its prior row partial state") &&
                expect(fused.storage != scalar.storage,
                       "partial-workspace difference is explicit and isolated before publication");
  return passed;
}

bool test_single_reduction_norm_arithmetic() {
  const std::array<double, 2U> projections{1.0, 2.0};
  double t = -1.0;
  double next_norm = -1.0;
  Status status = detail::single_reduction_fgmres_norm_for_test(
      Span<const double>{projections.data(), projections.size()}, 6.0, t,
      next_norm);
  bool passed = expect(status.code == StatusCode::ok && t == 1.0 &&
                           next_norm == 1.0,
                       "single-reduction norm accepts positive t exactly");

  t = -1.0;
  next_norm = -1.0;
  status = detail::single_reduction_fgmres_norm_for_test(
      Span<const double>{projections.data(), projections.size()}, 5.0, t,
      next_norm);
  passed &= expect(status.code == StatusCode::ok && t == 0.0 &&
                       next_norm == 0.0,
                   "single-reduction norm preserves exact zero happy breakdown");

  t = 0.0;
  next_norm = -1.0;
  status = detail::single_reduction_fgmres_norm_for_test(
      Span<const double>{projections.data(), projections.size()}, 4.0, t,
      next_norm);
  passed &= expect(status.code == StatusCode::numerical_failure &&
                       status.detail == 604U && t == -1.0 && next_norm == 0.0,
                   "strictly negative t is an unclamped norm breakdown");

  const std::array<double, 2U> nonfinite_projection{
      std::numeric_limits<double>::quiet_NaN(), 2.0};
  t = 0.0;
  next_norm = 0.0;
  status = detail::single_reduction_fgmres_norm_for_test(
      Span<const double>{nonfinite_projection.data(),
                         nonfinite_projection.size()},
      5.0, t, next_norm);
  passed &= expect(status.code == StatusCode::numerical_failure &&
                       status.detail == 605U && std::isnan(t) &&
                       std::isnan(next_norm),
                   "nonfinite projection is rejected before square root");

  t = 0.0;
  next_norm = 0.0;
  status = detail::single_reduction_fgmres_norm_for_test(
      Span<const double>{projections.data(), projections.size()},
      std::numeric_limits<double>::infinity(), t, next_norm);
  passed &= expect(status.code == StatusCode::numerical_failure &&
                       status.detail == 605U && std::isnan(t) &&
                       std::isnan(next_norm),
                   "nonfinite reduced norm is rejected before subtraction");
  return passed;
}

}  // namespace

int main() {
  return test_finite_bitwise_equivalence() &&
                 test_strip_boundaries_and_distinct_layouts() &&
                 test_nonfinite_failure_is_not_published() &&
                 test_nonfinite_active_strip_is_not_published() &&
                 test_single_reduction_norm_arithmetic()
             ? 0
             : 1;
}
