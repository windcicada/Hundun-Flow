// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_linear.hpp"

#include "field_view_interval_detail.hpp"
#include "solver_krylov_test_detail.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kLinearSolvePlan = 601U;
constexpr std::uint32_t kLinearSolveWorkspace = 602U;
constexpr std::uint32_t kLinearSolveBreakdown = 604U;
constexpr std::uint32_t kLinearSolveNonFinite = 605U;
constexpr std::uint32_t kLinearSolveConvergenceAudit = 606U;
// FGMRES owns 2 * restart + 8 vector slots.  The validated workspace
// contract stores that count in uint8_t, so restart and each Arnoldi update
// are bounded by (uint8_t_max - 8) / 2.
constexpr std::size_t kMaximumFgmresBasisUpdateCount =
    (std::numeric_limits<std::uint8_t>::max() - 8U) / 2U;
static_assert(kMaximumFgmresBasisUpdateCount == 123U);
constexpr std::size_t kFgmresBasisUpdateStripCells = 64U;
// Sanan--Schnepp--May Algorithm 14 uses a fixed private shift.  It is part of
// the registered method, not a case knob or a runtime spectral estimate.
constexpr double kSingleReductionFgmresShift = 1.0;
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
std::size_t g_forced_single_reduction_breakdown_count = 0U;
#endif

bool same_shape(Int3 left, Int3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

bool shape_contains(Int3 outer, Int3 inner) noexcept {
  return inner.x > 0 && inner.y > 0 && inner.z > 0 &&
         inner.x <= outer.x && inner.y <= outer.y && inner.z <= outer.z;
}

bool same_identity(LinearIdentity left, LinearIdentity right) noexcept {
  return left.symbolic == right.symbolic && left.numeric == right.numeric &&
         left.hierarchy == right.hierarchy &&
         left.workspace == right.workspace &&
         left.fingerprint == right.fingerprint;
}

template <class Left, class Right>
bool views_overlap(const BasicFieldView<Left>& left,
                   const BasicFieldView<Right>& right) noexcept {
  return detail::field_views_overlap(left, right);
}

bool valid_identity(LinearIdentity identity) noexcept {
  return identity.symbolic != 0U && identity.numeric != 0U &&
         identity.hierarchy != 0U && identity.workspace != 0U &&
         identity.fingerprint != 0U;
}

bool valid_preconditioner_status_scope(
    LinearPreconditionerStatusScope scope) noexcept {
  switch (scope) {
    case LinearPreconditionerStatusScope::rank_local:
    case LinearPreconditionerStatusScope::collective:
      return true;
  }
  return false;
}

bool matching_collective_operator_failure(
    Status status,
    LinearOperatorFailureProvenance provenance) noexcept {
  return !status && !provenance.status &&
         provenance.status.code == status.code &&
         provenance.status.detail == status.detail &&
         provenance.status_scope == LinearOperatorStatusScope::collective &&
         provenance.lowest_failing_rank >= 0;
}

bool valid_preconditioner_apply_lifecycle(
    LinearPreconditionerApplyLifecycle lifecycle) noexcept {
  switch (lifecycle) {
    case LinearPreconditionerApplyLifecycle::per_call_checked:
    case LinearPreconditionerApplyLifecycle::prepared_batch:
      return true;
  }
  return false;
}

std::uint64_t mix_contract(std::uint64_t hash, std::uint64_t value) noexcept {
  constexpr std::uint64_t prime = UINT64_C(1099511628211);
  hash ^= value;
  hash *= prime;
  return hash;
}

std::uint64_t double_bits(double value) noexcept {
  std::uint64_t result = 0U;
  static_assert(sizeof(result) == sizeof(value));
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

PlanFingerprint solve_contract_fingerprint(
    LinearAlgorithm algorithm, const LinearOperatorCertificate& op,
    const LinearPreconditionerCertificate& preconditioner,
    const LinearSolveInvocation& invocation,
    const LinearWorkspaceRequirements& workspace) noexcept {
  std::uint64_t hash = UINT64_C(1469598103934665603);
  hash = mix_contract(hash, static_cast<std::uint64_t>(algorithm));
  hash = mix_contract(hash, invocation.expected_identity.symbolic);
  hash = mix_contract(hash, invocation.expected_identity.hierarchy);
  hash = mix_contract(hash, invocation.expected_identity.numeric != 0U);
  hash = mix_contract(hash, invocation.expected_identity.workspace != 0U);
  hash = mix_contract(hash, invocation.expected_identity.fingerprint != 0U);
  hash = mix_contract(hash, op.collective_fingerprint);
  hash = mix_contract(hash, preconditioner.collective_fingerprint);
  hash = mix_contract(hash, static_cast<std::uint64_t>(op.operator_class));
  hash = mix_contract(
      hash, static_cast<std::uint64_t>(preconditioner.preconditioner_class));
  hash = mix_contract(
      hash, static_cast<std::uint64_t>(preconditioner.status_scope));
  hash = mix_contract(
      hash, static_cast<std::uint64_t>(preconditioner.apply_lifecycle));
  hash = mix_contract(hash, double_bits(invocation.control.absolute_tolerance));
  hash = mix_contract(hash, double_bits(invocation.control.relative_tolerance));
  hash = mix_contract(hash, invocation.control.maximum_iterations);
  hash = mix_contract(hash, invocation.control.true_residual_interval);
  hash = mix_contract(hash, invocation.control.restart);
  hash = mix_contract(
      hash, invocation.convergence_audit == nullptr
                ? 0U
                : invocation.convergence_audit->certificate()
                      .collective_fingerprint);
  hash = mix_contract(hash, workspace.maximum_restart);
  hash = mix_contract(hash, workspace.reduction_capacity);
  hash = mix_contract(hash, static_cast<std::uint64_t>(workspace.reduction_mode));
  hash = mix_contract(hash, workspace.execution_revision);
  return hash == 0U ? PlanFingerprint{1U} : hash;
}

template <class T>
bool valid_scalar_view(const BasicFieldView<T>& view) noexcept {
  if (view.base == nullptr || view.interior.x <= 0 || view.interior.y <= 0 ||
      view.interior.z <= 0 || view.ghosts.x < 0 || view.ghosts.y < 0 ||
      view.ghosts.z < 0 || view.components != 1U || view.stride_y == 0U ||
      view.stride_z == 0U || view.component_stride == 0U ||
      view.storage_identity == 0U || view.revision_domain == 0U) {
    return false;
  }
  const std::size_t padded_x =
      static_cast<std::size_t>(view.interior.x) +
      2U * static_cast<std::size_t>(view.ghosts.x);
  const std::size_t padded_y =
      static_cast<std::size_t>(view.interior.y) +
      2U * static_cast<std::size_t>(view.ghosts.y);
  if (view.stride_y < padded_x ||
      padded_y > std::numeric_limits<std::size_t>::max() / view.stride_y) {
    return false;
  }
  return view.stride_z >= view.stride_y * padded_y;
}

template <class Left, class Right>
Status local_dot(const BasicFieldView<Left>& left,
                 const BasicFieldView<Right>& right,
                 double& result) noexcept {
  double sum = 0.0;
  double correction = 0.0;
  for (std::int32_t z = 0; z < left.interior.z; ++z) {
    for (std::int32_t y = 0; y < left.interior.y; ++y) {
      for (std::int32_t x = 0; x < left.interior.x; ++x) {
        const double product = left.unchecked({x, y, z}, 0U) *
                               right.unchecked({x, y, z}, 0U);
        if (!std::isfinite(product)) {
          return {StatusCode::numerical_failure, kLinearSolveNonFinite};
        }
        const double next = sum + product;
        if (std::abs(sum) >= std::abs(product)) {
          correction += (sum - next) + product;
        } else {
          correction += (product - next) + sum;
        }
        sum = next;
      }
    }
  }
  result = sum + correction;
  return std::isfinite(result)
             ? Status{}
             : Status{StatusCode::numerical_failure,
                      kLinearSolveNonFinite};
}

template <class Left, class Right>
Status local_dot_pair(const BasicFieldView<Left>& left,
                      const BasicFieldView<Right>& right,
                      double& first_result, double& second_result) noexcept {
  double first_sum = 0.0;
  double first_correction = 0.0;
  double second_sum = 0.0;
  double second_correction = 0.0;
  for (std::int32_t z = 0; z < left.interior.z; ++z) {
    for (std::int32_t y = 0; y < left.interior.y; ++y) {
      for (std::int32_t x = 0; x < left.interior.x; ++x) {
        const double first_product =
            left.unchecked({x, y, z}, 0U) *
            right.unchecked({x, y, z}, 0U);
        if (!std::isfinite(first_product)) {
          return {StatusCode::numerical_failure, kLinearSolveNonFinite};
        }
        const double second_product =
            left.unchecked({x, y, z}, 0U) *
            left.unchecked({x, y, z}, 0U);
        if (!std::isfinite(second_product)) {
          return {StatusCode::numerical_failure, kLinearSolveNonFinite};
        }

        const double first_next = first_sum + first_product;
        if (std::abs(first_sum) >= std::abs(first_product)) {
          first_correction += (first_sum - first_next) + first_product;
        } else {
          first_correction += (first_product - first_next) + first_sum;
        }
        first_sum = first_next;

        const double second_next = second_sum + second_product;
        if (std::abs(second_sum) >= std::abs(second_product)) {
          second_correction +=
              (second_sum - second_next) + second_product;
        } else {
          second_correction +=
              (second_product - second_next) + second_sum;
        }
        second_sum = second_next;
      }
    }
  }
  const double first = first_sum + first_correction;
  const double second = second_sum + second_correction;
  if (!std::isfinite(first) || !std::isfinite(second)) {
    return {StatusCode::numerical_failure, kLinearSolveNonFinite};
  }
  first_result = first;
  second_result = second;
  return {};
}

template <class Source>
Status copy_field(const BasicFieldView<Source>& source,
                  FieldView destination) noexcept {
  for (std::int32_t z = 0; z < source.interior.z; ++z) {
    for (std::int32_t y = 0; y < source.interior.y; ++y) {
      for (std::int32_t x = 0; x < source.interior.x; ++x) {
        const double value = source.unchecked({x, y, z}, 0U);
        if (!std::isfinite(value)) {
          return {StatusCode::numerical_failure, kLinearSolveNonFinite};
        }
        destination.unchecked({x, y, z}, 0U) = value;
      }
    }
  }
  return {};
}

void fill_field(FieldView field, double value) noexcept {
  for (std::int32_t z = 0; z < field.interior.z; ++z) {
    for (std::int32_t y = 0; y < field.interior.y; ++y) {
      for (std::int32_t x = 0; x < field.interior.x; ++x) {
        field.unchecked({x, y, z}, 0U) = value;
      }
    }
  }
}

Status reset_to_zero_initial_guess(ConstFieldView rhs, FieldView solution,
                                   FieldView residual) noexcept {
  for (std::int32_t z = 0; z < rhs.interior.z; ++z) {
    for (std::int32_t y = 0; y < rhs.interior.y; ++y) {
      for (std::int32_t x = 0; x < rhs.interior.x; ++x) {
        const double value = rhs.unchecked({x, y, z}, 0U);
        solution.unchecked({x, y, z}, 0U) = 0.0;
        residual.unchecked({x, y, z}, 0U) = value;
      }
    }
  }
  return {};
}

Status axpby(FieldView destination, double left_scale,
             ConstFieldView left, double right_scale,
             ConstFieldView right) noexcept {
  for (std::int32_t z = 0; z < destination.interior.z; ++z) {
    for (std::int32_t y = 0; y < destination.interior.y; ++y) {
      for (std::int32_t x = 0; x < destination.interior.x; ++x) {
        const double value =
            left_scale * left.unchecked({x, y, z}, 0U) +
            right_scale * right.unchecked({x, y, z}, 0U);
        if (!std::isfinite(value)) {
          return {StatusCode::numerical_failure, kLinearSolveNonFinite};
        }
        destination.unchecked({x, y, z}, 0U) = value;
      }
    }
  }
  return {};
}

Status form_bicgstab_intermediate(FieldView destination,
                                  ConstFieldView residual,
                                  ConstFieldView image, double alpha,
                                  double& produced_local_scale) noexcept {
  double local_scale = 0.0;
  for (std::int32_t z = 0; z < destination.interior.z; ++z) {
    for (std::int32_t y = 0; y < destination.interior.y; ++y) {
      for (std::int32_t x = 0; x < destination.interior.x; ++x) {
        const Int3 cell{x, y, z};
        const double value =
            1.0 * residual.unchecked(cell, 0U) +
            -alpha * image.unchecked(cell, 0U);
        if (!std::isfinite(value)) {
          return {StatusCode::numerical_failure, kLinearSolveNonFinite};
        }
        destination.unchecked(cell, 0U) = value;
        local_scale = std::max(local_scale, std::abs(value));
      }
    }
  }
  produced_local_scale = local_scale;
  return {};
}

Status add_scaled(FieldView destination, double scale,
                  ConstFieldView source) noexcept {
  for (std::int32_t z = 0; z < destination.interior.z; ++z) {
    for (std::int32_t y = 0; y < destination.interior.y; ++y) {
      for (std::int32_t x = 0; x < destination.interior.x; ++x) {
        const double value = destination.unchecked({x, y, z}, 0U) +
                             scale * source.unchecked({x, y, z}, 0U);
        if (!std::isfinite(value)) {
          return {StatusCode::numerical_failure, kLinearSolveNonFinite};
        }
        destination.unchecked({x, y, z}, 0U) = value;
      }
    }
  }
  return {};
}

Status update_bicgstab_search_direction(FieldView direction,
                                        ConstFieldView residual,
                                        ConstFieldView image, double beta,
                                        double omega) noexcept {
  for (std::int32_t z = 0; z < direction.interior.z; ++z) {
    for (std::int32_t y = 0; y < direction.interior.y; ++y) {
      for (std::int32_t x = 0; x < direction.interior.x; ++x) {
        const Int3 cell{x, y, z};
        const double intermediate =
            1.0 * direction.unchecked(cell, 0U) +
            -omega * image.unchecked(cell, 0U);
        if (!std::isfinite(intermediate)) {
          return {StatusCode::numerical_failure, kLinearSolveNonFinite};
        }
        const double value = 1.0 * residual.unchecked(cell, 0U) +
                             beta * intermediate;
        if (!std::isfinite(value)) {
          return {StatusCode::numerical_failure, kLinearSolveNonFinite};
        }
        direction.unchecked(cell, 0U) = value;
      }
    }
  }
  return {};
}

Status update_bicgstab_solution_residual(
    FieldView solution, ConstFieldView first_correction,
    ConstFieldView second_correction, FieldView residual,
    ConstFieldView intermediate_residual, ConstFieldView image, double alpha,
    double omega, double* produced_local_scale = nullptr,
    const ConstFieldView* rho_left = nullptr, double* produced_local_rho = nullptr,
    bool* rho_available = nullptr) noexcept {
  double local_scale = 0.0;
  double rho_sum = 0.0;
  double rho_correction = 0.0;
  bool local_rho_available = rho_left != nullptr;
  if (rho_available != nullptr) {
    *rho_available = false;
  }
  for (std::int32_t z = 0; z < solution.interior.z; ++z) {
    for (std::int32_t y = 0; y < solution.interior.y; ++y) {
      for (std::int32_t x = 0; x < solution.interior.x; ++x) {
        const Int3 cell{x, y, z};
        const double solution_intermediate =
            solution.unchecked(cell, 0U) +
            alpha * first_correction.unchecked(cell, 0U);
        if (!std::isfinite(solution_intermediate)) {
          return {StatusCode::numerical_failure, kLinearSolveNonFinite};
        }
        const double next_solution =
            solution_intermediate +
            omega * second_correction.unchecked(cell, 0U);
        if (!std::isfinite(next_solution)) {
          return {StatusCode::numerical_failure, kLinearSolveNonFinite};
        }
        const double next_residual =
            1.0 * intermediate_residual.unchecked(cell, 0U) +
            -omega * image.unchecked(cell, 0U);
        if (!std::isfinite(next_residual)) {
          return {StatusCode::numerical_failure, kLinearSolveNonFinite};
        }
        if (local_rho_available) {
          const double product =
              rho_left->unchecked(cell, 0U) * next_residual;
          if (!std::isfinite(product)) {
            local_rho_available = false;
          } else {
            const double next = rho_sum + product;
            if (!std::isfinite(next)) {
              local_rho_available = false;
            } else {
              const double correction =
                  std::abs(rho_sum) >= std::abs(product)
                      ? (rho_sum - next) + product
                      : (product - next) + rho_sum;
              if (!std::isfinite(correction)) {
                local_rho_available = false;
              } else {
                rho_correction += correction;
                if (!std::isfinite(rho_correction)) {
                  local_rho_available = false;
                }
                rho_sum = next;
              }
            }
          }
        }
        solution.unchecked(cell, 0U) = next_solution;
        residual.unchecked(cell, 0U) = next_residual;
        local_scale = std::max(local_scale, std::abs(next_residual));
      }
    }
  }
  if (produced_local_scale != nullptr) {
    *produced_local_scale = local_scale;
  }
  if (local_rho_available && produced_local_rho != nullptr) {
    const double local_rho = rho_sum + rho_correction;
    if (std::isfinite(local_rho)) {
      *produced_local_rho = local_rho;
      if (rho_available != nullptr) {
        *rho_available = true;
      }
    }
  }
  return {};
}

Status add_scaled_basis(FieldView destination, const ConstFieldView* sources,
                        const double* scales, std::size_t count) noexcept {
  if (sources == nullptr || scales == nullptr || count == 0U ||
      count > kMaximumFgmresBasisUpdateCount) {
    return {StatusCode::invalid_plan, kLinearSolveWorkspace};
  }

  std::array<double, kFgmresBasisUpdateStripCells> strip;
  for (std::int32_t z = 0; z < destination.interior.z; ++z) {
    for (std::int32_t y = 0; y < destination.interior.y; ++y) {
      for (std::int32_t x_begin = 0; x_begin < destination.interior.x;
           x_begin += static_cast<std::int32_t>(
               kFgmresBasisUpdateStripCells)) {
        const std::size_t width = std::min(
            kFgmresBasisUpdateStripCells,
            static_cast<std::size_t>(destination.interior.x - x_begin));
        double* const destination_row =
            destination.base + static_cast<std::size_t>(z) *
                                   destination.stride_z +
            static_cast<std::size_t>(y) * destination.stride_y +
            static_cast<std::size_t>(x_begin);
        std::copy_n(destination_row, width, strip.data());

        for (std::size_t row = 0U; row < count; ++row) {
          const double* const source_row =
              sources[row].base + static_cast<std::size_t>(z) *
                                      sources[row].stride_z +
              static_cast<std::size_t>(y) * sources[row].stride_y +
              static_cast<std::size_t>(x_begin);
          const double scale = scales[row];
          std::uint32_t finite = 1U;
#if defined(__clang__)
#pragma clang loop vectorize(enable)
#endif
          for (std::size_t lane = 0U; lane < width; ++lane) {
            const double value =
                strip[lane] + scale * source_row[lane];
            strip[lane] = value;
            finite &= static_cast<std::uint32_t>(std::isfinite(value));
          }
          if (finite == 0U) {
            // This private strip is intentionally not published.  Earlier
            // completed strips may differ from the cell-major failure path,
            // but the caller reaches finish_failure before publishing a
            // solution.
            return {StatusCode::numerical_failure, kLinearSolveNonFinite};
          }
        }
        std::copy_n(strip.data(), width, destination_row);
      }
    }
  }
  return {};
}

Status single_reduction_fgmres_norm(const double* projections,
                                    std::size_t count, double bar_norm,
                                    double& t, double& next_norm) noexcept {
  if (projections == nullptr || count == 0U ||
      count > kMaximumFgmresBasisUpdateCount || !std::isfinite(bar_norm)) {
    t = std::numeric_limits<double>::quiet_NaN();
    next_norm = std::numeric_limits<double>::quiet_NaN();
    return {StatusCode::numerical_failure, kLinearSolveNonFinite};
  }
  t = bar_norm;
  for (std::size_t row = 0U; row < count; ++row) {
    const double projection = projections[row];
    if (!std::isfinite(projection)) {
      t = std::numeric_limits<double>::quiet_NaN();
      next_norm = std::numeric_limits<double>::quiet_NaN();
      return {StatusCode::numerical_failure, kLinearSolveNonFinite};
    }
    const double square = projection * projection;
    if (!std::isfinite(square)) {
      t = std::numeric_limits<double>::quiet_NaN();
      next_norm = std::numeric_limits<double>::quiet_NaN();
      return {StatusCode::numerical_failure, kLinearSolveNonFinite};
    }
    t = t - square;
    if (!std::isfinite(t)) {
      t = std::numeric_limits<double>::quiet_NaN();
      next_norm = std::numeric_limits<double>::quiet_NaN();
      return {StatusCode::numerical_failure, kLinearSolveNonFinite};
    }
  }
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  if (g_forced_single_reduction_breakdown_count == count) {
    t = -1.0;
  }
#endif
  if (t < 0.0) {
    next_norm = 0.0;
    return {StatusCode::numerical_failure, kLinearSolveBreakdown};
  }
  next_norm = std::sqrt(t);
  return std::isfinite(next_norm)
             ? Status{}
             : Status{StatusCode::numerical_failure,
                      kLinearSolveNonFinite};
}

Status scale_copy(FieldView destination, ConstFieldView source,
                  double scale) noexcept {
  for (std::int32_t z = 0; z < destination.interior.z; ++z) {
    for (std::int32_t y = 0; y < destination.interior.y; ++y) {
      for (std::int32_t x = 0; x < destination.interior.x; ++x) {
        const double value = scale * source.unchecked({x, y, z}, 0U);
        if (!std::isfinite(value)) {
          return {StatusCode::numerical_failure, kLinearSolveNonFinite};
        }
        destination.unchecked({x, y, z}, 0U) = value;
      }
    }
  }
  return {};
}

Status merge_status(Status primary, Status secondary) noexcept {
  return primary ? secondary : primary;
}

Status global_dot(ReductionEngine& reductions, ConstFieldView left,
                  ConstFieldView right, Status pending,
                  double& result) noexcept {
  double local = 0.0;
  const Status arithmetic = local_dot(left, right, local);
  pending = merge_status(pending, arithmetic);
  double global = 0.0;
  const Status reduced = reductions.checked_sum(
      {&local, 1U}, {&global, 1U}, pending);
  if (!reduced) {
    return reduced;
  }
  if (!std::isfinite(global)) {
    return {StatusCode::numerical_failure, kLinearSolveNonFinite};
  }
  result = global;
  return {};
}

Status global_norm(ReductionEngine& reductions, ConstFieldView field,
                   Status pending, double& norm) noexcept {
  double local_scale = 0.0;
  Status arithmetic{};
  for (std::int32_t z = 0; z < field.interior.z; ++z) {
    for (std::int32_t y = 0; y < field.interior.y; ++y) {
      for (std::int32_t x = 0; x < field.interior.x; ++x) {
        const double value = field.unchecked({x, y, z}, 0U);
        if (!std::isfinite(value)) {
          arithmetic = {StatusCode::numerical_failure,
                        kLinearSolveNonFinite};
        } else {
          local_scale = std::max(local_scale, std::abs(value));
        }
      }
    }
  }
  pending = merge_status(pending, arithmetic);
  double global_scale = 0.0;
  Status status = reductions.checked_max(
      {&local_scale, 1U}, {&global_scale, 1U}, pending);
  if (!status) {
    return status;
  }
  if (global_scale == 0.0) {
    norm = 0.0;
    return {};
  }
  double local_scaled_squares = 0.0;
  double correction = 0.0;
  for (std::int32_t z = 0; z < field.interior.z; ++z) {
    for (std::int32_t y = 0; y < field.interior.y; ++y) {
      for (std::int32_t x = 0; x < field.interior.x; ++x) {
        const double scaled = field.unchecked({x, y, z}, 0U) / global_scale;
        const double square = scaled * scaled;
        const double next = local_scaled_squares + square;
        if (std::abs(local_scaled_squares) >= std::abs(square)) {
          correction += (local_scaled_squares - next) + square;
        } else {
          correction += (square - next) + local_scaled_squares;
        }
        local_scaled_squares = next;
      }
    }
  }
  local_scaled_squares += correction;
  double global_scaled_squares = 0.0;
  status = reductions.checked_sum(
      {&local_scaled_squares, 1U}, {&global_scaled_squares, 1U});
  if (!status) {
    return status;
  }
  if (global_scaled_squares < 0.0 ||
      !std::isfinite(global_scaled_squares)) {
    return {StatusCode::numerical_failure, kLinearSolveNonFinite};
  }
  norm = global_scale * std::sqrt(global_scaled_squares);
  return std::isfinite(norm)
             ? Status{}
             : Status{StatusCode::numerical_failure,
                      kLinearSolveNonFinite};
}

Status stable_norm_from_global_scale(ConstFieldView field, double global_scale,
                                     double& norm) noexcept {
  if (global_scale == 0.0) {
    norm = 0.0;
    return {};
  }
  double local_scaled_squares = 0.0;
  double correction = 0.0;
  for (std::int32_t z = 0; z < field.interior.z; ++z) {
    for (std::int32_t y = 0; y < field.interior.y; ++y) {
      for (std::int32_t x = 0; x < field.interior.x; ++x) {
        const double scaled = field.unchecked({x, y, z}, 0U) / global_scale;
        const double square = scaled * scaled;
        const double next = local_scaled_squares + square;
        if (std::abs(local_scaled_squares) >= std::abs(square)) {
          correction += (local_scaled_squares - next) + square;
        } else {
          correction += (square - next) + local_scaled_squares;
        }
        local_scaled_squares = next;
      }
    }
  }
  local_scaled_squares += correction;
  norm = global_scale * std::sqrt(local_scaled_squares);
  return std::isfinite(norm)
             ? Status{}
             : Status{StatusCode::numerical_failure,
                      kLinearSolveNonFinite};
}

Status global_norm_from_local_scale(ReductionEngine& reductions,
                                    ConstFieldView field, double local_scale,
                                    Status pending, double& norm) noexcept {
  double global_scale = 0.0;
  const Status status = reductions.checked_max(
      {&local_scale, 1U}, {&global_scale, 1U}, pending);
  if (!status) {
    return status;
  }
  if (global_scale == 0.0) {
    norm = 0.0;
    return {};
  }
  double local_scaled_squares = 0.0;
  double correction = 0.0;
  for (std::int32_t z = 0; z < field.interior.z; ++z) {
    for (std::int32_t y = 0; y < field.interior.y; ++y) {
      for (std::int32_t x = 0; x < field.interior.x; ++x) {
        const double scaled = field.unchecked({x, y, z}, 0U) / global_scale;
        const double square = scaled * scaled;
        const double next = local_scaled_squares + square;
        if (std::abs(local_scaled_squares) >= std::abs(square)) {
          correction += (local_scaled_squares - next) + square;
        } else {
          correction += (square - next) + local_scaled_squares;
        }
        local_scaled_squares = next;
      }
    }
  }
  local_scaled_squares += correction;
  double global_scaled_squares = 0.0;
  const Status reduced = reductions.checked_sum(
      {&local_scaled_squares, 1U}, {&global_scaled_squares, 1U});
  if (!reduced) {
    return reduced;
  }
  if (global_scaled_squares < 0.0 ||
      !std::isfinite(global_scaled_squares)) {
    return {StatusCode::numerical_failure, kLinearSolveNonFinite};
  }
  norm = global_scale * std::sqrt(global_scaled_squares);
  return std::isfinite(norm)
             ? Status{}
             : Status{StatusCode::numerical_failure,
                      kLinearSolveNonFinite};
}

Status revise(SolverWorkspace& workspace, std::uint8_t slot,
              FieldView& view, Status pending = {}) noexcept {
  const Status revised = workspace.revise_vector(slot);
  if (revised) {
    const FieldView refreshed = workspace.vector(slot, view.interior);
    if (refreshed.base == nullptr) {
      return {StatusCode::invalid_plan, kLinearSolveWorkspace};
    }
    view = refreshed;
  }
  return merge_status(pending, revised);
}

struct OperatorApplyStatus {
  Status status{};
  bool collective_failure{};
};

OperatorApplyStatus apply_operator(const LinearOperator& linear_operator,
                                   FieldView input, FieldView& output,
                                   std::uint8_t output_slot,
                                   SolverWorkspace& workspace,
                                   LinearSolveResult& result) noexcept {
  const Status applied = linear_operator.apply(input, output);
  ++result.operator_applies;
  if (!applied) {
    const LinearOperatorFailureProvenance provenance =
        linear_operator.failure_provenance();
    if (matching_collective_operator_failure(applied, provenance)) {
      result.lowest_failing_rank = provenance.lowest_failing_rank;
      return {applied, true};
    }
  }
  return {revise(workspace, output_slot, output, applied), false};
}

Status apply_preconditioner(LinearPreconditioner& preconditioner,
                            ConstFieldView input, FieldView& output,
                            std::uint8_t output_slot,
                            std::uint32_t iteration,
                            SolverWorkspace& workspace,
                            LinearSolveResult& result,
                            const LinearPreconditionerBatchTicket* ticket =
                                nullptr) noexcept {
  const Status applied =
      ticket == nullptr
          ? preconditioner.apply(input, output, iteration)
          : preconditioner.apply_prepared(input, output, iteration, *ticket);
  ++result.preconditioner_applies;
  return revise(workspace, output_slot, output, applied);
}

void update_reduction_count(LinearSolveResult& result,
                            const ReductionEngine& reductions,
                            std::uint64_t initial_calls) noexcept {
  const std::uint64_t final_calls = reductions.counters().calls;
  const std::uint64_t total_calls =
      final_calls >= initial_calls ? final_calls - initial_calls : 0U;
  const std::uint64_t excluded_calls =
      result.recycle_reduction_calls <=
              std::numeric_limits<std::uint64_t>::max() -
                  result.recycle_capture_reduction_calls
          ? result.recycle_reduction_calls +
                result.recycle_capture_reduction_calls
          : total_calls;
  result.reduction_calls =
      total_calls >= excluded_calls
          ? total_calls - excluded_calls
          : 0U;
  result.lowest_failing_rank = reductions.lowest_failing_rank();
}

Status validate_iteration_account(const ResourceCounters* resources,
                                  std::uint32_t iterations) noexcept {
  if (resources == nullptr || iterations == 0U) {
    return {};
  }
  return iterations >
                 std::numeric_limits<std::uint64_t>::max() -
                     resources->linear_iterations
             ? Status{StatusCode::invalid_plan, kLinearSolveWorkspace}
             : Status{};
}

void publish_iteration_account(ResourceCounters* resources,
                               std::uint32_t iterations) noexcept {
  if (resources != nullptr) {
    resources->linear_iterations += iterations;
  }
}

Status collective_account(ResourceCounters* resources,
                          std::uint32_t iterations,
                          ReductionEngine& reductions) noexcept {
  const Status consensus = reductions.consensus(
      validate_iteration_account(resources, iterations));
  if (consensus) {
    publish_iteration_account(resources, iterations);
  }
  return consensus;
}

LinearSolveResult invalid_result(Status status,
                                 LinearTermination termination,
                                 ReductionEngine& reductions,
                                 std::uint64_t initial_calls) noexcept {
  LinearSolveResult result;
  result.status = status;
  result.termination = termination;
  update_reduction_count(result, reductions, initial_calls);
  return result;
}

Status validate_common(LinearAlgorithm algorithm,
                       const LinearOperator& linear_operator,
                       const LinearPreconditioner& preconditioner,
                       const LinearSolveInvocation& invocation,
                       const SolverWorkspace& workspace,
                       const ReductionEngine& reductions) noexcept {
  const LinearOperatorCertificate op = linear_operator.certificate();
  const LinearPreconditionerCertificate pc = preconditioner.certificate();
  const LinearWorkspaceRequirements& requirements = workspace.requirements();
  const FieldView workspace_vector = workspace.vector(
      0U, requirements.maximum_shape);
  const LinearSolveControl control = invocation.control;
  if (!valid_scalar_view(invocation.rhs) ||
      !valid_scalar_view(invocation.solution) ||
      !same_shape(invocation.rhs.interior, invocation.solution.interior) ||
      views_overlap(invocation.rhs, invocation.solution) ||
      views_overlap(invocation.rhs, workspace_vector) ||
      views_overlap(invocation.solution, workspace_vector) ||
      workspace.overlaps_storage(invocation.rhs) ||
      workspace.overlaps_storage(invocation.solution) ||
      !valid_identity(invocation.expected_identity) ||
      invocation.expected_identity.workspace != workspace.fingerprint() ||
      !same_identity(op.identity, invocation.expected_identity) ||
      !same_identity(pc.identity, invocation.expected_identity) ||
      op.collective_fingerprint == 0U || pc.collective_fingerprint == 0U ||
      !valid_preconditioner_status_scope(pc.status_scope) ||
      !valid_preconditioner_apply_lifecycle(pc.apply_lifecycle) ||
      (pc.apply_lifecycle ==
           LinearPreconditionerApplyLifecycle::prepared_batch &&
       (algorithm != LinearAlgorithm::fgmres &&
        algorithm != LinearAlgorithm::bicgstab)) ||
      (pc.apply_lifecycle ==
           LinearPreconditionerApplyLifecycle::prepared_batch &&
       pc.status_scope != LinearPreconditionerStatusScope::collective) ||
      !same_shape(op.local_shape, invocation.rhs.interior) ||
      requirements.algorithm != algorithm ||
      !shape_contains(requirements.maximum_shape, invocation.rhs.interior) ||
      workspace.fingerprint() == 0U || reductions.capacity() == 0U ||
      reductions.capacity() < requirements.reduction_capacity ||
      reductions.mode() != requirements.reduction_mode ||
      !std::isfinite(control.absolute_tolerance) ||
      !std::isfinite(control.relative_tolerance) ||
      control.absolute_tolerance < 0.0 || control.relative_tolerance < 0.0 ||
      (control.absolute_tolerance == 0.0 &&
       control.relative_tolerance == 0.0) ||
      control.maximum_iterations == 0U ||
      control.true_residual_interval == 0U ||
      (invocation.convergence_audit != nullptr &&
       ((algorithm != LinearAlgorithm::fgmres &&
         algorithm != LinearAlgorithm::bicgstab) ||
        invocation.convergence_audit->certificate().collective_fingerprint ==
            0U))) {
    return {StatusCode::invalid_plan, kLinearSolvePlan};
  }
  if (algorithm == LinearAlgorithm::pcg) {
    if (control.restart != 0U || op.operator_class != LinearOperatorClass::spd ||
        pc.preconditioner_class != LinearPreconditionerClass::fixed_spd) {
      return {StatusCode::invalid_plan, kLinearSolvePlan};
    }
  } else if (algorithm == LinearAlgorithm::fgmres) {
    if (control.restart == 0U ||
        control.restart > requirements.maximum_restart) {
      return {StatusCode::invalid_plan, kLinearSolvePlan};
    }
  } else if (control.restart != 0U ||
             pc.preconditioner_class !=
                 LinearPreconditionerClass::fixed_general) {
    return {StatusCode::invalid_plan, kLinearSolvePlan};
  }
  return {};
}

Status prepare_solve(LinearAlgorithm algorithm,
                     const LinearOperator& linear_operator,
                     LinearPreconditioner& preconditioner,
                     const LinearSolveInvocation& invocation,
                     const SolverWorkspace& workspace,
                     ReductionEngine& reductions,
                     LinearPreconditionerStatusScope* status_scope = nullptr,
                     LinearPreconditionerBatchTicket* batch_ticket = nullptr,
                     LinearPreconditionerApplyLifecycle* apply_lifecycle =
                         nullptr)
    noexcept {
  const LinearOperatorCertificate op = linear_operator.certificate();
  const LinearPreconditionerCertificate pc = preconditioner.certificate();
  if (status_scope != nullptr) {
    *status_scope = pc.status_scope;
  }
  if (apply_lifecycle != nullptr) {
    *apply_lifecycle = pc.apply_lifecycle;
  }
  const Status local_valid = validate_common(
      algorithm, linear_operator, preconditioner, invocation, workspace,
      reductions);
  const PlanFingerprint contract = solve_contract_fingerprint(
      algorithm, op, pc, invocation, workspace.requirements());
  // Agree on the callback contract before entering cold prepare.  This keeps
  // lifecycle mismatches and invalid values ahead of every callback while
  // retaining the existing two solve-preparation collectives.
  const Status contract_status = reductions.consensus_contract(contract);
  if (!contract_status) {
    return contract_status;
  }

  Status local = local_valid;
  if (local && pc.apply_lifecycle ==
                   LinearPreconditionerApplyLifecycle::prepared_batch) {
    if (batch_ticket == nullptr || (algorithm != LinearAlgorithm::fgmres &&
                                    algorithm != LinearAlgorithm::bicgstab)) {
      local = {StatusCode::invalid_plan, kLinearSolvePlan};
    } else {
      const bool bicgstab = algorithm == LinearAlgorithm::bicgstab;
      const std::size_t input_begin = bicgstab ? 3U : 1U;
      const std::size_t output_begin =
          bicgstab ? 5U
                   : static_cast<std::size_t>(invocation.control.restart) + 2U;
      const std::size_t slot_count =
          bicgstab ? 2U : static_cast<std::size_t>(invocation.control.restart);
      const std::size_t maximum_applications =
          bicgstab ? static_cast<std::size_t>(
                         invocation.control.maximum_iterations) *
                         2U
                   : invocation.control.maximum_iterations;
      if (slot_count == 0U || maximum_applications == 0U ||
          input_begin > std::numeric_limits<std::uint8_t>::max() ||
          output_begin > std::numeric_limits<std::uint8_t>::max() ||
          maximum_applications > std::numeric_limits<std::uint32_t>::max()) {
        local = {StatusCode::invalid_plan, kLinearSolveWorkspace};
      } else {
        const LinearPreconditionerBatchDescriptor descriptor{
            &workspace,
            invocation.rhs.interior,
            static_cast<std::uint8_t>(input_begin),
            static_cast<std::uint8_t>(output_begin),
            static_cast<std::uint8_t>(slot_count),
            static_cast<std::uint32_t>(maximum_applications)};
        *batch_ticket = {};
        local = preconditioner.prepare_batch(descriptor, *batch_ticket);
      }
    }
  }
  // Cold preparation is local/no-MPI and contributes to the second
  // solve-preparation collective; the contract collective above is what
  // prevents a lifecycle mismatch from entering any callback.
  const Status agreed = reductions.consensus(local);
  if (!agreed && batch_ticket != nullptr) {
    *batch_ticket = {};
  }
  return agreed;
}

Status publish_solution(ConstFieldView source, FieldView destination,
                        ResourceCounters* resources,
                        std::uint32_t iterations,
                        ReductionEngine& reductions) noexcept {
  const Status accounted =
      collective_account(resources, iterations, reductions);
  if (!accounted) {
    return accounted;
  }
  return copy_field(source, destination);
}

Status compute_true_residual(const LinearOperator& linear_operator,
                             ConstFieldView rhs, FieldView solution,
                             FieldView& residual, std::uint8_t residual_slot,
                             FieldView& operator_output,
                             std::uint8_t operator_slot,
                             SolverWorkspace& workspace,
                             ReductionEngine& reductions,
                             LinearSolveResult& result,
                             double& norm) noexcept {
  const OperatorApplyStatus applied = apply_operator(
      linear_operator, solution, operator_output, operator_slot, workspace,
      result);
  if (applied.collective_failure) {
    return applied.status;
  }
  Status pending = applied.status;
  const Status formed = axpby(residual, 1.0, rhs, -1.0,
                              as_const(operator_output));
  pending = revise(workspace, residual_slot, residual,
                   merge_status(pending, formed));
  return global_norm(reductions, as_const(residual), pending, norm);
}

Status compute_true_residual_norm_in_place(
    const LinearOperator& linear_operator, ConstFieldView rhs,
    FieldView solution, FieldView& operator_output,
    std::uint8_t operator_slot, SolverWorkspace& workspace,
    ReductionEngine& reductions, LinearSolveResult& result,
    double& norm) noexcept {
  const OperatorApplyStatus applied = apply_operator(
      linear_operator, solution, operator_output, operator_slot, workspace,
      result);
  if (applied.collective_failure) {
    return applied.status;
  }
  Status pending = applied.status;
  const Status formed = axpby(operator_output, 1.0, rhs, -1.0,
                              as_const(operator_output));
  pending = revise(workspace, operator_slot, operator_output,
                   merge_status(pending, formed));
  return global_norm(reductions, as_const(operator_output), pending, norm);
}

double convergence_limit(const LinearSolveInvocation& invocation,
                         double rhs_norm) noexcept {
  return std::max(invocation.control.absolute_tolerance,
                  invocation.control.relative_tolerance * rhs_norm);
}

Status audit_convergence(const LinearSolveInvocation& invocation,
                         ConstFieldView solution,
                         ConstFieldView true_residual,
                         ReductionEngine& reductions,
                         LinearSolveResult& result,
                         bool& accepted) noexcept {
  accepted = true;
  if (invocation.convergence_audit == nullptr) return {};
  if (result.convergence_audits ==
      std::numeric_limits<std::uint64_t>::max()) {
    return {StatusCode::invalid_plan, kLinearSolveConvergenceAudit};
  }
  LinearConvergenceAuditResult audit;
  const Status status = invocation.convergence_audit->evaluate(
      solution, true_residual, reductions, audit);
  if (!status) return status;
  const LinearConvergenceFailureProvenance& provenance =
      audit.failure_provenance;
  const bool valid_failure_provenance =
      !provenance.valid ||
      (!audit.accepted && audit.terminal_rejection &&
       provenance.global_index.x >= 0 && provenance.global_index.y >= 0 &&
       provenance.global_index.z >= 0 && provenance.owner_rank >= 0 &&
       provenance.pressure_activity <= 1U &&
       std::isfinite(provenance.predecessor_application_scale) &&
       provenance.predecessor_application_scale > 0.0 &&
       provenance.predecessor_application_scale <= 1.0 &&
       std::isfinite(provenance.predecessor_maximum_depletion) &&
       provenance.predecessor_maximum_depletion >= 0.0 &&
       std::isfinite(provenance.density) && provenance.density > 0.0 &&
       std::isfinite(provenance.psi) && provenance.psi >= 0.0 &&
       std::isfinite(provenance.raw_correction) &&
       std::isfinite(provenance.depletion) && provenance.depletion > 0.0 &&
       provenance.depletion == audit.maximum_depletion &&
       std::isfinite(provenance.bdf_storage) &&
       std::isfinite(provenance.flux_divergence) &&
       std::isfinite(provenance.reconstructed_rhs) &&
       std::isfinite(provenance.sealed_rhs) &&
       std::isfinite(provenance.rhs_absolute_mismatch) &&
       provenance.rhs_absolute_mismatch >= 0.0 &&
       std::isfinite(provenance.rhs_relative_mismatch) &&
       provenance.rhs_relative_mismatch >= 0.0);
  if (!std::isfinite(audit.metric) || audit.metric < 0.0 ||
      !std::isfinite(audit.limit) || !(audit.limit > 0.0) ||
      !std::isfinite(audit.application_scale) ||
      !(audit.application_scale > 0.0) || audit.application_scale > 1.0 ||
      !std::isfinite(audit.unscaled_metric) || audit.unscaled_metric < 0.0 ||
      !std::isfinite(audit.maximum_depletion) ||
      audit.maximum_depletion < 0.0 ||
      !std::isfinite(audit.operator_parity_error) ||
      audit.operator_parity_error < 0.0 ||
      !valid_failure_provenance ||
      (audit.accepted && audit.terminal_rejection) ||
      audit.accepted != (audit.metric <= audit.limit)) {
    return {StatusCode::invalid_plan, kLinearSolveConvergenceAudit};
  }
  ++result.convergence_audits;
  result.final_convergence_metric = audit.metric;
  result.convergence_limit = audit.limit;
  result.convergence_application_scale = audit.application_scale;
  result.convergence_unscaled_metric = audit.unscaled_metric;
  result.convergence_maximum_depletion = audit.maximum_depletion;
  result.convergence_operator_parity_error = audit.operator_parity_error;
  if (provenance.valid)
    result.convergence_failure_provenance = provenance;
  accepted = audit.accepted;
  if (!accepted) {
    if (result.convergence_rejections ==
        std::numeric_limits<std::uint64_t>::max()) {
      return {StatusCode::invalid_plan, kLinearSolveConvergenceAudit};
    }
    ++result.convergence_rejections;
  }
  if (audit.terminal_rejection) {
    return {StatusCode::rejected_step, kLinearSolveConvergenceAudit};
  }
  return {};
}

LinearSolveResult finish_failure(LinearSolveResult result, Status status,
                                 LinearTermination termination,
                                 ResourceCounters* resources,
                                 ReductionEngine& reductions,
                                 std::uint64_t initial_calls) noexcept {
  result.status = status;
  result.termination = termination;
  const int original_failing_rank =
      result.lowest_failing_rank >= 0
          ? result.lowest_failing_rank
          : reductions.lowest_failing_rank();
  if (resources != nullptr && result.iterations != 0U) {
    const Status accounted = collective_account(
        resources, result.iterations, reductions);
    if (!accounted) {
      result.status = accounted;
      result.termination = LinearTermination::invalid_plan;
    }
  }
  update_reduction_count(result, reductions, initial_calls);
  if (original_failing_rank >= 0) {
    result.lowest_failing_rank = original_failing_rank;
  }
  return result;
}

LinearSolveResult finish_success(LinearSolveResult result,
                                 ConstFieldView solution,
                                 const LinearSolveInvocation& invocation,
                                 ResourceCounters* resources,
                                 ReductionEngine& reductions,
                                 std::uint64_t initial_calls) noexcept {
  const Status published = publish_solution(solution, invocation.solution,
                                            resources, result.iterations,
                                            reductions);
  if (!published) {
    return finish_failure(result, published, LinearTermination::invalid_plan,
                          nullptr, reductions, initial_calls);
  }
  result.status = {};
  result.termination = LinearTermination::converged;
  update_reduction_count(result, reductions, initial_calls);
  return result;
}

}  // namespace

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
namespace detail {

Status fused_krylov_basis_update_for_test(
    FieldView destination, const ConstFieldView* sources,
    const double* scales, std::size_t count) noexcept {
  return add_scaled_basis(destination, sources, scales, count);
}

bool krylov_basis_update_inputs_disjoint_for_test(
    FieldView destination, const ConstFieldView* sources,
    std::size_t count) noexcept {
  if (destination.base == nullptr || sources == nullptr || count == 0U ||
      count > kMaximumFgmresBasisUpdateCount) {
    return false;
  }
  for (std::size_t row = 0U; row < count; ++row) {
    if (sources[row].base == nullptr ||
        views_overlap(as_const(destination), sources[row])) {
      return false;
    }
  }
  return true;
}

Status single_reduction_fgmres_norm_for_test(
    Span<const double> projections, double bar_norm, double& t,
    double& next_norm) noexcept {
  return single_reduction_fgmres_norm(projections.data, projections.size,
                                      bar_norm, t, next_norm);
}

void force_single_reduction_fgmres_breakdown_for_test(
    std::size_t projection_count) noexcept {
  g_forced_single_reduction_breakdown_count = projection_count;
}

}  // namespace detail
#endif

LinearSolveResult solve_pcg(const LinearOperator& linear_operator,
                            LinearPreconditioner& preconditioner,
                            const LinearSolveInvocation& invocation,
                            SolverWorkspace& workspace,
                            ReductionEngine& reductions,
                            ResourceCounters* resources) noexcept {
  const std::uint64_t initial_calls = reductions.counters().calls;
  const Status prepared = prepare_solve(LinearAlgorithm::pcg, linear_operator,
                                        preconditioner, invocation, workspace,
                                        reductions);
  if (!prepared) {
    return invalid_result(prepared, LinearTermination::invalid_plan,
                          reductions, initial_calls);
  }
  LinearSolveResult result;
  const Int3 shape = invocation.rhs.interior;
  FieldView x = workspace.vector(0U, shape);
  FieldView r = workspace.vector(1U, shape);
  FieldView z = workspace.vector(2U, shape);
  FieldView p = workspace.vector(3U, shape);
  FieldView ap = workspace.vector(4U, shape);

  double rhs_norm = 0.0;
  Status status = global_norm(reductions, invocation.rhs, {}, rhs_norm);
  if (!status) {
    return finish_failure(result, status, LinearTermination::non_finite,
                          resources, reductions, initial_calls);
  }
  if (rhs_norm == 0.0) {
    fill_field(x, 0.0);
    status = revise(workspace, 0U, x);
    status = reductions.consensus(status);
    status = merge_status(
        status, collective_account(resources, 0U, reductions));
    if (!status) {
      return finish_failure(result, status, LinearTermination::invalid_plan,
                            nullptr, reductions, initial_calls);
    }
    status = copy_field(as_const(x), invocation.solution);
    if (!status) {
      return finish_failure(result, status, LinearTermination::non_finite,
                            nullptr, reductions, initial_calls);
    }
    result.status = {};
    result.termination = LinearTermination::zero_rhs;
    result.initial_true_residual = 0.0;
    result.final_true_residual = 0.0;
    result.recursive_residual = 0.0;
    update_reduction_count(result, reductions, initial_calls);
    return result;
  }

  status = revise(workspace, 0U, x, copy_field(invocation.solution, x));
  status = reductions.consensus(status);
  if (!status) {
    return finish_failure(result, status, LinearTermination::invalid_plan,
                          resources, reductions, initial_calls);
  }
  status = compute_true_residual(linear_operator, invocation.rhs, x, r, 1U,
                                 ap, 4U, workspace, reductions, result,
                                 result.initial_true_residual);
  if (!status) {
    return finish_failure(result, status, LinearTermination::operator_failure,
                          resources, reductions, initial_calls);
  }
  result.final_true_residual = result.initial_true_residual;
  result.recursive_residual = result.initial_true_residual;
  const double tolerance = convergence_limit(invocation, rhs_norm);
  if (!std::isfinite(tolerance)) {
    return finish_failure(
        result, {StatusCode::numerical_failure, kLinearSolveNonFinite},
        LinearTermination::non_finite, resources, reductions, initial_calls);
  }
  if (result.initial_true_residual <= tolerance) {
    return finish_success(result, as_const(x), invocation, resources,
                          reductions, initial_calls);
  }

  status = apply_preconditioner(preconditioner, as_const(r), z, 2U, 0U,
                                workspace, result);
  double rho = 0.0;
  status = global_dot(reductions, as_const(r), as_const(z), status, rho);
  if (!status) {
    return finish_failure(result, status,
                          LinearTermination::preconditioner_failure,
                          resources, reductions, initial_calls);
  }
  if (!(rho > 0.0) || !std::isfinite(rho)) {
    return finish_failure(
        result, {StatusCode::numerical_failure, kLinearSolveBreakdown},
        LinearTermination::breakdown, resources, reductions, initial_calls);
  }
  status = revise(workspace, 3U, p, copy_field(as_const(z), p));
  status = reductions.consensus(status);
  if (!status) {
    return finish_failure(result, status, LinearTermination::non_finite,
                          resources, reductions, initial_calls);
  }

  while (result.iterations < invocation.control.maximum_iterations) {
    const OperatorApplyStatus operator_apply =
        apply_operator(linear_operator, p, ap, 4U, workspace, result);
    if (operator_apply.collective_failure) {
      return finish_failure(result, operator_apply.status,
                            LinearTermination::operator_failure, resources,
                            reductions, initial_calls);
    }
    double denominator = 0.0;
    status = global_dot(reductions, as_const(p), as_const(ap),
                        operator_apply.status, denominator);
    if (!status) {
      return finish_failure(result, status,
                            LinearTermination::operator_failure, resources,
                            reductions, initial_calls);
    }
    if (!(denominator > 0.0) || !std::isfinite(denominator)) {
      return finish_failure(
          result, {StatusCode::numerical_failure, kLinearSolveBreakdown},
          LinearTermination::breakdown, resources, reductions, initial_calls);
    }
    const double alpha = rho / denominator;
    if (!std::isfinite(alpha)) {
      return finish_failure(
          result, {StatusCode::numerical_failure, kLinearSolveNonFinite},
          LinearTermination::non_finite, resources, reductions, initial_calls);
    }
    status = revise(workspace, 0U, x, add_scaled(x, alpha, as_const(p)));
    status = revise(workspace, 1U, r,
                    merge_status(status,
                                 add_scaled(r, -alpha, as_const(ap))));
    ++result.iterations;
    double recursive = 0.0;
    status = global_norm(reductions, as_const(r), status, recursive);
    if (!status) {
      return finish_failure(result, status, LinearTermination::non_finite,
                            resources, reductions, initial_calls);
    }
    result.recursive_residual = recursive;
    const bool verify = recursive <= tolerance ||
                        result.iterations %
                                invocation.control.true_residual_interval ==
                            0U ||
                        result.iterations ==
                            invocation.control.maximum_iterations;
    if (verify) {
      status = compute_true_residual(
          linear_operator, invocation.rhs, x, r, 1U, ap, 4U, workspace,
          reductions, result, result.final_true_residual);
      if (!status) {
        return finish_failure(result, status,
                              LinearTermination::operator_failure, resources,
                              reductions, initial_calls);
      }
      result.recursive_residual = result.final_true_residual;
      if (result.final_true_residual <= tolerance) {
        return finish_success(result, as_const(x), invocation, resources,
                              reductions, initial_calls);
      }
      if (result.iterations == invocation.control.maximum_iterations) {
        return finish_failure(
            result, {StatusCode::rejected_step, kLinearSolveBreakdown},
            LinearTermination::maximum_iterations, resources, reductions,
            initial_calls);
      }
      status = apply_preconditioner(preconditioner, as_const(r), z, 2U,
                                    result.iterations, workspace, result);
      status = global_dot(reductions, as_const(r), as_const(z), status, rho);
      if (!status) {
        return finish_failure(result, status,
                              LinearTermination::preconditioner_failure,
                              resources, reductions, initial_calls);
      }
      if (!(rho > 0.0) || !std::isfinite(rho)) {
        return finish_failure(
            result, {StatusCode::numerical_failure, kLinearSolveBreakdown},
            LinearTermination::breakdown, resources, reductions,
            initial_calls);
      }
      status = revise(workspace, 3U, p, copy_field(as_const(z), p));
      status = reductions.consensus(status);
      if (!status) {
        return finish_failure(result, status, LinearTermination::non_finite,
                              resources, reductions, initial_calls);
      }
      continue;
    }

    status = apply_preconditioner(preconditioner, as_const(r), z, 2U,
                                  result.iterations, workspace, result);
    double rho_next = 0.0;
    status = global_dot(reductions, as_const(r), as_const(z), status,
                        rho_next);
    if (!status) {
      return finish_failure(result, status,
                            LinearTermination::preconditioner_failure,
                            resources, reductions, initial_calls);
    }
    if (!(rho_next > 0.0) || !std::isfinite(rho_next)) {
      return finish_failure(
          result, {StatusCode::numerical_failure, kLinearSolveBreakdown},
          LinearTermination::breakdown, resources, reductions, initial_calls);
    }
    const double beta = rho_next / rho;
    status = revise(workspace, 3U, p,
                    axpby(p, 1.0, as_const(z), beta, as_const(p)));
    status = reductions.consensus(status);
    if (!status) {
      return finish_failure(result, status, LinearTermination::non_finite,
                            resources, reductions, initial_calls);
    }
    rho = rho_next;
  }
  return finish_failure(result,
                        {StatusCode::rejected_step, kLinearSolveBreakdown},
                        LinearTermination::maximum_iterations, resources,
                        reductions, initial_calls);
}

LinearSolveResult solve_fgmres(const LinearOperator& linear_operator,
                               LinearPreconditioner& preconditioner,
                               const LinearSolveInvocation& invocation,
                               SolverWorkspace& workspace,
                               ReductionEngine& reductions,
                               ResourceCounters* resources) noexcept {
  const std::uint64_t initial_calls = reductions.counters().calls;
  LinearPreconditionerStatusScope preconditioner_status_scope{
      LinearPreconditionerStatusScope::rank_local};
  LinearPreconditionerApplyLifecycle preconditioner_apply_lifecycle{
      LinearPreconditionerApplyLifecycle::per_call_checked};
  LinearPreconditionerBatchTicket batch_ticket{};
  const Status prepared = prepare_solve(LinearAlgorithm::fgmres,
                                        linear_operator, preconditioner,
                                        invocation, workspace, reductions,
                                        &preconditioner_status_scope,
                                        &batch_ticket,
                                        &preconditioner_apply_lifecycle);
  if (!prepared) {
    return invalid_result(prepared, LinearTermination::invalid_plan,
                          reductions, initial_calls);
  }
  const bool prepared_batch =
      preconditioner_apply_lifecycle ==
      LinearPreconditionerApplyLifecycle::prepared_batch;
  LinearSolveResult result;
  const Int3 shape = invocation.rhs.interior;
  const std::uint32_t restart = invocation.control.restart;
  const std::uint8_t x_slot = 0U;
  const std::uint8_t v_begin = 1U;
  const std::uint8_t z_begin =
      static_cast<std::uint8_t>(v_begin + restart + 1U);
  const std::uint8_t ax_slot =
      static_cast<std::uint8_t>(2U * restart + 2U);
  FieldView x = workspace.vector(x_slot, shape);
  FieldView ax = workspace.vector(ax_slot, shape);

  double rhs_norm = 0.0;
  Status status = global_norm(reductions, invocation.rhs, {}, rhs_norm);
  if (!status) {
    return finish_failure(result, status, LinearTermination::non_finite,
                          resources, reductions, initial_calls);
  }
  if (rhs_norm == 0.0) {
    fill_field(x, 0.0);
    status = revise(workspace, x_slot, x);
    status = reductions.consensus(status);
    if (!status) {
      return finish_failure(result, status, LinearTermination::invalid_plan,
                            nullptr, reductions, initial_calls);
    }
    FieldView zero_residual = workspace.vector(v_begin, shape);
    fill_field(zero_residual, 0.0);
    status = revise(workspace, v_begin, zero_residual);
    status = reductions.consensus(status);
    bool audit_accepted = false;
    if (status)
      status = audit_convergence(invocation, as_const(x),
                                 as_const(zero_residual), reductions, result,
                                 audit_accepted);
    if (!status || !audit_accepted) {
      return finish_failure(
          result,
          status ? Status{StatusCode::rejected_step,
                          kLinearSolveConvergenceAudit}
                 : status,
          LinearTermination::convergence_audit_failure, resources,
          reductions, initial_calls);
    }
    const bool consume_projection = workspace.recycle_projection_pending();
    workspace.recycle_skip_projection(&result);
    if (consume_projection) {
      status = reductions.consensus({});
      if (!status) {
        return finish_failure(result, status,
                              LinearTermination::invalid_plan, nullptr,
                              reductions, initial_calls);
      }
    }
    status = collective_account(resources, 0U, reductions);
    if (!status) {
      return finish_failure(result, status, LinearTermination::invalid_plan,
                            nullptr, reductions, initial_calls);
    }
    status = copy_field(as_const(x), invocation.solution);
    if (!status) {
      return finish_failure(result, status, LinearTermination::non_finite,
                            nullptr, reductions, initial_calls);
    }
    result.status = {};
    result.termination = LinearTermination::zero_rhs;
    update_reduction_count(result, reductions, initial_calls);
    return result;
  }
  status = revise(workspace, x_slot, x,
                  copy_field(invocation.solution, x));
  status = reductions.consensus(status);
  if (!status) {
    return finish_failure(result, status, LinearTermination::invalid_plan,
                          resources, reductions, initial_calls);
  }
  FieldView v0 = workspace.vector(v_begin, shape);
  status = compute_true_residual(linear_operator, invocation.rhs, x, v0,
                                 v_begin, ax, ax_slot, workspace, reductions,
                                 result, result.initial_true_residual);
  if (!status) {
    return finish_failure(result, status, LinearTermination::operator_failure,
                          resources, reductions, initial_calls);
  }
  Status initial_guess_selection_status;
  const bool reset_initial_guess =
      workspace.recycle_capture_active() &&
      result.initial_true_residual > rhs_norm;
  if (reset_initial_guess) {
    initial_guess_selection_status =
        reset_to_zero_initial_guess(invocation.rhs, x, v0);
    initial_guess_selection_status = revise(
        workspace, x_slot, x, initial_guess_selection_status);
    initial_guess_selection_status = revise(
        workspace, v_begin, v0, initial_guess_selection_status);
    result.initial_true_residual = rhs_norm;
  }
  result.final_true_residual = result.initial_true_residual;
  result.recursive_residual = result.initial_true_residual;
  const double tolerance = convergence_limit(invocation, rhs_norm);
  if (!std::isfinite(tolerance)) {
    if (reset_initial_guess) {
      const Status selected =
          reductions.consensus(initial_guess_selection_status);
      if (!selected) {
        return finish_failure(result, selected,
                              LinearTermination::invalid_plan, resources,
                              reductions, initial_calls);
      }
    }
    return finish_failure(
        result, {StatusCode::numerical_failure, kLinearSolveNonFinite},
        LinearTermination::non_finite, resources, reductions, initial_calls);
  }
  bool convergence_audit_deferred = false;
  if (result.initial_true_residual <= tolerance) {
    if (reset_initial_guess) {
      status = reductions.consensus(initial_guess_selection_status);
      initial_guess_selection_status = {};
      if (!status) {
        return finish_failure(result, status,
                              LinearTermination::invalid_plan, resources,
                              reductions, initial_calls);
      }
    }
    const bool consume_projection = workspace.recycle_projection_pending();
    workspace.recycle_skip_projection(&result);
    if (consume_projection) {
      status = reductions.consensus({});
      if (!status) {
        return finish_failure(result, status,
                              LinearTermination::invalid_plan, resources,
                              reductions, initial_calls);
      }
    }
    bool audit_accepted = false;
    status = audit_convergence(invocation, as_const(x), as_const(v0),
                               reductions, result, audit_accepted);
    if (!status) {
      return finish_failure(
          result, status, LinearTermination::convergence_audit_failure,
          resources, reductions, initial_calls);
    }
    if (audit_accepted) {
      return finish_success(result, as_const(x), invocation, resources,
                            reductions, initial_calls);
    }
    if (result.initial_true_residual == 0.0) {
      return finish_failure(
          result,
          {StatusCode::rejected_step, kLinearSolveConvergenceAudit},
          LinearTermination::convergence_audit_failure, resources,
          reductions, initial_calls);
    }
    convergence_audit_deferred = true;
  }

  auto apply_recycle_projection = [&](FieldView& current_x,
                                      FieldView& current_residual,
                                      double baseline_residual) noexcept
      -> Status {
    if (!workspace.recycle_projection_pending()) return {};
    const std::uint64_t projection_reduction_begin =
        reductions.counters().calls;
    const std::size_t offered = workspace.recycle_correction_count();
    std::uint64_t operator_applies = 0U;
    std::size_t retained = 0U;
    bool projection_attempted = offered != 0U;
    bool projection_accepted = false;
    double projected_residual = 0.0;
    auto finish_projection = [&](Status pending) noexcept {
      const std::uint64_t final_calls = reductions.counters().calls;
      const std::uint64_t reduction_calls =
          final_calls >= projection_reduction_begin
              ? final_calls - projection_reduction_begin
              : 0U;
      workspace.recycle_set_projection_result(
          result, offered, static_cast<std::uint64_t>(retained),
          operator_applies, reduction_calls, projection_attempted,
          projection_accepted, projected_residual);
      const int projection_failure_rank =
          result.lowest_failing_rank >= 0
              ? result.lowest_failing_rank
              : reductions.lowest_failing_rank();
      workspace.recycle_skip_projection();
      // Projection ownership is single-use on every rank, including failure
      // and fallback.  Seal that state transition collectively without
      // counting it as checked numerical work.  Preserve the provenance from
      // the packet/operator that actually failed before this final agreement
      // can overwrite ReductionEngine's last-failure slot.
      const Status consumed = reductions.consensus(pending);
      if (projection_failure_rank >= 0) {
        result.lowest_failing_rank = projection_failure_rank;
      }
      return consumed;
    };

    if (offered == 0U) return finish_projection({});

    std::array<FieldView, kLinearRecycleMaximumDirections> images{};
    std::array<FieldView, kLinearRecycleMaximumDirections> directions{};
    std::array<std::uint8_t, kLinearRecycleMaximumDirections> image_slots{};
    std::array<double, kLinearRecycleMaximumDirections> pre_qr_norms{};
    std::array<std::size_t, kLinearRecycleMaximumDirections>
        retained_indices{};
    LinearSolveResult projection_operator_result{};
    double maximum_pre_qr_norm = 0.0;
    for (std::size_t index = 0U; index < offered; ++index) {
      directions[index] = workspace.recycle_correction_storage(index, shape);
      if (directions[index].base == nullptr) {
        return finish_projection(
            {StatusCode::invalid_plan, kLinearSolveWorkspace});
      }
    }
    for (std::size_t index = 0U; index < offered; ++index) {
      image_slots[index] =
          static_cast<std::uint8_t>(v_begin + 1U + index);
      images[index] = workspace.vector(image_slots[index], shape);
      bool aliases_direction = false;
      for (std::size_t direction = 0U; direction < offered; ++direction) {
        aliases_direction =
            detail::field_views_overlap(as_const(directions[direction]),
                                        as_const(images[index])) ||
            aliases_direction;
      }
      if (aliases_direction) {
        image_slots[index] = workspace.recycle_snapshot_slot();
        images[index] = workspace.recycle_snapshot(shape);
      }
      bool disjoint = images[index].base != nullptr;
      for (std::size_t direction = 0U; direction < offered; ++direction) {
        disjoint =
            !detail::field_views_overlap(as_const(directions[direction]),
                                         as_const(images[index])) &&
            disjoint;
      }
      for (std::size_t prior = 0U; prior < index; ++prior) {
        disjoint =
            !detail::field_views_overlap(as_const(images[prior]),
                                         as_const(images[index])) &&
            disjoint;
      }
      if (!disjoint) {
        return finish_projection(
            {StatusCode::invalid_plan, kLinearSolveWorkspace});
      }
      const OperatorApplyStatus applied = apply_operator(
          linear_operator, directions[index], images[index],
          image_slots[index], workspace, projection_operator_result);
      ++operator_applies;
      if (applied.collective_failure) {
        result.lowest_failing_rank =
            projection_operator_result.lowest_failing_rank;
        return finish_projection(applied.status);
      }
      const Status norm_status = global_norm(
          reductions, as_const(images[index]), applied.status,
          pre_qr_norms[index]);
      if (!norm_status) return finish_projection(norm_status);
      maximum_pre_qr_norm =
          std::max(maximum_pre_qr_norm, pre_qr_norms[index]);
    }

    if (!(maximum_pre_qr_norm > 0.0) ||
        !std::isfinite(maximum_pre_qr_norm)) {
      return finish_projection({});
    }
    const double rank_threshold =
        64.0 * std::numeric_limits<double>::epsilon() *
        maximum_pre_qr_norm;
    for (std::size_t index = 0U; index < offered; ++index) {
      if (!(pre_qr_norms[index] > rank_threshold)) continue;

      double image_norm = 0.0;
      if (retained != 0U) {
        Status orthogonalization_status{};
        for (std::size_t pass = 0U; pass < 2U; ++pass) {
          std::array<double, kLinearRecycleMaximumDirections> local_dots{};
          std::array<double, kLinearRecycleMaximumDirections> global_dots{};
          Status arithmetic{};
          for (std::size_t prior = 0U; prior < retained; ++prior) {
            const std::size_t prior_index = retained_indices[prior];
            double local = 0.0;
            arithmetic = merge_status(
                arithmetic,
                local_dot(as_const(images[prior_index]),
                          as_const(images[index]), local));
            local_dots[prior] = local;
          }
          const Status reduced = reductions.checked_sum(
              {local_dots.data(), retained}, {global_dots.data(), retained},
              merge_status(orthogonalization_status, arithmetic));
          if (!reduced) return finish_projection(reduced);
          orthogonalization_status = {};
          for (std::size_t prior = 0U; prior < retained; ++prior) {
            const std::size_t prior_index = retained_indices[prior];
            const double coefficient = global_dots[prior];
            if (!std::isfinite(coefficient)) {
              orthogonalization_status = merge_status(
                  orthogonalization_status,
                  {StatusCode::numerical_failure, kLinearSolveNonFinite});
              continue;
            }
            orthogonalization_status = merge_status(
                orthogonalization_status,
                add_scaled(images[index], -coefficient,
                           as_const(images[prior_index])));
            orthogonalization_status = merge_status(
                orthogonalization_status,
                add_scaled(directions[index], -coefficient,
                           as_const(directions[prior_index])));
          }
          orthogonalization_status = merge_status(
              orthogonalization_status,
              workspace.revise_vector(image_slots[index]));
          orthogonalization_status = merge_status(
              orthogonalization_status,
              workspace.revise_vector(
                  workspace.recycle_correction_logical_slot(index)));
        }
        const Status norm_status = global_norm(
            reductions, as_const(images[index]), orthogonalization_status,
            image_norm);
        if (!norm_status) return finish_projection(norm_status);
      } else {
        // No prior retained image means neither CGS pass has changed this
        // current A2*d image.  Reuse the bit-identical checked pre-QR norm
        // instead of repeating its max/sum packets.
        image_norm = pre_qr_norms[index];
      }
      if (!(image_norm > rank_threshold) || !std::isfinite(image_norm)) {
        continue;
      }
      const double inverse_norm = 1.0 / image_norm;
      Status normalized = scale_copy(images[index], as_const(images[index]),
                                     inverse_norm);
      normalized = merge_status(
          normalized,
          scale_copy(directions[index], as_const(directions[index]),
                     inverse_norm));
      normalized = merge_status(
          normalized, workspace.revise_vector(image_slots[index]));
      normalized = merge_status(
          normalized,
          workspace.revise_vector(
              workspace.recycle_correction_logical_slot(index)));
      normalized = reductions.consensus(normalized);
      if (!normalized) return finish_projection(normalized);
      retained_indices[retained] = index;
      ++retained;
    }

    if (retained == 0U) return finish_projection({});

    std::array<double, kLinearRecycleMaximumDirections> local_alpha{};
    std::array<double, kLinearRecycleMaximumDirections> global_alpha{};
    Status alpha_arithmetic{};
    for (std::size_t index = 0U; index < retained; ++index) {
      const std::size_t retained_index = retained_indices[index];
      double local = 0.0;
      alpha_arithmetic = merge_status(
          alpha_arithmetic,
          local_dot(as_const(images[retained_index]),
                    as_const(current_residual), local));
      local_alpha[index] = local;
    }
    Status alpha_status = reductions.checked_sum(
        {local_alpha.data(), retained}, {global_alpha.data(), retained},
        alpha_arithmetic);
    if (!alpha_status) return finish_projection(alpha_status);

    FieldView candidate = workspace.recycle_snapshot(shape);
    if (candidate.base == nullptr ||
        detail::field_views_overlap(as_const(candidate),
                                    as_const(current_residual))) {
      return finish_projection(
          {StatusCode::invalid_plan, kLinearSolveWorkspace});
    }
    Status candidate_status = copy_field(as_const(current_x), candidate);
    for (std::size_t index = 0U; index < retained; ++index) {
      const std::size_t retained_index = retained_indices[index];
      if (!std::isfinite(global_alpha[index])) {
        candidate_status = merge_status(
            candidate_status,
            {StatusCode::numerical_failure, kLinearSolveNonFinite});
      } else {
        candidate_status = merge_status(
            candidate_status,
            add_scaled(candidate, global_alpha[index],
                       as_const(directions[retained_index])));
      }
    }
    candidate_status = revise(workspace, workspace.recycle_snapshot_slot(),
                              candidate, candidate_status);
    candidate_status = reductions.consensus(candidate_status);
    if (!candidate_status) return finish_projection(candidate_status);

    FieldView candidate_operator = workspace.vector(ax_slot, shape);
    const OperatorApplyStatus operator_apply = apply_operator(
        linear_operator, candidate, candidate_operator, ax_slot, workspace,
        projection_operator_result);
    ++operator_applies;
    if (operator_apply.collective_failure) {
      result.lowest_failing_rank =
          projection_operator_result.lowest_failing_rank;
      return finish_projection(operator_apply.status);
    }
    Status residual_status = axpby(
        candidate_operator, 1.0, invocation.rhs, -1.0,
        as_const(candidate_operator));
    Status operator_status =
        merge_status(operator_apply.status, residual_status);
    operator_status = revise(workspace, ax_slot, candidate_operator,
                             operator_status);
    double candidate_norm = 0.0;
    operator_status = global_norm(reductions, as_const(candidate_operator),
                                  operator_status, candidate_norm);
    if (!operator_status) return finish_projection(operator_status);
    projected_residual = candidate_norm;
    if (std::isfinite(candidate_norm) && candidate_norm < baseline_residual) {
      Status admitted = copy_field(as_const(candidate), current_x);
      admitted = revise(workspace, x_slot, current_x, admitted);
      admitted = merge_status(
          admitted, copy_field(as_const(candidate_operator), current_residual));
      admitted = revise(workspace, v_begin, current_residual, admitted);
      admitted = reductions.consensus(admitted);
      if (!admitted) return finish_projection(admitted);
      result.final_true_residual = candidate_norm;
      result.recursive_residual = candidate_norm;
      projection_accepted = true;
    }
    return finish_projection({});
  };

  auto publish_recycle_cycle = [&](ConstFieldView solution) noexcept
      -> Status {
    if (!workspace.recycle_capture_active()) return {};
    const Status published =
        workspace.recycle_capture_cycle_publish(solution, reductions);
    if (published) {
      result.recycle_cycle_corrections =
          workspace.recycle_capture_cycle_corrections();
      result.recycle_capture_vector_passes =
          workspace.recycle_capture_vector_passes();
      result.recycle_capture_cycle_attempts =
          workspace.recycle_capture_cycle_attempts();
      result.recycle_capture_reduction_calls =
          workspace.recycle_capture_reduction_calls();
      result.recycle_capture_blocking_operations =
          workspace.recycle_capture_blocking_operations();
    }
    return published;
  };

  if (workspace.recycle_projection_pending()) {
    status = apply_recycle_projection(x, v0, result.final_true_residual);
    if (!status) {
      return finish_failure(
          result, status,
          status.detail == kLinearSolveNonFinite
              ? LinearTermination::non_finite
              : LinearTermination::operator_failure,
          resources, reductions, initial_calls);
    }
    if (result.recycle_projection_accepted &&
        result.final_true_residual <= tolerance) {
      bool audit_accepted = false;
      status = audit_convergence(invocation, as_const(x), as_const(v0),
                                 reductions, result, audit_accepted);
      if (!status) {
        return finish_failure(
            result, status, LinearTermination::convergence_audit_failure,
            resources, reductions, initial_calls);
      }
      if (audit_accepted) {
        return finish_success(result, as_const(x), invocation, resources,
                              reductions, initial_calls);
      }
      if (result.final_true_residual == 0.0) {
        return finish_failure(
            result,
            {StatusCode::rejected_step, kLinearSolveConvergenceAudit},
            LinearTermination::convergence_audit_failure, resources,
            reductions, initial_calls);
      }
      convergence_audit_deferred = true;
    }
  }

  const std::size_t rows = static_cast<std::size_t>(restart) + 1U;
  const std::size_t h_size = rows * static_cast<std::size_t>(restart);
  Span<double> h = workspace.scalars(0U, h_size);
  Span<double> cosine = workspace.scalars(h_size, restart);
  Span<double> sine = workspace.scalars(h_size + restart, restart);
  Span<double> g = workspace.scalars(h_size + 2U * restart, rows);
  Span<double> y = workspace.scalars(h_size + 2U * restart + rows,
                                     restart);
  Span<double> scratch = workspace.scalars(
      h_size + 3U * restart + rows, rows);
  if (h.data == nullptr || cosine.data == nullptr || sine.data == nullptr ||
      g.data == nullptr || y.data == nullptr || scratch.data == nullptr) {
    return finish_failure(
        result, {StatusCode::invalid_plan, kLinearSolveWorkspace},
        LinearTermination::invalid_plan, resources, reductions, initial_calls);
  }

  while (result.iterations < invocation.control.maximum_iterations) {
    if (workspace.recycle_capture_active()) {
      status = workspace.recycle_capture_cycle_start(
          as_const(x), reductions, initial_guess_selection_status);
      initial_guess_selection_status = {};
      if (!status) {
        return finish_failure(
            result, status,
            status.detail == kLinearSolveNonFinite
                ? LinearTermination::non_finite
                : LinearTermination::invalid_plan,
            resources, reductions, initial_calls);
      }
    }
    const double beta = result.final_true_residual;
    status = revise(workspace, v_begin, v0,
                    scale_copy(v0, as_const(v0), 1.0 / beta));
    status = reductions.consensus(status);
    if (!status) {
      return finish_failure(result, status, LinearTermination::non_finite,
                            resources, reductions, initial_calls);
    }
    for (std::size_t index = 0U; index < rows; ++index) {
      g.data[index] = 0.0;
    }
    g.data[0] = beta;
    bool restarted = false;
    // These fixed arrays are bound once per Arnoldi cycle and reused by the
    // two orthogonalization updates.  Only [0, dot_count) is read; leaving
    // the remaining capacity untouched avoids a full-capacity stack clear in
    // every column while retaining the validated uint8_t upper bound.
    std::array<ConstFieldView, kMaximumFgmresBasisUpdateCount> basis_views;
    std::array<double, kMaximumFgmresBasisUpdateCount> basis_scales;

    for (std::uint32_t column = 0U;
         column < restart &&
         result.iterations < invocation.control.maximum_iterations;
         ++column) {
      FieldView v = workspace.vector(
          static_cast<std::uint8_t>(v_begin + column), shape);
      const std::uint8_t z_slot =
          static_cast<std::uint8_t>(z_begin + column);
      FieldView z = workspace.vector(z_slot, shape);
      Status callback = apply_preconditioner(
          preconditioner, as_const(v), z, z_slot, result.iterations,
          workspace, result, prepared_batch ? &batch_ticket : nullptr);
      const Status preconditioner_consensus =
          preconditioner_status_scope ==
                  LinearPreconditionerStatusScope::collective
              ? callback
              : reductions.consensus(callback);
      if (!preconditioner_consensus) {
        return finish_failure(result, preconditioner_consensus,
                              LinearTermination::preconditioner_failure,
                              resources, reductions, initial_calls);
      }
      const OperatorApplyStatus operator_apply = apply_operator(
          linear_operator, z, ax, ax_slot, workspace, result);
      if (operator_apply.collective_failure) {
        return finish_failure(result, operator_apply.status,
                              LinearTermination::operator_failure, resources,
                              reductions, initial_calls);
      }
      callback = merge_status(
          operator_apply.status,
          revise(workspace, ax_slot, ax,
                 add_scaled(ax, -kSingleReductionFgmresShift,
                            as_const(v))));
      const std::size_t dot_count = static_cast<std::size_t>(column) + 1U;
      Status arithmetic{};
      for (std::uint32_t row = 0U; row <= column; ++row) {
        FieldView basis = workspace.vector(
            static_cast<std::uint8_t>(v_begin + row), shape);
        arithmetic = merge_status(
            arithmetic,
            local_dot(as_const(basis), as_const(ax), scratch.data[row]));
      }
      arithmetic = merge_status(
          arithmetic,
          local_dot(as_const(ax), as_const(ax), scratch.data[dot_count]));
      callback = merge_status(callback, arithmetic);
      double* const column_h =
          h.data + static_cast<std::size_t>(column) * rows;
      Status reduced = reductions.checked_sum(
          {scratch.data, dot_count + 1U},
          {column_h, dot_count + 1U}, callback);
      if (!reduced) {
        return finish_failure(
            result, reduced,
            reduced.detail == kLinearSolveNonFinite
                ? LinearTermination::non_finite
                : LinearTermination::operator_failure,
            resources, reductions, initial_calls);
      }
      const double bar_norm = column_h[dot_count];
      double t = 0.0;
      double next_norm = 0.0;
      for (std::uint32_t row = 0U; row <= column; ++row) {
        const FieldView basis = workspace.vector(
            static_cast<std::uint8_t>(v_begin + row), shape);
        basis_views[row] = as_const(basis);
        basis_scales[row] = -column_h[row];
      }
      status = single_reduction_fgmres_norm(
          column_h, dot_count, bar_norm, t, next_norm);
      const bool unsafe_recurrence =
          !status && status.detail == kLinearSolveBreakdown;
      if (!status && !unsafe_recurrence) {
        return finish_failure(result, status, LinearTermination::non_finite,
                              resources, reductions, initial_calls);
      }

      if (unsafe_recurrence && column != 0U) {
        ++result.iterations;
        const std::size_t valid_basis_count =
            static_cast<std::size_t>(column);
        for (std::size_t reverse = valid_basis_count; reverse > 0U;
             --reverse) {
          const std::size_t row = reverse - 1U;
          double value = g.data[row];
          for (std::size_t upper = row + 1U; upper < valid_basis_count;
               ++upper) {
            value -= h.data[upper * rows + row] * y.data[upper];
          }
          const double divisor = h.data[row * rows + row];
          if (divisor == 0.0 || !std::isfinite(divisor)) {
            return finish_failure(
                result, {StatusCode::numerical_failure, kLinearSolveBreakdown},
                LinearTermination::breakdown, resources, reductions,
                initial_calls);
          }
          y.data[row] = value / divisor;
          if (!std::isfinite(y.data[row])) {
            return finish_failure(
                result, {StatusCode::numerical_failure, kLinearSolveNonFinite},
                LinearTermination::non_finite, resources, reductions,
                initial_calls);
          }
        }
        status = {};
        for (std::size_t basis = 0U; basis < valid_basis_count; ++basis) {
          FieldView prior_z = workspace.vector(
              static_cast<std::uint8_t>(z_begin + basis), shape);
          status = merge_status(
              status, add_scaled(x, y.data[basis], as_const(prior_z)));
        }
        status = revise(workspace, x_slot, x, status);
        v0 = workspace.vector(v_begin, shape);
        ax = workspace.vector(ax_slot, shape);
        status = compute_true_residual(
            linear_operator, invocation.rhs, x, v0, v_begin, ax, ax_slot,
            workspace, reductions, result, result.final_true_residual);
        if (!status) {
          return finish_failure(result, status,
                                LinearTermination::operator_failure,
                                resources, reductions, initial_calls);
        }
        status = publish_recycle_cycle(as_const(x));
        if (!status) {
          return finish_failure(
              result, status,
              status.detail == kLinearSolveNonFinite
                  ? LinearTermination::non_finite
                  : LinearTermination::invalid_plan,
              resources, reductions, initial_calls);
        }
        result.recursive_residual = result.final_true_residual;
        if (result.final_true_residual <= tolerance) {
          bool audit_accepted = false;
          status = audit_convergence(invocation, as_const(x), as_const(v0),
                                     reductions, result, audit_accepted);
          if (!status) {
            return finish_failure(
                result, status,
                LinearTermination::convergence_audit_failure, resources,
                reductions, initial_calls);
          }
          if (audit_accepted) {
            return finish_success(result, as_const(x), invocation, resources,
                                  reductions, initial_calls);
          }
          if (result.final_true_residual == 0.0) {
            return finish_failure(
                result,
                {StatusCode::rejected_step, kLinearSolveConvergenceAudit},
                LinearTermination::convergence_audit_failure, resources,
                reductions, initial_calls);
          }
          convergence_audit_deferred = true;
        }
        if (result.iterations == invocation.control.maximum_iterations) {
          return finish_failure(
              result, {StatusCode::rejected_step, kLinearSolveBreakdown},
              LinearTermination::maximum_iterations, resources, reductions,
              initial_calls);
        }
        const double progress_guard =
            64.0 * std::numeric_limits<double>::epsilon() *
            std::max(beta, tolerance);
        if (!(result.final_true_residual + progress_guard < beta)) {
          return finish_failure(
              result, {StatusCode::numerical_failure, kLinearSolveBreakdown},
              LinearTermination::breakdown, resources, reductions,
              initial_calls);
        }
        if (result.norm_breakdown_restarts ==
            std::numeric_limits<std::uint64_t>::max()) {
          return finish_failure(
              result, {StatusCode::invalid_plan, kLinearSolveWorkspace},
              LinearTermination::invalid_plan, resources, reductions,
              initial_calls);
        }
        ++result.norm_breakdown_restarts;
        std::fill(h.data, h.data + h_size, 0.0);
        std::fill(cosine.data, cosine.data + restart, 0.0);
        std::fill(sine.data, sine.data + restart, 0.0);
        std::fill(g.data, g.data + rows, 0.0);
        std::fill(y.data, y.data + restart, 0.0);
        restarted = true;
        break;
      }

      bool explicitly_orthogonalized = false;
      if (unsafe_recurrence || next_norm == 0.0) {
        // The rearranged squared norm can become non-positive through
        // cancellation even when the Arnoldi residual is nonzero.  Form that
        // residual from the current operator vector and projections, apply a
        // second batched orthogonalization pass, then measure it stably
        // without repeating A/M application.
        Status recovery = add_scaled_basis(
            ax, basis_views.data(), basis_scales.data(), dot_count);
        recovery = revise(workspace, ax_slot, ax, recovery);
        Status correction_arithmetic{};
        for (std::uint32_t row = 0U; row <= column; ++row) {
          correction_arithmetic = merge_status(
              correction_arithmetic,
              local_dot(basis_views[row], as_const(ax), scratch.data[row]));
        }
        recovery = reductions.checked_sum(
            {scratch.data, dot_count}, {y.data, dot_count},
            merge_status(recovery, correction_arithmetic));
        if (!recovery) {
          return finish_failure(
              result, recovery,
              recovery.detail == kLinearSolveNonFinite
                  ? LinearTermination::non_finite
                  : LinearTermination::invalid_plan,
              resources, reductions, initial_calls);
        }
        Status correction_update{};
        for (std::uint32_t row = 0U; row <= column; ++row) {
          column_h[row] += y.data[row];
          if (!std::isfinite(column_h[row])) {
            correction_update = merge_status(
                correction_update,
                {StatusCode::numerical_failure, kLinearSolveNonFinite});
          }
          basis_scales[row] = -y.data[row];
        }
        correction_update = merge_status(
            correction_update,
            add_scaled_basis(ax, basis_views.data(), basis_scales.data(),
                             dot_count));
        recovery = revise(workspace, ax_slot, ax, correction_update);
        double explicit_norm = 0.0;
        recovery = global_norm(reductions, as_const(ax), recovery,
                               explicit_norm);
        if (!recovery) {
          return finish_failure(
              result, recovery,
              recovery.detail == kLinearSolveNonFinite
                  ? LinearTermination::non_finite
                  : LinearTermination::invalid_plan,
              resources, reductions, initial_calls);
        }
        next_norm = explicit_norm;
        explicitly_orthogonalized = true;
      }
      column_h[column] += kSingleReductionFgmresShift;
      if (!std::isfinite(column_h[column])) {
        return finish_failure(
            result, {StatusCode::numerical_failure, kLinearSolveNonFinite},
            LinearTermination::non_finite, resources, reductions,
            initial_calls);
      }
      column_h[column + 1U] = next_norm;

      for (std::uint32_t row = 0U; row < column; ++row) {
        double& upper =
            h.data[static_cast<std::size_t>(column) * rows + row];
        double& lower =
            h.data[static_cast<std::size_t>(column) * rows + row + 1U];
        const double rotated = cosine.data[row] * upper +
                               sine.data[row] * lower;
        lower = -sine.data[row] * upper + cosine.data[row] * lower;
        upper = rotated;
      }
      double& diagonal =
          h.data[static_cast<std::size_t>(column) * rows + column];
      double& subdiagonal =
          h.data[static_cast<std::size_t>(column) * rows + column + 1U];
      const double rotation_norm = std::hypot(diagonal, subdiagonal);
      if (!(rotation_norm > 0.0) || !std::isfinite(rotation_norm)) {
        return finish_failure(
            result, {StatusCode::numerical_failure, kLinearSolveBreakdown},
            LinearTermination::breakdown, resources, reductions,
            initial_calls);
      }
      cosine.data[column] = diagonal / rotation_norm;
      sine.data[column] = subdiagonal / rotation_norm;
      diagonal = rotation_norm;
      subdiagonal = 0.0;
      g.data[column + 1U] = -sine.data[column] * g.data[column];
      g.data[column] = cosine.data[column] * g.data[column];
      result.recursive_residual = std::abs(g.data[column + 1U]);
      ++result.iterations;

      const bool happy_breakdown = next_norm == 0.0;
      const bool cycle_end = happy_breakdown ||
                             column + 1U == restart ||
                             result.iterations ==
                                 invocation.control.maximum_iterations;
      const bool verify = happy_breakdown ||
                          (!convergence_audit_deferred &&
                           result.recursive_residual <= tolerance) ||
                          result.iterations %
                                  invocation.control.true_residual_interval ==
                              0U ||
                          cycle_end;
      if (!cycle_end) {
        Status basis_update_status{};
        if (!explicitly_orthogonalized) {
          basis_update_status = add_scaled_basis(
              ax, basis_views.data(), basis_scales.data(), dot_count);
          basis_update_status =
              revise(workspace, ax_slot, ax, basis_update_status);
        }
        FieldView next = workspace.vector(
            static_cast<std::uint8_t>(v_begin + column + 1U), shape);
        status = revise(workspace,
                        static_cast<std::uint8_t>(v_begin + column + 1U), next,
                        scale_copy(next, as_const(ax), 1.0 / next_norm));
        status = merge_status(status, basis_update_status);
        status = reductions.consensus(status);
        if (!status) {
          return finish_failure(result, status, LinearTermination::non_finite,
                                resources, reductions, initial_calls);
        }
      }
      if (!verify) {
        continue;
      }

      const std::size_t basis_count = static_cast<std::size_t>(column) + 1U;
      for (std::size_t reverse = basis_count; reverse > 0U; --reverse) {
        const std::size_t row = reverse - 1U;
        double value = g.data[row];
        for (std::size_t upper = row + 1U; upper < basis_count; ++upper) {
          value -= h.data[upper * rows + row] * y.data[upper];
        }
        const double divisor = h.data[row * rows + row];
        if (divisor == 0.0 || !std::isfinite(divisor)) {
          return finish_failure(
              result, {StatusCode::numerical_failure, kLinearSolveBreakdown},
              LinearTermination::breakdown, resources, reductions,
              initial_calls);
        }
        y.data[row] = value / divisor;
        if (!std::isfinite(y.data[row])) {
          return finish_failure(
              result, {StatusCode::numerical_failure, kLinearSolveNonFinite},
              LinearTermination::non_finite, resources, reductions,
              initial_calls);
        }
      }
      if (!cycle_end) {
        // A true-residual audit must not silently shorten the configured
        // Arnoldi cycle.  The final basis slot is still unused before the
        // cycle end, so form and audit the candidate there while preserving
        // x, V, and Z for the next column.
        const std::uint8_t trial_slot =
            static_cast<std::uint8_t>(v_begin + restart);
        FieldView trial = workspace.vector(trial_slot, shape);
        status = copy_field(as_const(x), trial);
        for (std::size_t basis = 0U; basis < basis_count; ++basis) {
          FieldView z = workspace.vector(
              static_cast<std::uint8_t>(z_begin + basis), shape);
          status = merge_status(
              status, add_scaled(trial, y.data[basis], as_const(z)));
        }
        status = revise(workspace, trial_slot, trial, status);
        ax = workspace.vector(ax_slot, shape);
        double audited_residual = 0.0;
        status = compute_true_residual_norm_in_place(
            linear_operator, invocation.rhs, trial, ax, ax_slot, workspace,
            reductions, result, audited_residual);
        if (!status) {
          return finish_failure(result, status,
                                LinearTermination::operator_failure,
                                resources, reductions, initial_calls);
        }
        result.final_true_residual = audited_residual;
        if (audited_residual <= tolerance) {
          bool audit_accepted = false;
          status = audit_convergence(invocation, as_const(trial),
                                     as_const(ax), reductions, result,
                                     audit_accepted);
          if (!status) {
            return finish_failure(
                result, status,
                LinearTermination::convergence_audit_failure, resources,
                reductions, initial_calls);
          }
          if (audit_accepted) {
            status = publish_recycle_cycle(as_const(trial));
            if (!status) {
              return finish_failure(
                  result, status,
                  status.detail == kLinearSolveNonFinite
                      ? LinearTermination::non_finite
                      : LinearTermination::invalid_plan,
                  resources, reductions, initial_calls);
            }
            return finish_success(result, as_const(trial), invocation,
                                  resources, reductions, initial_calls);
          }
          if (audited_residual == 0.0) {
            return finish_failure(
                result,
                {StatusCode::rejected_step, kLinearSolveConvergenceAudit},
                LinearTermination::convergence_audit_failure, resources,
                reductions, initial_calls);
          }
          convergence_audit_deferred = true;
        }
        continue;
      }

      status = {};
      for (std::size_t basis = 0U; basis < basis_count; ++basis) {
        FieldView z = workspace.vector(
            static_cast<std::uint8_t>(z_begin + basis), shape);
        status = merge_status(
            status, add_scaled(x, y.data[basis], as_const(z)));
      }
      status = revise(workspace, x_slot, x, status);
      v0 = workspace.vector(v_begin, shape);
      ax = workspace.vector(ax_slot, shape);
      status = compute_true_residual(
          linear_operator, invocation.rhs, x, v0, v_begin, ax, ax_slot,
          workspace, reductions, result, result.final_true_residual);
      if (!status) {
        return finish_failure(result, status,
                              LinearTermination::operator_failure, resources,
                              reductions, initial_calls);
      }
      status = publish_recycle_cycle(as_const(x));
      if (!status) {
        return finish_failure(
            result, status,
            status.detail == kLinearSolveNonFinite
                ? LinearTermination::non_finite
                : LinearTermination::invalid_plan,
            resources, reductions, initial_calls);
      }
      result.recursive_residual = result.final_true_residual;
      if (result.final_true_residual <= tolerance) {
        bool audit_accepted = false;
        status = audit_convergence(invocation, as_const(x), as_const(v0),
                                   reductions, result, audit_accepted);
        if (!status) {
          return finish_failure(
              result, status,
              LinearTermination::convergence_audit_failure, resources,
              reductions, initial_calls);
        }
        if (audit_accepted) {
          return finish_success(result, as_const(x), invocation, resources,
                                reductions, initial_calls);
        }
        if (result.final_true_residual == 0.0) {
          return finish_failure(
              result,
              {StatusCode::rejected_step, kLinearSolveConvergenceAudit},
              LinearTermination::convergence_audit_failure, resources,
              reductions, initial_calls);
        }
        convergence_audit_deferred = true;
      }
      if (result.iterations == invocation.control.maximum_iterations) {
        return finish_failure(
            result, {StatusCode::rejected_step, kLinearSolveBreakdown},
            LinearTermination::maximum_iterations, resources, reductions,
            initial_calls);
      }
      if (happy_breakdown) {
        const double progress_guard =
            64.0 * std::numeric_limits<double>::epsilon() *
            std::max(beta, tolerance);
        if (!(result.final_true_residual + progress_guard < beta)) {
          return finish_failure(
              result, {StatusCode::numerical_failure, kLinearSolveBreakdown},
              LinearTermination::breakdown, resources, reductions,
              initial_calls);
        }
        if (result.norm_breakdown_restarts ==
            std::numeric_limits<std::uint64_t>::max()) {
          return finish_failure(
              result, {StatusCode::invalid_plan, kLinearSolveWorkspace},
              LinearTermination::invalid_plan, resources, reductions,
              initial_calls);
        }
        ++result.norm_breakdown_restarts;
      }
      restarted = true;
      break;
    }
    if (!restarted) {
      break;
    }
  }
  return finish_failure(result,
                        {StatusCode::rejected_step, kLinearSolveBreakdown},
                        LinearTermination::maximum_iterations, resources,
                        reductions, initial_calls);
}

LinearSolveResult solve_bicgstab(const LinearOperator& linear_operator,
                                 LinearPreconditioner& preconditioner,
                                 const LinearSolveInvocation& invocation,
                                 SolverWorkspace& workspace,
                                 ReductionEngine& reductions,
                                 ResourceCounters* resources) noexcept {
  const std::uint64_t initial_calls = reductions.counters().calls;
  LinearPreconditionerStatusScope preconditioner_status_scope{
      LinearPreconditionerStatusScope::rank_local};
  LinearPreconditionerApplyLifecycle preconditioner_apply_lifecycle{
      LinearPreconditionerApplyLifecycle::per_call_checked};
  LinearPreconditionerBatchTicket batch_ticket{};
  const Status prepared = prepare_solve(
      LinearAlgorithm::bicgstab, linear_operator, preconditioner, invocation,
      workspace, reductions, &preconditioner_status_scope, &batch_ticket,
      &preconditioner_apply_lifecycle);
  if (!prepared) {
    return invalid_result(prepared, LinearTermination::invalid_plan,
                          reductions, initial_calls);
  }
  const bool prepared_batch =
      preconditioner_apply_lifecycle ==
      LinearPreconditionerApplyLifecycle::prepared_batch;
  LinearSolveResult result;
  const Int3 shape = invocation.rhs.interior;
  FieldView x = workspace.vector(0U, shape);
  FieldView r = workspace.vector(1U, shape);
  FieldView r_hat = workspace.vector(2U, shape);
  FieldView p = workspace.vector(3U, shape);
  FieldView s = workspace.vector(4U, shape);
  FieldView p_hat = workspace.vector(5U, shape);
  FieldView s_hat = workspace.vector(6U, shape);
  FieldView v = workspace.vector(7U, shape);
  FieldView t = workspace.vector(8U, shape);
  FieldView ax = workspace.vector(9U, shape);

  double rhs_norm = 0.0;
  Status status = global_norm(reductions, invocation.rhs, {}, rhs_norm);
  if (!status) {
    return finish_failure(result, status, LinearTermination::non_finite,
                          resources, reductions, initial_calls);
  }
  if (rhs_norm == 0.0) {
    fill_field(x, 0.0);
    status = revise(workspace, 0U, x);
    status = reductions.consensus(status);
    status = merge_status(
        status, collective_account(resources, 0U, reductions));
    if (!status) {
      return finish_failure(result, status, LinearTermination::invalid_plan,
                            nullptr, reductions, initial_calls);
    }
    if (invocation.convergence_audit != nullptr) {
      fill_field(r, 0.0);
      status = revise(workspace, 1U, r);
      status = reductions.consensus(status);
      bool audit_accepted = false;
      if (status) {
        status = audit_convergence(invocation, as_const(x), as_const(r),
                                   reductions, result, audit_accepted);
      }
      if (!status || !audit_accepted) {
        return finish_failure(result,
                              status ? Status{StatusCode::rejected_step,
                                              kLinearSolveConvergenceAudit}
                                     : status,
                              LinearTermination::convergence_audit_failure,
                              resources, reductions, initial_calls);
      }
    }
    status = copy_field(as_const(x), invocation.solution);
    if (!status) {
      return finish_failure(result, status, LinearTermination::non_finite,
                            nullptr, reductions, initial_calls);
    }
    result.status = {};
    result.termination = LinearTermination::zero_rhs;
    update_reduction_count(result, reductions, initial_calls);
    return result;
  }
  status = revise(workspace, 0U, x, copy_field(invocation.solution, x));
  status = reductions.consensus(status);
  if (!status) {
    return finish_failure(result, status, LinearTermination::invalid_plan,
                          resources, reductions, initial_calls);
  }
  status = compute_true_residual(linear_operator, invocation.rhs, x, r, 1U,
                                 ax, 9U, workspace, reductions, result,
                                 result.initial_true_residual);
  if (!status) {
    return finish_failure(result, status, LinearTermination::operator_failure,
                          resources, reductions, initial_calls);
  }
  result.final_true_residual = result.initial_true_residual;
  result.recursive_residual = result.initial_true_residual;
  const double tolerance = convergence_limit(invocation, rhs_norm);
  if (!std::isfinite(tolerance)) {
    return finish_failure(
        result, {StatusCode::numerical_failure, kLinearSolveNonFinite},
        LinearTermination::non_finite, resources, reductions, initial_calls);
  }
  if (result.initial_true_residual <= tolerance) {
    bool audit_accepted = false;
    status = audit_convergence(invocation, as_const(x), as_const(r), reductions,
                               result, audit_accepted);
    if (!status) {
      return finish_failure(result, status,
                            LinearTermination::convergence_audit_failure,
                            resources, reductions, initial_calls);
    }
    if (audit_accepted) {
      return finish_success(result, as_const(x), invocation, resources,
                            reductions, initial_calls);
    }
    if (result.initial_true_residual == 0.0) {
      return finish_failure(
          result, {StatusCode::rejected_step, kLinearSolveConvergenceAudit},
          LinearTermination::convergence_audit_failure, resources, reductions,
          initial_calls);
    }
  }
  status = revise(workspace, 2U, r_hat, copy_field(as_const(r), r_hat));
  status = reductions.consensus(status);
  if (!status) {
    return finish_failure(result, status, LinearTermination::non_finite,
                          resources, reductions, initial_calls);
  }

  double rho_previous = 1.0;
  double alpha = 1.0;
  double omega = 1.0;
  double cached_rho_local = 0.0;
  bool cached_rho_available = false;
  bool fresh = true;
  while (result.iterations < invocation.control.maximum_iterations) {
    double rho = 0.0;
    if (cached_rho_available) {
      cached_rho_available = false;
      double global = 0.0;
      status = reductions.checked_sum({&cached_rho_local, 1U}, {&global, 1U},
                                      {});
      if (status && !std::isfinite(global)) {
        status = {StatusCode::numerical_failure, kLinearSolveNonFinite};
      }
      if (status) {
        rho = global;
      }
    } else {
      status = global_dot(reductions, as_const(r_hat), as_const(r), {}, rho);
    }
    if (!status) {
      return finish_failure(result, status, LinearTermination::non_finite,
                            resources, reductions, initial_calls);
    }
    if (rho == 0.0 || !std::isfinite(rho)) {
      return finish_failure(
          result, {StatusCode::numerical_failure, kLinearSolveBreakdown},
          LinearTermination::breakdown, resources, reductions, initial_calls);
    }
    if (fresh) {
      status = revise(workspace, 3U, p, copy_field(as_const(r), p));
    } else {
      if (omega == 0.0 || !std::isfinite(omega)) {
        return finish_failure(
            result, {StatusCode::numerical_failure, kLinearSolveBreakdown},
            LinearTermination::breakdown, resources, reductions,
            initial_calls);
      }
      const double beta = (rho / rho_previous) * (alpha / omega);
      status = update_bicgstab_search_direction(
          p, as_const(r), as_const(v), beta, omega);
      status = revise(workspace, 3U, p, status);
    }
    status = reductions.consensus(status);
    if (!status) {
      return finish_failure(result, status, LinearTermination::non_finite,
                            resources, reductions, initial_calls);
    }

    const Status preconditioner_status = apply_preconditioner(
        preconditioner, as_const(p), p_hat, 5U, result.iterations * 2U,
        workspace, result, prepared_batch ? &batch_ticket : nullptr);
    const Status preconditioner_consensus =
        preconditioner_status_scope ==
                LinearPreconditionerStatusScope::collective
            ? preconditioner_status
            : reductions.consensus(preconditioner_status);
    if (!preconditioner_consensus) {
      return finish_failure(result, preconditioner_consensus,
                            LinearTermination::preconditioner_failure,
                            resources, reductions, initial_calls);
    }
    const OperatorApplyStatus operator_apply =
        apply_operator(linear_operator, p_hat, v, 7U, workspace, result);
    if (operator_apply.collective_failure) {
      return finish_failure(result, operator_apply.status,
                            LinearTermination::operator_failure, resources,
                            reductions, initial_calls);
    }
    double denominator = 0.0;
    status = global_dot(reductions, as_const(r_hat), as_const(v),
                        operator_apply.status, denominator);
    if (!status) {
      return finish_failure(
          result, status, LinearTermination::operator_failure, resources,
          reductions, initial_calls);
    }
    if (denominator == 0.0 || !std::isfinite(denominator)) {
      return finish_failure(
          result, {StatusCode::numerical_failure, kLinearSolveBreakdown},
          LinearTermination::breakdown, resources, reductions, initial_calls);
    }
    alpha = rho / denominator;
    double s_local_scale = 0.0;
    status = revise(
        workspace, 4U, s,
        form_bicgstab_intermediate(s, as_const(r), as_const(v), alpha,
                                   s_local_scale));
    double s_norm = 0.0;
    status = status
                 ? global_norm_from_local_scale(reductions, as_const(s),
                                                s_local_scale, {}, s_norm)
                 : global_norm(reductions, as_const(s), status, s_norm);
    if (!status) {
      return finish_failure(result, status, LinearTermination::non_finite,
                            resources, reductions, initial_calls);
    }
    if (s_norm <= tolerance) {
      status = revise(workspace, 0U, x,
                      add_scaled(x, alpha, as_const(p_hat)));
      ++result.iterations;
      status = compute_true_residual(
          linear_operator, invocation.rhs, x, r, 1U, ax, 9U, workspace,
          reductions, result, result.final_true_residual);
      if (!status) {
        return finish_failure(result, status,
                              LinearTermination::operator_failure, resources,
                              reductions, initial_calls);
      }
      result.recursive_residual = result.final_true_residual;
      if (result.final_true_residual <= tolerance) {
        bool audit_accepted = false;
        status = audit_convergence(invocation, as_const(x), as_const(r),
                                   reductions, result, audit_accepted);
        if (!status) {
          return finish_failure(result, status,
                                LinearTermination::convergence_audit_failure,
                                resources, reductions, initial_calls);
        }
        if (audit_accepted) {
          return finish_success(result, as_const(x), invocation, resources,
                                reductions, initial_calls);
        }
        if (result.final_true_residual == 0.0) {
          return finish_failure(
              result, {StatusCode::rejected_step, kLinearSolveConvergenceAudit},
              LinearTermination::convergence_audit_failure, resources,
              reductions, initial_calls);
        }
      }
      if (result.iterations == invocation.control.maximum_iterations) {
        return finish_failure(
            result, {StatusCode::rejected_step, kLinearSolveBreakdown},
            LinearTermination::maximum_iterations, resources, reductions,
            initial_calls);
      }
      status = revise(workspace, 2U, r_hat,
                      copy_field(as_const(r), r_hat));
      status = reductions.consensus(status);
      if (!status) {
        return finish_failure(result, status, LinearTermination::non_finite,
                              resources, reductions, initial_calls);
      }
      rho_previous = 1.0;
      alpha = 1.0;
      omega = 1.0;
      fresh = true;
      continue;
    }

    const Status s_preconditioner_status = apply_preconditioner(
        preconditioner, as_const(s), s_hat, 6U, result.iterations * 2U + 1U,
        workspace, result, prepared_batch ? &batch_ticket : nullptr);
    const Status s_preconditioner_consensus =
        preconditioner_status_scope ==
                LinearPreconditionerStatusScope::collective
            ? s_preconditioner_status
            : reductions.consensus(s_preconditioner_status);
    if (!s_preconditioner_consensus) {
      return finish_failure(result, s_preconditioner_consensus,
                            LinearTermination::preconditioner_failure,
                            resources, reductions, initial_calls);
    }
    const OperatorApplyStatus t_operator_apply = apply_operator(
        linear_operator, s_hat, t, 8U, workspace, result);
    if (t_operator_apply.collective_failure) {
      return finish_failure(result, t_operator_apply.status,
                            LinearTermination::operator_failure, resources,
                            reductions, initial_calls);
    }
    double local[2]{};
    double global[2]{};
    Status arithmetic =
        local_dot_pair(as_const(t), as_const(s), local[0], local[1]);
    arithmetic = merge_status(t_operator_apply.status, arithmetic);
    status = reductions.checked_sum({local, 2U}, {global, 2U}, arithmetic);
    if (!status) {
      return finish_failure(result, status, LinearTermination::operator_failure,
                            resources, reductions, initial_calls);
    }
    if (!(global[1] > 0.0) || !std::isfinite(global[0]) ||
        !std::isfinite(global[1])) {
      return finish_failure(
          result, {StatusCode::numerical_failure, kLinearSolveBreakdown},
          LinearTermination::breakdown, resources, reductions, initial_calls);
    }
    omega = global[0] / global[1];
    if (omega == 0.0 || !std::isfinite(omega)) {
      return finish_failure(
          result, {StatusCode::numerical_failure, kLinearSolveBreakdown},
          LinearTermination::breakdown, resources, reductions, initial_calls);
    }
    double r_local_scale = 0.0;
    double produced_rho_local = 0.0;
    bool produced_rho_available = false;
    const ConstFieldView rho_left = as_const(r_hat);
    status = update_bicgstab_solution_residual(
        x, as_const(p_hat), as_const(s_hat), r, as_const(s), as_const(t),
        alpha, omega, &r_local_scale, &rho_left, &produced_rho_local,
        &produced_rho_available);
    status = revise(workspace, 0U, x, status);
    status = revise(workspace, 1U, r, status);
    ++result.iterations;
    double recursive = 0.0;
    status = status
                 ? global_norm_from_local_scale(reductions, as_const(r),
                                                r_local_scale, {}, recursive)
                 : global_norm(reductions, as_const(r), status, recursive);
    if (!status) {
      return finish_failure(result, status, LinearTermination::non_finite,
                            resources, reductions, initial_calls);
    }
    result.recursive_residual = recursive;
    const bool verify = recursive <= tolerance ||
                        result.iterations %
                                invocation.control.true_residual_interval ==
                            0U ||
                        result.iterations ==
                            invocation.control.maximum_iterations;
    if (verify) {
      cached_rho_available = false;
      status = compute_true_residual(
          linear_operator, invocation.rhs, x, r, 1U, ax, 9U, workspace,
          reductions, result, result.final_true_residual);
      if (!status) {
        return finish_failure(result, status,
                              LinearTermination::operator_failure, resources,
                              reductions, initial_calls);
      }
      result.recursive_residual = result.final_true_residual;
      if (result.final_true_residual <= tolerance) {
        bool audit_accepted = false;
        status = audit_convergence(invocation, as_const(x), as_const(r),
                                   reductions, result, audit_accepted);
        if (!status) {
          return finish_failure(result, status,
                                LinearTermination::convergence_audit_failure,
                                resources, reductions, initial_calls);
        }
        if (audit_accepted) {
          return finish_success(result, as_const(x), invocation, resources,
                                reductions, initial_calls);
        }
        if (result.final_true_residual == 0.0) {
          return finish_failure(
              result, {StatusCode::rejected_step, kLinearSolveConvergenceAudit},
              LinearTermination::convergence_audit_failure, resources,
              reductions, initial_calls);
        }
      }
      if (result.iterations == invocation.control.maximum_iterations) {
        return finish_failure(
            result, {StatusCode::rejected_step, kLinearSolveBreakdown},
            LinearTermination::maximum_iterations, resources, reductions,
            initial_calls);
      }
      status = revise(workspace, 2U, r_hat,
                      copy_field(as_const(r), r_hat));
      status = reductions.consensus(status);
      if (!status) {
        return finish_failure(result, status, LinearTermination::non_finite,
                              resources, reductions, initial_calls);
      }
      rho_previous = 1.0;
      alpha = 1.0;
      omega = 1.0;
      fresh = true;
      continue;
    }
    cached_rho_local = produced_rho_local;
    cached_rho_available = produced_rho_available;
    rho_previous = rho;
    fresh = false;
  }
  return finish_failure(result,
                        {StatusCode::rejected_step, kLinearSolveBreakdown},
                        LinearTermination::maximum_iterations, resources,
                        reductions, initial_calls);
}

}  // namespace hundun::v04
