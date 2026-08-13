// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_linear.hpp"

#include "field_view_interval_detail.hpp"

#include <algorithm>
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
  hash = mix_contract(hash, invocation.expected_identity.numeric);
  hash = mix_contract(hash, invocation.expected_identity.hierarchy);
  hash = mix_contract(hash, invocation.expected_identity.workspace != 0U);
  hash = mix_contract(hash, invocation.expected_identity.fingerprint != 0U);
  hash = mix_contract(hash, op.collective_fingerprint);
  hash = mix_contract(hash, preconditioner.collective_fingerprint);
  hash = mix_contract(hash, static_cast<std::uint64_t>(op.operator_class));
  hash = mix_contract(
      hash, static_cast<std::uint64_t>(preconditioner.preconditioner_class));
  hash = mix_contract(hash, double_bits(invocation.control.absolute_tolerance));
  hash = mix_contract(hash, double_bits(invocation.control.relative_tolerance));
  hash = mix_contract(hash, invocation.control.maximum_iterations);
  hash = mix_contract(hash, invocation.control.true_residual_interval);
  hash = mix_contract(hash, invocation.control.restart);
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

Status apply_operator(const LinearOperator& linear_operator,
                      FieldView input, FieldView& output,
                      std::uint8_t output_slot, SolverWorkspace& workspace,
                      LinearSolveResult& result) noexcept {
  const Status applied = linear_operator.apply(input, output);
  ++result.operator_applies;
  return revise(workspace, output_slot, output, applied);
}

Status apply_preconditioner(LinearPreconditioner& preconditioner,
                            ConstFieldView input, FieldView& output,
                            std::uint8_t output_slot,
                            std::uint32_t iteration,
                            SolverWorkspace& workspace,
                            LinearSolveResult& result) noexcept {
  const Status applied = preconditioner.apply(input, output, iteration);
  ++result.preconditioner_applies;
  return revise(workspace, output_slot, output, applied);
}

void update_reduction_count(LinearSolveResult& result,
                            const ReductionEngine& reductions,
                            std::uint64_t initial_calls) noexcept {
  const std::uint64_t final_calls = reductions.counters().calls;
  result.reduction_calls =
      final_calls >= initial_calls ? final_calls - initial_calls : 0U;
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
      op.collective_fingerprint == 0U ||
      pc.collective_fingerprint == 0U ||
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
      control.true_residual_interval == 0U) {
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
             pc.preconditioner_class ==
                 LinearPreconditionerClass::flexible) {
    return {StatusCode::invalid_plan, kLinearSolvePlan};
  }
  return {};
}

Status prepare_solve(LinearAlgorithm algorithm,
                     const LinearOperator& linear_operator,
                     const LinearPreconditioner& preconditioner,
                     const LinearSolveInvocation& invocation,
                     const SolverWorkspace& workspace,
                     ReductionEngine& reductions) noexcept {
  const LinearOperatorCertificate op = linear_operator.certificate();
  const LinearPreconditionerCertificate pc = preconditioner.certificate();
  const Status valid = reductions.consensus(validate_common(
      algorithm, linear_operator, preconditioner, invocation, workspace,
      reductions));
  if (!valid) {
    return valid;
  }
  return reductions.consensus_contract(solve_contract_fingerprint(
      algorithm, op, pc, invocation, workspace.requirements()));
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
  Status pending = apply_operator(linear_operator, solution, operator_output,
                                  operator_slot, workspace, result);
  const Status formed = axpby(residual, 1.0, rhs, -1.0,
                              as_const(operator_output));
  pending = revise(workspace, residual_slot, residual,
                   merge_status(pending, formed));
  return global_norm(reductions, as_const(residual), pending, norm);
}

double convergence_limit(const LinearSolveInvocation& invocation,
                         double rhs_norm) noexcept {
  return std::max(invocation.control.absolute_tolerance,
                  invocation.control.relative_tolerance * rhs_norm);
}

LinearSolveResult finish_failure(LinearSolveResult result, Status status,
                                 LinearTermination termination,
                                 ResourceCounters* resources,
                                 ReductionEngine& reductions,
                                 std::uint64_t initial_calls) noexcept {
  result.status = status;
  result.termination = termination;
  const int original_failing_rank = reductions.lowest_failing_rank();
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
    status = apply_operator(linear_operator, p, ap, 4U, workspace, result);
    double denominator = 0.0;
    status = global_dot(reductions, as_const(p), as_const(ap), status,
                        denominator);
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
  const Status prepared = prepare_solve(LinearAlgorithm::fgmres,
                                        linear_operator, preconditioner,
                                        invocation, workspace, reductions);
  if (!prepared) {
    return invalid_result(prepared, LinearTermination::invalid_plan,
                          reductions, initial_calls);
  }
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
          workspace, result);
      const Status preconditioner_status = callback;
      const Status preconditioner_consensus =
          reductions.consensus(preconditioner_status);
      if (!preconditioner_consensus) {
        return finish_failure(result, preconditioner_consensus,
                              LinearTermination::preconditioner_failure,
                              resources, reductions, initial_calls);
      }
      const Status operator_status = apply_operator(
          linear_operator, z, ax, ax_slot, workspace, result);
      callback = operator_status;

      const std::size_t dot_count = static_cast<std::size_t>(column) + 1U;
      Status arithmetic{};
      for (std::uint32_t row = 0U; row <= column; ++row) {
        FieldView basis = workspace.vector(
            static_cast<std::uint8_t>(v_begin + row), shape);
        arithmetic = merge_status(
            arithmetic,
            local_dot(as_const(basis), as_const(ax), scratch.data[row]));
      }
      callback = merge_status(callback, arithmetic);
      Status reduced = reductions.checked_sum(
          {scratch.data, dot_count},
          {h.data + static_cast<std::size_t>(column) * rows, dot_count},
          callback);
      if (!reduced) {
        return finish_failure(
            result, reduced,
            !preconditioner_status
                ? LinearTermination::preconditioner_failure
                : (!operator_status ? LinearTermination::operator_failure
                                    : LinearTermination::non_finite),
            resources, reductions, initial_calls);
      }
      for (std::uint32_t row = 0U; row <= column; ++row) {
        FieldView basis = workspace.vector(
            static_cast<std::uint8_t>(v_begin + row), shape);
        status = add_scaled(
            ax,
            -h.data[static_cast<std::size_t>(column) * rows + row],
            as_const(basis));
        if (!status) {
          break;
        }
      }
      status = revise(workspace, ax_slot, ax, status);

      for (std::uint32_t row = 0U; row <= column; ++row) {
        FieldView basis = workspace.vector(
            static_cast<std::uint8_t>(v_begin + row), shape);
        arithmetic = merge_status(
            status,
            local_dot(as_const(basis), as_const(ax), scratch.data[row]));
        status = arithmetic;
      }
      reduced = reductions.checked_sum(
          {scratch.data, dot_count}, {y.data, dot_count}, status);
      if (!reduced) {
        return finish_failure(result, reduced, LinearTermination::non_finite,
                              resources, reductions, initial_calls);
      }
      for (std::uint32_t row = 0U; row <= column; ++row) {
        h.data[static_cast<std::size_t>(column) * rows + row] += y.data[row];
        FieldView basis = workspace.vector(
            static_cast<std::uint8_t>(v_begin + row), shape);
        status = add_scaled(ax, -y.data[row], as_const(basis));
        if (!status) {
          break;
        }
      }
      status = revise(workspace, ax_slot, ax, status);
      double next_norm = 0.0;
      status = global_norm(reductions, as_const(ax), status, next_norm);
      if (!status) {
        return finish_failure(result, status, LinearTermination::non_finite,
                              resources, reductions, initial_calls);
      }
      h.data[static_cast<std::size_t>(column) * rows + column + 1U] =
          next_norm;

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

      const bool happy_breakdown =
          next_norm <= std::numeric_limits<double>::epsilon() *
                           std::max(1.0, rotation_norm);
      const bool verify = happy_breakdown ||
                          result.recursive_residual <= tolerance ||
                          result.iterations %
                                  invocation.control.true_residual_interval ==
                              0U ||
                          column + 1U == restart ||
                          result.iterations ==
                              invocation.control.maximum_iterations;
      if (!verify) {
        FieldView next = workspace.vector(
            static_cast<std::uint8_t>(v_begin + column + 1U), shape);
        status = revise(workspace,
                        static_cast<std::uint8_t>(v_begin + column + 1U), next,
                        scale_copy(next, as_const(ax), 1.0 / next_norm));
        status = reductions.consensus(status);
        if (!status) {
          return finish_failure(result, status, LinearTermination::non_finite,
                                resources, reductions, initial_calls);
        }
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
      if (happy_breakdown) {
        return finish_failure(
            result, {StatusCode::numerical_failure, kLinearSolveBreakdown},
            LinearTermination::breakdown, resources, reductions,
            initial_calls);
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
  const Status prepared = prepare_solve(LinearAlgorithm::bicgstab,
                                        linear_operator, preconditioner,
                                        invocation, workspace, reductions);
  if (!prepared) {
    return invalid_result(prepared, LinearTermination::invalid_plan,
                          reductions, initial_calls);
  }
  LinearSolveResult result;
  const Int3 shape = invocation.rhs.interior;
  FieldView x = workspace.vector(0U, shape);
  FieldView r = workspace.vector(1U, shape);
  FieldView r_hat = workspace.vector(2U, shape);
  FieldView p = workspace.vector(3U, shape);
  FieldView v = workspace.vector(4U, shape);
  FieldView s = workspace.vector(5U, shape);
  FieldView p_hat = workspace.vector(6U, shape);
  FieldView s_hat = workspace.vector(7U, shape);
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
    return finish_success(result, as_const(x), invocation, resources,
                          reductions, initial_calls);
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
  bool fresh = true;
  while (result.iterations < invocation.control.maximum_iterations) {
    double rho = 0.0;
    status = global_dot(reductions, as_const(r_hat), as_const(r), {}, rho);
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
      status = axpby(p, 1.0, as_const(p), -omega, as_const(v));
      status = axpby(p, 1.0, as_const(r), beta, as_const(p));
      status = revise(workspace, 3U, p, status);
    }
    status = reductions.consensus(status);
    if (!status) {
      return finish_failure(result, status, LinearTermination::non_finite,
                            resources, reductions, initial_calls);
    }

    const Status preconditioner_status = apply_preconditioner(
        preconditioner, as_const(p), p_hat, 6U, result.iterations, workspace,
        result);
    const Status preconditioner_consensus =
        reductions.consensus(preconditioner_status);
    if (!preconditioner_consensus) {
      return finish_failure(result, preconditioner_consensus,
                            LinearTermination::preconditioner_failure,
                            resources, reductions, initial_calls);
    }
    const Status operator_status = apply_operator(
        linear_operator, p_hat, v, 4U, workspace, result);
    double denominator = 0.0;
    status = global_dot(reductions, as_const(r_hat), as_const(v),
                        operator_status, denominator);
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
    status = revise(workspace, 5U, s,
                    axpby(s, 1.0, as_const(r), -alpha, as_const(v)));
    double s_norm = 0.0;
    status = global_norm(reductions, as_const(s), status, s_norm);
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
        return finish_success(result, as_const(x), invocation, resources,
                              reductions, initial_calls);
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
        preconditioner, as_const(s), s_hat, 7U, result.iterations, workspace,
        result);
    const Status s_preconditioner_consensus =
        reductions.consensus(s_preconditioner_status);
    if (!s_preconditioner_consensus) {
      return finish_failure(result, s_preconditioner_consensus,
                            LinearTermination::preconditioner_failure,
                            resources, reductions, initial_calls);
    }
    const Status t_operator_status = apply_operator(
        linear_operator, s_hat, t, 8U, workspace, result);
    double local[2]{};
    double global[2]{};
    Status arithmetic = local_dot(as_const(t), as_const(s), local[0]);
    arithmetic = merge_status(
        arithmetic, local_dot(as_const(t), as_const(t), local[1]));
    arithmetic = merge_status(t_operator_status, arithmetic);
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
    status = add_scaled(x, alpha, as_const(p_hat));
    status = merge_status(status, add_scaled(x, omega, as_const(s_hat)));
    status = revise(workspace, 0U, x, status);
    status = revise(workspace, 1U, r,
                    merge_status(status,
                                 axpby(r, 1.0, as_const(s), -omega,
                                       as_const(t))));
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
          linear_operator, invocation.rhs, x, r, 1U, ax, 9U, workspace,
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
    rho_previous = rho;
    fresh = false;
  }
  return finish_failure(result,
                        {StatusCode::rejected_step, kLinearSolveBreakdown},
                        LinearTermination::maximum_iterations, resources,
                        reductions, initial_calls);
}

}  // namespace hundun::v04
