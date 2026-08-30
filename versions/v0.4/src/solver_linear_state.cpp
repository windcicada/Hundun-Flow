// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_linear.hpp"

#include "field_view_interval_detail.hpp"

#include "solver_linear_state_detail.hpp"

#include <cstddef>
#include <cstdint>
#include <cmath>
#include <limits>
#include <string>

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kLinearSymbolic = 501U;
constexpr std::uint32_t kLinearNumeric = 502U;
constexpr std::uint32_t kLinearHierarchy = 503U;
constexpr std::uint32_t kLinearWorkspace = 504U;
constexpr std::uint32_t kLinearCounter = 505U;
constexpr std::uint32_t kLinearRegistry = 506U;

bool valid_location(LinearLocation location) noexcept {
  return location == LinearLocation::cell ||
         location == LinearLocation::face_x ||
         location == LinearLocation::face_y ||
         location == LinearLocation::face_z;
}

bool valid_backend(LinearBackend backend) noexcept {
  return backend == LinearBackend::native_cartesian ||
         backend == LinearBackend::hypre_struct;
}

bool valid_algorithm(LinearAlgorithm algorithm) noexcept {
  return algorithm == LinearAlgorithm::pcg ||
         algorithm == LinearAlgorithm::fgmres ||
         algorithm == LinearAlgorithm::bicgstab;
}

bool valid_reduction(ReductionMode mode) noexcept {
  return mode == ReductionMode::reproducible_tree ||
         mode == ReductionMode::mpi_allreduce;
}

bool valid_shape(Int3 shape) noexcept {
  return shape.x > 0 && shape.y > 0 && shape.z > 0;
}

bool valid_symbolic_spec(SymbolicSpec spec) noexcept {
  return spec.operator_kind != 0U && valid_location(spec.location) &&
         spec.topology != 0U && spec.boundary_layout != 0U &&
         spec.partition != 0U && spec.eb_interface != 0U &&
         spec.stencil != 0U && valid_backend(spec.backend);
}

bool valid_coefficients(CoefficientRevisions revisions) noexcept {
  return revisions.diagonal != 0U && revisions.off_diagonal != 0U &&
         revisions.time != 0U && revisions.material_transport != 0U &&
         revisions.numeric_boundary != 0U && revisions.constraint != 0U;
}

bool valid_policy(HierarchyPolicyIdentity policy) noexcept {
  return policy.coarsening != 0U && policy.transfer_smoother != 0U &&
         policy.policy_epoch != 0U && valid_backend(policy.backend);
}

RevisionToken next_generation(RevisionToken current) noexcept {
  return current == 0U ? RevisionToken{1U} : current + 1U;
}

bool valid_view_shape(FieldView view) noexcept {
  if (view.base == nullptr || !valid_shape(view.interior) ||
      view.ghosts.x < 0 || view.ghosts.y < 0 || view.ghosts.z < 0 ||
      view.components == 0U || view.revision == 0U ||
      view.storage_identity == 0U || view.revision_domain == 0U ||
      view.stride_y == 0U || view.stride_z == 0U ||
      view.component_stride == 0U) {
    return false;
  }
  std::size_t padded_x = 0U;
  std::size_t padded_y = 0U;
  std::size_t padded_z = 0U;
  std::size_t two_ghosts = 0U;
  std::size_t minimum_stride_z = 0U;
  std::size_t minimum_component_stride = 0U;
  return detail::checked_linear_multiply(
             static_cast<std::size_t>(view.ghosts.x), 2U, two_ghosts) &&
         detail::checked_linear_add(static_cast<std::size_t>(view.interior.x),
                                    two_ghosts, padded_x) &&
         detail::checked_linear_multiply(
             static_cast<std::size_t>(view.ghosts.y), 2U, two_ghosts) &&
         detail::checked_linear_add(static_cast<std::size_t>(view.interior.y),
                                    two_ghosts, padded_y) &&
         detail::checked_linear_multiply(
             static_cast<std::size_t>(view.ghosts.z), 2U, two_ghosts) &&
         detail::checked_linear_add(static_cast<std::size_t>(view.interior.z),
                                    two_ghosts, padded_z) &&
         view.stride_y >= padded_x &&
         detail::checked_linear_multiply(view.stride_y, padded_y,
                                         minimum_stride_z) &&
         view.stride_z >= minimum_stride_z &&
         detail::checked_linear_multiply(view.stride_z, padded_z,
                                         minimum_component_stride) &&
         view.component_stride >= minimum_component_stride;
}

bool sufficient_workspace_views(
    const LinearWorkspaceRequirements& requirements, FieldView vector_bundle,
    FieldView scalar_buffer) noexcept {
  if (!valid_view_shape(vector_bundle) || !valid_view_shape(scalar_buffer) ||
      detail::field_views_overlap(as_const(vector_bundle),
                                  as_const(scalar_buffer)) ||
      vector_bundle.storage_identity != scalar_buffer.storage_identity ||
      vector_bundle.revision_domain != scalar_buffer.revision_domain ||
      vector_bundle.field == scalar_buffer.field ||
      vector_bundle.replica != 0U || scalar_buffer.replica != 0U ||
      vector_bundle.components < requirements.vector_slots ||
      vector_bundle.ghosts.x < requirements.ghost_width ||
      vector_bundle.ghosts.y < requirements.ghost_width ||
      vector_bundle.ghosts.z < requirements.ghost_width ||
      vector_bundle.interior.x < requirements.maximum_shape.x ||
      vector_bundle.interior.y < requirements.maximum_shape.y ||
      vector_bundle.interior.z < requirements.maximum_shape.z ||
      scalar_buffer.components != 1U || scalar_buffer.ghosts.x != 0 ||
      scalar_buffer.ghosts.y != 0 || scalar_buffer.ghosts.z != 0 ||
      scalar_buffer.interior.y != 1 || scalar_buffer.interior.z != 1 ||
      static_cast<std::size_t>(scalar_buffer.interior.x) <
          requirements.scalar_doubles) {
    return false;
  }
  return true;
}

bool same_shape(Int3 left, Int3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

bool shape_contains(Int3 outer, Int3 inner) noexcept {
  return valid_shape(inner) && inner.x <= outer.x && inner.y <= outer.y &&
         inner.z <= outer.z;
}

Status copy_recycle_field(ConstFieldView source, FieldView destination) noexcept {
  for (std::int32_t z = 0; z < destination.interior.z; ++z) {
    for (std::int32_t y = 0; y < destination.interior.y; ++y) {
      for (std::int32_t x = 0; x < destination.interior.x; ++x) {
        const double value = source.unchecked({x, y, z}, 0U);
        if (!std::isfinite(value)) {
          return {StatusCode::numerical_failure, kLinearWorkspace};
        }
        destination.unchecked({x, y, z}, 0U) = value;
      }
    }
  }
  return {};
}

Status difference_recycle_field(ConstFieldView solution,
                                ConstFieldView snapshot,
                                FieldView destination,
                                bool& nonzero) noexcept {
  nonzero = false;
  for (std::int32_t z = 0; z < destination.interior.z; ++z) {
    for (std::int32_t y = 0; y < destination.interior.y; ++y) {
      for (std::int32_t x = 0; x < destination.interior.x; ++x) {
        const double value = solution.unchecked({x, y, z}, 0U) -
                             snapshot.unchecked({x, y, z}, 0U);
        if (!std::isfinite(value)) {
          return {StatusCode::numerical_failure, kLinearWorkspace};
        }
        nonzero = nonzero || value != 0.0;
        destination.unchecked({x, y, z}, 0U) = value;
      }
    }
  }
  return {};
}

}  // namespace

namespace detail {

bool same_symbolic_spec(SymbolicSpec left, SymbolicSpec right) noexcept {
  return left.operator_kind == right.operator_kind &&
         left.location == right.location && left.topology == right.topology &&
         left.boundary_layout == right.boundary_layout &&
         left.partition == right.partition &&
         left.eb_interface == right.eb_interface &&
         left.stencil == right.stencil && left.backend == right.backend;
}

bool same_coefficients(CoefficientRevisions left,
                       CoefficientRevisions right) noexcept {
  return left.diagonal == right.diagonal &&
         left.off_diagonal == right.off_diagonal && left.time == right.time &&
         left.material_transport == right.material_transport &&
         left.numeric_boundary == right.numeric_boundary &&
         left.constraint == right.constraint;
}

bool same_hierarchy_policy(HierarchyPolicyIdentity left,
                           HierarchyPolicyIdentity right) noexcept {
  return left.coarsening == right.coarsening &&
         left.transfer_smoother == right.transfer_smoother &&
         left.policy_epoch == right.policy_epoch && left.backend == right.backend;
}

bool same_workspace_requirements(const LinearWorkspaceRequirements& left,
                                 const LinearWorkspaceRequirements& right)
    noexcept {
  return left.algorithm == right.algorithm &&
         left.maximum_shape.x == right.maximum_shape.x &&
         left.maximum_shape.y == right.maximum_shape.y &&
         left.maximum_shape.z == right.maximum_shape.z &&
         left.ghost_width == right.ghost_width &&
         left.maximum_restart == right.maximum_restart &&
         left.reduction_mode == right.reduction_mode &&
         left.execution_revision == right.execution_revision &&
         left.vector_slots == right.vector_slots &&
         left.scalar_doubles == right.scalar_doubles &&
         left.reduction_capacity == right.reduction_capacity &&
         left.fingerprint == right.fingerprint;
}

bool valid_workspace_requirements(
    const LinearWorkspaceRequirements& requirements) noexcept {
  LinearWorkspaceRequirements expected;
  return make_linear_workspace_requirements(
             requirements.algorithm, requirements.maximum_shape,
             requirements.ghost_width, requirements.maximum_restart,
             requirements.reduction_mode, requirements.execution_revision,
             expected) &&
         same_workspace_requirements(requirements, expected);
}

bool same_field_binding(FieldView left, FieldView right) noexcept {
  return left.base == right.base && left.interior.x == right.interior.x &&
         left.interior.y == right.interior.y &&
         left.interior.z == right.interior.z &&
         left.ghosts.x == right.ghosts.x &&
         left.ghosts.y == right.ghosts.y &&
         left.ghosts.z == right.ghosts.z &&
         left.components == right.components &&
         left.stride_y == right.stride_y && left.stride_z == right.stride_z &&
         left.component_stride == right.component_stride &&
         left.replica == right.replica && left.field == right.field &&
         left.revision == right.revision &&
         left.storage_identity == right.storage_identity &&
         left.revision_domain == right.revision_domain;
}

PlanFingerprint symbolic_fingerprint(SymbolicSpec spec) noexcept {
  std::uint64_t hash = kLinearFnvOffset;
  hash = linear_hash_mix(hash, spec.operator_kind);
  hash = linear_hash_mix(hash, static_cast<std::uint64_t>(spec.location));
  hash = linear_hash_mix(hash, spec.topology);
  hash = linear_hash_mix(hash, spec.boundary_layout);
  hash = linear_hash_mix(hash, spec.partition);
  hash = linear_hash_mix(hash, spec.eb_interface);
  hash = linear_hash_mix(hash, spec.stencil);
  hash = linear_hash_mix(hash, static_cast<std::uint64_t>(spec.backend));
  return finish_linear_hash(hash);
}

PlanFingerprint numeric_fingerprint(PlanFingerprint symbolic,
                                    CoefficientRevisions revisions,
                                    PlanFingerprint content) noexcept {
  std::uint64_t hash = linear_hash_mix(kLinearFnvOffset, symbolic);
  hash = linear_hash_mix(hash, revisions.diagonal);
  hash = linear_hash_mix(hash, revisions.off_diagonal);
  hash = linear_hash_mix(hash, revisions.time);
  hash = linear_hash_mix(hash, revisions.material_transport);
  hash = linear_hash_mix(hash, revisions.numeric_boundary);
  hash = linear_hash_mix(hash, revisions.constraint);
  hash = linear_hash_mix(hash, content);
  return finish_linear_hash(hash);
}

PlanFingerprint hierarchy_fingerprint(
    PlanFingerprint symbolic, PlanFingerprint numeric,
    HierarchyPolicyIdentity policy, PlanFingerprint content) noexcept {
  std::uint64_t hash = linear_hash_mix(kLinearFnvOffset, symbolic);
  hash = linear_hash_mix(hash, numeric);
  hash = linear_hash_mix(hash, policy.coarsening);
  hash = linear_hash_mix(hash, policy.transfer_smoother);
  hash = linear_hash_mix(hash, policy.policy_epoch);
  hash = linear_hash_mix(hash, static_cast<std::uint64_t>(policy.backend));
  hash = linear_hash_mix(hash, content);
  return finish_linear_hash(hash);
}

PlanFingerprint workspace_fingerprint(
    const LinearWorkspaceRequirements& requirements, FieldView vector_bundle,
    FieldView scalar_buffer) noexcept {
  std::uint64_t hash = linear_hash_mix(kLinearFnvOffset,
                                       requirements.fingerprint);
  hash = linear_hash_mix(hash, vector_bundle.storage_identity);
  hash = linear_hash_mix(hash, vector_bundle.revision_domain);
  hash = linear_hash_mix(hash, vector_bundle.field);
  hash = linear_hash_mix(hash, scalar_buffer.field);
  hash = linear_hash_mix(hash, vector_bundle.revision);
  hash = linear_hash_mix(hash, scalar_buffer.revision);
  hash = linear_hash_mix(
      hash, static_cast<std::uint64_t>(
                reinterpret_cast<std::uintptr_t>(vector_bundle.base)));
  hash = linear_hash_mix(
      hash, static_cast<std::uint64_t>(
                reinterpret_cast<std::uintptr_t>(scalar_buffer.base)));
  return finish_linear_hash(hash);
}

bool checked_linear_add(std::size_t left, std::size_t right,
                        std::size_t& out) noexcept {
  if (right > std::numeric_limits<std::size_t>::max() - left) {
    return false;
  }
  out = left + right;
  return true;
}

bool checked_linear_multiply(std::size_t left, std::size_t right,
                             std::size_t& out) noexcept {
  if (left != 0U && right > std::numeric_limits<std::size_t>::max() / left) {
    return false;
  }
  out = left * right;
  return true;
}

bool checked_counter_increment(std::uint64_t current,
                               std::uint64_t& out) noexcept {
  if (current == std::numeric_limits<std::uint64_t>::max()) {
    return false;
  }
  out = current + 1U;
  return true;
}

}  // namespace detail

Status SymbolicPlan::compile(SymbolicSpec spec, SymbolicPlan& out,
                             LinearLifecycleCounters* counters) noexcept {
  if (!valid_symbolic_spec(spec)) {
    return {StatusCode::invalid_plan, kLinearSymbolic};
  }
  if (out.valid() && detail::same_symbolic_spec(out.spec_, spec)) {
    return {};
  }
  const RevisionToken generation = next_generation(out.generation_);
  std::uint64_t counter = 0U;
  if (generation == 0U ||
      (counters != nullptr &&
       !detail::checked_counter_increment(counters->symbolic_builds,
                                          counter))) {
    return {StatusCode::invalid_plan, kLinearCounter};
  }
  SymbolicPlan candidate;
  candidate.spec_ = spec;
  candidate.fingerprint_ = detail::symbolic_fingerprint(spec);
  candidate.generation_ = generation;
  out = candidate;
  if (counters != nullptr) {
    counters->symbolic_builds = counter;
  }
  return {};
}

Status NumericState::refill(const SymbolicPlan& symbolic,
                            CoefficientRevisions revisions,
                            PlanFingerprint content_fingerprint,
                            LinearLifecycleCounters* counters) noexcept {
  if (!symbolic.valid() || !valid_coefficients(revisions) ||
      content_fingerprint == 0U) {
    return {StatusCode::invalid_plan, kLinearNumeric};
  }
  const bool same_authority =
      symbolic_fingerprint_ == symbolic.fingerprint_ &&
      symbolic_generation_ == symbolic.generation_;
  if (fingerprint_ != 0U && same_authority &&
      detail::same_coefficients(revisions_, revisions)) {
    return content_fingerprint_ == content_fingerprint
               ? Status{}
               : Status{StatusCode::invalid_plan, kLinearNumeric};
  }
  const RevisionToken generation = next_generation(generation_);
  std::uint64_t counter = 0U;
  if (generation == 0U ||
      (counters != nullptr &&
       !detail::checked_counter_increment(counters->numeric_refills,
                                          counter))) {
    return {StatusCode::invalid_plan, kLinearCounter};
  }
  NumericState candidate;
  candidate.revisions_ = revisions;
  candidate.symbolic_fingerprint_ = symbolic.fingerprint_;
  candidate.symbolic_generation_ = symbolic.generation_;
  candidate.content_fingerprint_ = content_fingerprint;
  candidate.fingerprint_ = detail::numeric_fingerprint(
      symbolic.fingerprint_, revisions, content_fingerprint);
  candidate.generation_ = generation;
  *this = candidate;
  if (counters != nullptr) {
    counters->numeric_refills = counter;
  }
  return {};
}

bool NumericState::valid_for(const SymbolicPlan& symbolic) const noexcept {
  return fingerprint_ != 0U && symbolic.valid() &&
         symbolic_fingerprint_ == symbolic.fingerprint_ &&
         symbolic_generation_ == symbolic.generation_;
}

Status HierarchyState::update(const SymbolicPlan& symbolic,
                              const NumericState& numeric,
                              HierarchyPolicyIdentity policy,
                              HierarchyUpdate decision,
                              PlanFingerprint content_fingerprint,
                              LinearLifecycleCounters* counters) noexcept {
  if (!numeric.valid_for(symbolic) || !valid_policy(policy) ||
      content_fingerprint == 0U) {
    return {StatusCode::invalid_plan, kLinearHierarchy};
  }
  const bool exists = fingerprint_ != 0U;
  const bool same_symbolic = symbolic_fingerprint_ == symbolic.fingerprint_ &&
                             symbolic_generation_ == symbolic.generation_;
  const bool same_numeric = numeric_fingerprint_ == numeric.fingerprint_ &&
                            numeric_generation_ == numeric.generation_;
  const bool same_policy = detail::same_hierarchy_policy(policy_, policy);
  const bool same_content = content_fingerprint_ == content_fingerprint;
  if (decision == HierarchyUpdate::retain) {
    return exists && same_symbolic && same_numeric && same_policy &&
                   same_content
               ? Status{}
               : Status{StatusCode::invalid_plan, kLinearHierarchy};
  }
  if (decision == HierarchyUpdate::refresh_numeric &&
      (!exists || !same_symbolic || !same_policy)) {
    return {StatusCode::invalid_plan, kLinearHierarchy};
  }
  if (decision != HierarchyUpdate::refresh_numeric &&
      decision != HierarchyUpdate::rebuild) {
    return {StatusCode::invalid_plan, kLinearHierarchy};
  }
  if (decision == HierarchyUpdate::refresh_numeric && same_numeric) {
    return same_content
               ? Status{}
               : Status{StatusCode::invalid_plan, kLinearHierarchy};
  }
  if (decision == HierarchyUpdate::rebuild && exists && same_symbolic &&
      same_numeric && same_policy) {
    return same_content
               ? Status{}
               : Status{StatusCode::invalid_plan, kLinearHierarchy};
  }
  const RevisionToken generation = next_generation(generation_);
  std::uint64_t counter = 0U;
  const std::uint64_t current =
      decision == HierarchyUpdate::rebuild
          ? (counters == nullptr ? 0U : counters->hierarchy_rebuilds)
          : (counters == nullptr ? 0U : counters->hierarchy_refreshes);
  if (generation == 0U ||
      (counters != nullptr &&
       !detail::checked_counter_increment(current, counter))) {
    return {StatusCode::invalid_plan, kLinearCounter};
  }
  HierarchyState candidate;
  candidate.policy_ = policy;
  candidate.symbolic_fingerprint_ = symbolic.fingerprint_;
  candidate.symbolic_generation_ = symbolic.generation_;
  candidate.numeric_fingerprint_ = numeric.fingerprint_;
  candidate.numeric_generation_ = numeric.generation_;
  candidate.content_fingerprint_ = content_fingerprint;
  candidate.fingerprint_ = detail::hierarchy_fingerprint(
      symbolic.fingerprint_, numeric.fingerprint_, policy,
      content_fingerprint);
  candidate.generation_ = generation;
  *this = candidate;
  if (counters != nullptr) {
    if (decision == HierarchyUpdate::rebuild) {
      counters->hierarchy_rebuilds = counter;
    } else {
      counters->hierarchy_refreshes = counter;
    }
  }
  return {};
}

bool HierarchyState::valid_for(const SymbolicPlan& symbolic,
                               const NumericState& numeric) const noexcept {
  return fingerprint_ != 0U && numeric.valid_for(symbolic) &&
         symbolic_fingerprint_ == symbolic.fingerprint_ &&
         symbolic_generation_ == symbolic.generation_ &&
         numeric_fingerprint_ == numeric.fingerprint_ &&
         numeric_generation_ == numeric.generation_;
}

Status make_linear_workspace_requirements(
    LinearAlgorithm algorithm, Int3 maximum_shape, std::uint8_t ghost_width,
    std::uint32_t maximum_restart, ReductionMode reduction_mode,
    RevisionToken execution_revision,
    LinearWorkspaceRequirements& out) noexcept {
  if (!valid_algorithm(algorithm) || !valid_reduction(reduction_mode) ||
      !valid_shape(maximum_shape) || ghost_width > 16U ||
      execution_revision == 0U ||
      (algorithm == LinearAlgorithm::fgmres && maximum_restart == 0U) ||
      (algorithm != LinearAlgorithm::fgmres && maximum_restart != 0U)) {
    return {StatusCode::invalid_plan, kLinearWorkspace};
  }
  std::size_t local_cells = 0U;
  if (!detail::checked_linear_multiply(
          static_cast<std::size_t>(maximum_shape.x),
          static_cast<std::size_t>(maximum_shape.y), local_cells) ||
      !detail::checked_linear_multiply(
          local_cells, static_cast<std::size_t>(maximum_shape.z),
          local_cells)) {
    return {StatusCode::invalid_plan, kLinearWorkspace};
  }
  std::size_t vector_slots =
      algorithm == LinearAlgorithm::pcg ? 5U : 10U;
  std::size_t scalar_doubles =
      algorithm == LinearAlgorithm::pcg ? 16U : 32U;
  std::size_t reduction_capacity = 2U;
  if (algorithm == LinearAlgorithm::fgmres) {
    if (!detail::checked_linear_multiply(
            static_cast<std::size_t>(maximum_restart), 2U, vector_slots) ||
        !detail::checked_linear_add(vector_slots, 8U, vector_slots) ||
        vector_slots > std::numeric_limits<std::uint8_t>::max()) {
      return {StatusCode::invalid_plan, kLinearWorkspace};
    }
    std::size_t restart_plus_one = 0U;
    std::size_t hessenberg = 0U;
    std::size_t vectors = 0U;
    if (!detail::checked_linear_add(
            static_cast<std::size_t>(maximum_restart), 1U,
            restart_plus_one) ||
        !detail::checked_linear_multiply(
            restart_plus_one, static_cast<std::size_t>(maximum_restart),
            hessenberg) ||
        !detail::checked_linear_multiply(
            static_cast<std::size_t>(maximum_restart), 4U, vectors) ||
        !detail::checked_linear_add(hessenberg, vectors, scalar_doubles) ||
        !detail::checked_linear_add(scalar_doubles, restart_plus_one,
                                    scalar_doubles) ||
        !detail::checked_linear_add(scalar_doubles, 16U,
                                    scalar_doubles)) {
      return {StatusCode::invalid_plan, kLinearWorkspace};
    }
    reduction_capacity =
        restart_plus_one < kLinearRecycleMaximumDirections
            ? kLinearRecycleMaximumDirections
            : restart_plus_one;
  }
  std::size_t maximum_vector_doubles = 0U;
  if (!detail::checked_linear_multiply(local_cells, vector_slots,
                                       maximum_vector_doubles) ||
      maximum_vector_doubles >
          static_cast<std::size_t>(std::numeric_limits<std::ptrdiff_t>::max()) /
              sizeof(double) ||
      scalar_doubles >
          static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max())) {
    return {StatusCode::invalid_plan, kLinearWorkspace};
  }
  std::uint64_t hash = detail::kLinearFnvOffset;
  hash = detail::linear_hash_mix(hash, static_cast<std::uint64_t>(algorithm));
  hash = detail::linear_hash_mix(hash,
                                 static_cast<std::uint64_t>(maximum_shape.x));
  hash = detail::linear_hash_mix(hash,
                                 static_cast<std::uint64_t>(maximum_shape.y));
  hash = detail::linear_hash_mix(hash,
                                 static_cast<std::uint64_t>(maximum_shape.z));
  hash = detail::linear_hash_mix(hash, ghost_width);
  hash = detail::linear_hash_mix(hash, maximum_restart);
  hash = detail::linear_hash_mix(hash,
                                 static_cast<std::uint64_t>(reduction_mode));
  hash = detail::linear_hash_mix(hash, execution_revision);
  hash = detail::linear_hash_mix(hash, reduction_capacity);
  LinearWorkspaceRequirements candidate;
  candidate.algorithm = algorithm;
  candidate.maximum_shape = maximum_shape;
  candidate.ghost_width = ghost_width;
  candidate.maximum_restart = maximum_restart;
  candidate.reduction_mode = reduction_mode;
  candidate.execution_revision = execution_revision;
  candidate.vector_slots = static_cast<std::uint8_t>(vector_slots);
  candidate.scalar_doubles = scalar_doubles;
  candidate.reduction_capacity = reduction_capacity;
  candidate.fingerprint = detail::finish_linear_hash(hash);
  out = candidate;
  return {};
}

Status register_linear_workspace(
    FieldRegistry& registry, std::string_view stable_prefix,
    const LinearWorkspaceRequirements& requirements,
    LinearWorkspaceFieldIds& out) noexcept {
  if (!detail::valid_workspace_requirements(requirements) ||
      stable_prefix.empty()) {
    return {StatusCode::invalid_plan, kLinearRegistry};
  }
  try {
    std::string vector_name(stable_prefix);
    std::string scalar_name(stable_prefix);
    vector_name += ".vectors";
    scalar_name += ".scalars";
    FieldRegistry candidate = registry;
    FieldId vectors = 0U;
    FieldId scalars = 0U;
    const Status declared_vectors = candidate.declare_field(
        vector_name, requirements.vector_slots, requirements.ghost_width,
        vectors);
    if (!declared_vectors) {
      return declared_vectors;
    }
    const Status declared_scalars =
        candidate.declare_field(scalar_name, 1U, 0U, scalars);
    if (!declared_scalars) {
      return declared_scalars;
    }
    registry = std::move(candidate);
    out = LinearWorkspaceFieldIds{vectors, scalars};
    return {};
  } catch (const std::bad_alloc&) {
    return {StatusCode::allocation_failure, 0U};
  } catch (...) {
    return {StatusCode::invalid_plan, kLinearRegistry};
  }
}

Status SolverWorkspace::bind(
    const LinearWorkspaceRequirements& requirements, FieldView vector_bundle,
    FieldView scalar_buffer, SolverWorkspace& out,
    LinearLifecycleCounters* counters) noexcept {
  if (!detail::valid_workspace_requirements(requirements) ||
      !sufficient_workspace_views(requirements, vector_bundle,
                                  scalar_buffer)) {
    return {StatusCode::invalid_plan, kLinearWorkspace};
  }
  if (out.fingerprint_ != 0U) {
    if (detail::same_workspace_requirements(out.requirements_, requirements) &&
        detail::same_field_binding(out.vector_bundle_, vector_bundle) &&
        detail::same_field_binding(out.scalar_buffer_, scalar_buffer)) {
      return {};
    }
    return {StatusCode::invalid_plan, kLinearWorkspace};
  }
  std::uint64_t counter = 0U;
  if (counters != nullptr &&
      !detail::checked_counter_increment(counters->workspace_bindings,
                                         counter)) {
    return {StatusCode::invalid_plan, kLinearCounter};
  }
  SolverWorkspace candidate;
  candidate.requirements_ = requirements;
  candidate.vector_bundle_ = vector_bundle;
  candidate.scalar_buffer_ = scalar_buffer;
  for (std::size_t slot = 0U; slot < requirements.vector_slots; ++slot) {
    candidate.vector_revisions_[slot] = static_cast<RevisionToken>(slot) + 1U;
  }
  candidate.next_vector_revision_ =
      static_cast<RevisionToken>(requirements.vector_slots) + 1U;
  candidate.fingerprint_ = detail::workspace_fingerprint(
      requirements, vector_bundle, scalar_buffer);
  candidate.binding_identity_ = 1U;
  candidate.recycle_ = {};
  out = candidate;
  if (counters != nullptr) {
    counters->workspace_bindings = counter;
  }
  return {};
}

bool SolverWorkspace::valid_for(
    const LinearWorkspaceRequirements& requirements) const noexcept {
  return fingerprint_ != 0U &&
         detail::same_workspace_requirements(requirements_, requirements);
}

FieldView SolverWorkspace::vector(std::uint8_t slot,
                                  Int3 active_shape) const noexcept {
  if (fingerprint_ == 0U || slot >= requirements_.vector_slots ||
      !valid_shape(active_shape) ||
      active_shape.x > requirements_.maximum_shape.x ||
      active_shape.y > requirements_.maximum_shape.y ||
      active_shape.z > requirements_.maximum_shape.z) {
    return {};
  }
  FieldView result = vector_bundle_;
  result.base += static_cast<std::size_t>(slot) * result.component_stride;
  result.interior = active_shape;
  result.components = 1U;
  result.revision = vector_revisions_[slot];
  return result;
}

Status SolverWorkspace::revise_vector(std::uint8_t slot) noexcept {
  if (fingerprint_ == 0U || slot >= requirements_.vector_slots ||
      vector_revisions_[slot] == 0U ||
      next_vector_revision_ == 0U ||
      next_vector_revision_ ==
          std::numeric_limits<RevisionToken>::max()) {
    return {StatusCode::invalid_plan, kLinearWorkspace};
  }
  vector_revisions_[slot] = next_vector_revision_++;
  return {};
}

bool SolverWorkspace::overlaps_storage(ConstFieldView view) const noexcept {
  if (fingerprint_ == 0U || view.base == nullptr) {
    return false;
  }
  return detail::field_views_overlap(view, as_const(vector_bundle_)) ||
         detail::field_views_overlap(view, as_const(scalar_buffer_));
}

bool SolverWorkspace::overlaps_storage(FieldView view) const noexcept {
  return overlaps_storage(as_const(view));
}

Span<double> SolverWorkspace::scalars(std::size_t offset,
                                      std::size_t count) const noexcept {
  std::size_t end = 0U;
  if (fingerprint_ == 0U || count == 0U ||
      !detail::checked_linear_add(offset, count, end) ||
      end > requirements_.scalar_doubles) {
    return {};
  }
  return {scalar_buffer_.base + offset, count};
}

Status SolverWorkspace::recycle_begin_capture(
    Int3 shape, PlanFingerprint source_identity) noexcept {
  if (recycle_.capture_active || recycle_.projection_pending ||
      requirements_.algorithm != LinearAlgorithm::fgmres ||
      !valid_shape(shape) || source_identity == 0U ||
      !shape_contains(requirements_.maximum_shape, shape) ||
      requirements_.maximum_restart == 0U ||
      requirements_.vector_slots <
          2U * static_cast<std::size_t>(requirements_.maximum_restart) +
              8U) {
    return {StatusCode::invalid_plan, kLinearWorkspace};
  }
  recycle_ = {};
  recycle_.capture_active = true;
  recycle_.shape = shape;
  recycle_.source_identity = source_identity;
  recycle_.snapshot_slot = static_cast<std::uint8_t>(
      2U * static_cast<std::size_t>(requirements_.maximum_restart) + 3U);
  return {};
}

Status SolverWorkspace::recycle_begin_projection(
    Int3 shape, PlanFingerprint current_identity) noexcept {
  if (!recycle_.capture_active || recycle_.cycle_active ||
      recycle_.projection_pending ||
      !valid_shape(shape) || !same_shape(shape, recycle_.shape) ||
      current_identity == 0U || recycle_.source_identity == 0U) {
    return {StatusCode::invalid_plan, kLinearWorkspace};
  }
  recycle_.current_identity = current_identity;
  recycle_.capture_active = false;
  recycle_.capture_vector_passes = 0U;
  recycle_.capture_cycle_attempts = 0U;
  recycle_.capture_reduction_calls = 0U;
  recycle_.capture_blocking_operations = 0U;
  recycle_.published_cycle_corrections = 0U;
  recycle_.projection_pending = true;
  return {};
}

void SolverWorkspace::recycle_clear() noexcept { recycle_ = {}; }

std::uint8_t SolverWorkspace::recycle_snapshot_slot() const noexcept {
  if (requirements_.algorithm != LinearAlgorithm::fgmres ||
      requirements_.maximum_restart == 0U) {
    return 0U;
  }
  if (recycle_.snapshot_slot != 0U) return recycle_.snapshot_slot;
  return static_cast<std::uint8_t>(
      2U * static_cast<std::size_t>(requirements_.maximum_restart) + 3U);
}

std::uint8_t SolverWorkspace::recycle_pool_slot(std::size_t index) const
    noexcept {
  if (index >= kLinearRecycleMaximumDirections + 1U ||
      requirements_.algorithm != LinearAlgorithm::fgmres ||
      requirements_.maximum_restart == 0U) {
    return 0U;
  }
  return static_cast<std::uint8_t>(
      2U * static_cast<std::size_t>(requirements_.maximum_restart) + 3U +
      index);
}

std::uint8_t SolverWorkspace::recycle_correction_logical_slot(
    std::size_t index) const noexcept {
  if ((!recycle_.capture_active && !recycle_.projection_pending) ||
      index >= recycle_.correction_count) {
    return 0U;
  }
  const std::size_t logical =
      (static_cast<std::size_t>(recycle_.oldest_correction) + index) %
      kLinearRecycleMaximumDirections;
  return recycle_.correction_order[logical];
}

FieldView SolverWorkspace::recycle_snapshot(Int3 shape) const noexcept {
  return vector(recycle_snapshot_slot(), shape);
}

FieldView SolverWorkspace::recycle_correction_storage(
    std::size_t index, Int3 shape) const noexcept {
  if ((!recycle_.capture_active && !recycle_.projection_pending) ||
      index >= recycle_.correction_count) {
    return {};
  }
  const std::size_t logical =
      (static_cast<std::size_t>(recycle_.oldest_correction) + index) %
      kLinearRecycleMaximumDirections;
  return vector(recycle_.correction_order[logical], shape);
}

ConstFieldView SolverWorkspace::recycle_correction(
    std::size_t index, Int3 shape) const noexcept {
  if ((!recycle_.capture_active && !recycle_.projection_pending) ||
      index >= recycle_.correction_count) {
    return {};
  }
  const std::size_t logical =
      (static_cast<std::size_t>(recycle_.oldest_correction) + index) %
      kLinearRecycleMaximumDirections;
  return as_const(vector(recycle_.correction_order[logical], shape));
}

Status SolverWorkspace::recycle_capture_cycle_start(
    ConstFieldView solution, ReductionEngine& reductions,
    Status prerequisite) noexcept {
  Status local_status{};
  FieldView snapshot{};
  if (!recycle_.capture_active || recycle_.cycle_active ||
      !same_shape(solution.interior, recycle_.shape) ||
      solution.base == nullptr) {
    local_status = {StatusCode::invalid_plan, kLinearWorkspace};
  } else {
    snapshot = recycle_snapshot(recycle_.shape);
    if (snapshot.base == nullptr ||
        detail::field_views_overlap(solution, as_const(snapshot))) {
      local_status = {StatusCode::invalid_plan, kLinearWorkspace};
    } else {
      local_status = copy_recycle_field(solution, snapshot);
      // A finite check can fail after a partial local copy.  Every rank with
      // the valid session/view still advances the same private slot revision;
      // otherwise a rank-selective failure would poison the next Halo use of
      // this workspace even after the session is discarded.
      const Status revised = revise_vector(recycle_snapshot_slot());
      if (local_status) local_status = revised;
    }
  }
  if (!prerequisite) local_status = prerequisite;
  const Status agreed_revised = reductions.consensus(local_status);
  if (!agreed_revised) {
    recycle_clear();
    return agreed_revised;
  }
  if (recycle_.capture_blocking_operations ==
      std::numeric_limits<std::uint64_t>::max()) {
    recycle_clear();
    return {StatusCode::invalid_plan, kLinearWorkspace};
  }
  ++recycle_.capture_blocking_operations;
  recycle_.cycle_active = true;
  return {};
}

Status SolverWorkspace::recycle_capture_cycle_publish(
    ConstFieldView solution, ReductionEngine& reductions) noexcept {
  const std::uint64_t reduction_begin = reductions.counters().calls;
  auto account_capture_reductions = [&]() noexcept -> Status {
    const std::uint64_t final_calls = reductions.counters().calls;
    if (final_calls < reduction_begin ||
        final_calls - reduction_begin >
            std::numeric_limits<std::uint64_t>::max() -
                recycle_.capture_reduction_calls) {
      return {StatusCode::invalid_plan, kLinearWorkspace};
    }
    recycle_.capture_reduction_calls += final_calls - reduction_begin;
    return {};
  };
  if (!recycle_.capture_active || !recycle_.cycle_active ||
      !same_shape(solution.interior, recycle_.shape) ||
      solution.base == nullptr || recycle_.capture_vector_passes >
                                      std::numeric_limits<std::uint64_t>::max() -
                                          2U) {
    recycle_clear();
    return {StatusCode::invalid_plan, kLinearWorkspace};
  }
  FieldView snapshot = recycle_snapshot(recycle_.shape);
  Status local_status{};
  if (snapshot.base == nullptr ||
      detail::field_views_overlap(solution, as_const(snapshot))) {
    local_status = {StatusCode::invalid_plan, kLinearWorkspace};
  }
  bool local_nonzero = false;
  if (local_status) {
    // The snapshot is also the only free pool slot.  Compute the correction
    // in place so every publish uses exactly one difference pass, including
    // a zero correction; the retained ring remains untouched until the
    // global nonzero decision below.
    local_status = difference_recycle_field(
        solution, as_const(snapshot), snapshot, local_nonzero);
    // As at cycle start, the in-place difference can fail after a partial
    // write.  Keep the private vector revision progression collective even
    // when the numerical status is not successful.
    const Status revised = revise_vector(recycle_.snapshot_slot);
    if (local_status) local_status = revised;
  }
  const double local_nonzero_value = local_nonzero ? 1.0 : 0.0;
  double global_nonzero = 0.0;
  const Status nonzero_status = reductions.checked_max(
      {&local_nonzero_value, 1U}, {&global_nonzero, 1U}, local_status);
  if (!nonzero_status) {
    recycle_clear();
    return nonzero_status;
  }
  if (recycle_.capture_cycle_attempts ==
      std::numeric_limits<std::uint64_t>::max()) {
    recycle_clear();
    return {StatusCode::invalid_plan, kLinearWorkspace};
  }
  ++recycle_.capture_cycle_attempts;
  if (recycle_.capture_blocking_operations ==
      std::numeric_limits<std::uint64_t>::max()) {
    recycle_clear();
    return {StatusCode::invalid_plan, kLinearWorkspace};
  }
  ++recycle_.capture_blocking_operations;
  if (recycle_.capture_vector_passes >
      std::numeric_limits<std::uint64_t>::max() - 2U) {
    recycle_clear();
    return {StatusCode::invalid_plan, kLinearWorkspace};
  }
  recycle_.capture_vector_passes += 2U;
  if (!(global_nonzero > 0.0) || !std::isfinite(global_nonzero)) {
    recycle_.cycle_active = false;
    const Status accounted = account_capture_reductions();
    if (!accounted) recycle_clear();
    return accounted;
  }
  const std::uint8_t published_slot = recycle_.snapshot_slot;
  std::uint8_t next_snapshot_slot = 0U;
  if (recycle_.correction_count < kLinearRecycleMaximumDirections) {
    for (std::size_t pool = 0U; pool < kLinearRecycleMaximumDirections + 1U;
         ++pool) {
      const std::uint8_t candidate = recycle_pool_slot(pool);
      bool used = candidate == published_slot;
      for (std::size_t index = 0U; index < recycle_.correction_count;
           ++index) {
        used = used || recycle_.correction_order[index] == candidate;
      }
      if (!used) {
        next_snapshot_slot = candidate;
        break;
      }
    }
    if (next_snapshot_slot == 0U) {
      (void)account_capture_reductions();
      recycle_clear();
      return {StatusCode::invalid_plan, kLinearWorkspace};
    }
  } else {
    next_snapshot_slot =
        recycle_.correction_order[recycle_.oldest_correction];
  }
  const Status accounted = account_capture_reductions();
  if (!accounted) {
    recycle_clear();
    return accounted;
  }
  if (recycle_.published_cycle_corrections ==
      std::numeric_limits<std::uint64_t>::max()) {
    recycle_clear();
    return {StatusCode::invalid_plan, kLinearWorkspace};
  }
  ++recycle_.published_cycle_corrections;
  if (recycle_.correction_count < kLinearRecycleMaximumDirections) {
    recycle_.correction_order[recycle_.correction_count] =
        published_slot;
    ++recycle_.correction_count;
  } else {
    recycle_.correction_order[recycle_.oldest_correction] =
        published_slot;
    recycle_.oldest_correction = static_cast<std::uint8_t>(
        (static_cast<std::size_t>(recycle_.oldest_correction) + 1U) %
        kLinearRecycleMaximumDirections);
  }
  recycle_.snapshot_slot = next_snapshot_slot;
  recycle_.cycle_active = false;
  return {};
}

void SolverWorkspace::recycle_set_projection_result(
    LinearSolveResult& result, std::uint64_t offered,
    std::uint64_t retained, std::uint64_t operator_applies,
    std::uint64_t reduction_calls, bool attempted, bool accepted,
    double projected_residual) noexcept {
  result.recycle_offered_directions = offered;
  result.recycle_retained_directions = retained;
  result.recycle_operator_applies = operator_applies;
  result.recycle_reduction_calls = reduction_calls;
  result.recycle_projection_attempted = attempted;
  result.recycle_projection_accepted = accepted;
  result.recycle_projected_true_residual = projected_residual;
}

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
Status SolverWorkspace::recycle_begin_capture_for_test(
    Int3 shape, PlanFingerprint source_identity) noexcept {
  return recycle_begin_capture(shape, source_identity);
}

Status SolverWorkspace::recycle_begin_projection_for_test(
    Int3 shape, PlanFingerprint current_identity) noexcept {
  return recycle_begin_projection(shape, current_identity);
}

Status SolverWorkspace::recycle_capture_cycle_start_for_test(
    ConstFieldView solution, ReductionEngine& reductions,
    Status prerequisite) noexcept {
  return recycle_capture_cycle_start(solution, reductions, prerequisite);
}

Status SolverWorkspace::recycle_capture_cycle_publish_for_test(
    ConstFieldView solution, ReductionEngine& reductions) noexcept {
  return recycle_capture_cycle_publish(solution, reductions);
}

void SolverWorkspace::recycle_clear_for_test() noexcept { recycle_clear(); }

std::size_t SolverWorkspace::recycle_correction_count_for_test()
    const noexcept {
  return recycle_correction_count();
}

std::uint64_t SolverWorkspace::recycle_capture_vector_passes_for_test()
    const noexcept {
  return recycle_capture_vector_passes();
}

std::uint64_t SolverWorkspace::recycle_capture_cycle_attempts_for_test()
    const noexcept {
  return recycle_capture_cycle_attempts();
}

std::uint64_t SolverWorkspace::recycle_capture_reduction_calls_for_test()
    const noexcept {
  return recycle_capture_reduction_calls();
}

std::uint64_t SolverWorkspace::recycle_capture_blocking_operations_for_test()
    const noexcept {
  return recycle_capture_blocking_operations();
}

std::uint64_t SolverWorkspace::recycle_capture_cycle_corrections_for_test()
    const noexcept {
  return recycle_capture_cycle_corrections();
}

std::uint8_t SolverWorkspace::recycle_snapshot_slot_for_test() const noexcept {
  return recycle_snapshot_slot();
}

std::uint8_t SolverWorkspace::recycle_correction_physical_slot_for_test(
    std::size_t index) const noexcept {
  return recycle_correction_logical_slot(index);
}

ConstFieldView SolverWorkspace::recycle_correction_for_test(
    std::size_t index, Int3 shape) const noexcept {
  return recycle_correction(index, shape);
}
#endif

LinearIdentity compose_linear_identity(
    const SymbolicPlan& symbolic, const NumericState& numeric,
    const HierarchyState& hierarchy,
    const SolverWorkspace& workspace) noexcept {
  if (!numeric.valid_for(symbolic) ||
      !hierarchy.valid_for(symbolic, numeric) || workspace.fingerprint_ == 0U) {
    return {};
  }
  std::uint64_t hash = detail::kLinearFnvOffset;
  hash = detail::linear_hash_mix(hash, symbolic.fingerprint_);
  hash = detail::linear_hash_mix(hash, numeric.fingerprint_);
  hash = detail::linear_hash_mix(hash, hierarchy.fingerprint_);
  hash = detail::linear_hash_mix(hash, workspace.fingerprint_);
  return {symbolic.fingerprint_, numeric.fingerprint_, hierarchy.fingerprint_,
          workspace.fingerprint_, detail::finish_linear_hash(hash)};
}

}  // namespace hundun::v04
