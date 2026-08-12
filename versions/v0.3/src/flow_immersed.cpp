// SPDX-License-Identifier: Apache-2.0

#include "hundun/flow_immersed.hpp"

#include "hundun/diag_immersed_module.hpp"

#include "flow_fixed_step_detail.hpp"
#include "flow_checkpoint_v3_detail.hpp"
#include "flow_immersed_density_detail.hpp"
#include "rt_mpi_error_detail.hpp"
#include "flow_immersed_piso_detail.hpp"
#include "flow_immersed_wale_detail.hpp"
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
#include "flow_immersed_access_detail.hpp"
#endif

#include "fvm_immersed_boundary_authority_detail.hpp"
#include "hundun/fvm_cell_centered.hpp"
#include "hundun/fvm_immersed_operator.hpp"
#include "hundun/fvm_immersed_reconstruction.hpp"
#include "hundun/ib_domain.hpp"
#include "hundun/ib_ghost_stencil_plan.hpp"
#include "hundun/ib_local_flow_pattern.hpp"
#include "hundun/lin_ghosted_vector.hpp"
#include "hundun/lin_ghosted_vector_halo.hpp"
#include "hundun/lin_preconditioners.hpp"
#include "hundun/lin_restarted_gmres.hpp"
#include "hundun/mesh_geometry.hpp"
#include "hundun/mesh_topology.hpp"
#include "hundun/rt_collective_status.hpp"
#include "hundun/rt_error.hpp"
#include "hundun/rt_field_access_plan.hpp"
#include "hundun/rt_field_registry.hpp"
#include "hundun/rt_field_storage.hpp"
#include "hundun/rt_halo_exchange.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_structured_decomposition.hpp"
#include "ib_wall_force_detail.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <functional>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <string>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace hundun::flow {
namespace {

constexpr std::uint64_t kFnvOffset = 1469598103934665603ULL;
constexpr std::uint64_t kFnvPrime = 1099511628211ULL;
constexpr runtime::PhaseId kStatePhase = 1800U;
constexpr runtime::ActorId kStateActor = 1800U;
constexpr std::uint32_t kCompactPressureRestartLength = 12U;

std::uint64_t hash_u64(std::uint64_t hash, std::uint64_t value) noexcept {
  for (std::size_t byte = 0U; byte < sizeof(value); ++byte) {
    hash ^= (value >> (byte * 8U)) & 0xffU;
    hash *= kFnvPrime;
  }
  return hash;
}

std::uint64_t fp64_bits(double value) noexcept;

std::uint64_t diagnostic_report_seal(
    const ImmersedFlowStepAttemptReport &attempt) noexcept {
  std::uint64_t hash = kFnvOffset;
  if (const auto *ideal =
          std::get_if<IdealGasStepAttemptReport>(&attempt.base)) {
    if (!detail::DensityClosureAdapter::report_authenticated(*ideal))
      return 0U;
    hash = hash_u64(hash, 0x494445414c474153ULL);
    hash = hash_u64(hash,
                    detail::DensityClosureAdapter::report_seal(*ideal));
  } else if (const auto *material =
          std::get_if<MaterialDensityStepAttemptReport>(&attempt.base)) {
    if (!detail::DensityClosureBridge::report_authenticated(*material))
      return 0U;
    hash = hash_u64(hash, 0x4d4154455249414cULL);
    hash = hash_u64(hash,
                    detail::DensityClosureBridge::report_seal(*material));
  } else {
    if (!std::holds_alternative<StepAttemptReport>(attempt.base))
      return 0U;
    const auto &report = std::get<StepAttemptReport>(attempt.base);
    hash = hash_u64(hash, static_cast<std::uint64_t>(report.disposition));
    hash = hash_u64(hash, static_cast<std::uint64_t>(report.reason));
    hash = hash_u64(
        hash, static_cast<std::uint64_t>(report.lowest_failing_rank + 1));
    hash = hash_u64(hash, report.pressure_corrector_count);
    hash = hash_u64(hash, fp64_bits(report.attempted_dt_s));
    hash = hash_u64(hash, fp64_bits(report.suggested_dt_s));
    const auto add_solve = [&](const linear::SolveReport &solve) {
      hash = hash_u64(hash, static_cast<std::uint64_t>(solve.reason));
      hash = hash_u64(hash, solve.iterations);
      hash = hash_u64(hash, fp64_bits(solve.initial_residual));
      hash = hash_u64(hash, fp64_bits(solve.recursive_residual));
      hash = hash_u64(hash, fp64_bits(solve.final_residual));
      hash = hash_u64(hash, solve.matvec_count);
      hash = hash_u64(hash, solve.preconditioner_apply_count);
      hash = hash_u64(hash, solve.global_reduction_count);
      hash = hash_u64(
          hash, static_cast<std::uint64_t>(solve.lowest_failing_rank + 1));
    };
    for (const auto &solve : report.momentum.components)
      add_solve(solve);
    for (const auto &solve : report.pressure)
      add_solve(solve);
    hash = hash_u64(hash, fp64_bits(report.final_continuity_normalized_l2));
    hash = hash_u64(hash, fp64_bits(report.final_pressure_residual_l2));
    for (const auto value : report.final_momentum_normalized_l2)
      hash = hash_u64(hash, fp64_bits(value));
    hash = hash_u64(hash, report.final_transport_normalized_l2.size());
    for (const auto value : report.final_transport_normalized_l2)
      hash = hash_u64(hash, fp64_bits(value));
    hash = hash_u64(
        hash, fp64_bits(report.final_mass_relative_conservation_defect));
    for (const auto value : report.final_momentum_relative_conservation_defect)
      hash = hash_u64(hash, fp64_bits(value));
    hash = hash_u64(
        hash, report.final_transport_relative_conservation_defect.size());
    for (const auto value : report.final_transport_relative_conservation_defect)
      hash = hash_u64(hash, fp64_bits(value));
    hash =
        hash_u64(hash, report.final_backflow_evidence.has_value() ? 1U : 0U);
    if (report.final_backflow_evidence.has_value()) {
      const auto &evidence = *report.final_backflow_evidence;
      hash = hash_u64(hash, evidence.patch_id);
      hash = hash_u64(hash, evidence.step);
      hash = hash_u64(hash, fp64_bits(evidence.time_s));
      hash = hash_u64(
          hash, fp64_bits(evidence.minimum_outward_mass_flux_kg_per_s));
      hash = hash_u64(hash, evidence.global_face_id);
      hash = hash_u64(
          hash, static_cast<std::uint64_t>(evidence.lowest_failing_rank + 1));
    }
  }
  if (attempt.linear_solve_failure.has_value()) {
    const auto &failure = *attempt.linear_solve_failure;
    hash = hash_u64(hash, 0x4c494e4641494c55ULL);
    hash = hash_u64(hash, static_cast<std::uint64_t>(failure.phase));
    hash = hash_u64(hash, failure.pressure_corrector_index);
    hash = hash_u64(hash, failure.component_index);
    hash = hash_u64(hash, static_cast<std::uint64_t>(failure.solve.reason));
    hash = hash_u64(hash, failure.solve.iterations);
    hash = hash_u64(hash, fp64_bits(failure.solve.initial_residual));
    hash = hash_u64(hash, fp64_bits(failure.solve.recursive_residual));
    hash = hash_u64(hash, fp64_bits(failure.solve.final_residual));
    hash = hash_u64(hash, failure.solve.matvec_count);
    hash = hash_u64(hash, failure.solve.preconditioner_apply_count);
    hash = hash_u64(hash, failure.solve.global_reduction_count);
    hash = hash_u64(
        hash,
        static_cast<std::uint64_t>(failure.solve.lowest_failing_rank + 1));
    hash = hash_u64(hash, fp64_bits(failure.independent_residual_l2));
    hash = hash_u64(hash, fp64_bits(failure.rhs_l2));
    hash = hash_u64(hash, fp64_bits(failure.acceptance_threshold));
  }
  hash = hash_u64(hash, attempt.force.has_value() ? 1U : 0U);
  if (attempt.force.has_value()) {
    const auto add = [&](const immersed::ForceComponents &force) {
      for (const auto value : {force.pressure_N, force.total_N,
                               force.viscous_N}) {
        hash = hash_u64(hash, fp64_bits(value.x));
        hash = hash_u64(hash, fp64_bits(value.y));
        hash = hash_u64(hash, fp64_bits(value.z));
      }
    };
    add(attempt.force->operator_force);
    add(attempt.force->budget_reaction);
    add(attempt.force->surface_traction);
    add(attempt.force->consistency);
  }
  hash = hash_u64(hash, attempt.wale.has_value() ? 1U : 0U);
  if (attempt.wale.has_value()) {
    hash = hash_u64(hash, attempt.wale->identity.value);
    hash = hash_u64(hash,
                    fp64_bits(attempt.wale->minimum_nu_t_m2_per_s));
    hash = hash_u64(hash,
                    fp64_bits(attempt.wale->maximum_nu_t_m2_per_s));
    hash = hash_u64(hash, fp64_bits(attempt.wale->l2_nu_t_m2_per_s));
    hash = hash_u64(hash, attempt.wale->exact_zero_count);
    hash = hash_u64(hash, attempt.wale->owned_active_count);
  }
  return hash == 0U ? 1U : hash;
}

std::uint64_t diagnostic_snapshot_seal(
    const ImmersedFlowStepAttemptReport &attempt,
    config::DensityModel density_model,
    const std::optional<immersed::WallForceSample> &wall_force,
    const std::optional<finite_volume::ImmersedOperatorReport>
        &local_flow_pattern) noexcept {
  std::uint64_t hash = diagnostic_report_seal(attempt);
  if (hash == 0U)
    return 0U;
  hash = hash_u64(hash, static_cast<std::uint64_t>(density_model));
  hash = hash_u64(hash, wall_force.has_value() ? 1U : 0U);
  if (wall_force.has_value()) {
    const auto add_real3 = [&](runtime::Real3 value) {
      hash = hash_u64(hash, fp64_bits(value.x));
      hash = hash_u64(hash, fp64_bits(value.y));
      hash = hash_u64(hash, fp64_bits(value.z));
    };
    for (const auto value : {
             wall_force->surface_traction.pressure_N,
             wall_force->surface_traction.total_N,
             wall_force->surface_traction.viscous_N,
             wall_force->moment_about_global_origin.pressure_N_m,
             wall_force->moment_about_global_origin.total_N_m,
             wall_force->moment_about_global_origin.viscous_N_m,
             wall_force->area_vector_closure_m2,
         })
      add_real3(value);
    hash = hash_u64(hash, wall_force->quadrature_point_count);
    hash = hash_u64(
        hash, static_cast<std::uint64_t>(wall_force->lowest_failing_rank + 1));
  }
  hash = hash_u64(hash, local_flow_pattern.has_value() ? 1U : 0U);
  if (local_flow_pattern.has_value()) {
    hash = hash_u64(hash, local_flow_pattern->row_fingerprint);
    hash = hash_u64(hash, local_flow_pattern->replacement_group_count);
    hash = hash_u64(hash,
                    local_flow_pattern->algebraic_occurrence_count);
    hash = hash_u64(
        hash, fp64_bits(local_flow_pattern->replacement_coefficient_l2));
    hash = hash_u64(hash, local_flow_pattern->limiting_case_status);
  }
  return hash == 0U ? 1U : hash;
}

bool finite(runtime::Real3 value) noexcept {
  return std::isfinite(value.x) && std::isfinite(value.y) &&
         std::isfinite(value.z);
}

double norm_squared(runtime::Real3 value) noexcept {
  return value.x * value.x + value.y * value.y + value.z * value.z;
}

double dot(runtime::Real3 left, runtime::Real3 right) noexcept {
  return left.x * right.x + left.y * right.y + left.z * right.z;
}

runtime::Real3 multiply(double scale, runtime::Real3 value) noexcept {
  return {scale * value.x, scale * value.y, scale * value.z};
}

runtime::Real3 subtract(runtime::Real3 left, runtime::Real3 right) noexcept {
  return {left.x - right.x, left.y - right.y, left.z - right.z};
}

ForceAttemptReport assemble_candidate_force_report_from_budget(
    const immersed::ForceComponents &budget_reaction,
    const immersed::ForceComponents &surface_traction) {
  const auto negate = [](runtime::Real3 value) noexcept {
    return runtime::Real3{-value.x, -value.y, -value.z};
  };
  ForceAttemptReport result;
  result.budget_reaction = budget_reaction;
  result.operator_force.pressure_N =
      negate(result.budget_reaction.pressure_N);
  result.operator_force.viscous_N =
      negate(result.budget_reaction.viscous_N);
  result.operator_force.total_N = negate(result.budget_reaction.total_N);
  result.surface_traction = surface_traction;
  result.consistency.pressure_N =
      subtract(result.operator_force.pressure_N,
               result.surface_traction.pressure_N);
  result.consistency.viscous_N =
      subtract(result.operator_force.viscous_N,
               result.surface_traction.viscous_N);
  result.consistency.total_N =
      subtract(result.operator_force.total_N,
               result.surface_traction.total_N);
  return result;
}

double component(runtime::Real3 value, std::size_t index) noexcept {
  return index == 0U ? value.x : index == 1U ? value.y : value.z;
}

std::size_t checked_bytes(std::size_t count) {
  if (count > std::numeric_limits<std::size_t>::max() / sizeof(double))
    throw runtime::Error("immersed-flow vector byte count overflows");
  return count * sizeof(double);
}

std::uint64_t fp64_bits(double value) noexcept {
  std::uint64_t bits{};
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

std::uint64_t geometry_identity(const mesh::MeshTopology &topology,
                                const mesh::MeshGeometry &geometry) {
  std::uint64_t hash = kFnvOffset;
  hash = hash_u64(hash, static_cast<std::uint64_t>(geometry.mapping_kind()));
  const auto extent = topology.global_extent();
  hash = hash_u64(hash, static_cast<std::uint64_t>(extent.x));
  hash = hash_u64(hash, static_cast<std::uint64_t>(extent.y));
  hash = hash_u64(hash, static_cast<std::uint64_t>(extent.z));
  for (int k = 0; k <= extent.z; ++k) {
    for (int j = 0; j <= extent.y; ++j) {
      for (int i = 0; i <= extent.x; ++i) {
        const auto vertex = geometry.vertex_position_m({i, j, k});
        hash = hash_u64(hash, fp64_bits(vertex.x));
        hash = hash_u64(hash, fp64_bits(vertex.y));
        hash = hash_u64(hash, fp64_bits(vertex.z));
      }
    }
  }
  return hash;
}

struct ActiveConnection final {
  std::size_t row{};
  mesh::LocalFaceId face{};
  mesh::GlobalCellId neighbour{};
  double diagonal_coefficient{};
  double neighbour_coefficient{};
  double mass_flux_coefficient{};
  bool pressure_reference{};
  bool local_neighbour{};
  std::size_t value_offset{std::numeric_limits<std::size_t>::max()};
  int neighbour_rank{-1};
};

struct ExactPredictorResponse final {
  struct Work final {
    std::uint64_t response_count{};
    std::uint64_t solve_count{};
    std::uint64_t iteration_count{};
    std::uint64_t matvec_count{};
    std::uint64_t preconditioner_apply_count{};
    std::uint64_t global_reduction_count{};
  } work;
  std::vector<double> divergence_per_volume;
  std::array<std::vector<double>, 3> velocity_increment;
  std::vector<double> face_mass_flux_increment;
  std::vector<double> face_velocity_increment;
};
static_assert(
    std::is_nothrow_copy_assignable_v<ExactPredictorResponse::Work>);
static_assert(std::is_nothrow_copy_assignable_v<les::WaleSummary>);

std::uint64_t checked_exact_predictor_work_sum(std::uint64_t left,
                                               std::uint64_t right) {
  if (right > std::numeric_limits<std::uint64_t>::max() - left)
    throw runtime::Error("immersed-flow exact predictor work counter would wrap");
  return left + right;
}

void add_exact_predictor_work(ExactPredictorResponse::Work &aggregate,
                              const ExactPredictorResponse::Work &increment) {
  aggregate.response_count = checked_exact_predictor_work_sum(
      aggregate.response_count, increment.response_count);
  aggregate.solve_count = checked_exact_predictor_work_sum(
      aggregate.solve_count, increment.solve_count);
  aggregate.iteration_count = checked_exact_predictor_work_sum(
      aggregate.iteration_count, increment.iteration_count);
  aggregate.matvec_count = checked_exact_predictor_work_sum(
      aggregate.matvec_count, increment.matvec_count);
  aggregate.preconditioner_apply_count = checked_exact_predictor_work_sum(
      aggregate.preconditioner_apply_count,
      increment.preconditioner_apply_count);
  aggregate.global_reduction_count = checked_exact_predictor_work_sum(
      aggregate.global_reduction_count, increment.global_reduction_count);
}

void add_exact_predictor_solve(ExactPredictorResponse::Work &aggregate,
                               const linear::SolveReport &solve) {
  ExactPredictorResponse::Work increment;
  increment.solve_count = 1U;
  increment.iteration_count = solve.iterations;
  increment.matvec_count = solve.matvec_count;
  increment.preconditioner_apply_count = solve.preconditioner_apply_count;
  increment.global_reduction_count = solve.global_reduction_count;
  add_exact_predictor_work(aggregate, increment);
}

void resize_exact_predictor_response(ExactPredictorResponse &response,
                                     std::size_t cell_count,
                                     std::size_t face_count) {
  response.divergence_per_volume.resize(cell_count);
  for (auto &component : response.velocity_increment)
    component.resize(cell_count);
  response.face_mass_flux_increment.resize(face_count);
  if (face_count > std::numeric_limits<std::size_t>::max() / 3U)
    throw runtime::Error(
        "immersed-flow exact predictor face response size overflows");
  response.face_velocity_increment.resize(face_count * 3U);
}

bool exact_predictor_response_layout_matches(
    const ExactPredictorResponse &response, std::size_t cell_count,
    std::size_t face_count) noexcept {
  return response.divergence_per_volume.size() == cell_count &&
         std::all_of(response.velocity_increment.begin(),
                     response.velocity_increment.end(),
                     [cell_count](const auto &component) {
                       return component.size() == cell_count;
                     }) &&
         response.face_mass_flux_increment.size() == face_count &&
         face_count <= std::numeric_limits<std::size_t>::max() / 3U &&
         response.face_velocity_increment.size() == face_count * 3U;
}

void copy_exact_predictor_response(const ExactPredictorResponse &source,
                                   ExactPredictorResponse &destination) {
  const std::size_t cell_count = destination.divergence_per_volume.size();
  const std::size_t face_count = destination.face_mass_flux_increment.size();
  if (!exact_predictor_response_layout_matches(source, cell_count,
                                               face_count) ||
      !exact_predictor_response_layout_matches(destination, cell_count,
                                               face_count))
    throw runtime::Error(
        "immersed-flow exact predictor response layout is invalid");
  destination.work = source.work;
  std::copy(source.divergence_per_volume.begin(),
            source.divergence_per_volume.end(),
            destination.divergence_per_volume.begin());
  for (std::size_t component = 0U; component < 3U; ++component)
    std::copy(source.velocity_increment[component].begin(),
              source.velocity_increment[component].end(),
              destination.velocity_increment[component].begin());
  std::copy(source.face_mass_flux_increment.begin(),
            source.face_mass_flux_increment.end(),
            destination.face_mass_flux_increment.begin());
  std::copy(source.face_velocity_increment.begin(),
            source.face_velocity_increment.end(),
            destination.face_velocity_increment.begin());
}

void subtract_exact_predictor_response(
    const ExactPredictorResponse &source,
    const ExactPredictorResponse &affine_baseline,
    ExactPredictorResponse &destination) {
  const std::size_t cell_count = destination.divergence_per_volume.size();
  const std::size_t face_count = destination.face_mass_flux_increment.size();
  if (!exact_predictor_response_layout_matches(source, cell_count,
                                               face_count) ||
      !exact_predictor_response_layout_matches(affine_baseline, cell_count,
                                               face_count) ||
      !exact_predictor_response_layout_matches(destination, cell_count,
                                               face_count))
    throw runtime::Error(
        "immersed-flow exact predictor baseline layout is invalid");
  destination.work = source.work;
  for (std::size_t row = 0U; row < cell_count; ++row) {
    destination.divergence_per_volume[row] =
        source.divergence_per_volume[row] -
        affine_baseline.divergence_per_volume[row];
    for (std::size_t component = 0U; component < 3U; ++component)
      destination.velocity_increment[component][row] =
          source.velocity_increment[component][row] -
          affine_baseline.velocity_increment[component][row];
  }
  for (std::size_t face = 0U; face < face_count; ++face) {
    destination.face_mass_flux_increment[face] =
        source.face_mass_flux_increment[face] -
        affine_baseline.face_mass_flux_increment[face];
    for (std::size_t component = 0U; component < 3U; ++component) {
      const std::size_t offset = face * 3U + component;
      destination.face_velocity_increment[offset] =
          source.face_velocity_increment[offset] -
          affine_baseline.face_velocity_increment[offset];
    }
  }
}

class ActivePressureTwoLevelPreconditioner;

class ActivePressureOperator final : public linear::LinearOperator {
public:
  ActivePressureOperator(const mesh::MeshTopology &topology,
                         const mesh::MeshGeometry &geometry,
                         const runtime::StructuredDecomposition &decomposition,
                         const immersed::ActiveCellLayout &active,
                         const immersed::ActiveBoundaryLayout &active_boundary,
                         const boundary::BoundaryRegistry &boundaries,
                         const runtime::MpiContext &mpi,
                         execution::ExecutionContext &context)
      : topology_(&topology), geometry_(&geometry), mpi_(&mpi),
        context_(&context), decomposition_(&decomposition),
        layout_(active.owned_active_count(), active.ordered_global_ids()),
        exchange_values_(context, layout_),
        active_halo_(linear::GhostedVectorHalo::create(decomposition, topology,
                                                       context, layout_)),
        diagonal_(active.owned_active_count(), 0.0),
        owned_ids_(
            active.ordered_global_ids().begin(),
            active.ordered_global_ids().begin() +
                static_cast<std::ptrdiff_t>(active.owned_active_count())) {
    active_sqrt_volumes_.resize(layout_.local_count());
    active_offsets_.reserve(layout_.local_count());
    for (std::size_t index = 0U; index < layout_.local_count(); ++index) {
      const auto id = layout_.global_ids()[index];
      const auto local = topology_->find_local_cell(id);
      if (!local.has_value() ||
          (index < layout_.owned_count()) !=
              (topology_->cell_ownership(*local) ==
               mesh::EntityOwnership::owned))
        throw runtime::Error(
            "immersed-flow active pressure local layout is invalid");
      active_offsets_.emplace(id, index);
      active_sqrt_volumes_[index] =
          std::sqrt(geometry_->cell_volume_m3(*local));
    }
    if (std::any_of(active_sqrt_volumes_.begin(), active_sqrt_volumes_.end(),
                    [](double value) {
                      return !(value > 0.0) || !std::isfinite(value);
                    }))
      throw runtime::Error(
          "immersed-flow active pressure volume scaling is invalid");
    for (auto &values : active_momentum_diagonal_)
      values.resize(layout_.local_count());
    face_velocity_mobility_.resize(topology_->local_face_count() * 3U);
    face_mass_flux_coefficient_.resize(topology_->local_face_count());
    build_connections(active, active_boundary, boundaries);
    double reference_present = has_pressure_reference_ ? 1.0 : 0.0;
    mpi_->allreduce_fp64_in_place(&reference_present, 1U,
                                  runtime::Fp64ReductionOperation::maximum);
    has_pressure_reference_ = reference_present == 1.0;
    if (has_pressure_reference_ != active_boundary.has_pressure_reference())
      throw runtime::Error(
          "immersed-flow active pressure reference disagrees with its layout");
    bind_connections();
  }

  linear::VectorLayout domain_layout() const override { return layout_; }
  linear::VectorLayout range_layout() const override { return layout_; }
  const execution::ExecutionContext &context() const override {
    return *context_;
  }
  std::uint64_t revision() const override { return revision_; }

  void prepare_face_coefficients(
      double rho_ref,
      const std::array<std::vector<double>, 3> &momentum_diagonal) {
    if (!(rho_ref > 0.0) || !std::isfinite(rho_ref))
      throw runtime::Error("immersed-flow pressure mobility input is invalid");
    prepare_face_coefficients_impl(nullptr, rho_ref, momentum_diagonal);
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    test_density_ = rho_ref;
#endif
  }

  void prepare_face_coefficients(
      const std::vector<double> &face_density,
      const std::array<std::vector<double>, 3> &momentum_diagonal) {
    prepare_face_coefficients_impl(&face_density, 0.0, momentum_diagonal);
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    test_density_ = 0.0;
#endif
  }

  void prepare_face_coefficients_impl(
      const std::vector<double> *face_density, double rho_ref,
      const std::array<std::vector<double>, 3> &momentum_diagonal) {
    const bool local_layout_ok =
        momentum_diagonal[0].size() == owned_ids_.size() &&
        momentum_diagonal[1].size() == owned_ids_.size() &&
        momentum_diagonal[2].size() == owned_ids_.size() &&
        (face_density == nullptr ||
         face_density->size() == topology_->local_face_count());
    const auto layout_status = runtime::collective_status(
        *mpi_, local_layout_ok,
        "immersed-flow pressure mobility diagonal layout is invalid");
    if (!layout_status.ok)
      throw runtime::Error("immersed-flow pressure mobility input is invalid");
    for (std::size_t component_index = 0U; component_index < 3U;
         ++component_index) {
      auto owned = exchange_values_.owned_view();
      std::copy(momentum_diagonal[component_index].begin(),
                momentum_diagonal[component_index].end(), owned.data());
      active_halo_.exchange(exchange_values_);
      const auto local =
          static_cast<const linear::GhostedVector &>(exchange_values_)
              .local_view();
      std::copy(local.data(), local.data() + local.size(),
                active_momentum_diagonal_[component_index].begin());
    }
    for (const auto &values : active_momentum_diagonal_) {
      if (std::any_of(values.begin(), values.end(), [](double value) {
            return !(value > 0.0) || !std::isfinite(value);
          }))
        throw runtime::Error("immersed-flow pressure mobility diagonal is invalid");
    }
    std::fill(face_velocity_mobility_.begin(), face_velocity_mobility_.end(),
              0.0);
    std::fill(face_mass_flux_coefficient_.begin(),
              face_mass_flux_coefficient_.end(), 0.0);
    for (mesh::LocalFaceId face = 0U; face < topology_->local_face_count();
         ++face) {
      const auto owner = topology_->owner(face);
      const auto neighbour = topology_->neighbour(face);
      const auto owner_offset =
          active_offset_if_present(topology_->global_cell_id(owner));
      const auto neighbour_offset =
          neighbour.has_value()
              ? active_offset_if_present(topology_->global_cell_id(*neighbour))
              : kMissingGlobalOffset;
      if (owner_offset == kMissingGlobalOffset &&
          neighbour_offset == kMissingGlobalOffset)
        continue;
      const auto area =
          geometry_->face_area_vector_m2(face, mesh::FaceSide::owner);
      const double area_magnitude = std::sqrt(norm_squared(area));
      if (!(area_magnitude > 0.0) || !std::isfinite(area_magnitude))
        throw runtime::Error("immersed-flow pressure mobility face is invalid");
      const auto unit_normal = multiply(1.0 / area_magnitude, area);
      double owner_weight = 1.0;
      double neighbour_weight = 0.0;
      if (owner_offset != kMissingGlobalOffset &&
          neighbour_offset != kMissingGlobalOffset) {
        const auto owner_to_face = subtract(geometry_->face_center_m(face),
                                            geometry_->cell_center_m(owner));
        const auto face_to_neighbour =
            subtract(geometry_->face_displacement_m(face), owner_to_face);
        const double owner_distance = std::sqrt(norm_squared(owner_to_face));
        const double neighbour_distance =
            std::sqrt(norm_squared(face_to_neighbour));
        const double total_distance = owner_distance + neighbour_distance;
        if (!(total_distance > 0.0) || !std::isfinite(total_distance))
          throw runtime::Error("immersed-flow pressure mobility weights are invalid");
        owner_weight = neighbour_distance / total_distance;
        neighbour_weight = owner_distance / total_distance;
      } else if (owner_offset == kMissingGlobalOffset) {
        owner_weight = 0.0;
        neighbour_weight = 1.0;
      }
      double normal_mobility = 0.0;
      for (std::size_t component_index = 0U; component_index < 3U;
           ++component_index) {
        double rate = 0.0;
        if (owner_offset != kMissingGlobalOffset) {
          rate += owner_weight *
                  active_momentum_diagonal_[component_index][owner_offset] /
                  geometry_->cell_volume_m3(owner);
        }
        if (neighbour_offset != kMissingGlobalOffset) {
          rate += neighbour_weight *
                  active_momentum_diagonal_[component_index][neighbour_offset] /
                  geometry_->cell_volume_m3(*neighbour);
        }
        if (!(rate > 0.0) || !std::isfinite(rate))
          throw runtime::Error("immersed-flow pressure mobility rate is invalid");
        const double mobility = 1.0 / rate;
        face_velocity_mobility_[static_cast<std::size_t>(face) * 3U +
                                component_index] = mobility;
        const double normal_component = component(unit_normal, component_index);
        normal_mobility += normal_component * normal_component * mobility;
      }
      const double rho_face =
          face_density == nullptr ? rho_ref : (*face_density)[face];
      const double coefficient = rho_face * normal_mobility;
      if (!(coefficient > 0.0) || !std::isfinite(coefficient))
        throw runtime::Error("immersed-flow pressure face coefficient is invalid");
      face_mass_flux_coefficient_[face] = coefficient;
    }
    for (auto &connection : connections_) {
      connection.mass_flux_coefficient =
          face_mass_flux_coefficient_[connection.face];
      if (!(connection.mass_flux_coefficient > 0.0) ||
          !std::isfinite(connection.mass_flux_coefficient))
        throw runtime::Error(
            "immersed-flow pressure connection coefficient is invalid");
    }
    coefficients_prepared_ = true;
  }

  double face_velocity_mobility(mesh::LocalFaceId face,
                                std::size_t component_index) const {
    if (!coefficients_prepared_ || face >= topology_->local_face_count() ||
        component_index >= 3U)
      throw runtime::Error("immersed-flow pressure face mobility is unavailable");
    const double value =
        face_velocity_mobility_[static_cast<std::size_t>(face) * 3U +
                                component_index];
    if (!(value > 0.0) || !std::isfinite(value))
      throw runtime::Error("immersed-flow pressure face mobility is invalid");
    return value;
  }

  double face_mass_flux_coefficient(mesh::LocalFaceId face) const {
    if (!coefficients_prepared_ || face >= face_mass_flux_coefficient_.size())
      throw runtime::Error("immersed-flow pressure face coefficient is unavailable");
    const double value = face_mass_flux_coefficient_[face];
    if (!(value > 0.0) || !std::isfinite(value))
      throw runtime::Error("immersed-flow pressure face coefficient is invalid");
    return value;
  }

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  double cell_velocity_mobility(mesh::GlobalCellId cell,
                                std::size_t component_index) const {
    if (!coefficients_prepared_ || component_index >= 3U)
      throw runtime::Error("immersed-flow pressure cell mobility is unavailable");
    const auto offset = active_offset_if_present(cell);
    if (offset == kMissingGlobalOffset)
      throw runtime::Error("immersed-flow pressure cell mobility is unavailable");
    const double value = active_sqrt_volumes_[offset] *
                         active_sqrt_volumes_[offset] /
                         active_momentum_diagonal_[component_index][offset];
    if (!(value > 0.0) || !std::isfinite(value))
      throw runtime::Error("immersed-flow pressure cell mobility is invalid");
    return value;
  }

  double density() const noexcept { return test_density_; }
#endif

  runtime::HaloPerformanceCounters performance_counters() const noexcept {
    return active_halo_.performance_counters();
  }

  void replace(std::uint64_t dependency_revision) {
    if (!coefficients_prepared_ || dependency_revision == 0U) {
      throw runtime::Error("immersed-flow pressure coefficient is invalid");
    }
    revision_ = dependency_revision;
    std::fill(diagonal_.begin(), diagonal_.end(), 0.0);
    for (const auto &connection : connections_) {
      const double value =
          connection.mass_flux_coefficient * connection.diagonal_coefficient;
      if (!(value > 0.0) || !std::isfinite(value) ||
          !std::isfinite(diagonal_[connection.row] + value))
        throw runtime::Error("immersed-flow pressure diagonal is invalid");
      diagonal_[connection.row] += value;
    }
    if (std::any_of(diagonal_.begin(), diagonal_.end(), [](double value) {
          return !(value > 0.0) || !std::isfinite(value);
        })) {
      throw runtime::Error(
          "immersed-flow active pressure region has an unconstrained row");
    }
  }

  execution::ExecutionEvent
  apply(execution::VectorView<const double> x,
        execution::VectorView<double> y) const override {
    if (active_)
      throw runtime::Error("immersed-flow pressure operator is already active");
    struct Guard final {
      explicit Guard(bool &value) : value_(&value) { value = true; }
      ~Guard() { *value_ = false; }
      bool *value_;
    } guard(active_);
    if (x.size() != layout_.owned_count() ||
        y.size() != layout_.owned_count() || x.stride() != 1U ||
        y.stride() != 1U || !y.writable() ||
        x.backend_identity() != context_->backend_identity() ||
        y.backend_identity() != context_->backend_identity()) {
      throw runtime::Error("immersed-flow pressure operator view is incompatible");
    }
    const double *const x_values = x.data();
    double *const y_values = y.data();
    if (y.size() != 0U)
      std::fill(y_values, y_values + y.size(), 0.0);
    auto exchange_owned = exchange_values_.owned_view();
    std::copy(x_values, x_values + x.size(), exchange_owned.data());
    active_halo_.begin(exchange_values_);
    last_apply_schedule_ = 1U;
    const auto accumulate = [&](const ActiveConnection &connection,
                                double neighbour_value) {
      y_values[connection.row] +=
          connection.mass_flux_coefficient *
          (connection.diagonal_coefficient * x_values[connection.row] -
           connection.neighbour_coefficient * neighbour_value);
    };
    for (const auto &connection : connections_) {
      if (connection.pressure_reference) {
        accumulate(connection, 0.0);
        continue;
      }
      if (connection.local_neighbour)
        accumulate(connection, x_values[connection.value_offset]);
    }
    last_apply_schedule_ = last_apply_schedule_ * 10U + 2U;
    active_halo_.wait(exchange_values_);
    last_apply_schedule_ = last_apply_schedule_ * 10U + 3U;
    const auto exchanged_values =
        static_cast<const linear::GhostedVector &>(exchange_values_)
            .local_view();
    for (const auto &connection : connections_) {
      if (connection.pressure_reference || connection.local_neighbour)
        continue;
      if (connection.value_offset == kMissingGlobalOffset)
        throw runtime::Error("immersed-flow pressure neighbour value is unavailable");
      accumulate(connection, exchanged_values[connection.value_offset]);
    }
    last_apply_schedule_ = last_apply_schedule_ * 10U + 4U;
    for (std::size_t row = 0U; row < y.size(); ++row)
      if (!std::isfinite(y_values[row]))
        throw runtime::Error("immersed-flow pressure result is non-finite");
    return execution::ExecutionEvent::completed();
  }

  bool has_diagonal() const override { return true; }
  execution::ExecutionEvent
  diagonal(execution::VectorView<double> output) const override {
    if (output.size() != diagonal_.size() || output.stride() != 1U ||
        !output.writable() ||
        output.backend_identity() != context_->backend_identity()) {
      throw runtime::Error("immersed-flow pressure diagonal view is incompatible");
    }
    double *const values = output.data();
    for (std::size_t index = 0U; index < diagonal_.size(); ++index)
      values[index] = diagonal_[index];
    return execution::ExecutionEvent::completed();
  }

  double exchanged_value(mesh::GlobalCellId id) const {
    const auto offset = active_offset_if_present(id);
    if (offset == kMissingGlobalOffset)
      throw runtime::Error("immersed-flow exchanged pressure value is unavailable");
    const auto values =
        static_cast<const linear::GhostedVector &>(exchange_values_)
            .local_view();
    return values[offset] / active_sqrt_volumes_[offset];
  }

  void exchange_scaled_values(
      execution::VectorView<const double> scaled_values) const {
    if (scaled_values.size() != layout_.owned_count() ||
        scaled_values.stride() != 1U ||
        scaled_values.backend_identity() != context_->backend_identity())
      throw runtime::Error(
          "immersed-flow scaled pressure exchange view is incompatible");
    auto owned = exchange_values_.owned_view();
    std::copy(scaled_values.data(), scaled_values.data() + scaled_values.size(),
              owned.data());
    active_halo_.exchange(exchange_values_);
  }

  bool has_pressure_reference() const noexcept {
    return has_pressure_reference_;
  }

  std::uint64_t last_apply_schedule() const noexcept {
    return last_apply_schedule_;
  }

private:
  friend class ActivePressureTwoLevelPreconditioner;

  void build_connections(const immersed::ActiveCellLayout &active,
                         const immersed::ActiveBoundaryLayout &active_boundary,
                         const boundary::BoundaryRegistry &boundaries) {
    for (mesh::LocalFaceId face = 0U; face < topology_->local_face_count();
         ++face) {
      const auto neighbour = topology_->neighbour(face);
      const auto add = [&](mesh::LocalCellId cell, mesh::LocalCellId other) {
        if (topology_->cell_ownership(cell) != mesh::EntityOwnership::owned)
          return;
        const auto row = active.active_index(cell);
        const auto other_index = active.active_index(other);
        if (!row.has_value() || !other_index.has_value() ||
            *row >= active.owned_active_count())
          return;
        const auto area =
            geometry_->face_area_vector_m2(face, mesh::FaceSide::owner);
        const auto displacement = geometry_->face_displacement_m(face);
        const double projection = dot(area, displacement);
        const double factor = dot(area, area) / projection;
        const double cell_volume = geometry_->cell_volume_m3(cell);
        const double other_volume = geometry_->cell_volume_m3(other);
        const double diagonal_coefficient = factor / cell_volume;
        const double neighbour_coefficient =
            factor / std::sqrt(cell_volume * other_volume);
        if (!(projection > 0.0) || !(diagonal_coefficient > 0.0) ||
            !(neighbour_coefficient > 0.0) ||
            !std::isfinite(diagonal_coefficient) ||
            !std::isfinite(neighbour_coefficient))
          throw runtime::Error("immersed-flow pressure face metric is invalid");
        connections_.push_back({*row, face, topology_->global_cell_id(other),
                                diagonal_coefficient, neighbour_coefficient,
                                0.0, false, false, kMissingGlobalOffset});
      };
      if (!neighbour.has_value()) {
        const auto patch = topology_->patch_id(face);
        if (!patch.has_value() || boundaries.patch(*patch).pressure_rule() !=
                                      boundary::PressureRule::prescribed_value)
          continue;
        const auto &active_faces = active_boundary.patch_faces(*patch);
        if (std::find(active_faces.begin(), active_faces.end(),
                      topology_->global_face_id(face)) == active_faces.end())
          continue;
        const auto owner = topology_->owner(face);
        const auto row = active.active_index(owner);
        if (!row.has_value() || *row >= active.owned_active_count())
          continue;
        const auto area =
            geometry_->face_area_vector_m2(face, mesh::FaceSide::owner);
        const auto displacement = geometry_->face_displacement_m(face);
        const double projection = dot(area, displacement);
        const double factor = dot(area, area) / projection;
        const double coefficient = factor / geometry_->cell_volume_m3(owner);
        if (!(coefficient > 0.0) || !std::isfinite(coefficient))
          throw runtime::Error("immersed-flow pressure reference metric is invalid");
        connections_.push_back({*row, face, topology_->global_cell_id(owner),
                                coefficient, 0.0, 0.0, true, false,
                                kMissingGlobalOffset});
        has_pressure_reference_ = true;
        continue;
      }
      add(topology_->owner(face), *neighbour);
      if (!topology_->periodic_pair(face).has_value())
        add(*neighbour, topology_->owner(face));
    }
    std::sort(connections_.begin(), connections_.end(),
              [](const auto &left, const auto &right) {
                return std::tuple{left.row, left.neighbour,
                                  left.pressure_reference} <
                       std::tuple{right.row, right.neighbour,
                                  right.pressure_reference};
              });
  }

  void bind_connections() {
    for (auto &connection : connections_) {
      if (connection.pressure_reference)
        continue;
      const auto offset = active_offset_if_present(connection.neighbour);
      if (offset == kMissingGlobalOffset)
        throw runtime::Error(
            "immersed-flow pressure connection neighbour is unavailable");
      connection.local_neighbour = offset < layout_.owned_count();
      connection.value_offset = offset;
      connection.neighbour_rank = connection.local_neighbour
                                      ? mpi_->rank()
                                      : owning_rank(connection.neighbour);
    }
  }

  static int owner_coordinate(int extent, int partitions, int coordinate) {
    if (extent <= 0 || partitions <= 0 || coordinate < 0 ||
        coordinate >= extent)
      throw runtime::Error(
          "immersed-flow pressure owner coordinate is invalid");
    const int quotient = extent / partitions;
    const int remainder = extent % partitions;
    const int expanded = (quotient + 1) * remainder;
    if (coordinate < expanded)
      return coordinate / (quotient + 1);
    if (quotient == 0)
      throw runtime::Error(
          "immersed-flow pressure decomposition has an empty partition");
    return remainder + (coordinate - expanded) / quotient;
  }

  int owning_rank(mesh::GlobalCellId cell) const {
    const auto local = topology_->find_local_cell(cell);
    if (!local.has_value())
      throw runtime::Error(
          "immersed-flow pressure neighbour owner is unavailable");
    const auto global_cell = topology_->global_cell(*local);
    const auto global_extent = decomposition_->global_extent();
    const auto process_grid = decomposition_->process_grid();
    std::array<int, 3> coordinates{
        owner_coordinate(global_extent.x, process_grid.x, global_cell.x),
        owner_coordinate(global_extent.y, process_grid.y, global_cell.y),
        owner_coordinate(global_extent.z, process_grid.z, global_cell.z)};
    int rank = MPI_PROC_NULL;
    runtime::detail::check_mpi(
        MPI_Cart_rank(decomposition_->comm(), coordinates.data(), &rank),
        "MPI_Cart_rank immersed-flow pressure neighbour owner");
    if (rank == MPI_PROC_NULL || rank < 0 || rank >= mpi_->size())
      throw runtime::Error("immersed-flow pressure neighbour owner is invalid");
    return rank;
  }

  std::size_t active_offset_if_present(mesh::GlobalCellId cell) const noexcept {
    const auto found = active_offsets_.find(cell);
    return found == active_offsets_.end() ? kMissingGlobalOffset
                                          : found->second;
  }

  static constexpr std::size_t kMissingGlobalOffset =
      std::numeric_limits<std::size_t>::max();
  const mesh::MeshTopology *topology_;
  const mesh::MeshGeometry *geometry_;
  const runtime::MpiContext *mpi_;
  execution::ExecutionContext *context_;
  const runtime::StructuredDecomposition *decomposition_;
  linear::VectorLayout layout_;
  mutable linear::GhostedVector exchange_values_;
  mutable linear::GhostedVectorHalo active_halo_;
  std::vector<ActiveConnection> connections_;
  std::vector<double> diagonal_;
  std::vector<mesh::GlobalCellId> owned_ids_;
  std::vector<double> active_sqrt_volumes_;
  std::unordered_map<mesh::GlobalCellId, std::size_t> active_offsets_;
  std::array<std::vector<double>, 3> active_momentum_diagonal_;
  std::vector<double> face_velocity_mobility_;
  std::vector<double> face_mass_flux_coefficient_;
  std::uint64_t revision_{1U};
  bool coefficients_prepared_{};
  mutable bool active_{};
  bool has_pressure_reference_{};
  mutable std::uint64_t last_apply_schedule_{};
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  double test_density_{};
#endif
};

// Two-level compact-pressure approximation: rank-local zero-fill incomplete
// Cholesky for the symmetric positive local block plus one normalized constant
// coarse basis per MPI rank. The global compact operator is not required to be
// symmetric, so the coarse matrix uses pivoted LU and the consumer is FGMRES.
class ActivePressureTwoLevelPreconditioner final
    : public linear::Preconditioner {
public:
  explicit ActivePressureTwoLevelPreconditioner(
      ActivePressureOperator &linear_operator)
      : linear_operator_(&linear_operator),
        context_(&linear_operator.context()) {
    build_pattern();
  }

  void update(const linear::LinearOperator &linear_operator,
              std::uint64_t revision) override {
    if (&linear_operator != linear_operator_ || revision == 0U ||
        linear_operator_->revision() != revision ||
        linear_operator_->domain_layout() != linear_operator_->range_layout() ||
        &linear_operator_->context() != context_)
      throw runtime::Error(
          "immersed-flow block IC(0) preconditioner input is invalid");
    prepared_ = false;
    const std::size_t count = linear_operator_->layout_.owned_count();
    if (linear_operator_->diagonal_.size() != count ||
        row_offsets_.size() != count + 1U)
      throw runtime::Error(
          "immersed-flow block IC(0) preconditioner layout changed");

    std::copy(linear_operator_->diagonal_.begin(),
              linear_operator_->diagonal_.end(), pivot_.begin());
    std::fill(lower_matrix_.begin(), lower_matrix_.end(), 0.0);
    std::fill(upper_matrix_.begin(), upper_matrix_.end(), 0.0);
    for (const auto &connection : linear_operator_->connections_) {
      if (connection.pressure_reference || !connection.local_neighbour)
        continue;
      const std::size_t row = connection.row;
      const std::size_t column = connection.value_offset;
      if (row >= count || column >= count)
        throw runtime::Error(
            "immersed-flow block IC(0) connection is out of range");
      const double value =
          -connection.mass_flux_coefficient * connection.neighbour_coefficient;
      if (!std::isfinite(value))
        throw runtime::Error(
            "immersed-flow block IC(0) coefficient is non-finite");
      if (row == column) {
        pivot_[row] += value;
        continue;
      }
      const std::size_t lower_row = std::max(row, column);
      const std::size_t lower_column = std::min(row, column);
      const std::size_t slot = lower_slot(lower_row, lower_column);
      if (row > column)
        lower_matrix_[slot] += value;
      else
        upper_matrix_[slot] += value;
    }

    constexpr double symmetry_factor = 16384.0;
    for (std::size_t slot = 0U; slot < lower_matrix_.size(); ++slot) {
      const double lower = lower_matrix_[slot];
      const double upper = upper_matrix_[slot];
      const double scale = std::max({1.0, std::abs(lower), std::abs(upper)});
      if (!std::isfinite(lower) || !std::isfinite(upper) ||
          std::abs(lower - upper) >
              symmetry_factor * std::numeric_limits<double>::epsilon() * scale)
        throw runtime::Error(
            "immersed-flow compact pressure block is not symmetric");
      lower_factor_[slot] = 0.5 * (lower + upper);
    }

    for (std::size_t row = 0U; row < count; ++row) {
      const std::size_t begin = row_offsets_[row];
      const std::size_t end = row_offsets_[row + 1U];
      for (std::size_t slot = begin; slot < end; ++slot) {
        const std::size_t column = columns_[slot];
        double value = lower_factor_[slot];
        std::size_t left = begin;
        std::size_t right = row_offsets_[column];
        const std::size_t right_end = row_offsets_[column + 1U];
        while (left < slot && right < right_end) {
          const std::size_t left_column = columns_[left];
          const std::size_t right_column = columns_[right];
          if (left_column < right_column) {
            ++left;
          } else if (right_column < left_column) {
            ++right;
          } else {
            value -= lower_factor_[left] * pivot_[left_column] *
                     lower_factor_[right];
            ++left;
            ++right;
          }
        }
        if (!(pivot_[column] > 0.0) || !std::isfinite(pivot_[column]) ||
            !std::isfinite(value / pivot_[column]))
          throw runtime::Error(
              "immersed-flow block IC(0) factor pivot is invalid");
        lower_factor_[slot] = value / pivot_[column];
      }
      double diagonal = pivot_[row];
      for (std::size_t slot = begin; slot < end; ++slot) {
        const std::size_t column = columns_[slot];
        diagonal -= lower_factor_[slot] * lower_factor_[slot] * pivot_[column];
      }
      if (!(diagonal > 0.0) || !std::isfinite(diagonal))
        throw runtime::Error(
            "immersed-flow block IC(0) factor is not positive definite");
      pivot_[row] = diagonal;
    }
    update_coarse_factor();
    revision_ = revision;
    prepared_ = true;
  }

  execution::ExecutionEvent
  apply(execution::VectorView<const double> residual,
        execution::VectorView<double> correction) const override {
    const std::size_t count = linear_operator_->layout_.owned_count();
    if (!prepared_ || linear_operator_->revision() != revision_ ||
        residual.size() != count || correction.size() != count ||
        residual.stride() != 1U || correction.stride() != 1U ||
        !correction.writable() ||
        residual.backend_identity() != context_->backend_identity() ||
        correction.backend_identity() != context_->backend_identity() ||
        (residual.allocation_identity() == correction.allocation_identity() &&
         residual.offset_bytes() != correction.offset_bytes()))
      throw runtime::Error(
          "immersed-flow block IC(0) preconditioner view is incompatible");
    const double *const input = residual.data();
    double *const output = correction.data();
    for (std::size_t row = 0U; row < count; ++row) {
      double value = input[row];
      for (std::size_t slot = row_offsets_[row]; slot < row_offsets_[row + 1U];
           ++slot)
        value -= lower_factor_[slot] * output[columns_[slot]];
      if (!std::isfinite(value))
        throw runtime::Error(
            "immersed-flow block IC(0) forward solve is non-finite");
      output[row] = value;
    }
    for (std::size_t row = 0U; row < count; ++row) {
      output[row] /= pivot_[row];
      if (!std::isfinite(output[row]))
        throw runtime::Error(
            "immersed-flow block IC(0) diagonal solve is non-finite");
    }
    for (std::size_t reverse = count; reverse-- > 0U;) {
      const double value = output[reverse];
      for (std::size_t slot = row_offsets_[reverse];
           slot < row_offsets_[reverse + 1U]; ++slot)
        output[columns_[slot]] -= lower_factor_[slot] * value;
    }
    if (std::any_of(output, output + count,
                    [](double value) { return !std::isfinite(value); }))
      throw runtime::Error(
          "immersed-flow block IC(0) backward solve is non-finite");
    apply_coarse_correction(residual, output);
    return execution::ExecutionEvent::completed();
  }

private:
  void build_pattern() {
    const std::size_t count = linear_operator_->layout_.owned_count();
    std::vector<std::vector<std::size_t>> rows(count);
    for (const auto &connection : linear_operator_->connections_) {
      if (connection.pressure_reference || !connection.local_neighbour)
        continue;
      const std::size_t row = connection.row;
      const std::size_t column = connection.value_offset;
      if (row >= count || column >= count)
        throw runtime::Error(
            "immersed-flow block IC(0) pattern is out of range");
      if (row != column)
        rows[std::max(row, column)].push_back(std::min(row, column));
    }
    row_offsets_.resize(count + 1U, 0U);
    for (std::size_t row = 0U; row < count; ++row) {
      auto &columns = rows[row];
      std::sort(columns.begin(), columns.end());
      columns.erase(std::unique(columns.begin(), columns.end()), columns.end());
      if (columns.size() >
          std::numeric_limits<std::size_t>::max() - row_offsets_[row])
        throw runtime::Error(
            "immersed-flow block IC(0) pattern size overflows");
      row_offsets_[row + 1U] = row_offsets_[row] + columns.size();
      columns_.insert(columns_.end(), columns.begin(), columns.end());
    }
    lower_matrix_.resize(columns_.size());
    upper_matrix_.resize(columns_.size());
    lower_factor_.resize(columns_.size());
    pivot_.resize(count);
    const std::size_t ranks =
        static_cast<std::size_t>(linear_operator_->mpi_->size());
    if (ranks == 0U || ranks > std::numeric_limits<std::size_t>::max() / ranks)
      throw runtime::Error(
          "immersed-flow pressure coarse rank count is invalid");
    coarse_factor_.resize(ranks * ranks);
    coarse_pivots_.resize(ranks);
    coarse_inverse_sqrt_counts_.resize(ranks);
    coarse_rhs_.resize(ranks);
  }

  void update_coarse_factor() {
    const std::size_t ranks =
        static_cast<std::size_t>(linear_operator_->mpi_->size());
    const std::size_t rank =
        static_cast<std::size_t>(linear_operator_->mpi_->rank());
    const std::uint64_t local_count =
        static_cast<std::uint64_t>(linear_operator_->layout_.owned_count());
    std::vector<std::uint64_t> counts(ranks);
    runtime::detail::check_mpi(
        MPI_Allgather(&local_count, 1, MPI_UINT64_T, counts.data(), 1,
                      MPI_UINT64_T, linear_operator_->mpi_->comm()),
        "MPI_Allgather immersed-flow pressure coarse counts");
    for (std::size_t peer = 0U; peer < ranks; ++peer) {
      if (counts[peer] == 0U)
        throw runtime::Error(
            "immersed-flow pressure coarse rank has no active cells");
      coarse_inverse_sqrt_counts_[peer] =
          1.0 / std::sqrt(static_cast<double>(counts[peer]));
    }

    std::fill(coarse_factor_.begin(), coarse_factor_.end(), 0.0);
    const double row_scale = coarse_inverse_sqrt_counts_[rank];
    for (const double diagonal : linear_operator_->diagonal_)
      coarse_factor_[rank * ranks + rank] += diagonal * row_scale * row_scale;
    for (const auto &connection : linear_operator_->connections_) {
      if (connection.pressure_reference)
        continue;
      if (connection.neighbour_rank < 0 ||
          connection.neighbour_rank >= linear_operator_->mpi_->size())
        throw runtime::Error(
            "immersed-flow pressure coarse neighbour rank is invalid");
      const std::size_t column_rank =
          static_cast<std::size_t>(connection.neighbour_rank);
      const double value =
          -connection.mass_flux_coefficient * connection.neighbour_coefficient;
      coarse_factor_[rank * ranks + column_rank] +=
          value * row_scale * coarse_inverse_sqrt_counts_[column_rank];
    }
    linear_operator_->mpi_->allreduce_fp64_in_place(
        coarse_factor_.data(), coarse_factor_.size(),
        runtime::Fp64ReductionOperation::sum);

    coarse_pinned_ = !linear_operator_->has_pressure_reference_;
    if (coarse_pinned_) {
      for (std::size_t peer = 0U; peer < ranks; ++peer) {
        coarse_factor_[peer] = 0.0;
        coarse_factor_[peer * ranks] = 0.0;
      }
      coarse_factor_[0] = 1.0;
    }
    for (std::size_t column = 0U; column < ranks; ++column) {
      std::size_t selected = column;
      double magnitude = std::abs(coarse_factor_[column * ranks + column]);
      for (std::size_t row = column + 1U; row < ranks; ++row) {
        const double candidate = std::abs(coarse_factor_[row * ranks + column]);
        if (candidate > magnitude) {
          magnitude = candidate;
          selected = row;
        }
      }
      if (!(magnitude > 0.0) || !std::isfinite(magnitude))
        throw runtime::Error(
            "immersed-flow pressure coarse matrix is singular");
      coarse_pivots_[column] = selected;
      if (selected != column)
        for (std::size_t entry = 0U; entry < ranks; ++entry)
          std::swap(coarse_factor_[column * ranks + entry],
                    coarse_factor_[selected * ranks + entry]);
      const double pivot = coarse_factor_[column * ranks + column];
      for (std::size_t row = column + 1U; row < ranks; ++row) {
        coarse_factor_[row * ranks + column] /= pivot;
        const double multiplier = coarse_factor_[row * ranks + column];
        for (std::size_t entry = column + 1U; entry < ranks; ++entry)
          coarse_factor_[row * ranks + entry] -=
              multiplier * coarse_factor_[column * ranks + entry];
      }
    }
    if (std::any_of(coarse_factor_.begin(), coarse_factor_.end(),
                    [](double value) { return !std::isfinite(value); }))
      throw runtime::Error(
          "immersed-flow pressure coarse factor is non-finite");
  }

  void apply_coarse_correction(execution::VectorView<const double> residual,
                               double *output) const {
    const std::size_t ranks = coarse_rhs_.size();
    const std::size_t rank =
        static_cast<std::size_t>(linear_operator_->mpi_->rank());
    std::fill(coarse_rhs_.begin(), coarse_rhs_.end(), 0.0);
    double local = 0.0;
    for (std::size_t row = 0U; row < residual.size(); ++row)
      local += residual[row];
    coarse_rhs_[rank] = local * coarse_inverse_sqrt_counts_[rank];
    linear_operator_->mpi_->allreduce_fp64_in_place(
        coarse_rhs_.data(), coarse_rhs_.size(),
        runtime::Fp64ReductionOperation::sum);
    if (coarse_pinned_)
      coarse_rhs_[0] = 0.0;
    for (std::size_t column = 0U; column < ranks; ++column)
      if (coarse_pivots_[column] != column)
        std::swap(coarse_rhs_[column], coarse_rhs_[coarse_pivots_[column]]);
    for (std::size_t row = 0U; row < ranks; ++row)
      for (std::size_t column = 0U; column < row; ++column)
        coarse_rhs_[row] -=
            coarse_factor_[row * ranks + column] * coarse_rhs_[column];
    for (std::size_t reverse = ranks; reverse-- > 0U;) {
      for (std::size_t column = reverse + 1U; column < ranks; ++column)
        coarse_rhs_[reverse] -=
            coarse_factor_[reverse * ranks + column] * coarse_rhs_[column];
      coarse_rhs_[reverse] /= coarse_factor_[reverse * ranks + reverse];
    }
    const double coarse_value =
        coarse_rhs_[rank] * coarse_inverse_sqrt_counts_[rank];
    if (!std::isfinite(coarse_value))
      throw runtime::Error(
          "immersed-flow pressure coarse correction is non-finite");
    for (std::size_t row = 0U; row < residual.size(); ++row)
      output[row] += coarse_value;
  }

  std::size_t lower_slot(std::size_t row, std::size_t column) const {
    const auto begin =
        columns_.begin() + static_cast<std::ptrdiff_t>(row_offsets_[row]);
    const auto end =
        columns_.begin() + static_cast<std::ptrdiff_t>(row_offsets_[row + 1U]);
    const auto found = std::lower_bound(begin, end, column);
    if (found == end || *found != column)
      throw runtime::Error(
          "immersed-flow block IC(0) pattern entry is unavailable");
    return static_cast<std::size_t>(std::distance(columns_.begin(), found));
  }

  ActivePressureOperator *linear_operator_{};
  const execution::ExecutionContext *context_{};
  std::vector<std::size_t> row_offsets_;
  std::vector<std::size_t> columns_;
  std::vector<double> lower_matrix_;
  std::vector<double> upper_matrix_;
  std::vector<double> lower_factor_;
  std::vector<double> pivot_;
  std::vector<double> coarse_factor_;
  std::vector<std::size_t> coarse_pivots_;
  std::vector<double> coarse_inverse_sqrt_counts_;
  mutable std::vector<double> coarse_rhs_;
  bool coarse_pinned_{};
  std::uint64_t revision_{};
  bool prepared_{};
};

class ExactPredictorSchurOperator final : public linear::LinearOperator {
public:
  using Evaluator = std::function<const ExactPredictorResponse &(
      const std::vector<double> &)>;

  ExactPredictorSchurOperator(
      const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
      const immersed::ActiveCellLayout &active,
      ActivePressureOperator &compact_operator,
      execution::ExecutionContext &context)
      : context_(&context), compact_operator_(&compact_operator),
        layout_(active.owned_active_count(), active.ordered_global_ids()),
        sqrt_volumes_(active.owned_active_count()),
        pressure_workspace_(active.owned_active_count()),
        cached_scaled_input_(active.owned_active_count()) {
    if (compact_operator_->domain_layout() != layout_ ||
        compact_operator_->range_layout() != layout_ ||
        &compact_operator_->context() != context_)
      throw runtime::Error(
          "immersed-flow exact predictor Schur layout is invalid");
    for (std::size_t row = 0U; row < active.owned_active_count(); ++row) {
      const auto local =
          topology.find_local_cell(active.ordered_global_ids()[row]);
      if (!local.has_value())
        throw runtime::Error(
            "immersed-flow exact predictor Schur volume is unavailable");
      sqrt_volumes_[row] = std::sqrt(geometry.cell_volume_m3(*local));
    }
    resize_exact_predictor_response(cached_response_,
                                    active.owned_active_count(),
                                    topology.local_face_count());
  }

  linear::VectorLayout domain_layout() const override { return layout_; }
  linear::VectorLayout range_layout() const override { return layout_; }
  const execution::ExecutionContext &context() const override {
    return *context_;
  }
  std::uint64_t revision() const override { return revision_; }

  void replace(std::uint64_t revision, Evaluator evaluator) {
    if (revision == 0U || compact_operator_->revision() != revision ||
        !evaluator)
      throw runtime::Error(
          "immersed-flow exact predictor Schur input is invalid");
    cache_valid_ = false;
    accumulated_work_ = {};
    revision_ = revision;
    evaluator_ = std::move(evaluator);
  }

  execution::ExecutionEvent
  apply(execution::VectorView<const double> x,
        execution::VectorView<double> y) const override {
    if (!evaluator_ || x.size() != sqrt_volumes_.size() ||
        y.size() != sqrt_volumes_.size() || x.stride() != 1U ||
        y.stride() != 1U || !y.writable() ||
        x.backend_identity() != context_->backend_identity() ||
        y.backend_identity() != context_->backend_identity())
      throw runtime::Error(
          "immersed-flow exact predictor Schur view is incompatible");
    cache_valid_ = false;
    std::uint64_t input_fingerprint = hash_u64(kFnvOffset, revision_);
    input_fingerprint =
        hash_u64(input_fingerprint, static_cast<std::uint64_t>(x.size()));
    for (std::size_t row = 0U; row < x.size(); ++row) {
      cached_scaled_input_[row] = x[row];
      input_fingerprint = hash_u64(input_fingerprint, fp64_bits(x[row]));
      pressure_workspace_[row] = x[row] / sqrt_volumes_[row];
    }
    last_apply_schedule_ = 1U;
    const auto &response = evaluator_(pressure_workspace_);
    last_apply_schedule_ = last_apply_schedule_ * 10U + 2U;
    if (!exact_predictor_response_layout_matches(
            response, y.size(),
            cached_response_.face_mass_flux_increment.size()))
      throw runtime::Error(
          "immersed-flow exact predictor Schur response is invalid");
    for (std::size_t row = 0U; row < y.size(); ++row) {
      y[row] = response.divergence_per_volume[row] * sqrt_volumes_[row];
      if (!std::isfinite(y[row]))
        throw runtime::Error(
            "immersed-flow exact predictor Schur result is non-finite");
    }
    copy_exact_predictor_response(response, cached_response_);
    add_exact_predictor_work(accumulated_work_, response.work);
    cached_input_fingerprint_ = input_fingerprint;
    cached_revision_ = revision_;
    cache_valid_ = true;
    last_apply_schedule_ = last_apply_schedule_ * 10U + 3U;
    last_apply_schedule_ = last_apply_schedule_ * 10U + 4U;
    return execution::ExecutionEvent::completed();
  }

  bool has_diagonal() const override { return false; }
  execution::ExecutionEvent
  diagonal(execution::VectorView<double>) const override {
    throw runtime::Error(
        "immersed-flow exact predictor Schur has no diagonal capability");
  }

  const ExactPredictorResponse &cached_response(
      execution::VectorView<const double> scaled_input) const {
    if (!cache_matches(scaled_input))
      throw runtime::Error(
          "immersed-flow exact predictor Schur cached response is stale");
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    ++cached_response_consumption_count_;
#endif
    return cached_response_;
  }

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  std::uint64_t cached_response_consumption_count() const noexcept {
    return cached_response_consumption_count_;
  }

  bool cached_input_mutation_rejected() const noexcept {
    if (!cache_valid_ || cached_scaled_input_.empty())
      return false;
    std::uint64_t mutated_fingerprint = hash_u64(kFnvOffset, revision_);
    mutated_fingerprint = hash_u64(
        mutated_fingerprint,
        static_cast<std::uint64_t>(cached_scaled_input_.size()));
    bool exact_bits_differ = false;
    for (std::size_t row = 0U; row < cached_scaled_input_.size(); ++row) {
      std::uint64_t bits = fp64_bits(cached_scaled_input_[row]);
      if (row == 0U) {
        bits ^= 1U;
        exact_bits_differ =
            bits != fp64_bits(cached_scaled_input_[row]);
      }
      mutated_fingerprint = hash_u64(mutated_fingerprint, bits);
    }
    return exact_bits_differ &&
           mutated_fingerprint != cached_input_fingerprint_;
  }
#endif

  bool has_pressure_reference() const noexcept {
    return compact_operator_->has_pressure_reference();
  }

  const ExactPredictorResponse::Work &accumulated_work() const noexcept {
    return accumulated_work_;
  }

  std::uint64_t last_apply_schedule() const noexcept {
    return last_apply_schedule_;
  }

private:
  bool cache_matches(
      execution::VectorView<const double> scaled_input) const noexcept {
    if (!cache_valid_ || cached_revision_ != revision_ ||
        scaled_input.size() != cached_scaled_input_.size() ||
        scaled_input.stride() != 1U ||
        scaled_input.backend_identity() != context_->backend_identity())
      return false;
    std::uint64_t fingerprint = hash_u64(kFnvOffset, revision_);
    fingerprint = hash_u64(
        fingerprint, static_cast<std::uint64_t>(scaled_input.size()));
    for (std::size_t row = 0U; row < scaled_input.size(); ++row) {
      if (fp64_bits(scaled_input[row]) !=
          fp64_bits(cached_scaled_input_[row]))
        return false;
      fingerprint = hash_u64(fingerprint, fp64_bits(scaled_input[row]));
    }
    return fingerprint == cached_input_fingerprint_;
  }

  execution::ExecutionContext *context_{};
  ActivePressureOperator *compact_operator_{};
  linear::VectorLayout layout_;
  std::vector<double> sqrt_volumes_;
  mutable std::vector<double> pressure_workspace_;
  mutable std::vector<double> cached_scaled_input_;
  mutable ExactPredictorResponse cached_response_;
  mutable ExactPredictorResponse::Work accumulated_work_;
  std::uint64_t revision_{1U};
  Evaluator evaluator_;
  mutable std::uint64_t cached_input_fingerprint_{};
  mutable std::uint64_t cached_revision_{};
  mutable bool cache_valid_{};
  mutable std::uint64_t last_apply_schedule_{};
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  mutable std::uint64_t cached_response_consumption_count_{};
#endif
};

class CompactSchurPreconditioner final : public linear::Preconditioner {
public:
  CompactSchurPreconditioner(ExactPredictorSchurOperator &exact_operator,
                             ActivePressureOperator &compact_operator,
                             const linear::LinearSolver &solver,
                             linear::Preconditioner &backend) noexcept
      : exact_operator_(&exact_operator), compact_operator_(&compact_operator),
        solver_(&solver), backend_(&backend) {}

  void update(const linear::LinearOperator &linear_operator,
              std::uint64_t revision) override {
    if (&linear_operator != exact_operator_ ||
        exact_operator_->revision() != revision ||
        compact_operator_->revision() != revision ||
        exact_operator_->domain_layout() !=
            compact_operator_->domain_layout() ||
        exact_operator_->range_layout() !=
            compact_operator_->range_layout() ||
        &exact_operator_->context() != &compact_operator_->context())
      throw runtime::Error(
          "immersed-flow compact Schur preconditioner input is invalid");
    backend_->update(*compact_operator_, revision);
    last_failure_.reset();
    revision_ = revision;
    prepared_ = true;
  }

  execution::ExecutionEvent
  apply(execution::VectorView<const double> residual,
        execution::VectorView<double> correction) const override {
    if (!prepared_ || exact_operator_->revision() != revision_ ||
        compact_operator_->revision() != revision_)
      throw runtime::Error(
          "immersed-flow compact Schur preconditioner revision is stale");
    for (std::size_t row = 0U; row < correction.size(); ++row)
      correction[row] = 0.0;
    // This is an algorithmic inexact-inverse contract for the compact
    // Poisson approximation, not an outer scientific acceptance threshold.
    // Outer FGMRES may consume the varying correction produced by the inner
    // Krylov solve. A rank-local IC(0) block supplies a bounded-memory compact
    // approximation without asserting global symmetry across MPI ranks.
    linear::SolveControl control;
    control.atol = 0.0;
    control.rtol = kInnerRelativeTolerance;
    control.max_iterations = kInnerMaximumIterations;
    control.residual_recompute_interval = 20U;
    const auto report = solver_->solve(*compact_operator_, *backend_, residual,
                                       correction, control);
    const bool usable =
        report.reason == linear::SolveTerminationReason::converged ||
        report.reason ==
            linear::SolveTerminationReason::zero_right_hand_side ||
        (report.reason ==
             linear::SolveTerminationReason::maximum_iterations &&
         std::isfinite(report.initial_residual) &&
         std::isfinite(report.final_residual) &&
         report.final_residual < report.initial_residual);
    if (!usable) {
      last_failure_ = report;
      throw runtime::Error(
          "immersed-flow compact Schur approximate solve failed");
    }
    return execution::ExecutionEvent::completed();
  }

  const std::optional<linear::SolveReport> &last_failure() const noexcept {
    return last_failure_;
  }

private:
  static constexpr double kInnerRelativeTolerance = 1.0e-2;
  static constexpr std::uint64_t kInnerMaximumIterations = 200U;
  ExactPredictorSchurOperator *exact_operator_{};
  ActivePressureOperator *compact_operator_{};
  const linear::LinearSolver *solver_{};
  linear::Preconditioner *backend_{};
  mutable std::optional<linear::SolveReport> last_failure_;
  std::uint64_t revision_{};
  bool prepared_{};
};

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
std::unique_ptr<ActivePressureOperator> make_unit_mobility_pressure_probe(
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    const runtime::StructuredDecomposition &decomposition,
    const immersed::ImmersedDomain &domain,
    const boundary::BoundaryRegistry &boundaries,
    const runtime::MpiContext &mpi, execution::ExecutionContext &execution,
    const std::vector<mesh::LocalCellId> &owned_active_cells) {
  auto probe = std::make_unique<ActivePressureOperator>(
      topology, geometry, decomposition, domain.active_cells(),
      domain.active_boundaries(), boundaries, mpi, execution);
  std::array<std::vector<double>, 3> diagonal;
  for (auto &component_diagonal : diagonal)
    component_diagonal.resize(owned_active_cells.size());
  for (std::size_t row = 0U; row < owned_active_cells.size(); ++row)
    for (auto &component_diagonal : diagonal)
      component_diagonal[row] =
          geometry.cell_volume_m3(owned_active_cells[row]);
  probe->prepare_face_coefficients(1.0, diagonal);
  probe->replace(1U);
  return probe;
}

std::array<double, 3> pressure_operator_algebra_probe(
    const linear::LinearOperator &pressure_operator,
    const std::vector<mesh::GlobalCellId> &active_ids,
    execution::ExecutionContext &execution, const runtime::MpiContext &mpi) {
  const auto layout = pressure_operator.domain_layout();
  const std::size_t count = layout.owned_count();
  if (pressure_operator.range_layout() != layout ||
      layout.global_ids().size() != active_ids.size() ||
      !std::equal(layout.global_ids().begin(), layout.global_ids().end(),
                  active_ids.begin(), active_ids.end()))
    throw runtime::Error(
        "immersed-flow pressure algebra layout is incompatible");
  execution::Buffer x(execution, checked_bytes(count));
  execution::Buffer y(execution, checked_bytes(count));
  execution::Buffer ax(execution, checked_bytes(count));
  execution::Buffer ay(execution, checked_bytes(count));
  auto xv = x.view(0U, count);
  auto yv = y.view(0U, count);
  for (std::size_t row = 0U; row < count; ++row) {
    const double phase = static_cast<double>((active_ids[row] % 104729U) + 1U);
    xv[row] = std::sin(phase * 0.017);
    yv[row] = std::cos(phase * 0.013);
  }
  pressure_operator
      .apply(static_cast<const execution::Buffer &>(x).view(0U, count),
             ax.view(0U, count))
      .wait();
  pressure_operator
      .apply(static_cast<const execution::Buffer &>(y).view(0U, count),
             ay.view(0U, count))
      .wait();
  const auto axv = static_cast<const execution::Buffer &>(ax).view(0U, count);
  const auto ayv = static_cast<const execution::Buffer &>(ay).view(0U, count);
  double totals[6]{};
  for (std::size_t row = 0U; row < count; ++row) {
    totals[0] += xv[row] * ayv[row];
    totals[1] += yv[row] * axv[row];
    totals[2] += xv[row] * axv[row];
    totals[3] += yv[row] * ayv[row];
    totals[4] += xv[row] * xv[row];
    totals[5] += yv[row] * yv[row];
  }
  mpi.allreduce_fp64_in_place(totals, 6U, runtime::Fp64ReductionOperation::sum);
  const double symmetry_scale =
      std::max({1.0, std::abs(totals[0]), std::abs(totals[1])});
  return {std::abs(totals[0] - totals[1]) / symmetry_scale,
          totals[2] / totals[4], totals[3] / totals[5]};
}

template <class PressureOperator>
std::vector<double> evaluate_pressure_operator_values(
    PressureOperator &pressure_operator,
    const mesh::MeshGeometry &geometry, execution::ExecutionContext &execution,
    const std::vector<mesh::LocalCellId> &owned_active_cells,
    const std::vector<double> &owned_active_pressure_pa) {
  const std::size_t count = owned_active_cells.size();
  execution::Buffer x(execution, checked_bytes(count));
  execution::Buffer ax(execution, checked_bytes(count));
  auto xv = x.view(0U, count);
  for (std::size_t row = 0U; row < count; ++row)
    xv[row] = owned_active_pressure_pa[row] *
              std::sqrt(geometry.cell_volume_m3(owned_active_cells[row]));
  pressure_operator
      .apply(static_cast<const execution::Buffer &>(x).view(0U, count),
             ax.view(0U, count))
      .wait();
  const auto axv = static_cast<const execution::Buffer &>(ax).view(0U, count);
  std::vector<double> result(count);
  for (std::size_t row = 0U; row < count; ++row)
    result[row] =
        axv[row] / std::sqrt(geometry.cell_volume_m3(owned_active_cells[row]));
  return result;
}

detail::ImmersedFlowPressureFluxIdentityData make_pressure_flux_identity_report(
    const mesh::MeshTopology &topology, const runtime::MpiContext &mpi,
    const std::vector<mesh::LocalCellId> &owned_active_cells,
    const std::vector<double> &operator_values,
    const std::vector<double> &routed_values) {
  const std::size_t count = owned_active_cells.size();
  if (operator_values.size() != count || routed_values.size() != count)
    throw runtime::Error("immersed-flow pressure identity layout is invalid");
  struct Wire final {
    std::uint64_t global_cell_id{};
    double operator_value{};
    double routed_value{};
  };
  static_assert(std::is_trivially_copyable_v<Wire>);
  std::vector<Wire> local(count);
  for (std::size_t row = 0U; row < count; ++row) {
    local[row] = {topology.global_cell_id(owned_active_cells[row]),
                  operator_values[row], routed_values[row]};
  }
  if (local.size() >
      static_cast<std::size_t>(std::numeric_limits<int>::max()) / sizeof(Wire))
    throw runtime::Error("immersed-flow pressure identity exceeds MPI int");
  const int local_bytes = static_cast<int>(local.size() * sizeof(Wire));
  std::vector<int> counts(static_cast<std::size_t>(mpi.size()));
  runtime::detail::check_mpi(
      MPI_Allgather(&local_bytes, 1, MPI_INT, counts.data(), 1, MPI_INT,
                    mpi.comm()),
      "MPI_Allgather(immersed-flow pressure identity counts)");
  std::vector<int> offsets(counts.size(), 0);
  std::size_t total_bytes = 0U;
  for (std::size_t rank = 0U; rank < counts.size(); ++rank) {
    if (counts[rank] < 0 ||
        total_bytes > static_cast<std::size_t>(std::numeric_limits<int>::max()))
      throw runtime::Error("immersed-flow pressure identity payload is unsupported");
    offsets[rank] = static_cast<int>(total_bytes);
    total_bytes += static_cast<std::size_t>(counts[rank]);
  }
  if (total_bytes % sizeof(Wire) != 0U ||
      total_bytes > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    throw runtime::Error("immersed-flow pressure identity payload is invalid");
  std::vector<Wire> gathered(total_bytes / sizeof(Wire));
  runtime::detail::check_mpi(
      MPI_Allgatherv(local.data(), local_bytes, MPI_BYTE, gathered.data(),
                     counts.data(), offsets.data(), MPI_BYTE, mpi.comm()),
      "MPI_Allgatherv(immersed-flow pressure identity values)");
  std::sort(gathered.begin(), gathered.end(),
            [](const auto &left, const auto &right) {
              return left.global_cell_id < right.global_cell_id;
            });
  if (gathered.empty() ||
      std::adjacent_find(gathered.begin(), gathered.end(),
                         [](const auto &left, const auto &right) {
                           return left.global_cell_id == right.global_cell_id;
                         }) != gathered.end())
    throw runtime::Error("immersed-flow pressure identity global layout is invalid");

  detail::ImmersedFlowPressureFluxIdentityData report;
  report.active_global_cell_ids.reserve(gathered.size());
  report.operator_residual_per_volume.reserve(gathered.size());
  report.routed_flux_divergence_per_volume.reserve(gathered.size());
  double maximum_difference = -1.0;
  double operator_square = 0.0;
  double routed_square = 0.0;
  double difference_square = 0.0;
  for (const auto &cell : gathered) {
    if (!std::isfinite(cell.operator_value) ||
        !std::isfinite(cell.routed_value))
      throw runtime::Error("immersed-flow pressure identity is non-finite");
    const double difference = cell.operator_value - cell.routed_value;
    const double absolute_difference = std::abs(difference);
    report.active_global_cell_ids.push_back(cell.global_cell_id);
    report.operator_residual_per_volume.push_back(cell.operator_value);
    report.routed_flux_divergence_per_volume.push_back(cell.routed_value);
    operator_square += cell.operator_value * cell.operator_value;
    routed_square += cell.routed_value * cell.routed_value;
    difference_square += difference * difference;
    report.operator_linf =
        std::max(report.operator_linf, std::abs(cell.operator_value));
    report.routed_linf =
        std::max(report.routed_linf, std::abs(cell.routed_value));
    report.difference_linf =
        std::max(report.difference_linf, absolute_difference);
    if (absolute_difference > maximum_difference) {
      maximum_difference = absolute_difference;
      report.maximum_difference_global_cell_id = cell.global_cell_id;
      report.maximum_difference_operator_value = cell.operator_value;
      report.maximum_difference_routed_value = cell.routed_value;
    }
  }
  report.operator_l2 = std::sqrt(operator_square);
  report.routed_l2 = std::sqrt(routed_square);
  report.difference_l2 = std::sqrt(difference_square);
  return report;
}
#endif

struct ActiveViscousConnection final {
  std::size_t row{};
  mesh::GlobalCellId neighbour{};
  mesh::LocalFaceId face{};
  double coefficient_per_viscosity{};
  double effective_dynamic_viscosity_pa_s{};
  bool local_neighbour{};
  std::size_t value_offset{std::numeric_limits<std::size_t>::max()};
};

struct ActiveWallViscousCoefficient final {
  mesh::GlobalCellId cell{};
  double coefficient_per_viscosity{};
};

struct ActiveWallViscousContribution final {
  std::size_t row{};
  double coefficient_per_viscosity{};
};

class ActiveMomentumOperator final : public linear::LinearOperator {
public:
  ActiveMomentumOperator(
      const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
      const runtime::StructuredDecomposition &decomposition,
      const immersed::ActiveCellLayout &active,
      std::vector<ActiveWallViscousCoefficient> wall_coefficients,
      execution::ExecutionContext &context)
      : context_(&context),
        layout_(active.owned_active_count(), active.ordered_global_ids()),
        exchange_values_(context, layout_),
        active_halo_(linear::GhostedVectorHalo::create(
            decomposition, topology, context, layout_)),
        active_count_(active.ordered_global_ids().size()),
        face_count_(topology.local_face_count()),
        time_diagonal_(active.owned_active_count(), 1.0),
        spatial_diagonal_(active.owned_active_count(), 0.0),
        candidate_spatial_diagonal_(active.owned_active_count(), 0.0),
        diagonal_(active.owned_active_count(), 1.0),
        uniform_active_viscosity_(active.ordered_global_ids().size(), 0.0),
        uniform_face_viscosity_(topology.local_face_count(), 0.0) {
    active_offsets_.reserve(active.ordered_global_ids().size());
    for (std::size_t index = 0U; index < active.ordered_global_ids().size();
         ++index)
      active_offsets_.emplace(active.ordered_global_ids()[index], index);
    connections_ = build_connections(topology, geometry, active);
    for (const auto &wall : wall_coefficients) {
      const auto found = active_offsets_.find(wall.cell);
      if (found == active_offsets_.end() ||
          found->second >= layout_.owned_count() ||
          !(wall.coefficient_per_viscosity > 0.0) ||
          !std::isfinite(wall.coefficient_per_viscosity))
        throw runtime::Error(
            "immersed-flow wall viscous reference coefficient is invalid");
      wall_coefficients_.push_back(
          {found->second, wall.coefficient_per_viscosity});
    }
    candidate_connection_viscosity_.resize(connections_.size(), 0.0);
    bind_connections();
  }

  void replace(const std::vector<double> &time_diagonal,
               double dynamic_viscosity_pa_s) {
    if (!(dynamic_viscosity_pa_s >= 0.0) ||
        !std::isfinite(dynamic_viscosity_pa_s))
      throw runtime::Error("immersed-flow momentum diagonal is invalid");
    std::fill(uniform_active_viscosity_.begin(),
              uniform_active_viscosity_.end(), dynamic_viscosity_pa_s);
    std::fill(uniform_face_viscosity_.begin(), uniform_face_viscosity_.end(),
              dynamic_viscosity_pa_s);
    replace_impl(time_diagonal, uniform_active_viscosity_,
                 uniform_face_viscosity_, dynamic_viscosity_pa_s > 0.0);
  }

  void replace(
      const std::vector<double> &time_diagonal,
      const std::vector<double> &effective_dynamic_viscosity_by_active_cell,
      const std::vector<double> &effective_dynamic_viscosity_by_face) {
    replace_impl(time_diagonal, effective_dynamic_viscosity_by_active_cell,
                 effective_dynamic_viscosity_by_face, true);
  }

  linear::VectorLayout domain_layout() const override { return layout_; }
  linear::VectorLayout range_layout() const override { return layout_; }
  const execution::ExecutionContext &context() const override {
    return *context_;
  }
  std::uint64_t revision() const override { return revision_; }
  execution::ExecutionEvent
  apply(execution::VectorView<const double> x,
        execution::VectorView<double> y) const override {
    return apply_impl(x, y, false);
  }

  execution::ExecutionEvent
  apply_spatial(execution::VectorView<const double> x,
                execution::VectorView<double> y) const {
    return apply_impl(x, y, true);
  }

  bool has_diagonal() const override { return true; }
  const std::vector<double> &assembled_diagonal() const noexcept {
    return diagonal_;
  }
  runtime::HaloPerformanceCounters performance_counters() const noexcept {
    return active_halo_.performance_counters();
  }
  execution::ExecutionEvent
  diagonal(execution::VectorView<double> output) const override {
    if (output.size() != diagonal_.size() || output.stride() != 1U ||
        !output.writable() ||
        output.backend_identity() != context_->backend_identity())
      throw runtime::Error("immersed-flow momentum diagonal view is incompatible");
    for (std::size_t index = 0U; index < diagonal_.size(); ++index)
      output[index] = diagonal_[index];
    return execution::ExecutionEvent::completed();
  }

private:
  void replace_impl(
      const std::vector<double> &time_diagonal,
      const std::vector<double> &effective_dynamic_viscosity_by_active_cell,
      const std::vector<double> &effective_dynamic_viscosity_by_face,
      bool viscous_exchange_required) {
    if (time_diagonal.size() != diagonal_.size() ||
        effective_dynamic_viscosity_by_active_cell.size() != active_count_ ||
        effective_dynamic_viscosity_by_face.size() != face_count_ ||
        std::any_of(time_diagonal.begin(), time_diagonal.end(),
                    [](double value) {
                      return !(value > 0.0) || !std::isfinite(value);
                    }) ||
        std::any_of(effective_dynamic_viscosity_by_active_cell.begin(),
                    effective_dynamic_viscosity_by_active_cell.end(),
                    [](double value) {
                      return !(value >= 0.0) || !std::isfinite(value);
                    }) ||
        std::any_of(effective_dynamic_viscosity_by_face.begin(),
                    effective_dynamic_viscosity_by_face.end(),
                    [](double value) {
                      return !(value >= 0.0) || !std::isfinite(value);
                    }))
      throw runtime::Error("immersed-flow momentum diagonal is invalid");
    if (revision_ == std::numeric_limits<std::uint64_t>::max())
      throw runtime::Error("immersed-flow momentum revision would wrap");

    std::fill(candidate_spatial_diagonal_.begin(),
              candidate_spatial_diagonal_.end(), 0.0);
    for (std::size_t index = 0U; index < connections_.size(); ++index) {
      const auto &connection = connections_[index];
      const double viscosity =
          effective_dynamic_viscosity_by_face[connection.face];
      const double contribution =
          viscosity * connection.coefficient_per_viscosity;
      if (!(viscosity >= 0.0) || !std::isfinite(viscosity) ||
          !(contribution >= 0.0) || !std::isfinite(contribution) ||
          !std::isfinite(candidate_spatial_diagonal_[connection.row] +
                         contribution))
        throw runtime::Error("immersed-flow momentum diagonal is invalid");
      candidate_connection_viscosity_[index] = viscosity;
      candidate_spatial_diagonal_[connection.row] += contribution;
    }
    for (const auto &wall : wall_coefficients_) {
      const double contribution =
          effective_dynamic_viscosity_by_active_cell[wall.row] *
          wall.coefficient_per_viscosity;
      if (!(contribution >= 0.0) || !std::isfinite(contribution) ||
          !std::isfinite(candidate_spatial_diagonal_[wall.row] +
                         contribution))
        throw runtime::Error("immersed-flow momentum diagonal is invalid");
      candidate_spatial_diagonal_[wall.row] += contribution;
    }
    for (std::size_t row = 0U; row < diagonal_.size(); ++row) {
      const double assembled =
          time_diagonal[row] + candidate_spatial_diagonal_[row];
      if (!(assembled > 0.0) || !std::isfinite(assembled))
        throw runtime::Error("immersed-flow momentum diagonal is invalid");
    }

    time_diagonal_ = time_diagonal;
    spatial_diagonal_.swap(candidate_spatial_diagonal_);
    for (std::size_t index = 0U; index < connections_.size(); ++index)
      connections_[index].effective_dynamic_viscosity_pa_s =
          candidate_connection_viscosity_[index];
    for (std::size_t row = 0U; row < diagonal_.size(); ++row)
      diagonal_[row] = time_diagonal_[row] + spatial_diagonal_[row];
    viscous_exchange_required_ = viscous_exchange_required;
    ++revision_;
  }
  execution::ExecutionEvent apply_impl(execution::VectorView<const double> x,
                                       execution::VectorView<double> y,
                                       bool spatial_only) const {
    validate_views(x, y);
    if (active_)
      throw runtime::Error("immersed-flow momentum operator is already active");
    struct Guard final {
      explicit Guard(bool &active) : active_(&active) { active = true; }
      ~Guard() { *active_ = false; }
      bool *active_;
    } guard(active_);
    for (std::size_t row = 0U; row < diagonal_.size(); ++row) {
      const double coefficient = spatial_diagonal_[row] +
                                 (spatial_only ? 0.0 : time_diagonal_[row]);
      y[row] = coefficient * x[row];
    }
    if (!viscous_exchange_required_)
      return execution::ExecutionEvent::completed();
    auto exchange_owned = exchange_values_.owned_view();
    std::copy(x.data(), x.data() + x.size(), exchange_owned.data());
    active_halo_.begin(exchange_values_);
    const auto subtract_neighbour = [&](const ActiveViscousConnection &item,
                                        double value) {
      y[item.row] -= item.effective_dynamic_viscosity_pa_s *
                     item.coefficient_per_viscosity * value;
    };
    for (const auto &connection : connections_)
      if (connection.local_neighbour)
        subtract_neighbour(connection, x[connection.value_offset]);
    active_halo_.wait(exchange_values_);
    const auto local =
        static_cast<const linear::GhostedVector &>(exchange_values_)
            .local_view();
    for (const auto &connection : connections_) {
      if (connection.local_neighbour)
        continue;
      if (connection.value_offset == kMissingGlobalOffset)
        throw runtime::Error(
            "immersed-flow momentum reference neighbour is unavailable");
      subtract_neighbour(connection, local[connection.value_offset]);
    }
    for (std::size_t row = 0U; row < y.size(); ++row)
      if (!std::isfinite(y[row]))
        throw runtime::Error("immersed-flow momentum result is non-finite");
    return execution::ExecutionEvent::completed();
  }

  void validate_views(execution::VectorView<const double> x,
                      execution::VectorView<double> y) const {
    if (x.size() != diagonal_.size() || y.size() != diagonal_.size() ||
        x.stride() != 1U || y.stride() != 1U || !y.writable() ||
        x.backend_identity() != context_->backend_identity() ||
        y.backend_identity() != context_->backend_identity())
      throw runtime::Error("immersed-flow momentum view is incompatible");
  }

  std::vector<ActiveViscousConnection>
  build_connections(const mesh::MeshTopology &topology,
                    const mesh::MeshGeometry &geometry,
                    const immersed::ActiveCellLayout &active) const {
    std::vector<ActiveViscousConnection> connections;
    const auto add = [&](mesh::LocalFaceId face, mesh::LocalCellId cell,
                         mesh::LocalCellId other) {
      if (topology.cell_ownership(cell) != mesh::EntityOwnership::owned)
        return;
      const auto row = active.active_index(cell);
      const auto neighbour = active.active_index(other);
      if (!row.has_value() || !neighbour.has_value() ||
          *row >= active.owned_active_count())
        return;
      const auto area =
          geometry.face_area_vector_m2(face, mesh::FaceSide::owner);
      const auto displacement = geometry.face_displacement_m(face);
      const double projection = dot(area, displacement);
      const double coefficient = dot(area, area) / projection;
      if (!(projection > 0.0) || !(coefficient > 0.0) ||
          !std::isfinite(coefficient))
        throw runtime::Error(
            "immersed-flow viscous reference face metric is invalid");
      connections.push_back({*row, topology.global_cell_id(other), face,
                             coefficient, 0.0, false,
                             kMissingGlobalOffset});
    };
    for (mesh::LocalFaceId face = 0U; face < topology.local_face_count();
         ++face) {
      const auto neighbour = topology.neighbour(face);
      if (!neighbour.has_value())
        continue;
      add(face, topology.owner(face), *neighbour);
      if (!topology.periodic_pair(face).has_value())
        add(face, *neighbour, topology.owner(face));
    }
    std::sort(connections.begin(), connections.end(),
              [](const auto &left, const auto &right) {
                return std::tuple{left.row, left.neighbour, left.face,
                                  left.coefficient_per_viscosity} <
                       std::tuple{right.row, right.neighbour, right.face,
                                  right.coefficient_per_viscosity};
              });
    return connections;
  }

  void bind_connections() {
    for (auto &connection : connections_) {
      const auto found = active_offsets_.find(connection.neighbour);
      if (found == active_offsets_.end())
        throw runtime::Error(
            "immersed-flow momentum connection neighbour is unavailable");
      connection.local_neighbour = found->second < layout_.owned_count();
      connection.value_offset = found->second;
    }
  }

  static constexpr std::size_t kMissingGlobalOffset =
      std::numeric_limits<std::size_t>::max();
  execution::ExecutionContext *context_;
  linear::VectorLayout layout_;
  mutable linear::GhostedVector exchange_values_;
  mutable linear::GhostedVectorHalo active_halo_;
  std::size_t active_count_{};
  std::size_t face_count_{};
  std::vector<ActiveViscousConnection> connections_;
  std::vector<ActiveWallViscousContribution> wall_coefficients_;
  std::vector<double> time_diagonal_;
  std::vector<double> spatial_diagonal_;
  std::vector<double> candidate_spatial_diagonal_;
  std::vector<double> candidate_connection_viscosity_;
  std::vector<double> diagonal_;
  std::vector<double> uniform_active_viscosity_;
  std::vector<double> uniform_face_viscosity_;
  std::unordered_map<mesh::GlobalCellId, std::size_t> active_offsets_;
  bool viscous_exchange_required_{};
  mutable bool active_{};
  std::uint64_t revision_{1U};
};

enum class AttemptFailureStage : std::uint8_t {
  none,
  after_corrector_1,
  after_corrector_2,
  final_momentum_residual,
  after_final_transport,
  final_wall_penetration,
  final_continuity_residual,
  final_pressure_residual,
  final_force_reconstruction,
  before_commit
};

enum class WallInputFault : std::uint8_t {
  none,
  non_positive_density,
  non_finite_density,
  non_positive_coefficient,
  non_finite_coefficient,
  stale_preconditioner_revision
};

constexpr runtime::PhaseId kImmersedFlowScratchPhase = 1830U;
constexpr runtime::ActorId kImmersedFlowScratchActor = 1830U;

runtime::FieldDescriptor
scratch_cell(std::string name, std::uint32_t components, int ghost_width) {
  return {std::move(name),
          "1",
          "stage3_flow",
          runtime::FunctionSpace::cell_average,
          runtime::ScalarType::float64,
          components,
          ghost_width,
          false,
          runtime::RestartPolicy::transient,
          runtime::OutputPolicy::never};
}

runtime::FieldDescriptor scratch_face(std::string name,
                                      std::uint32_t components) {
  return {std::move(name),
          "1",
          "stage3_flow",
          runtime::FunctionSpace::face_value,
          runtime::ScalarType::float64,
          components,
          0,
          false,
          runtime::RestartPolicy::transient,
          runtime::OutputPolicy::never};
}

class ImmersedFlowScratch final {
public:
  ImmersedFlowScratch(runtime::FieldLayoutSet layout, int ghost_width) {
    lagged_velocity = registry.declare_field(
        scratch_cell("stage3_lagged_velocity", 3U, ghost_width));
    velocity_gradient = registry.declare_field(
        scratch_cell("stage3_velocity_gradient", 9U, ghost_width));
    pressure_gradient = registry.declare_field(
        scratch_cell("stage3_pressure_gradient", 3U, ghost_width));
    momentum_face =
        registry.declare_field(scratch_face("stage3_momentum_face", 3U));
    dynamic_viscosity =
        registry.declare_field(scratch_face("stage3_dynamic_viscosity", 1U));
    dynamic_viscosity_cell = registry.declare_field(
        scratch_cell("stage3_dynamic_viscosity_cell", 1U, ghost_width));
    momentum_residual =
        registry.declare_field(scratch_cell("stage3_momentum_residual", 3U, 0));
    pressure_probe_mass_flux = finite_volume::declare_face_mass_flux(registry);
    registry.freeze();
    access = std::make_unique<runtime::FieldAccessPlan>(registry);
    for (runtime::FieldId id = 0U;
         id < static_cast<runtime::FieldId>(registry.size()); ++id)
      access->declare_access(kImmersedFlowScratchPhase, kImmersedFlowScratchActor, id,
                             runtime::AccessMode::read_write);
    access->freeze();
    storage = std::make_unique<runtime::FieldStorage>(registry, layout);
  }

  runtime::FieldRegistry registry;
  std::unique_ptr<runtime::FieldAccessPlan> access;
  std::unique_ptr<runtime::FieldStorage> storage;
  runtime::FieldId lagged_velocity{};
  runtime::FieldId velocity_gradient{};
  runtime::FieldId pressure_gradient{};
  runtime::FieldId momentum_face{};
  runtime::FieldId dynamic_viscosity{};
  runtime::FieldId dynamic_viscosity_cell{};
  runtime::FieldId momentum_residual{};
  runtime::FieldId pressure_probe_mass_flux{};
};

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
std::atomic<AttemptFailureStage> injected_failure_stage{
    AttemptFailureStage::none};
std::atomic<int> injected_failure_rank{-1};
std::atomic<WallInputFault> injected_wall_failure{WallInputFault::none};
std::atomic<int> injected_wall_failure_rank{-1};
std::atomic<double> injected_momentum_time_diagonal_scale{1.0};
std::atomic<int> injected_linear_failure_phase{-1};
std::atomic<std::uint32_t> injected_linear_failure_component{
    std::numeric_limits<std::uint32_t>::max()};
#endif

} // namespace

detail::ImmersedWallPressureCondition
detail::make_immersed_wall_pressure_condition(
    const ImmersedWallPressureInput &input) {
  if (!(input.rho_wall_kg_per_m3 > 0.0) ||
      !std::isfinite(input.rho_wall_kg_per_m3)) {
    throw runtime::Error(
        "immersed pressure wall density must be positive and finite");
  }
  if (!(input.effective_transformed_measure_m2 > 0.0) ||
      !std::isfinite(input.effective_transformed_measure_m2)) {
    throw runtime::Error(
        "immersed pressure transformed measure must be positive and finite");
  }
  if (!finite(input.solid_to_fluid_unit_normal) ||
      !finite(input.momentum_velocity_correction_m3_s_per_kg) ||
      !std::isfinite(input.current_normal_gradient_pa_per_m) ||
      !std::isfinite(input.predictor_mass_flux_kg_per_s)) {
    throw runtime::Error("immersed pressure wall input is non-finite");
  }
  const double normal_length_squared =
      norm_squared(input.solid_to_fluid_unit_normal);
  const double normal_tolerance = 64.0 * std::numeric_limits<double>::epsilon();
  if (std::abs(normal_length_squared - 1.0) > normal_tolerance) {
    throw runtime::Error("immersed pressure wall normal is not unit length");
  }
  const auto &normal = input.solid_to_fluid_unit_normal;
  const auto &correction = input.momentum_velocity_correction_m3_s_per_kg;
  const double normal_mobility = normal.x * normal.x * correction.x +
                                 normal.y * normal.y * correction.y +
                                 normal.z * normal.z * correction.z;
  if (!(normal_mobility > 0.0) || !std::isfinite(normal_mobility)) {
    throw runtime::Error(
        "immersed pressure wall velocity correction must be positive and "
        "finite");
  }
  const double coefficient = input.rho_wall_kg_per_m3 *
                             input.effective_transformed_measure_m2 *
                             normal_mobility;
  if (!(coefficient > 0.0) || !std::isfinite(coefficient)) {
    throw runtime::Error(
        "immersed pressure wall correction coefficient must be positive and "
        "finite");
  }
  const double gradient = input.predictor_mass_flux_kg_per_s / coefficient;
  if (!std::isfinite(gradient)) {
    throw runtime::Error(
        "immersed pressure wall correction gradient is non-finite");
  }
  const double corrected_gradient =
      input.current_normal_gradient_pa_per_m + gradient;
  if (!std::isfinite(corrected_gradient))
    throw runtime::Error(
        "immersed pressure corrected wall gradient is non-finite");
  return {input.link,  input.predictor_mass_flux_kg_per_s,
          coefficient, input.current_normal_gradient_pa_per_m,
          gradient,    corrected_gradient,
          0.0};
}

std::uint64_t detail::make_immersed_pressure_revision(
    std::uint64_t previous_revision, std::uint64_t density_fingerprint,
    std::uint64_t momentum_diagonal_fingerprint,
    std::uint64_t geometry_fingerprint, std::uint64_t active_layout_fingerprint,
    std::uint64_t wall_condition_fingerprint) {
  std::uint64_t hash = kFnvOffset;
  hash = hash_u64(hash, previous_revision);
  hash = hash_u64(hash, density_fingerprint);
  hash = hash_u64(hash, momentum_diagonal_fingerprint);
  hash = hash_u64(hash, geometry_fingerprint);
  hash = hash_u64(hash, active_layout_fingerprint);
  hash = hash_u64(hash, wall_condition_fingerprint);
  if (hash == 0U)
    throw runtime::Error("immersed pressure dependency identity is invalid");
  return hash;
}

struct FixedStepImmersedFlow::Impl final {
  Impl(const runtime::StructuredDecomposition &decomposition,
       const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
       const boundary::BoundaryRegistry &boundaries,
       const immersed::ImmersedDomain *domain,
       const immersed::GhostStencilPlan *ghost_plan,
       const immersed::WallQuadraturePlan *wall_quadrature,
       const immersed::LocalFlowPatternTransform *transform,
       const les::WaleModel *wale, ImmersedFlowDensitySetup density_setup,
       const runtime::MpiContext &mpi,
       execution::ExecutionContext &execution_context,
       runtime::HaloExchange &halo, const linear::LinearSolver &momentum_solver,
       std::array<linear::Preconditioner *, 3> momentum_preconditioners,
       const linear::LinearSolver &pressure_solver,
       linear::Preconditioner &pressure_preconditioner)
      : mpi(&mpi), domain(domain), ghost_plan(ghost_plan),
        wall_quadrature(wall_quadrature), transform(transform), wale(wale),
        density_setup(std::move(density_setup)),
        decomposition(&decomposition), topology(&topology), geometry(&geometry),
        boundaries(&boundaries),
        halo(&halo), execution(&execution_context),
        momentum_solver(&momentum_solver), pressure_solver(&pressure_solver),
        momentum_preconditioners(momentum_preconditioners),
        pressure_preconditioner(&pressure_preconditioner),
        geometry_fingerprint(geometry_identity(topology, geometry)) {
    const FlowFieldIds canonical_fields{};
    const bool canonical_field_identity =
        this->density_setup.fields.density == canonical_fields.density &&
        this->density_setup.fields.velocity == canonical_fields.velocity &&
        this->density_setup.fields.mechanical_pressure ==
            canonical_fields.mechanical_pressure &&
        this->density_setup.fields.face_velocity ==
            canonical_fields.face_velocity &&
        this->density_setup.fields.face_mass_flux ==
            canonical_fields.face_mass_flux &&
        this->density_setup.fields.transported_cell_fields.empty();
    const bool density_setup_valid =
        (this->density_setup.model == config::DensityModel::constant &&
         this->density_setup.registry == nullptr && canonical_field_identity &&
         !this->density_setup.material_transport.has_value() &&
         this->density_setup.ideal_gas_closure == nullptr) ||
        (this->density_setup.model == config::DensityModel::material &&
         this->density_setup.registry != nullptr &&
         this->density_setup.material_transport.has_value() &&
         this->density_setup.ideal_gas_closure == nullptr) ||
        (this->density_setup.model == config::DensityModel::ideal_gas &&
         this->density_setup.registry != nullptr &&
         this->density_setup.material_transport.has_value() &&
         this->density_setup.ideal_gas_closure != nullptr);
    const auto density_setup_status = runtime::collective_status(
        mpi, density_setup_valid,
        "immersed-flow density setup is invalid");
    if (!density_setup_status.ok)
      throw runtime::Error(density_setup_status.message +
                           " (lowest failing rank " +
                           std::to_string(density_setup_status.failing_rank) +
                           ")");
    const bool any_immersed = domain != nullptr || ghost_plan != nullptr ||
                              wall_quadrature != nullptr ||
                              transform != nullptr;
    const bool complete_immersed =
        domain != nullptr && ghost_plan != nullptr && transform != nullptr;
    if (any_immersed && !complete_immersed) {
      throw runtime::Error(
          "immersed-flow immersed dependencies must be supplied as one set");
    }
    if (complete_immersed) {
      if (wall_quadrature != nullptr) {
        const auto wall_halo_status = runtime::collective_status(
            mpi,
            wall_quadrature->maximum_halo_reach() > 0U &&
                wall_quadrature->maximum_halo_reach() <=
                    static_cast<std::uint32_t>(halo.ghost_width()),
            "immersed-flow wall quadrature Halo width is insufficient");
        if (!wall_halo_status.ok)
          throw runtime::Error(
              wall_halo_status.message + " (lowest failing rank " +
              std::to_string(wall_halo_status.failing_rank) + ")");
      }
      physical_boundary_fvm.emplace(
          finite_volume::CellCenteredFvmOperators::create(topology, geometry));
      reconstruction.emplace(finite_volume::ImmersedReconstruction::create(
          topology, geometry, boundaries, *domain, *ghost_plan, decomposition,
          mpi, halo));
      immersed_operator.emplace(finite_volume::ImmersedOperatorAdapter::create(
          topology, geometry, *domain, *ghost_plan, *transform,
          *reconstruction));
      scratch = std::make_unique<ImmersedFlowScratch>(
          runtime::FieldLayoutSet{decomposition.local_extent(),
                                  topology.local_face_count()},
          halo.ghost_width());
      initialize_wall_links();
      initialize_pressure_authority_storage();
      initialize_active_algebra();
      if (this->density_setup.model == config::DensityModel::material) {
        density_adapter.emplace<detail::ImmersedMaterialDensityAdapter>(
            detail::ImmersedMaterialDensityAdapter::create(
                *this->density_setup.registry, decomposition, topology,
                geometry, boundaries, *domain, mpi, halo,
                this->density_setup.fields,
                *this->density_setup.material_transport));
      } else if (this->density_setup.model ==
                 config::DensityModel::ideal_gas) {
        density_adapter.emplace<detail::ImmersedIdealGasDensityAdapter>(
            detail::ImmersedIdealGasDensityAdapter::create(
                *this->density_setup.registry, decomposition, topology,
                geometry, boundaries, *domain, mpi, halo,
                this->density_setup.fields,
                *this->density_setup.material_transport,
                *this->density_setup.ideal_gas_closure));
      }
      if (wall_quadrature != nullptr) {
        wall_force_integrator.emplace(
            immersed::WallForceIntegrator::create(*wall_quadrature, mpi));
      }
    } else {
      if (this->density_setup.model != config::DensityModel::constant) {
        if (this->density_setup.model == config::DensityModel::ideal_gas) {
          const bool closure_matches =
              detail::DensityClosureAdapter::matches_body_fitted(
                  *this->density_setup.ideal_gas_closure, topology, geometry,
                  boundaries, mpi, *this->density_setup.registry,
                  this->density_setup.fields);
          const auto closure_status = runtime::collective_status(
              mpi, closure_matches,
              "body-fitted ideal-gas closure collaborator does not match");
          if (!closure_status.ok)
            throw runtime::Error(
                closure_status.message + " (lowest failing rank " +
                std::to_string(closure_status.failing_rank) + ")");
        }
        body_fitted_material.emplace(FixedStepMaterialDensityFlow::create(
            decomposition, topology, geometry, boundaries, mpi,
            execution_context, halo, momentum_solver,
            momentum_preconditioners, pressure_solver,
            pressure_preconditioner, *this->density_setup.registry,
            this->density_setup.fields,
            *this->density_setup.material_transport));
        if (this->density_setup.model == config::DensityModel::ideal_gas) {
          body_fitted_closure_hooks.emplace(
              detail::DensityClosureAdapter::bind(
                  *this->density_setup.ideal_gas_closure,
                  this->density_setup.material_transport->enthalpy_density,
                  0.0));
        }
      } else {
        body_fitted.emplace(FixedStepConstantDensityFlow::create(
            decomposition, topology, geometry, boundaries, mpi,
            execution_context, halo, momentum_solver,
            momentum_preconditioners, pressure_solver,
            pressure_preconditioner));
      }
    }
  }

  struct WallRuntime final {
    immersed::ImmersedLink link;
    mesh::LocalCellId fluid{};
    mesh::LocalFaceId face{};
    double effective_measure_m2{};
  };

  void initialize_pressure_authority_storage() {
    const auto make_storage = [&]() {
      std::vector<finite_volume::detail::ImmersedWallNormalGradient> result(
          wall_links.size());
      for (std::size_t index = 0U; index < wall_links.size(); ++index)
        result[index] = {wall_links[index].link.id, 0.0};
      return result;
    };
    history_pressure_wall_gradients = make_storage();
    committed_pressure_wall_gradients = make_storage();
    pending_pressure_wall_gradients = make_storage();
  }

  std::vector<ActiveWallViscousCoefficient>
  active_viscous_wall_coefficients() const {
    const auto &active = domain->active_cells();
    std::vector<ActiveWallViscousCoefficient> result;
    result.reserve(wall_links.size() + 6U);
    for (const auto &wall : wall_links) {
      const double distance = std::sqrt(norm_squared(subtract(
          geometry->cell_center_m(wall.fluid), wall.link.wall_intercept_m)));
      const double coefficient = wall.effective_measure_m2 / distance;
      if (!(distance > 0.0) || !(coefficient > 0.0) ||
          !std::isfinite(coefficient))
        throw runtime::Error(
            "immersed-flow wall viscous reference metric is invalid");
      result.push_back({topology->global_cell_id(wall.fluid), coefficient});
    }
    for (mesh::LocalFaceId face = 0U; face < topology->local_face_count();
         ++face) {
      if (topology->neighbour(face).has_value())
        continue;
      const auto patch = topology->patch_id(face);
      if (!patch.has_value())
        continue;
      const auto rule = boundaries->patch(*patch).velocity_rule();
      if (rule != boundary::VelocityRule::prescribed_zero &&
          rule != boundary::VelocityRule::prescribed_inlet)
        continue;
      const auto owner = topology->owner(face);
      const auto row = active.active_index(owner);
      if (!row.has_value() || *row >= active.owned_active_count())
        continue;
      const auto area =
          geometry->face_area_vector_m2(face, mesh::FaceSide::owner);
      const double area_magnitude = std::sqrt(norm_squared(area));
      const double distance = std::sqrt(norm_squared(subtract(
          geometry->face_center_m(face), geometry->cell_center_m(owner))));
      const double coefficient = area_magnitude / distance;
      if (!(area_magnitude > 0.0) || !(distance > 0.0) ||
          !(coefficient > 0.0) || !std::isfinite(coefficient))
        throw runtime::Error(
            "immersed-flow physical-wall viscous reference metric is invalid");
      result.push_back({topology->global_cell_id(owner), coefficient});
    }
    return result;
  }

  void initialize_active_algebra() {
    const auto &active = domain->active_cells();
    const auto wall_coefficients = active_viscous_wall_coefficients();
    pressure_operator = std::make_unique<ActivePressureOperator>(
        *topology, *geometry, *decomposition, active,
        domain->active_boundaries(), *boundaries, *mpi, *execution);
    exact_pressure_operator =
        std::make_unique<ExactPredictorSchurOperator>(
            *topology, *geometry, active, *pressure_operator, *execution);
    compact_pressure_preconditioner =
        std::make_unique<ActivePressureTwoLevelPreconditioner>(
            *pressure_operator);
    // The compact solve is nested inside the outer pressure FGMRES. A
    // separate short-restart instance bounds simultaneous Krylov storage;
    // it does not change either solve's residual or iteration contract.
    compact_pressure_solver = std::make_unique<linear::RestartedFgmresSolver>(
        *execution, *mpi, kCompactPressureRestartLength);
    compact_schur_preconditioner = std::make_unique<CompactSchurPreconditioner>(
        *exact_pressure_operator, *pressure_operator, *compact_pressure_solver,
        *compact_pressure_preconditioner);
    for (std::size_t component_index = 0U; component_index < 3U;
         ++component_index) {
      momentum_operators[component_index] =
          std::make_unique<ActiveMomentumOperator>(
              *topology, *geometry, *decomposition, active, wall_coefficients,
              *execution);
      momentum_rhs[component_index] = std::make_unique<execution::Buffer>(
          *execution, checked_bytes(active.owned_active_count()));
      momentum_solution[component_index] = std::make_unique<execution::Buffer>(
          *execution, checked_bytes(active.owned_active_count()));
    }
    pressure_rhs = std::make_unique<execution::Buffer>(
        *execution, checked_bytes(active.owned_active_count()));
    pressure_solution = std::make_unique<execution::Buffer>(
        *execution, checked_bytes(active.owned_active_count()));
    pressure_residual = std::make_unique<execution::Buffer>(
        *execution, checked_bytes(active.owned_active_count()));
    owned_active_cells.reserve(active.owned_active_count());
    for (std::size_t index = 0U; index < active.owned_active_count(); ++index) {
      const auto local =
          topology->find_local_cell(active.ordered_global_ids()[index]);
      if (!local.has_value() ||
          topology->cell_ownership(*local) != mesh::EntityOwnership::owned)
        throw runtime::Error(
            "immersed-flow owned active layout does not map to owned cells");
      owned_active_cells.push_back(*local);
    }
    attempt_continuity_divergence.resize(active.owned_active_count());
  }


  void initialize_wall_links() {
    struct Wire final {
      std::uint64_t id{};
      std::uint64_t fluid{};
      std::uint64_t solid{};
      std::uint64_t triangle{};
      double wall[3]{};
      double normal[3]{};
      double fraction{};
    };
    std::vector<Wire> local;
    local.reserve(domain->links().size());
    for (const auto &link : domain->links()) {
      local.push_back(
          {link.id,
           link.fluid_cell,
           link.solid_cell,
           link.triangle,
           {link.wall_intercept_m.x, link.wall_intercept_m.y,
            link.wall_intercept_m.z},
           {link.solid_to_fluid_normal.x, link.solid_to_fluid_normal.y,
            link.solid_to_fluid_normal.z},
           link.fluid_to_wall_fraction});
    }
    if (local.size() >
        static_cast<std::size_t>(std::numeric_limits<int>::max()) /
            sizeof(Wire))
      throw runtime::Error("immersed-flow wall-link payload exceeds MPI int");
    const int local_bytes = static_cast<int>(local.size() * sizeof(Wire));
    std::vector<int> counts(static_cast<std::size_t>(mpi->size()));
    runtime::detail::check_mpi(MPI_Allgather(&local_bytes, 1, MPI_INT,
                                             counts.data(), 1, MPI_INT,
                                             mpi->comm()),
                               "MPI_Allgather(immersed-flow wall-link counts)");
    std::vector<int> offsets(counts.size(), 0);
    std::size_t total_bytes = 0U;
    for (std::size_t rank = 0U; rank < counts.size(); ++rank) {
      if (counts[rank] < 0 ||
          total_bytes >
              static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw runtime::Error("immersed-flow wall-link payload is unsupported");
      offsets[rank] = static_cast<int>(total_bytes);
      total_bytes += static_cast<std::size_t>(counts[rank]);
    }
    if (total_bytes % sizeof(Wire) != 0U ||
        total_bytes > static_cast<std::size_t>(std::numeric_limits<int>::max()))
      throw runtime::Error("immersed-flow wall-link payload size is invalid");
    std::vector<Wire> global(total_bytes / sizeof(Wire));
    runtime::detail::check_mpi(
        MPI_Allgatherv(local.data(), local_bytes, MPI_BYTE, global.data(),
                       counts.data(), offsets.data(), MPI_BYTE, mpi->comm()),
        "MPI_Allgatherv(immersed-flow wall links)");
    std::sort(
        global.begin(), global.end(),
        [](const auto &left, const auto &right) { return left.id < right.id; });
    if (std::adjacent_find(global.begin(), global.end(),
                           [](const auto &left, const auto &right) {
                             return left.id == right.id;
                           }) != global.end())
      throw runtime::Error("immersed-flow wall-link IDs are duplicated");
    struct FacePair final {
      mesh::GlobalCellId first{};
      mesh::GlobalCellId second{};
      mesh::LocalFaceId face{};
    };
    std::vector<FacePair> face_pairs;
    face_pairs.reserve(topology->local_face_count());
    for (mesh::LocalFaceId face = 0U; face < topology->local_face_count();
         ++face) {
      const auto neighbour = topology->neighbour(face);
      if (!neighbour.has_value())
        continue;
      const auto owner_id = topology->global_cell_id(topology->owner(face));
      const auto neighbour_id = topology->global_cell_id(*neighbour);
      face_pairs.push_back({std::min(owner_id, neighbour_id),
                            std::max(owner_id, neighbour_id), face});
    }
    std::sort(face_pairs.begin(), face_pairs.end(),
              [](const auto &left, const auto &right) {
                return std::tuple{left.first, left.second, left.face} <
                       std::tuple{right.first, right.second, right.face};
              });
    for (const auto &wire : global) {
      const auto fluid = topology->find_local_cell(wire.fluid);
      if (!fluid.has_value() ||
          topology->cell_ownership(*fluid) != mesh::EntityOwnership::owned)
        continue;
      const auto solid = topology->find_local_cell(wire.solid);
      if (!solid.has_value())
        throw runtime::Error("immersed-flow wall solid cell is unavailable");
      const auto key = std::pair{std::min(wire.fluid, wire.solid),
                                 std::max(wire.fluid, wire.solid)};
      const auto found = std::lower_bound(
          face_pairs.begin(), face_pairs.end(), key,
          [](const FacePair &candidate, const auto &target) {
            return std::pair{candidate.first, candidate.second} < target;
          });
      if (found == face_pairs.end() ||
          std::pair{found->first, found->second} != key)
        throw runtime::Error("immersed-flow wall link face is unavailable");
      const auto link_face = found->face;
      const runtime::Real3 normal{wire.normal[0], wire.normal[1],
                                  wire.normal[2]};
      const auto fluid_center = geometry->cell_center_m(*fluid);
      const auto solid_center = geometry->cell_center_m(*solid);
      const auto face_center = geometry->face_center_m(link_face);
      const double center_distance =
          std::sqrt(norm_squared(subtract(solid_center, fluid_center)));
      const double owner_to_face =
          std::sqrt(norm_squared(subtract(face_center, fluid_center)));
      if (!(center_distance > 0.0) || !(owner_to_face > 0.0) ||
          !(owner_to_face < center_distance) ||
          !std::isfinite(center_distance) || !std::isfinite(owner_to_face))
        throw runtime::Error("immersed-flow wall interpolation distance is invalid");
      auto area_from_fluid =
          geometry->face_area_vector_m2(link_face, mesh::FaceSide::owner);
      if (topology->owner(link_face) != *fluid)
        area_from_fluid = multiply(-1.0, area_from_fluid);
      const double effective_measure = -dot(area_from_fluid, normal);
      if (!(effective_measure > 0.0) || !std::isfinite(effective_measure))
        throw runtime::Error("immersed-flow signed wall measure is invalid");
      wall_links.push_back({{wire.id,
                             wire.fluid,
                             wire.solid,
                             wire.triangle,
                             {wire.wall[0], wire.wall[1], wire.wall[2]},
                             normal,
                             wire.fraction},
                            *fluid,
                            link_face,
                            effective_measure});
    }
  }

  struct CellIndex final {
    int i{};
    int j{};
    int k{};
  };

  CellIndex field_index(mesh::LocalCellId cell) const {
    const auto global = topology->global_cell(cell);
    const auto box = topology->owned_global_box();
    const auto extent = topology->global_extent();
    const auto local =
        runtime::Int3{box.end.x - box.begin.x, box.end.y - box.begin.y,
                      box.end.z - box.begin.z};
    const auto axis = [](int coordinate, int begin, int end, int global_size,
                         int local_size) {
      int result = coordinate - begin;
      if (result >= local_size && coordinate == global_size - 1 && begin == 0)
        result = -1;
      if (result < 0 && coordinate == 0 && end == global_size)
        result = local_size;
      return result;
    };
    return {axis(global.x, box.begin.x, box.end.x, extent.x, local.x),
            axis(global.y, box.begin.y, box.end.y, extent.y, local.y),
            axis(global.z, box.begin.z, box.end.z, extent.z, local.z)};
  }

  template <class Function>
  detail::MaterialDensityStageResult
  prepare_density_authority_collectively(Function &&function) const {
    MaterialTransportFailureReason local_failure{
        MaterialTransportFailureReason::none};
    try {
      function();
    } catch (const runtime::detail::MpiOperationError &) {
      throw;
    } catch (const detail::ImmersedDensityAuthorityFailure &failure) {
      local_failure = failure.reason;
    } catch (...) {
      local_failure = MaterialTransportFailureReason::invalid_input;
    }
    const auto status = runtime::collective_status(
        *mpi, local_failure == MaterialTransportFailureReason::none,
        "immersed material density authority is invalid");
    if (status.ok)
      return {};
    int payload{};
    if (mpi->rank() == status.failing_rank)
      payload = static_cast<int>(local_failure);
    runtime::detail::check_mpi(
        MPI_Bcast(&payload, 1, MPI_INT, status.failing_rank, mpi->comm()),
        "MPI_Bcast(immersed material density failure)");
    if (payload < static_cast<int>(MaterialTransportFailureReason::none) ||
        payload > static_cast<int>(
                      MaterialTransportFailureReason::collective_operation))
      throw runtime::Error(
          "immersed material density failure reason is invalid");
    const auto selected =
        static_cast<MaterialTransportFailureReason>(payload);
    if (selected == MaterialTransportFailureReason::invalid_input)
      throw runtime::Error(
          "immersed material density authority construction failed");
    if (selected != MaterialTransportFailureReason::non_finite_state &&
        selected != MaterialTransportFailureReason::non_positive_density)
      throw runtime::Error(
          "immersed material density failure is not recoverable");
    return {selected, status.failing_rank};
  }

  detail::ImmersedDensityAttemptAuthority density_authority(
      FlowState &state, runtime::FieldStorage &layer) const {
    if (domain == nullptr || ghost_plan == nullptr)
      throw runtime::Error("immersed density authority is unavailable");
    const auto density = layer.acquire_read<double>(
        state.solver_access_plan(), kStatePhase, kStateActor,
        state.fields().density);
    detail::ImmersedDensityAttemptAuthority result;
    const auto &active = domain->active_cells();
    result.owned_active_density.resize(active.owned_active_count());
    result.face_density.assign(topology->local_face_count(), 0.0);
    result.wall_density.reserve(wall_links.size());
    result.fingerprint = kFnvOffset;
    for (std::size_t row = 0U; row < active.owned_active_count(); ++row) {
      const auto local =
          topology->find_local_cell(active.ordered_global_ids()[row]);
      if (!local.has_value())
        throw runtime::Error("immersed active density cell is unavailable");
      const auto index = field_index(*local);
      const double value = density(index.i, index.j, index.k, 0);
      if (!std::isfinite(value))
        throw detail::ImmersedDensityAuthorityFailure{
            MaterialTransportFailureReason::non_finite_state};
      if (!(value > 0.0))
        throw detail::ImmersedDensityAuthorityFailure{
            MaterialTransportFailureReason::non_positive_density};
      result.owned_active_density[row] = value;
      result.fingerprint = hash_u64(result.fingerprint, fp64_bits(value));
    }
    for (const auto &wall : wall_links) {
      auto value = detail::reconstruct_immersed_wall_density(
          wall.link, ghost_plan->reconstruction(wall.link.id), density);
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
      if (mpi->rank() ==
          injected_wall_failure_rank.load(std::memory_order_relaxed)) {
        const auto failure =
            injected_wall_failure.load(std::memory_order_relaxed);
        if (failure == WallInputFault::non_positive_density)
          value.rho_wall_kg_per_m3 = 0.0;
        else if (failure == WallInputFault::non_finite_density)
          value.rho_wall_kg_per_m3 =
              std::numeric_limits<double>::quiet_NaN();
      }
#endif
      result.wall_density.push_back(value);
      result.fingerprint =
          hash_u64(result.fingerprint, fp64_bits(value.rho_wall_kg_per_m3));
      result.fingerprint = hash_u64(
          result.fingerprint,
          fp64_bits(value.normal_derivative_kg_per_m4));
    }
    for (mesh::LocalFaceId face = 0U; face < topology->local_face_count();
         ++face) {
      const auto owner = topology->owner(face);
      const auto neighbour = topology->neighbour(face);
      const bool owner_active =
          domain->region(owner) == immersed::CellRegion::fluid;
      const bool neighbour_active =
          neighbour.has_value() &&
          domain->region(*neighbour) == immersed::CellRegion::fluid;
      double value = 0.0;
      if (owner_active && neighbour_active) {
        const auto owner_index = field_index(owner);
        const auto neighbour_index = field_index(*neighbour);
        const auto owner_to_face = subtract(geometry->face_center_m(face),
                                            geometry->cell_center_m(owner));
        const auto face_to_neighbour =
            subtract(geometry->face_displacement_m(face), owner_to_face);
        const double owner_distance = std::sqrt(norm_squared(owner_to_face));
        const double neighbour_distance =
            std::sqrt(norm_squared(face_to_neighbour));
        const double total_distance = owner_distance + neighbour_distance;
        if (!(total_distance > 0.0) || !std::isfinite(total_distance))
          throw runtime::Error("immersed face density metric is invalid");
        value = (neighbour_distance *
                     density(owner_index.i, owner_index.j, owner_index.k, 0) +
                 owner_distance * density(neighbour_index.i,
                                          neighbour_index.j,
                                          neighbour_index.k, 0)) /
                total_distance;
      } else if (!neighbour.has_value() && owner_active) {
        const auto owner_index = field_index(owner);
        value = density(owner_index.i, owner_index.j, owner_index.k, 0);
      } else if (owner_active != neighbour_active) {
        const auto wall = std::find_if(
            wall_links.begin(), wall_links.end(),
            [face](const WallRuntime &candidate) {
              return candidate.face == face;
            });
        if (wall == wall_links.end())
          throw runtime::Error(
              "immersed interface face has no density authority row");
        const auto offset =
            static_cast<std::size_t>(std::distance(wall_links.begin(), wall));
        value = result.wall_density[offset].rho_wall_kg_per_m3;
      }
      if ((owner_active || neighbour_active) && !std::isfinite(value))
        throw detail::ImmersedDensityAuthorityFailure{
            MaterialTransportFailureReason::non_finite_state};
      if ((owner_active || neighbour_active) && !(value > 0.0))
        throw detail::ImmersedDensityAuthorityFailure{
            MaterialTransportFailureReason::non_positive_density};
      result.face_density[face] = value;
      result.fingerprint = hash_u64(result.fingerprint, fp64_bits(value));
    }
    if (result.fingerprint == 0U)
      result.fingerprint = 1U;
    return result;
  }

  detail::ImmersedWaleAttemptAuthority prepare_wale_authority(
      FlowState &state, double molecular_dynamic_viscosity_pa_s,
      const MomentumTimeStencil &stencil,
      runtime::FieldView<const double> rho_attempt) {
    if (wale == nullptr || domain == nullptr || reconstruction == std::nullopt ||
        physical_boundary_fvm == std::nullopt || scratch == nullptr)
      throw runtime::Error("immersed-flow WALE authority is unavailable");
    if (pending_wale_evaluation_count != 0U)
      throw runtime::Error(
          "immersed-flow WALE was evaluated more than once in one attempt");
    ++pending_wale_evaluation_count;
    const auto &active = domain->active_cells();
    const auto &access = state.solver_access_plan();
    const auto fields = state.fields();
    auto &committed = state.solver_layer(FlowLayer::committed);
    auto &history = state.solver_layer(FlowLayer::history);
    const auto velocity_n = committed.acquire_read<double>(
        access, kStatePhase, kStateActor, fields.velocity);
    const auto velocity_nm1 = history.acquire_read<double>(
        access, kStatePhase, kStateActor, fields.velocity);

    auto lagged = scratch->storage->acquire_write<double>(
        *scratch->access, kImmersedFlowScratchPhase, kImmersedFlowScratchActor,
        scratch->lagged_velocity);
    const auto extent = lagged.interior_extent();
    const int ghost = lagged.ghost_width();
    for (int k = -ghost; k < extent.z + ghost; ++k)
      for (int j = -ghost; j < extent.y + ghost; ++j)
        for (int i = -ghost; i < extent.x + ghost; ++i)
          for (int component = 0; component < 3; ++component)
            lagged(i, j, k, component) =
                std::numeric_limits<double>::quiet_NaN();
    const double ratio = stencil.order == MomentumTimeOrder::bdf2
                             ? stencil.dt_s / stencil.previous_dt_s
                             : 0.0;
    for (const auto global_cell : active.ordered_global_ids()) {
      const auto local = topology->find_local_cell(global_cell);
      if (!local.has_value() ||
          domain->region(*local) != immersed::CellRegion::fluid)
        throw runtime::Error(
            "immersed-flow WALE active cell is unavailable");
      const auto index = field_index(*local);
      for (int component = 0; component < 3; ++component) {
        const double current =
            velocity_n(index.i, index.j, index.k, component);
        const double previous =
            velocity_nm1(index.i, index.j, index.k, component);
        const double value = stencil.order == MomentumTimeOrder::bdf2
                                 ? current + ratio * (current - previous)
                                 : current;
        if (!std::isfinite(value))
          throw runtime::Error(
              "immersed-flow WALE lagged velocity is non-finite");
        lagged(index.i, index.j, index.k, component) = value;
      }
    }

    const finite_volume::ReconstructionFieldBinding lagged_binding{
        *scratch->storage, *scratch->access, kImmersedFlowScratchPhase,
        kImmersedFlowScratchActor, scratch->lagged_velocity};
    const finite_volume::ReconstructionFieldBinding gradient_binding{
        *scratch->storage, *scratch->access, kImmersedFlowScratchPhase,
        kImmersedFlowScratchActor, scratch->velocity_gradient};
    reconstruction->compute_gradient(
        finite_volume::GradientScheme::weighted_least_squares,
        finite_volume::FiniteVolumeQuantity::velocity(), lagged_binding,
        gradient_binding);
    halo->exchange(*scratch->storage, scratch->velocity_gradient);
    const auto gradient = scratch->storage->acquire_read<double>(
        *scratch->access, kImmersedFlowScratchPhase, kImmersedFlowScratchActor,
        scratch->velocity_gradient);
    std::uint64_t committed_hash = kFnvOffset;
    std::uint64_t history_hash = kFnvOffset;
    std::uint64_t gradient_hash = kFnvOffset;
    std::uint64_t density_hash = kFnvOffset;
    for (const auto global_cell : active.ordered_global_ids()) {
      const auto local = topology->find_local_cell(global_cell);
      if (!local.has_value())
        throw runtime::Error(
            "immersed-flow WALE fingerprint cell is unavailable");
      const auto index = field_index(*local);
      for (int component = 0; component < 3; ++component) {
        committed_hash = hash_u64(
            committed_hash,
            fp64_bits(velocity_n(index.i, index.j, index.k, component)));
        history_hash = hash_u64(
            history_hash,
            fp64_bits(velocity_nm1(index.i, index.j, index.k, component)));
      }
      for (int component = 0; component < 9; ++component)
        gradient_hash = hash_u64(
            gradient_hash,
            fp64_bits(gradient(index.i, index.j, index.k, component)));
      density_hash = hash_u64(
          density_hash,
          fp64_bits(rho_attempt(index.i, index.j, index.k, 0)));
    }

    auto coefficients = wale->evaluate(
        {state.metadata().step + 1U,
         stencil.dt_s,
         stencil.order == MomentumTimeOrder::bdf2
             ? les::WaleTimeOrder::bdf2
             : les::WaleTimeOrder::backward_euler,
         committed_hash,
         history_hash,
         gradient_hash,
         density_hash,
         gradient,
         rho_attempt});
    if (coefficients.owned_active_count() != active.owned_active_count() ||
        coefficients.local_active_count() !=
            active.ordered_global_ids().size())
      throw runtime::Error(
          "immersed-flow WALE coefficient layout is not active-cell order");

    const auto mu_sgs = coefficients.mu_sgs_pa_s();
    std::vector<double> effective_by_active(mu_sgs.size());
    auto cell_viscosity = scratch->storage->acquire_write<double>(
        *scratch->access, kImmersedFlowScratchPhase, kImmersedFlowScratchActor,
        scratch->dynamic_viscosity_cell);
    const auto cell_extent = cell_viscosity.interior_extent();
    const int cell_ghost = cell_viscosity.ghost_width();
    for (int k = -cell_ghost; k < cell_extent.z + cell_ghost; ++k)
      for (int j = -cell_ghost; j < cell_extent.y + cell_ghost; ++j)
        for (int i = -cell_ghost; i < cell_extent.x + cell_ghost; ++i)
          cell_viscosity(i, j, k, 0) = molecular_dynamic_viscosity_pa_s;
    for (std::size_t row = 0U; row < effective_by_active.size(); ++row) {
      const double effective =
          molecular_dynamic_viscosity_pa_s + mu_sgs[row];
      if (!(effective >= molecular_dynamic_viscosity_pa_s) ||
          !std::isfinite(effective))
        throw runtime::Error(
            "immersed-flow WALE effective viscosity is invalid");
      effective_by_active[row] = effective;
      const auto local =
          topology->find_local_cell(active.ordered_global_ids()[row]);
      if (!local.has_value())
        throw runtime::Error(
            "immersed-flow WALE viscosity cell is unavailable");
      const auto index = field_index(*local);
      cell_viscosity(index.i, index.j, index.k, 0) = effective;
    }

    const auto cell_viscosity_read = scratch->storage->acquire_read<double>(
        *scratch->access, kImmersedFlowScratchPhase, kImmersedFlowScratchActor,
        scratch->dynamic_viscosity_cell);
    auto face_viscosity = scratch->storage->acquire_face_write<double>(
        *scratch->access, kImmersedFlowScratchPhase, kImmersedFlowScratchActor,
        scratch->dynamic_viscosity);
    physical_boundary_fvm->interpolate_cell_scalar_to_faces(
        cell_viscosity_read, face_viscosity);
    for (mesh::LocalFaceId face = 0U; face < topology->local_face_count();
         ++face) {
      const auto owner = active.active_index(topology->owner(face));
      const auto neighbour_cell = topology->neighbour(face);
      const auto neighbour = neighbour_cell.has_value()
                                 ? active.active_index(*neighbour_cell)
                                 : std::optional<std::size_t>{};
      if (neighbour_cell.has_value() &&
          owner.has_value() != neighbour.has_value()) {
        face_viscosity(face, 0) =
            effective_by_active[owner.has_value() ? *owner : *neighbour];
      } else if (!neighbour_cell.has_value() &&
                 !topology->periodic_pair(face).has_value() &&
                 owner.has_value()) {
        face_viscosity(face, 0) = effective_by_active[*owner];
      }
      if (!(face_viscosity(face, 0) >= 0.0) ||
          !std::isfinite(face_viscosity(face, 0)))
        throw runtime::Error(
            "immersed-flow WALE face viscosity is invalid");
    }
    std::vector<double> effective_by_face(topology->local_face_count());
    for (mesh::LocalFaceId face = 0U; face < topology->local_face_count();
         ++face)
      effective_by_face[face] = face_viscosity(face, 0);

    std::uint64_t wall_fingerprint = kFnvOffset;
    for (const auto &wall : wall_links) {
      const auto row = active.active_index(wall.fluid);
      if (!row.has_value())
        throw runtime::Error(
            "immersed-flow WALE wall viscosity row is unavailable");
      wall_fingerprint = hash_u64(wall_fingerprint, wall.link.id);
      wall_fingerprint =
          hash_u64(wall_fingerprint, fp64_bits(effective_by_active[*row]));
    }
    if (wall_fingerprint == 0U)
      wall_fingerprint = 1U;
    return detail::ImmersedWaleAttemptAuthority(
        std::move(coefficients), std::move(effective_by_active),
        std::move(effective_by_face), wall_fingerprint);
  }

  static bool solve_succeeded(linear::SolveTerminationReason reason) noexcept {
    return reason == linear::SolveTerminationReason::converged ||
           reason == linear::SolveTerminationReason::zero_right_hand_side;
  }

  static StepFailureReason material_step_failure(
      MaterialTransportFailureReason reason) noexcept {
    switch (reason) {
    case MaterialTransportFailureReason::none:
      return StepFailureReason::none;
    case MaterialTransportFailureReason::invalid_input:
      return StepFailureReason::invalid_input;
    case MaterialTransportFailureReason::collective_operation:
      return StepFailureReason::collective_operation;
    case MaterialTransportFailureReason::non_finite_state:
    case MaterialTransportFailureReason::non_positive_density:
      return StepFailureReason::transport_failure;
    case MaterialTransportFailureReason::final_density_residual:
      return StepFailureReason::final_continuity_residual;
    case MaterialTransportFailureReason::final_transport_residual:
      return StepFailureReason::final_transport_residual;
    case MaterialTransportFailureReason::final_conservation_defect:
      return StepFailureReason::final_conservation_defect;
    }
    return StepFailureReason::invalid_input;
  }

  static runtime::Real3 add(runtime::Real3 left,
                            runtime::Real3 right) noexcept {
    return {left.x + right.x, left.y + right.y, left.z + right.z};
  }

  static bool finite_force(const immersed::ForceComponents &force) noexcept {
    return finite(force.pressure_N) && finite(force.viscous_N) &&
           finite(force.total_N);
  }

  runtime::CollectiveStatus injected_status(AttemptFailureStage stage) const {
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    const bool armed =
        injected_failure_stage.load(std::memory_order_relaxed) == stage;
    const int configured =
        injected_failure_rank.load(std::memory_order_relaxed);
    const int selected =
        configured >= 0 ? configured : (mpi->size() > 1 ? 1 : 0);
    return runtime::collective_status(*mpi, !armed || mpi->rank() != selected,
                                      "injected immersed-flow transaction failure");
#else
    static_cast<void>(stage);
    return {true, -1, {}};
#endif
  }

  double global_l2(execution::VectorView<const double> values) const {
    double sum = 0.0;
    for (std::size_t index = 0U; index < values.size(); ++index)
      sum += values[index] * values[index];
    mpi->allreduce_fp64_in_place(&sum, 1U,
                                 runtime::Fp64ReductionOperation::sum);
    return std::sqrt(sum);
  }

  linear::SolveControl relative_equation_control(
      const linear::SolveControl &control,
      execution::VectorView<const double> right_hand_side) const {
    linear::SolveControl result = control;
    const double right_hand_side_l2 = global_l2(right_hand_side);
    if (std::isfinite(control.atol) && control.atol >= 0.0 &&
        std::isfinite(control.rtol) && control.rtol > 0.0 &&
        std::isfinite(right_hand_side_l2))
      result.atol =
          std::min(control.atol, control.rtol * right_hand_side_l2);
    return result;
  }

  bool assess_final_momentum(
      FlowState &state, double rho_ref, const MomentumTimeStencil &stencil,
      const std::vector<double> &residual_n,
      const std::vector<double> &residual_nm1,
      const std::vector<double> &implicit_reference_n,
      const std::vector<double> &implicit_reference_nm1,
      const std::vector<double> &implicit_reference_trial,
      const std::vector<double> &pressure_n,
      const std::vector<double> &pressure_nm1,
      const std::vector<finite_volume::detail::ImmersedWallNormalGradient>
          &final_wall_gradients,
      StepAttemptReport &report) const {
    auto &trial = state.solver_layer(FlowLayer::trial);
    auto &committed = state.solver_layer(FlowLayer::committed);
    auto &history = state.solver_layer(FlowLayer::history);
    if (implicit_reference_n.size() != residual_n.size() ||
        implicit_reference_trial.size() != residual_n.size() ||
        (stencil.order == MomentumTimeOrder::bdf2 &&
         implicit_reference_nm1.size() != residual_nm1.size()))
      throw runtime::Error(
          "immersed-flow implicit viscous reference layout is invalid");
    const auto &access = state.solver_access_plan();
    const auto fields = state.fields();
    std::vector<double> pressure_final = assemble_immersed_pressure_residual(
        state, trial, &final_wall_gradients);
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    pending_final_momentum_pressure_residual = pressure_final;
    pending_final_pressure_wall_gradients = final_wall_gradients;
#endif
    const auto velocity = trial.acquire_read<double>(
        access, kStatePhase, kStateActor, fields.velocity);
    const auto velocity_n = committed.acquire_read<double>(
        access, kStatePhase, kStateActor, fields.velocity);
    const auto velocity_nm1 = history.acquire_read<double>(
        access, kStatePhase, kStateActor, fields.velocity);
    const auto density = trial.acquire_read<double>(
        access, kStatePhase, kStateActor, fields.density);
    const auto density_n = committed.acquire_read<double>(
        access, kStatePhase, kStateActor, fields.density);
    const auto density_nm1 = history.acquire_read<double>(
        access, kStatePhase, kStateActor, fields.density);
    const bool variable_density =
        density_setup.model != config::DensityModel::constant;
    std::array<double, 9> norm_sums{};
    std::array<double, 6> conservation_sums{};
    for (std::size_t row = 0U; row < owned_active_cells.size(); ++row) {
      const auto cell = owned_active_cells[row];
      const auto index = field_index(cell);
      const double volume_over_dt =
          geometry->cell_volume_m3(cell) / stencil.dt_s;
      const double rho_trial =
          variable_density
              ? density(index.i, index.j, index.k, 0)
              : rho_ref;
      const double rho_current =
          variable_density
              ? density_n(index.i, index.j, index.k, 0)
              : rho_ref;
      const double rho_history =
          variable_density
              ? density_nm1(index.i, index.j, index.k, 0)
              : rho_ref;
      if (!(rho_trial > 0.0) || !(rho_current > 0.0) ||
          !(rho_history > 0.0) || !std::isfinite(rho_trial) ||
          !std::isfinite(rho_current) || !std::isfinite(rho_history))
        throw runtime::Error(
            "immersed final momentum density authority is invalid");
      for (std::size_t component_index = 0U; component_index < 3U;
           ++component_index) {
        const double trial_term =
            volume_over_dt * stencil.alpha0 * rho_trial *
            velocity(index.i, index.j, index.k,
                     static_cast<int>(component_index));
        const double current_term =
            volume_over_dt * stencil.alpha1 * rho_current *
            velocity_n(index.i, index.j, index.k,
                       static_cast<int>(component_index));
        const double history_term =
            volume_over_dt * stencil.alpha2 * rho_history *
            velocity_nm1(index.i, index.j, index.k,
                         static_cast<int>(component_index));
        const auto offset = row * 3U + component_index;
        const double spatial_term =
            stencil.order == MomentumTimeOrder::backward_euler
                ? residual_n[offset] - pressure_n[offset] -
                      implicit_reference_n[offset] + pressure_final[offset] +
                      implicit_reference_trial[offset]
                : 2.0 * (residual_n[offset] - pressure_n[offset] -
                         implicit_reference_n[offset]) -
                      (residual_nm1[offset] - pressure_nm1[offset] -
                       implicit_reference_nm1[offset]) +
                      pressure_final[offset] + implicit_reference_trial[offset];
        double source_term = 0.0;
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
        if (!manufactured_body_source.empty())
          source_term = -geometry->cell_volume_m3(cell) *
                        manufactured_body_source[offset];
#endif
        const double raw = trial_term + current_term + history_term +
                           spatial_term + source_term;
        const double scale = std::abs(trial_term) + std::abs(current_term) +
                             std::abs(history_term) + std::abs(spatial_term) +
                             std::abs(source_term);
        const double cancellation_reference =
            std::abs(residual_n[offset]) + std::abs(pressure_n[offset]) +
            std::abs(implicit_reference_n[offset]) +
            std::abs(implicit_reference_trial[offset]) +
            std::abs(pressure_final[offset]) +
            (stencil.order == MomentumTimeOrder::bdf2
                 ? std::abs(residual_nm1[offset]) +
                       std::abs(pressure_nm1[offset]) +
                       std::abs(implicit_reference_nm1[offset])
                 : 0.0);
        norm_sums[component_index * 3U] += raw * raw;
        norm_sums[component_index * 3U + 1U] += scale * scale;
        norm_sums[component_index * 3U + 2U] +=
            cancellation_reference * cancellation_reference;
        if (variable_density) {
          conservation_sums[component_index * 2U] += raw;
          conservation_sums[component_index * 2U + 1U] += scale;
        }
      }
    }
    mpi->allreduce_fp64_in_place(norm_sums.data(), norm_sums.size(),
                                 runtime::Fp64ReductionOperation::sum);
    if (variable_density)
      mpi->allreduce_fp64_in_place(
          conservation_sums.data(), conservation_sums.size(),
          runtime::Fp64ReductionOperation::sum);
    bool accepted = true;
    for (std::size_t component_index = 0U; component_index < 3U;
         ++component_index) {
      const double numerator = std::sqrt(norm_sums[component_index * 3U]);
      const double denominator =
          std::sqrt(norm_sums[component_index * 3U + 1U]);
      const double cancellation_reference =
          std::sqrt(norm_sums[component_index * 3U + 2U]);
      const double roundoff_limit = 512.0 *
                                    std::numeric_limits<double>::epsilon() *
                                    cancellation_reference;
      const double normalized =
          numerator <= roundoff_limit
              ? 0.0
              : (denominator > 0.0
                     ? numerator / denominator
                     : (numerator == 0.0
                            ? 0.0
                            : std::numeric_limits<double>::infinity()));
      report.final_momentum_normalized_l2[component_index] = normalized;
      if (variable_density) {
        const double conservation_numerator =
            std::abs(conservation_sums[component_index * 2U]);
        const double conservation_denominator =
            conservation_sums[component_index * 2U + 1U];
        const double conservation_roundoff =
            512.0 * std::numeric_limits<double>::epsilon() *
            conservation_denominator;
        report.final_momentum_relative_conservation_defect[component_index] =
            conservation_numerator <= conservation_roundoff
                ? 0.0
                : conservation_numerator /
                      std::max(conservation_denominator,
                               std::numeric_limits<double>::min());
      }
      accepted = accepted && std::isfinite(normalized) && normalized <= 1.0e-9;
    }
    return accepted;
  }

  void project_active_pressure_rhs(execution::VectorView<double> values) const {
    double totals[2]{};
    for (std::size_t row = 0U; row < owned_active_cells.size(); ++row) {
      const double volume = geometry->cell_volume_m3(owned_active_cells[row]);
      const double null_value = std::sqrt(volume);
      totals[0] += values[row] * null_value;
      totals[1] += volume;
    }
    mpi->allreduce_fp64_in_place(totals, 2U,
                                 runtime::Fp64ReductionOperation::sum);
    if (!(totals[1] > 0.0) || !std::isfinite(totals[0]) ||
        !std::isfinite(totals[1]))
      throw runtime::Error("immersed-flow pressure RHS projection is invalid");
    const double component = totals[0] / totals[1];
    for (std::size_t row = 0U; row < values.size(); ++row) {
      const double volume = geometry->cell_volume_m3(owned_active_cells[row]);
      values[row] -= std::sqrt(volume) * component;
    }
  }

  void normalize_active_pressure_solution(
      execution::VectorView<double> values) const {
    double totals[2]{};
    for (std::size_t row = 0U; row < owned_active_cells.size(); ++row) {
      const double volume = geometry->cell_volume_m3(owned_active_cells[row]);
      totals[0] += values[row] * std::sqrt(volume);
      totals[1] += volume;
    }
    mpi->allreduce_fp64_in_place(totals, 2U,
                                 runtime::Fp64ReductionOperation::sum);
    if (!(totals[1] > 0.0) || !std::isfinite(totals[0]) ||
        !std::isfinite(totals[1]))
      throw runtime::Error("immersed-flow pressure normalization is invalid");
    const double mean = totals[0] / totals[1];
    for (std::size_t row = 0U; row < values.size(); ++row) {
      const double volume = geometry->cell_volume_m3(owned_active_cells[row]);
      values[row] -= std::sqrt(volume) * mean;
    }
  }

  void require_canonical_inactive(FlowState &state,
                                  runtime::FieldStorage &layer,
                                  double rho_ref) const {
    const auto &access = state.solver_access_plan();
    const auto fields = state.fields();
    const auto density = layer.acquire_read<double>(
        access, kStatePhase, kStateActor, fields.density);
    const auto velocity = layer.acquire_read<double>(
        access, kStatePhase, kStateActor, fields.velocity);
    const auto pressure = layer.acquire_read<double>(
        access, kStatePhase, kStateActor, fields.mechanical_pressure);
    const auto positive_zero = [](double value) {
      return value == 0.0 && !std::signbit(value);
    };
    for (mesh::LocalCellId cell = 0U; cell < topology->owned_cell_count();
         ++cell) {
      const auto index = field_index(cell);
      if (domain->region(cell) == immersed::CellRegion::fluid) {
        if (density_setup.model == config::DensityModel::constant &&
            density(index.i, index.j, index.k, 0) != rho_ref)
          throw runtime::Error(
              "immersed-flow active constant density does not match physics");
        if (density_setup.model != config::DensityModel::constant &&
            (!(density(index.i, index.j, index.k, 0) > 0.0) ||
             !std::isfinite(density(index.i, index.j, index.k, 0))))
          throw runtime::Error(
              "immersed-flow active material density is invalid");
        continue;
      }
      if (!positive_zero(density(index.i, index.j, index.k, 0)) ||
          !positive_zero(pressure(index.i, index.j, index.k, 0)))
        throw runtime::Error(
            "immersed-flow inactive scalar slot is not canonical positive zero");
      for (int component_index = 0; component_index < 3; ++component_index)
        if (!positive_zero(
                velocity(index.i, index.j, index.k, component_index)))
          throw runtime::Error(
              "immersed-flow inactive velocity slot is not canonical positive zero");
      for (const auto field : fields.transported_cell_fields) {
        const auto transported = layer.acquire_read<double>(
            access, kStatePhase, kStateActor, field);
        if (!positive_zero(transported(index.i, index.j, index.k, 0)))
          throw runtime::Error(
              "immersed-flow inactive transport slot is not canonical positive zero");
      }
    }
    const auto face_velocity = layer.acquire_face_read<double>(
        access, kStatePhase, kStateActor, fields.face_velocity);
    const auto face_flux = layer.acquire_face_read<double>(
        access, kStatePhase, kStateActor, fields.face_mass_flux);
    for (mesh::LocalFaceId face = 0U; face < topology->local_face_count();
         ++face) {
      const auto neighbour = topology->neighbour(face);
      const bool owner_active =
          domain->region(topology->owner(face)) == immersed::CellRegion::fluid;
      const bool neighbour_active =
          neighbour.has_value() &&
          domain->region(*neighbour) == immersed::CellRegion::fluid;
      if (owner_active && (!neighbour.has_value() || neighbour_active))
        continue;
      if (!positive_zero(face_flux(face, 0)))
        throw runtime::Error(
            "immersed-flow inactive face flux is not canonical positive zero");
      for (int component_index = 0; component_index < 3; ++component_index)
        if (!positive_zero(face_velocity(face, component_index)))
          throw runtime::Error(
              "immersed-flow inactive face velocity is not canonical positive zero");
    }
  }

  std::vector<double> assemble_immersed_momentum_residual(
      FlowState &state, runtime::FieldStorage &accepted,
      double dynamic_viscosity,
      const std::vector<finite_volume::detail::ImmersedWallNormalGradient>
          *wall_normal_gradients = nullptr
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
      ,
      const finite_volume::ImmersedOperatorAdapter *operator_override = nullptr
#endif
  ) const {
    const auto &state_access = state.solver_access_plan();
    const auto fields = state.fields();
    const finite_volume::ReconstructionFieldBinding velocity_binding{
        accepted, state_access, kStatePhase, kStateActor, fields.velocity};
    const finite_volume::ReconstructionFieldBinding gradient_binding{
        *scratch->storage, *scratch->access, kImmersedFlowScratchPhase,
        kImmersedFlowScratchActor, scratch->velocity_gradient};
    reconstruction->compute_gradient(
        finite_volume::GradientScheme::weighted_least_squares,
        finite_volume::FiniteVolumeQuantity::velocity(), velocity_binding,
        gradient_binding);
    halo->exchange(*scratch->storage, scratch->velocity_gradient);
    const auto mass_flux = finite_volume::FaceMassFlux::acquire(
        state.solver_registry(), accepted, state_access, kStatePhase,
        kStateActor, fields.face_mass_flux, *topology);
    const finite_volume::ReconstructionFieldBinding face_binding{
        *scratch->storage, *scratch->access, kImmersedFlowScratchPhase,
        kImmersedFlowScratchActor, scratch->momentum_face};
    reconstruction->reconstruct_momentum_faces(mass_flux, velocity_binding,
                                               face_binding);
    auto viscosity = scratch->storage->acquire_face_write<double>(
        *scratch->access, kImmersedFlowScratchPhase, kImmersedFlowScratchActor,
        scratch->dynamic_viscosity);
    if (active_wale_authority != nullptr) {
      const auto &effective =
          active_wale_authority->effective_dynamic_viscosity_by_face();
      if (effective.size() != topology->local_face_count())
        throw runtime::Error(
            "immersed-flow WALE face-viscosity layout is invalid");
      for (mesh::LocalFaceId face = 0U; face < topology->local_face_count();
           ++face)
        viscosity(face, 0) = effective[face];
    } else {
      for (mesh::LocalFaceId face = 0U; face < topology->local_face_count();
           ++face)
        viscosity(face, 0) = dynamic_viscosity;
    }
    auto residual = scratch->storage->acquire_write<double>(
        *scratch->access, kImmersedFlowScratchPhase, kImmersedFlowScratchActor,
        scratch->momentum_residual);
    const auto extent = residual.interior_extent();
    for (int k = 0; k < extent.z; ++k)
      for (int j = 0; j < extent.y; ++j)
        for (int i = 0; i < extent.x; ++i)
          for (int component_index = 0; component_index < 3; ++component_index)
            residual(i, j, k, component_index) = 0.0;
    const auto face_velocity = scratch->storage->acquire_face_read<double>(
        *scratch->access, kImmersedFlowScratchPhase, kImmersedFlowScratchActor,
        scratch->momentum_face);
    const auto velocity = accepted.acquire_read<double>(
        state_access, kStatePhase, kStateActor, fields.velocity);
    const auto pressure = accepted.acquire_read<double>(
        state_access, kStatePhase, kStateActor, fields.mechanical_pressure);
    const auto gradient = scratch->storage->acquire_read<double>(
        *scratch->access, kImmersedFlowScratchPhase, kImmersedFlowScratchActor,
        scratch->velocity_gradient);
    const auto viscosity_read = scratch->storage->acquire_face_read<double>(
        *scratch->access, kImmersedFlowScratchPhase, kImmersedFlowScratchActor,
        scratch->dynamic_viscosity);
    const auto &selected_operator =
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
        operator_override != nullptr ? *operator_override : *immersed_operator;
#else
        *immersed_operator;
#endif
    if (wall_normal_gradients == nullptr) {
      selected_operator.accumulate_momentum(mass_flux, face_velocity, velocity,
                                            pressure, gradient, viscosity_read,
                                            residual);
    } else {
      finite_volume::detail::accumulate_momentum_with_wall_normal_constraints(
          selected_operator, mass_flux, face_velocity, velocity, pressure,
          *wall_normal_gradients, gradient, viscosity_read, residual);
    }
    physical_boundary_fvm->physical_boundary_momentum_contributions(
        *boundaries, mass_flux, face_velocity, velocity, gradient,
        viscosity_read, physical_boundary_momentum);
    physical_boundary_fvm->physical_boundary_pressure_contributions(
        *boundaries, pressure, physical_boundary_pressure);
    if (physical_boundary_momentum.size() != physical_boundary_pressure.size())
      throw runtime::Error("immersed-flow physical momentum boundary sets differ");
    for (std::size_t boundary_index = 0U;
         boundary_index < physical_boundary_momentum.size(); ++boundary_index) {
      const auto &momentum = physical_boundary_momentum[boundary_index];
      const auto &pressure_part = physical_boundary_pressure[boundary_index];
      if (momentum.global_face_id != pressure_part.global_face_id)
        throw runtime::Error(
            "immersed-flow physical momentum boundary identities differ");
      const auto face = topology->find_local_face(momentum.global_face_id);
      if (!face.has_value())
        throw runtime::Error(
            "immersed-flow physical momentum boundary face is unavailable");
      const auto owner = topology->owner(*face);
      const auto row = domain->active_cells().active_index(owner);
      if (!row.has_value())
        continue;
      if (*row >= owned_active_cells.size())
        throw runtime::Error(
            "immersed-flow physical momentum boundary row is not owned");
      const auto index = field_index(owner);
      for (std::size_t component_index = 0U; component_index < 3U;
           ++component_index)
        residual(index.i, index.j, index.k,
                 static_cast<int>(component_index)) +=
            momentum.convective[component_index] +
            momentum.viscous[component_index] +
            pressure_part.pressure[component_index];
    }
    std::vector<double> result(owned_active_cells.size() * 3U);
    for (std::size_t row = 0U; row < owned_active_cells.size(); ++row) {
      const auto index = field_index(owned_active_cells[row]);
      for (std::size_t component_index = 0U; component_index < 3U;
           ++component_index)
        result[row * 3U + component_index] = residual(
            index.i, index.j, index.k, static_cast<int>(component_index));
    }
    return result;
  }

  std::vector<double>
  apply_implicit_viscous_reference(FlowState &state,
                                   runtime::FieldStorage &layer) const {
    const auto &access = state.solver_access_plan();
    const auto fields = state.fields();
    const auto velocity = layer.acquire_read<double>(
        access, kStatePhase, kStateActor, fields.velocity);
    const std::size_t count = owned_active_cells.size();
    std::vector<double> result(count * 3U);
    for (std::size_t component_index = 0U; component_index < 3U;
         ++component_index) {
      auto input = momentum_solution[component_index]->view(0U, count);
      auto output = momentum_rhs[component_index]->view(0U, count);
      for (std::size_t row = 0U; row < count; ++row) {
        const auto index = field_index(owned_active_cells[row]);
        input[row] = velocity(index.i, index.j, index.k,
                              static_cast<int>(component_index));
      }
      momentum_operators[component_index]
          ->apply_spatial(static_cast<const execution::Buffer &>(
                              *momentum_solution[component_index])
                              .view(0U, count),
                          output)
          .wait();
      for (std::size_t row = 0U; row < count; ++row)
        result[row * 3U + component_index] = output[row];
    }
    return result;
  }

  struct FinalForceEvidence final {
    std::optional<ForceAttemptReport> report;
    std::optional<immersed::WallForceSample> sample;
    int failing_rank{-1};
  };

  FinalForceEvidence collect_final_force(
      FlowState &state, double dynamic_viscosity,
      const std::vector<finite_volume::detail::ImmersedWallNormalGradient>
          &final_wall_gradients) const {
    if (!wall_force_integrator.has_value())
      return {};

    const auto &access = state.solver_access_plan();
    const auto fields = state.fields();
    auto &trial = state.solver_layer(FlowLayer::trial);
    halo->exchange(trial, fields.velocity);
    halo->exchange(trial, fields.mechanical_pressure);
    {
      const finite_volume::detail::ForceAuthorityEvaluationScope
          force_authority;
      static_cast<void>(assemble_immersed_momentum_residual(
          state, trial, dynamic_viscosity, &final_wall_gradients));
    }

    const auto local_budget_reaction =
        immersed_operator->report().budget_reaction_N;
    std::array<double, 6> reduced_budget_reaction{
        local_budget_reaction.pressure[0], local_budget_reaction.pressure[1],
        local_budget_reaction.pressure[2], local_budget_reaction.viscous[0],
        local_budget_reaction.viscous[1],  local_budget_reaction.viscous[2]};
    const bool local_budget_reaction_ok =
        std::all_of(reduced_budget_reaction.begin(),
                    reduced_budget_reaction.end(),
                    [](double value) { return std::isfinite(value); });
    const auto budget_reaction_status = runtime::collective_status(
        *mpi, local_budget_reaction_ok,
        "immersed-flow final immersed budget reaction is non-finite");
    if (!budget_reaction_status.ok)
      return {std::nullopt, std::nullopt,
              budget_reaction_status.failing_rank};
    runtime::detail::check_mpi(
        MPI_Allreduce(MPI_IN_PLACE, reduced_budget_reaction.data(),
                      static_cast<int>(reduced_budget_reaction.size()),
                      MPI_DOUBLE, MPI_SUM, mpi->comm()),
        "MPI_Allreduce(immersed-flow final immersed budget reaction)");

    auto cell_viscosity = scratch->storage->acquire_write<double>(
        *scratch->access, kImmersedFlowScratchPhase, kImmersedFlowScratchActor,
        scratch->dynamic_viscosity_cell);
    const auto extent = cell_viscosity.interior_extent();
    const int ghost = cell_viscosity.ghost_width();
    for (int k = -ghost; k < extent.z + ghost; ++k)
      for (int j = -ghost; j < extent.y + ghost; ++j)
        for (int i = -ghost; i < extent.x + ghost; ++i)
          cell_viscosity(i, j, k, 0) = dynamic_viscosity;
    if (active_wale_authority != nullptr) {
      const auto &active = domain->active_cells();
      const auto &effective =
          active_wale_authority->effective_dynamic_viscosity_by_active_cell();
      if (effective.size() != active.ordered_global_ids().size())
        throw runtime::Error(
            "immersed-flow WALE wall-viscosity layout is invalid");
      for (std::size_t row = 0U; row < effective.size(); ++row) {
        const auto local =
            topology->find_local_cell(active.ordered_global_ids()[row]);
        if (!local.has_value())
          throw runtime::Error(
              "immersed-flow WALE wall-viscosity cell is unavailable");
        const auto index = field_index(*local);
        cell_viscosity(index.i, index.j, index.k, 0) = effective[row];
      }
    }

    const auto pressure = trial.acquire_read<double>(
        access, kStatePhase, kStateActor, fields.mechanical_pressure);
    const auto velocity = trial.acquire_read<double>(
        access, kStatePhase, kStateActor, fields.velocity);
    const auto gradient = scratch->storage->acquire_read<double>(
        *scratch->access, kImmersedFlowScratchPhase, kImmersedFlowScratchActor,
        scratch->velocity_gradient);
    const auto viscosity = scratch->storage->acquire_read<double>(
        *scratch->access, kImmersedFlowScratchPhase, kImmersedFlowScratchActor,
        scratch->dynamic_viscosity_cell);
    if (active_wale_authority != nullptr) {
      std::uint64_t actual_wall_fingerprint = kFnvOffset;
      for (const auto &wall : wall_links) {
        const auto index = field_index(wall.fluid);
        const double value = viscosity(index.i, index.j, index.k, 0);
        actual_wall_fingerprint =
            hash_u64(actual_wall_fingerprint, wall.link.id);
        actual_wall_fingerprint =
            hash_u64(actual_wall_fingerprint, fp64_bits(value));
      }
      if (actual_wall_fingerprint == 0U)
        actual_wall_fingerprint = 1U;
      const auto wall_viscosity_status = runtime::collective_status(
          *mpi,
          actual_wall_fingerprint ==
              active_wale_authority
                  ->wall_effective_viscosity_fingerprint(),
          "immersed-flow WALE wall viscosity lost attempt authority");
      if (!wall_viscosity_status.ok)
        return {std::nullopt, std::nullopt,
                wall_viscosity_status.failing_rank};
      pending_wall_effective_viscosity_fingerprint = actual_wall_fingerprint;
    }
    std::vector<immersed::detail::WallPressureNormalGradient> force_gradients;
    force_gradients.reserve(final_wall_gradients.size());
    for (const auto &wall_gradient : final_wall_gradients)
      force_gradients.push_back({wall_gradient.link, wall_gradient.value});
    const auto surface =
        immersed::detail::integrate_with_wall_pressure_authority(
            *wall_force_integrator, pressure, velocity, gradient, viscosity,
            force_gradients);
    if (surface.lowest_failing_rank >= 0)
      return {std::nullopt, std::nullopt, surface.lowest_failing_rank};

    immersed::ForceComponents budget_reaction;
    budget_reaction.pressure_N = {
        reduced_budget_reaction[0], reduced_budget_reaction[1],
        reduced_budget_reaction[2]};
    budget_reaction.viscous_N = {
        reduced_budget_reaction[3], reduced_budget_reaction[4],
        reduced_budget_reaction[5]};
    budget_reaction.total_N =
        add(budget_reaction.pressure_N, budget_reaction.viscous_N);
    const auto result = assemble_candidate_force_report_from_budget(
        budget_reaction, surface.surface_traction);
    const bool local_result_ok = finite_force(result.operator_force) &&
                                 finite_force(result.budget_reaction) &&
                                 finite_force(result.surface_traction) &&
                                 finite_force(result.consistency);
    const auto result_status = runtime::collective_status(
        *mpi, local_result_ok,
        "immersed-flow final immersed force evidence is non-finite");
    if (!result_status.ok)
      return {std::nullopt, std::nullopt, result_status.failing_rank};
    return {result, surface, -1};
  }

  void compute_immersed_pressure_gradient(
      FlowState &state, runtime::FieldStorage &source,
      const std::vector<finite_volume::detail::ImmersedWallNormalGradient>
          *wall_normal_gradients = nullptr
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
      ,
      const finite_volume::ImmersedOperatorAdapter *operator_override = nullptr
#endif
  ) const {
    const auto &state_access = state.solver_access_plan();
    const auto fields = state.fields();
    auto probe_flux = scratch->storage->acquire_face_write<double>(
        *scratch->access, kImmersedFlowScratchPhase, kImmersedFlowScratchActor,
        scratch->pressure_probe_mass_flux);
    auto viscosity = scratch->storage->acquire_face_write<double>(
        *scratch->access, kImmersedFlowScratchPhase, kImmersedFlowScratchActor,
        scratch->dynamic_viscosity);
    for (mesh::LocalFaceId face = 0U; face < topology->local_face_count();
         ++face) {
      probe_flux(face, 0) = 0.0;
      viscosity(face, 0) = 0.0;
    }
    auto pressure_residual_view = scratch->storage->acquire_write<double>(
        *scratch->access, kImmersedFlowScratchPhase, kImmersedFlowScratchActor,
        scratch->momentum_residual);
    const auto extent = pressure_residual_view.interior_extent();
    for (int k = 0; k < extent.z; ++k)
      for (int j = 0; j < extent.y; ++j)
        for (int i = 0; i < extent.x; ++i)
          for (int component_index = 0; component_index < 3; ++component_index)
            pressure_residual_view(i, j, k, component_index) = 0.0;
    const auto mass_flux = finite_volume::FaceMassFlux::acquire(
        scratch->registry, *scratch->storage, *scratch->access,
        kImmersedFlowScratchPhase, kImmersedFlowScratchActor,
        scratch->pressure_probe_mass_flux, *topology);
    const auto face_velocity = scratch->storage->acquire_face_read<double>(
        *scratch->access, kImmersedFlowScratchPhase, kImmersedFlowScratchActor,
        scratch->momentum_face);
    const auto velocity = source.acquire_read<double>(
        state_access, kStatePhase, kStateActor, fields.velocity);
    const auto pressure = source.acquire_read<double>(
        state_access, kStatePhase, kStateActor, fields.mechanical_pressure);
    const auto velocity_gradient = scratch->storage->acquire_read<double>(
        *scratch->access, kImmersedFlowScratchPhase, kImmersedFlowScratchActor,
        scratch->velocity_gradient);
    const auto viscosity_read = scratch->storage->acquire_face_read<double>(
        *scratch->access, kImmersedFlowScratchPhase, kImmersedFlowScratchActor,
        scratch->dynamic_viscosity);
    const auto &selected_operator =
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
        operator_override != nullptr ? *operator_override : *immersed_operator;
#else
        *immersed_operator;
#endif
    if (wall_normal_gradients == nullptr) {
      selected_operator.accumulate_momentum(
          mass_flux, face_velocity, velocity, pressure, velocity_gradient,
          viscosity_read, pressure_residual_view);
    } else {
      finite_volume::detail::accumulate_momentum_with_wall_normal_constraints(
          selected_operator, mass_flux, face_velocity, velocity, pressure,
          *wall_normal_gradients, velocity_gradient, viscosity_read,
          pressure_residual_view);
    }
    physical_boundary_fvm->physical_boundary_pressure_contributions(
        *boundaries, pressure, physical_boundary_pressure);
    for (const auto &pressure_part : physical_boundary_pressure) {
      const auto face = topology->find_local_face(pressure_part.global_face_id);
      if (!face.has_value())
        throw runtime::Error(
            "immersed-flow physical pressure boundary face is unavailable");
      const auto owner = topology->owner(*face);
      const auto row = domain->active_cells().active_index(owner);
      if (!row.has_value())
        continue;
      if (*row >= owned_active_cells.size())
        throw runtime::Error(
            "immersed-flow physical pressure boundary row is not owned");
      const auto index = field_index(owner);
      for (std::size_t component_index = 0U; component_index < 3U;
           ++component_index)
        pressure_residual_view(index.i, index.j, index.k,
                               static_cast<int>(component_index)) +=
            pressure_part.pressure[component_index];
    }
    auto gradient = scratch->storage->acquire_write<double>(
        *scratch->access, kImmersedFlowScratchPhase, kImmersedFlowScratchActor,
        scratch->pressure_gradient);
    for (mesh::LocalCellId cell = 0U; cell < topology->owned_cell_count();
         ++cell) {
      const auto index = field_index(cell);
      const double inverse_volume = 1.0 / geometry->cell_volume_m3(cell);
      for (int component_index = 0; component_index < 3; ++component_index)
        gradient(index.i, index.j, index.k, component_index) =
            pressure_residual_view(index.i, index.j, index.k, component_index) *
            inverse_volume;
    }
    halo->exchange(*scratch->storage, scratch->pressure_gradient);
  }

  std::vector<double> assemble_immersed_pressure_residual(
      FlowState &state, runtime::FieldStorage &source,
      const std::vector<finite_volume::detail::ImmersedWallNormalGradient>
          *wall_normal_gradients = nullptr
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
      ,
      const finite_volume::ImmersedOperatorAdapter *operator_override = nullptr
#endif
  ) const {
    compute_immersed_pressure_gradient(state, source, wall_normal_gradients
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
                                       ,
                                       operator_override
#endif
    );
    const auto gradient = scratch->storage->acquire_read<double>(
        *scratch->access, kImmersedFlowScratchPhase, kImmersedFlowScratchActor,
        scratch->pressure_gradient);
    std::vector<double> result(owned_active_cells.size() * 3U);
    for (std::size_t row = 0U; row < owned_active_cells.size(); ++row) {
      const auto cell = owned_active_cells[row];
      const auto index = field_index(cell);
      const double volume = geometry->cell_volume_m3(cell);
      for (std::size_t component_index = 0U; component_index < 3U;
           ++component_index)
        result[row * 3U + component_index] =
            volume * gradient(index.i, index.j, index.k,
                              static_cast<int>(component_index));
    }
    return result;
  }

  enum class WallPressureGradientState { current, corrected };

  std::vector<finite_volume::detail::ImmersedWallNormalGradient>
  wall_pressure_gradients(
      const std::vector<detail::ImmersedWallPressureCondition> &walls,
      WallPressureGradientState selected) const {
    if (walls.size() != wall_links.size())
      throw runtime::Error("immersed-flow wall pressure-gradient layout is invalid");
    std::vector<finite_volume::detail::ImmersedWallNormalGradient> result;
    result.reserve(walls.size());
    for (std::size_t index = 0U; index < walls.size(); ++index) {
      if (walls[index].link != wall_links[index].link.id)
        throw runtime::Error(
            "immersed-flow wall pressure-gradient link is inconsistent");
      const double value =
          selected == WallPressureGradientState::current
              ? walls[index].current_normal_gradient_pa_per_m
              : walls[index].corrected_normal_gradient_pa_per_m;
      if (!std::isfinite(value))
        throw runtime::Error(
            "immersed-flow wall pressure-gradient value is non-finite");
      result.push_back({walls[index].link, value});
    }
    return result;
  }

  std::vector<finite_volume::detail::ImmersedWallNormalGradient>
  reconstructed_wall_pressure_gradients(
      const FlowState &state, const runtime::FieldStorage &source) const {
    const auto pressure = source.acquire_read<double>(
        state.solver_access_plan(), kStatePhase, kStateActor,
        state.fields().mechanical_pressure);
    std::vector<finite_volume::detail::ImmersedWallNormalGradient> result;
    result.reserve(wall_links.size());
    for (const auto &wall : wall_links) {
      const auto &row_reconstruction =
          finite_volume::detail::ImmersedBoundaryAuthorityAccess::
              row_reconstruction(*immersed_operator, wall.link.id);
      const auto gradient = row_reconstruction.gradient(
          wall.link.wall_intercept_m, pressure, 0U);
      const double value = dot(gradient, wall.link.solid_to_fluid_normal);
      if (!std::isfinite(value))
        throw runtime::Error(
            "immersed-flow reconstructed wall pressure gradient is non-finite");
      result.push_back({wall.link.id, value});
    }
    return result;
  }

  std::vector<detail::ImmersedWallPressureCondition> assemble_face_predictor(
      FlowState &state, double rho_ref, const MomentumTimeStencil &stencil,
      const std::array<std::vector<double>, 3> &diagonal,
      const std::vector<finite_volume::detail::ImmersedWallNormalGradient>
          *current_wall_pressure_gradients = nullptr,
      const detail::ImmersedDensityAttemptAuthority *attempt_density = nullptr,
      const detail::ImmersedDensityAttemptAuthority *committed_density = nullptr,
      const detail::ImmersedDensityAttemptAuthority *history_density = nullptr)
      const {
    const bool variable_density = attempt_density != nullptr;
    if (variable_density != (committed_density != nullptr) ||
        variable_density != (history_density != nullptr))
      throw runtime::Error(
          "immersed predictor density authorities are incomplete");
    if (variable_density) {
      for (const auto *authority :
           {attempt_density, committed_density, history_density}) {
        if (authority->face_density.size() != topology->local_face_count() ||
            authority->wall_density.size() != wall_links.size() ||
            authority->owned_active_density.size() !=
                owned_active_cells.size())
          throw runtime::Error(
              "immersed predictor density authority layout is invalid");
        for (std::size_t index = 0U; index < wall_links.size(); ++index)
          if (authority->wall_density[index].link != wall_links[index].link.id)
            throw runtime::Error(
                "immersed predictor wall-density ordering is invalid");
      }
    }
    const auto rho_face = [&](const auto *authority,
                              mesh::LocalFaceId face) {
      return authority == nullptr ? rho_ref : authority->face_density[face];
    };
    const auto rho_wall = [&](const auto *authority, std::size_t wall) {
      return authority == nullptr
                 ? rho_ref
                 : authority->wall_density[wall].rho_wall_kg_per_m3;
    };
    if (current_wall_pressure_gradients != nullptr) {
      if (current_wall_pressure_gradients->size() != wall_links.size())
        throw runtime::Error(
            "immersed-flow current wall pressure-gradient layout is invalid");
      for (std::size_t index = 0U; index < wall_links.size(); ++index) {
        const auto &gradient = (*current_wall_pressure_gradients)[index];
        if (gradient.link != wall_links[index].link.id ||
            !std::isfinite(gradient.value))
          throw runtime::Error(
              "immersed-flow current wall pressure-gradient authority is invalid");
      }
    }
    const auto &access = state.solver_access_plan();
    const auto fields = state.fields();
    auto &trial = state.solver_layer(FlowLayer::trial);
    auto &committed = state.solver_layer(FlowLayer::committed);
    auto &history = state.solver_layer(FlowLayer::history);
    const auto velocity = trial.acquire_read<double>(
        access, kStatePhase, kStateActor, fields.velocity);
    const auto velocity_n = committed.acquire_read<double>(
        access, kStatePhase, kStateActor, fields.velocity);
    const auto velocity_nm1 = history.acquire_read<double>(
        access, kStatePhase, kStateActor, fields.velocity);
    const auto face_n = committed.acquire_face_read<double>(
        access, kStatePhase, kStateActor, fields.face_velocity);
    const auto face_nm1 = history.acquire_face_read<double>(
        access, kStatePhase, kStateActor, fields.face_velocity);
    halo->exchange(trial, fields.velocity);
    halo->exchange(trial, fields.mechanical_pressure);
    compute_immersed_pressure_gradient(state, trial,
                                       current_wall_pressure_gradients);
    const auto pressure = trial.acquire_read<double>(
        access, kStatePhase, kStateActor, fields.mechanical_pressure);
    const auto pressure_gradient = scratch->storage->acquire_read<double>(
        *scratch->access, kImmersedFlowScratchPhase, kImmersedFlowScratchActor,
        scratch->pressure_gradient);
    auto face_velocity = trial.acquire_face_write<double>(
        access, kStatePhase, kStateActor, fields.face_velocity);
    auto face_flux = trial.acquire_face_write<double>(
        access, kStatePhase, kStateActor, fields.face_mass_flux);

    for (mesh::LocalFaceId face = 0U; face < topology->local_face_count();
         ++face) {
      const auto neighbour = topology->neighbour(face);
      const auto owner = topology->owner(face);
      const auto owner_active = domain->active_cells().active_index(owner);
      const auto neighbour_active =
          neighbour.has_value()
              ? domain->active_cells().active_index(*neighbour)
              : std::optional<std::size_t>{};
      if (!neighbour.has_value() && owner_active.has_value()) {
        const auto patch = topology->patch_id(face);
        if (!patch.has_value())
          throw runtime::Error("immersed-flow physical face has no boundary patch");
        const auto owner_index = field_index(owner);
        const runtime::Real3 interior{
            velocity(owner_index.i, owner_index.j, owner_index.k, 0),
            velocity(owner_index.i, owner_index.j, owner_index.k, 1),
            velocity(owner_index.i, owner_index.j, owner_index.k, 2)};
        const auto area =
            geometry->face_area_vector_m2(face, mesh::FaceSide::owner);
        const auto values =
            boundaries->evaluate_velocity(*patch, interior, area);
        runtime::Real3 candidate = values.face;
        if (boundaries->patch(*patch).pressure_rule() ==
            boundary::PressureRule::prescribed_value) {
          const auto displacement = geometry->face_displacement_m(face);
          const double area_magnitude = std::sqrt(norm_squared(area));
          if (!(area_magnitude > 0.0) || !std::isfinite(area_magnitude))
            throw runtime::Error(
                "immersed-flow pressure-outlet face area is invalid");
          const auto unit_normal = multiply(1.0 / area_magnitude, area);
          const double normal_distance = dot(displacement, unit_normal);
          if (!(normal_distance > 0.0) || !std::isfinite(normal_distance))
            throw runtime::Error(
                "immersed-flow pressure-outlet normal distance is invalid");
          const auto pressure_values = boundaries->evaluate_pressure(
              *patch, pressure(owner_index.i, owner_index.j, owner_index.k, 0));
          const runtime::Real3 owner_gradient{
              pressure_gradient(owner_index.i, owner_index.j, owner_index.k, 0),
              pressure_gradient(owner_index.i, owner_index.j, owner_index.k, 1),
              pressure_gradient(owner_index.i, owner_index.j, owner_index.k,
                                2)};
          const double pressure_defect =
              (pressure_values.face -
               pressure(owner_index.i, owner_index.j, owner_index.k, 0) -
               dot(owner_gradient, displacement)) /
              normal_distance;
          if (!std::isfinite(pressure_defect))
            throw runtime::Error(
                "immersed-flow pressure-outlet interpolation defect is invalid");
          const auto row = *owner_active;
          if (row >= diagonal[0].size() || row >= diagonal[1].size() ||
              row >= diagonal[2].size())
            throw runtime::Error(
                "immersed-flow pressure-outlet active row is not owned");
          for (std::size_t component_index = 0U; component_index < 3U;
               ++component_index) {
            const auto component_offset = static_cast<int>(component_index);
            const double discrepancy_n =
                face_n(face, component_offset) -
                velocity_n(owner_index.i, owner_index.j, owner_index.k,
                           component_offset);
            double discrepancy_nm1 = 0.0;
            if (stencil.order == MomentumTimeOrder::bdf2)
              discrepancy_nm1 = face_nm1(face, component_offset) -
                                velocity_nm1(owner_index.i, owner_index.j,
                                             owner_index.k, component_offset);
            const double history_contribution =
                pressure_operator->face_velocity_mobility(face,
                                                          component_index) *
                (((-stencil.alpha1) * rho_face(committed_density, face) *
                  discrepancy_n) +
                 ((-stencil.alpha2) * rho_face(history_density, face) *
                  discrepancy_nm1)) /
                stencil.dt_s;
            if (component_index == 0U)
              candidate.x += history_contribution;
            else if (component_index == 1U)
              candidate.y += history_contribution;
            else
              candidate.z += history_contribution;
          }
          candidate.x -= pressure_operator->face_velocity_mobility(face, 0U) *
                         pressure_defect * unit_normal.x;
          candidate.y -= pressure_operator->face_velocity_mobility(face, 1U) *
                         pressure_defect * unit_normal.y;
          candidate.z -= pressure_operator->face_velocity_mobility(face, 2U) *
                         pressure_defect * unit_normal.z;
          if (!std::isfinite(candidate.x) || !std::isfinite(candidate.y) ||
              !std::isfinite(candidate.z))
            throw runtime::Error(
                "immersed-flow pressure-outlet predictor is non-finite");
        }
        face_velocity(face, 0) = candidate.x;
        face_velocity(face, 1) = candidate.y;
        face_velocity(face, 2) = candidate.z;
        const auto rule = boundaries->patch(*patch).mass_flux_rule();
        if (rule == boundary::MassFluxRule::identically_zero) {
          face_flux(face, 0) = 0.0;
        } else if (rule == boundary::MassFluxRule::prescribed_inlet_state ||
                   rule == boundary::MassFluxRule::outflow_only) {
          const double value =
              rho_face(attempt_density, face) * dot(candidate, area);
          if (!std::isfinite(value))
            throw runtime::Error(
                "immersed-flow physical boundary mass flux is non-finite");
          face_flux(face, 0) = value == 0.0 ? 0.0 : value;
        } else {
          throw runtime::Error(
              "immersed-flow physical face has an invalid mass-flux rule");
        }
        continue;
      }
      if (!owner_active.has_value() || !neighbour_active.has_value()) {
        face_flux(face, 0) = 0.0;
        for (int component_index = 0; component_index < 3; ++component_index)
          face_velocity(face, component_index) = 0.0;
        continue;
      }
      const auto owner_index = field_index(owner);
      const auto neighbour_index = field_index(*neighbour);
      const auto displacement = geometry->face_displacement_m(face);
      const auto area =
          geometry->face_area_vector_m2(face, mesh::FaceSide::owner);
      const double area_magnitude = std::sqrt(norm_squared(area));
      const runtime::Real3 unit_normal = multiply(1.0 / area_magnitude, area);
      runtime::Real3 pressure_gradient_face{};
      pressure_gradient_face.x =
          0.5 *
          (pressure_gradient(owner_index.i, owner_index.j, owner_index.k, 0) +
           pressure_gradient(neighbour_index.i, neighbour_index.j,
                             neighbour_index.k, 0));
      pressure_gradient_face.y =
          0.5 *
          (pressure_gradient(owner_index.i, owner_index.j, owner_index.k, 1) +
           pressure_gradient(neighbour_index.i, neighbour_index.j,
                             neighbour_index.k, 1));
      pressure_gradient_face.z =
          0.5 *
          (pressure_gradient(owner_index.i, owner_index.j, owner_index.k, 2) +
           pressure_gradient(neighbour_index.i, neighbour_index.j,
                             neighbour_index.k, 2));
      const double normal_distance = dot(displacement, unit_normal);
      const double pressure_defect =
          (pressure(neighbour_index.i, neighbour_index.j, neighbour_index.k,
                    0) -
           pressure(owner_index.i, owner_index.j, owner_index.k, 0) -
           dot(pressure_gradient_face, displacement)) /
          normal_distance;
      runtime::Real3 candidate{};
      for (std::size_t component_index = 0U; component_index < 3U;
           ++component_index) {
        const double predictor =
            0.5 *
            (velocity(owner_index.i, owner_index.j, owner_index.k,
                      static_cast<int>(component_index)) +
             velocity(neighbour_index.i, neighbour_index.j, neighbour_index.k,
                      static_cast<int>(component_index)));
        const double cell_n =
            0.5 *
            (velocity_n(owner_index.i, owner_index.j, owner_index.k,
                        static_cast<int>(component_index)) +
             velocity_n(neighbour_index.i, neighbour_index.j, neighbour_index.k,
                        static_cast<int>(component_index)));
        const double discrepancy_n =
            face_n(face, static_cast<int>(component_index)) - cell_n;
        double discrepancy_nm1 = 0.0;
        if (stencil.order == MomentumTimeOrder::bdf2) {
          const double cell_nm1 =
              0.5 * (velocity_nm1(owner_index.i, owner_index.j, owner_index.k,
                                  static_cast<int>(component_index)) +
                     velocity_nm1(neighbour_index.i, neighbour_index.j,
                                  neighbour_index.k,
                                  static_cast<int>(component_index)));
          discrepancy_nm1 =
              face_nm1(face, static_cast<int>(component_index)) - cell_nm1;
        }
        const double mobility =
            pressure_operator->face_velocity_mobility(face, component_index);
        const double value = predictor +
                             (-mobility * pressure_defect *
                              component(unit_normal, component_index)) +
                             mobility *
                                 (((-stencil.alpha1) *
                                   rho_face(committed_density, face) *
                                   discrepancy_n) +
                                  ((-stencil.alpha2) *
                                   rho_face(history_density, face) *
                                   discrepancy_nm1)) /
                                 stencil.dt_s;
        if (component_index == 0U)
          candidate.x = value;
        else if (component_index == 1U)
          candidate.y = value;
        else
          candidate.z = value;
        face_velocity(face, static_cast<int>(component_index)) = value;
      }
      double flux = rho_face(attempt_density, face) * dot(candidate, area);
      if (flux == 0.0)
        flux = 0.0;
      face_flux(face, 0) = flux;
    }

    {
      auto pressure_free_predictor = scratch->storage->acquire_write<double>(
          *scratch->access, kImmersedFlowScratchPhase, kImmersedFlowScratchActor,
          scratch->velocity_gradient);
      const auto extent = pressure_free_predictor.interior_extent();
      for (int k = 0; k < extent.z; ++k)
        for (int j = 0; j < extent.y; ++j)
          for (int i = 0; i < extent.x; ++i)
            for (int component_index = 0; component_index < 9;
                 ++component_index)
              pressure_free_predictor(i, j, k, component_index) = 0.0;
      for (std::size_t row = 0U; row < owned_active_cells.size(); ++row) {
        const auto cell = owned_active_cells[row];
        const auto index = field_index(cell);
        const double volume = geometry->cell_volume_m3(cell);
        for (std::size_t component_index = 0U; component_index < 3U;
             ++component_index) {
          const double mobility = volume / diagonal[component_index][row];
          pressure_free_predictor(index.i, index.j, index.k,
                                  static_cast<int>(component_index)) =
              velocity(index.i, index.j, index.k,
                       static_cast<int>(component_index)) +
              mobility * pressure_gradient(index.i, index.j, index.k,
                                           static_cast<int>(component_index));
        }
      }
    }
    halo->exchange(*scratch->storage, scratch->velocity_gradient);
    const auto pressure_free_predictor = scratch->storage->acquire_read<double>(
        *scratch->access, kImmersedFlowScratchPhase, kImmersedFlowScratchActor,
        scratch->velocity_gradient);

    std::vector<detail::ImmersedWallPressureCondition> conditions;
    conditions.reserve(wall_links.size());
    for (std::size_t wall_index = 0U; wall_index < wall_links.size();
         ++wall_index) {
      const auto &wall = wall_links[wall_index];
      const auto active_index = domain->active_cells().active_index(wall.fluid);
      if (!active_index.has_value())
        throw runtime::Error("immersed-flow wall fluid cell is inactive");
      runtime::Real3 predictor{};
      runtime::Real3 correction{};
      const auto &wall_reconstruction =
          finite_volume::detail::ImmersedBoundaryAuthorityAccess::
              row_reconstruction(*immersed_operator, wall.link.id);
      double wall_normal_gradient = 0.0;
      if (current_wall_pressure_gradients == nullptr) {
        const auto wall_gradient = wall_reconstruction.gradient(
            wall.link.wall_intercept_m, pressure, 0U);
        wall_normal_gradient =
            dot(wall_gradient, wall.link.solid_to_fluid_normal);
      } else {
        wall_normal_gradient =
            (*current_wall_pressure_gradients)[wall_index].value;
      }
      if (!std::isfinite(wall_normal_gradient))
        throw runtime::Error(
            "immersed-flow current wall pressure gradient is non-finite");
      for (std::size_t component_index = 0U; component_index < 3U;
           ++component_index) {
        double value =
            wall_reconstruction.value(wall.link.wall_intercept_m,
                                      pressure_free_predictor, component_index);
        const double normal_component =
            component(wall.link.solid_to_fluid_normal, component_index);
        const double reconstructed_velocity_n = wall_reconstruction.value(
            wall.link.wall_intercept_m, velocity_n, component_index);
        const double discrepancy_n =
            face_n(wall.face, static_cast<int>(component_index)) -
            reconstructed_velocity_n;
        double discrepancy_nm1 = 0.0;
        if (stencil.order == MomentumTimeOrder::bdf2) {
          const double reconstructed_velocity_nm1 = wall_reconstruction.value(
              wall.link.wall_intercept_m, velocity_nm1, component_index);
          discrepancy_nm1 =
              face_nm1(wall.face, static_cast<int>(component_index)) -
              reconstructed_velocity_nm1;
        }
        const double mobility = geometry->cell_volume_m3(wall.fluid) /
                                diagonal[component_index][*active_index];
        const double history_contribution =
            mobility *
            (((-stencil.alpha1) * rho_wall(committed_density, wall_index) *
              discrepancy_n) +
             ((-stencil.alpha2) * rho_wall(history_density, wall_index) *
              discrepancy_nm1)) /
            stencil.dt_s;
        value += history_contribution;
        const double wall_pressure =
            -mobility * wall_normal_gradient * normal_component;
        value += wall_pressure;
        if (component_index == 0U) {
          predictor.x = value;
          correction.x = mobility;
        } else if (component_index == 1U) {
          predictor.y = value;
          correction.y = mobility;
        } else {
          predictor.z = value;
          correction.z = mobility;
        }
      }
      const double wall_density = rho_wall(attempt_density, wall_index);
      const double mdot_star = wall_density * wall.effective_measure_m2 *
                               dot(predictor, wall.link.solid_to_fluid_normal);
      detail::ImmersedWallPressureInput input{wall.link.id,
                                              wall_density,
                                              wall.effective_measure_m2,
                                              wall.link.solid_to_fluid_normal,
                                              correction,
                                              wall_normal_gradient,
                                              mdot_star};
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
      if (mpi->rank() ==
          injected_wall_failure_rank.load(std::memory_order_relaxed)) {
        switch (injected_wall_failure.load(std::memory_order_relaxed)) {
        case WallInputFault::non_positive_density:
          input.rho_wall_kg_per_m3 = 0.0;
          break;
        case WallInputFault::non_finite_density:
          input.rho_wall_kg_per_m3 = std::numeric_limits<double>::quiet_NaN();
          break;
        case WallInputFault::non_positive_coefficient:
          input.effective_transformed_measure_m2 = 0.0;
          break;
        case WallInputFault::non_finite_coefficient:
          input.effective_transformed_measure_m2 =
              std::numeric_limits<double>::infinity();
          break;
        case WallInputFault::none:
        case WallInputFault::stale_preconditioner_revision:
          break;
        }
      }
#endif
      conditions.push_back(
          detail::make_immersed_wall_pressure_condition(input));
      face_flux(wall.face, 0) = 0.0;
      for (int component_index = 0; component_index < 3; ++component_index)
        face_velocity(wall.face, component_index) = 0.0;
    }
    return conditions;
  }

  void ensure_exact_predictor_workspace(const FlowState &template_state,
                                        double rho_ref) {
    if (exact_predictor_probe_state.has_value()) {
      const auto fields = template_state.fields();
      if (fields.density != exact_predictor_fields.density ||
          fields.velocity != exact_predictor_fields.velocity ||
          fields.mechanical_pressure !=
              exact_predictor_fields.mechanical_pressure ||
          fields.face_velocity != exact_predictor_fields.face_velocity ||
          fields.face_mass_flux != exact_predictor_fields.face_mass_flux ||
          fields.transported_cell_fields !=
              exact_predictor_fields.transported_cell_fields)
        throw runtime::Error(
            "immersed-flow exact predictor field layout changed");
      if (!exact_predictor_response_layout_matches(
              exact_predictor_response_workspace,
              owned_active_cells.size(), topology->local_face_count()) ||
          !exact_predictor_response_layout_matches(
              exact_affine_response_workspace, owned_active_cells.size(),
              topology->local_face_count()) ||
          !exact_predictor_response_layout_matches(
              exact_homogeneous_baseline_workspace,
              owned_active_cells.size(), topology->local_face_count()) ||
          !exact_predictor_response_layout_matches(
              exact_homogeneous_response_workspace,
              owned_active_cells.size(),
              topology->local_face_count()))
        throw runtime::Error(
            "immersed-flow exact predictor response workspace changed");
      return;
    }
    const auto box = topology->owned_global_box();
    const runtime::FieldLayoutSet layout{
        {box.end.x - box.begin.x, box.end.y - box.begin.y,
         box.end.z - box.begin.z},
        topology->local_face_count()};
    FlowLayerValues zeros;
    zeros.density.resize(topology->owned_cell_count(), rho_ref);
    zeros.velocity.resize(topology->owned_cell_count() * 3U, 0.0);
    zeros.mechanical_pressure.resize(topology->owned_cell_count(), 0.0);
    zeros.face_velocity.resize(topology->local_face_count() * 3U, 0.0);
    zeros.face_mass_flux.resize(topology->local_face_count(), 0.0);
    zeros.transported_cell_fields.resize(
        template_state.fields().transported_cell_fields.size());
    for (auto &transported : zeros.transported_cell_fields)
      transported.resize(topology->owned_cell_count(), 0.0);
    ExactPredictorResponse candidate_response;
    resize_exact_predictor_response(candidate_response,
                                    owned_active_cells.size(),
                                    topology->local_face_count());
    ExactPredictorResponse candidate_affine_response;
    resize_exact_predictor_response(candidate_affine_response,
                                    owned_active_cells.size(),
                                    topology->local_face_count());
    ExactPredictorResponse candidate_homogeneous_baseline;
    resize_exact_predictor_response(candidate_homogeneous_baseline,
                                    owned_active_cells.size(),
                                    topology->local_face_count());
    ExactPredictorResponse candidate_homogeneous_response;
    resize_exact_predictor_response(candidate_homogeneous_response,
                                    owned_active_cells.size(),
                                    topology->local_face_count());
    auto candidate_fields = template_state.fields();
    auto candidate_state = FlowState::create(
        template_state.solver_registry(), layout, candidate_fields,
        template_state.metadata());
    candidate_state.seed_accepted_layers(zeros, zeros);
    static_assert(std::is_nothrow_move_assignable_v<FlowFieldIds>);
    static_assert(std::is_nothrow_move_constructible_v<FlowState>);
    exact_predictor_fields = std::move(candidate_fields);
    exact_predictor_probe_state.emplace(std::move(candidate_state));
    exact_predictor_response_workspace = std::move(candidate_response);
    exact_affine_response_workspace =
        std::move(candidate_affine_response);
    exact_homogeneous_baseline_workspace =
        std::move(candidate_homogeneous_baseline);
    exact_homogeneous_response_workspace =
        std::move(candidate_homogeneous_response);
    ++exact_predictor_probe_creation_count;
  }

  const ExactPredictorResponse &exact_predictor_response(
      double rho_ref,
      const MomentumTimeStencil &stencil,
      const std::array<std::vector<double>, 3> &diagonal,
      const std::vector<double> &owned_pressure_pa,
      const std::vector<finite_volume::detail::ImmersedWallNormalGradient>
          &wall_gradients,
      const linear::SolveControl &momentum_control,
      ImmersedLinearSolvePhase solve_phase,
      std::uint32_t pressure_corrector_index,
      const detail::ImmersedDensityAttemptAuthority *attempt_density = nullptr,
      const detail::ImmersedDensityAttemptAuthority *committed_density =
          nullptr,
      const detail::ImmersedDensityAttemptAuthority *history_density =
          nullptr) {
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    const auto allocation_events_before =
        execution::allocation_counters().allocation_events;
#endif
    const std::size_t count = owned_active_cells.size();
    if (!(rho_ref > 0.0) || !std::isfinite(rho_ref) ||
        owned_pressure_pa.size() != count)
      throw runtime::Error(
          "immersed-flow exact predictor response input is invalid");
    for (const auto &component_diagonal : diagonal)
      if (component_diagonal.size() != count)
        throw runtime::Error(
            "immersed-flow exact predictor response diagonal is invalid");

    const std::array<std::size_t, 6> response_capacities_before{
        exact_predictor_response_workspace.divergence_per_volume.capacity(),
        exact_predictor_response_workspace.velocity_increment[0].capacity(),
        exact_predictor_response_workspace.velocity_increment[1].capacity(),
        exact_predictor_response_workspace.velocity_increment[2].capacity(),
        exact_predictor_response_workspace.face_mass_flux_increment.capacity(),
        exact_predictor_response_workspace.face_velocity_increment.capacity()};

    if (!exact_predictor_probe_state.has_value())
      throw runtime::Error(
          "immersed-flow exact predictor workspace is unavailable");
    auto &probe_state = *exact_predictor_probe_state;
    auto &trial = probe_state.solver_layer(FlowLayer::trial);
    const auto &access = probe_state.solver_access_plan();
    {
      auto pressure = trial.acquire_write<double>(
          access, kStatePhase, kStateActor,
          exact_predictor_fields.mechanical_pressure);
      const auto extent = pressure.interior_extent();
      for (int k = 0; k < extent.z; ++k)
        for (int j = 0; j < extent.y; ++j)
          for (int i = 0; i < extent.x; ++i)
            pressure(i, j, k, 0) = 0.0;
      for (std::size_t row = 0U; row < count; ++row) {
        const auto index = field_index(owned_active_cells[row]);
        pressure(index.i, index.j, index.k, 0) = owned_pressure_pa[row];
      }
    }
    halo->exchange(trial, exact_predictor_fields.mechanical_pressure);
    const auto pressure_response = assemble_immersed_pressure_residual(
        probe_state, trial, &wall_gradients);
    if (pressure_response.size() != count * 3U)
      throw runtime::Error(
          "immersed-flow exact predictor response pressure layout is invalid");

    auto &response = exact_predictor_response_workspace;
    response.work = {};
    response.work.response_count = 1U;
    for (std::size_t component_index = 0U; component_index < 3U;
         ++component_index) {
      auto rhs_values = momentum_rhs[component_index]->view(0U, count);
      auto solution_values = momentum_solution[component_index]->view(0U,
                                                                      count);
      for (std::size_t row = 0U; row < count; ++row) {
        rhs_values[row] = -pressure_response[row * 3U + component_index];
        solution_values[row] = 0.0;
      }
      momentum_preconditioners[component_index]->update(
          *momentum_operators[component_index],
          momentum_operators[component_index]->revision());
      const auto rhs_const =
          static_cast<const execution::Buffer &>(
              *momentum_rhs[component_index])
              .view(0U, count);
      auto inner_control = momentum_control;
      inner_control.atol = std::min(inner_control.atol, 1.0e-15);
      inner_control.rtol = std::min(inner_control.rtol, 1.0e-14);
      inner_control.max_iterations =
          std::max<std::uint64_t>(inner_control.max_iterations, 500U);
      const auto control =
          relative_equation_control(inner_control, rhs_const);
      auto solve = momentum_solver->solve(
          *momentum_operators[component_index],
          *momentum_preconditioners[component_index], rhs_const,
          solution_values, control);
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
      if (injected_linear_failure_phase.load(std::memory_order_relaxed) ==
              static_cast<int>(solve_phase) &&
          injected_linear_failure_component.load(std::memory_order_relaxed) ==
              component_index) {
        solve.reason = linear::SolveTerminationReason::maximum_iterations;
        solve.iterations = control.max_iterations;
        solve.recursive_residual =
            std::isfinite(solve.recursive_residual)
                ? std::max(1.0, solve.recursive_residual)
                : 1.0;
        solve.final_residual = solve.recursive_residual;
      }
#endif
      add_exact_predictor_solve(response.work, solve);
      const auto solve_status = runtime::collective_status(
          *mpi, solve_succeeded(solve.reason),
          "immersed-flow exact predictor response momentum solve failed");
      if (!solve_status.ok) {
        solve.lowest_failing_rank = solve_status.failing_rank;
        last_linear_solve_failure = ImmersedLinearSolveFailure{
            solve_phase, pressure_corrector_index,
            static_cast<std::uint32_t>(component_index), solve, 0.0, 0.0,
            0.0};
        throw std::pair{solve_status.failing_rank,
                        StepFailureReason::pressure_linear_solve};
      }
      const auto solved_values =
          static_cast<const execution::Buffer &>(
              *momentum_solution[component_index])
              .view(0U, count);
      response.velocity_increment[component_index].assign(
          solved_values.data(), solved_values.data() + count);
    }
    {
      auto velocity = trial.acquire_write<double>(
          access, kStatePhase, kStateActor, exact_predictor_fields.velocity);
      const auto extent = velocity.interior_extent();
      for (int k = 0; k < extent.z; ++k)
        for (int j = 0; j < extent.y; ++j)
          for (int i = 0; i < extent.x; ++i)
            for (int component_index = 0; component_index < 3;
                 ++component_index)
              velocity(i, j, k, component_index) = 0.0;
      for (std::size_t row = 0U; row < count; ++row) {
        const auto index = field_index(owned_active_cells[row]);
        for (std::size_t component_index = 0U; component_index < 3U;
             ++component_index)
          velocity(index.i, index.j, index.k,
                   static_cast<int>(component_index)) =
              response.velocity_increment[component_index][row];
      }
    }
    halo->exchange(trial, exact_predictor_fields.velocity);
    static_cast<void>(assemble_face_predictor(
        probe_state, rho_ref, stencil, diagonal, &wall_gradients,
        attempt_density, committed_density, history_density));
    const auto face_flux = trial.acquire_face_read<double>(
        access, kStatePhase, kStateActor,
        exact_predictor_fields.face_mass_flux);
    const auto face_velocity = trial.acquire_face_read<double>(
        access, kStatePhase, kStateActor, exact_predictor_fields.face_velocity);
    response.face_mass_flux_increment.resize(topology->local_face_count());
    response.face_velocity_increment.resize(topology->local_face_count() * 3U);
    response.divergence_per_volume.assign(count, 0.0);
    for (mesh::LocalFaceId face = 0U; face < topology->local_face_count();
         ++face) {
      const double value = face_flux(face, 0);
      response.face_mass_flux_increment[face] = value;
      for (std::size_t component_index = 0U; component_index < 3U;
           ++component_index)
        response.face_velocity_increment[face * 3U + component_index] =
            face_velocity(face, static_cast<int>(component_index));
      const auto owner_row =
          domain->active_cells().active_index(topology->owner(face));
      if (owner_row.has_value() && *owner_row < count)
        response.divergence_per_volume[*owner_row] += value;
      const auto neighbour = topology->neighbour(face);
      if (!neighbour.has_value() || topology->periodic_pair(face).has_value())
        continue;
      const auto neighbour_row =
          domain->active_cells().active_index(*neighbour);
      if (neighbour_row.has_value() && *neighbour_row < count)
        response.divergence_per_volume[*neighbour_row] -= value;
    }
    for (std::size_t row = 0U; row < count; ++row)
      response.divergence_per_volume[row] /=
          geometry->cell_volume_m3(owned_active_cells[row]);
    const std::array<std::size_t, 6> response_capacities_after{
        response.divergence_per_volume.capacity(),
        response.velocity_increment[0].capacity(),
        response.velocity_increment[1].capacity(),
        response.velocity_increment[2].capacity(),
        response.face_mass_flux_increment.capacity(),
        response.face_velocity_increment.capacity()};
    if (response_capacities_after != response_capacities_before)
      ++exact_predictor_response_capacity_growth_count;
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    const auto allocation_events_after =
        execution::allocation_counters().allocation_events;
    if (allocation_events_after < allocation_events_before)
      throw runtime::Error(
          "immersed-flow exact predictor allocation counter regressed");
    last_exact_predictor_response_allocation_event_count =
        allocation_events_after - allocation_events_before;
#endif
    return response;
  }

  bool assess_final_wall_penetration(
      FlowState &state,
      const std::vector<detail::ImmersedWallPressureCondition> &walls) const {
    if (walls.size() != wall_links.size())
      return false;
    const auto &access = state.solver_access_plan();
    const auto fields = state.fields();
    const auto &trial = state.solver_layer(FlowLayer::trial);
    const auto face_flux = trial.acquire_face_read<double>(
        access, kStatePhase, kStateActor, fields.face_mass_flux);
    for (std::size_t index = 0U; index < walls.size(); ++index) {
      const auto &condition = walls[index];
      const auto &wall = wall_links[index];
      const double derived = condition.predictor_mass_flux_kg_per_s -
                             condition.correction_coefficient *
                                 condition.correction_normal_gradient_pa_per_m;
      const double tolerance =
          64.0 * std::numeric_limits<double>::epsilon() *
          std::max(1.0, std::abs(condition.predictor_mass_flux_kg_per_s));
      if (condition.link != wall.link.id || !std::isfinite(derived) ||
          std::abs(derived) > tolerance ||
          condition.corrected_mass_flux_kg_per_s != 0.0 ||
          std::signbit(condition.corrected_mass_flux_kg_per_s) ||
          face_flux(wall.face, 0) != 0.0 ||
          std::signbit(face_flux(wall.face, 0))) {
        return false;
      }
    }
    return true;
  }

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  detail::ImmersedFlowSpatialEnergyData
  spatial_energy_terms(const FlowState &state, FlowLayer selected,
                       double rho_ref, double dynamic_viscosity) const {
    if (!(rho_ref > 0.0) || !std::isfinite(rho_ref) ||
        !(dynamic_viscosity >= 0.0) || !std::isfinite(dynamic_viscosity))
      throw runtime::Error("immersed-flow spatial-energy inputs are invalid");
    const auto values = state.snapshot(selected);
    const auto box = topology->owned_global_box();
    const runtime::Int3 local_extent{box.end.x - box.begin.x,
                                     box.end.y - box.begin.y,
                                     box.end.z - box.begin.z};
    auto probe_state = FlowState::create(
        state.solver_registry(), {local_extent, topology->local_face_count()},
        state.fields(), state.metadata());
    probe_state.seed_accepted_layers(values, values);
    auto &source = probe_state.solver_layer(FlowLayer::committed);
    const auto fields = probe_state.fields();
    halo->exchange(source, fields.velocity);
    halo->exchange(source, fields.mechanical_pressure);

    auto probe_operator = finite_volume::ImmersedOperatorAdapter::create(
        *topology, *geometry, *domain, *ghost_plan, *transform,
        *reconstruction);

    const auto total = assemble_immersed_momentum_residual(
        probe_state, source, dynamic_viscosity, nullptr, &probe_operator);
    const auto inviscid = assemble_immersed_momentum_residual(
        probe_state, source, 0.0, nullptr, &probe_operator);
    const auto pressure = assemble_immersed_pressure_residual(
        probe_state, source, nullptr, &probe_operator);
    const auto implicit_viscous_reference =
        apply_implicit_viscous_reference(probe_state, source);
    const std::size_t count = owned_active_cells.size();
    if (count > std::numeric_limits<std::size_t>::max() / 3U ||
        total.size() != count * 3U || inviscid.size() != total.size() ||
        pressure.size() != total.size() ||
        implicit_viscous_reference.size() != total.size())
      throw runtime::Error("immersed-flow spatial-energy layout is invalid");

    const auto &access = probe_state.solver_access_plan();
    const auto velocity = source.acquire_read<double>(
        access, kStatePhase, kStateActor, fields.velocity);
    const auto pressure_field = source.acquire_read<double>(
        access, kStatePhase, kStateActor, fields.mechanical_pressure);
    const auto face_flux = source.acquire_face_read<double>(
        access, kStatePhase, kStateActor, fields.face_mass_flux);
    std::vector<double> mass_divergence(count, 0.0);
    const auto accumulate_divergence = [&](mesh::LocalCellId cell,
                                           double increment) {
      const auto row = domain->active_cells().active_index(cell);
      if (!row.has_value() || *row >= count)
        return;
      const double candidate = mass_divergence[*row] + increment;
      if (!std::isfinite(candidate))
        throw runtime::Error(
            "immersed-flow spatial-energy mass divergence is non-finite");
      mass_divergence[*row] = candidate;
    };
    double wall_flux_linf = 0.0;
    double physical_flux_linf = 0.0;
    for (mesh::LocalFaceId face = 0U; face < topology->local_face_count();
         ++face) {
      const double value = face_flux(face, 0);
      if (!std::isfinite(value))
        throw runtime::Error("immersed-flow spatial-energy face flux is non-finite");
      accumulate_divergence(topology->owner(face), value);
      const auto neighbour = topology->neighbour(face);
      if (!neighbour.has_value()) {
        physical_flux_linf = std::max(physical_flux_linf, std::abs(value));
        continue;
      }
      const bool owner_active = domain->active_cells()
                                    .active_index(topology->owner(face))
                                    .has_value();
      const bool neighbour_active =
          domain->active_cells().active_index(*neighbour).has_value();
      if (owner_active != neighbour_active)
        wall_flux_linf = std::max(wall_flux_linf, std::abs(value));
      if (!topology->periodic_pair(face).has_value())
        accumulate_divergence(*neighbour, -value);
    }
    double boundary_fluxes[2]{wall_flux_linf, physical_flux_linf};
    mpi->allreduce_fp64_in_place(boundary_fluxes, 2U,
                                 runtime::Fp64ReductionOperation::maximum);

    struct Wire final {
      std::uint64_t global_cell_id{};
      double velocity[3]{};
      double pressure{};
      double mass_divergence{};
      double total[3]{};
      double convective[3]{};
      double pressure_residual[3]{};
      double viscous[3]{};
      double implicit_viscous_reference[3]{};
      double volume{};
    };
    static_assert(std::is_trivially_copyable_v<Wire>);
    std::vector<Wire> local(count);
    for (std::size_t row = 0U; row < count; ++row) {
      const auto cell = owned_active_cells[row];
      const auto index = field_index(cell);
      auto &wire = local[row];
      wire.global_cell_id = topology->global_cell_id(cell);
      wire.pressure = pressure_field(index.i, index.j, index.k, 0);
      wire.mass_divergence = mass_divergence[row];
      wire.volume = geometry->cell_volume_m3(cell);
      for (std::size_t component_index = 0U; component_index < 3U;
           ++component_index) {
        const auto offset = row * 3U + component_index;
        wire.velocity[component_index] = velocity(
            index.i, index.j, index.k, static_cast<int>(component_index));
        wire.total[component_index] = total[offset];
        wire.pressure_residual[component_index] = pressure[offset];
        wire.convective[component_index] = inviscid[offset] - pressure[offset];
        wire.viscous[component_index] = total[offset] - inviscid[offset];
        wire.implicit_viscous_reference[component_index] =
            implicit_viscous_reference[offset];
      }
    }
    if (local.size() >
        static_cast<std::size_t>(std::numeric_limits<int>::max()) /
            sizeof(Wire))
      throw runtime::Error("immersed-flow spatial-energy payload exceeds MPI int");
    const int local_bytes = static_cast<int>(local.size() * sizeof(Wire));
    std::vector<int> counts(static_cast<std::size_t>(mpi->size()));
    runtime::detail::check_mpi(MPI_Allgather(&local_bytes, 1, MPI_INT,
                                             counts.data(), 1, MPI_INT,
                                             mpi->comm()),
                               "MPI_Allgather(immersed-flow spatial-energy counts)");
    std::vector<int> offsets(counts.size(), 0);
    std::size_t total_bytes = 0U;
    for (std::size_t rank = 0U; rank < counts.size(); ++rank) {
      if (counts[rank] < 0 ||
          static_cast<std::size_t>(counts[rank]) >
              static_cast<std::size_t>(std::numeric_limits<int>::max()) -
                  total_bytes)
        throw runtime::Error("immersed-flow spatial-energy payload is unsupported");
      offsets[rank] = static_cast<int>(total_bytes);
      total_bytes += static_cast<std::size_t>(counts[rank]);
    }
    if (total_bytes % sizeof(Wire) != 0U)
      throw runtime::Error("immersed-flow spatial-energy payload size is invalid");
    std::vector<Wire> gathered(total_bytes / sizeof(Wire));
    runtime::detail::check_mpi(
        MPI_Allgatherv(local.data(), local_bytes, MPI_BYTE, gathered.data(),
                       counts.data(), offsets.data(), MPI_BYTE, mpi->comm()),
        "MPI_Allgatherv(immersed-flow spatial-energy values)");
    std::sort(gathered.begin(), gathered.end(),
              [](const auto &left, const auto &right) {
                return left.global_cell_id < right.global_cell_id;
              });
    if (gathered.empty() ||
        std::adjacent_find(gathered.begin(), gathered.end(),
                           [](const auto &left, const auto &right) {
                             return left.global_cell_id == right.global_cell_id;
                           }) != gathered.end())
      throw runtime::Error("immersed-flow spatial-energy layout is invalid");

    detail::ImmersedFlowSpatialEnergyData report;
    report.active_global_cell_ids.reserve(gathered.size());
    report.cell_volume_m3.reserve(gathered.size());
    report.velocity_m_per_s.reserve(gathered.size() * 3U);
    report.pressure_pa.reserve(gathered.size());
    report.mass_divergence_kg_per_s.reserve(gathered.size());
    report.total_residual_N.reserve(gathered.size() * 3U);
    report.convective_residual_N.reserve(gathered.size() * 3U);
    report.pressure_residual_N.reserve(gathered.size() * 3U);
    report.viscous_residual_N.reserve(gathered.size() * 3U);
    report.implicit_viscous_reference_residual_N.reserve(gathered.size() * 3U);
    double closure_square = 0.0;
    double divergence_square = 0.0;
    double maximum_closure = -1.0;
    for (const auto &cell : gathered) {
      if (!(cell.volume > 0.0) || !std::isfinite(cell.volume) ||
          !std::isfinite(cell.pressure) || !std::isfinite(cell.mass_divergence))
        throw runtime::Error("immersed-flow spatial-energy value is invalid");
      report.active_global_cell_ids.push_back(cell.global_cell_id);
      report.cell_volume_m3.push_back(cell.volume);
      report.pressure_pa.push_back(cell.pressure);
      report.mass_divergence_kg_per_s.push_back(cell.mass_divergence);
      divergence_square += cell.mass_divergence * cell.mass_divergence;
      report.mass_divergence_linf_kg_per_s = std::max(
          report.mass_divergence_linf_kg_per_s, std::abs(cell.mass_divergence));
      double speed_squared = 0.0;
      for (std::size_t component_index = 0U; component_index < 3U;
           ++component_index) {
        const double u = cell.velocity[component_index];
        const double total_value = cell.total[component_index];
        const double convective = cell.convective[component_index];
        const double pressure_value = cell.pressure_residual[component_index];
        const double viscous = cell.viscous[component_index];
        const double implicit_reference_value =
            cell.implicit_viscous_reference[component_index];
        const double parts = convective + pressure_value + viscous;
        const double closure = total_value - parts;
        if (!std::isfinite(u) || !std::isfinite(total_value) ||
            !std::isfinite(convective) || !std::isfinite(pressure_value) ||
            !std::isfinite(viscous) ||
            !std::isfinite(implicit_reference_value) || !std::isfinite(closure))
          throw runtime::Error("immersed-flow spatial-energy term is non-finite");
        report.velocity_m_per_s.push_back(u);
        report.total_residual_N.push_back(total_value);
        report.convective_residual_N.push_back(convective);
        report.pressure_residual_N.push_back(pressure_value);
        report.viscous_residual_N.push_back(viscous);
        report.implicit_viscous_reference_residual_N.push_back(
            implicit_reference_value);
        speed_squared += u * u;
        report.total_power_W += u * total_value;
        report.convective_power_W += u * convective;
        report.pressure_power_W += u * pressure_value;
        report.viscous_power_W += u * viscous;
        report.implicit_viscous_reference_power_W +=
            u * implicit_reference_value;
        closure_square += closure * closure;
        report.residual_closure_linf_N =
            std::max(report.residual_closure_linf_N, std::abs(closure));
        if (std::abs(closure) > maximum_closure) {
          maximum_closure = std::abs(closure);
          report.maximum_closure_global_cell_id = cell.global_cell_id;
          report.maximum_closure_component =
              static_cast<std::uint32_t>(component_index);
          report.maximum_closure_total_value_N = total_value;
          report.maximum_closure_parts_value_N = parts;
        }
      }
      report.kinetic_energy_J += 0.5 * rho_ref * cell.volume * speed_squared;
      report.centered_convective_power_W +=
          0.5 * speed_squared * cell.mass_divergence;
      report.pressure_continuity_power_W +=
          cell.pressure * cell.mass_divergence / rho_ref;
    }
    report.reconstruction_power_W =
        report.convective_power_W - report.centered_convective_power_W;
    report.pressure_adjoint_defect_W =
        report.pressure_power_W + report.pressure_continuity_power_W;
    report.residual_closure_l2_N = std::sqrt(closure_square);
    report.mass_divergence_l2_kg_per_s = std::sqrt(divergence_square);
    report.stationary_wall_flux_linf_kg_per_s = boundary_fluxes[0];
    report.physical_boundary_flux_linf_kg_per_s = boundary_fluxes[1];
    return report;
  }

  detail::ImmersedFlowPressureCorrectionRecord make_cell_pressure_correction_record(
      const std::vector<double> &exact_velocity_change,
      const std::vector<double> &momentum_operator_velocity_change,
      const std::vector<double> &lfp_pressure_residual_change,
      const std::vector<double> &pressure_before_pa,
      const std::vector<double> &pressure_correction_pa,
      const std::vector<double> &pressure_after_pa,
      const std::vector<double> &pressure_rhs_per_volume,
      const std::vector<double> &compact_pressure_action_per_volume,
      const std::vector<double> &interface_pressure_action_per_volume,
      const std::vector<double> &hybrid_pressure_action_per_volume,
      const std::vector<double> &affine_wall_source_per_volume,
      const std::vector<double> &compact_predictor_defect_per_volume,
      const std::vector<double> &hybrid_predictor_defect_per_volume,
      const std::vector<detail::ImmersedWallPressureCondition> &walls) const {
    const std::size_t local_count = owned_active_cells.size();
    if (local_count > std::numeric_limits<std::size_t>::max() / 3U ||
        exact_velocity_change.size() != local_count * 3U ||
        momentum_operator_velocity_change.size() != local_count * 3U ||
        lfp_pressure_residual_change.size() != local_count * 3U ||
        pressure_before_pa.size() != local_count ||
        pressure_correction_pa.size() != local_count ||
        pressure_after_pa.size() != local_count ||
        pressure_rhs_per_volume.size() != local_count ||
        compact_pressure_action_per_volume.size() != local_count ||
        interface_pressure_action_per_volume.size() != local_count ||
        hybrid_pressure_action_per_volume.size() != local_count ||
        affine_wall_source_per_volume.size() != local_count ||
        compact_predictor_defect_per_volume.size() != local_count ||
        hybrid_predictor_defect_per_volume.size() != local_count)
      throw runtime::Error(
          "immersed-flow cell pressure-correction capture is invalid");
    struct Wire final {
      std::uint64_t global_cell_id{};
      double velocity_change[3]{};
      double momentum[3]{};
      double pressure[3]{};
      double pressure_before{};
      double pressure_correction{};
      double pressure_after{};
      double pressure_rhs{};
      double compact_pressure_action{};
      double interface_pressure_action{};
      double hybrid_pressure_action{};
      double affine_wall_source{};
      double compact_predictor_defect{};
      double hybrid_predictor_defect{};
    };
    static_assert(std::is_trivially_copyable_v<Wire>);
    std::vector<Wire> local(local_count);
    for (std::size_t row = 0U; row < local_count; ++row) {
      local[row].global_cell_id =
          topology->global_cell_id(owned_active_cells[row]);
      for (std::size_t component_index = 0U; component_index < 3U;
           ++component_index) {
        const auto offset = row * 3U + component_index;
        local[row].velocity_change[component_index] =
            exact_velocity_change[offset];
        local[row].momentum[component_index] =
            momentum_operator_velocity_change[offset];
        local[row].pressure[component_index] =
            lfp_pressure_residual_change[offset];
      }
      local[row].pressure_before = pressure_before_pa[row];
      local[row].pressure_correction = pressure_correction_pa[row];
      local[row].pressure_after = pressure_after_pa[row];
      local[row].pressure_rhs = pressure_rhs_per_volume[row];
      local[row].compact_pressure_action =
          compact_pressure_action_per_volume[row];
      local[row].interface_pressure_action =
          interface_pressure_action_per_volume[row];
      local[row].hybrid_pressure_action =
          hybrid_pressure_action_per_volume[row];
      local[row].affine_wall_source = affine_wall_source_per_volume[row];
      local[row].compact_predictor_defect =
          compact_predictor_defect_per_volume[row];
      local[row].hybrid_predictor_defect =
          hybrid_predictor_defect_per_volume[row];
    }
    if (local.size() >
        static_cast<std::size_t>(std::numeric_limits<int>::max()) /
            sizeof(Wire))
      throw runtime::Error(
          "immersed-flow cell pressure-correction capture exceeds MPI int");
    const int local_bytes = static_cast<int>(local.size() * sizeof(Wire));
    std::vector<int> counts(static_cast<std::size_t>(mpi->size()));
    runtime::detail::check_mpi(
        MPI_Allgather(&local_bytes, 1, MPI_INT, counts.data(), 1, MPI_INT,
                      mpi->comm()),
        "MPI_Allgather(immersed-flow cell pressure-correction counts)");
    std::vector<int> offsets(counts.size(), 0);
    std::size_t total_bytes = 0U;
    for (std::size_t rank = 0U; rank < counts.size(); ++rank) {
      if (counts[rank] < 0 ||
          total_bytes >
              static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw runtime::Error(
            "immersed-flow cell pressure-correction payload is unsupported");
      offsets[rank] = static_cast<int>(total_bytes);
      total_bytes += static_cast<std::size_t>(counts[rank]);
    }
    if (total_bytes % sizeof(Wire) != 0U ||
        total_bytes > static_cast<std::size_t>(std::numeric_limits<int>::max()))
      throw runtime::Error(
          "immersed-flow cell pressure-correction payload size is invalid");
    std::vector<Wire> gathered(total_bytes / sizeof(Wire));
    runtime::detail::check_mpi(
        MPI_Allgatherv(local.data(), local_bytes, MPI_BYTE, gathered.data(),
                       counts.data(), offsets.data(), MPI_BYTE, mpi->comm()),
        "MPI_Allgatherv(immersed-flow cell pressure-correction values)");
    std::sort(gathered.begin(), gathered.end(),
              [](const auto &left, const auto &right) {
                return left.global_cell_id < right.global_cell_id;
              });
    if (gathered.empty() ||
        std::adjacent_find(gathered.begin(), gathered.end(),
                           [](const auto &left, const auto &right) {
                             return left.global_cell_id == right.global_cell_id;
                           }) != gathered.end())
      throw runtime::Error(
          "immersed-flow cell pressure-correction layout is invalid");

    detail::ImmersedFlowPressureCorrectionRecord record;
    record.active_global_cell_ids.reserve(gathered.size());
    record.exact_velocity_change.reserve(gathered.size() * 3U);
    record.momentum_operator_velocity_change.reserve(gathered.size() * 3U);
    record.lfp_pressure_residual_change.reserve(gathered.size() * 3U);
    record.pressure_before_pa.reserve(gathered.size());
    record.pressure_correction_pa.reserve(gathered.size());
    record.pressure_after_pa.reserve(gathered.size());
    record.pressure_rhs_per_volume.reserve(gathered.size());
    record.compact_pressure_action_per_volume.reserve(gathered.size());
    record.interface_pressure_action_per_volume.reserve(gathered.size());
    record.hybrid_pressure_action_per_volume.reserve(gathered.size());
    record.affine_wall_source_per_volume.reserve(gathered.size());
    record.compact_predictor_defect_per_volume.reserve(gathered.size());
    record.hybrid_predictor_defect_per_volume.reserve(gathered.size());
    double momentum_square = 0.0;
    double pressure_square = 0.0;
    double closure_square = 0.0;
    double maximum_closure = -1.0;
    for (const auto &cell : gathered) {
      if (!std::isfinite(cell.pressure_before) ||
          !std::isfinite(cell.pressure_correction) ||
          !std::isfinite(cell.pressure_after) ||
          !std::isfinite(cell.pressure_rhs) ||
          !std::isfinite(cell.compact_pressure_action) ||
          !std::isfinite(cell.interface_pressure_action) ||
          !std::isfinite(cell.hybrid_pressure_action) ||
          !std::isfinite(cell.affine_wall_source) ||
          !std::isfinite(cell.compact_predictor_defect) ||
          !std::isfinite(cell.hybrid_predictor_defect))
        throw runtime::Error("immersed-flow pressure-mode capture is non-finite");
      record.active_global_cell_ids.push_back(cell.global_cell_id);
      record.pressure_before_pa.push_back(cell.pressure_before);
      record.pressure_correction_pa.push_back(cell.pressure_correction);
      record.pressure_after_pa.push_back(cell.pressure_after);
      record.pressure_rhs_per_volume.push_back(cell.pressure_rhs);
      record.compact_pressure_action_per_volume.push_back(
          cell.compact_pressure_action);
      record.interface_pressure_action_per_volume.push_back(
          cell.interface_pressure_action);
      record.hybrid_pressure_action_per_volume.push_back(
          cell.hybrid_pressure_action);
      record.affine_wall_source_per_volume.push_back(cell.affine_wall_source);
      record.compact_predictor_defect_per_volume.push_back(
          cell.compact_predictor_defect);
      record.hybrid_predictor_defect_per_volume.push_back(
          cell.hybrid_predictor_defect);
      for (std::size_t component_index = 0U; component_index < 3U;
           ++component_index) {
        const double velocity_change = cell.velocity_change[component_index];
        const double momentum = cell.momentum[component_index];
        const double pressure = cell.pressure[component_index];
        const double closure = momentum + pressure;
        if (!std::isfinite(velocity_change) || !std::isfinite(momentum) ||
            !std::isfinite(pressure) ||
            !std::isfinite(closure))
          throw runtime::Error(
              "immersed-flow cell pressure-correction capture is non-finite");
        record.exact_velocity_change.push_back(velocity_change);
        record.momentum_operator_velocity_change.push_back(momentum);
        record.lfp_pressure_residual_change.push_back(pressure);
        momentum_square += momentum * momentum;
        pressure_square += pressure * pressure;
        closure_square += closure * closure;
        record.momentum_linf =
            std::max(record.momentum_linf, std::abs(momentum));
        record.pressure_linf =
            std::max(record.pressure_linf, std::abs(pressure));
        record.closure_linf = std::max(record.closure_linf, std::abs(closure));
        if (std::abs(closure) > maximum_closure) {
          maximum_closure = std::abs(closure);
          record.maximum_closure_global_cell_id = cell.global_cell_id;
          record.maximum_closure_component =
              static_cast<std::uint32_t>(component_index);
          record.maximum_closure_momentum_value = momentum;
          record.maximum_closure_pressure_value = pressure;
        }
      }
    }
    record.momentum_l2 = std::sqrt(momentum_square);
    record.pressure_l2 = std::sqrt(pressure_square);
    record.closure_l2 = std::sqrt(closure_square);
    double wall_sums[2]{};
    double wall_maxima[2]{};
    for (const auto &wall : walls) {
      wall_sums[0] +=
          wall.predictor_mass_flux_kg_per_s * wall.predictor_mass_flux_kg_per_s;
      wall_sums[1] += wall.correction_normal_gradient_pa_per_m *
                      wall.correction_normal_gradient_pa_per_m;
      wall_maxima[0] =
          std::max(wall_maxima[0], std::abs(wall.predictor_mass_flux_kg_per_s));
      wall_maxima[1] = std::max(
          wall_maxima[1], std::abs(wall.correction_normal_gradient_pa_per_m));
    }
    mpi->allreduce_fp64_in_place(wall_sums, 2U,
                                 runtime::Fp64ReductionOperation::sum);
    mpi->allreduce_fp64_in_place(wall_maxima, 2U,
                                 runtime::Fp64ReductionOperation::maximum);
    record.wall_predictor_mass_flux_l2_kg_per_s = std::sqrt(wall_sums[0]);
    record.wall_predictor_mass_flux_linf_kg_per_s = wall_maxima[0];
    record.wall_correction_gradient_l2_pa_per_m = std::sqrt(wall_sums[1]);
    record.wall_correction_gradient_linf_pa_per_m = wall_maxima[1];
    return record;
  }
#endif

  PressureCorrectionReport correct_active_pressure(
      FlowState &state, double rho_ref, const MomentumTimeStencil &stencil,
      const std::array<std::vector<double>, 3> &diagonal,
      const std::vector<detail::ImmersedWallPressureCondition> &walls,
      const linear::SolveControl &control,
      const linear::SolveControl &momentum_control,
      std::uint32_t pressure_corrector_index,
      const detail::ImmersedDensityAttemptAuthority *attempt_density = nullptr,
      const detail::ImmersedDensityAttemptAuthority *committed_density =
          nullptr,
      const detail::ImmersedDensityAttemptAuthority *history_density =
          nullptr) {
    PressureCorrectionReport report{};
    const std::size_t count = owned_active_cells.size();
    const bool variable_density = attempt_density != nullptr;
    if (variable_density != (committed_density != nullptr) ||
        variable_density != (history_density != nullptr))
      throw runtime::Error(
          "immersed pressure density authorities are incomplete");
    if (variable_density) {
      for (const auto *authority :
           {attempt_density, committed_density, history_density})
        if (authority->owned_active_density.size() != count ||
            authority->face_density.size() != topology->local_face_count() ||
            authority->wall_density.size() != wall_links.size())
          throw runtime::Error(
              "immersed pressure density authority layout is invalid");
    }
    auto rhs = pressure_rhs->view(0U, count);
    auto correction = pressure_solution->view(0U, count);
    if (count != 0U) {
      std::fill(rhs.data(), rhs.data() + count, 0.0);
      std::fill(correction.data(), correction.data() + count, 0.0);
    }
    const auto &access = state.solver_access_plan();
    const auto fields = state.fields();
    auto &trial = state.solver_layer(FlowLayer::trial);
    const auto flux = trial.acquire_face_read<double>(
        access, kStatePhase, kStateActor, fields.face_mass_flux);
    if (variable_density) {
      for (std::size_t row = 0U; row < count; ++row) {
        const double temporal =
            (stencil.alpha0 * attempt_density->owned_active_density[row] +
             stencil.alpha1 * committed_density->owned_active_density[row] +
             stencil.alpha2 * history_density->owned_active_density[row]) /
            stencil.dt_s;
        const double volume =
            geometry->cell_volume_m3(owned_active_cells[row]);
        if (!std::isfinite(temporal))
          throw runtime::Error(
              "immersed material pressure temporal residual is non-finite");
        rhs[row] -= temporal * std::sqrt(volume);
      }
    }
    for (mesh::LocalFaceId face = 0U; face < topology->local_face_count();
         ++face) {
      const auto neighbour = topology->neighbour(face);
      if (!neighbour.has_value()) {
        const auto owner_row =
            domain->active_cells().active_index(topology->owner(face));
        if (owner_row.has_value() && *owner_row < count) {
          const double volume = geometry->cell_volume_m3(topology->owner(face));
          rhs[*owner_row] -= flux(face, 0) / std::sqrt(volume);
        }
        continue;
      }
      const auto owner_index =
          domain->active_cells().active_index(topology->owner(face));
      const auto neighbour_index =
          domain->active_cells().active_index(*neighbour);
      if (!owner_index.has_value() || !neighbour_index.has_value())
        continue;
      const double value = flux(face, 0);
      if (*owner_index < count) {
        const double volume = geometry->cell_volume_m3(topology->owner(face));
        rhs[*owner_index] -= value / std::sqrt(volume);
      }
      if (!topology->periodic_pair(face).has_value() &&
          *neighbour_index < count) {
        const double volume = geometry->cell_volume_m3(*neighbour);
        rhs[*neighbour_index] += value / std::sqrt(volume);
      }
    }
    const auto current_wall_gradients =
        wall_pressure_gradients(walls, WallPressureGradientState::current);
    const auto corrected_wall_gradients =
        wall_pressure_gradients(walls, WallPressureGradientState::corrected);
    auto wall_gradient_increments = corrected_wall_gradients;
    for (std::size_t index = 0U; index < wall_gradient_increments.size();
         ++index) {
      if (wall_gradient_increments[index].link !=
          current_wall_gradients[index].link)
        throw runtime::Error(
            "immersed-flow exact pressure wall-gradient increment layout is invalid");
      wall_gradient_increments[index].value -=
          current_wall_gradients[index].value;
    }

    std::uint64_t diagonal_fingerprint = kFnvOffset;
    for (const auto &component_values : diagonal)
      for (const double value : component_values)
        diagonal_fingerprint = hash_u64(diagonal_fingerprint, fp64_bits(value));
    std::uint64_t local_wall_condition_fingerprint = 0U;
    for (const auto &wall : walls) {
      std::uint64_t fingerprint = kFnvOffset;
      fingerprint = hash_u64(fingerprint, wall.link);
      fingerprint =
          hash_u64(fingerprint, fp64_bits(wall.correction_coefficient));
      fingerprint = hash_u64(fingerprint,
                             fp64_bits(wall.current_normal_gradient_pa_per_m));
      fingerprint = hash_u64(
          fingerprint, fp64_bits(wall.correction_normal_gradient_pa_per_m));
      fingerprint = hash_u64(
          fingerprint, fp64_bits(wall.corrected_normal_gradient_pa_per_m));
      local_wall_condition_fingerprint ^= fingerprint;
    }
    std::uint64_t wall_condition_fingerprint{};
    runtime::detail::check_mpi(
        MPI_Allreduce(&local_wall_condition_fingerprint,
                      &wall_condition_fingerprint, 1, MPI_UINT64_T, MPI_BXOR,
                      mpi->comm()),
        "MPI_Allreduce(immersed-flow wall-condition fingerprint)");
    const std::uint64_t density_fingerprint =
        variable_density ? attempt_density->fingerprint : fp64_bits(rho_ref);
    const std::uint64_t next_dependency_identity =
        detail::make_immersed_pressure_revision(
            0U, density_fingerprint, diagonal_fingerprint,
            geometry_fingerprint,
            domain->active_cells().fingerprint(), wall_condition_fingerprint);
    double dependency_changed =
        next_dependency_identity != pressure_dependency_identity ? 1.0 : 0.0;
    mpi->allreduce_fp64_in_place(&dependency_changed, 1U,
                                 runtime::Fp64ReductionOperation::maximum);
    if (dependency_changed == 1.0) {
      if (pressure_revision == std::numeric_limits<std::uint64_t>::max())
        throw runtime::Error("immersed-flow pressure revision would wrap");
      ++pressure_revision;
    }
    pressure_dependency_identity = next_dependency_identity;
    pressure_operator->replace(pressure_revision);
    ensure_exact_predictor_workspace(state, rho_ref);
    std::vector<finite_volume::detail::ImmersedWallNormalGradient>
        zero_wall_gradients;
    zero_wall_gradients.reserve(wall_links.size());
    for (const auto &wall : wall_links)
      zero_wall_gradients.push_back({wall.link.id, 0.0});
    std::vector<double> zero_pressure(count, 0.0);
    const auto &evaluated_homogeneous_baseline = exact_predictor_response(
        rho_ref, stencil, diagonal, zero_pressure, zero_wall_gradients,
        momentum_control, ImmersedLinearSolvePhase::pressure_affine_momentum,
        pressure_corrector_index, attempt_density, committed_density,
        history_density);
    copy_exact_predictor_response(evaluated_homogeneous_baseline,
                                  exact_homogeneous_baseline_workspace);
    exact_pressure_operator->replace(
        pressure_revision,
        [this, rho_ref, stencil, diagonal, zero_wall_gradients,
         momentum_control, pressure_corrector_index, attempt_density,
         committed_density, history_density](const std::vector<double> &pressure)
            -> const ExactPredictorResponse & {
          const auto &raw_response = exact_predictor_response(
              rho_ref, stencil, diagonal, pressure, zero_wall_gradients,
              momentum_control,
              ImmersedLinearSolvePhase::pressure_homogeneous_momentum,
              pressure_corrector_index, attempt_density, committed_density,
              history_density);
          subtract_exact_predictor_response(
              raw_response, exact_homogeneous_baseline_workspace,
              exact_homogeneous_response_workspace);
          return exact_homogeneous_response_workspace;
        });
    const auto &evaluated_affine_response = exact_predictor_response(
        rho_ref, stencil, diagonal, zero_pressure, wall_gradient_increments,
        momentum_control,
        ImmersedLinearSolvePhase::pressure_affine_momentum,
        pressure_corrector_index, attempt_density, committed_density,
        history_density);
    copy_exact_predictor_response(evaluated_affine_response,
                                  exact_affine_response_workspace);
    const auto &exact_affine_response = exact_affine_response_workspace;
    auto affine_source = pressure_residual->view(0U, count);
    for (std::size_t row = 0U; row < count; ++row) {
      const double sqrt_volume =
          std::sqrt(geometry->cell_volume_m3(owned_active_cells[row]));
      affine_source[row] =
          exact_affine_response.divergence_per_volume[row] * sqrt_volume;
      rhs[row] -= affine_source[row];
    }
    if (!exact_pressure_operator->has_pressure_reference())
      project_active_pressure_rhs(rhs);
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    std::vector<double> pressure_rhs_per_volume(count);
    std::vector<double> cell_pressure_before_correction(count);
    {
      const auto pressure_before = trial.acquire_read<double>(
          access, kStatePhase, kStateActor, fields.mechanical_pressure);
      for (std::size_t row = 0U; row < count; ++row) {
        const auto index = field_index(owned_active_cells[row]);
        cell_pressure_before_correction[row] =
            pressure_before(index.i, index.j, index.k, 0);
      }
    }
    const auto compact_pressure_action_per_volume =
        evaluate_pressure_operator_values(
            *pressure_operator, *geometry, *execution, owned_active_cells,
            cell_pressure_before_correction);
    const auto hybrid_pressure_action_per_volume =
        evaluate_pressure_operator_values(
            *exact_pressure_operator, *geometry, *execution,
            owned_active_cells,
            cell_pressure_before_correction);
    std::vector<double> interface_pressure_action_per_volume(count);
    std::vector<double> affine_wall_source_per_volume(count);
    std::vector<double> compact_predictor_defect_per_volume(count);
    std::vector<double> hybrid_predictor_defect_per_volume(count);
    for (std::size_t row = 0U; row < count; ++row) {
      const double sqrt_volume =
          std::sqrt(geometry->cell_volume_m3(owned_active_cells[row]));
      pressure_rhs_per_volume[row] = rhs[row] / sqrt_volume;
      interface_pressure_action_per_volume[row] =
          hybrid_pressure_action_per_volume[row] -
          compact_pressure_action_per_volume[row];
      affine_wall_source_per_volume[row] = affine_source[row] / sqrt_volume;
      compact_predictor_defect_per_volume[row] =
          pressure_rhs_per_volume[row] +
          compact_pressure_action_per_volume[row];
      hybrid_predictor_defect_per_volume[row] =
          compact_predictor_defect_per_volume[row] +
          interface_pressure_action_per_volume[row];
    }
#endif
    bool preconditioner_ok = true;
    try {
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
      if (injected_wall_failure.load(std::memory_order_relaxed) ==
              WallInputFault::stale_preconditioner_revision &&
          mpi->rank() ==
              injected_wall_failure_rank.load(std::memory_order_relaxed)) {
        compact_schur_preconditioner->update(*exact_pressure_operator,
                                             pressure_revision - 1U);
      }
#endif
    } catch (...) {
      preconditioner_ok = false;
    }
    const auto preconditioner_status = runtime::collective_status(
        *mpi, preconditioner_ok,
        "immersed-flow pressure preconditioner revision is stale");
    if (!preconditioner_status.ok) {
      report.disposition = PressureCorrectionDisposition::non_retryable_failure;
      report.reason = StepFailureReason::invalid_input;
      report.lowest_failing_rank = preconditioner_status.failing_rank;
      report.accepted = false;
      return report;
    }
    compact_schur_preconditioner->update(*exact_pressure_operator,
                                         pressure_revision);
    double volume_scales[2]{};
    for (const auto cell : owned_active_cells) {
      const double sqrt_volume = std::sqrt(geometry->cell_volume_m3(cell));
      volume_scales[0] = std::max(volume_scales[0], sqrt_volume);
      volume_scales[1] = std::max(volume_scales[1], 1.0 / sqrt_volume);
    }
    mpi->allreduce_fp64_in_place(volume_scales, 2U,
                                 runtime::Fp64ReductionOperation::maximum);
    if (!(volume_scales[0] > 0.0) || !(volume_scales[1] > 0.0) ||
        !std::isfinite(volume_scales[0]) || !std::isfinite(volume_scales[1])) {
      throw runtime::Error("immersed-flow pressure volume scaling is invalid");
    }
    const double minimum_sqrt_volume = 1.0 / volume_scales[1];
    linear::SolveControl transformed_control = control;
    transformed_control.atol *= minimum_sqrt_volume;
    transformed_control.rtol *= minimum_sqrt_volume / volume_scales[0];
    report.solve = pressure_solver->solve(
        *exact_pressure_operator, *compact_schur_preconditioner,
        static_cast<const execution::Buffer &>(*pressure_rhs).view(0U, count),
        correction, transformed_control);
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    if (injected_linear_failure_phase.load(std::memory_order_relaxed) ==
        static_cast<int>(ImmersedLinearSolvePhase::pressure_outer)) {
      report.solve.reason =
          linear::SolveTerminationReason::maximum_iterations;
      report.solve.iterations = transformed_control.max_iterations;
      report.solve.recursive_residual =
          std::isfinite(report.solve.recursive_residual)
              ? std::max(1.0, report.solve.recursive_residual)
              : 1.0;
      report.solve.final_residual = report.solve.recursive_residual;
    }
#endif
    if (!exact_pressure_operator->has_pressure_reference())
      normalize_active_pressure_solution(correction);
    auto residual = pressure_residual->view(0U, count);
    exact_pressure_operator->apply(correction, residual).wait();
    const auto &exact_homogeneous_response =
        exact_pressure_operator->cached_response(
            static_cast<const execution::Buffer &>(*pressure_solution)
                .view(0U, count));
    for (std::size_t row = 0U; row < count; ++row) {
      const double sqrt_volume =
          std::sqrt(geometry->cell_volume_m3(owned_active_cells[row]));
      residual[row] = (rhs[row] - residual[row]) / sqrt_volume;
    }
    report.independent_residual_l2 =
        global_l2(static_cast<const execution::Buffer &>(*pressure_residual)
                      .view(0U, count));
    for (std::size_t row = 0U; row < count; ++row) {
      const double sqrt_volume =
          std::sqrt(geometry->cell_volume_m3(owned_active_cells[row]));
      residual[row] = rhs[row] / sqrt_volume;
    }
    report.rhs_l2 =
        global_l2(static_cast<const execution::Buffer &>(*pressure_residual)
                      .view(0U, count));
    const double threshold =
        std::max(control.atol, control.rtol * report.rhs_l2);
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    if (injected_linear_failure_phase.load(std::memory_order_relaxed) ==
        static_cast<int>(
            ImmersedLinearSolvePhase::pressure_independent_residual))
      report.independent_residual_l2 =
          std::max(1.0, std::nextafter(threshold,
                                       std::numeric_limits<double>::infinity()));
#endif
    report.accepted = solve_succeeded(report.solve.reason) &&
                      std::isfinite(report.independent_residual_l2) &&
                      report.independent_residual_l2 <= threshold;
    const auto accepted = runtime::collective_status(
        *mpi, report.accepted,
        "immersed-flow pressure correction failed its residual contract");
    if (!accepted.ok) {
      const auto compact_failure = compact_schur_preconditioner->last_failure();
      const auto phase =
          compact_failure.has_value()
              ? ImmersedLinearSolvePhase::pressure_compact_preconditioner
          : solve_succeeded(report.solve.reason)
              ? ImmersedLinearSolvePhase::pressure_independent_residual
              : ImmersedLinearSolvePhase::pressure_outer;
      auto diagnostic_solve =
          compact_failure.has_value() ? *compact_failure : report.solve;
      diagnostic_solve.lowest_failing_rank = accepted.failing_rank;
      report.solve.lowest_failing_rank = accepted.failing_rank;
      last_linear_solve_failure =
          ImmersedLinearSolveFailure{phase,
                                     pressure_corrector_index,
                                     std::numeric_limits<std::uint32_t>::max(),
                                     diagnostic_solve,
                                     report.independent_residual_l2,
                                     report.rhs_l2,
                                     threshold};
      report.disposition = PressureCorrectionDisposition::recoverable_failure;
      report.reason = StepFailureReason::pressure_linear_solve;
      report.lowest_failing_rank = accepted.failing_rank;
      report.accepted = false;
      return report;
    }
    add_exact_predictor_work(pending_exact_predictor_work,
                             exact_homogeneous_baseline_workspace.work);
    add_exact_predictor_work(pending_exact_predictor_work,
                             exact_affine_response.work);
    add_exact_predictor_work(pending_exact_predictor_work,
                             exact_pressure_operator->accumulated_work());
    const auto lfp_pressure_before_correction =
        assemble_immersed_pressure_residual(state, trial,
                                            &current_wall_gradients);
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    if (!pending_previous_lfp_pressure_after_correction.empty()) {
      if (pending_previous_lfp_pressure_after_correction.size() !=
          lfp_pressure_before_correction.size())
        throw runtime::Error(
            "immersed-flow inter-corrector pressure authority layout is invalid");
      double difference_square = 0.0;
      for (std::size_t offset = 0U;
           offset < lfp_pressure_before_correction.size(); ++offset) {
        const double difference =
            lfp_pressure_before_correction[offset] -
            pending_previous_lfp_pressure_after_correction[offset];
        difference_square += difference * difference;
      }
      mpi->allreduce_fp64_in_place(&difference_square, 1U,
                                   runtime::Fp64ReductionOperation::sum);
      pending_inter_corrector_authority_difference_l2 =
          std::sqrt(difference_square);
    }
#endif
    for (std::size_t row = 0U; row < count; ++row)
      residual[row] = correction[row];
    for (std::size_t row = 0U; row < count; ++row) {
      const double sqrt_volume =
          std::sqrt(geometry->cell_volume_m3(owned_active_cells[row]));
      correction[row] /= sqrt_volume;
    }
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    const std::vector<double> pressure_correction_pa(correction.data(),
                                                     correction.data() + count);
#endif

    auto pressure = trial.acquire_write<double>(
        access, kStatePhase, kStateActor, fields.mechanical_pressure);
    for (std::size_t row = 0U; row < count; ++row) {
      const auto index = field_index(owned_active_cells[row]);
      pressure(index.i, index.j, index.k, 0) += correction[row];
    }
    if (!exact_pressure_operator->has_pressure_reference()) {
      double pressure_totals[2]{};
      for (std::size_t row = 0U; row < count; ++row) {
        const auto cell = owned_active_cells[row];
        const auto index = field_index(cell);
        const double volume = geometry->cell_volume_m3(cell);
        pressure_totals[0] += pressure(index.i, index.j, index.k, 0) * volume;
        pressure_totals[1] += volume;
      }
      mpi->allreduce_fp64_in_place(pressure_totals, 2U,
                                   runtime::Fp64ReductionOperation::sum);
      const double pressure_mean = pressure_totals[0] / pressure_totals[1];
      for (const auto cell : owned_active_cells) {
        const auto index = field_index(cell);
        pressure(index.i, index.j, index.k, 0) -= pressure_mean;
      }
    }
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    std::vector<double> cell_pressure_after_correction(count);
    for (std::size_t row = 0U; row < count; ++row) {
      const auto index = field_index(owned_active_cells[row]);
      cell_pressure_after_correction[row] =
          pressure(index.i, index.j, index.k, 0);
    }
#endif
    halo->exchange(trial, fields.mechanical_pressure);
    const auto lfp_pressure_after_correction =
        assemble_immersed_pressure_residual(state, trial,
                                            &corrected_wall_gradients);
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    pending_previous_lfp_pressure_after_correction =
        lfp_pressure_after_correction;
#endif
    if (lfp_pressure_before_correction.size() != count * 3U ||
        lfp_pressure_after_correction.size() != count * 3U)
      throw runtime::Error("immersed-flow LFP pressure-correction layout is invalid");
    std::vector<double> lfp_pressure_residual_change(count * 3U);
    for (std::size_t offset = 0U; offset < lfp_pressure_residual_change.size();
         ++offset)
      lfp_pressure_residual_change[offset] =
          lfp_pressure_after_correction[offset] -
          lfp_pressure_before_correction[offset];
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    std::vector<double> exact_velocity_change(count * 3U, 0.0);
    std::vector<double> momentum_operator_velocity_change(count * 3U, 0.0);
#endif
    for (std::size_t component_index = 0U; component_index < 3U;
         ++component_index) {
      auto momentum_change_rhs = momentum_rhs[component_index]->view(0U, count);
      auto velocity_change =
          momentum_solution[component_index]->view(0U, count);
      for (std::size_t row = 0U; row < count; ++row)
        velocity_change[row] =
            exact_homogeneous_response
                .velocity_increment[component_index][row] +
            exact_affine_response.velocity_increment[component_index][row];
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
      momentum_operators[component_index]
          ->apply(static_cast<const execution::Buffer &>(
                      *momentum_solution[component_index])
                      .view(0U, count),
                  momentum_change_rhs)
          .wait();
      for (std::size_t row = 0U; row < count; ++row)
        momentum_operator_velocity_change[row * 3U + component_index] =
            momentum_change_rhs[row];
      for (std::size_t row = 0U; row < count; ++row)
        exact_velocity_change[row * 3U + component_index] =
            velocity_change[row];
#endif
    }
    {
      auto velocity = trial.acquire_write<double>(access, kStatePhase,
                                                  kStateActor, fields.velocity);
      for (std::size_t row = 0U; row < count; ++row) {
        const auto index = field_index(owned_active_cells[row]);
        for (std::size_t component_index = 0U; component_index < 3U;
             ++component_index) {
          const double increment =
              momentum_solution[component_index]->view(0U, count)[row];
          if (!std::isfinite(increment))
            throw runtime::Error(
                "immersed-flow LFP pressure-correction velocity is non-finite");
          velocity(index.i, index.j, index.k,
                   static_cast<int>(component_index)) += increment;
        }
      }
    }
    halo->exchange(trial, fields.velocity);
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    double applied_wall_gradient_norm = 0.0;
    for (const auto &wall : walls)
      applied_wall_gradient_norm += wall.correction_normal_gradient_pa_per_m *
                                    wall.correction_normal_gradient_pa_per_m;
    mpi->allreduce_fp64_in_place(&applied_wall_gradient_norm, 1U,
                                 runtime::Fp64ReductionOperation::sum);
    last_wall_pressure_gradient_application_norm =
        std::max(last_wall_pressure_gradient_application_norm,
                 std::sqrt(applied_wall_gradient_norm));
#endif
    auto face_flux = trial.acquire_face_write<double>(
        access, kStatePhase, kStateActor, fields.face_mass_flux);
    auto face_velocity = trial.acquire_face_write<double>(
        access, kStatePhase, kStateActor, fields.face_velocity);
    for (mesh::LocalFaceId face = 0U; face < topology->local_face_count();
         ++face) {
      face_flux(face, 0) +=
          exact_homogeneous_response.face_mass_flux_increment[face] +
          exact_affine_response.face_mass_flux_increment[face];
      for (std::size_t component_index = 0U; component_index < 3U;
           ++component_index)
        face_velocity(face, static_cast<int>(component_index)) +=
            exact_homogeneous_response
                .face_velocity_increment[face * 3U + component_index] +
            exact_affine_response
                .face_velocity_increment[face * 3U + component_index];
    }
    for (std::size_t wall_index = 0U; wall_index < wall_links.size();
         ++wall_index) {
      const auto &wall = wall_links[wall_index];
      face_flux(wall.face, 0) = 0.0;
      for (int component_index = 0; component_index < 3; ++component_index)
        face_velocity(wall.face, component_index) = 0.0;
      const auto row = domain->active_cells().active_index(wall.fluid);
      if (!row.has_value() || *row >= count)
        throw runtime::Error("immersed-flow wall correction row is invalid");
    }
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    pending_cell_pressure_correction_records.push_back(
        make_cell_pressure_correction_record(
            exact_velocity_change, momentum_operator_velocity_change,
            lfp_pressure_residual_change,
            cell_pressure_before_correction, pressure_correction_pa,
            cell_pressure_after_correction, pressure_rhs_per_volume,
            compact_pressure_action_per_volume,
            interface_pressure_action_per_volume,
            hybrid_pressure_action_per_volume,
            affine_wall_source_per_volume,
            compact_predictor_defect_per_volume,
            hybrid_predictor_defect_per_volume, walls));
#endif
    report.disposition = PressureCorrectionDisposition::accepted;
    report.reason = StepFailureReason::none;
    report.lowest_failing_rank = -1;
    report.accepted = true;
    return report;
  }

  StepAttemptReport
  attempt_immersed(FlowState &state, const ImmersedFlowPhysics &physics,
                   const MomentumTimeStencil &stencil,
                   const linear::SolveControl &momentum_control,
                   const linear::SolveControl &pressure_control,
                   std::optional<ForceAttemptReport> *accepted_force,
                   std::optional<les::WaleSummary> *accepted_wale,
                   std::optional<MaterialDensityTransportReport>
                       *accepted_material,
                   std::optional<MaterialDensityTransportReport>
                       *accepted_post_closure_material,
                   std::optional<IdealGasClosureReport> *accepted_closure,
                   MaterialTransportFailureReason *material_failure,
                   std::uint64_t *material_attempt_identity) {
    StepAttemptReport report{};
    report.attempted_dt_s = stencil.dt_s;
    if (accepted_force != nullptr)
      accepted_force->reset();
    if (accepted_wale != nullptr)
      accepted_wale->reset();
    if (accepted_material != nullptr)
      accepted_material->reset();
    if (accepted_post_closure_material != nullptr)
      accepted_post_closure_material->reset();
    if (accepted_closure != nullptr)
      accepted_closure->reset();
    struct AcceptedMaterialGuard final {
      std::optional<MaterialDensityTransportReport> *output{};
      ~AcceptedMaterialGuard() noexcept {
        if (output != nullptr)
          output->reset();
      }
    } accepted_material_guard{accepted_material};
    if (material_failure != nullptr)
      *material_failure = MaterialTransportFailureReason::none;
    if (material_attempt_identity != nullptr)
      *material_attempt_identity = 0U;
    auto *material_adapter =
        std::get_if<detail::ImmersedMaterialDensityAdapter>(&density_adapter);
    auto *ideal_adapter =
        std::get_if<detail::ImmersedIdealGasDensityAdapter>(&density_adapter);
    const bool variable_density_adapter =
        material_adapter != nullptr || ideal_adapter != nullptr;
    const auto material_field_count = [&] {
      return material_adapter != nullptr
                 ? material_adapter->material_field_count()
                 : ideal_adapter != nullptr
                       ? ideal_adapter->material_field_count()
                       : 0U;
    };
    if (variable_density_adapter) {
      const auto fields = static_cast<std::size_t>(
          material_field_count());
      report.final_transport_normalized_l2.assign(fields, 0.0);
      report.final_transport_relative_conservation_defect.assign(fields, 0.0);
    }
    struct DensityAdapterGuard final {
      detail::ImmersedMaterialDensityAdapter *material{};
      detail::ImmersedIdealGasDensityAdapter *ideal{};
      ~DensityAdapterGuard() noexcept {
        if (material != nullptr)
          material->rollback();
        if (ideal != nullptr)
          ideal->rollback();
      }
    } density_adapter_guard{material_adapter, ideal_adapter};
    last_wale_evaluation_count = 0U;
    last_wale_coefficient_identity = {};
    last_wall_effective_viscosity_fingerprint = 0U;
    pending_diagnostic_wall_force_sample.reset();
    pending_diagnostic_local_flow_pattern.reset();
    pending_wale_evaluation_count = 0U;
    pending_wall_effective_viscosity_fingerprint = 0U;
    last_linear_solve_failure.reset();
    std::fill(attempt_continuity_divergence.begin(),
              attempt_continuity_divergence.end(), 0.0);
    pending_exact_predictor_work = {};
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    last_wall_pressure_gradient_application_norm = 0.0;
    pending_cell_pressure_correction_records.clear();
    pending_previous_lfp_pressure_after_correction.clear();
    pending_inter_corrector_authority_difference_l2 = 0.0;
    pending_inter_step_authority_difference_l2 =
        std::numeric_limits<double>::quiet_NaN();
    pending_bdf2_history_authority_difference_l2 =
        std::numeric_limits<double>::quiet_NaN();
    pending_final_momentum_pressure_residual.clear();
    pending_final_pressure_wall_gradients.clear();
#endif
    bool attempt_active = false;
    try {
      if (variable_density_adapter) {
        bool local_prepare_ok = true;
        try {
          if (material_adapter != nullptr)
            material_adapter->prepare_attempt();
          else
            ideal_adapter->prepare_attempt();
        } catch (...) {
          local_prepare_ok = false;
        }
        if (material_attempt_identity != nullptr)
          *material_attempt_identity =
              material_adapter != nullptr
                  ? material_adapter->attempt_identity()
                  : ideal_adapter->attempt_identity();
        const auto prepare_status = runtime::collective_status(
            *mpi, local_prepare_ok,
            "immersed variable-density workspace preparation failed");
        if (!prepare_status.ok) {
          report.disposition =
              StepAttemptDisposition::non_retryable_failure;
          report.reason = StepFailureReason::invalid_input;
          report.lowest_failing_rank = prepare_status.failing_rank;
          report.suggested_dt_s = stencil.dt_s;
          return report;
        }
      }
      bool local_input_ok = true;
      try {
        const auto expected = make_momentum_time_stencil(
            stencil.order, stencil.dt_s, stencil.previous_dt_s);
        const bool ordinary_density_contract =
            density_setup.model != config::DensityModel::ideal_gas &&
            !physics.cp_J_per_kg_K.has_value() &&
            !physics.gas_constant_J_per_kg_K.has_value() &&
            !physics.thermodynamic_pressure_pa.has_value();
        const bool ideal_density_contract =
            density_setup.model == config::DensityModel::ideal_gas &&
            ideal_adapter != nullptr && physics.cp_J_per_kg_K.has_value() &&
            physics.gas_constant_J_per_kg_K.has_value() &&
            physics.thermodynamic_pressure_pa.has_value() &&
            fp64_bits(*physics.cp_J_per_kg_K) ==
                fp64_bits(ideal_adapter->cp_J_per_kg_K()) &&
            fp64_bits(*physics.gas_constant_J_per_kg_K) ==
                fp64_bits(ideal_adapter->gas_constant_J_per_kg_K()) &&
            fp64_bits(*physics.thermodynamic_pressure_pa) ==
                fp64_bits(ideal_adapter->configured_pressure_pa()) &&
            *physics.thermodynamic_pressure_pa > 0.0 &&
            std::isfinite(*physics.thermodynamic_pressure_pa);
        local_input_ok =
            expected.alpha0 == stencil.alpha0 &&
            expected.alpha1 == stencil.alpha1 &&
            expected.alpha2 == stencil.alpha2 &&
            physics.density_model == density_setup.model &&
            (density_setup.model == config::DensityModel::constant ||
             (density_setup.model == config::DensityModel::material &&
              domain != nullptr) ||
             (density_setup.model == config::DensityModel::ideal_gas &&
              domain != nullptr)) &&
            (ordinary_density_contract || ideal_density_contract) &&
            physics.rho_ref_kg_per_m3 > 0.0 &&
            std::isfinite(physics.rho_ref_kg_per_m3) &&
            physics.dynamic_viscosity_pa_s >= 0.0 &&
            std::isfinite(physics.dynamic_viscosity_pa_s) &&
            active_wale_authority == nullptr &&
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
            (manufactured_body_source.empty() ||
             (manufactured_body_source.size() ==
                  owned_active_cells.size() * 3U &&
              std::all_of(
                  manufactured_body_source.begin(),
                  manufactured_body_source.end(),
                  [](double value) { return std::isfinite(value); }))) &&
#endif
            !state.attempt_active();
      } catch (...) {
        local_input_ok = false;
      }
      const auto input_status = runtime::collective_status(
          *mpi, local_input_ok, "immersed-flow fixed-step input is invalid");
      if (!input_status.ok) {
        report.disposition = StepAttemptDisposition::non_retryable_failure;
        report.reason = StepFailureReason::invalid_input;
        report.lowest_failing_rank = input_status.failing_rank;
        report.suggested_dt_s = 0.0;
        return report;
      }
      synchronize_pressure_authority_binding(state);
      struct PressureAuthorityBindingGuard final {
        Impl *flow;
        const FlowState *state;
        ~PressureAuthorityBindingGuard() noexcept {
          flow->refresh_pressure_authority_binding(*state);
        }
      } pressure_authority_binding_guard{this, &state};
      bool local_begin_ok = true;
      try {
        state.begin_attempt();
        attempt_active = true;
      } catch (...) {
        local_begin_ok = false;
      }
      const auto begin_status = runtime::collective_status(
          *mpi, local_begin_ok, "immersed-flow attempt transaction could not begin");
      if (!begin_status.ok) {
        if (attempt_active)
          state.rollback_attempt();
        attempt_active = false;
        report.disposition = StepAttemptDisposition::non_retryable_failure;
        report.reason = StepFailureReason::invalid_input;
        report.lowest_failing_rank = begin_status.failing_rank;
        report.suggested_dt_s = 0.0;
        return report;
      }
      if (ideal_adapter != nullptr) {
        bool local_closure_begin_ok = true;
        int preflight_failure_rank = -1;
        try {
          ideal_adapter->begin_closure(state);
        } catch (const detail::DensityClosurePreflightFailure &failure) {
          local_closure_begin_ok = false;
          preflight_failure_rank = failure.failing_rank();
        } catch (const runtime::detail::MpiOperationError &) {
          throw;
        } catch (...) {
          local_closure_begin_ok = false;
        }
        if (preflight_failure_rank >= 0) {
          state.rollback_attempt();
          attempt_active = false;
          report.disposition =
              StepAttemptDisposition::non_retryable_failure;
          report.reason = StepFailureReason::invalid_input;
          report.lowest_failing_rank = preflight_failure_rank;
          report.suggested_dt_s = stencil.dt_s;
          return report;
        }
        const auto closure_begin_status = runtime::collective_status(
            *mpi, local_closure_begin_ok,
            "immersed ideal-gas closure transaction could not begin");
        if (!closure_begin_status.ok) {
          state.rollback_attempt();
          attempt_active = false;
          report.disposition =
              StepAttemptDisposition::non_retryable_failure;
          report.reason = StepFailureReason::invalid_input;
          report.lowest_failing_rank = closure_begin_status.failing_rank;
          report.suggested_dt_s = stencil.dt_s;
          return report;
        }
      }
      const auto &access = state.solver_access_plan();
      const auto fields = state.fields();
      auto &trial = state.solver_layer(FlowLayer::trial);
      auto &committed = state.solver_layer(FlowLayer::committed);
      auto &history = state.solver_layer(FlowLayer::history);
      bool inactive_state_ok = true;
      try {
        require_canonical_inactive(state, history, physics.rho_ref_kg_per_m3);
        require_canonical_inactive(state, committed, physics.rho_ref_kg_per_m3);
        require_canonical_inactive(state, trial, physics.rho_ref_kg_per_m3);
      } catch (...) {
        inactive_state_ok = false;
      }
      const auto inactive_status = runtime::collective_status(
          *mpi, inactive_state_ok, "immersed-flow inactive state invariant failed");
      if (!inactive_status.ok) {
        state.rollback_attempt();
        attempt_active = false;
        report.disposition = StepAttemptDisposition::non_retryable_failure;
        report.reason = StepFailureReason::invalid_input;
        report.lowest_failing_rank = inactive_status.failing_rank;
        report.suggested_dt_s = 0.0;
        return report;
      }
      for (auto *accepted : {&committed, &history}) {
        halo->exchange(*accepted, fields.density);
        halo->exchange(*accepted, fields.velocity);
        halo->exchange(*accepted, fields.mechanical_pressure);
        for (const auto field : fields.transported_cell_fields)
          halo->exchange(*accepted, field);
      }
      if (variable_density_adapter) {
        {
          const auto current = committed.acquire_face_read<double>(
              access, kStatePhase, kStateActor, fields.face_mass_flux);
          const auto previous = history.acquire_face_read<double>(
              access, kStatePhase, kStateActor, fields.face_mass_flux);
          auto predictor = trial.acquire_face_write<double>(
              access, kStatePhase, kStateActor, fields.face_mass_flux);
          for (mesh::LocalFaceId face = 0U;
               face < topology->local_face_count(); ++face) {
            const auto neighbour = topology->neighbour(face);
            const bool owner_active =
                domain->region(topology->owner(face)) ==
                immersed::CellRegion::fluid;
            const bool neighbour_active =
                neighbour.has_value() &&
                domain->region(*neighbour) == immersed::CellRegion::fluid;
            const bool physical_active =
                !neighbour.has_value() && owner_active;
            double value = 0.0;
            if ((owner_active && neighbour_active) || physical_active) {
              value = stencil.order == MomentumTimeOrder::backward_euler
                          ? current(face, 0)
                          : 2.0 * current(face, 0) - previous(face, 0);
            }
            if (!std::isfinite(value))
              throw runtime::Error(
                  "immersed material predictor mass flux is non-finite");
            predictor(face, 0) = value == 0.0 ? 0.0 : value;
          }
        }
        const auto stage = material_adapter != nullptr
                               ? material_adapter->stage_predictor_transport(
                                     state, stencil)
                               : ideal_adapter->stage_predictor_transport(
                                     state, stencil);
        if (stage.reason != MaterialTransportFailureReason::none) {
          if (material_failure != nullptr)
            *material_failure = stage.reason;
          state.rollback_attempt();
          attempt_active = false;
          report.reason = material_step_failure(stage.reason);
          report.lowest_failing_rank = stage.lowest_failing_rank;
          const bool non_retryable =
              stage.reason == MaterialTransportFailureReason::invalid_input ||
              stage.reason ==
                  MaterialTransportFailureReason::collective_operation;
          report.disposition =
              non_retryable
                  ? StepAttemptDisposition::non_retryable_failure
                  : StepAttemptDisposition::recoverable_failure;
          report.suggested_dt_s =
              non_retryable ? stencil.dt_s : 0.5 * stencil.dt_s;
          return report;
        }
        if (ideal_adapter != nullptr) {
          const auto closure_stage = ideal_adapter->evaluate(
              state, detail::DensityClosureStage::predictor);
          if (accepted_closure != nullptr)
            *accepted_closure = ideal_adapter->take_closure_report();
          if (!closure_stage.accepted) {
            state.rollback_attempt();
            attempt_active = false;
            report.reason = StepFailureReason::density_closure_failure;
            report.lowest_failing_rank =
                closure_stage.lowest_failing_rank;
            report.disposition =
                closure_stage.recoverable
                    ? StepAttemptDisposition::recoverable_failure
                    : StepAttemptDisposition::non_retryable_failure;
            report.suggested_dt_s = closure_stage.recoverable
                                        ? 0.5 * stencil.dt_s
                                        : stencil.dt_s;
            return report;
          }
        }
        halo->exchange(trial, fields.density);
        for (const auto field : fields.transported_cell_fields)
          halo->exchange(trial, field);
        if (ideal_adapter != nullptr)
          ideal_adapter->after_halo(
              detail::DensityClosureStage::predictor);
      }
      std::optional<detail::ImmersedDensityAttemptAuthority>
          material_density_trial;
      std::optional<detail::ImmersedDensityAttemptAuthority>
          material_density_committed;
      std::optional<detail::ImmersedDensityAttemptAuthority>
          material_density_history;
      const auto require_density_authority = [&](auto &&function) {
        const auto density_stage = prepare_density_authority_collectively(
            std::forward<decltype(function)>(function));
        if (density_stage.reason == MaterialTransportFailureReason::none)
          return;
        if (material_failure != nullptr)
          *material_failure = density_stage.reason;
        throw std::pair{density_stage.lowest_failing_rank,
                        material_step_failure(density_stage.reason)};
      };
      if (variable_density_adapter) {
        require_density_authority([&] {
          material_density_trial.emplace(density_authority(state, trial));
          material_density_committed.emplace(
              density_authority(state, committed));
          material_density_history.emplace(density_authority(state, history));
        });
      }
      std::optional<detail::ImmersedWaleAttemptAuthority> wale_authority;
      struct WaleAuthorityBindingGuard final {
        Impl *flow;
        ~WaleAuthorityBindingGuard() noexcept {
          flow->active_wale_authority = nullptr;
        }
      } wale_authority_binding_guard{this};
      if (wale != nullptr) {
        const auto rho_attempt =
            (!variable_density_adapter ? committed : trial)
                .acquire_read<double>(access, kStatePhase, kStateActor,
                                      fields.density);
        wale_authority.emplace(prepare_wale_authority(
            state, physics.dynamic_viscosity_pa_s, stencil, rho_attempt));
        active_wale_authority = &*wale_authority;
      }
      const std::size_t count = owned_active_cells.size();
      std::array<std::vector<double>, 3> diagonal;
      for (auto &values : diagonal)
        values.resize(count);
      bool local_momentum_preparation_ok = true;
      try {
        for (std::size_t component_index = 0U; component_index < 3U;
             ++component_index) {
          for (std::size_t row = 0U; row < count; ++row) {
            const auto cell = owned_active_cells[row];
            const double volume = geometry->cell_volume_m3(cell);
            double diagonal_scale = 1.0;
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
            diagonal_scale = injected_momentum_time_diagonal_scale.load(
                std::memory_order_relaxed);
            if (!(diagonal_scale > 0.0) || !std::isfinite(diagonal_scale))
              throw runtime::Error(
                  "immersed-flow momentum diagonal test scale is invalid");
#endif
            const double rho_attempt =
                material_density_trial.has_value()
                    ? material_density_trial->owned_active_density[row]
                    : physics.rho_ref_kg_per_m3;
            diagonal[component_index][row] =
                diagonal_scale * stencil.alpha0 * rho_attempt * volume /
                stencil.dt_s;
          }
          if (active_wale_authority == nullptr) {
            momentum_operators[component_index]->replace(
                diagonal[component_index], physics.dynamic_viscosity_pa_s);
          } else {
            momentum_operators[component_index]->replace(
                diagonal[component_index],
                active_wale_authority
                    ->effective_dynamic_viscosity_by_active_cell(),
                active_wale_authority
                    ->effective_dynamic_viscosity_by_face());
          }
          diagonal[component_index] =
              momentum_operators[component_index]->assembled_diagonal();
        }
      } catch (const runtime::detail::MpiOperationError &) {
        throw;
      } catch (...) {
        local_momentum_preparation_ok = false;
      }
      const auto momentum_preparation_status = runtime::collective_status(
          *mpi, local_momentum_preparation_ok,
          "immersed-flow momentum equation preparation failed");
      if (!momentum_preparation_status.ok) {
        state.rollback_attempt();
        attempt_active = false;
        report.disposition = StepAttemptDisposition::non_retryable_failure;
        report.reason = StepFailureReason::invalid_input;
        report.lowest_failing_rank = momentum_preparation_status.failing_rank;
        report.suggested_dt_s = 0.0;
        return report;
      }
      std::vector<finite_volume::detail::ImmersedWallNormalGradient>
          reconstructed_committed_pressure_authority;
      if (!committed_pressure_authority_available)
        reconstructed_committed_pressure_authority =
            reconstructed_wall_pressure_gradients(state, committed);
      const auto *committed_pressure_authority =
          committed_pressure_authority_available
              ? &committed_pressure_wall_gradients
              : &reconstructed_committed_pressure_authority;
      std::vector<finite_volume::detail::ImmersedWallNormalGradient>
          reconstructed_history_pressure_authority;
      if (stencil.order == MomentumTimeOrder::bdf2 &&
          !history_pressure_authority_available)
        reconstructed_history_pressure_authority =
            reconstructed_wall_pressure_gradients(state, history);
      const auto *history_pressure_authority =
          history_pressure_authority_available
              ? &history_pressure_wall_gradients
              : &reconstructed_history_pressure_authority;
      const auto residual_n = assemble_immersed_momentum_residual(
          state, committed, physics.dynamic_viscosity_pa_s,
          committed_pressure_authority);
      const auto implicit_reference_n =
          apply_implicit_viscous_reference(state, committed);
      const auto pressure_n = assemble_immersed_pressure_residual(
          state, committed, committed_pressure_authority);
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
      if (committed_final_momentum_pressure_residual_step ==
          state.metadata().step) {
        if (pressure_n.size() !=
            committed_final_momentum_pressure_residual.size())
          throw runtime::Error(
              "immersed-flow inter-step pressure authority layout is invalid");
        double difference_square = 0.0;
        for (std::size_t offset = 0U; offset < pressure_n.size(); ++offset) {
          const double difference =
              pressure_n[offset] -
              committed_final_momentum_pressure_residual[offset];
          difference_square += difference * difference;
        }
        mpi->allreduce_fp64_in_place(&difference_square, 1U,
                                     runtime::Fp64ReductionOperation::sum);
        pending_inter_step_authority_difference_l2 =
            std::sqrt(difference_square);
      }
#endif
      std::vector<double> residual_nm1;
      std::vector<double> implicit_reference_nm1;
      std::vector<double> pressure_nm1;
      if (stencil.order == MomentumTimeOrder::bdf2) {
        residual_nm1 = assemble_immersed_momentum_residual(
            state, history, physics.dynamic_viscosity_pa_s,
            history_pressure_authority);
        implicit_reference_nm1 =
            apply_implicit_viscous_reference(state, history);
        pressure_nm1 = assemble_immersed_pressure_residual(
            state, history, history_pressure_authority);
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
        if (history_final_momentum_pressure_residual_step !=
                std::numeric_limits<std::uint64_t>::max() &&
            history_final_momentum_pressure_residual_step + 1U ==
                state.metadata().step) {
          if (pressure_nm1.size() !=
              history_final_momentum_pressure_residual.size())
            throw runtime::Error(
                "immersed-flow BDF2 history pressure authority layout is invalid");
          double difference_square = 0.0;
          for (std::size_t offset = 0U; offset < pressure_nm1.size();
               ++offset) {
            const double difference =
                pressure_nm1[offset] -
                history_final_momentum_pressure_residual[offset];
            difference_square += difference * difference;
          }
          mpi->allreduce_fp64_in_place(&difference_square, 1U,
                                       runtime::Fp64ReductionOperation::sum);
          pending_bdf2_history_authority_difference_l2 =
              std::sqrt(difference_square);
        }
#endif
      }
      const auto velocity_n = committed.acquire_read<double>(
          access, kStatePhase, kStateActor, fields.velocity);
      const auto velocity_nm1 = history.acquire_read<double>(
          access, kStatePhase, kStateActor, fields.velocity);
      for (std::size_t component_index = 0U; component_index < 3U;
           ++component_index) {
        auto rhs = momentum_rhs[component_index]->view(0U, count);
        auto solution = momentum_solution[component_index]->view(0U, count);
        for (std::size_t row = 0U; row < count; ++row) {
          const auto cell = owned_active_cells[row];
          const auto index = field_index(cell);
          const double volume = geometry->cell_volume_m3(cell);
          const auto offset = row * 3U + component_index;
          const double predictor_spatial_residual =
              stencil.order == MomentumTimeOrder::backward_euler
                  ? residual_n[offset] - implicit_reference_n[offset]
                  : 2.0 * (residual_n[offset] - pressure_n[offset] -
                           implicit_reference_n[offset]) -
                        (residual_nm1[offset] - pressure_nm1[offset] -
                         implicit_reference_nm1[offset]) +
                        pressure_n[offset];
          const double rho_n =
              material_density_committed.has_value()
                  ? material_density_committed->owned_active_density[row]
                  : physics.rho_ref_kg_per_m3;
          const double rho_nm1 =
              material_density_history.has_value()
                  ? material_density_history->owned_active_density[row]
                  : physics.rho_ref_kg_per_m3;
          rhs[row] =
              -(volume / stencil.dt_s) *
                  (stencil.alpha1 * rho_n *
                       velocity_n(index.i, index.j, index.k,
                                  static_cast<int>(component_index)) +
                   stencil.alpha2 * rho_nm1 *
                       velocity_nm1(index.i, index.j, index.k,
                                    static_cast<int>(component_index))) -
                     predictor_spatial_residual;
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
          if (!manufactured_body_source.empty())
            rhs[row] += volume * manufactured_body_source[offset];
#endif
          solution[row] = velocity_n(index.i, index.j, index.k,
                                     static_cast<int>(component_index));
        }
        momentum_preconditioners[component_index]->update(
            *momentum_operators[component_index],
            momentum_operators[component_index]->revision());
        const auto momentum_rhs_const =
            static_cast<const execution::Buffer &>(
                *momentum_rhs[component_index])
                .view(0U, count);
        auto strict_momentum_control = momentum_control;
        strict_momentum_control.atol =
            std::min(strict_momentum_control.atol, 1.0e-15);
        strict_momentum_control.rtol =
            std::min(strict_momentum_control.rtol, 1.0e-14);
        strict_momentum_control.max_iterations =
            std::max<std::uint64_t>(strict_momentum_control.max_iterations,
                                    500U);
        const auto predictor_control = relative_equation_control(
            strict_momentum_control, momentum_rhs_const);
        report.momentum.components[component_index] =
            momentum_solver->solve(*momentum_operators[component_index],
                                   *momentum_preconditioners[component_index],
                                   momentum_rhs_const, solution,
                                   predictor_control);
        const auto momentum_status = runtime::collective_status(
            *mpi,
            solve_succeeded(report.momentum.components[component_index].reason),
            "immersed-flow momentum solve did not converge");
        if (!momentum_status.ok) {
          auto &solve = report.momentum.components[component_index];
          solve.lowest_failing_rank = momentum_status.failing_rank;
          last_linear_solve_failure = ImmersedLinearSolveFailure{
              ImmersedLinearSolvePhase::predictor_momentum,
              0U,
              static_cast<std::uint32_t>(component_index),
              solve,
              0.0,
              0.0,
              0.0};
          throw std::pair{momentum_status.failing_rank,
                          StepFailureReason::momentum_linear_solve};
        }
      }
      auto velocity = trial.acquire_write<double>(access, kStatePhase,
                                                  kStateActor, fields.velocity);
      for (std::size_t row = 0U; row < count; ++row) {
        const auto index = field_index(owned_active_cells[row]);
        for (std::size_t component_index = 0U; component_index < 3U;
             ++component_index)
          velocity(index.i, index.j, index.k,
                   static_cast<int>(component_index)) =
              momentum_solution[component_index]->view(0U, count)[row];
      }
      halo->exchange(trial, fields.velocity);
      const auto prepare_pressure_density = [&]() {
        if (material_density_trial.has_value())
          pressure_operator->prepare_face_coefficients(
              material_density_trial->face_density, diagonal);
        else
          pressure_operator->prepare_face_coefficients(
              physics.rho_ref_kg_per_m3, diagonal);
      };
      prepare_pressure_density();

      std::vector<detail::ImmersedWallPressureCondition> walls;
      const auto assemble_walls_collectively = [&](const auto *current) {
        bool local_ok = true;
        try {
          walls = assemble_face_predictor(state, physics.rho_ref_kg_per_m3,
                                          stencil, diagonal, current,
                                          material_density_trial
                                              ? &*material_density_trial
                                              : nullptr,
                                          material_density_committed
                                              ? &*material_density_committed
                                              : nullptr,
                                          material_density_history
                                              ? &*material_density_history
                                              : nullptr);
        } catch (...) {
          local_ok = false;
        }
        const auto status = runtime::collective_status(
            *mpi, local_ok, "immersed-flow wall pressure input is invalid");
        if (!status.ok) {
          throw std::pair{status.failing_rank,
                          StepFailureReason::non_finite_trial};
        }
      };
      assemble_walls_collectively(committed_pressure_authority);
      const auto first =
          correct_active_pressure(state, physics.rho_ref_kg_per_m3, stencil,
                                  diagonal, walls, pressure_control,
                                  momentum_control, 1U,
                                  material_density_trial
                                      ? &*material_density_trial
                                      : nullptr,
                                  material_density_committed
                                      ? &*material_density_committed
                                      : nullptr,
                                  material_density_history
                                      ? &*material_density_history
                                      : nullptr);
      report.pressure[0] = first.solve;
      if (!first.accepted) {
        state.rollback_attempt();
        attempt_active = false;
        report.disposition =
            first.disposition ==
                    PressureCorrectionDisposition::recoverable_failure
                ? StepAttemptDisposition::recoverable_failure
                : StepAttemptDisposition::non_retryable_failure;
        report.reason = first.reason;
        report.lowest_failing_rank = first.lowest_failing_rank;
        return report;
      }
      report.pressure_corrector_count = 1U;
      const auto failure_one =
          injected_status(AttemptFailureStage::after_corrector_1);
      if (!failure_one.ok)
        throw std::pair{failure_one.failing_rank,
                        StepFailureReason::non_finite_trial};

      if (variable_density_adapter) {
        const auto provisional =
            material_adapter != nullptr
                ? material_adapter->stage_after_corrector_one(state, stencil)
                : ideal_adapter->stage_after_corrector_one(state, stencil);
        if (provisional.reason != MaterialTransportFailureReason::none) {
          if (material_failure != nullptr)
            *material_failure = provisional.reason;
          throw std::pair{provisional.lowest_failing_rank,
                          material_step_failure(provisional.reason)};
        }
        if (ideal_adapter != nullptr) {
          const auto closure_stage = ideal_adapter->evaluate(
              state, detail::DensityClosureStage::provisional);
          if (accepted_closure != nullptr)
            *accepted_closure = ideal_adapter->take_closure_report();
          if (!closure_stage.accepted) {
            state.rollback_attempt();
            attempt_active = false;
            report.reason = StepFailureReason::density_closure_failure;
            report.lowest_failing_rank =
                closure_stage.lowest_failing_rank;
            report.disposition =
                closure_stage.recoverable
                    ? StepAttemptDisposition::recoverable_failure
                    : StepAttemptDisposition::non_retryable_failure;
            report.suggested_dt_s = closure_stage.recoverable
                                        ? 0.5 * stencil.dt_s
                                        : stencil.dt_s;
            return report;
          }
        }
        halo->exchange(trial, fields.density);
        for (const auto field : fields.transported_cell_fields)
          halo->exchange(trial, field);
        if (ideal_adapter != nullptr)
          ideal_adapter->after_halo(
              detail::DensityClosureStage::provisional);
        require_density_authority([&] {
          material_density_trial.emplace(density_authority(state, trial));
        });
        for (std::size_t component_index = 0U; component_index < 3U;
             ++component_index) {
          for (std::size_t row = 0U; row < count; ++row) {
            const double volume =
                geometry->cell_volume_m3(owned_active_cells[row]);
            diagonal[component_index][row] =
                stencil.alpha0 *
                material_density_trial->owned_active_density[row] * volume /
                stencil.dt_s;
          }
          if (active_wale_authority == nullptr) {
            momentum_operators[component_index]->replace(
                diagonal[component_index], physics.dynamic_viscosity_pa_s);
          } else {
            momentum_operators[component_index]->replace(
                diagonal[component_index],
                active_wale_authority
                    ->effective_dynamic_viscosity_by_active_cell(),
                active_wale_authority
                    ->effective_dynamic_viscosity_by_face());
          }
          diagonal[component_index] =
              momentum_operators[component_index]->assembled_diagonal();
        }
        prepare_pressure_density();
      }

      const auto first_corrected_wall_gradients =
          wall_pressure_gradients(walls, WallPressureGradientState::corrected);
      assemble_walls_collectively(&first_corrected_wall_gradients);
      const auto second =
          correct_active_pressure(state, physics.rho_ref_kg_per_m3, stencil,
                                  diagonal, walls, pressure_control,
                                  momentum_control, 2U,
                                  material_density_trial
                                      ? &*material_density_trial
                                      : nullptr,
                                  material_density_committed
                                      ? &*material_density_committed
                                      : nullptr,
                                  material_density_history
                                      ? &*material_density_history
                                      : nullptr);
      report.pressure[1] = second.solve;
      if (!second.accepted) {
        state.rollback_attempt();
        attempt_active = false;
        report.disposition =
            second.disposition ==
                    PressureCorrectionDisposition::recoverable_failure
                ? StepAttemptDisposition::recoverable_failure
                : StepAttemptDisposition::non_retryable_failure;
        report.reason = second.reason;
        report.lowest_failing_rank = second.lowest_failing_rank;
        return report;
      }
      report.pressure_corrector_count = 2U;
      last_corrector_count = 2U;
      const auto failure_two =
          injected_status(AttemptFailureStage::after_corrector_2);
      if (!failure_two.ok)
        throw std::pair{failure_two.failing_rank,
                        StepFailureReason::non_finite_trial};
      if (variable_density_adapter) {
        const auto final =
            material_adapter != nullptr
                ? material_adapter->finalize_from_corrector_two_flux(state,
                                                                      stencil)
                : ideal_adapter->finalize_from_corrector_two_flux(state,
                                                                   stencil);
        if (material_failure != nullptr)
          *material_failure = final.reason;
        if (final.reason != MaterialTransportFailureReason::none)
          throw std::pair{final.lowest_failing_rank,
                          material_step_failure(final.reason)};
        const auto final_report =
            material_adapter != nullptr
                ? material_adapter->take_final_report()
                : ideal_adapter->take_final_report();
        if (!final_report.has_value())
          throw runtime::Error(
              "immersed material final report is unavailable");
        report.final_transport_normalized_l2 =
            final_report->transport_normalized_l2();
        report.final_transport_relative_conservation_defect =
            final_report->transport_relative_conservation_defect();
        report.final_mass_relative_conservation_defect =
            final_report->mass_relative_conservation_defect();
        if (accepted_material != nullptr)
          *accepted_material = *final_report;
        if (ideal_adapter != nullptr)
          accepted_material_guard.output = nullptr;
        if (ideal_adapter != nullptr) {
          const auto closure_stage = ideal_adapter->evaluate(
              state, detail::DensityClosureStage::final);
          if (accepted_closure != nullptr)
            *accepted_closure = ideal_adapter->take_closure_report();
          if (!closure_stage.accepted) {
            state.rollback_attempt();
            attempt_active = false;
            report.reason = StepFailureReason::density_closure_failure;
            report.lowest_failing_rank =
                closure_stage.lowest_failing_rank;
            report.disposition =
                closure_stage.recoverable
                    ? StepAttemptDisposition::recoverable_failure
                    : StepAttemptDisposition::non_retryable_failure;
            report.suggested_dt_s = closure_stage.recoverable
                                        ? 0.5 * stencil.dt_s
                                        : stencil.dt_s;
            return report;
          }
        }
        halo->exchange(trial, fields.density);
        for (const auto field : fields.transported_cell_fields)
          halo->exchange(trial, field);
        if (ideal_adapter != nullptr) {
          ideal_adapter->after_halo(detail::DensityClosureStage::final);
          const auto post =
              ideal_adapter->assess_after_closure(state, stencil);
          const auto post_report =
              ideal_adapter->take_post_closure_report();
          if (!post_report)
            throw runtime::Error(
                "immersed ideal-gas post-closure report is unavailable");
          if (accepted_post_closure_material != nullptr)
            *accepted_post_closure_material = *post_report;
          report.final_transport_normalized_l2 =
              post_report->transport_normalized_l2();
          report.final_transport_relative_conservation_defect =
              post_report->transport_relative_conservation_defect();
          report.final_mass_relative_conservation_defect =
              post_report->mass_relative_conservation_defect();
          if (post.reason != MaterialTransportFailureReason::none)
            throw std::pair{post.lowest_failing_rank,
                            material_step_failure(post.reason)};
        }
      }
      const auto implicit_reference_trial =
          apply_implicit_viscous_reference(state, trial);
      const auto final_wall_gradients =
          wall_pressure_gradients(walls, WallPressureGradientState::corrected);
      if (!assess_final_momentum(state, physics.rho_ref_kg_per_m3, stencil,
                                 residual_n, residual_nm1, implicit_reference_n,
                                 implicit_reference_nm1,
                                 implicit_reference_trial, pressure_n,
                                 pressure_nm1, final_wall_gradients, report))
        throw std::pair{0, StepFailureReason::final_momentum_residual};
      const auto forced_momentum =
          injected_status(AttemptFailureStage::final_momentum_residual);
      if (!forced_momentum.ok)
        throw std::pair{forced_momentum.failing_rank,
                        StepFailureReason::final_momentum_residual};
      const auto failure_transport =
          injected_status(AttemptFailureStage::after_final_transport);
      if (!failure_transport.ok) {
        if (ideal_adapter != nullptr) {
          state.rollback_attempt();
          attempt_active = false;
          report.disposition =
              StepAttemptDisposition::non_retryable_failure;
          report.reason = StepFailureReason::invalid_input;
          report.lowest_failing_rank = failure_transport.failing_rank;
          report.suggested_dt_s = stencil.dt_s;
          return report;
        }
        if (material_adapter != nullptr && material_failure != nullptr)
          *material_failure = MaterialTransportFailureReason::non_finite_state;
        throw std::pair{failure_transport.failing_rank,
                        StepFailureReason::transport_failure};
      }

      bool final_wall_ok = assess_final_wall_penetration(state, walls);
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
      const int configured_failure_rank =
          injected_failure_rank.load(std::memory_order_relaxed);
      const int selected_failure_rank = configured_failure_rank >= 0
                                            ? configured_failure_rank
                                            : (mpi->size() > 1 ? 1 : 0);
      if (injected_failure_stage.load(std::memory_order_relaxed) ==
              AttemptFailureStage::final_wall_penetration &&
          mpi->rank() == selected_failure_rank) {
        final_wall_ok = false;
      }
#endif
      const auto final_wall_status = runtime::collective_status(
          *mpi, final_wall_ok,
          "immersed-flow final wall penetration residual is invalid");
      if (!final_wall_status.ok)
        throw std::pair{final_wall_status.failing_rank,
                        StepFailureReason::final_continuity_residual};

      const auto final_flux = trial.acquire_face_read<double>(
          access, kStatePhase, kStateActor, fields.face_mass_flux);
      std::fill(attempt_continuity_divergence.begin(),
                attempt_continuity_divergence.end(), 0.0);
      std::vector<double> continuity_absolute;
      if (variable_density_adapter)
        continuity_absolute.assign(count, 0.0);
      for (mesh::LocalFaceId face = 0U; face < topology->local_face_count();
           ++face) {
        const double value = final_flux(face, 0);
        const auto owner_index =
            domain->active_cells().active_index(topology->owner(face));
        if (owner_index.has_value() && *owner_index < count) {
          attempt_continuity_divergence[*owner_index] += value;
          if (variable_density_adapter)
            continuity_absolute[*owner_index] += std::abs(value);
        }
        const auto neighbour = topology->neighbour(face);
        if (!neighbour.has_value() || topology->periodic_pair(face).has_value())
          continue;
        const auto neighbour_index =
            domain->active_cells().active_index(*neighbour);
        if (neighbour_index.has_value() && *neighbour_index < count) {
          attempt_continuity_divergence[*neighbour_index] -= value;
          if (variable_density_adapter)
            continuity_absolute[*neighbour_index] += std::abs(value);
        }
      }
      double continuity_parts[2]{};
      if (variable_density_adapter) {
        const auto rho_final = trial.acquire_read<double>(
            access, kStatePhase, kStateActor, fields.density);
        const auto rho_n = committed.acquire_read<double>(
            access, kStatePhase, kStateActor, fields.density);
        const auto rho_nm1 = history.acquire_read<double>(
            access, kStatePhase, kStateActor, fields.density);
        for (std::size_t row = 0U; row < count; ++row) {
          const auto cell = owned_active_cells[row];
          const auto index = field_index(cell);
          const double volume = geometry->cell_volume_m3(cell);
          const double temporal =
              (stencil.alpha0 * rho_final(index.i, index.j, index.k, 0) +
               stencil.alpha1 * rho_n(index.i, index.j, index.k, 0) +
               stencil.alpha2 * rho_nm1(index.i, index.j, index.k, 0)) *
              volume / stencil.dt_s;
          const double raw =
              temporal + attempt_continuity_divergence[row];
          const double scale =
              (std::abs(stencil.alpha0 *
                        rho_final(index.i, index.j, index.k, 0)) +
               std::abs(stencil.alpha1 *
                        rho_n(index.i, index.j, index.k, 0)) +
               std::abs(stencil.alpha2 *
                        rho_nm1(index.i, index.j, index.k, 0))) *
                  volume / stencil.dt_s +
              continuity_absolute[row];
          if (!std::isfinite(raw) || !std::isfinite(scale))
            throw runtime::Error(
                "immersed material final continuity authority is invalid");
          continuity_parts[0] += raw * raw;
          continuity_parts[1] += scale * scale;
        }
      } else {
        for (std::size_t row = 0U; row < count; ++row) {
          const auto cell = owned_active_cells[row];
          const double volume = geometry->cell_volume_m3(cell);
          const double normalized =
              attempt_continuity_divergence[row] / volume;
          continuity_parts[0] += normalized * normalized * volume;
          continuity_parts[1] += volume;
        }
      }
      mpi->allreduce_fp64_in_place(continuity_parts, 2U,
                                   runtime::Fp64ReductionOperation::sum);
      if (variable_density_adapter) {
        const double numerator = std::sqrt(continuity_parts[0]);
        const double denominator = std::sqrt(continuity_parts[1]);
        report.final_continuity_normalized_l2 =
            denominator == 0.0
                ? (numerator == 0.0
                       ? 0.0
                       : std::numeric_limits<double>::infinity())
                : numerator / denominator;
      } else {
        report.final_continuity_normalized_l2 =
            std::sqrt(continuity_parts[0] / continuity_parts[1]);
      }
      report.final_pressure_residual_l2 = second.independent_residual_l2;
      const auto forced_continuity =
          injected_status(AttemptFailureStage::final_continuity_residual);
      if (!forced_continuity.ok)
        throw std::pair{forced_continuity.failing_rank,
                        StepFailureReason::final_continuity_residual};
      const auto forced_pressure =
          injected_status(AttemptFailureStage::final_pressure_residual);
      if (!forced_pressure.ok)
        throw std::pair{forced_pressure.failing_rank,
                        StepFailureReason::final_pressure_residual};
      const double final_pressure_limit = std::max(
          pressure_control.atol, pressure_control.rtol * second.rhs_l2);
      if (!std::isfinite(report.final_continuity_normalized_l2) ||
          report.final_continuity_normalized_l2 > 1.0e-10 ||
          !std::isfinite(report.final_pressure_residual_l2)) {
        throw std::pair{0, StepFailureReason::final_continuity_residual};
      }
      if (report.final_pressure_residual_l2 > final_pressure_limit)
        throw std::pair{0, StepFailureReason::final_pressure_residual};
      const auto metadata = state.metadata();
      if (ideal_adapter != nullptr && boundaries->open_domain())
        ideal_adapter->before_outlet(state);
      const auto backflow = boundaries->assess_final_pressure_outlet_flux(
          *topology, *mpi, final_flux, metadata.step + 1U,
          metadata.time_s + stencil.dt_s);
      report.final_backflow_evidence = backflow.evidence;
      if (backflow.decision == boundary::FinalFluxDecision::outlet_backflow) {
        const int failing_rank = backflow.evidence.has_value()
                                     ? backflow.evidence->lowest_failing_rank
                                     : -1;
        throw std::pair{failing_rank, StepFailureReason::boundary_backflow};
      }
      if (wall_force_integrator.has_value()) {
        auto force = collect_final_force(state, physics.dynamic_viscosity_pa_s,
                                         final_wall_gradients);
        if (!force.report.has_value() || !force.sample.has_value())
          throw std::pair{force.failing_rank,
                          StepFailureReason::final_conservation_defect};
        const auto forced_force =
            injected_status(AttemptFailureStage::final_force_reconstruction);
        if (!forced_force.ok)
          throw std::pair{forced_force.failing_rank,
                          StepFailureReason::final_conservation_defect};
        pending_diagnostic_wall_force_sample = *force.sample;
        pending_diagnostic_local_flow_pattern = immersed_operator->report();
        if (accepted_force != nullptr)
          *accepted_force = std::move(*force.report);
      }
      bool local_commit_preparation_ok = true;
      try {
        if (final_wall_gradients.size() !=
            pending_pressure_wall_gradients.size())
          throw runtime::Error(
              "immersed-flow final pressure authority layout is invalid");
        for (std::size_t index = 0U; index < final_wall_gradients.size();
             ++index) {
          if (final_wall_gradients[index].link !=
              pending_pressure_wall_gradients[index].link)
            throw runtime::Error(
                "immersed-flow final pressure authority ordering is invalid");
        }
        if (active_wale_authority != nullptr &&
            pending_wale_evaluation_count !=
                active_wale_authority->evaluation_count())
          throw runtime::Error(
              "immersed-flow WALE evaluation count lost attempt authority");
        const AcceptedStepMetadata accepted_metadata{
            metadata.step + 1U, metadata.time_s + stencil.dt_s, stencil.dt_s,
            metadata.dt_s, stencil.order};
        if (material_adapter != nullptr) {
          material_adapter->prepare_commit();
          state.prepare_commit_attempt(accepted_metadata);
        } else if (ideal_adapter != nullptr) {
          const int before_prepare_rank =
              ideal_adapter->before_prepare(state, accepted_metadata);
          if (before_prepare_rank >= 0 &&
              mpi->rank() == before_prepare_rank)
            throw runtime::Error(
                "immersed ideal-gas pre-commit hook failed");
          state.prepare_commit_attempt(accepted_metadata);
          const int closure_prepare_rank = ideal_adapter->prepare_commit();
          if (closure_prepare_rank >= 0 &&
              mpi->rank() == closure_prepare_rank)
            throw runtime::Error(
                "immersed ideal-gas closure preparation failed");
        } else {
          state.prepare_commit_attempt(accepted_metadata);
        }
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
        if (injected_failure_stage.load(std::memory_order_relaxed) ==
            AttemptFailureStage::before_commit) {
          const int configured_commit_failure_rank =
              injected_failure_rank.load(std::memory_order_relaxed);
          const int selected_commit_failure_rank =
              configured_commit_failure_rank >= 0
                  ? configured_commit_failure_rank
                  : (mpi->size() > 1 ? 1 : 0);
          if (mpi->rank() == selected_commit_failure_rank)
            throw runtime::Error(
                "injected immersed-flow commit preparation failure");
        }
#endif
      } catch (...) {
        local_commit_preparation_ok = false;
      }
      const auto commit_preparation_status = runtime::collective_status(
          *mpi, local_commit_preparation_ok,
          "immersed-flow collective commit preparation failed");
      if (!commit_preparation_status.ok) {
        state.rollback_attempt();
        attempt_active = false;
        report.final_continuity_normalized_l2 = 0.0;
        report.final_pressure_residual_l2 = 0.0;
        if (material_adapter != nullptr) {
          report.final_momentum_normalized_l2.fill(0.0);
          report.final_momentum_relative_conservation_defect.fill(0.0);
        }
        report.disposition = StepAttemptDisposition::non_retryable_failure;
        report.reason = StepFailureReason::invalid_input;
        report.lowest_failing_rank =
            commit_preparation_status.failing_rank;
        report.suggested_dt_s = stencil.dt_s;
        return report;
      }
      for (std::size_t index = 0U; index < final_wall_gradients.size();
           ++index)
        pending_pressure_wall_gradients[index].value =
            final_wall_gradients[index].value;
      // The collective-ready boundary is above.  Everything below is a
      // no-throw publication of already prepared storage or scalar authority.
      state.publish_commit_attempt();
      if (material_adapter != nullptr) {
        material_adapter->publish();
        density_adapter_guard.material = nullptr;
        accepted_material_guard.output = nullptr;
      } else if (ideal_adapter != nullptr) {
        ideal_adapter->publish();
        density_adapter_guard.ideal = nullptr;
        accepted_material_guard.output = nullptr;
      }
      static_assert(noexcept(history_pressure_wall_gradients.swap(
          committed_pressure_wall_gradients)));
      static_assert(noexcept(committed_pressure_wall_gradients.swap(
          pending_pressure_wall_gradients)));
      history_pressure_wall_gradients.swap(committed_pressure_wall_gradients);
      committed_pressure_wall_gradients.swap(pending_pressure_wall_gradients);
      history_pressure_authority_available =
          committed_pressure_authority_available;
      committed_pressure_authority_available = true;
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
      static_assert(noexcept(committed_final_momentum_pressure_residual.swap(
          pending_final_momentum_pressure_residual)));
      history_final_momentum_pressure_residual.swap(
          committed_final_momentum_pressure_residual);
      committed_final_momentum_pressure_residual.swap(
          pending_final_momentum_pressure_residual);
      history_final_pressure_wall_gradients.swap(
          committed_final_pressure_wall_gradients);
      committed_final_pressure_wall_gradients.swap(
          pending_final_pressure_wall_gradients);
      history_final_momentum_pressure_residual_step =
          committed_final_momentum_pressure_residual_step;
      committed_final_momentum_pressure_residual_step = metadata.step + 1U;
#endif
      committed_exact_predictor_work = pending_exact_predictor_work;
      if (active_wale_authority != nullptr) {
        last_wale_evaluation_count = pending_wale_evaluation_count;
        last_wale_coefficient_identity =
            active_wale_authority->summary().identity;
        last_wall_effective_viscosity_fingerprint =
            pending_wall_effective_viscosity_fingerprint == 0U
                ? active_wale_authority
                      ->wall_effective_viscosity_fingerprint()
                : pending_wall_effective_viscosity_fingerprint;
        if (accepted_wale != nullptr)
          *accepted_wale = active_wale_authority->summary();
      }
      static_assert(std::is_nothrow_copy_assignable_v<
                    std::optional<immersed::WallForceSample>>);
      committed_diagnostic_wall_force_sample =
          pending_diagnostic_wall_force_sample;
      static_assert(std::is_nothrow_copy_assignable_v<std::optional<
                    finite_volume::ImmersedOperatorReport>>);
      committed_diagnostic_local_flow_pattern =
          pending_diagnostic_local_flow_pattern;
      attempt_active = false;
      report.disposition = StepAttemptDisposition::committed;
      report.reason = StepFailureReason::none;
      report.lowest_failing_rank = -1;
      report.suggested_dt_s = stencil.dt_s;
      return report;
    } catch (const std::pair<int, StepFailureReason> &failure) {
      if (attempt_active)
        state.rollback_attempt();
      report.disposition = StepAttemptDisposition::recoverable_failure;
      report.reason = failure.second;
      report.lowest_failing_rank = failure.first;
      report.suggested_dt_s = 0.5 * stencil.dt_s;
      return report;
    } catch (const runtime::detail::MpiOperationError &) {
      if (attempt_active)
        state.rollback_attempt();
      if (variable_density_adapter)
        throw;
      report.disposition = StepAttemptDisposition::non_retryable_failure;
      report.reason = StepFailureReason::collective_operation;
      report.lowest_failing_rank = -1;
      return report;
    } catch (const std::exception &error) {
      static_cast<void>(error);
      if (attempt_active)
        state.rollback_attempt();
      if (variable_density_adapter)
        throw;
      report.disposition = StepAttemptDisposition::non_retryable_failure;
      report.reason = StepFailureReason::invalid_input;
      report.lowest_failing_rank = -1;
      return report;
    } catch (...) {
      if (attempt_active)
        state.rollback_attempt();
      if (variable_density_adapter)
        throw;
      report.disposition = StepAttemptDisposition::non_retryable_failure;
      report.reason = StepFailureReason::invalid_input;
      report.lowest_failing_rank = -1;
      return report;
    }
  }

  void synchronize_pressure_authority_binding(const FlowState &state) noexcept {
    const auto metadata = state.metadata();
    const bool matches =
        pressure_authority_state_instance ==
            detail::FlowStateSolverAccess::instance_identity(state) &&
        pressure_authority_state_mutation ==
            state.diagnostic_mutation_identity() &&
        pressure_authority_state_step == metadata.step;
    if (!matches) {
      history_pressure_authority_available = false;
      committed_pressure_authority_available = false;
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
      history_final_momentum_pressure_residual.clear();
      committed_final_momentum_pressure_residual.clear();
      pending_final_momentum_pressure_residual.clear();
      history_final_pressure_wall_gradients.clear();
      committed_final_pressure_wall_gradients.clear();
      pending_final_pressure_wall_gradients.clear();
      history_final_momentum_pressure_residual_step =
          std::numeric_limits<std::uint64_t>::max();
      committed_final_momentum_pressure_residual_step =
          std::numeric_limits<std::uint64_t>::max();
#endif
    }
    pressure_authority_state_instance =
        detail::FlowStateSolverAccess::instance_identity(state);
    pressure_authority_state_mutation = state.diagnostic_mutation_identity();
    pressure_authority_state_step = metadata.step;
  }

  void refresh_pressure_authority_binding(const FlowState &state) noexcept {
    pressure_authority_state_instance =
        detail::FlowStateSolverAccess::instance_identity(state);
    pressure_authority_state_mutation = state.diagnostic_mutation_identity();
    pressure_authority_state_step = state.metadata().step;
  }

  const runtime::MpiContext *mpi{};
  const immersed::ImmersedDomain *domain{};
  const immersed::GhostStencilPlan *ghost_plan{};
  const immersed::WallQuadraturePlan *wall_quadrature{};
  const immersed::LocalFlowPatternTransform *transform{};
  const les::WaleModel *wale{};
  ImmersedFlowDensitySetup density_setup;
  detail::ImmersedDensityAdapter density_adapter;
  const runtime::StructuredDecomposition *decomposition{};
  const mesh::MeshTopology *topology{};
  const mesh::MeshGeometry *geometry{};
  const boundary::BoundaryRegistry *boundaries{};
  runtime::HaloExchange *halo{};
  execution::ExecutionContext *execution{};
  const linear::LinearSolver *momentum_solver{};
  const linear::LinearSolver *pressure_solver{};
  std::array<linear::Preconditioner *, 3> momentum_preconditioners{};
  linear::Preconditioner *pressure_preconditioner{};
  std::uint64_t geometry_fingerprint{};
  mutable diagnostics::Stage3PerformanceCounters performance;
  std::optional<FixedStepConstantDensityFlow> body_fitted;
  std::optional<FixedStepMaterialDensityFlow> body_fitted_material;
  std::optional<detail::DensityClosureHooks> body_fitted_closure_hooks;
  std::optional<finite_volume::CellCenteredFvmOperators> physical_boundary_fvm;
  std::optional<finite_volume::ImmersedReconstruction> reconstruction;
  std::optional<immersed::WallForceIntegrator> wall_force_integrator;
  std::optional<finite_volume::ImmersedOperatorAdapter> immersed_operator;
  std::unique_ptr<ImmersedFlowScratch> scratch;
  std::array<std::unique_ptr<ActiveMomentumOperator>, 3> momentum_operators;
  std::array<std::unique_ptr<execution::Buffer>, 3> momentum_rhs;
  std::array<std::unique_ptr<execution::Buffer>, 3> momentum_solution;
  std::unique_ptr<ActivePressureOperator> pressure_operator;
  std::unique_ptr<ExactPredictorSchurOperator> exact_pressure_operator;
  std::unique_ptr<ActivePressureTwoLevelPreconditioner>
      compact_pressure_preconditioner;
  std::unique_ptr<linear::RestartedFgmresSolver> compact_pressure_solver;
  std::unique_ptr<CompactSchurPreconditioner> compact_schur_preconditioner;
  std::unique_ptr<execution::Buffer> pressure_rhs;
  std::unique_ptr<execution::Buffer> pressure_solution;
  std::unique_ptr<execution::Buffer> pressure_residual;
  std::optional<FlowState> exact_predictor_probe_state;
  FlowFieldIds exact_predictor_fields;
  ExactPredictorResponse exact_predictor_response_workspace;
  ExactPredictorResponse exact_affine_response_workspace;
  ExactPredictorResponse exact_homogeneous_baseline_workspace;
  ExactPredictorResponse exact_homogeneous_response_workspace;
  ExactPredictorResponse::Work pending_exact_predictor_work;
  ExactPredictorResponse::Work committed_exact_predictor_work;
  std::optional<ImmersedLinearSolveFailure> last_linear_solve_failure;
  std::uint64_t exact_predictor_probe_creation_count{};
  std::uint64_t exact_predictor_response_capacity_growth_count{};
  std::vector<mesh::LocalCellId> owned_active_cells;
  std::vector<WallRuntime> wall_links;
  mutable std::vector<finite_volume::PhysicalBoundaryMomentumContribution>
      physical_boundary_momentum;
  mutable std::vector<finite_volume::PhysicalBoundaryPressureContribution>
      physical_boundary_pressure;
  std::vector<double> attempt_continuity_divergence;
  std::uint64_t pressure_revision{};
  std::uint64_t pressure_dependency_identity{};
  std::uint32_t last_corrector_count{};
  const detail::ImmersedWaleAttemptAuthority *active_wale_authority{};
  std::uint64_t pending_wale_evaluation_count{};
  mutable std::uint64_t pending_wall_effective_viscosity_fingerprint{};
  std::uint64_t last_wale_evaluation_count{};
  les::WaleCoefficientIdentity last_wale_coefficient_identity{};
  std::uint64_t last_wall_effective_viscosity_fingerprint{};
  const FlowState *last_diagnostic_state{};
  std::uint64_t last_diagnostic_state_identity{};
  std::uint64_t last_diagnostic_report_seal{};
  std::uint64_t diagnostic_attempt_generation{};
  double last_diagnostic_rho_ref_kg_per_m3{};
  config::DensityModel last_diagnostic_density_model{
      config::DensityModel::constant};
  std::optional<immersed::WallForceSample>
      pending_diagnostic_wall_force_sample;
  std::optional<immersed::WallForceSample>
      committed_diagnostic_wall_force_sample;
  std::optional<finite_volume::ImmersedOperatorReport>
      pending_diagnostic_local_flow_pattern;
  std::optional<finite_volume::ImmersedOperatorReport>
      committed_diagnostic_local_flow_pattern;
  std::vector<finite_volume::detail::ImmersedWallNormalGradient>
      history_pressure_wall_gradients;
  std::vector<finite_volume::detail::ImmersedWallNormalGradient>
      committed_pressure_wall_gradients;
  std::vector<finite_volume::detail::ImmersedWallNormalGradient>
      pending_pressure_wall_gradients;
  std::uint64_t pressure_authority_state_instance{};
  std::uint64_t pressure_authority_state_mutation{};
  std::uint64_t pressure_authority_state_step{
      std::numeric_limits<std::uint64_t>::max()};
  bool history_pressure_authority_available{};
  bool committed_pressure_authority_available{};
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  std::vector<double> manufactured_body_source;
  std::uint64_t last_exact_predictor_response_allocation_event_count{};
  double last_wall_pressure_gradient_application_norm{};
  std::vector<detail::ImmersedFlowPressureCorrectionRecord>
      pending_cell_pressure_correction_records;
  std::vector<double> pending_previous_lfp_pressure_after_correction;
  double pending_inter_corrector_authority_difference_l2{};
  double pending_inter_step_authority_difference_l2{
      std::numeric_limits<double>::quiet_NaN()};
  double pending_bdf2_history_authority_difference_l2{
      std::numeric_limits<double>::quiet_NaN()};
  mutable std::vector<double> pending_final_momentum_pressure_residual;
  mutable std::vector<finite_volume::detail::ImmersedWallNormalGradient>
      pending_final_pressure_wall_gradients;
  std::vector<double> history_final_momentum_pressure_residual;
  std::vector<finite_volume::detail::ImmersedWallNormalGradient>
      history_final_pressure_wall_gradients;
  std::vector<double> committed_final_momentum_pressure_residual;
  std::vector<finite_volume::detail::ImmersedWallNormalGradient>
      committed_final_pressure_wall_gradients;
  std::uint64_t committed_final_momentum_pressure_residual_step{
      std::numeric_limits<std::uint64_t>::max()};
  std::uint64_t history_final_momentum_pressure_residual_step{
      std::numeric_limits<std::uint64_t>::max()};
#endif
};

FixedStepImmersedFlow FixedStepImmersedFlow::create(
    const runtime::StructuredDecomposition &decomposition,
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    const boundary::BoundaryRegistry &boundaries,
    const immersed::ImmersedDomain *domain,
    const immersed::GhostStencilPlan *ghost_plan,
    const immersed::WallQuadraturePlan *wall_quadrature,
    const immersed::LocalFlowPatternTransform *transform,
    const les::WaleModel *wale, const runtime::MpiContext &mpi,
    execution::ExecutionContext &execution_context, runtime::HaloExchange &halo,
    const linear::LinearSolver &momentum_solver,
    std::array<linear::Preconditioner *, 3> momentum_preconditioners,
    const linear::LinearSolver &pressure_solver,
    linear::Preconditioner &pressure_preconditioner) {
  return create(decomposition, topology, geometry, boundaries, domain,
                ghost_plan, wall_quadrature, transform, wale,
                ImmersedFlowDensitySetup{}, mpi, execution_context, halo,
                momentum_solver, momentum_preconditioners, pressure_solver,
                pressure_preconditioner);
}

FixedStepImmersedFlow FixedStepImmersedFlow::create(
    const runtime::StructuredDecomposition &decomposition,
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    const boundary::BoundaryRegistry &boundaries,
    const immersed::ImmersedDomain *domain,
    const immersed::GhostStencilPlan *ghost_plan,
    const immersed::WallQuadraturePlan *wall_quadrature,
    const immersed::LocalFlowPatternTransform *transform,
    const les::WaleModel *wale, ImmersedFlowDensitySetup density_setup,
    const runtime::MpiContext &mpi,
    execution::ExecutionContext &execution_context, runtime::HaloExchange &halo,
    const linear::LinearSolver &momentum_solver,
    std::array<linear::Preconditioner *, 3> momentum_preconditioners,
    const linear::LinearSolver &pressure_solver,
    linear::Preconditioner &pressure_preconditioner) {
  return FixedStepImmersedFlow(std::make_unique<Impl>(
      decomposition, topology, geometry, boundaries, domain, ghost_plan,
      wall_quadrature, transform, wale, std::move(density_setup), mpi,
      execution_context, halo,
      momentum_solver, momentum_preconditioners, pressure_solver,
      pressure_preconditioner));
}

FixedStepImmersedFlow::FixedStepImmersedFlow(std::unique_ptr<Impl> impl) noexcept
    : impl_(std::move(impl)) {}
FixedStepImmersedFlow::~FixedStepImmersedFlow() noexcept = default;
FixedStepImmersedFlow::FixedStepImmersedFlow(FixedStepImmersedFlow &&) noexcept =
    default;

struct detail::ImmersedFlowCheckpointPreparedRestore::Impl final {
  std::vector<finite_volume::detail::ImmersedWallNormalGradient> history;
  std::vector<finite_volume::detail::ImmersedWallNormalGradient> committed;
  std::vector<finite_volume::detail::ImmersedWallNormalGradient> pending;
  bool history_available{};
  bool committed_available{};
};

detail::ImmersedFlowCheckpointPreparedRestore::
    ImmersedFlowCheckpointPreparedRestore(std::unique_ptr<Impl> impl) noexcept
    : impl_(std::move(impl)) {}
detail::ImmersedFlowCheckpointPreparedRestore::
    ~ImmersedFlowCheckpointPreparedRestore() noexcept = default;
detail::ImmersedFlowCheckpointPreparedRestore::
    ImmersedFlowCheckpointPreparedRestore(
        ImmersedFlowCheckpointPreparedRestore &&) noexcept = default;

detail::CheckpointV3RankAuthority detail::ImmersedFlowCheckpointAccess::snapshot(
    const FixedStepImmersedFlow &flow, std::int32_t rank,
    runtime::Box3 owned_box) {
  if (!flow.impl_ || !flow.impl_->domain) {
    throw runtime::Error("immersed-flow checkpoint source is unavailable");
  }
  const auto encode = [](bool available, const auto &gradients) {
    std::vector<detail::CheckpointV3AuthorityGradient> result;
    if (!available) {
      return result;
    }
    result.reserve(gradients.size());
    for (const auto &gradient : gradients) {
      result.push_back({gradient.link, gradient.value});
    }
    return result;
  };
  return {rank,
          owned_box,
          {},
          {},
          0U,
          0U,
          0U,
          flow.impl_->history_pressure_authority_available,
          flow.impl_->committed_pressure_authority_available,
          encode(flow.impl_->history_pressure_authority_available,
                 flow.impl_->history_pressure_wall_gradients),
          encode(flow.impl_->committed_pressure_authority_available,
                 flow.impl_->committed_pressure_wall_gradients)};
}

detail::ImmersedFlowCheckpointPreparedRestore
detail::ImmersedFlowCheckpointAccess::prepare_restore(
    const FixedStepImmersedFlow &flow,
    const CheckpointV3RankAuthority &authority) {
  if (!flow.impl_ || !flow.impl_->domain) {
    throw runtime::Error("immersed-flow checkpoint target is unavailable");
  }
  const auto &links = flow.impl_->wall_links;
  const auto validate = [&links](bool available, const auto &rows) {
    if (!available) {
      return rows.empty();
    }
    if (rows.size() != links.size()) {
      return false;
    }
    for (std::size_t index = 0U; index < links.size(); ++index) {
      if (rows[index].link != links[index].link.id ||
          !std::isfinite(rows[index].value)) {
        return false;
      }
    }
    return true;
  };
  if (!validate(authority.history_available, authority.history) ||
      !validate(authority.committed_available, authority.committed)) {
    throw runtime::Error(
        "immersed-flow checkpoint authority layout is incompatible");
  }
  auto prepared = std::make_unique<ImmersedFlowCheckpointPreparedRestore::Impl>();
  const auto materialize = [&links](bool available, const auto &rows) {
    std::vector<finite_volume::detail::ImmersedWallNormalGradient> result;
    result.reserve(links.size());
    for (std::size_t index = 0U; index < links.size(); ++index) {
      result.push_back(
          {links[index].link.id, available ? rows[index].value : 0.0});
    }
    return result;
  };
  prepared->history = materialize(authority.history_available, authority.history);
  prepared->committed =
      materialize(authority.committed_available, authority.committed);
  prepared->pending = materialize(false, authority.history);
  prepared->history_available = authority.history_available;
  prepared->committed_available = authority.committed_available;
  return ImmersedFlowCheckpointPreparedRestore(std::move(prepared));
}

void detail::ImmersedFlowCheckpointAccess::publish_restore(
    FixedStepImmersedFlow &flow,
    ImmersedFlowCheckpointPreparedRestore &&prepared,
    const FlowState &state) noexcept {
  if (!flow.impl_ || !prepared.impl_) {
    return;
  }
  static_assert(noexcept(flow.impl_->history_pressure_wall_gradients.swap(
      prepared.impl_->history)));
  flow.impl_->history_pressure_wall_gradients.swap(prepared.impl_->history);
  flow.impl_->committed_pressure_wall_gradients.swap(
      prepared.impl_->committed);
  flow.impl_->pending_pressure_wall_gradients.swap(prepared.impl_->pending);
  flow.impl_->history_pressure_authority_available =
      prepared.impl_->history_available;
  flow.impl_->committed_pressure_authority_available =
      prepared.impl_->committed_available;
  flow.impl_->pressure_revision = 0U;
  flow.impl_->pressure_dependency_identity = 0U;
  flow.impl_->last_diagnostic_state = nullptr;
  flow.impl_->last_diagnostic_state_identity = 0U;
  flow.impl_->last_diagnostic_report_seal = 0U;
  flow.impl_->pending_diagnostic_wall_force_sample.reset();
  flow.impl_->committed_diagnostic_wall_force_sample.reset();
  flow.impl_->pending_diagnostic_local_flow_pattern.reset();
  flow.impl_->committed_diagnostic_local_flow_pattern.reset();
  ++flow.impl_->diagnostic_attempt_generation;
  flow.impl_->refresh_pressure_authority_binding(state);
}

ImmersedFlowStepAttemptReport FixedStepImmersedFlow::attempt(
    FlowState &state, const ImmersedFlowPhysics &physics,
    const MomentumTimeStencil &stencil,
    const linear::SolveControl &momentum_control,
    const linear::SolveControl &pressure_control) const {
  ImmersedFlowStepAttemptReport result{};
  if (!impl_) {
    StepAttemptReport failure{};
    failure.disposition = StepAttemptDisposition::non_retryable_failure;
    failure.reason = StepFailureReason::invalid_input;
    result.base = failure;
    return result;
  }
  if (impl_->domain != nullptr) {
    impl_->last_diagnostic_state = nullptr;
    impl_->last_diagnostic_state_identity = 0U;
    impl_->last_diagnostic_report_seal = 0U;
    impl_->pending_diagnostic_wall_force_sample.reset();
    impl_->committed_diagnostic_wall_force_sample.reset();
    impl_->pending_diagnostic_local_flow_pattern.reset();
    impl_->committed_diagnostic_local_flow_pattern.reset();
    ++impl_->diagnostic_attempt_generation;
  }
  const auto record_diagnostic_attempt = [&](const auto &attempt_report) {
    if (impl_->domain == nullptr)
      return;
    impl_->last_diagnostic_state = &state;
    impl_->last_diagnostic_state_identity =
        state.diagnostic_mutation_identity();
    impl_->last_diagnostic_report_seal =
        diagnostic_snapshot_seal(
            attempt_report, physics.density_model,
            impl_->committed_diagnostic_wall_force_sample,
            impl_->committed_diagnostic_local_flow_pattern);
    impl_->last_diagnostic_rho_ref_kg_per_m3 = physics.rho_ref_kg_per_m3;
    impl_->last_diagnostic_density_model = physics.density_model;
  };
  const bool ordinary_density =
      (physics.density_model == config::DensityModel::constant ||
       physics.density_model == config::DensityModel::material) &&
      !physics.cp_J_per_kg_K.has_value() &&
      !physics.gas_constant_J_per_kg_K.has_value() &&
      !physics.thermodynamic_pressure_pa.has_value();
  const bool ideal_density =
      physics.density_model == config::DensityModel::ideal_gas &&
      impl_->density_setup.ideal_gas_closure != nullptr &&
      physics.cp_J_per_kg_K.has_value() &&
      physics.gas_constant_J_per_kg_K.has_value() &&
      physics.thermodynamic_pressure_pa.has_value() &&
      fp64_bits(*physics.cp_J_per_kg_K) ==
          fp64_bits(detail::DensityClosureAdapter::cp_J_per_kg_K(
              *impl_->density_setup.ideal_gas_closure)) &&
      fp64_bits(*physics.gas_constant_J_per_kg_K) ==
          fp64_bits(detail::DensityClosureAdapter::gas_constant_J_per_kg_K(
              *impl_->density_setup.ideal_gas_closure)) &&
      fp64_bits(*physics.thermodynamic_pressure_pa) ==
          fp64_bits(detail::DensityClosureAdapter::configured_pressure_pa(
              *impl_->density_setup.ideal_gas_closure));
  const bool local_supported =
      physics.density_model == impl_->density_setup.model &&
      (ordinary_density || ideal_density);
  runtime::CollectiveStatus supported{};
  try {
    supported = runtime::collective_status(
        *impl_->mpi, local_supported,
        "immersed-flow fixed-step physics is not implemented");
  } catch (const runtime::detail::MpiOperationError &) {
    StepAttemptReport failure{};
    failure.disposition = StepAttemptDisposition::non_retryable_failure;
    failure.reason = StepFailureReason::collective_operation;
    failure.lowest_failing_rank = -1;
    result.base = failure;
    record_diagnostic_attempt(result);
    return result;
  }
  if (!supported.ok) {
    StepAttemptReport failure{};
    failure.disposition = StepAttemptDisposition::non_retryable_failure;
    failure.reason = StepFailureReason::invalid_input;
    failure.lowest_failing_rank = supported.failing_rank;
    result.base = failure;
    record_diagnostic_attempt(result);
    return result;
  }
  if (impl_->domain != nullptr) {
    std::optional<MaterialDensityTransportReport> material_report;
    std::optional<MaterialDensityTransportReport> post_closure_report;
    std::optional<IdealGasClosureReport> closure_report;
    MaterialTransportFailureReason material_failure{
        MaterialTransportFailureReason::none};
    std::uint64_t material_attempt_identity{};
    auto report = impl_->attempt_immersed(
        state, physics, stencil, momentum_control, pressure_control,
        &result.force, &result.wale, &material_report,
        &post_closure_report, &closure_report, &material_failure,
        &material_attempt_identity);
    result.linear_solve_failure = impl_->last_linear_solve_failure;
    diagnostics::Stage3PerformanceCounters increment;
    const auto operator_work = impl_->immersed_operator->report();
    increment.step_ghost_constraints =
        operator_work.algebraic_occurrence_count;
    increment.step_lfp_transforms = operator_work.replacement_group_count;
    increment.step_immersed_rows = operator_work.algebraic_occurrence_count;
    if (report.pressure_corrector_count != 0U) {
      if (impl_->wall_links.size() >
          std::numeric_limits<std::uint64_t>::max() /
              report.pressure_corrector_count)
        throw runtime::Error(
            "immersed pressure wall work counter would overflow");
      increment.step_pressure_wall_constraints =
          static_cast<std::uint64_t>(impl_->wall_links.size()) *
          report.pressure_corrector_count;
    }
    if (result.force.has_value()) {
      increment.step_wall_quadrature_evaluations =
          impl_->wall_quadrature->local_points().size();
      increment.step_force_reductions = 3U;
    }
    const auto add_work = [](std::uint64_t &target, std::uint64_t value) {
      if (value > std::numeric_limits<std::uint64_t>::max() - target)
        throw runtime::Error(
            "immersed-flow performance counter would overflow");
      target += value;
    };
    add_work(impl_->performance.step_ghost_constraints,
             increment.step_ghost_constraints);
    add_work(impl_->performance.step_lfp_transforms,
             increment.step_lfp_transforms);
    add_work(impl_->performance.step_immersed_rows,
             increment.step_immersed_rows);
    add_work(impl_->performance.step_pressure_wall_constraints,
             increment.step_pressure_wall_constraints);
    add_work(impl_->performance.step_wall_quadrature_evaluations,
             increment.step_wall_quadrature_evaluations);
    add_work(impl_->performance.step_force_reductions,
             increment.step_force_reductions);
    if (impl_->density_setup.model == config::DensityModel::material) {
      const auto *adapter =
          std::get_if<detail::ImmersedMaterialDensityAdapter>(
              &impl_->density_adapter);
      if (adapter == nullptr)
        throw runtime::Error("immersed material adapter is unavailable");
      result.base = detail::DensityClosureBridge::make_material_report(
          std::move(report), std::move(material_report), material_failure,
          adapter->material_field_count(),
          adapter->shared_face_mass_flux_field(), material_attempt_identity);
    } else if (impl_->density_setup.model ==
               config::DensityModel::ideal_gas) {
      const auto *adapter =
          std::get_if<detail::ImmersedIdealGasDensityAdapter>(
              &impl_->density_adapter);
      if (adapter == nullptr)
        throw runtime::Error("immersed ideal-gas adapter is unavailable");
      auto material =
          detail::DensityClosureBridge::make_material_closure_report(
              std::move(report), std::move(material_report),
              std::move(post_closure_report), material_failure,
              adapter->material_field_count(),
              adapter->shared_face_mass_flux_field(),
              material_attempt_identity);
      result.base = detail::DensityClosureAdapter::make_ideal_report(
          std::move(material), std::move(closure_report),
          material_attempt_identity);
    } else {
      result.base = std::move(report);
    }
    const auto committed = std::visit(
        [](const auto &base) {
          using Base = std::decay_t<decltype(base)>;
          if constexpr (std::is_same_v<Base, StepAttemptReport>)
            return base.disposition == StepAttemptDisposition::committed;
          else if constexpr (
              std::is_same_v<Base, MaterialDensityStepAttemptReport>)
            return base.flow().disposition ==
                   StepAttemptDisposition::committed;
          else
            return base.flow().flow().disposition ==
                   StepAttemptDisposition::committed;
        },
        result.base);
    if (!committed) {
      result.force.reset();
      result.wale.reset();
    }
    record_diagnostic_attempt(result);
    return result;
  }
  if (impl_->density_setup.model != config::DensityModel::constant) {
    if (!impl_->body_fitted_material.has_value()) {
      StepAttemptReport failure{};
      failure.disposition = StepAttemptDisposition::non_retryable_failure;
      failure.reason = StepFailureReason::invalid_input;
      result.base = failure;
      return result;
    }
    impl_->last_wale_evaluation_count = 0U;
    impl_->last_wale_coefficient_identity = {};
    les::WaleSummary wale_summary;
    auto report = detail::DensityClosureBridge::attempt_with_optional_wale(
        *impl_->body_fitted_material, state,
        physics.dynamic_viscosity_pa_s, stencil, momentum_control,
        pressure_control,
        impl_->density_setup.model == config::DensityModel::ideal_gas
            ? &*impl_->body_fitted_closure_hooks
            : nullptr,
        impl_->wale,
        impl_->wale == nullptr ? nullptr : &wale_summary);
    impl_->last_wale_evaluation_count =
        detail::DensityClosureBridge::wale_evaluation_count(
            *impl_->body_fitted_material);
    impl_->last_corrector_count = report.flow().pressure_corrector_count;
    if (report.flow().disposition == StepAttemptDisposition::committed &&
        impl_->wale != nullptr) {
      result.wale = wale_summary;
      impl_->last_wale_coefficient_identity = wale_summary.identity;
    }
    if (impl_->density_setup.model == config::DensityModel::ideal_gas) {
      const auto attempt_identity = report.attempt_identity();
      std::optional<IdealGasClosureReport> closure_report;
      auto latest = detail::DensityClosureAdapter::latest_report(
          *impl_->density_setup.ideal_gas_closure);
      if (latest.attempt_identity() == attempt_identity)
        closure_report = std::move(latest);
      result.base = detail::DensityClosureAdapter::make_ideal_report(
          std::move(report), std::move(closure_report), attempt_identity);
    } else {
      result.base = std::move(report);
    }
    return result;
  }
  if (!impl_->body_fitted.has_value()) {
    StepAttemptReport failure{};
    failure.disposition = StepAttemptDisposition::non_retryable_failure;
    failure.reason = StepFailureReason::invalid_input;
    result.base = failure;
    return result;
  }
  impl_->last_wale_evaluation_count = 0U;
  impl_->last_wale_coefficient_identity = {};
  StepAttemptReport report;
  if (impl_->wale == nullptr) {
    report = impl_->body_fitted->attempt(
        state, physics.rho_ref_kg_per_m3, physics.dynamic_viscosity_pa_s,
        stencil, momentum_control, pressure_control);
  } else {
    les::WaleSummary wale_summary;
    report = impl_->body_fitted->attempt_with_wale(
        state, physics.rho_ref_kg_per_m3, physics.dynamic_viscosity_pa_s,
        stencil, momentum_control, pressure_control, *impl_->wale,
        wale_summary);
    if (report.disposition == StepAttemptDisposition::committed) {
      result.wale = wale_summary;
      impl_->last_wale_evaluation_count = 1U;
      impl_->last_wale_coefficient_identity = wale_summary.identity;
    }
  }
  impl_->last_corrector_count = report.pressure_corrector_count;
  result.base = std::move(report);
  return result;
}

ImmersedFlowDiagnosticSource
FixedStepImmersedFlow::diagnostic_source(
    const FlowState &state,
    const ImmersedFlowStepAttemptReport &report) const {
  if (!impl_ || !impl_->domain || !impl_->ghost_plan ||
      impl_->last_diagnostic_state != &state || state.attempt_active() ||
      state.diagnostic_mutation_identity() == 0U ||
      state.diagnostic_mutation_identity() !=
          impl_->last_diagnostic_state_identity ||
      diagnostic_snapshot_seal(
          report, impl_->last_diagnostic_density_model,
          impl_->committed_diagnostic_wall_force_sample,
          impl_->committed_diagnostic_local_flow_pattern) == 0U ||
      diagnostic_snapshot_seal(
          report, impl_->last_diagnostic_density_model,
          impl_->committed_diagnostic_wall_force_sample,
          impl_->committed_diagnostic_local_flow_pattern) !=
          impl_->last_diagnostic_report_seal) {
    throw runtime::Error("immersed-flow diagnostic source is stale");
  }
  ImmersedFlowDiagnosticSource source;
  source.attempt_generation_ = &impl_->diagnostic_attempt_generation;
  source.expected_attempt_generation_ = impl_->diagnostic_attempt_generation;
  source.report_ = report;
  source.rank_ = impl_->mpi->rank();
  source.step_ = state.metadata().step;
  source.time_s_ = state.metadata().time_s;

  const auto &access = state.solver_access_plan();
  const auto &committed = state.layer(FlowLayer::committed);
  const auto flux = committed.acquire_face_read<double>(
      access, kStatePhase, kStateActor, state.fields().face_mass_flux);
  double penetration_sum = 0.0;
  for (const auto &wall : impl_->wall_links) {
    const double denominator =
        impl_->last_diagnostic_rho_ref_kg_per_m3 * wall.effective_measure_m2;
    if (!(denominator > 0.0) || !std::isfinite(denominator)) {
      source.maximum_wall_penetration_m_per_s_ =
          std::numeric_limits<double>::quiet_NaN();
      source.mean_wall_penetration_m_per_s_ =
          std::numeric_limits<double>::quiet_NaN();
      break;
    }
    const double penetration = std::abs(flux(wall.face, 0) / denominator);
    if (!std::isfinite(penetration)) {
      source.maximum_wall_penetration_m_per_s_ =
          std::numeric_limits<double>::quiet_NaN();
      source.mean_wall_penetration_m_per_s_ =
          std::numeric_limits<double>::quiet_NaN();
      break;
    }
    source.maximum_wall_penetration_m_per_s_ =
        std::max(source.maximum_wall_penetration_m_per_s_, penetration);
    penetration_sum += penetration;
  }
  if (!impl_->wall_links.empty() &&
      std::isfinite(source.maximum_wall_penetration_m_per_s_))
    source.mean_wall_penetration_m_per_s_ =
        penetration_sum / static_cast<double>(impl_->wall_links.size());
  source.classified_cell_count_ = impl_->topology->owned_cell_count();
  source.active_cell_count_ = impl_->domain->active_cells().owned_active_count();
  source.immersed_link_count_ = impl_->wall_links.size();
  for (const auto &wall : impl_->wall_links) {
    source.donor_reference_count_ +=
        impl_->ghost_plan->velocity_constraint(wall.link.id, 0U).donors.size();
    source.donor_reference_count_ +=
        impl_->ghost_plan->zero_normal_constraint(wall.link.id).donors.size();
    source.donor_reference_count_ +=
        impl_->ghost_plan->density_extrapolation(wall.link.id).donors.size();
  }
  source.wall_quadrature_point_count_ =
      impl_->wall_quadrature == nullptr
          ? 0U
          : impl_->wall_quadrature->local_points().size();
  source.density_model_ = impl_->last_diagnostic_density_model;
  source.wall_force_available_ =
      impl_->committed_diagnostic_wall_force_sample.has_value() &&
      report.force.has_value();
  if (source.wall_force_available_)
    source.wall_force_sample_ = *impl_->committed_diagnostic_wall_force_sample;
  source.local_flow_pattern_available_ =
      impl_->committed_diagnostic_local_flow_pattern.has_value();
  if (source.local_flow_pattern_available_) {
    const auto &snapshot = *impl_->committed_diagnostic_local_flow_pattern;
    source.local_flow_pattern_algorithm_fingerprint_ =
        impl_->transform->algorithm_fingerprint();
    source.local_flow_pattern_row_fingerprint_ = snapshot.row_fingerprint;
    source.local_flow_pattern_replacement_group_count_ =
        snapshot.replacement_group_count;
    source.local_flow_pattern_algebraic_occurrence_count_ =
        snapshot.algebraic_occurrence_count;
    source.local_flow_pattern_replacement_coefficient_l2_ =
        snapshot.replacement_coefficient_l2;
    source.local_flow_pattern_limiting_case_status_ =
        snapshot.limiting_case_status;
  }
  source.snapshot_seal_ = impl_->last_diagnostic_report_seal;
  source.reduction_counters_ = impl_->mpi->fp64_reduction_counters();
  source.halo_counters_ = impl_->halo->performance_counters();
  source.allocation_counters_ = execution::allocation_counters();
  return source;
}

diagnostics::Stage3PerformanceCounters
FixedStepImmersedFlow::performance_counters() const noexcept {
  return impl_ ? impl_->performance
               : diagnostics::Stage3PerformanceCounters{};
}

void ImmersedFlowDiagnosticSource::validate() const {
  if (attempt_generation_ == nullptr ||
      *attempt_generation_ != expected_attempt_generation_)
    throw runtime::Error("immersed-flow diagnostic source is stale");
}

const ImmersedFlowStepAttemptReport &
ImmersedFlowDiagnosticSource::report() const {
  validate();
  return report_;
}
int ImmersedFlowDiagnosticSource::rank() const {
  validate();
  return rank_;
}
std::uint64_t ImmersedFlowDiagnosticSource::committed_step() const {
  validate();
  return step_;
}
double ImmersedFlowDiagnosticSource::committed_time_s() const {
  validate();
  return time_s_;
}
double ImmersedFlowDiagnosticSource::maximum_wall_penetration_m_per_s() const {
  validate();
  return maximum_wall_penetration_m_per_s_;
}
double ImmersedFlowDiagnosticSource::mean_wall_penetration_m_per_s() const {
  validate();
  return mean_wall_penetration_m_per_s_;
}
std::uint64_t ImmersedFlowDiagnosticSource::classified_cell_count() const {
  validate();
  return classified_cell_count_;
}
std::uint64_t ImmersedFlowDiagnosticSource::active_cell_count() const {
  validate();
  return active_cell_count_;
}
std::uint64_t ImmersedFlowDiagnosticSource::immersed_link_count() const {
  validate();
  return immersed_link_count_;
}
std::uint64_t ImmersedFlowDiagnosticSource::donor_reference_count() const {
  validate();
  return donor_reference_count_;
}
std::uint64_t
ImmersedFlowDiagnosticSource::wall_quadrature_point_count() const {
  validate();
  return wall_quadrature_point_count_;
}
config::DensityModel ImmersedFlowDiagnosticSource::density_model() const {
  validate();
  return density_model_;
}
bool ImmersedFlowDiagnosticSource::wall_force_available() const {
  validate();
  return wall_force_available_;
}
const immersed::WallForceSample &
ImmersedFlowDiagnosticSource::wall_force_sample() const {
  validate();
  if (!wall_force_available_)
    throw runtime::Error("immersed-flow wall-force snapshot is unavailable");
  return wall_force_sample_;
}
bool ImmersedFlowDiagnosticSource::local_flow_pattern_available() const {
  validate();
  return local_flow_pattern_available_;
}
std::uint64_t
ImmersedFlowDiagnosticSource::local_flow_pattern_algorithm_fingerprint() const {
  validate();
  return local_flow_pattern_algorithm_fingerprint_;
}
std::uint64_t
ImmersedFlowDiagnosticSource::local_flow_pattern_row_fingerprint() const {
  validate();
  return local_flow_pattern_row_fingerprint_;
}
std::uint64_t ImmersedFlowDiagnosticSource::
    local_flow_pattern_replacement_group_count() const {
  validate();
  return local_flow_pattern_replacement_group_count_;
}
std::uint64_t ImmersedFlowDiagnosticSource::
    local_flow_pattern_algebraic_occurrence_count() const {
  validate();
  return local_flow_pattern_algebraic_occurrence_count_;
}
double ImmersedFlowDiagnosticSource::
    local_flow_pattern_replacement_coefficient_l2() const {
  validate();
  return local_flow_pattern_replacement_coefficient_l2_;
}
std::uint64_t ImmersedFlowDiagnosticSource::
    local_flow_pattern_limiting_case_status() const {
  validate();
  return local_flow_pattern_limiting_case_status_;
}
std::uint64_t ImmersedFlowDiagnosticSource::snapshot_seal() const {
  validate();
  return snapshot_seal_;
}
runtime::Fp64ReductionCounters
ImmersedFlowDiagnosticSource::reduction_counters() const {
  validate();
  return reduction_counters_;
}
runtime::HaloPerformanceCounters
ImmersedFlowDiagnosticSource::halo_counters() const {
  validate();
  return halo_counters_;
}
execution::AllocationCounters
ImmersedFlowDiagnosticSource::allocation_counters() const {
  validate();
  return allocation_counters_;
}

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
void detail::ImmersedFlowAccess::set_linear_solve_failure(
    ImmersedLinearSolvePhase phase, std::uint32_t component) noexcept {
  injected_linear_failure_component.store(component,
                                           std::memory_order_relaxed);
  injected_linear_failure_phase.store(static_cast<int>(phase),
                                      std::memory_order_relaxed);
}

void detail::ImmersedFlowAccess::clear_linear_solve_failure() noexcept {
  injected_linear_failure_phase.store(-1, std::memory_order_relaxed);
  injected_linear_failure_component.store(
      std::numeric_limits<std::uint32_t>::max(), std::memory_order_relaxed);
}

void detail::ImmersedFlowAccess::set_failure_stage(
    std::uint8_t stage, int failing_rank) noexcept {
  injected_failure_rank.store(failing_rank, std::memory_order_relaxed);
  injected_failure_stage.store(static_cast<AttemptFailureStage>(stage),
                               std::memory_order_relaxed);
}

void detail::ImmersedFlowAccess::clear_failure_stage() noexcept {
  injected_failure_stage.store(AttemptFailureStage::none,
                               std::memory_order_relaxed);
  injected_failure_rank.store(-1, std::memory_order_relaxed);
}

void detail::ImmersedFlowAccess::set_wall_input_failure(
    std::uint8_t failure, int failing_rank) noexcept {
  injected_wall_failure_rank.store(failing_rank, std::memory_order_relaxed);
  injected_wall_failure.store(static_cast<WallInputFault>(failure),
                              std::memory_order_relaxed);
}

void detail::ImmersedFlowAccess::clear_wall_input_failure() noexcept {
  injected_wall_failure.store(WallInputFault::none,
                              std::memory_order_relaxed);
  injected_wall_failure_rank.store(-1, std::memory_order_relaxed);
}

ForceAttemptReport
detail::ImmersedFlowAccess::assemble_force_attempt_report_from_budget(
    const immersed::ForceComponents &budget_reaction,
    const immersed::ForceComponents &surface_traction) {
  return assemble_candidate_force_report_from_budget(budget_reaction,
                                                     surface_traction);
}

void detail::ImmersedFlowAccess::set_momentum_time_diagonal_scale(
    double scale) noexcept {
  injected_momentum_time_diagonal_scale.store(scale, std::memory_order_relaxed);
}

void detail::ImmersedFlowAccess::clear_momentum_time_diagonal_scale() noexcept {
  injected_momentum_time_diagonal_scale.store(1.0, std::memory_order_relaxed);
}

void detail::ImmersedFlowAccess::set_manufactured_body_source(
    FixedStepImmersedFlow &flow, std::vector<double> owned_active_source) {
  if (!flow.impl_)
    throw runtime::Error("immersed-flow test flow is empty");
  flow.impl_->manufactured_body_source = std::move(owned_active_source);
}

void detail::ImmersedFlowAccess::clear_manufactured_body_source(
    FixedStepImmersedFlow &flow) noexcept {
  if (flow.impl_)
    flow.impl_->manufactured_body_source.clear();
}

std::uint64_t detail::ImmersedFlowAccess::pressure_revision(
    const FixedStepImmersedFlow &flow) noexcept {
  return flow.impl_ ? flow.impl_->pressure_revision : 0U;
}

bool detail::ImmersedFlowAccess::has_active_pressure_reference(
    const FixedStepImmersedFlow &flow) noexcept {
  return flow.impl_ && flow.impl_->exact_pressure_operator &&
         flow.impl_->exact_pressure_operator->has_pressure_reference();
}

std::uint64_t detail::ImmersedFlowAccess::pressure_apply_schedule(
    const FixedStepImmersedFlow &flow) noexcept {
  return flow.impl_ && flow.impl_->exact_pressure_operator
             ? flow.impl_->exact_pressure_operator->last_apply_schedule()
             : 0U;
}

std::uint32_t detail::ImmersedFlowAccess::last_corrector_count(
    const FixedStepImmersedFlow &flow) noexcept {
  return flow.impl_ ? flow.impl_->last_corrector_count : 0U;
}

std::uint64_t detail::ImmersedFlowAccess::wale_evaluation_count(
    const FixedStepImmersedFlow &flow) noexcept {
  return flow.impl_ ? flow.impl_->last_wale_evaluation_count : 0U;
}

les::WaleCoefficientIdentity
detail::ImmersedFlowAccess::wale_coefficient_identity(
    const FixedStepImmersedFlow &flow) noexcept {
  return flow.impl_ ? flow.impl_->last_wale_coefficient_identity
                    : les::WaleCoefficientIdentity{};
}

std::uint64_t
detail::ImmersedFlowAccess::wall_effective_viscosity_fingerprint(
    const FixedStepImmersedFlow &flow) noexcept {
  return flow.impl_ ? flow.impl_->last_wall_effective_viscosity_fingerprint
                    : 0U;
}

double detail::ImmersedFlowAccess::last_wall_pressure_gradient_application_norm(
    const FixedStepImmersedFlow &flow) noexcept {
  return flow.impl_ ? flow.impl_->last_wall_pressure_gradient_application_norm
                    : 0.0;
}

std::vector<detail::ImmersedFlowWallMeasure>
detail::ImmersedFlowAccess::wall_effective_measures(
    const FixedStepImmersedFlow &flow) {
  std::vector<detail::ImmersedFlowWallMeasure> result;
  if (!flow.impl_)
    return result;
  result.reserve(flow.impl_->wall_links.size());
  for (const auto &wall : flow.impl_->wall_links)
    result.push_back({wall.link.id, wall.effective_measure_m2});
  return result;
}

std::uint64_t
detail::ImmersedFlowAccess::immersed_operator_structure_fingerprint(
    const FixedStepImmersedFlow &flow) {
#ifdef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS
  if (!flow.impl_ || !flow.impl_->immersed_operator.has_value())
    return 0U;
  struct RowWire final {
    std::uint64_t active_cell{};
    std::uint64_t fingerprint{};
  };
  static_assert(sizeof(RowWire) == 2U * sizeof(std::uint64_t));
  std::vector<RowWire> local;
  for (const auto &row :
       finite_volume::detail::ImmersedBoundaryAuthorityAccess::rows(
           *flow.impl_->immersed_operator)) {
    std::uint64_t fingerprint = kFnvOffset;
    fingerprint = hash_u64(fingerprint, row.active_cell);
    fingerprint = hash_u64(fingerprint, row.row_replacement_fingerprint);
    fingerprint = hash_u64(fingerprint, row.replacement_group_count);
    fingerprint =
        hash_u64(fingerprint, static_cast<std::uint64_t>(row.links.size()));
    for (const auto &link : row.links) {
      fingerprint = hash_u64(fingerprint, link.id);
      fingerprint =
          hash_u64(fingerprint, static_cast<std::uint64_t>(link.occurrence));
      fingerprint = hash_u64(fingerprint, link.replacement_fingerprint);
    }
    local.push_back({row.active_cell, fingerprint});
  }
  if (local.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()) /
                         sizeof(RowWire))
    throw runtime::Error(
        "immersed-flow operator fingerprint payload exceeds MPI int");
  const int local_bytes = static_cast<int>(local.size() * sizeof(RowWire));
  std::vector<int> counts(static_cast<std::size_t>(flow.impl_->mpi->size()));
  runtime::detail::check_mpi(
      MPI_Allgather(&local_bytes, 1, MPI_INT, counts.data(), 1, MPI_INT,
                    flow.impl_->mpi->comm()),
      "MPI_Allgather(immersed-flow operator fingerprint counts)");
  std::vector<int> offsets(counts.size(), 0);
  std::size_t total_bytes = 0U;
  for (std::size_t rank = 0U; rank < counts.size(); ++rank) {
    if (counts[rank] < 0 ||
        total_bytes > static_cast<std::size_t>(std::numeric_limits<int>::max()))
      throw runtime::Error(
          "immersed-flow operator fingerprint payload is unsupported");
    offsets[rank] = static_cast<int>(total_bytes);
    total_bytes += static_cast<std::size_t>(counts[rank]);
  }
  if (total_bytes % sizeof(RowWire) != 0U ||
      total_bytes > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    throw runtime::Error(
        "immersed-flow operator fingerprint payload size is invalid");
  std::vector<RowWire> gathered(total_bytes / sizeof(RowWire));
  runtime::detail::check_mpi(MPI_Allgatherv(local.data(), local_bytes, MPI_BYTE,
                                            gathered.data(), counts.data(),
                                            offsets.data(), MPI_BYTE,
                                            flow.impl_->mpi->comm()),
                             "MPI_Allgatherv(immersed-flow operator fingerprints)");
  std::sort(gathered.begin(), gathered.end(),
            [](const auto &left, const auto &right) {
              return left.active_cell < right.active_cell;
            });
  if (std::adjacent_find(gathered.begin(), gathered.end(),
                         [](const auto &left, const auto &right) {
                           return left.active_cell == right.active_cell;
                         }) != gathered.end())
    throw runtime::Error("immersed-flow operator fingerprint row is duplicated");
  std::uint64_t result = kFnvOffset;
  result = hash_u64(result, static_cast<std::uint64_t>(gathered.size()));
  for (const auto &row : gathered) {
    result = hash_u64(result, row.active_cell);
    result = hash_u64(result, row.fingerprint);
  }
  return result;
#else
  static_cast<void>(flow);
  return 0U;
#endif
}

const finite_volume::ImmersedOperatorAdapter &
detail::ImmersedFlowAccess::immersed_operator(
    const FixedStepImmersedFlow &flow) {
#ifdef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS
  if (!flow.impl_ || !flow.impl_->immersed_operator.has_value())
    throw runtime::Error("immersed-flow immersed operator is unavailable");
  return *flow.impl_->immersed_operator;
#else
  static_cast<void>(flow);
  throw runtime::Error("immersed-flow immersed operator test access is unavailable");
#endif
}

std::vector<detail::ImmersedFlowWallGradient>
detail::ImmersedFlowAccess::final_pressure_wall_gradients(
    const FixedStepImmersedFlow &flow) {
  if (!flow.impl_)
    throw runtime::Error("immersed-flow test flow is empty");
  std::vector<detail::ImmersedFlowWallGradient> result;
  result.reserve(flow.impl_->committed_final_pressure_wall_gradients.size());
  for (const auto &gradient :
       flow.impl_->committed_final_pressure_wall_gradients)
    result.push_back({gradient.link, gradient.value});
  return result;
}

std::array<double, 3> detail::ImmersedFlowAccess::active_pressure_algebra_probe(
    const FixedStepImmersedFlow &flow) {
  if (!flow.impl_ || !flow.impl_->exact_pressure_operator ||
      !flow.impl_->domain ||
      !flow.impl_->ghost_plan || !flow.impl_->topology ||
      !flow.impl_->geometry || !flow.impl_->boundaries ||
      !flow.impl_->execution || !flow.impl_->mpi)
    throw runtime::Error("immersed-flow pressure algebra probe is unavailable");
  auto probe = make_unit_mobility_pressure_probe(
      *flow.impl_->topology, *flow.impl_->geometry, *flow.impl_->decomposition,
      *flow.impl_->domain, *flow.impl_->boundaries, *flow.impl_->mpi,
      *flow.impl_->execution, flow.impl_->owned_active_cells);
  const auto &ids = flow.impl_->domain->active_cells().ordered_global_ids();
  return pressure_operator_algebra_probe(*probe, ids, *flow.impl_->execution,
                                         *flow.impl_->mpi);
}

std::array<double, 3>
detail::ImmersedFlowAccess::active_pressure_current_algebra_probe(
    const FixedStepImmersedFlow &flow) {
  if (!flow.impl_ || !flow.impl_->pressure_operator || !flow.impl_->domain ||
      !flow.impl_->execution || !flow.impl_->mpi ||
      flow.impl_->pressure_operator->revision() == 0U)
    throw runtime::Error(
        "immersed-flow current pressure algebra probe is unavailable");
  const auto &ids = flow.impl_->domain->active_cells().ordered_global_ids();
  return pressure_operator_algebra_probe(*flow.impl_->pressure_operator, ids,
                                         *flow.impl_->execution,
                                         *flow.impl_->mpi);
}

detail::ImmersedFlowActiveExchangeData
detail::ImmersedFlowAccess::active_exchange_performance(
    const FixedStepImmersedFlow &flow) noexcept {
  detail::ImmersedFlowActiveExchangeData result;
  if (!flow.impl_ || !flow.impl_->domain || !flow.impl_->pressure_operator)
    return result;
  const auto &active = flow.impl_->domain->active_cells();
  result.owned_active_count = active.owned_active_count();
  result.ghost_active_count =
      active.ordered_global_ids().size() - active.owned_active_count();
  result.pressure = flow.impl_->pressure_operator->performance_counters();
  for (std::size_t component = 0U; component < result.momentum.size();
       ++component)
    if (flow.impl_->momentum_operators[component])
      result.momentum[component] =
          flow.impl_->momentum_operators[component]->performance_counters();
  return result;
}

std::array<double, 9>
detail::ImmersedFlowAccess::active_momentum_reference_probe(
    const FixedStepImmersedFlow &flow, double rho_ref_kg_per_m3,
    double dynamic_viscosity_pa_s, const MomentumTimeStencil &stencil) {
  if (!flow.impl_ || !flow.impl_->domain || !flow.impl_->topology ||
      !flow.impl_->geometry || !flow.impl_->execution || !flow.impl_->mpi ||
      !(rho_ref_kg_per_m3 > 0.0) || !std::isfinite(rho_ref_kg_per_m3) ||
      !(dynamic_viscosity_pa_s >= 0.0) ||
      !std::isfinite(dynamic_viscosity_pa_s))
    throw runtime::Error(
        "immersed-flow implicit viscous reference probe is unavailable");
  const std::size_t count = flow.impl_->owned_active_cells.size();
  std::vector<double> time_diagonal(count);
  for (std::size_t row = 0U; row < count; ++row)
    time_diagonal[row] = stencil.alpha0 * rho_ref_kg_per_m3 *
                         flow.impl_->geometry->cell_volume_m3(
                             flow.impl_->owned_active_cells[row]) /
                         stencil.dt_s;
  auto wall_coefficients = flow.impl_->active_viscous_wall_coefficients();
  ActiveMomentumOperator probe(
      *flow.impl_->topology, *flow.impl_->geometry, *flow.impl_->decomposition,
      flow.impl_->domain->active_cells(), wall_coefficients,
      *flow.impl_->execution);
  probe.replace(time_diagonal, dynamic_viscosity_pa_s);
  execution::Buffer x(*flow.impl_->execution, checked_bytes(count));
  execution::Buffer y(*flow.impl_->execution, checked_bytes(count));
  execution::Buffer ax(*flow.impl_->execution, checked_bytes(count));
  execution::Buffer ay(*flow.impl_->execution, checked_bytes(count));
  execution::Buffer trial(*flow.impl_->execution, checked_bytes(count));
  execution::Buffer atrial(*flow.impl_->execution, checked_bytes(count));
  execution::Buffer diagonal(*flow.impl_->execution, checked_bytes(count));
  auto xv = x.view(0U, count);
  auto yv = y.view(0U, count);
  auto trial_values = trial.view(0U, count);
  const auto &ids = flow.impl_->domain->active_cells().ordered_global_ids();
  for (std::size_t row = 0U; row < count; ++row) {
    const double phase = static_cast<double>((ids[row] % 104729U) + 1U);
    xv[row] = std::sin(phase * 0.017);
    yv[row] = std::cos(phase * 0.013);
    trial_values[row] = stencil.order == MomentumTimeOrder::backward_euler
                            ? xv[row]
                            : 2.0 * xv[row] - yv[row];
  }
  std::unordered_map<mesh::GlobalCellId, std::size_t> owned_rows;
  owned_rows.reserve(count);
  for (std::size_t row = 0U; row < count; ++row)
    owned_rows.emplace(ids[row], row);
  double local_wall_reference_sum = 0.0;
  double local_wall_reference_scale = 0.0;
  for (const auto &wall : wall_coefficients) {
    const auto found = owned_rows.find(wall.cell);
    if (found == owned_rows.end())
      throw runtime::Error(
          "immersed-flow implicit viscous reference wall row is unavailable");
    const double value = dynamic_viscosity_pa_s *
                         wall.coefficient_per_viscosity * xv[found->second];
    local_wall_reference_sum += value;
    local_wall_reference_scale += std::abs(value);
  }
  probe
      .apply(static_cast<const execution::Buffer &>(x).view(0U, count),
             ax.view(0U, count))
      .wait();
  probe
      .apply(static_cast<const execution::Buffer &>(y).view(0U, count),
             ay.view(0U, count))
      .wait();
  probe
      .apply(static_cast<const execution::Buffer &>(trial).view(0U, count),
             atrial.view(0U, count))
      .wait();
  probe.diagonal(diagonal.view(0U, count)).wait();
  const auto axv = static_cast<const execution::Buffer &>(ax).view(0U, count);
  const auto ayv = static_cast<const execution::Buffer &>(ay).view(0U, count);
  const auto atrial_values =
      static_cast<const execution::Buffer &>(atrial).view(0U, count);
  const auto dv =
      static_cast<const execution::Buffer &>(diagonal).view(0U, count);
  double totals[6]{};
  double minimum_normalized_diagonal = std::numeric_limits<double>::infinity();
  double maximum_spatial_increment = 0.0;
  double maximum_deferred_cancellation_error = 0.0;
  double maximum_omitted_history_residual = 0.0;
  double off_diagonal_action_squared = 0.0;
  double conservation_totals[4]{0.0, local_wall_reference_sum,
                                local_wall_reference_scale, 0.0};
  for (std::size_t row = 0U; row < count; ++row) {
    totals[0] += xv[row] * ayv[row];
    totals[1] += yv[row] * axv[row];
    totals[2] += xv[row] * axv[row];
    totals[3] += yv[row] * ayv[row];
    totals[4] += xv[row] * xv[row];
    totals[5] += yv[row] * yv[row];
    const double volume = flow.impl_->geometry->cell_volume_m3(
        flow.impl_->owned_active_cells[row]);
    minimum_normalized_diagonal =
        std::min(minimum_normalized_diagonal, dv[row] / volume);
    maximum_spatial_increment = std::max(
        maximum_spatial_increment, (dv[row] - time_diagonal[row]) / volume);
    const double off_diagonal_action = axv[row] - dv[row] * xv[row];
    off_diagonal_action_squared += off_diagonal_action * off_diagonal_action;
    const bool backward_euler =
        stencil.order == MomentumTimeOrder::backward_euler;
    const double trial_term =
        atrial_values[row] - time_diagonal[row] * trial_values[row];
    const double current_term = axv[row] - time_diagonal[row] * xv[row];
    const double history_term = ayv[row] - time_diagonal[row] * yv[row];
    const double full = backward_euler
                            ? trial_term - current_term
                            : trial_term - 2.0 * current_term + history_term;
    const double scale =
        std::max(1.0, std::abs(trial_term) + 2.0 * std::abs(current_term) +
                          std::abs(history_term));
    maximum_deferred_cancellation_error =
        std::max(maximum_deferred_cancellation_error, std::abs(full) / scale);
    const double omitted_history =
        backward_euler ? trial_term : trial_term - 2.0 * current_term;
    maximum_omitted_history_residual = std::max(
        maximum_omitted_history_residual, std::abs(omitted_history) / volume);
    conservation_totals[0] += current_term;
    conservation_totals[2] += std::abs(current_term);
    conservation_totals[3] += xv[row] * current_term;
  }
  flow.impl_->mpi->allreduce_fp64_in_place(
      totals, 6U, runtime::Fp64ReductionOperation::sum);
  flow.impl_->mpi->allreduce_fp64_in_place(
      &off_diagonal_action_squared, 1U, runtime::Fp64ReductionOperation::sum);
  flow.impl_->mpi->allreduce_fp64_in_place(
      conservation_totals, 4U, runtime::Fp64ReductionOperation::sum);
  runtime::detail::check_mpi(
      MPI_Allreduce(MPI_IN_PLACE, &minimum_normalized_diagonal, 1, MPI_DOUBLE,
                    MPI_MIN, flow.impl_->mpi->comm()),
      "MPI_Allreduce(immersed-flow momentum minimum diagonal)");
  flow.impl_->mpi->allreduce_fp64_in_place(
      &maximum_spatial_increment, 1U, runtime::Fp64ReductionOperation::maximum);
  double cancellation_evidence[2]{maximum_deferred_cancellation_error,
                                  maximum_omitted_history_residual};
  flow.impl_->mpi->allreduce_fp64_in_place(
      cancellation_evidence, 2U, runtime::Fp64ReductionOperation::maximum);
  const double symmetry_scale =
      std::max({1.0, std::abs(totals[0]), std::abs(totals[1])});
  const double conservation_defect =
      conservation_totals[2] > 0.0
          ? std::abs(conservation_totals[0] - conservation_totals[1]) /
                conservation_totals[2]
          : (conservation_totals[0] == conservation_totals[1]
                 ? 0.0
                 : std::numeric_limits<double>::infinity());
  return {std::abs(totals[0] - totals[1]) / symmetry_scale,
          totals[2] / totals[4],
          minimum_normalized_diagonal,
          maximum_spatial_increment,
          cancellation_evidence[0],
          cancellation_evidence[1],
          std::sqrt(off_diagonal_action_squared /
                    std::max(totals[4], std::numeric_limits<double>::min())),
          conservation_defect,
          conservation_totals[3]};
}

std::vector<double> detail::ImmersedFlowAccess::active_pressure_operator_values(
    const FixedStepImmersedFlow &flow,
    std::vector<double> owned_active_pressure_pa,
    std::vector<detail::ImmersedFlowWallGradient> wall_gradients) {
  if (!flow.impl_ || !flow.impl_->exact_pressure_operator ||
      !flow.impl_->pressure_operator || !flow.impl_->domain ||
      !flow.impl_->ghost_plan || !flow.impl_->topology ||
      !flow.impl_->geometry || !flow.impl_->boundaries ||
      !flow.impl_->execution || !flow.impl_->mpi)
    throw runtime::Error("immersed-flow pressure operator probe is unavailable");
  const std::size_t count =
      flow.impl_->domain->active_cells().owned_active_count();
  if (owned_active_pressure_pa.size() != count)
    throw runtime::Error("immersed-flow pressure operator probe size is invalid");
  std::sort(wall_gradients.begin(), wall_gradients.end(),
            [](const auto &left, const auto &right) {
              return left.link < right.link;
            });
  if (std::adjacent_find(wall_gradients.begin(), wall_gradients.end(),
                         [](const auto &left, const auto &right) {
                           return left.link == right.link;
                         }) != wall_gradients.end())
    throw runtime::Error(
        "immersed-flow pressure operator probe gradient is duplicated");
  for (const auto &gradient : wall_gradients)
    if (!std::isfinite(gradient.normal_gradient_pa_per_m) ||
        gradient.normal_gradient_pa_per_m != 0.0)
      throw runtime::Error(
          "immersed-flow homogeneous pressure operator probe gradient is invalid");
  auto result = evaluate_pressure_operator_values(
      *flow.impl_->exact_pressure_operator, *flow.impl_->geometry,
      *flow.impl_->execution, flow.impl_->owned_active_cells,
      owned_active_pressure_pa);
  return result;
}

detail::ImmersedFlowPredictorWorkspaceData
detail::ImmersedFlowAccess::exact_predictor_workspace(
    const FixedStepImmersedFlow &flow) noexcept {
  if (!flow.impl_)
    return {};
  return {flow.impl_->exact_predictor_probe_creation_count,
          flow.impl_->exact_predictor_response_capacity_growth_count,
          flow.impl_->exact_pressure_operator
              ? flow.impl_->exact_pressure_operator
                    ->cached_response_consumption_count()
              : 0U,
          flow.impl_->exact_pressure_operator &&
              flow.impl_->exact_pressure_operator
                  ->cached_input_mutation_rejected(),
          flow.impl_->last_exact_predictor_response_allocation_event_count,
          flow.impl_->committed_exact_predictor_work.response_count,
          flow.impl_->committed_exact_predictor_work.solve_count,
          flow.impl_->committed_exact_predictor_work.iteration_count,
          flow.impl_->committed_exact_predictor_work.matvec_count,
          flow.impl_->committed_exact_predictor_work
              .preconditioner_apply_count,
          flow.impl_->committed_exact_predictor_work.global_reduction_count};
}

detail::ImmersedFlowOperatorDiagonalData
detail::ImmersedFlowAccess::exact_predictor_diagonal_contract(
    const FixedStepImmersedFlow &flow) {
  if (!flow.impl_ || !flow.impl_->exact_pressure_operator ||
      !flow.impl_->execution || !flow.impl_->mpi)
    throw runtime::Error(
        "immersed-flow exact predictor diagonal probe is unavailable");
  const auto &linear_operator = *flow.impl_->exact_pressure_operator;
  detail::ImmersedFlowOperatorDiagonalData report;
  report.advertised = linear_operator.has_diagonal();
  if (!report.advertised)
    return report;

  const auto layout = linear_operator.domain_layout();
  const std::size_t count = layout.owned_count();
  mesh::GlobalCellId local_minimum =
      std::numeric_limits<mesh::GlobalCellId>::max();
  if (count != 0U)
    local_minimum =
        *std::min_element(layout.global_ids().begin(),
                          std::next(layout.global_ids().begin(),
                                    static_cast<std::ptrdiff_t>(count)));
  mesh::GlobalCellId selected = local_minimum;
  runtime::detail::check_mpi(
      MPI_Allreduce(MPI_IN_PLACE, &selected, 1, MPI_UINT64_T, MPI_MIN,
                    flow.impl_->mpi->comm()),
      "MPI_Allreduce(immersed-flow exact predictor diagonal probe row)");
  if (selected == std::numeric_limits<mesh::GlobalCellId>::max())
    throw runtime::Error(
        "immersed-flow exact predictor diagonal probe layout is empty");

  execution::Buffer basis(*flow.impl_->execution, checked_bytes(count));
  execution::Buffer action(*flow.impl_->execution, checked_bytes(count));
  execution::Buffer declared(*flow.impl_->execution, checked_bytes(count));
  auto basis_values = basis.view(0U, count);
  for (std::size_t row = 0U; row < count; ++row)
    basis_values[row] = 0.0;
  std::optional<std::size_t> selected_row;
  for (std::size_t row = 0U; row < count; ++row) {
    if (layout.global_ids()[row] == selected) {
      selected_row = row;
      basis_values[row] = 1.0;
      break;
    }
  }
  linear_operator.diagonal(declared.view(0U, count)).wait();
  linear_operator
      .apply(static_cast<const execution::Buffer &>(basis).view(0U, count),
             action.view(0U, count))
      .wait();

  double values[3]{};
  if (selected_row.has_value()) {
    const auto declared_values =
        static_cast<const execution::Buffer &>(declared).view(0U, count);
    const auto action_values =
        static_cast<const execution::Buffer &>(action).view(0U, count);
    values[0] = declared_values[*selected_row];
    values[1] = action_values[*selected_row];
    values[2] = 1.0;
  }
  flow.impl_->mpi->allreduce_fp64_in_place(
      values, 3U, runtime::Fp64ReductionOperation::sum);
  if (values[2] != 1.0 || !std::isfinite(values[0]) ||
      !std::isfinite(values[1]))
    throw runtime::Error(
        "immersed-flow exact predictor diagonal probe result is invalid");
  report.declared_value = values[0];
  report.applied_basis_value = values[1];
  report.relative_difference =
      std::abs(values[0] - values[1]) /
      std::max({1.0, std::abs(values[0]), std::abs(values[1])});
  return report;
}

detail::ImmersedFlowPressureFluxIdentityData
detail::ImmersedFlowAccess::exact_predictor_schur_identity(
    FixedStepImmersedFlow &flow, const FlowState &state,
    double rho_ref_kg_per_m3, const MomentumTimeStencil &stencil,
    std::vector<double> owned_active_pressure_pa,
    std::vector<detail::ImmersedFlowWallGradient> wall_gradients) {
  if (!flow.impl_ || !flow.impl_->exact_pressure_operator ||
      !flow.impl_->pressure_operator || !flow.impl_->domain ||
      !flow.impl_->topology || !flow.impl_->geometry ||
      !flow.impl_->boundaries || !flow.impl_->execution || !flow.impl_->mpi ||
      !flow.impl_->halo || !flow.impl_->ghost_plan ||
      !(rho_ref_kg_per_m3 > 0.0) || !std::isfinite(rho_ref_kg_per_m3))
    throw runtime::Error(
        "immersed-flow exact predictor Schur probe is unavailable");
  const std::size_t count = flow.impl_->owned_active_cells.size();
  if (owned_active_pressure_pa.size() != count)
    throw runtime::Error(
        "immersed-flow exact predictor Schur pressure layout is invalid");
  std::sort(wall_gradients.begin(), wall_gradients.end(),
            [](const auto &left, const auto &right) {
              return left.link < right.link;
            });
  if (std::adjacent_find(wall_gradients.begin(), wall_gradients.end(),
                         [](const auto &left, const auto &right) {
                           return left.link == right.link;
                         }) != wall_gradients.end())
    throw runtime::Error(
        "immersed-flow exact predictor Schur wall gradient is duplicated");
  std::vector<finite_volume::detail::ImmersedWallNormalGradient> gradients;
  gradients.reserve(wall_gradients.size());
  for (const auto &gradient : wall_gradients) {
    if (!std::isfinite(gradient.normal_gradient_pa_per_m))
      throw runtime::Error(
          "immersed-flow exact predictor Schur wall gradient is non-finite");
    gradients.push_back(
        {gradient.link, gradient.normal_gradient_pa_per_m});
  }

  auto candidate_values = evaluate_pressure_operator_values(
      *flow.impl_->exact_pressure_operator, *flow.impl_->geometry,
      *flow.impl_->execution, flow.impl_->owned_active_cells,
      owned_active_pressure_pa);

  const auto box = flow.impl_->topology->owned_global_box();
  const runtime::FieldLayoutSet layout{
      {box.end.x - box.begin.x, box.end.y - box.begin.y,
       box.end.z - box.begin.z},
      flow.impl_->topology->local_face_count()};
  FlowLayerValues zeros;
  zeros.density.resize(flow.impl_->topology->owned_cell_count(),
                       rho_ref_kg_per_m3);
  zeros.velocity.resize(flow.impl_->topology->owned_cell_count() * 3U, 0.0);
  zeros.mechanical_pressure.resize(flow.impl_->topology->owned_cell_count(),
                                   0.0);
  zeros.face_velocity.resize(flow.impl_->topology->local_face_count() * 3U,
                             0.0);
  zeros.face_mass_flux.resize(flow.impl_->topology->local_face_count(), 0.0);
  auto probe_state = FlowState::create(
      detail::FlowStateSolverAccess::registry(state), layout, state.fields(),
      state.metadata());
  probe_state.seed_accepted_layers(zeros, zeros);
  auto &trial = detail::FlowStateSolverAccess::layer(probe_state,
                                                     FlowLayer::trial);
  const auto &probe_access =
      detail::FlowStateSolverAccess::access(probe_state);
  {
    auto pressure = trial.acquire_write<double>(
        probe_access, kStatePhase, kStateActor,
        state.fields().mechanical_pressure);
    for (std::size_t row = 0U; row < count; ++row) {
      const auto index =
          flow.impl_->field_index(flow.impl_->owned_active_cells[row]);
      pressure(index.i, index.j, index.k, 0) =
          owned_active_pressure_pa[row];
    }
  }
  flow.impl_->halo->exchange(trial, state.fields().mechanical_pressure);
  const auto pressure_residual = flow.impl_->assemble_immersed_pressure_residual(
      probe_state, trial, &gradients);
  if (pressure_residual.size() != count * 3U)
    throw runtime::Error(
        "immersed-flow exact predictor Schur pressure residual is invalid");

  std::array<std::vector<double>, 3> diagonal;
  std::array<std::vector<double>, 3> velocity_increment;
  for (std::size_t component_index = 0U; component_index < 3U;
       ++component_index) {
    execution::Buffer diagonal_buffer(*flow.impl_->execution,
                                      checked_bytes(count));
    flow.impl_->momentum_operators[component_index]
        ->diagonal(diagonal_buffer.view(0U, count))
        .wait();
    const auto diagonal_values =
        static_cast<const execution::Buffer &>(diagonal_buffer)
            .view(0U, count);
    diagonal[component_index].assign(diagonal_values.data(),
                                     diagonal_values.data() + count);
    execution::Buffer rhs_buffer(*flow.impl_->execution, checked_bytes(count));
    execution::Buffer solution_buffer(*flow.impl_->execution,
                                      checked_bytes(count));
    auto rhs_values = rhs_buffer.view(0U, count);
    auto solution_values = solution_buffer.view(0U, count);
    for (std::size_t row = 0U; row < count; ++row) {
      rhs_values[row] = -pressure_residual[row * 3U + component_index];
      solution_values[row] = 0.0;
    }
    linear::JacobiPreconditioner preconditioner(*flow.impl_->execution);
    preconditioner.update(*flow.impl_->momentum_operators[component_index],
                          flow.impl_->momentum_operators[component_index]
                              ->revision());
    linear::SolveControl control;
    control.atol = 1.0e-15;
    control.rtol = 1.0e-14;
    control.max_iterations = 500U;
    const auto solve = flow.impl_->momentum_solver->solve(
        *flow.impl_->momentum_operators[component_index], preconditioner,
        static_cast<const execution::Buffer &>(rhs_buffer).view(0U, count),
        solution_values, control);
    if (!flow.impl_->solve_succeeded(solve.reason))
      throw runtime::Error(
          "immersed-flow exact predictor Schur momentum solve failed");
    const auto solved_values =
        static_cast<const execution::Buffer &>(solution_buffer)
            .view(0U, count);
    velocity_increment[component_index].assign(
        solved_values.data(), solved_values.data() + count);
  }
  linear::SolveControl exact_control;
  exact_control.atol = 1.0e-15;
  exact_control.rtol = 1.0e-14;
  exact_control.max_iterations = 500U;
  flow.impl_->ensure_exact_predictor_workspace(state, rho_ref_kg_per_m3);
  const auto candidate_affine = flow.impl_->exact_predictor_response(
      rho_ref_kg_per_m3, stencil, diagonal,
      std::vector<double>(count, 0.0), gradients, exact_control,
      ImmersedLinearSolvePhase::pressure_affine_momentum, 1U);
  for (std::size_t row = 0U; row < count; ++row)
    candidate_values[row] += candidate_affine.divergence_per_volume[row];
  {
    auto velocity = trial.acquire_write<double>(
        probe_access, kStatePhase, kStateActor, state.fields().velocity);
    for (std::size_t row = 0U; row < count; ++row) {
      const auto index =
          flow.impl_->field_index(flow.impl_->owned_active_cells[row]);
      for (std::size_t component_index = 0U; component_index < 3U;
           ++component_index)
        velocity(index.i, index.j, index.k,
                 static_cast<int>(component_index)) =
            velocity_increment[component_index][row];
    }
  }
  flow.impl_->halo->exchange(trial, state.fields().velocity);
  static_cast<void>(flow.impl_->assemble_face_predictor(
      probe_state, rho_ref_kg_per_m3, stencil, diagonal, &gradients));
  const auto face_flux = trial.acquire_face_read<double>(
      probe_access, kStatePhase, kStateActor, state.fields().face_mass_flux);
  std::vector<double> literal_values(count, 0.0);
  for (mesh::LocalFaceId face = 0U;
       face < flow.impl_->topology->local_face_count(); ++face) {
    const double value = face_flux(face, 0);
    const auto owner_row = flow.impl_->domain->active_cells().active_index(
        flow.impl_->topology->owner(face));
    if (owner_row.has_value() && *owner_row < count)
      literal_values[*owner_row] += value;
    const auto neighbour = flow.impl_->topology->neighbour(face);
    if (!neighbour.has_value() ||
        flow.impl_->topology->periodic_pair(face).has_value())
      continue;
    const auto neighbour_row =
        flow.impl_->domain->active_cells().active_index(*neighbour);
    if (neighbour_row.has_value() && *neighbour_row < count)
      literal_values[*neighbour_row] -= value;
  }
  for (std::size_t row = 0U; row < count; ++row)
    literal_values[row] /= flow.impl_->geometry->cell_volume_m3(
        flow.impl_->owned_active_cells[row]);
  return make_pressure_flux_identity_report(
      *flow.impl_->topology, *flow.impl_->mpi, flow.impl_->owned_active_cells,
      candidate_values, literal_values);
}

detail::ImmersedFlowSpatialEnergyData
detail::ImmersedFlowAccess::spatial_energy_terms(
    FixedStepImmersedFlow &flow, const FlowState &state, FlowLayer selected,
    double rho_ref_kg_per_m3, double dynamic_viscosity_pa_s) {
  if (!flow.impl_ || !flow.impl_->domain || !flow.impl_->topology ||
      !flow.impl_->geometry || !flow.impl_->halo || !flow.impl_->mpi ||
      !flow.impl_->scratch || !flow.impl_->scratch->storage ||
      !flow.impl_->scratch->access)
    throw runtime::Error("immersed-flow spatial-energy probe is unavailable");
  return flow.impl_->spatial_energy_terms(state, selected, rho_ref_kg_per_m3,
                                          dynamic_viscosity_pa_s);
}

detail::ImmersedFlowExactMomentumResidualData
detail::ImmersedFlowAccess::exact_momentum_residual_terms(
    FixedStepImmersedFlow &flow, const FlowState &state,
    const FlowLayerValues &exact_trial, const MomentumTimeStencil &stencil,
    double rho_ref_kg_per_m3, double dynamic_viscosity_pa_s,
    const std::vector<double> &owned_active_source_N_per_m3,
    std::array<std::vector<detail::ImmersedFlowWallGradient>, 3>
        wall_pressure_gradients) {
  if (!flow.impl_ || !flow.impl_->domain || !flow.impl_->topology ||
      !flow.impl_->geometry || !flow.impl_->halo || !flow.impl_->mpi ||
      !flow.impl_->execution || !flow.impl_->scratch ||
      !flow.impl_->scratch->storage || !flow.impl_->scratch->access ||
      !(rho_ref_kg_per_m3 > 0.0) || !std::isfinite(rho_ref_kg_per_m3) ||
      !(dynamic_viscosity_pa_s >= 0.0) ||
      !std::isfinite(dynamic_viscosity_pa_s))
    throw runtime::Error(
        "immersed-flow exact momentum residual probe is unavailable");
  const auto expected_stencil = make_momentum_time_stencil(
      stencil.order, stencil.dt_s, stencil.previous_dt_s);
  if (expected_stencil.alpha0 != stencil.alpha0 ||
      expected_stencil.alpha1 != stencil.alpha1 ||
      expected_stencil.alpha2 != stencil.alpha2)
    throw runtime::Error(
        "immersed-flow exact momentum residual stencil is invalid");

  const std::size_t count = flow.impl_->owned_active_cells.size();
  if (count > std::numeric_limits<std::size_t>::max() / 3U ||
      owned_active_source_N_per_m3.size() != count * 3U ||
      std::any_of(owned_active_source_N_per_m3.begin(),
                  owned_active_source_N_per_m3.end(),
                  [](double value) { return !std::isfinite(value); }))
    throw runtime::Error(
        "immersed-flow exact momentum residual source layout is invalid");

  const auto convert_gradients = [&](std::vector<detail::ImmersedFlowWallGradient>
                                         values) {
    std::sort(values.begin(), values.end(),
              [](const auto &left, const auto &right) {
                return left.link < right.link;
              });
    if (values.size() != flow.impl_->wall_links.size() ||
        std::adjacent_find(values.begin(), values.end(),
                           [](const auto &left, const auto &right) {
                             return left.link == right.link;
                           }) != values.end())
      throw runtime::Error(
          "immersed-flow exact momentum wall-gradient layout is invalid");
    std::vector<finite_volume::detail::ImmersedWallNormalGradient> result;
    result.reserve(values.size());
    for (std::size_t index = 0U; index < values.size(); ++index) {
      if (values[index].link != flow.impl_->wall_links[index].link.id ||
          !std::isfinite(values[index].normal_gradient_pa_per_m))
        throw runtime::Error(
            "immersed-flow exact momentum wall-gradient authority is invalid");
      result.push_back(
          {values[index].link, values[index].normal_gradient_pa_per_m});
    }
    return result;
  };
  std::array<std::vector<finite_volume::detail::ImmersedWallNormalGradient>, 3>
      gradients{convert_gradients(std::move(wall_pressure_gradients[0])),
                convert_gradients(std::move(wall_pressure_gradients[1])),
                convert_gradients(std::move(wall_pressure_gradients[2]))};

  const auto box = flow.impl_->topology->owned_global_box();
  const runtime::FieldLayoutSet layout{
      {box.end.x - box.begin.x, box.end.y - box.begin.y,
       box.end.z - box.begin.z},
      flow.impl_->topology->local_face_count()};
  const auto make_probe_state = [&](const FlowLayerValues &values) {
    auto result = FlowState::create(detail::FlowStateSolverAccess::registry(
                                        state),
                                    layout, state.fields(), state.metadata());
    result.seed_accepted_layers(values, values);
    auto &layer = detail::FlowStateSolverAccess::layer(
        result, FlowLayer::committed);
    for (const auto field : {state.fields().density, state.fields().velocity,
                             state.fields().mechanical_pressure})
      flow.impl_->halo->exchange(layer, field);
    return result;
  };
  auto history_state = make_probe_state(state.snapshot(FlowLayer::history));
  auto current_state = make_probe_state(state.snapshot(FlowLayer::committed));
  auto trial_state = make_probe_state(exact_trial);

  auto wall_coefficients = flow.impl_->active_viscous_wall_coefficients();
  ActiveMomentumOperator reference_operator(
      *flow.impl_->topology, *flow.impl_->geometry, *flow.impl_->decomposition,
      flow.impl_->domain->active_cells(), std::move(wall_coefficients),
      *flow.impl_->execution);
  std::vector<double> time_diagonal(count);
  for (std::size_t row = 0U; row < count; ++row)
    time_diagonal[row] =
        stencil.alpha0 * rho_ref_kg_per_m3 *
        flow.impl_->geometry->cell_volume_m3(
            flow.impl_->owned_active_cells[row]) /
        stencil.dt_s;
  reference_operator.replace(time_diagonal, dynamic_viscosity_pa_s);

  struct SpatialTerms final {
    std::vector<double> convective;
    std::vector<double> viscous;
    std::vector<double> pressure;
    std::vector<double> background_pressure;
    std::vector<double> pressure_wall_defect;
    std::vector<double> implicit_reference;
  };
  const auto evaluate_spatial = [&](FlowState &probe_state,
                                    const auto &wall_gradients) {
    auto probe_operator = finite_volume::ImmersedOperatorAdapter::create(
        *flow.impl_->topology, *flow.impl_->geometry, *flow.impl_->domain,
        *flow.impl_->ghost_plan, *flow.impl_->transform,
        *flow.impl_->reconstruction);
    auto &layer = detail::FlowStateSolverAccess::layer(
        probe_state, FlowLayer::committed);
    const auto total = flow.impl_->assemble_immersed_momentum_residual(
        probe_state, layer, dynamic_viscosity_pa_s, &wall_gradients,
        &probe_operator);
    const auto inviscid = flow.impl_->assemble_immersed_momentum_residual(
        probe_state, layer, 0.0, &wall_gradients, &probe_operator);
    const auto pressure = flow.impl_->assemble_immersed_pressure_residual(
        probe_state, layer, &wall_gradients, &probe_operator);
    if (total.size() != count * 3U || inviscid.size() != total.size() ||
        pressure.size() != total.size())
      throw runtime::Error(
          "immersed-flow exact momentum spatial layout is invalid");

    SpatialTerms result;
    result.convective.resize(total.size());
    result.viscous.resize(total.size());
    result.pressure = pressure;
    result.background_pressure = pressure;
    result.pressure_wall_defect.assign(total.size(), 0.0);
    result.implicit_reference.resize(total.size());
#ifdef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS
    const auto boundary_rows = finite_volume::detail::
        ImmersedBoundaryAuthorityAccess::last_boundary_row_evaluations(
            probe_operator);
    for (const auto &boundary_row : boundary_rows) {
      const auto local =
          flow.impl_->topology->find_local_cell(boundary_row.active_cell);
      if (!local.has_value())
        throw runtime::Error(
            "immersed-flow exact momentum pressure row is unavailable");
      const auto active =
          flow.impl_->domain->active_cells().active_index(*local);
      if (!active.has_value() || *active >= count)
        throw runtime::Error(
            "immersed-flow exact momentum pressure row is not owned");
      for (std::size_t component_index = 0U; component_index < 3U;
           ++component_index) {
        const auto offset = *active * 3U + component_index;
        result.background_pressure[offset] =
            boundary_row.residual_before_wall[component_index] +
            boundary_row.background_contribution.pressure[component_index];
        result.pressure_wall_defect[offset] =
            result.pressure[offset] - result.background_pressure[offset];
      }
    }
#else
    throw runtime::Error(
        "immersed-flow exact momentum pressure-row capture is unavailable");
#endif
    const auto &access = detail::FlowStateSolverAccess::access(probe_state);
    const auto velocity = layer.acquire_read<double>(
        access, kStatePhase, kStateActor, probe_state.fields().velocity);
    for (std::size_t component_index = 0U; component_index < 3U;
         ++component_index) {
      execution::Buffer input(*flow.impl_->execution, checked_bytes(count));
      execution::Buffer output(*flow.impl_->execution, checked_bytes(count));
      auto input_view = input.view(0U, count);
      for (std::size_t row = 0U; row < count; ++row) {
        const auto index = flow.impl_->field_index(
            flow.impl_->owned_active_cells[row]);
        input_view[row] = velocity(
            index.i, index.j, index.k, static_cast<int>(component_index));
      }
      reference_operator
          .apply_spatial(static_cast<const execution::Buffer &>(input).view(
                             0U, count),
                         output.view(0U, count))
          .wait();
      const auto output_view =
          static_cast<const execution::Buffer &>(output).view(0U, count);
      for (std::size_t row = 0U; row < count; ++row) {
        const auto offset = row * 3U + component_index;
        result.convective[offset] = inviscid[offset] - pressure[offset];
        result.viscous[offset] = total[offset] - inviscid[offset];
        result.implicit_reference[offset] = output_view[row];
      }
    }
    return result;
  };

  const auto history_terms = evaluate_spatial(history_state, gradients[0]);
  const auto current_terms = evaluate_spatial(current_state, gradients[1]);
  const auto trial_terms = evaluate_spatial(trial_state, gradients[2]);
  const auto history_values = state.snapshot(FlowLayer::history);
  const auto current_values = state.snapshot(FlowLayer::committed);

  std::unordered_set<mesh::GlobalCellId> interface_cells;
  interface_cells.reserve(flow.impl_->domain->links().size());
  for (const auto &link : flow.impl_->domain->links())
    interface_cells.insert(link.fluid_cell);

  struct Wire final {
    std::uint64_t global_cell_id{};
    std::uint8_t interface_row{};
    std::uint8_t padding[7]{};
    double volume{};
    double term[7][3]{};
    double background_pressure[3]{};
    double pressure_wall_defect[3]{};
  };
  static_assert(std::is_trivially_copyable_v<Wire>);
  std::vector<Wire> local(count);
  for (std::size_t row = 0U; row < count; ++row) {
    const auto cell = flow.impl_->owned_active_cells[row];
    const auto gid = flow.impl_->topology->global_cell_id(cell);
    auto &wire = local[row];
    wire.global_cell_id = gid;
    wire.interface_row = interface_cells.count(gid) != 0U ? 1U : 0U;
    wire.volume = flow.impl_->geometry->cell_volume_m3(cell);
    for (std::size_t component_index = 0U; component_index < 3U;
         ++component_index) {
      const auto offset = row * 3U + component_index;
      const auto state_offset = static_cast<std::size_t>(cell) * 3U +
                                component_index;
      const double time = rho_ref_kg_per_m3 * wire.volume / stencil.dt_s *
                          (stencil.alpha0 * exact_trial.velocity[state_offset] +
                           stencil.alpha1 * current_values.velocity[state_offset] +
                           stencil.alpha2 * history_values.velocity[state_offset]);
      const double convective =
          stencil.order == MomentumTimeOrder::backward_euler
              ? current_terms.convective[offset]
              : 2.0 * current_terms.convective[offset] -
                    history_terms.convective[offset];
      const double viscous_remainder =
          stencil.order == MomentumTimeOrder::backward_euler
              ? current_terms.viscous[offset] -
                    current_terms.implicit_reference[offset]
              : 2.0 * (current_terms.viscous[offset] -
                       current_terms.implicit_reference[offset]) -
                    (history_terms.viscous[offset] -
                     history_terms.implicit_reference[offset]);
      const double pressure = trial_terms.pressure[offset];
      const double implicit_reference =
          trial_terms.implicit_reference[offset];
      const double source =
          -wire.volume * owned_active_source_N_per_m3[offset];
      const double total = time + convective + viscous_remainder + pressure +
                           implicit_reference + source;
      wire.term[0][component_index] = time;
      wire.term[1][component_index] = convective;
      wire.term[2][component_index] = viscous_remainder;
      wire.term[3][component_index] = pressure;
      wire.term[4][component_index] = implicit_reference;
      wire.term[5][component_index] = source;
      wire.term[6][component_index] = total;
      wire.background_pressure[component_index] =
          trial_terms.background_pressure[offset];
      wire.pressure_wall_defect[component_index] =
          trial_terms.pressure_wall_defect[offset];
    }
  }

  if (local.size() >
      static_cast<std::size_t>(std::numeric_limits<int>::max()) / sizeof(Wire))
    throw runtime::Error(
        "immersed-flow exact momentum residual payload exceeds MPI int");
  const int local_bytes = static_cast<int>(local.size() * sizeof(Wire));
  std::vector<int> counts(static_cast<std::size_t>(flow.impl_->mpi->size()));
  runtime::detail::check_mpi(
      MPI_Allgather(&local_bytes, 1, MPI_INT, counts.data(), 1, MPI_INT,
                    flow.impl_->mpi->comm()),
      "MPI_Allgather(immersed-flow exact momentum residual counts)");
  std::vector<int> offsets(counts.size(), 0);
  std::size_t total_bytes = 0U;
  for (std::size_t rank = 0U; rank < counts.size(); ++rank) {
    if (counts[rank] < 0 ||
        total_bytes > static_cast<std::size_t>(std::numeric_limits<int>::max()))
      throw runtime::Error(
          "immersed-flow exact momentum residual payload is unsupported");
    offsets[rank] = static_cast<int>(total_bytes);
    total_bytes += static_cast<std::size_t>(counts[rank]);
  }
  if (total_bytes % sizeof(Wire) != 0U ||
      total_bytes > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    throw runtime::Error(
        "immersed-flow exact momentum residual payload size is invalid");
  std::vector<Wire> gathered(total_bytes / sizeof(Wire));
  runtime::detail::check_mpi(
      MPI_Allgatherv(local.data(), local_bytes, MPI_BYTE, gathered.data(),
                     counts.data(), offsets.data(), MPI_BYTE,
                     flow.impl_->mpi->comm()),
      "MPI_Allgatherv(immersed-flow exact momentum residual values)");
  std::sort(gathered.begin(), gathered.end(),
            [](const auto &left, const auto &right) {
              return left.global_cell_id < right.global_cell_id;
            });
  if (gathered.empty() ||
      std::adjacent_find(gathered.begin(), gathered.end(),
                         [](const auto &left, const auto &right) {
                           return left.global_cell_id == right.global_cell_id;
                         }) != gathered.end())
    throw runtime::Error(
        "immersed-flow exact momentum residual global layout is invalid");

  detail::ImmersedFlowExactMomentumResidualData report;
  report.active_global_cell_ids.reserve(gathered.size());
  report.immersed_interface_row.reserve(gathered.size());
  report.cell_volume_m3.reserve(gathered.size());
  std::array<std::vector<double> *, 7> outputs{
      &report.time_residual_N,
      &report.convective_residual_N,
      &report.viscous_remainder_residual_N,
      &report.pressure_residual_N,
      &report.implicit_viscous_reference_residual_N,
      &report.source_residual_N,
      &report.total_residual_N};
  for (auto *values : outputs)
    values->reserve(gathered.size() * 3U);
  report.background_pressure_residual_N.reserve(gathered.size() * 3U);
  report.pressure_wall_defect_residual_N.reserve(gathered.size() * 3U);
  for (const auto &wire : gathered) {
    if ((wire.interface_row != 0U && wire.interface_row != 1U) ||
        !(wire.volume > 0.0) || !std::isfinite(wire.volume))
      throw runtime::Error(
          "immersed-flow exact momentum residual value is invalid");
    report.active_global_cell_ids.push_back(wire.global_cell_id);
    report.immersed_interface_row.push_back(wire.interface_row);
    report.cell_volume_m3.push_back(wire.volume);
    for (std::size_t component_index = 0U; component_index < 3U;
         ++component_index) {
      for (std::size_t term_index = 0U; term_index < outputs.size();
           ++term_index) {
        if (!std::isfinite(wire.term[term_index][component_index]))
          throw runtime::Error(
              "immersed-flow exact momentum residual term is non-finite");
        outputs[term_index]->push_back(
            wire.term[term_index][component_index]);
      }
      if (!std::isfinite(wire.background_pressure[component_index]) ||
          !std::isfinite(wire.pressure_wall_defect[component_index]))
        throw runtime::Error(
            "immersed-flow exact momentum pressure split is non-finite");
      report.background_pressure_residual_N.push_back(
          wire.background_pressure[component_index]);
      report.pressure_wall_defect_residual_N.push_back(
          wire.pressure_wall_defect[component_index]);
      const double reordered =
          (wire.term[0][component_index] + wire.term[5][component_index]) +
          (wire.term[1][component_index] + wire.term[3][component_index]) +
          (wire.term[2][component_index] + wire.term[4][component_index]);
      report.pointwise_split_closure_linf_N = std::max(
          report.pointwise_split_closure_linf_N,
          std::abs(wire.term[6][component_index] - reordered));
    }
  }
  return report;
}

bool detail::ImmersedFlowAccess::
    exact_momentum_residual_report_is_self_consistent(
        const detail::ImmersedFlowExactMomentumResidualData &report) noexcept {
  const std::size_t count = report.active_global_cell_ids.size();
  const std::size_t vector_size = count * 3U;
  if (count == 0U || count > std::numeric_limits<std::size_t>::max() / 3U ||
      report.immersed_interface_row.size() != count ||
      report.cell_volume_m3.size() != count ||
      report.time_residual_N.size() != vector_size ||
      report.convective_residual_N.size() != vector_size ||
      report.viscous_remainder_residual_N.size() != vector_size ||
      report.pressure_residual_N.size() != vector_size ||
      report.background_pressure_residual_N.size() != vector_size ||
      report.pressure_wall_defect_residual_N.size() != vector_size ||
      report.implicit_viscous_reference_residual_N.size() != vector_size ||
      report.source_residual_N.size() != vector_size ||
      report.total_residual_N.size() != vector_size ||
      !std::is_sorted(report.active_global_cell_ids.begin(),
                      report.active_global_cell_ids.end()) ||
      std::adjacent_find(report.active_global_cell_ids.begin(),
                         report.active_global_cell_ids.end()) !=
          report.active_global_cell_ids.end())
    return false;
  double closure_linf = 0.0;
  for (std::size_t row = 0U; row < count; ++row) {
    if ((report.immersed_interface_row[row] != 0U &&
         report.immersed_interface_row[row] != 1U) ||
        !(report.cell_volume_m3[row] > 0.0) ||
        !std::isfinite(report.cell_volume_m3[row]))
      return false;
    for (std::size_t component_index = 0U; component_index < 3U;
         ++component_index) {
      const auto offset = row * 3U + component_index;
      const std::array<double, 7> values{
          report.time_residual_N[offset],
          report.convective_residual_N[offset],
          report.viscous_remainder_residual_N[offset],
          report.pressure_residual_N[offset],
          report.implicit_viscous_reference_residual_N[offset],
          report.source_residual_N[offset], report.total_residual_N[offset]};
      const double background_pressure =
          report.background_pressure_residual_N[offset];
      const double pressure_wall_defect =
          report.pressure_wall_defect_residual_N[offset];
      if (std::any_of(values.begin(), values.end(),
                      [](double value) { return !std::isfinite(value); }) ||
          !std::isfinite(background_pressure) ||
          !std::isfinite(pressure_wall_defect))
        return false;
      const double pressure_split_scale =
          std::max({1.0, std::abs(values[3]), std::abs(background_pressure),
                    std::abs(pressure_wall_defect)});
      if (std::abs(values[3] -
                   (background_pressure + pressure_wall_defect)) >
          512.0 * std::numeric_limits<double>::epsilon() *
              pressure_split_scale)
        return false;
      const double reordered = (values[0] + values[5]) +
                               (values[1] + values[3]) +
                               (values[2] + values[4]);
      const double scale = std::accumulate(
          values.begin(), values.begin() + 6, 1.0,
          [](double result, double value) { return result + std::abs(value); });
      const double closure = std::abs(values[6] - reordered);
      if (closure > 2048.0 * std::numeric_limits<double>::epsilon() * scale)
        return false;
      closure_linf = std::max(closure_linf, closure);
    }
  }
  const double closure_scale =
      std::max({1.0, std::abs(closure_linf),
                std::abs(report.pointwise_split_closure_linf_N)});
  return std::isfinite(report.pointwise_split_closure_linf_N) &&
         std::abs(closure_linf - report.pointwise_split_closure_linf_N) <=
             64.0 * std::numeric_limits<double>::epsilon() * closure_scale;
}

detail::ImmersedFlowPressureCorrectionData
detail::ImmersedFlowAccess::cell_pressure_correction_authority(
    const FixedStepImmersedFlow &flow) {
  if (!flow.impl_ ||
      flow.impl_->pending_cell_pressure_correction_records.size() != 2U)
    throw runtime::Error(
        "immersed-flow cell pressure-correction capture is unavailable");
  return {flow.impl_->pending_cell_pressure_correction_records,
          flow.impl_->pending_inter_corrector_authority_difference_l2};
}

std::vector<detail::ImmersedFlowOutletPredictorValue>
detail::ImmersedFlowAccess::physical_outlet_predictor_values(
    FixedStepImmersedFlow &flow, FlowState &state, double rho_ref_kg_per_m3,
    const MomentumTimeStencil &stencil, double momentum_diagonal) {
  if (!flow.impl_ || !flow.impl_->pressure_operator || !flow.impl_->domain ||
      !flow.impl_->topology || !flow.impl_->geometry ||
      !flow.impl_->boundaries || !(momentum_diagonal > 0.0) ||
      !std::isfinite(momentum_diagonal))
    throw runtime::Error("immersed-flow outlet predictor probe is unavailable");
  const std::size_t count =
      flow.impl_->domain->active_cells().owned_active_count();
  std::array<std::vector<double>, 3> diagonal{
      std::vector<double>(count, momentum_diagonal),
      std::vector<double>(count, momentum_diagonal),
      std::vector<double>(count, momentum_diagonal)};
  flow.impl_->pressure_operator->prepare_face_coefficients(rho_ref_kg_per_m3,
                                                           diagonal);
  static_cast<void>(flow.impl_->assemble_face_predictor(
      state, rho_ref_kg_per_m3, stencil, diagonal));
  const auto trial = state.snapshot(FlowLayer::trial);
  std::vector<detail::ImmersedFlowOutletPredictorValue> result;
  for (mesh::LocalFaceId face = 0U;
       face < flow.impl_->topology->local_face_count(); ++face) {
    if (flow.impl_->topology->neighbour(face).has_value())
      continue;
    const auto patch = flow.impl_->topology->patch_id(face);
    const auto active = flow.impl_->domain->active_cells().active_index(
        flow.impl_->topology->owner(face));
    if (!patch.has_value() || !active.has_value() || *active >= count ||
        flow.impl_->boundaries->patch(*patch).pressure_rule() !=
            boundary::PressureRule::prescribed_value)
      continue;
    result.push_back(
        {flow.impl_->topology->global_face_id(face),
         {trial.face_velocity[static_cast<std::size_t>(face) * 3U],
          trial.face_velocity[static_cast<std::size_t>(face) * 3U + 1U],
          trial.face_velocity[static_cast<std::size_t>(face) * 3U + 2U]},
         trial.face_mass_flux[face]});
  }
  std::sort(result.begin(), result.end(),
            [](const auto &left, const auto &right) {
              return left.global_face_id < right.global_face_id;
            });
  return result;
}

std::vector<detail::ImmersedFlowWallPredictorValue>
detail::ImmersedFlowAccess::wall_predictor_values(
    FixedStepImmersedFlow &flow, FlowState &state, double rho_ref_kg_per_m3,
    const MomentumTimeStencil &stencil, double momentum_diagonal) {
  if (!flow.impl_ || !flow.impl_->pressure_operator || !flow.impl_->domain ||
      !flow.impl_->halo || !(momentum_diagonal > 0.0) ||
      !std::isfinite(momentum_diagonal))
    throw runtime::Error("immersed-flow wall predictor probe is unavailable");
  const auto fields = state.fields();
  for (auto layer : {FlowLayer::committed, FlowLayer::history}) {
    auto &storage = detail::FlowStateSolverAccess::layer(state, layer);
    flow.impl_->halo->exchange(storage, fields.velocity);
    flow.impl_->halo->exchange(storage, fields.mechanical_pressure);
  }
  std::optional<detail::ImmersedDensityAttemptAuthority> trial_density;
  std::optional<detail::ImmersedDensityAttemptAuthority> committed_density;
  std::optional<detail::ImmersedDensityAttemptAuthority> history_density;
  if (flow.impl_->density_setup.model == config::DensityModel::material) {
    for (auto layer : {FlowLayer::trial, FlowLayer::committed,
                       FlowLayer::history}) {
      auto &storage = detail::FlowStateSolverAccess::layer(state, layer);
      flow.impl_->halo->exchange(storage, fields.density);
    }
    const auto density_stage =
        flow.impl_->prepare_density_authority_collectively([&] {
          trial_density.emplace(flow.impl_->density_authority(
              state, detail::FlowStateSolverAccess::layer(state,
                                                          FlowLayer::trial)));
          committed_density.emplace(flow.impl_->density_authority(
              state, detail::FlowStateSolverAccess::layer(
                         state, FlowLayer::committed)));
          history_density.emplace(flow.impl_->density_authority(
              state, detail::FlowStateSolverAccess::layer(
                         state, FlowLayer::history)));
        });
    if (density_stage.reason != MaterialTransportFailureReason::none)
      throw runtime::Error(
          "immersed-flow wall predictor density probe is invalid");
  }
  const std::size_t count =
      flow.impl_->domain->active_cells().owned_active_count();
  std::array<std::vector<double>, 3> diagonal{
      std::vector<double>(count, momentum_diagonal),
      std::vector<double>(count, momentum_diagonal),
      std::vector<double>(count, momentum_diagonal)};
  if (trial_density.has_value())
    flow.impl_->pressure_operator->prepare_face_coefficients(
        trial_density->face_density, diagonal);
  else
    flow.impl_->pressure_operator->prepare_face_coefficients(rho_ref_kg_per_m3,
                                                             diagonal);
  const auto conditions = flow.impl_->assemble_face_predictor(
      state, rho_ref_kg_per_m3, stencil, diagonal, nullptr,
      trial_density ? &*trial_density : nullptr,
      committed_density ? &*committed_density : nullptr,
      history_density ? &*history_density : nullptr);
  std::vector<detail::ImmersedFlowWallPredictorValue> result;
  result.reserve(conditions.size());
  for (std::size_t index = 0U; index < conditions.size(); ++index) {
    const auto &condition = conditions[index];
    const double wall_density =
        trial_density ? trial_density->wall_density[index].rho_wall_kg_per_m3
                      : rho_ref_kg_per_m3;
    const double wall_density_normal_derivative =
        trial_density
            ? trial_density->wall_density[index]
                  .normal_derivative_kg_per_m4
            : 0.0;
    result.push_back({condition.link, condition.predictor_mass_flux_kg_per_s,
                      wall_density, wall_density_normal_derivative});
  }
  std::sort(result.begin(), result.end(),
            [](const auto &left, const auto &right) {
              return left.link < right.link;
            });
  return result;
}

detail::ImmersedFlowFinalPressureResidualData
detail::ImmersedFlowAccess::final_momentum_pressure_residual_routes(
    FixedStepImmersedFlow &flow, FlowState &state) {
  if (!flow.impl_ || !flow.impl_->domain || !flow.impl_->topology ||
      !flow.impl_->mpi ||
      flow.impl_->committed_final_momentum_pressure_residual_step !=
          state.metadata().step)
    throw runtime::Error(
        "immersed-flow final momentum pressure residual capture is unavailable");
  const std::size_t local_count = flow.impl_->owned_active_cells.size();
  if (local_count > std::numeric_limits<std::size_t>::max() / 3U ||
      flow.impl_->committed_final_momentum_pressure_residual.size() !=
          local_count * 3U)
    throw runtime::Error(
        "immersed-flow final momentum pressure residual capture is invalid");

  auto &committed =
      detail::FlowStateSolverAccess::layer(state, FlowLayer::committed);
  const auto direct = flow.impl_->assemble_immersed_pressure_residual(
      state, committed, &flow.impl_->committed_final_pressure_wall_gradients);
  if (direct.size() != local_count * 3U)
    throw runtime::Error(
        "immersed-flow final momentum pressure residual reassembly is invalid");

  struct ResidualWire final {
    std::uint64_t global_cell_id{};
    double used[3]{};
    double direct[3]{};
  };
  static_assert(std::is_trivially_copyable_v<ResidualWire>);
  std::vector<ResidualWire> local(local_count);
  for (std::size_t row = 0U; row < local_count; ++row) {
    local[row].global_cell_id = flow.impl_->topology->global_cell_id(
        flow.impl_->owned_active_cells[row]);
    for (std::size_t component = 0U; component < 3U; ++component) {
      const auto offset = row * 3U + component;
      local[row].used[component] =
          flow.impl_->committed_final_momentum_pressure_residual[offset];
      local[row].direct[component] = direct[offset];
    }
  }
  if (local.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()) /
                         sizeof(ResidualWire))
    throw runtime::Error(
        "immersed-flow final momentum pressure residual payload exceeds MPI int");
  const int local_bytes = static_cast<int>(local.size() * sizeof(ResidualWire));
  std::vector<int> counts(static_cast<std::size_t>(flow.impl_->mpi->size()));
  runtime::detail::check_mpi(
      MPI_Allgather(&local_bytes, 1, MPI_INT, counts.data(), 1, MPI_INT,
                    flow.impl_->mpi->comm()),
      "MPI_Allgather(immersed-flow final pressure residual counts)");
  std::vector<int> offsets(counts.size(), 0);
  std::size_t total_bytes = 0U;
  for (std::size_t rank = 0U; rank < counts.size(); ++rank) {
    if (counts[rank] < 0 ||
        total_bytes > static_cast<std::size_t>(std::numeric_limits<int>::max()))
      throw runtime::Error(
          "immersed-flow final momentum pressure residual payload is unsupported");
    offsets[rank] = static_cast<int>(total_bytes);
    total_bytes += static_cast<std::size_t>(counts[rank]);
  }
  if (total_bytes % sizeof(ResidualWire) != 0U ||
      total_bytes > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    throw runtime::Error(
        "immersed-flow final momentum pressure residual payload size is invalid");
  std::vector<ResidualWire> gathered(total_bytes / sizeof(ResidualWire));
  runtime::detail::check_mpi(
      MPI_Allgatherv(local.data(), local_bytes, MPI_BYTE, gathered.data(),
                     counts.data(), offsets.data(), MPI_BYTE,
                     flow.impl_->mpi->comm()),
      "MPI_Allgatherv(immersed-flow final pressure residual values)");
  std::sort(gathered.begin(), gathered.end(),
            [](const auto &left, const auto &right) {
              return left.global_cell_id < right.global_cell_id;
            });
  if (gathered.empty() ||
      std::adjacent_find(gathered.begin(), gathered.end(),
                         [](const auto &left, const auto &right) {
                           return left.global_cell_id == right.global_cell_id;
                         }) != gathered.end())
    throw runtime::Error(
        "immersed-flow final momentum pressure residual global layout is invalid");

  detail::ImmersedFlowFinalPressureResidualData report;
  report.active_global_cell_ids.reserve(gathered.size());
  report.used_residual.reserve(gathered.size() * 3U);
  report.direct_residual.reserve(gathered.size() * 3U);
  double maximum_difference = -1.0;
  double used_square = 0.0;
  double direct_square = 0.0;
  double difference_square = 0.0;
  for (const auto &cell : gathered) {
    report.active_global_cell_ids.push_back(cell.global_cell_id);
    for (std::size_t component = 0U; component < 3U; ++component) {
      const double used = cell.used[component];
      const double reassembled = cell.direct[component];
      if (!std::isfinite(used) || !std::isfinite(reassembled))
        throw runtime::Error(
            "immersed-flow final momentum pressure residual is non-finite");
      const double difference = used - reassembled;
      const double absolute_difference = std::abs(difference);
      report.used_residual.push_back(used);
      report.direct_residual.push_back(reassembled);
      used_square += used * used;
      direct_square += reassembled * reassembled;
      difference_square += difference * difference;
      report.used_linf = std::max(report.used_linf, std::abs(used));
      report.direct_linf = std::max(report.direct_linf, std::abs(reassembled));
      report.difference_linf =
          std::max(report.difference_linf, absolute_difference);
      if (absolute_difference > maximum_difference) {
        maximum_difference = absolute_difference;
        report.maximum_difference_global_cell_id = cell.global_cell_id;
        report.maximum_difference_component =
            static_cast<std::uint32_t>(component);
        report.maximum_difference_used_value = used;
        report.maximum_difference_direct_value = reassembled;
      }
    }
  }
  report.used_l2 = std::sqrt(used_square);
  report.direct_l2 = std::sqrt(direct_square);
  report.difference_l2 = std::sqrt(difference_square);
  if (!std::isfinite(report.used_l2) || !std::isfinite(report.direct_l2) ||
      !std::isfinite(report.difference_l2))
    throw runtime::Error(
        "immersed-flow final momentum pressure residual norms are non-finite");
  return report;
}

double detail::ImmersedFlowAccess::inter_step_pressure_authority_difference_l2(
    const FixedStepImmersedFlow &flow) noexcept {
  return flow.impl_ == nullptr
             ? std::numeric_limits<double>::quiet_NaN()
             : flow.impl_->pending_inter_step_authority_difference_l2;
}

double
detail::ImmersedFlowAccess::bdf2_history_pressure_authority_difference_l2(
    const FixedStepImmersedFlow &flow) noexcept {
  return flow.impl_ == nullptr
             ? std::numeric_limits<double>::quiet_NaN()
             : flow.impl_->pending_bdf2_history_authority_difference_l2;
}

std::array<std::uint64_t, 2>
detail::ImmersedFlowAccess::pressure_authority_fingerprints(
    const FixedStepImmersedFlow &flow) noexcept {
  if (flow.impl_ == nullptr)
    return {};
  const auto fingerprint = [](bool available, const auto &gradients) {
    std::uint64_t hash = hash_u64(kFnvOffset, available ? 1U : 0U);
    hash = hash_u64(
        hash, available ? static_cast<std::uint64_t>(gradients.size()) : 0U);
    if (!available)
      return hash;
    for (const auto &gradient : gradients) {
      hash = hash_u64(hash, gradient.link);
      hash = hash_u64(hash, fp64_bits(gradient.value));
    }
    return hash;
  };
  return {fingerprint(flow.impl_->history_pressure_authority_available,
                      flow.impl_->history_pressure_wall_gradients),
          fingerprint(flow.impl_->committed_pressure_authority_available,
                      flow.impl_->committed_pressure_wall_gradients)};
}

detail::ImmersedFlowPressureAuthorityData
detail::ImmersedFlowAccess::pressure_authority_state(
    const FixedStepImmersedFlow &flow) {
  if (flow.impl_ == nullptr)
    return {};
  const auto snapshot = [](const auto &gradients) {
    std::vector<detail::ImmersedFlowWallGradient> result;
    result.reserve(gradients.size());
    for (const auto &gradient : gradients)
      result.push_back({gradient.link, gradient.value});
    return result;
  };
  return {snapshot(flow.impl_->history_pressure_wall_gradients),
          snapshot(flow.impl_->committed_pressure_wall_gradients),
          snapshot(flow.impl_->pending_pressure_wall_gradients),
          flow.impl_->history_pressure_authority_available,
          flow.impl_->committed_pressure_authority_available};
}

#endif

} // namespace hundun::flow
