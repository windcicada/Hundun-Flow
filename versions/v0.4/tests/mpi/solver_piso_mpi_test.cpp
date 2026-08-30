// SPDX-License-Identifier: Apache-2.0

#include "../support/ibm_force_fixture.hpp"
#include "../support/candidate_boundary_fixture.hpp"
#include "../support/piso_fixture.hpp"
#include "parallel_halo_detail.hpp"
#include "solver_piso_detail.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string_view>
#include <utility>
#include <vector>

namespace {

using namespace hundun::v04;
using namespace hundun::v04::test;

bool expect(bool condition, int rank, std::string_view description) {
  if (!condition) {
    std::cerr << "rank " << rank << " FAIL: " << description << '\n';
  }
  return condition;
}

bool all_true(bool local, MPI_Comm communicator) {
  const int input = local ? 1 : 0;
  int result = 0;
  return MPI_Allreduce(&input, &result, 1, MPI_INT, MPI_MIN, communicator) ==
             MPI_SUCCESS &&
         result != 0;
}

std::vector<double> face_values(ConstFaceFluxView flux) {
  std::vector<double> values;
  const std::array<ConstFaceFieldView, 3U> faces{flux.x, flux.y, flux.z};
  for (ConstFaceFieldView face : faces) {
    for (std::int32_t z = 0; z < face.extents.z; ++z)
      for (std::int32_t y = 0; y < face.extents.y; ++y)
        for (std::int32_t x = 0; x < face.extents.x; ++x)
          values.push_back(face.unchecked({x, y, z}));
  }
  return values;
}

double periodic_value(std::int32_t index, std::int32_t count, double phase) {
  constexpr double two_pi = 6.283185307179586476925286766559;
  return std::sin(two_pi * (static_cast<double>(index) + 0.5) /
                      static_cast<double>(count) +
                  phase);
}

double cell_volume(const PeriodicPisoFixture& fixture, Int3 local) {
  const Int3 global{fixture.patch.begin.x + local.x,
                    fixture.patch.begin.y + local.y,
                    fixture.patch.begin.z + local.z};
  return fixture.geometry.x().widths().data[global.x] *
         fixture.geometry.y().widths().data[global.y] *
         fixture.geometry.z().widths().data[global.z];
}

double time_step_for_bdf(BdfCoefficients bdf) {
  if (bdf.order == 1U) return 1.0 / bdf.a0;
  const double ratio = 1.0 / (std::sqrt(-bdf.a1 / bdf.a2) - 1.0);
  return (1.0 + ratio) / -bdf.a1;
}

struct LinearResources {
  OwnedField vectors;
  OwnedField scalars;
  SolverWorkspace workspace;
  ReductionEngine reductions;
  HaloEngine operator_halo;
  PressureLinearOperator linear_operator;
  MgWorkspaceRequirements mg_requirements;
  OwnedField mg_vectors;
  MgWorkspace mg_workspace;
  HaloEngine mg_halo;
  std::vector<HaloEngine> coarse_halos;
  std::vector<HaloEngine*> coarse_halo_pointers;
  NativeCartesianMgPlan mg;
  MgPlanCounters mg_counters;
  NativeCartesianMgSpec mg_spec;
};

class NthHaloFailureOperator final : public LinearOperator {
 public:
  NthHaloFailureOperator(LinearOperator& delegated,
                         detail::HaloFailurePoint point, int failing_rank,
                         std::uint32_t failing_apply) noexcept
      : delegated_(&delegated),
        point_(point),
        failing_rank_(failing_rank),
        failing_apply_(failing_apply) {}

  LinearOperatorCertificate certificate() const noexcept override {
    return delegated_->certificate();
  }

  Status apply(FieldView input, FieldView output) const noexcept override {
    ++calls_;
    if (calls_ == failing_apply_) {
      detail::set_halo_failure_for_test(point_, failing_rank_);
    }
    return delegated_->apply(input, output);
  }

  LinearOperatorFailureProvenance failure_provenance() const
      noexcept override {
    return delegated_->failure_provenance();
  }

  std::uint32_t calls() const noexcept { return calls_; }

 private:
  LinearOperator* delegated_{};
  detail::HaloFailurePoint point_{detail::HaloFailurePoint::none};
  int failing_rank_{-1};
  std::uint32_t failing_apply_{};
  mutable std::uint32_t calls_{};
};

class NthPostApplyFailureOperator final : public LinearOperator {
 public:
  NthPostApplyFailureOperator(LinearOperator& delegated, int rank,
                              int failing_rank,
                              std::uint32_t failing_apply) noexcept
      : delegated_(&delegated),
        rank_(rank),
        failing_rank_(failing_rank),
        failing_apply_(failing_apply) {}

  LinearOperatorCertificate certificate() const noexcept override {
    return delegated_->certificate();
  }

  Status apply(FieldView input, FieldView output) const noexcept override {
    ++calls_;
    const Status delegated = delegated_->apply(input, output);
    if (!delegated) return delegated;
    return rank_ == failing_rank_ && calls_ == failing_apply_
               ? Status{StatusCode::numerical_failure, 8520U}
               : Status{};
  }

 private:
  LinearOperator* delegated_{};
  int rank_{};
  int failing_rank_{-1};
  std::uint32_t failing_apply_{};
  mutable std::uint32_t calls_{};
};

bool exact_halo_delta(HaloRuntimeCounters before, HaloRuntimeCounters after,
                      HaloPlanStats plan,
                      std::uint64_t exchanges) noexcept {
  return after.begin_calls == before.begin_calls + exchanges &&
         after.finish_calls == before.finish_calls + exchanges &&
         after.messages_started ==
             before.messages_started +
                 exchanges * plan.maximum_messages_per_exchange &&
         after.bytes_packed ==
             before.bytes_packed +
                 exchanges * plan.maximum_bytes_per_exchange &&
         after.bytes_unpacked ==
             before.bytes_unpacked +
                 exchanges * plan.receive_capacity_doubles * sizeof(double);
}

bool initialize_linear_resources(
    MPI_Comm world, PeriodicPisoFixture& fixture,
    const PressureCorrectionCertificate& pressure,
    PressureCorrectionSystemView system, LinearIdentity& identity,
    MgCoefficientIdentity coefficients, LinearResources& resources) {
  LinearWorkspaceRequirements requirements;
  if (!make_linear_workspace_requirements(
          LinearAlgorithm::fgmres, fixture.patch.cells, 1U,
          fixture.piso.pressure_solve().restart,
          ReductionMode::mpi_allreduce, 7001U, requirements)) {
    return false;
  }
  constexpr StorageIdentity krylov_storage = 7101U;
  resources.vectors = make_field(70U, fixture.patch.cells,
                                 requirements.vector_slots, 1U, 7002U,
                                 krylov_storage);
  resources.scalars = make_field(
      71U,
      {static_cast<std::int32_t>(requirements.scalar_doubles), 1, 1},
      1U, 0U, 7003U, krylov_storage);
  if (!SolverWorkspace::bind(requirements, resources.vectors.view,
                             resources.scalars.view, resources.workspace) ||
      !ReductionEngine::compile(world, ReductionMode::mpi_allreduce,
                                requirements.reduction_capacity,
                                resources.reductions)) {
    return false;
  }
  identity.workspace = resources.workspace.fingerprint();
  const std::array<HaloFieldSpec, 1U> operator_fields{{
      {resources.vectors.view.field, 1U, 1U}}};
  if (!resources.operator_halo.reserve(
          world, fixture.patch,
          {operator_fields.data(), operator_fields.size()},
          fixture.boundary.halo_topology()) ||
      !fixture.coupler.bind_pressure_operator(
          {world, &resources.operator_halo, 7004U,
           resources.vectors.view.field},
          system, resources.linear_operator)) {
    return false;
  }

  resources.mg_spec.communicator = world;
  resources.mg_spec.geometry = &fixture.geometry;
  resources.mg_spec.patch = fixture.patch;
  resources.mg_spec.boundaries = {
      MgBoundaryKind::periodic, MgBoundaryKind::periodic,
      MgBoundaryKind::periodic, MgBoundaryKind::periodic,
      MgBoundaryKind::periodic, MgBoundaryKind::periodic};
  resources.mg_spec.null_space = MgNullSpace::none;
  resources.mg_spec.identity = identity;
  resources.mg_spec.coefficients = coefficients;
  if (!make_mg_workspace_requirements(
          world, fixture.geometry, fixture.patch, resources.mg_spec.policy,
          7005U, resources.mg_requirements)) {
    return false;
  }
  constexpr StorageIdentity mg_storage = 7102U;
  resources.mg_vectors = make_field(
      80U, resources.mg_requirements.arena_shape, 1U, 1U, 7006U,
      mg_storage);
  if (!MgWorkspace::bind(resources.mg_requirements,
                         resources.mg_vectors.view,
                         resources.mg_workspace)) {
    return false;
  }
  const std::array<HaloFieldSpec, 1U> mg_fields{{
      {resources.mg_vectors.view.field, 1U, 1U}}};
  if (!resources.mg_halo.reserve(
          world, resources.mg_requirements.levels[0U].patch,
          {mg_fields.data(), mg_fields.size()},
          fixture.boundary.halo_topology())) {
    return false;
  }
  resources.coarse_halos.resize(resources.mg_requirements.level_count - 1U);
  resources.coarse_halo_pointers.resize(resources.coarse_halos.size());
  for (std::size_t level = 1U;
       level < resources.mg_requirements.level_count; ++level) {
    if (!resources.coarse_halos[level - 1U].reserve(
            world, resources.mg_requirements.levels[level].patch,
            {mg_fields.data(), mg_fields.size()},
            fixture.boundary.halo_topology())) {
      return false;
    }
    resources.coarse_halo_pointers[level - 1U] =
        &resources.coarse_halos[level - 1U];
  }
  const MgRuntimeServices services{
      &resources.mg_halo, &resources.reductions, &resources.mg_workspace,
      {resources.coarse_halo_pointers.data(),
       resources.coarse_halo_pointers.size()}};
  return static_cast<bool>(fixture.coupler.compile_native_pressure_mg(
      pressure, resources.mg_spec, services, system, resources.mg,
      &resources.mg_counters));
}

struct OwnedCoefficientFace {
  std::vector<double> storage;
  FaceFieldView view{};
};

OwnedCoefficientFace make_coefficient_face(CartesianAxis axis, Int3 cells,
                                            StorageIdentity identity) {
  OwnedCoefficientFace result;
  Int3 extents = cells;
  if (axis == CartesianAxis::x) ++extents.x;
  if (axis == CartesianAxis::y) ++extents.y;
  if (axis == CartesianAxis::z) ++extents.z;
  const std::size_t stride_y = static_cast<std::size_t>(extents.x);
  const std::size_t stride_z = stride_y * extents.y;
  result.storage.assign(stride_z * extents.z, 1.0);
  result.view = {result.storage.data(), extents, stride_y, stride_z, axis,
                 identity, 8501U};
  return result;
}

Status seed_cycle_projection(
    Int3 cells, LinearIdentity identity,
    const std::array<std::vector<double>,
                     kLinearRecycleMaximumDirections>& cycle_directions,
    std::size_t cycle_direction_count, OwnedField& caller,
    SolverWorkspace& workspace, ReductionEngine& reductions) {
  fill(caller, 0.0);
  Status seeded =
      workspace.recycle_begin_capture_for_test(cells, identity.fingerprint);
  for (std::size_t direction = 0U;
       direction < cycle_direction_count && seeded; ++direction) {
    if (cycle_directions[direction].size() !=
        static_cast<std::size_t>(cells.x * cells.y * cells.z)) {
      return {StatusCode::invalid_plan, 8510U};
    }
    seeded = workspace.recycle_capture_cycle_start_for_test(
        as_const(caller.view), reductions);
    for (std::int32_t z = 0; z < cells.z && seeded; ++z) {
      for (std::int32_t y = 0; y < cells.y; ++y) {
        for (std::int32_t x = 0; x < cells.x; ++x) {
          const std::size_t local =
              (static_cast<std::size_t>(z) * cells.y + y) * cells.x + x;
          caller.view.unchecked({x, y, z}, 0U) +=
              cycle_directions[direction][local];
        }
      }
    }
    ++caller.view.revision;
    if (seeded) {
      seeded = workspace.recycle_capture_cycle_publish_for_test(
          as_const(caller.view), reductions);
    }
  }
  fill(caller, 0.0);
  if (seeded) {
    seeded = workspace.recycle_begin_projection_for_test(
        cells, identity.fingerprint);
  }
  return seeded;
}

bool test_real_pressure_halo_failure_provenance(
    MPI_Comm world, int rank, Int3 cells, ConstFieldView rhs,
    LinearIdentity identity, LinearSolveControl control,
    LinearOperator& exact_operator, LinearResources& linear,
    const std::array<std::vector<double>,
                     kLinearRecycleMaximumDirections>& cycle_directions,
    std::size_t cycle_direction_count, bool full_solve_retry) {
  int size = 0;
  MPI_Comm_size(world, &size);
  const int failing_rank = size - 1;
  const std::array failure_cases{
      std::pair{detail::HaloFailurePoint::start,
                detail::halo_detail_start_failure},
      std::pair{detail::HaloFailurePoint::completion,
                detail::halo_detail_completion_failure}};
  bool passed = true;

  for (const auto [failure_point, expected_detail] : failure_cases) {
    OwnedField caller =
        make_field(91U, cells, 1U, 0U, 8502U, 8503U);
    const auto seed_projection = [&]() {
      return static_cast<bool>(seed_cycle_projection(
          cells, identity, cycle_directions, cycle_direction_count, caller,
          linear.workspace, linear.reductions));
    };
    bool projection_seeded = seed_projection();
    if (full_solve_retry &&
        failure_point == detail::HaloFailurePoint::start) {
      const LinearSolveInvocation clean_invocation{
          rhs, caller.view, identity, control, nullptr};
      const LinearSolveResult clean = solve_fgmres(
          exact_operator, linear.mg, clean_invocation, linear.workspace,
          linear.reductions);
      passed &= expect(
          projection_seeded && clean.status.code == StatusCode::ok &&
              clean.termination == LinearTermination::converged,
          rank,
          "reconstructed C1 directions reproduce a clean C2 solve before fault injection");
      projection_seeded = seed_projection();
    }
    const std::vector<double> caller_before = caller.storage;

    NthHaloFailureOperator fault_operator(exact_operator, failure_point,
                                          failing_rank, 2U);
    const LinearSolveInvocation invocation{
        rhs, caller.view, identity, control, nullptr};
    const LinearSolveResult failed = solve_fgmres(
        fault_operator, linear.mg, invocation, linear.workspace,
        linear.reductions);
    const LinearOperatorFailureProvenance provenance =
        fault_operator.failure_provenance();
    const int halo_lowest_failing_rank =
        linear.operator_halo.lowest_failing_rank();
    detail::clear_halo_failure_for_test();

    const std::array<std::uint64_t, 2U> local_status{
        static_cast<std::uint64_t>(failed.status.code),
        failed.status.detail};
    std::array<std::uint64_t, 2U> minimum_status{};
    std::array<std::uint64_t, 2U> maximum_status{};
    const int minimum_result = MPI_Allreduce(
        local_status.data(), minimum_status.data(),
        static_cast<int>(local_status.size()), MPI_UINT64_T, MPI_MIN, world);
    const int maximum_result = MPI_Allreduce(
        local_status.data(), maximum_status.data(),
        static_cast<int>(local_status.size()), MPI_UINT64_T, MPI_MAX, world);

    passed &= expect(
        minimum_result == MPI_SUCCESS && maximum_result == MPI_SUCCESS &&
            minimum_status == maximum_status &&
            failed.status.code == StatusCode::mpi_failure &&
            failed.status.detail == expected_detail &&
            failed.termination == LinearTermination::operator_failure &&
            failed.lowest_failing_rank == failing_rank &&
            failed.operator_applies == 1U &&
            projection_seeded &&
            failed.recycle_offered_directions == cycle_direction_count &&
            failed.recycle_projection_attempted &&
            !failed.recycle_projection_accepted &&
            failed.recycle_operator_applies == 1U &&
            caller.storage == caller_before,
        rank,
        "real current-A2 pressure Halo fault reaches the solver gate with true provenance and no caller publication");
    passed &= expect(
        provenance.status.code == failed.status.code &&
            provenance.status.detail == failed.status.detail &&
            provenance.status_scope ==
                LinearOperatorStatusScope::collective &&
            provenance.lowest_failing_rank == failing_rank &&
            halo_lowest_failing_rank == failing_rank &&
            !linear.operator_halo.active() && linear.operator_halo.ready(),
        rank,
        "pressure adapter snapshots exact collective Halo provenance");

    const Status stale_projection =
        linear.workspace.recycle_begin_projection_for_test(
            cells, identity.fingerprint);
    const bool retry_seeded = full_solve_retry ? seed_projection() : true;
    const HaloRuntimeCounters retry_halo_before =
        linear.operator_halo.runtime_counters();
    LinearSolveResult retry;
    Status direct_retry{};
    if (full_solve_retry) {
      retry = solve_fgmres(exact_operator, linear.mg, invocation,
                           linear.workspace, linear.reductions);
    } else {
      FieldView input = linear.workspace.vector(0U, cells);
      FieldView output = linear.workspace.vector(1U, cells);
      direct_retry = exact_operator.apply(input, output);
    }
    const HaloRuntimeCounters retry_halo_after =
        linear.operator_halo.runtime_counters();
    const LinearOperatorFailureProvenance cleared =
        exact_operator.failure_provenance();
    const std::uint64_t retry_operator_exchanges =
        full_solve_retry
            ? retry.operator_applies + retry.recycle_operator_applies
            : 1U;
    const bool retry_matrix =
        stale_projection.code == StatusCode::invalid_plan && retry_seeded &&
            (full_solve_retry
                 ? retry.status.code == StatusCode::ok &&
                       retry.termination == LinearTermination::converged
                 : static_cast<bool>(direct_retry)) &&
            exact_halo_delta(retry_halo_before, retry_halo_after,
                             linear.operator_halo.plan_stats(),
                             retry_operator_exchanges) &&
            cleared.status &&
            cleared.status_scope ==
                LinearOperatorStatusScope::rank_local &&
            cleared.lowest_failing_rank == -1 &&
            linear.operator_halo.lowest_failing_rank() == -1 &&
            !linear.operator_halo.active() && linear.operator_halo.ready();
    if (!retry_matrix) {
      const HaloPlanStats plan = linear.operator_halo.plan_stats();
      std::cerr
          << "rank " << rank << " C2 Halo retry status="
          << static_cast<unsigned>(full_solve_retry ? retry.status.code
                                                    : direct_retry.code)
          << " detail="
          << (full_solve_retry ? retry.status.detail : direct_retry.detail)
          << " term=" << static_cast<unsigned>(retry.termination)
          << " stale/seed="
          << static_cast<unsigned>(stale_projection.code) << '/'
          << retry_seeded << " ordinary/recycle applies="
          << retry.operator_applies << '/' << retry.recycle_operator_applies
          << " offered/retained/accepted="
          << retry.recycle_offered_directions << '/'
          << retry.recycle_retained_directions << '/'
          << retry.recycle_projection_accepted
          << " initial/projected/final=" << retry.initial_true_residual << '/'
          << retry.recycle_projected_true_residual << '/'
          << retry.final_true_residual
          << " halo begin/finish/messages/pack/unpack delta="
          << (retry_halo_after.begin_calls - retry_halo_before.begin_calls)
          << '/'
          << (retry_halo_after.finish_calls - retry_halo_before.finish_calls)
          << '/'
          << (retry_halo_after.messages_started -
              retry_halo_before.messages_started)
          << '/'
          << (retry_halo_after.bytes_packed - retry_halo_before.bytes_packed)
          << '/'
          << (retry_halo_after.bytes_unpacked -
              retry_halo_before.bytes_unpacked)
          << " per=" << plan.maximum_messages_per_exchange << '/'
          << plan.maximum_bytes_per_exchange << '/'
          << plan.receive_capacity_doubles * sizeof(double) << '\n';
    }
    passed &= expect(
        retry_matrix, rank,
        "consumed C2 Halo fault immediately retries the same pressure/Halo resources with exact messages and bytes");
  }
  return all_true(passed, world);
}

bool test_real_pressure_projection_failure_lifecycle(
    MPI_Comm world, int rank, Int3 cells, ConstFieldView rhs,
    LinearIdentity identity, LinearSolveControl control,
    PressureLinearOperator& exact_operator, LinearResources& linear,
    const std::array<std::vector<double>,
                     kLinearRecycleMaximumDirections>& cycle_directions,
    std::size_t cycle_direction_count) {
  int size = 0;
  MPI_Comm_size(world, &size);
  const int failing_rank = size - 1;
  bool passed = true;

  const auto retry_projection = [&](OwnedField& caller) {
    const bool seeded = static_cast<bool>(seed_cycle_projection(
        cells, identity, cycle_directions, cycle_direction_count, caller,
        linear.workspace, linear.reductions));
    const HaloRuntimeCounters before =
        linear.operator_halo.runtime_counters();
    const LinearSolveInvocation invocation{
        rhs, caller.view, identity, control, nullptr};
    const LinearSolveResult retry = solve_fgmres(
        exact_operator, linear.mg, invocation, linear.workspace,
        linear.reductions);
    const HaloRuntimeCounters after =
        linear.operator_halo.runtime_counters();
    return seeded && retry.status.code == StatusCode::ok &&
           retry.termination == LinearTermination::converged &&
           exact_halo_delta(before, after, linear.operator_halo.plan_stats(),
                            retry.operator_applies +
                                retry.recycle_operator_applies);
  };

  {
    OwnedField caller = make_field(92U, cells, 1U, 0U, 8521U, 8522U);
    const bool seeded = static_cast<bool>(seed_cycle_projection(
        cells, identity, cycle_directions, cycle_direction_count, caller,
        linear.workspace, linear.reductions));
    const std::vector<double> caller_before = caller.storage;
    NthPostApplyFailureOperator failing_operator(
        exact_operator, rank, failing_rank, 2U);
    const LinearSolveInvocation invocation{
        rhs, caller.view, identity, control, nullptr};
    const LinearSolveResult failed = solve_fgmres(
        failing_operator, linear.mg, invocation, linear.workspace,
        linear.reductions);
    const Status stale = linear.workspace.recycle_begin_projection_for_test(
        cells, identity.fingerprint);
    const bool failure_matrix =
        seeded && failed.status.code == StatusCode::numerical_failure &&
        failed.status.detail == 8520U &&
        failed.termination == LinearTermination::operator_failure &&
        failed.lowest_failing_rank == failing_rank &&
        failed.operator_applies == 1U &&
        failed.recycle_offered_directions == cycle_direction_count &&
        failed.recycle_projection_attempted &&
        !failed.recycle_projection_accepted &&
        failed.recycle_operator_applies == 1U &&
        caller.storage == caller_before &&
        stale.code == StatusCode::invalid_plan;
    const bool retried = retry_projection(caller);
    passed &= expect(
        failure_matrix && retried, rank,
        "rank-selective real current-A2 failure publishes no pressure correction, consumes once and immediately retries");
  }

  {
    OwnedField caller = make_field(93U, cells, 1U, 0U, 8523U, 8524U);
    const bool seeded = static_cast<bool>(seed_cycle_projection(
        cells, identity, cycle_directions, cycle_direction_count, caller,
        linear.workspace, linear.reductions));
    const std::vector<double> caller_before = caller.storage;
    const Status armed =
        linear.reductions.arm_checked_sum_fault_for_test(7U, failing_rank);
    const LinearSolveInvocation invocation{
        rhs, caller.view, identity, control, nullptr};
    const LinearSolveResult failed = solve_fgmres(
        exact_operator, linear.mg, invocation, linear.workspace,
        linear.reductions);
    linear.reductions.clear_checked_sum_fault_for_test();
    const Status stale = linear.workspace.recycle_begin_projection_for_test(
        cells, identity.fingerprint);
    const bool failure_matrix =
        seeded && static_cast<bool>(armed) &&
        failed.status.code == StatusCode::numerical_failure &&
        failed.status.detail != 0U &&
        failed.termination == LinearTermination::operator_failure &&
        failed.lowest_failing_rank == failing_rank &&
        failed.recycle_offered_directions == cycle_direction_count &&
        failed.recycle_projection_attempted &&
        !failed.recycle_projection_accepted &&
        failed.recycle_operator_applies == cycle_direction_count &&
        caller.storage == caller_before &&
        stale.code == StatusCode::invalid_plan;
    const bool retried = retry_projection(caller);
    passed &= expect(
        failure_matrix && retried, rank,
        "rank-selective real C2 QR reduction failure publishes no pressure correction, consumes once and immediately retries");
  }

  return all_true(passed, world);
}

bool test_negative_density_refresh_collective(MPI_Comm world, int rank) {
  int size = 0;
  MPI_Comm_size(world, &size);
  if (size != 2) return true;

  // A rank-local prerequisite must enter the same refresh Halo gate as the
  // ranks that are locally ready.  Keep this probe on an independent fixture
  // so the failed attempt cannot alter the density-recovery oracle below.
  {
    PeriodicPisoFixture prerequisite_fixture;
    bool prerequisite_passed =
        expect(prerequisite_fixture.initialize(10, world), rank,
               "rank-selective prerequisite fixture compiles");
    if (!all_true(prerequisite_passed, world)) return false;
    const BdfCoefficients prerequisite_bdf{10.0, -15.0, 5.0, 2U};
    const Status prerequisite =
        rank == size - 1 ? Status{StatusCode::numerical_failure, 1591U}
                         : Status{};
    const HaloRuntimeCounters prerequisite_halo_before =
        prerequisite_fixture.halo.runtime_counters();
    PisoIntermediateCertificate prerequisite_certificate{};
    const Status prerequisite_failed = prerequisite_fixture.coupler.refresh(
        prerequisite_fixture.intermediate_input(prerequisite_bdf, 7311U),
        prerequisite_certificate, prerequisite);
    const HaloRuntimeCounters prerequisite_halo_after =
        prerequisite_fixture.halo.runtime_counters();
    const std::array<std::uint64_t, 2U> local_prerequisite_status{
        static_cast<std::uint64_t>(prerequisite_failed.code),
        prerequisite_failed.detail};
    std::array<std::uint64_t, 2U> minimum_prerequisite_status{};
    std::array<std::uint64_t, 2U> maximum_prerequisite_status{};
    const int minimum_prerequisite_result = MPI_Allreduce(
        local_prerequisite_status.data(), minimum_prerequisite_status.data(),
        static_cast<int>(local_prerequisite_status.size()), MPI_UINT64_T,
        MPI_MIN, world);
    const int maximum_prerequisite_result = MPI_Allreduce(
        local_prerequisite_status.data(), maximum_prerequisite_status.data(),
        static_cast<int>(local_prerequisite_status.size()), MPI_UINT64_T,
        MPI_MAX, world);
    prerequisite_passed &= expect(
        minimum_prerequisite_result == MPI_SUCCESS &&
            maximum_prerequisite_result == MPI_SUCCESS &&
            minimum_prerequisite_status == maximum_prerequisite_status &&
            prerequisite_failed.code == StatusCode::numerical_failure &&
            prerequisite_failed.detail == 1591U,
        rank,
        "rank-selective prerequisite publishes the same numerical failure/detail on both ranks");
    prerequisite_passed &= expect(
        !prerequisite_certificate.valid(), rank,
        "rank-selective prerequisite publishes no certificate");
    prerequisite_passed &= expect(
        prerequisite_halo_after.begin_calls ==
                prerequisite_halo_before.begin_calls + 1U &&
            prerequisite_halo_after.finish_calls ==
                prerequisite_halo_before.finish_calls &&
            prerequisite_halo_after.messages_started ==
                prerequisite_halo_before.messages_started &&
            prerequisite_halo_after.bytes_packed ==
                prerequisite_halo_before.bytes_packed &&
            prerequisite_halo_after.bytes_unpacked ==
                prerequisite_halo_before.bytes_unpacked,
        rank,
        "rank-selective prerequisite begins one collective halo attempt without finish/messages/bytes");
    if (!all_true(prerequisite_passed, world)) return false;

    PisoIntermediateCertificate prerequisite_retry_certificate{};
    const Status prerequisite_retry = prerequisite_fixture.coupler.refresh(
        prerequisite_fixture.intermediate_input(prerequisite_bdf, 7312U),
        prerequisite_retry_certificate);
    prerequisite_passed &= expect(
        static_cast<bool>(prerequisite_retry) &&
            prerequisite_retry_certificate.valid(),
        rank, "restored prerequisite refreshes clean corrector one");
    if (!all_true(prerequisite_passed, world)) return false;
  }

  PeriodicPisoFixture fixture;
  bool passed = expect(fixture.initialize(9, world), rank,
                       "negative-density refresh fixture compiles");
  if (!all_true(passed, world)) return false;

  const Int3 bad_cell{0, 0, 0};
  if (rank == size - 1) {
    fixture.density.view.unchecked(bad_cell, 0U) = -1.0;
  }
  const BdfCoefficients bdf{10.0, -15.0, 5.0, 2U};
  const HaloRuntimeCounters failure_halo_before =
      fixture.halo.runtime_counters();
  PisoIntermediateCertificate failed_certificate{};
  const Status failed = fixture.coupler.refresh(
      fixture.intermediate_input(bdf, 7301U), failed_certificate);
  const HaloRuntimeCounters failure_halo_after =
      fixture.halo.runtime_counters();

  const std::array<std::uint64_t, 2U> local_status{
      static_cast<std::uint64_t>(failed.code), failed.detail};
  std::array<std::uint64_t, 2U> minimum_status{};
  std::array<std::uint64_t, 2U> maximum_status{};
  const int minimum_result = MPI_Allreduce(
      local_status.data(), minimum_status.data(),
      static_cast<int>(local_status.size()), MPI_UINT64_T, MPI_MIN, world);
  const int maximum_result = MPI_Allreduce(
      local_status.data(), maximum_status.data(),
      static_cast<int>(local_status.size()), MPI_UINT64_T, MPI_MAX, world);
  passed &= expect(
      minimum_result == MPI_SUCCESS && maximum_result == MPI_SUCCESS &&
          minimum_status == maximum_status &&
          failed.code == StatusCode::numerical_failure && failed.detail == 1504U,
      rank,
      "negative owned density publishes the same numerical failure/detail on both ranks");
  passed &= expect(!failed_certificate.valid(), rank,
                   "failed density refresh publishes no certificate");
  passed &= expect(
      failure_halo_after.begin_calls == failure_halo_before.begin_calls + 1U &&
          failure_halo_after.finish_calls == failure_halo_before.finish_calls &&
          failure_halo_after.messages_started ==
              failure_halo_before.messages_started &&
          failure_halo_after.bytes_packed == failure_halo_before.bytes_packed &&
          failure_halo_after.bytes_unpacked ==
              failure_halo_before.bytes_unpacked,
      rank,
      "failed refresh begins one collective halo attempt without finish/messages/bytes");
  if (!all_true(passed, world)) return false;

  fixture.density.view.unchecked(bad_cell, 0U) = 1.0;
  PisoIntermediateCertificate clean_certificate{};
  const Status clean = fixture.coupler.refresh(
      fixture.intermediate_input(bdf, 7302U), clean_certificate);
  passed &= expect(static_cast<bool>(clean) && clean_certificate.valid(), rank,
                   "restored density refreshes clean corrector one");
  return all_true(passed, world);
}

bool test_post_halo_thermophysical_revalidation_collective(
    MPI_Comm world, int rank) {
  int size = 0;
  MPI_Comm_size(world, &size);
  if (size < 2) return true;

  PeriodicPisoFixture fixture;
  bool passed = expect(fixture.initialize(12, world), rank,
                       "post-halo revalidation fixture compiles");
  if (!all_true(passed, world)) return false;
  const BdfCoefficients bdf{10.0, -15.0, 5.0, 2U};
  const int failing_rank = size - 1;
  detail::arm_piso_post_halo_revalidation_failure_once_for_test(failing_rank);
  const HaloRuntimeCounters before = fixture.halo.runtime_counters();
  PisoIntermediateCertificate rejected_certificate;
  const Status rejected = fixture.coupler.refresh(
      fixture.intermediate_input(bdf, 7351U), rejected_certificate);
  detail::clear_piso_post_halo_revalidation_failure_for_test();
  const HaloRuntimeCounters after = fixture.halo.runtime_counters();

  const std::array<std::uint64_t, 2U> local_status{
      static_cast<std::uint64_t>(rejected.code), rejected.detail};
  std::array<std::uint64_t, 2U> minimum_status{};
  std::array<std::uint64_t, 2U> maximum_status{};
  const int minimum_result = MPI_Allreduce(
      local_status.data(), minimum_status.data(),
      static_cast<int>(local_status.size()), MPI_UINT64_T, MPI_MIN, world);
  const int maximum_result = MPI_Allreduce(
      local_status.data(), maximum_status.data(),
      static_cast<int>(local_status.size()), MPI_UINT64_T, MPI_MAX, world);
  passed &= expect(
      minimum_result == MPI_SUCCESS && maximum_result == MPI_SUCCESS &&
          minimum_status == maximum_status &&
          rejected.code == StatusCode::invalid_plan &&
          rejected.detail == 1503U && !rejected_certificate.valid(),
      rank,
      "single-rank post-halo thermophysical mismatch fails closed identically on every rank");
  passed &= expect(after.begin_calls == before.begin_calls + 1U &&
                       after.finish_calls == before.finish_calls + 1U &&
                       after.bytes_unpacked > before.bytes_unpacked,
                   rank,
                   "post-halo mismatch reaches completed exchange before collective rejection");
  if (!all_true(passed, world)) return false;

  // Periodic/inter-rank ghosts are outside the physical-ghost digest.  Their
  // post-issuance contents may differ until the next halo exchange.
  if (rank == failing_rank)
    fixture.density.view.unchecked({-1, 0, 0}, 0U) = 9.0;
  PisoIntermediateCertificate retried_certificate;
  const Status retried = fixture.coupler.refresh(
      fixture.intermediate_input(bdf, 7352U), retried_certificate);
  passed &= expect(static_cast<bool>(retried) && retried_certificate.valid(),
                   rank,
                   "rank-local periodic/MPI ghost mutation does not invalidate physical authority");
  return all_true(passed, world);
}

bool test_frozen_momentum_candidate_collective(MPI_Comm world, int rank) {
  int size = 0;
  MPI_Comm_size(world, &size);
  PeriodicPisoFixture fixture;
  bool passed = expect(fixture.initialize(9, world, true), rank,
                       "distributed three-species frozen-candidate fixture compiles");
  if (!all_true(passed, world)) return false;
  const Int3 cells = fixture.patch.cells;
  const BdfCoefficients bdf{10.0, -15.0, 5.0, 2U};
  OwnedField base_pressure = make_field(2U, cells, 1U, 0U, 7377U, 8377U);
  OwnedField base_enthalpy = make_field(3U, cells, 1U, 0U, 7378U, 8378U);
  OwnedField base_temperature =
      make_field(4U, cells, 1U, 0U, 7379U, 8379U);
  std::array<OwnedField, 2U> semantic_species{
      make_field(8U, cells, 1U, 0U, 73791U, 83791U),
      make_field(9U, cells, 1U, 0U, 73792U, 83792U)};
  std::array<OwnedField, 2U> raw_species{
      make_field(230U, cells, 1U, 0U, 73793U, 83793U),
      make_field(231U, cells, 1U, 0U, 73794U, 83794U)};
  fill(semantic_species[0U], 0.2);
  fill(semantic_species[1U], 0.3);
  fill(raw_species[0U], 0.2);
  fill(raw_species[1U], 0.3);
  const std::array<double, 2U> base_composition{{0.2, 0.3}};
  const std::array<ConstFieldView, 2U> semantic_species_views{{
      as_const(semantic_species[0U].view),
      as_const(semantic_species[1U].view)}};
  const std::array<ConstFieldView, 2U> raw_species_views{{
      as_const(raw_species[0U].view), as_const(raw_species[1U].view)}};
  fill(base_pressure, 0.0);
  double base_enthalpy_value = 0.0;
  double base_cp = 0.0;
  double base_gas_constant = 0.0;
  Status status = fixture.thermodynamics.mixture_enthalpy(
      300.0, {base_composition.data(), base_composition.size()},
      base_enthalpy_value, base_cp, base_gas_constant);
  passed &= expect(static_cast<bool>(status), rank,
                   "distributed frozen-candidate base enthalpy closes at 300 K");
  if (!all_true(passed, world)) return false;
  fill(base_enthalpy, base_enthalpy_value);
  fill(base_temperature, 0.0);
  PisoIntermediateInput intermediate_input =
      fixture.intermediate_input(bdf, 7381U);
  constexpr double absolute_pressure_reference = 101325.0;
  intermediate_input.thermophysical_boundary.binding.pressure_reference =
      absolute_pressure_reference;
  status = {};
  for (std::int32_t z = 0; z < cells.z && status; ++z)
    for (std::int32_t y = 0; y < cells.y && status; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        ThermoState thermo;
        status = fixture.thermodynamics.evaluate_from_reference_pressure(
            absolute_pressure_reference,
            base_pressure.view.unchecked(cell, 0U),
            base_enthalpy.view.unchecked(cell, 0U),
            {base_composition.data(), base_composition.size()},
            {fixture.velocity.view.unchecked(cell, 0U),
             fixture.velocity.view.unchecked(cell, 1U),
             fixture.velocity.view.unchecked(cell, 2U)},
            thermo);
        if (!status) break;
        fixture.density.view.unchecked(cell, 0U) = thermo.rho;
        base_temperature.view.unchecked(cell, 0U) = thermo.temperature;
      }
  passed &= expect(static_cast<bool>(status), rank,
                   "distributed frozen-candidate base p/h evaluates through EOS");
  if (!all_true(passed, world)) return false;
  PisoIntermediateCertificate intermediate;
  status = fixture.coupler.refresh(intermediate_input, intermediate);
  passed &= expect(static_cast<bool>(status) && intermediate.valid(), rank,
                   "distributed frozen-candidate C1 refreshes");
  if (!all_true(passed, world)) return false;

  OwnedField accepted = make_field(50U, cells, 1U, 0U, 7382U, 8382U);
  OwnedField previous = make_field(51U, cells, 1U, 0U, 7383U, 8383U);
  OwnedField drho_dp = make_field(6U, cells, 1U, 0U, 7384U, 8384U);
  OwnedField diagonal = make_field(53U, cells, 1U, 0U, 7385U, 8385U);
  OwnedField rhs = make_field(54U, cells, 1U, 0U, 7386U, 8386U);
  fill(accepted, 1.0);
  fill(previous, 1.0);
  for (std::int32_t z = 0; z < cells.z && status; ++z)
    for (std::int32_t y = 0; y < cells.y && status; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        ThermoState thermo;
        status = fixture.thermodynamics.evaluate_from_reference_pressure(
            absolute_pressure_reference,
            base_pressure.view.unchecked(cell, 0U),
            base_enthalpy.view.unchecked(cell, 0U),
            {base_composition.data(), base_composition.size()},
            {fixture.velocity.view.unchecked(cell, 0U),
             fixture.velocity.view.unchecked(cell, 1U),
             fixture.velocity.view.unchecked(cell, 2U)},
            thermo, base_temperature.view.unchecked(cell, 0U));
        if (!status) break;
        drho_dp.view.unchecked(cell, 0U) = thermo.drho_dp_hY;
      }
  passed &= expect(static_cast<bool>(status), rank,
                   "distributed frozen-candidate compressibility replays from EOS");
  if (!all_true(passed, world)) return false;
  const PressureCorrectionInput pressure_input{
      intermediate,
      intermediate_input.pressure_reference,
      as_const(fixture.density.view),
      as_const(accepted.view),
      as_const(previous.view),
      as_const(drho_dp.view),
      bdf,
      intermediate_input.momentum.time,
      intermediate_input.momentum.geometry,
      intermediate_input.numeric_boundary};
  PressureCorrectionCertificate pressure;
  status = fixture.coupler.assemble_pressure_system(
      pressure_input, {diagonal.view, rhs.view}, pressure);
  passed &= expect(static_cast<bool>(status) && pressure.valid(), rank,
                   "distributed frozen-candidate pressure system assembles");
  if (!all_true(passed, world)) return false;

  PisoFrozenMomentumStageAuthority authority;
  status = fixture.coupler.make_frozen_momentum_stage_authority(
      intermediate, pressure, authority);
  passed &= expect(static_cast<bool>(status) && authority.valid(), rank,
                   "distributed frozen-candidate authority freezes same-target state");
  if (!all_true(passed, world)) return false;
  OwnedField raw_dp = make_field(90U, cells, 1U, 1U, 7387U, 8387U);
  OwnedField scaled_dp = make_field(200U, cells, 1U, 1U, 7388U, 8388U);
  OwnedField candidate_velocity =
      make_field(201U, cells, 3U, 1U, 7389U, 8389U);
  OwnedField candidate_density =
      make_field(202U, cells, 1U, 1U, 7390U, 8390U);
  fill(raw_dp, 0.0);
  fill(scaled_dp, -1.0);
  fill(candidate_velocity, -2.0);
  const auto restore_candidate_density = [&]() noexcept {
    for (std::int32_t z = -1; z <= cells.z; ++z)
      for (std::int32_t y = -1; y <= cells.y; ++y)
        for (std::int32_t x = -1; x <= cells.x; ++x) {
          const unsigned outside =
              static_cast<unsigned>(x < 0 || x >= cells.x) +
              static_cast<unsigned>(y < 0 || y >= cells.y) +
              static_cast<unsigned>(z < 0 || z >= cells.z);
          if (outside <= 1U)
            candidate_density.view.unchecked({x, y, z}, 0U) =
                fixture.density.view.unchecked({x, y, z}, 0U);
        }
  };
  restore_candidate_density();
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 local{x, y, z};
        const Int3 global{fixture.patch.begin.x + x,
                          fixture.patch.begin.y + y,
                          fixture.patch.begin.z + z};
        raw_dp.view.unchecked(local, 0U) =
            50.0 * periodic_value(global.x,
                                  fixture.geometry.global_cells().x, 0.3) +
            20.0 * periodic_value(global.y,
                                  fixture.geometry.global_cells().y, -0.2);
      }
  const std::array<HaloFieldSpec, 1U> candidate_contract{{
      {scaled_dp.view.field, 1U, 1U}}};
  HaloEngine candidate_halo;
  HaloEngine foreign_halo;
  passed &= expect(
      static_cast<bool>(candidate_halo.reserve(
          world, fixture.patch,
          {candidate_contract.data(), candidate_contract.size()},
          fixture.boundary.halo_topology())) &&
          static_cast<bool>(foreign_halo.reserve(
              world, fixture.patch,
              {candidate_contract.data(), candidate_contract.size()},
              fixture.boundary.halo_topology())),
      rank, "distributed independent candidate halos reserve");
  FaceFluxStorage candidate_flux_storage;
  FaceFluxView candidate_flux;
  passed &= expect(
      static_cast<bool>(FaceFluxStorage::allocate_workspace(
          cells, 1U, candidate_flux_storage)) &&
          static_cast<bool>(candidate_flux_storage.workspace_view(
              0U, 7391U, candidate_flux)),
      rank, "distributed candidate flux scratch allocates");
  if (!all_true(passed, world)) return false;

  PisoFrozenMomentumPressureStageCertificate pressure_stage;
  status = fixture.coupler.form_frozen_momentum_scaled_pressure(
      authority, as_const(raw_dp.view), candidate_halo, 0.5, scaled_dp.view,
      pressure_stage);
  passed &= expect(static_cast<bool>(status) && pressure_stage.valid() &&
                       scaled_dp.view.revision == 7388U,
                   rank,
                   "distributed form preserves caller-issued write revision");
  PisoFrozenMomentumVelocityStageCertificate skipped;
  const Status skipped_status =
      fixture.coupler.stage_frozen_momentum_velocity(
          authority, pressure_stage, candidate_halo, as_const(scaled_dp.view),
          candidate_velocity.view, skipped);
  passed &= expect(skipped_status.code == StatusCode::invalid_plan &&
                       !skipped.valid(),
                   rank,
                   "distributed velocity rejects skipped correction halo");

  std::array<FieldView, 1U> foreign_fields{scaled_dp.view};
  HaloTicket foreign_ticket;
  status = foreign_halo.begin(7392U,
                              {foreign_fields.data(), foreign_fields.size()},
                              foreign_ticket);
  if (status)
    status = foreign_halo.finish(
        foreign_ticket, {foreign_fields.data(), foreign_fields.size()});
  scaled_dp.view = foreign_fields[0U];
  PisoFrozenMomentumVelocityStageCertificate foreign;
  const Status foreign_status =
      fixture.coupler.stage_frozen_momentum_velocity(
          authority, pressure_stage, foreign_halo, as_const(scaled_dp.view),
          candidate_velocity.view, foreign);
  passed &= expect(static_cast<bool>(status) &&
                       foreign_status.code == StatusCode::invalid_plan &&
                       !foreign.valid(),
                   rank,
                   "distributed velocity rejects a foreign compatible halo");

  std::array<FieldView, 1U> candidate_fields{scaled_dp.view};
  HaloTicket candidate_ticket;
  status = candidate_halo.begin(
      7393U, {candidate_fields.data(), candidate_fields.size()},
      candidate_ticket);
  if (status)
    status = candidate_halo.finish(
        candidate_ticket, {candidate_fields.data(), candidate_fields.size()});
  scaled_dp.view = candidate_fields[0U];
  PisoFrozenMomentumVelocityStageCertificate velocity_stage;
  if (status)
    status = fixture.coupler.stage_frozen_momentum_velocity(
        authority, pressure_stage, candidate_halo, as_const(scaled_dp.view),
        candidate_velocity.view, velocity_stage);
  passed &= expect(static_cast<bool>(status) && velocity_stage.valid(), rank,
                   "distributed bound halo stages candidate velocity");
  PisoFrozenMomentumPressureStageCertificate live_halo_rejected;
  const Status live_halo_status =
      fixture.coupler.form_frozen_momentum_scaled_pressure(
          authority, as_const(raw_dp.view), fixture.correction_halo, 0.5,
          scaled_dp.view, live_halo_rejected);
  passed &= expect(live_halo_status.code == StatusCode::invalid_plan &&
                       !live_halo_rejected.valid(),
                   rank, "distributed form rejects the live correction halo");
  if (!all_true(passed, world)) return false;

  if (rank == size - 1)
    candidate_density.view.unchecked({0, 0, 0}, 0U) = -1.0;
  PisoFrozenMomentumFluxStageCertificate bad_flux;
  const Status bad_status = fixture.coupler.stage_frozen_momentum_flux(
      authority, velocity_stage, as_const(candidate_density.view),
      candidate_flux, bad_flux);
  const std::array<std::uint64_t, 2U> local_status{
      static_cast<std::uint64_t>(bad_status.code), bad_status.detail};
  std::array<std::uint64_t, 2U> minimum_status{};
  std::array<std::uint64_t, 2U> maximum_status{};
  MPI_Allreduce(local_status.data(), minimum_status.data(), 2, MPI_UINT64_T,
                MPI_MIN, world);
  MPI_Allreduce(local_status.data(), maximum_status.data(), 2, MPI_UINT64_T,
                MPI_MAX, world);
  passed &= expect(
      minimum_status == maximum_status &&
          bad_status.code == StatusCode::numerical_failure &&
          !bad_flux.valid(),
      rank,
      "one-rank nonpositive candidate fails identically with no flux certificate");
  candidate_density.view.unchecked({0, 0, 0}, 0U) =
      fixture.density.view.unchecked({0, 0, 0}, 0U);
  PisoFrozenMomentumFluxStageCertificate recovered_flux;
  status = fixture.coupler.stage_frozen_momentum_flux(
      authority, velocity_stage, as_const(candidate_density.view),
      candidate_flux, recovered_flux);
  PisoFrozenMomentumStageAuthority replay;
  const Status replay_status =
      fixture.coupler.make_frozen_momentum_stage_authority(
          intermediate, pressure, replay);
  passed &= expect(static_cast<bool>(status) && recovered_flux.valid() &&
                       static_cast<bool>(replay_status) && replay.valid(),
                   rank,
                   "collective candidate failure consumes no current pressure/intermediate authority");
  const RevisionToken current_revision = scaled_dp.view.revision;
  ++scaled_dp.view.revision;
  PisoFrozenMomentumVelocityStageCertificate stale_revision;
  const Status stale_status =
      fixture.coupler.stage_frozen_momentum_velocity(
          authority, pressure_stage, candidate_halo, as_const(scaled_dp.view),
          candidate_velocity.view, stale_revision);
  scaled_dp.view.revision = current_revision;
  passed &= expect(stale_status.code == StatusCode::invalid_plan &&
                       !stale_revision.valid(),
                   rank,
                   "distributed velocity rejects stale caller revision/ghost proof");

  ++scaled_dp.view.revision;
  ++candidate_velocity.view.revision;
  ++candidate_density.view.revision;
  ++candidate_flux.revision;
  restore_candidate_density();
  PisoFrozenMomentumPressureStageCertificate baseline_pressure_stage;
  status = fixture.coupler.form_frozen_momentum_scaled_pressure(
      authority, as_const(raw_dp.view), candidate_halo, 0.0, scaled_dp.view,
      baseline_pressure_stage);
  std::array<FieldView, 1U> baseline_halo_fields{scaled_dp.view};
  HaloTicket baseline_halo_ticket;
  if (status)
    status = candidate_halo.begin(
        7394U, {baseline_halo_fields.data(), baseline_halo_fields.size()},
        baseline_halo_ticket);
  if (status)
    status = candidate_halo.finish(
        baseline_halo_ticket,
        {baseline_halo_fields.data(), baseline_halo_fields.size()});
  scaled_dp.view = baseline_halo_fields[0U];
  PisoFrozenMomentumVelocityStageCertificate baseline_velocity_stage;
  if (status)
    status = fixture.coupler.stage_frozen_momentum_velocity(
        authority, baseline_pressure_stage, candidate_halo,
        as_const(scaled_dp.view), candidate_velocity.view,
        baseline_velocity_stage);
  PisoFrozenMomentumFluxStageCertificate baseline_flux_stage;
  if (status)
    status = fixture.coupler.stage_frozen_momentum_flux(
        authority, baseline_velocity_stage, as_const(candidate_density.view),
        candidate_flux, baseline_flux_stage);
  passed &= expect(static_cast<bool>(status) &&
                       baseline_pressure_stage.valid() &&
                       baseline_velocity_stage.valid() &&
                       baseline_flux_stage.valid(),
                   rank,
                   "distributed alpha-zero baseline stages exact U and base flux");

  OwnedField raw_dh = make_field(93U, cells, 1U, 0U, 7395U, 8395U);
  OwnedField scaled_dh = make_field(207U, cells, 1U, 0U, 7396U, 8396U);
  OwnedField candidate_pressure =
      make_field(208U, cells, 1U, 0U, 7397U, 8397U);
  OwnedField candidate_enthalpy =
      make_field(203U, cells, 1U, 0U, 7398U, 8398U);
  OwnedField candidate_temperature =
      make_field(204U, cells, 1U, 0U, 7399U, 8399U);
  fill(raw_dh, 0.0);
  fill(scaled_dh, 0.0);
  fill(candidate_enthalpy, base_enthalpy_value);
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x)
        candidate_temperature.view.unchecked({x, y, z}, 0U) =
            base_temperature.view.unchecked({x, y, z}, 0U);
  PisoExactEosClosureIdentity closure;
  closure.thermodynamics =
      intermediate_input.pressure_reference.thermodynamics;
  closure.pressure_reference =
      intermediate_input.pressure_reference.pressure_reference;
  closure.composition = exact_composition_identity_for_test(
      closure.thermodynamics,
      {raw_species_views.data(), raw_species_views.size()}, cells);
  closure.pressure_state =
      make_piso_field_revision_identity(as_const(base_pressure.view));
  closure.pressure_correction =
      make_piso_field_revision_identity(as_const(scaled_dp.view));
  closure.enthalpy_state =
      make_piso_field_revision_identity(as_const(base_enthalpy.view));
  closure.enthalpy_correction =
      make_piso_field_revision_identity(as_const(scaled_dh.view));
  closure.candidate_enthalpy =
      make_piso_field_revision_identity(as_const(candidate_enthalpy.view));
  closure.candidate_density =
      make_piso_field_revision_identity(as_const(candidate_density.view));
  closure.candidate_temperature =
      make_piso_field_revision_identity(as_const(candidate_temperature.view));
  closure.closure = 7401U;
  PisoExactThermodynamicCandidateView thermodynamic;
  thermodynamic.enthalpy = as_const(candidate_enthalpy.view);
  thermodynamic.density = as_const(candidate_density.view);
  thermodynamic.temperature = as_const(candidate_temperature.view);
  thermodynamic.closure = closure;
  thermodynamic.pressure_compressibility = as_const(drho_dp.view);
  thermodynamic.independent_species =
      {semantic_species_views.data(), semantic_species_views.size()};
  ReductionEngine stationary_reductions;
  status = ReductionEngine::compile(world, ReductionMode::mpi_allreduce, 2U,
                                    stationary_reductions);
  if (status)
    status = fixture.prepare_closed_gauge(
        pressure, intermediate_input.pressure_reference,
        absolute_pressure_reference, as_const(base_pressure.view),
        as_const(scaled_dp.view), thermodynamic.pressure_compressibility,
        closure.closure, stationary_reductions, thermodynamic.closed_gauge);
  if (status) {
    for (std::int32_t z = 0; z < cells.z; ++z)
      for (std::int32_t y = 0; y < cells.y; ++y)
        for (std::int32_t x = 0; x < cells.x; ++x)
          candidate_pressure.view.unchecked({x, y, z}, 0U) =
              base_pressure.view.unchecked({x, y, z}, 0U) +
              scaled_dp.view.unchecked({x, y, z}, 0U) -
              thermodynamic.closed_gauge.shift;
  }
  const PisoCoupledStateView baseline_state{
      fixture.velocity.view, base_pressure.view, base_enthalpy.view,
      fixture.density.view, base_temperature.view};
  PisoFrozenMomentumExactCandidateInput baseline_exact_input{
      as_const(raw_dh.view),
      as_const(scaled_dp.view),
      as_const(scaled_dh.view),
      baseline_state,
      as_const(candidate_pressure.view),
      thermodynamic,
      as_const(candidate_velocity.view),
      as_const(candidate_flux),
      {},
      {raw_species_views.data(), raw_species_views.size()}};
  PisoFrozenMomentumExactCandidateCertificate exact_baseline;
  if (status)
    status = fixture.coupler.certify_frozen_momentum_exact_baseline(
        authority, baseline_pressure_stage, baseline_velocity_stage,
        baseline_flux_stage, baseline_exact_input,
        stationary_reductions, exact_baseline);
  passed &= expect(static_cast<bool>(status) && exact_baseline.valid(), rank,
                   "distributed physical three-species alpha-zero exact baseline certifies");
  if (!all_true(passed, world)) return false;

  std::array<OwnedField, 2U> replay_raw_species{
      make_field(232U, cells, 1U, 0U, 7402U, 8402U),
      make_field(233U, cells, 1U, 0U, 7403U, 8403U)};
  fill(replay_raw_species[0U], 0.2);
  fill(replay_raw_species[1U], 0.3);
  const std::array<ConstFieldView, 2U> replay_raw_views{{
      as_const(replay_raw_species[0U].view),
      as_const(replay_raw_species[1U].view)}};
  PisoFrozenMomentumExactCandidateInput replay_input = baseline_exact_input;
  replay_input.independent_species =
      {replay_raw_views.data(), replay_raw_views.size()};
  PisoFrozenMomentumExactCandidateCertificate replayed_baseline;
  Status replayed_status =
      fixture.coupler.certify_frozen_momentum_exact_baseline(
          authority, baseline_pressure_stage, baseline_velocity_stage,
          baseline_flux_stage, replay_input, stationary_reductions,
          replayed_baseline);
  passed &= expect(
      static_cast<bool>(replayed_status) && replayed_baseline.valid() &&
          replayed_baseline.canonical_lineage() ==
              exact_baseline.canonical_lineage() &&
          replayed_baseline.candidate_state_provenance() ==
              exact_baseline.candidate_state_provenance() &&
          replayed_baseline.scratch_binding() !=
              exact_baseline.scratch_binding(),
      rank,
      "same physical periodic composition replays canonical provenance across distinct scratch IDs while retaining a distinct binding");
  if (!all_true(passed, world)) return false;

  const auto rejects_exact = [&](PisoFrozenMomentumExactCandidateInput input) {
    PisoFrozenMomentumExactCandidateCertificate rejected = exact_baseline;
    const Status rejected_status =
        fixture.coupler.certify_frozen_momentum_exact_baseline(
            authority, baseline_pressure_stage, baseline_velocity_stage,
            baseline_flux_stage, input, stationary_reductions, rejected);
    return rejected_status.code == StatusCode::invalid_plan &&
           !rejected.valid();
  };
  PisoFrozenMomentumExactCandidateInput poisoned_input = baseline_exact_input;
  if (rank == size - 1) {
    poisoned_input.independent_species =
        {raw_species_views.data(), raw_species_views.size() - 1U};
  }
  passed &= expect(rejects_exact(poisoned_input), rank,
                   "one-rank wrong periodic composition count rejects collectively and clears a reused certificate");
  if (!all_true(passed, world)) return false;

  std::array<ConstFieldView, 2U> reordered_semantic = semantic_species_views;
  if (rank == size - 1) std::swap(reordered_semantic[0U], reordered_semantic[1U]);
  poisoned_input = baseline_exact_input;
  poisoned_input.thermodynamic.independent_species =
      {reordered_semantic.data(), reordered_semantic.size()};
  passed &= expect(rejects_exact(poisoned_input), rank,
                   "one-rank periodic semantic species reorder rejects collectively");
  if (!all_true(passed, world)) return false;

  std::array<ConstFieldView, 2U> aliased_raw_species = raw_species_views;
  if (rank == size - 1) {
    aliased_raw_species[0U] = as_const(candidate_pressure.view);
    aliased_raw_species[0U].field = raw_species_views[0U].field;
  }
  poisoned_input = baseline_exact_input;
  poisoned_input.independent_species =
      {aliased_raw_species.data(), aliased_raw_species.size()};
  passed &= expect(
      rejects_exact(poisoned_input), rank,
      "one-rank raw species aliasing candidate pressure is collectively rejected before exact issuance");
  if (!all_true(passed, world)) return false;

  const Int3 poison_cell{0, 0, 0};
  const double raw_species_before =
      raw_species[0U].view.unchecked(poison_cell, 0U);
  if (rank == size - 1)
    raw_species[0U].view.unchecked(poison_cell, 0U) =
        std::nextafter(raw_species_before,
                       std::numeric_limits<double>::infinity());
  passed &= expect(rejects_exact(baseline_exact_input), rank,
                   "one-rank no-revision periodic species mutation rejects collectively");
  raw_species[0U].view.unchecked(poison_cell, 0U) = raw_species_before;
  if (!all_true(passed, world)) return false;

  const double temperature_before =
      candidate_temperature.view.unchecked(poison_cell, 0U);
  if (rank == size - 1)
    candidate_temperature.view.unchecked(poison_cell, 0U) =
        std::nextafter(temperature_before,
                       std::numeric_limits<double>::infinity());
  passed &= expect(rejects_exact(baseline_exact_input), rank,
                   "one-rank no-revision periodic temperature mutation rejects exact EOS");
  candidate_temperature.view.unchecked(poison_cell, 0U) = temperature_before;
  if (!all_true(passed, world)) return false;

  const double density_before =
      candidate_density.view.unchecked(poison_cell, 0U);
  if (rank == size - 1)
    candidate_density.view.unchecked(poison_cell, 0U) =
        std::nextafter(density_before,
                       std::numeric_limits<double>::infinity());
  passed &= expect(rejects_exact(baseline_exact_input), rank,
                   "one-rank no-revision periodic density mutation invalidates the staged flux replay");
  candidate_density.view.unchecked(poison_cell, 0U) = density_before;
  if (!all_true(passed, world)) return false;

  const double compressibility_before =
      drho_dp.view.unchecked(poison_cell, 0U);
  if (rank == size - 1)
    drho_dp.view.unchecked(poison_cell, 0U) =
        std::nextafter(compressibility_before,
                       std::numeric_limits<double>::infinity());
  passed &= expect(rejects_exact(baseline_exact_input), rank,
                   "one-rank no-revision periodic compressibility mutation rejects the closed-gauge EOS replay");
  drho_dp.view.unchecked(poison_cell, 0U) = compressibility_before;
  if (!all_true(passed, world)) return false;

  if (size > 1) {
    fixture.coupler.set_frozen_stationary_tolerances_for_test(
        rank == size - 1 ? 2.0e-10 : 1.0e-10, 1.0e-10);
    PressureEnergyStationaryCertificate mismatched_stationary;
    const Status mismatch_status =
        fixture.coupler.certify_frozen_momentum_stationary(
            authority, replayed_baseline, 0.0, 0.0, stationary_reductions,
            mismatched_stationary);
    const std::array<std::uint64_t, 2U> local_mismatch{
        static_cast<std::uint64_t>(mismatch_status.code),
        mismatch_status.detail};
    std::array<std::uint64_t, 2U> minimum_mismatch{};
    std::array<std::uint64_t, 2U> maximum_mismatch{};
    MPI_Allreduce(local_mismatch.data(), minimum_mismatch.data(), 2,
                  MPI_UINT64_T, MPI_MIN, world);
    MPI_Allreduce(local_mismatch.data(), maximum_mismatch.data(), 2,
                  MPI_UINT64_T, MPI_MAX, world);
    passed &= expect(
        minimum_mismatch == maximum_mismatch && !mismatch_status &&
            !mismatched_stationary.valid(),
        rank,
        "rank-local stationary tolerance mismatch fails the shared consensus contract");
  }
  fixture.coupler.set_frozen_stationary_tolerances_for_test(1.0e-10,
                                                             1.0e-10);
  ReductionEngine restored_reductions;
  status = ReductionEngine::compile(world, ReductionMode::mpi_allreduce, 7U,
                                    restored_reductions);
  PressureEnergyStationaryCertificate stationary;
  if (status)
    status = fixture.coupler.certify_frozen_momentum_stationary(
        authority, replayed_baseline, 0.0, 0.0, restored_reductions,
        stationary);
  passed &= expect(static_cast<bool>(status) && stationary.valid(), rank,
                   "restored stationary tolerances reuse the exact authority");

  OwnedField positive_scaled_dp =
      make_field(200U, cells, 1U, 1U, 7501U, 8501U);
  OwnedField positive_scaled_dh =
      make_field(207U, cells, 1U, 0U, 7502U, 8502U);
  OwnedField positive_pressure =
      make_field(208U, cells, 1U, 0U, 7503U, 8503U);
  OwnedField positive_enthalpy =
      make_field(203U, cells, 1U, 0U, 7504U, 8504U);
  OwnedField positive_density =
      make_field(202U, cells, 1U, 1U, 7505U, 8505U);
  OwnedField positive_temperature =
      make_field(204U, cells, 1U, 0U, 7506U, 8506U);
  OwnedField positive_compressibility =
      make_field(6U, cells, 1U, 0U, 7507U, 8507U);
  OwnedField positive_velocity =
      make_field(201U, cells, 3U, 1U, 7508U, 8508U);
  fill(positive_scaled_dp, 0.0);
  fill(positive_scaled_dh, 0.0);
  fill(positive_pressure, 0.0);
  fill(positive_enthalpy, base_enthalpy_value);
  fill(positive_density, 0.0);
  fill(positive_temperature, 0.0);
  fill(positive_compressibility, 0.0);
  fill(positive_velocity, 0.0);
  const std::array<HaloFieldSpec, 1U> positive_correction_contract{{
      {positive_scaled_dp.view.field, 1U, 1U}}};
  HaloEngine positive_correction_halo;
  status = positive_correction_halo.reserve(
      world, fixture.patch,
      {positive_correction_contract.data(),
       positive_correction_contract.size()},
      fixture.boundary.halo_topology());
  PisoFrozenMomentumPressureStageCertificate positive_pressure_stage;
  if (status)
    status = fixture.coupler.form_frozen_momentum_scaled_pressure(
        authority, as_const(raw_dp.view), positive_correction_halo, 0.25,
        positive_scaled_dp.view, positive_pressure_stage);
  std::array<FieldView, 1U> positive_correction_fields{
      positive_scaled_dp.view};
  HaloTicket positive_correction_ticket;
  if (status)
    status = positive_correction_halo.begin(
        7509U,
        {positive_correction_fields.data(), positive_correction_fields.size()},
        positive_correction_ticket);
  if (status)
    status = positive_correction_halo.finish(
        positive_correction_ticket,
        {positive_correction_fields.data(), positive_correction_fields.size()});
  positive_scaled_dp.view = positive_correction_fields[0U];
  PisoFrozenMomentumVelocityStageCertificate positive_velocity_stage;
  if (status)
    status = fixture.coupler.stage_frozen_momentum_velocity(
        authority, positive_pressure_stage, positive_correction_halo,
        as_const(positive_scaled_dp.view), positive_velocity.view,
        positive_velocity_stage);
  for (std::int32_t z = -1; z <= cells.z && status; ++z)
    for (std::int32_t y = -1; y <= cells.y && status; ++y)
      for (std::int32_t x = -1; x <= cells.x; ++x) {
        const unsigned outside =
            static_cast<unsigned>(x < 0 || x >= cells.x) +
            static_cast<unsigned>(y < 0 || y >= cells.y) +
            static_cast<unsigned>(z < 0 || z >= cells.z);
        if (outside > 1U) continue;
        const Int3 cell{x, y, z};
        const Int3 wrapped{(x % cells.x + cells.x) % cells.x,
                           (y % cells.y + cells.y) % cells.y,
                           (z % cells.z + cells.z) % cells.z};
        ThermoState replayed;
        status = fixture.thermodynamics.evaluate_from_reference_pressure(
            absolute_pressure_reference,
            base_pressure.view.unchecked(wrapped, 0U) +
                positive_scaled_dp.view.unchecked(cell, 0U),
            positive_enthalpy.view.unchecked(wrapped, 0U),
            {base_composition.data(), base_composition.size()},
            {positive_velocity.view.unchecked(cell, 0U),
             positive_velocity.view.unchecked(cell, 1U),
             positive_velocity.view.unchecked(cell, 2U)},
            replayed, base_temperature.view.unchecked(wrapped, 0U));
        if (!status) break;
        positive_density.view.unchecked(cell, 0U) = replayed.rho;
        if (outside == 0U) {
          positive_temperature.view.unchecked(wrapped, 0U) =
              replayed.temperature;
          positive_compressibility.view.unchecked(wrapped, 0U) =
              replayed.drho_dp_hY;
        }
      }
  PisoExactEosClosureIdentity positive_closure;
  positive_closure.thermodynamics = closure.thermodynamics;
  positive_closure.pressure_reference = closure.pressure_reference;
  positive_closure.composition = closure.composition;
  positive_closure.pressure_state = closure.pressure_state;
  positive_closure.pressure_correction = make_piso_field_revision_identity(
      as_const(positive_scaled_dp.view));
  positive_closure.enthalpy_state = closure.enthalpy_state;
  positive_closure.enthalpy_correction = make_piso_field_revision_identity(
      as_const(positive_scaled_dh.view));
  positive_closure.candidate_enthalpy = make_piso_field_revision_identity(
      as_const(positive_enthalpy.view));
  positive_closure.candidate_density = make_piso_field_revision_identity(
      as_const(positive_density.view));
  positive_closure.candidate_temperature = make_piso_field_revision_identity(
      as_const(positive_temperature.view));
  positive_closure.closure = 7510U;
  PisoExactThermodynamicCandidateView positive_thermodynamic;
  positive_thermodynamic.enthalpy = as_const(positive_enthalpy.view);
  positive_thermodynamic.density = as_const(positive_density.view);
  positive_thermodynamic.temperature = as_const(positive_temperature.view);
  positive_thermodynamic.closure = positive_closure;
  positive_thermodynamic.pressure_compressibility =
      as_const(positive_compressibility.view);
  positive_thermodynamic.independent_species =
      {semantic_species_views.data(), semantic_species_views.size()};
  if (status)
    status = fixture.prepare_closed_gauge(
        pressure, intermediate_input.pressure_reference,
        absolute_pressure_reference, as_const(base_pressure.view),
        as_const(positive_scaled_dp.view),
        positive_thermodynamic.pressure_compressibility,
        positive_closure.closure, restored_reductions,
        positive_thermodynamic.closed_gauge);
  if (status) {
    for (std::int32_t z = 0; z < cells.z; ++z)
      for (std::int32_t y = 0; y < cells.y; ++y)
        for (std::int32_t x = 0; x < cells.x; ++x)
          positive_pressure.view.unchecked({x, y, z}, 0U) =
              base_pressure.view.unchecked({x, y, z}, 0U) +
              positive_scaled_dp.view.unchecked({x, y, z}, 0U) -
              positive_thermodynamic.closed_gauge.shift;
  }
  FaceFluxStorage positive_flux_storage;
  FaceFluxView positive_flux;
  if (status)
    status = FaceFluxStorage::allocate_workspace(cells, 1U,
                                                 positive_flux_storage);
  if (status)
    status = positive_flux_storage.workspace_view(0U, 7511U,
                                                  positive_flux);
  PisoFrozenMomentumFluxStageCertificate positive_flux_stage;
  if (status)
    status = fixture.coupler.stage_frozen_momentum_flux(
        authority, positive_velocity_stage, as_const(positive_density.view),
        positive_flux, positive_flux_stage);
  const PisoFrozenMomentumExactCandidateInput positive_input{
      as_const(raw_dh.view),
      as_const(positive_scaled_dp.view),
      as_const(positive_scaled_dh.view),
      baseline_state,
      as_const(positive_pressure.view),
      positive_thermodynamic,
      as_const(positive_velocity.view),
      as_const(positive_flux),
      {},
      {raw_species_views.data(), raw_species_views.size()}};
  const double positive_species_before =
      raw_species[1U].view.unchecked(poison_cell, 0U);
  if (rank == size - 1)
    raw_species[1U].view.unchecked(poison_cell, 0U) =
        std::nextafter(positive_species_before,
                       std::numeric_limits<double>::infinity());
  PisoFrozenMomentumExactCandidateCertificate rejected_positive =
      replayed_baseline;
  Status rejected_positive_status = status;
  if (rejected_positive_status)
    rejected_positive_status =
        fixture.coupler.certify_frozen_momentum_exact_candidate(
            authority, replayed_baseline, positive_pressure_stage,
            positive_velocity_stage, positive_flux_stage, positive_input,
            restored_reductions, rejected_positive);
  passed &= expect(rejected_positive_status.code == StatusCode::invalid_plan &&
                       !rejected_positive.valid(),
                   rank,
                   "one-rank positive-alpha periodic composition mutation rejects collectively and clears a reused certificate");
  raw_species[1U].view.unchecked(poison_cell, 0U) = positive_species_before;
  PisoFrozenMomentumExactCandidateCertificate positive_exact;
  if (status)
    status = fixture.coupler.certify_frozen_momentum_exact_candidate(
        authority, replayed_baseline, positive_pressure_stage,
        positive_velocity_stage, positive_flux_stage, positive_input,
        restored_reductions, positive_exact);
  passed &= expect(static_cast<bool>(status) && positive_exact.valid() &&
                       positive_exact.alpha() == 0.25 &&
                       positive_exact.baseline_state_provenance() ==
                           replayed_baseline.candidate_state_provenance() &&
                       positive_exact.baseline_mass_flux_provenance() ==
                           replayed_baseline.candidate_mass_flux_provenance(),
                   rank,
                   "distributed physical three-species positive-alpha exact candidate certifies against the alpha-zero baseline");

  FieldRegistry publication_registry;
  FieldSchema publication_schema;
  std::array<ArenaFieldRequest, 10U> publication_requests{};
  for (std::size_t field = 0U; field < publication_requests.size(); ++field) {
    FieldId declared = 0U;
    const std::string name =
        "periodic.terminal.composition." + std::to_string(field);
    if (status)
      status = publication_registry.declare_field(name, 1U, 0U, declared);
    if (status && declared != field)
      status = {StatusCode::invalid_plan, 1U};
    publication_requests[field] = {
        declared, {1, 1, 1}, {0U},
        field == 8U || field == 9U ? FieldLifetime::state_layer
                                   : FieldLifetime::step_scratch};
  }
  if (status) status = publication_registry.freeze(publication_schema);
  ArenaLayout publication_layout;
  StateLayers publication_layers;
  AttemptTransaction publication_transaction;
  FaceFluxStorage publication_flux_storage;
  FinalFaceFluxAuthority publication_flux_authority;
  FinalFaceFluxWriter publication_flux_writer;
  if (status)
    status = ArenaLayout::compile(
        publication_schema,
        {publication_requests.data(), publication_requests.size()},
        publication_layout);
  if (status)
    status = StateLayers::allocate(publication_layout, publication_layers);
  if (status)
    status = AttemptTransaction::create(
        publication_layers.field_count(), 1U,
        publication_layers.field_count(), publication_transaction);
  if (status)
    status = FaceFluxStorage::allocate_final(cells,
                                             publication_flux_storage);
  if (status)
    status = publication_flux_authority.claim(
        fixture.piso.pressure_stage(), fixture.piso.final_flux_slot(),
        publication_transaction, publication_flux_writer);
  if (status) status = publication_transaction.begin(publication_layers);
  if (status) status = publication_transaction.revise_trial(8U);
  if (status) status = publication_transaction.revise_trial(9U);
  semantic_species[0U].view.revision =
      publication_transaction.trial_revision(8U);
  semantic_species[1U].view.revision =
      publication_transaction.trial_revision(9U);
  const std::array<ConstFieldView, 2U> publication_species_views{{
      as_const(semantic_species[0U].view),
      as_const(semantic_species[1U].view)}};
  passed &= expect(static_cast<bool>(status), rank,
                   "periodic three-species publication transaction reserves semantic revisions");
  if (!all_true(passed, world)) return false;

  PisoExactThermodynamicCandidateView publication_thermodynamic =
      thermodynamic;
  publication_thermodynamic.independent_species = {
      publication_species_views.data(), publication_species_views.size()};
  PisoFrozenMomentumExactCandidateInput publication_baseline_input =
      replay_input;
  publication_baseline_input.thermodynamic = publication_thermodynamic;
  PisoFrozenMomentumExactCandidateCertificate publication_baseline;
  status = fixture.coupler.certify_frozen_momentum_exact_baseline(
      authority, baseline_pressure_stage, baseline_velocity_stage,
      baseline_flux_stage, publication_baseline_input, restored_reductions,
      publication_baseline);
  PressureEnergyStationaryCertificate publication_stationary;
  if (status)
    status = fixture.coupler.certify_frozen_momentum_stationary(
        authority, publication_baseline, 0.0, 0.0, restored_reductions,
        publication_stationary);
  PisoStateCorrectionCertificate corrected_one;
  if (status)
    status = fixture.coupler.commit_frozen_momentum_stationary_trial_state(
        authority, publication_baseline, publication_stationary,
        baseline_state, fixture.trial_flux, restored_reductions,
        corrected_one);
  passed &= expect(static_cast<bool>(status) && corrected_one.valid() &&
                       corrected_one.corrector == 1U,
                   rank,
                   "periodic three-species exact baseline commits C1 through stationary authority");
  if (!all_true(passed, world)) return false;

  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        accepted.view.unchecked(cell, 0U) =
            fixture.density.view.unchecked(cell, 0U);
        previous.view.unchecked(cell, 0U) =
            fixture.density.view.unchecked(cell, 0U);
      }

  PisoIntermediateInput c2_input = intermediate_input;
  c2_input.corrector = 2U;
  c2_input.temporal_reference = {};
  c2_input.committed_face_history = {};
  c2_input.prior_corrector = corrected_one.state;
  c2_input.pressure_reference = corrected_one.output_pressure_reference;
  c2_input.thermophysical_boundary.binding.pressure_reference =
      publication_thermodynamic.closed_gauge.next_pressure_reference;
  c2_input.thermophysical_boundary.binding.independent_species = {
      publication_species_views.data(), publication_species_views.size()};
  c2_input.trial_velocity = as_const(fixture.velocity.view);
  c2_input.density = fixture.density.view;
  ConstFaceFluxView c2_trial = as_const(fixture.trial_flux);
  c2_trial.revision = corrected_one.state;
  c2_input.trial_flux = c2_trial;
  ++c2_input.momentum.state;
  PisoIntermediateCertificate c2_intermediate;
  status = fixture.coupler.refresh(c2_input, c2_intermediate);
  const PressureCorrectionInput c2_pressure_input{
      c2_intermediate,
      c2_input.pressure_reference,
      as_const(fixture.density.view),
      as_const(accepted.view),
      as_const(previous.view),
      as_const(drho_dp.view),
      bdf,
      c2_input.momentum.time,
      c2_input.momentum.geometry,
      c2_input.numeric_boundary};
  PressureCorrectionCertificate c2_pressure;
  if (status)
    status = fixture.coupler.assemble_pressure_system(
        c2_pressure_input, {diagonal.view, rhs.view}, c2_pressure);
  PisoFrozenMomentumStageAuthority c2_authority;
  if (status)
    status = fixture.coupler.make_frozen_momentum_stage_authority(
        c2_intermediate, c2_pressure, c2_authority);
  passed &= expect(static_cast<bool>(status) && c2_authority.valid() &&
                       c2_authority.corrector() == 2U,
                   rank,
                   "periodic three-species corrected C1 lineage issues C2 frozen authority");
  if (!all_true(passed, world)) return false;

  ++scaled_dp.view.revision;
  ++scaled_dh.view.revision;
  ++candidate_pressure.view.revision;
  ++candidate_enthalpy.view.revision;
  ++candidate_velocity.view.revision;
  ++candidate_density.view.revision;
  ++candidate_temperature.view.revision;
  ++candidate_flux.revision;
  fill(scaled_dh, 0.0);
  fill(candidate_enthalpy, base_enthalpy_value);
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x)
        candidate_temperature.view.unchecked({x, y, z}, 0U) =
            base_temperature.view.unchecked({x, y, z}, 0U);
  restore_candidate_density();
  PisoFrozenMomentumPressureStageCertificate c2_pressure_stage;
  status = fixture.coupler.form_frozen_momentum_scaled_pressure(
      c2_authority, as_const(raw_dp.view), candidate_halo, 0.0,
      scaled_dp.view, c2_pressure_stage);
  std::array<FieldView, 1U> c2_halo_fields{scaled_dp.view};
  HaloTicket c2_halo_ticket;
  if (status)
    status = candidate_halo.begin(
        7520U, {c2_halo_fields.data(), c2_halo_fields.size()},
        c2_halo_ticket);
  if (status)
    status = candidate_halo.finish(
        c2_halo_ticket,
        {c2_halo_fields.data(), c2_halo_fields.size()});
  scaled_dp.view = c2_halo_fields[0U];
  PisoFrozenMomentumVelocityStageCertificate c2_velocity_stage;
  if (status)
    status = fixture.coupler.stage_frozen_momentum_velocity(
        c2_authority, c2_pressure_stage, candidate_halo,
        as_const(scaled_dp.view), candidate_velocity.view,
        c2_velocity_stage);
  PisoFrozenMomentumFluxStageCertificate c2_flux_stage;
  if (status)
    status = fixture.coupler.stage_frozen_momentum_flux(
        c2_authority, c2_velocity_stage, as_const(candidate_density.view),
        candidate_flux, c2_flux_stage);

  PisoExactEosClosureIdentity c2_closure = closure;
  c2_closure.pressure_reference =
      corrected_one.output_pressure_reference.pressure_reference;
  c2_closure.pressure_state =
      make_piso_field_revision_identity(as_const(base_pressure.view));
  c2_closure.pressure_correction =
      make_piso_field_revision_identity(as_const(scaled_dp.view));
  c2_closure.enthalpy_state =
      make_piso_field_revision_identity(as_const(base_enthalpy.view));
  c2_closure.enthalpy_correction =
      make_piso_field_revision_identity(as_const(scaled_dh.view));
  c2_closure.candidate_enthalpy =
      make_piso_field_revision_identity(as_const(candidate_enthalpy.view));
  c2_closure.candidate_density =
      make_piso_field_revision_identity(as_const(candidate_density.view));
  c2_closure.candidate_temperature =
      make_piso_field_revision_identity(as_const(candidate_temperature.view));
  c2_closure.closure = 7521U;
  PisoExactThermodynamicCandidateView c2_thermodynamic;
  c2_thermodynamic.enthalpy = as_const(candidate_enthalpy.view);
  c2_thermodynamic.density = as_const(candidate_density.view);
  c2_thermodynamic.temperature = as_const(candidate_temperature.view);
  c2_thermodynamic.closure = c2_closure;
  c2_thermodynamic.pressure_compressibility = as_const(drho_dp.view);
  c2_thermodynamic.independent_species = {
      publication_species_views.data(), publication_species_views.size()};
  if (status)
    status = fixture.prepare_closed_gauge(
        c2_pressure, corrected_one.output_pressure_reference,
        publication_thermodynamic.closed_gauge.next_pressure_reference,
        as_const(base_pressure.view), as_const(scaled_dp.view),
        c2_thermodynamic.pressure_compressibility, c2_closure.closure,
        restored_reductions, c2_thermodynamic.closed_gauge);
  if (status) {
    for (std::int32_t z = 0; z < cells.z; ++z)
      for (std::int32_t y = 0; y < cells.y; ++y)
        for (std::int32_t x = 0; x < cells.x; ++x)
          candidate_pressure.view.unchecked({x, y, z}, 0U) =
              base_pressure.view.unchecked({x, y, z}, 0U) +
              scaled_dp.view.unchecked({x, y, z}, 0U) -
              c2_thermodynamic.closed_gauge.shift;
  }
  const PisoFrozenMomentumExactCandidateInput c2_exact_input{
      as_const(raw_dh.view),
      as_const(scaled_dp.view),
      as_const(scaled_dh.view),
      baseline_state,
      as_const(candidate_pressure.view),
      c2_thermodynamic,
      as_const(candidate_velocity.view),
      as_const(candidate_flux),
      {},
      {replay_raw_views.data(), replay_raw_views.size()}};
  PisoFrozenMomentumExactCandidateCertificate c2_exact;
  if (status)
    status = fixture.coupler.certify_frozen_momentum_exact_baseline(
        c2_authority, c2_pressure_stage, c2_velocity_stage, c2_flux_stage,
        c2_exact_input, restored_reductions, c2_exact);
  PressureEnergyStationaryCertificate c2_stationary;
  if (status)
    status = fixture.coupler.certify_frozen_momentum_stationary(
        c2_authority, c2_exact, 0.0, 0.0, restored_reductions,
        c2_stationary);
  PendingFaceFluxView pending;
  if (status)
    status = publication_flux_writer.begin_pending(
        publication_transaction, publication_flux_storage, pending);
  PisoStateCorrectionCertificate corrected_two;
  if (status)
    status = fixture.coupler.commit_frozen_momentum_stationary_pending_state(
        c2_authority, c2_exact, c2_stationary, baseline_state, pending,
        restored_reductions, corrected_two);
  passed &= expect(static_cast<bool>(status) && corrected_two.valid() &&
                       corrected_two.corrector == 2U && pending.valid(),
                   rank,
                   "periodic three-species exact C2 establishes opaque pending flux authority");
  if (!all_true(passed, world)) return false;

  PisoTerminalAuditInput audit;
  audit.correction = corrected_two;
  audit.pressure_reference = corrected_two.output_pressure_reference;
  audit.density = as_const(fixture.density.view);
  audit.eos_density = as_const(fixture.density.view);
  audit.density_accepted = as_const(accepted.view);
  audit.density_previous = as_const(previous.view);
  audit.pressure_perturbation = as_const(base_pressure.view);
  audit.drho_dp_h_y = c2_thermodynamic.pressure_compressibility;
  audit.bdf = bdf;
  audit.step_dt = time_step_for_bdf(bdf);
  audit.convective_cfl_limit = 1.0e6;
  audit.thermophysical_boundary.binding.independent_species = {
      publication_species_views.data(), publication_species_views.size()};
  double local_closed_mass = 0.0;
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        local_closed_mass += fixture.density.view.unchecked(cell, 0U) *
                             cell_volume(fixture, cell);
      }
  MPI_Allreduce(&local_closed_mass, &audit.closed_mass_target, 1, MPI_DOUBLE,
                MPI_SUM, world);

  if (size > 1) {
    ReductionEngine foreign_reductions;
    const Status compiled_foreign = ReductionEngine::compile(
        MPI_COMM_SELF, ReductionMode::mpi_allreduce, 7U,
        foreign_reductions);
    PisoAttemptReport foreign_reduction_report;
    PisoTerminalCertificate foreign_reduction_terminal;
    const Status foreign_reduction_status =
        compiled_foreign
            ? fixture.coupler.audit_pending_final(
                  audit, pending, foreign_reductions,
                  foreign_reduction_report, foreign_reduction_terminal)
            : compiled_foreign;
    passed &= expect(
        foreign_reduction_status.code == StatusCode::invalid_plan &&
            !foreign_reduction_terminal.valid() && pending.valid(),
        rank,
        "terminal audit rejects a ReductionEngine from a foreign communicator before its first reduction");
  }

  std::array<OwnedField, 2U> foreign_species{
      make_field(8U, cells, 1U, 0U, publication_species_views[0U].revision,
                 8522U),
      make_field(9U, cells, 1U, 0U, publication_species_views[1U].revision,
                 8523U)};
  fill(foreign_species[0U], 0.2);
  fill(foreign_species[1U], 0.3);
  const std::array<ConstFieldView, 2U> foreign_species_views{{
      as_const(foreign_species[0U].view),
      as_const(foreign_species[1U].view)}};
  PisoTerminalAuditInput foreign_audit = audit;
  foreign_audit.thermophysical_boundary.binding.independent_species = {
      foreign_species_views.data(), foreign_species_views.size()};
  PisoAttemptReport foreign_report;
  PisoTerminalCertificate foreign_terminal;
  const Status foreign_composition_status =
      fixture.coupler.audit_pending_final(
      foreign_audit, pending, restored_reductions, foreign_report,
      foreign_terminal);
  passed &= expect(foreign_composition_status.code ==
                           StatusCode::invalid_plan &&
                       !foreign_terminal.valid() && pending.valid(),
                   rank,
                   "equal-valued foreign periodic species storage cannot replace the exact commit binding at terminal audit");

  const double semantic_before =
      semantic_species[0U].view.unchecked(poison_cell, 0U);
  if (rank == size - 1)
    semantic_species[0U].view.unchecked(poison_cell, 0U) =
        std::nextafter(semantic_before,
                       std::numeric_limits<double>::infinity());
  PisoAttemptReport washed_report;
  PisoTerminalCertificate washed_terminal;
  const Status washed_status = fixture.coupler.audit_pending_final(
      foreign_audit, pending, restored_reductions, washed_report,
      washed_terminal);
  passed &= expect(washed_status.code == StatusCode::invalid_plan &&
                       !washed_terminal.valid() && pending.valid(),
                   rank,
                   "foreign old-valued storage cannot wash an unrevisioned live periodic composition mutation");
  semantic_species[0U].view.unchecked(poison_cell, 0U) = semantic_before;

  PisoAttemptReport terminal_report;
  PisoTerminalCertificate terminal;
  status = fixture.coupler.audit_pending_final(
      audit, pending, restored_reductions, terminal_report, terminal);
  passed &= expect(static_cast<bool>(status) && terminal.valid(), rank,
                   "restored semantic periodic composition re-audits against the exact C2 binding");
  if (!all_true(passed, world)) return false;

  const std::array<RevisionDependency, 2U> publication_dependencies{{
      {AttemptTransaction::field_revision_source(8U),
       publication_species_views[0U].revision},
      {AttemptTransaction::field_revision_source(9U),
       publication_species_views[1U].revision}}};
  if (rank == size - 1)
    semantic_species[1U].view.unchecked(poison_cell, 0U) =
        std::nextafter(
            semantic_species[1U].view.unchecked(poison_cell, 0U),
            std::numeric_limits<double>::infinity());
  const Status stale_terminal_status = fixture.coupler.publish_pending_final(
      terminal,
      {publication_dependencies.data(), publication_dependencies.size()},
      {publication_species_views.data(), publication_species_views.size()},
      restored_reductions, publication_flux_writer, pending);
  passed &= expect(stale_terminal_status.code == StatusCode::invalid_plan &&
                       pending.valid(),
                   rank,
                   "post-terminal unrevisioned live species mutation rejects publication without consuming pending flux");
  semantic_species[1U].view.unchecked(poison_cell, 0U) = 0.3;

  status = fixture.coupler.audit_pending_final(
      audit, pending, restored_reductions, terminal_report, terminal);
  if (status)
    status = fixture.coupler.publish_pending_final(
        terminal,
        {publication_dependencies.data(), publication_dependencies.size()},
        {publication_species_views.data(), publication_species_views.size()},
        restored_reductions, publication_flux_writer, pending);
  if (status)
    status = publication_transaction.collective_finish(world, Status{});
  ConstFaceFluxView committed_flux;
  passed &= expect(
      static_cast<bool>(status) && publication_transaction.committed() &&
          publication_flux_writer.committed(publication_flux_storage,
                                             committed_flux) &&
          committed_flux.revision == corrected_two.face_flux,
      rank,
      "restored periodic composition terminal certificate publishes the exact final flux collectively");
  return all_true(passed, world);
}

bool test_pressure_correction_fraction_red(MPI_Comm world, int rank) {
  int size = 0;
  MPI_Comm_size(world, &size);
  if (size != 2) return true;

  // Keep this pressure-fraction probe independent from the collective
  // refresh-failure fixture.  The highest rank owns one finite correction that
  // would deplete rho by exactly two without a common affine pressure scale.
  PeriodicPisoFixture fixture;
  const bool initialized = fixture.initialize(9, world);
  if (!all_true(initialized, world)) return false;
  const Int3 cells = fixture.patch.cells;
  fill(fixture.density, 1.0);
  fill(fixture.velocity, 0.0);
  fill(fixture.momentum_diagonal, 1.0);
  fill(fixture.momentum_rhs, 0.0);

  const BdfCoefficients bdf{10.0, -15.0, 5.0, 2U};
  PisoIntermediateCertificate intermediate;
  Status status = fixture.coupler.refresh(
      fixture.intermediate_input(bdf, 7401U), intermediate);
  if (!all_true(static_cast<bool>(status) && intermediate.valid(), world))
    return false;

  OwnedField accepted = make_field(50U, cells, 1U, 0U, 7402U, 8402U);
  OwnedField previous = make_field(51U, cells, 1U, 0U, 7403U, 8403U);
  OwnedField drho_dp = make_field(52U, cells, 1U, 0U, 7404U, 8404U);
  OwnedField pressure_diagonal =
      make_field(53U, cells, 1U, 0U, 7405U, 8405U);
  OwnedField pressure_rhs = make_field(54U, cells, 1U, 0U, 7406U, 8406U);
  fill(accepted, 0.9);
  fill(previous, 0.7);
  fill(drho_dp, 1.0);
  const PressureCorrectionSystemView pressure_system{
      pressure_diagonal.view, pressure_rhs.view};
  PisoIntermediateInput intermediate_input =
      fixture.intermediate_input(bdf, 7401U);
  PressureCorrectionInput pressure_input;
  pressure_input.intermediate = intermediate;
  pressure_input.pressure_reference = intermediate_input.pressure_reference;
  pressure_input.density_trial = as_const(fixture.density.view);
  pressure_input.density_accepted = as_const(accepted.view);
  pressure_input.density_previous = as_const(previous.view);
  pressure_input.drho_dp_h_y = as_const(drho_dp.view);
  pressure_input.bdf = bdf;
  pressure_input.time = intermediate_input.momentum.time;
  pressure_input.geometry = intermediate_input.momentum.geometry;
  pressure_input.numeric_boundary = intermediate_input.numeric_boundary;
  PressureCorrectionCertificate pressure;
  status = fixture.coupler.assemble_pressure_system(
      pressure_input, pressure_system, pressure);
  if (!all_true(static_cast<bool>(status) && pressure.valid(), world))
    return false;

  OwnedField pressure_perturbation =
      make_field(2U, cells, 1U, 0U, 7407U, 8407U);
  OwnedField correction = make_field(90U, cells, 1U, 1U, 7408U, 8408U);
  fill(pressure_perturbation, 0.0);
  fill(correction, 0.0);
  std::array<std::int32_t, 3U> bad_global_values{};
  if (rank == size - 1) {
    bad_global_values = {fixture.patch.begin.x, fixture.patch.begin.y,
                         fixture.patch.begin.z};
  }
  if (MPI_Bcast(bad_global_values.data(),
                static_cast<int>(bad_global_values.size()), MPI_INT,
                size - 1, world) != MPI_SUCCESS)
    return false;
  const Int3 bad_global{bad_global_values[0U], bad_global_values[1U],
                        bad_global_values[2U]};
  constexpr double raw_pressure_correction = -2.0;
  if (rank == size - 1)
    correction.view.unchecked({0, 0, 0}, 0U) = raw_pressure_correction;

  const Int3 global_cells = fixture.geometry.global_cells();
  const auto wrap = [](std::int32_t value, std::int32_t extent) noexcept {
    value %= extent;
    return value < 0 ? value + extent : value;
  };
  const auto raw_value = [&](Int3 global) noexcept {
    return wrap(global.x, global_cells.x) == wrap(bad_global.x,
                                                  global_cells.x) &&
                   wrap(global.y, global_cells.y) == wrap(bad_global.y,
                                                           global_cells.y) &&
                   wrap(global.z, global_cells.z) == wrap(bad_global.z,
                                                           global_cells.z)
               ? raw_pressure_correction
               : 0.0;
  };
  const auto axis_global = [](Int3 value, CartesianAxis axis) noexcept {
    return axis == CartesianAxis::x
               ? value.x
               : (axis == CartesianAxis::y ? value.y : value.z);
  };
  const auto axis_width = [&](CartesianAxis axis,
                              std::int32_t global) noexcept {
    const std::int32_t wrapped =
        wrap(global, axis == CartesianAxis::x
                       ? global_cells.x
                       : (axis == CartesianAxis::y ? global_cells.y
                                                   : global_cells.z));
    return axis == CartesianAxis::x
               ? fixture.geometry.x().widths().data[wrapped]
               : (axis == CartesianAxis::y
                      ? fixture.geometry.y().widths().data[wrapped]
                      : fixture.geometry.z().widths().data[wrapped]);
  };
  const auto raw_gradient = [&](Int3 local,
                                CartesianAxis axis) noexcept {
    const Int3 global{fixture.patch.begin.x + local.x,
                      fixture.patch.begin.y + local.y,
                      fixture.patch.begin.z + local.z};
    Int3 minus = global;
    Int3 plus = global;
    if (axis == CartesianAxis::x) {
      --minus.x;
      ++plus.x;
    } else if (axis == CartesianAxis::y) {
      --minus.y;
      ++plus.y;
    } else {
      --minus.z;
      ++plus.z;
    }
    return 0.5 / axis_width(axis, axis_global(global, axis)) *
           (raw_value(plus) - raw_value(minus));
  };
  const auto raw_face_jump = [&](CartesianAxis axis, Int3 face) noexcept {
    Int3 right{fixture.patch.begin.x + face.x,
               fixture.patch.begin.y + face.y,
               fixture.patch.begin.z + face.z};
    Int3 left = right;
    if (axis == CartesianAxis::x)
      --left.x;
    else if (axis == CartesianAxis::y)
      --left.y;
    else
      --left.z;
    return raw_value(right) - raw_value(left);
  };

  const HaloRuntimeCounters halo_before = fixture.correction_halo.runtime_counters();
  const std::vector<double> density_before = fixture.density.storage;
  const std::vector<double> velocity_before = fixture.velocity.storage;
  const std::vector<double> pressure_before = pressure_perturbation.storage;
  const std::vector<double> flux_before =
      face_values(as_const(fixture.trial_flux));
  ReductionEngine reductions;
  const Status reduction_status = ReductionEngine::compile(
      world, ReductionMode::mpi_allreduce, 1U, reductions);
  if (!all_true(static_cast<bool>(reduction_status), world)) return false;
  PisoStateCorrectionCertificate certificate{};
  status = fixture.coupler.correct_trial_state(
      pressure, correction.view,
      {fixture.velocity.view, pressure_perturbation.view,
       fixture.density.view, as_const(drho_dp.view)},
      fixture.trial_flux, reductions, certificate);
  const HaloRuntimeCounters halo_after = fixture.correction_halo.runtime_counters();
  const bool affine_closed_fail_closed =
      status.code == StatusCode::invalid_plan && !certificate.valid() &&
      fixture.density.storage == density_before &&
      fixture.velocity.storage == velocity_before &&
      pressure_perturbation.storage == pressure_before &&
      face_values(as_const(fixture.trial_flux)) == flux_before &&
      halo_after.begin_calls == halo_before.begin_calls &&
      halo_after.finish_calls == halo_before.finish_calls &&
      (rank != size - 1 ||
       correction.view.unchecked({0, 0, 0}, 0U) ==
           raw_pressure_correction);
  return all_true(
      expect(affine_closed_fail_closed, rank,
             "closed-mass affine correction cannot bypass exact EOS/gauge authority and rejects atomically"),
      world);

}

bool test_exact_eos_correction_collective_transaction(MPI_Comm world,
                                                      int rank) {
  int size = 0;
  MPI_Comm_size(world, &size);

  PeriodicPisoFixture fixture;
  bool passed = expect(fixture.initialize(9, world), rank,
                       "exact-EOS correction fixture compiles");
  if (!all_true(passed, world)) return false;
  const Int3 cells = fixture.patch.cells;
  fill(fixture.density, 1.0);
  fill(fixture.velocity, 0.0);
  fill(fixture.momentum_diagonal, 1.0);
  fill(fixture.momentum_rhs, 0.0);

  const BdfCoefficients bdf{10.0, -15.0, 5.0, 2U};
  OwnedField accepted = make_field(50U, cells, 1U, 0U, 8941U, 9941U);
  OwnedField previous = make_field(51U, cells, 1U, 0U, 8942U, 9942U);
  OwnedField drho_dp = make_field(6U, cells, 1U, 0U, 8943U, 9943U);
  OwnedField pressure_diagonal =
      make_field(53U, cells, 1U, 0U, 8944U, 9944U);
  OwnedField pressure_rhs =
      make_field(54U, cells, 1U, 0U, 8945U, 9945U);
  OwnedField pressure = make_field(2U, cells, 1U, 1U, 8946U, 9946U);
  OwnedField correction = make_field(90U, cells, 1U, 1U, 8947U, 9947U);
  OwnedField enthalpy = make_field(3U, cells, 1U, 1U, 8948U, 9948U);
  OwnedField temperature = make_field(4U, cells, 1U, 1U, 8949U, 9949U);
  OwnedField dh = make_field(91U, cells, 1U, 0U, 8950U, 9950U);
  OwnedField candidate_enthalpy =
      make_field(92U, cells, 1U, 0U, 8951U, 9951U);
  OwnedField candidate_density =
      make_field(93U, cells, 1U, 0U, 8952U, 9952U);
  OwnedField candidate_temperature =
      make_field(94U, cells, 1U, 0U, 8953U, 9953U);
  OwnedField candidate_compressibility =
      make_field(95U, cells, 1U, 0U, 8954U, 9954U);
  fill(accepted, 1.0);
  fill(previous, 1.0);
  fill(drho_dp, 1.0);
  fill(pressure, 0.0);
  fill(correction, 0.0);
  fill(enthalpy, 0.0);
  fill(temperature, 300.0);
  fill(dh, 2.0);
  fill(candidate_enthalpy, 0.0);
  fill(candidate_density, 0.0);
  fill(candidate_temperature, 0.0);
  fill(candidate_compressibility, 0.0);

  const PressureCorrectionSystemView pressure_system{
      pressure_diagonal.view, pressure_rhs.view};
  PisoIntermediateInput intermediate_input;
  PressureCorrectionInput pressure_input;
  PressureCorrectionCertificate pressure_certificate;
  constexpr double pressure_reference_value = 101325.0;
  ReductionEngine reductions;
  passed &= expect(
      static_cast<bool>(ReductionEngine::compile(
          world, ReductionMode::mpi_allreduce, 2U, reductions)),
      rank, "exact-EOS correction reductions compile");
  if (!all_true(passed, world)) return false;
  auto populate_exact_candidate = [&](bool limiting_density) {
    for (std::int32_t z = 0; z < cells.z; ++z) {
      for (std::int32_t y = 0; y < cells.y; ++y) {
        for (std::int32_t x = 0; x < cells.x; ++x) {
          const Int3 cell{x, y, z};
          const bool limiting = limiting_density && rank == size - 1 &&
                                x == 0 && y == 0 && z == 0;
          const double target_density = limiting ? 0.5 : 1.0;
          double reference_h = 0.0;
          double cp = 0.0;
          double gas_constant = 0.0;
          Status local = fixture.thermodynamics.mixture_enthalpy(
              300.0, {}, reference_h, cp, gas_constant);
          const double absolute_pressure =
              pressure_reference_value + pressure.view.unchecked(cell, 0U) +
              correction.view.unchecked(cell, 0U);
          const double target_temperature =
              absolute_pressure / (gas_constant * target_density);
          double target_enthalpy = 0.0;
          if (local) {
            local = fixture.thermodynamics.mixture_enthalpy(
                target_temperature, {}, target_enthalpy, cp, gas_constant);
          }
          ThermoState replayed;
          if (local) {
            local = fixture.thermodynamics.evaluate_from_reference_pressure(
                pressure_reference_value,
                pressure.view.unchecked(cell, 0U) +
                    correction.view.unchecked(cell, 0U),
                target_enthalpy, {}, {}, replayed, target_temperature);
          }
          if (!local) return local;
          candidate_enthalpy.view.unchecked(cell, 0U) = target_enthalpy;
          enthalpy.view.unchecked(cell, 0U) =
              target_enthalpy - dh.view.unchecked(cell, 0U);
          candidate_density.view.unchecked(cell, 0U) = replayed.rho;
          candidate_temperature.view.unchecked(cell, 0U) =
              replayed.temperature;
          candidate_compressibility.view.unchecked(cell, 0U) =
              replayed.drho_dp_hY;
        }
      }
    }
    return Status{};
  };
  passed &= expect(static_cast<bool>(populate_exact_candidate(false)), rank,
                   "distributed exact-EOS fixture starts from a physical candidate");
  if (!all_true(passed, world)) return false;
  auto rebuild_c1 = [&](RevisionToken time,
                        double absolute_pressure_reference =
                            pressure_reference_value) {
    intermediate_input = fixture.intermediate_input(bdf, time);
    intermediate_input.thermophysical_boundary.binding.pressure_reference =
        absolute_pressure_reference;
    PisoIntermediateCertificate intermediate;
    Status status = fixture.coupler.refresh(intermediate_input, intermediate);
    if (!status) return status;
    pressure_input.intermediate = intermediate;
    pressure_input.pressure_reference =
        intermediate_input.pressure_reference;
    pressure_input.density_trial = as_const(fixture.density.view);
    pressure_input.density_accepted = as_const(accepted.view);
    pressure_input.density_previous = as_const(previous.view);
    pressure_input.drho_dp_h_y = as_const(drho_dp.view);
    pressure_input.bdf = bdf;
    pressure_input.time = intermediate_input.momentum.time;
    pressure_input.geometry = intermediate_input.momentum.geometry;
    pressure_input.numeric_boundary = intermediate_input.numeric_boundary;
    return fixture.coupler.assemble_pressure_system(
        pressure_input, pressure_system, pressure_certificate);
  };
  auto exact_candidate = [&]() {
    PisoExactEosClosureIdentity closure;
    closure.thermodynamics =
        pressure_input.pressure_reference.thermodynamics;
    closure.pressure_reference =
        pressure_input.pressure_reference.pressure_reference;
    closure.composition = exact_composition_identity_for_test(
        closure.thermodynamics, {}, cells);
    closure.pressure_state =
        make_piso_field_revision_identity(as_const(pressure.view));
    closure.pressure_correction =
        make_piso_field_revision_identity(as_const(correction.view));
    closure.enthalpy_state =
        make_piso_field_revision_identity(as_const(enthalpy.view));
    closure.enthalpy_correction =
        make_piso_field_revision_identity(as_const(dh.view));
    closure.candidate_enthalpy = make_piso_field_revision_identity(
        as_const(candidate_enthalpy.view));
    closure.candidate_density = make_piso_field_revision_identity(
        as_const(candidate_density.view));
    closure.candidate_temperature = make_piso_field_revision_identity(
        as_const(candidate_temperature.view));
    closure.closure = 8955U;
    PisoExactThermodynamicCandidateView candidate;
    candidate.enthalpy = as_const(candidate_enthalpy.view);
    candidate.density = as_const(candidate_density.view);
    candidate.temperature = as_const(candidate_temperature.view);
    candidate.closure = closure;
    candidate.pressure_compressibility =
        as_const(candidate_compressibility.view);
    candidate.pressure_compressibility.field = drho_dp.view.field;
    return candidate;
  };
  auto certified_candidate = [&](Status& gauge_status) {
    PisoExactThermodynamicCandidateView candidate = exact_candidate();
    gauge_status = fixture.prepare_closed_gauge(
        pressure_certificate, pressure_input.pressure_reference,
        pressure_reference_value, as_const(pressure.view),
        as_const(correction.view), candidate.pressure_compressibility,
        candidate.closure.closure, reductions, candidate.closed_gauge);
    return candidate;
  };
  const PisoCoupledStateView state{fixture.velocity.view, pressure.view,
                                   enthalpy.view, fixture.density.view,
                                   temperature.view};

  auto unchanged = [&](const std::vector<double>& velocity_before,
                       const std::vector<double>& pressure_before,
                       const std::vector<double>& enthalpy_before,
                       const std::vector<double>& density_before,
                       const std::vector<double>& temperature_before,
                       const std::vector<double>& flux_before) {
    return fixture.velocity.storage == velocity_before &&
           pressure.storage == pressure_before &&
           enthalpy.storage == enthalpy_before &&
           fixture.density.storage == density_before &&
           temperature.storage == temperature_before &&
           face_values(as_const(fixture.trial_flux)) == flux_before;
  };
  auto reject = [&](PisoExactThermodynamicCandidateView candidate,
                    StatusCode expected) {
    const std::vector<double> velocity_before = fixture.velocity.storage;
    const std::vector<double> pressure_before = pressure.storage;
    const std::vector<double> enthalpy_before = enthalpy.storage;
    const std::vector<double> density_before = fixture.density.storage;
    const std::vector<double> temperature_before = temperature.storage;
    const std::vector<double> flux_before =
        face_values(as_const(fixture.trial_flux));
    PisoStateCorrectionCertificate certificate;
    const Status status = fixture.coupler.correct_coupled_trial_state(
        pressure_certificate, correction.view, as_const(dh.view), state,
        candidate, fixture.trial_flux, reductions, certificate);
    const bool atomic = unchanged(velocity_before, pressure_before,
                                  enthalpy_before, density_before,
                                  temperature_before, flux_before);
    return status.code == expected && !certificate.valid() && atomic;
  };

  Status status = rebuild_c1(8956U);
  Status gauge_status;
  PisoExactThermodynamicCandidateView candidate =
      certified_candidate(gauge_status);
  if (rank == size - 1)
    candidate_density.view.unchecked({0, 0, 0}, 0U) = -1.0;
  // The gauge authority was prepared before the EOS candidate mutation, so
  // the coupler must still reject the now-stale transaction collectively.
  passed &= expect(static_cast<bool>(status) && static_cast<bool>(gauge_status) &&
                       reject(candidate,
                              StatusCode::numerical_failure),
                   rank,
                   "one-rank non-positive exact-EOS candidate rejects on every rank without partial publication");
  if (!all_true(passed, world)) return false;

  passed &= expect(static_cast<bool>(populate_exact_candidate(false)), rank,
                   "exact-EOS fixture restores its physical density after poison");
  if (!all_true(passed, world)) return false;
  fill(pressure, 0.0);
  if (rank == size - 1) {
    pressure.view.unchecked({0, 0, 0}, 0U) =
        -pressure_reference_value;
  }
  fill(correction, 0.0);
  status = rebuild_c1(8957U);
  candidate = exact_candidate();
  gauge_status = fixture.prepare_closed_gauge(
      pressure_certificate, pressure_input.pressure_reference,
      pressure_reference_value, as_const(pressure.view),
      as_const(correction.view), candidate.pressure_compressibility,
      candidate.closure.closure, reductions, candidate.closed_gauge);
  passed &= expect(static_cast<bool>(status) &&
                       gauge_status.code == StatusCode::numerical_failure,
                   rank,
                   "one-rank zero exact-EOS absolute pressure rejects during collective gauge preparation");
  if (!all_true(passed, world)) return false;

  fill(pressure, 0.0);
  if (rank == size - 1) {
    pressure.view.unchecked({0, 0, 0}, 0U) =
        std::numeric_limits<double>::max();
  }
  status = rebuild_c1(8958U, std::numeric_limits<double>::max());
  candidate = exact_candidate();
  gauge_status = fixture.prepare_closed_gauge(
      pressure_certificate, pressure_input.pressure_reference,
      std::numeric_limits<double>::max(), as_const(pressure.view),
      as_const(correction.view), candidate.pressure_compressibility,
      candidate.closure.closure, reductions, candidate.closed_gauge);
  passed &= expect(static_cast<bool>(status) &&
                       gauge_status.code == StatusCode::numerical_failure,
                   rank,
                   "one-rank infinite exact-EOS absolute pressure rejects during collective gauge preparation");
  if (!all_true(passed, world)) return false;

  fill(pressure, 0.0);
  passed &= expect(static_cast<bool>(populate_exact_candidate(false)), rank,
                   "exact-EOS fixture restores its physical state after pressure poison");
  if (!all_true(passed, world)) return false;
  status = rebuild_c1(8959U);
  PisoExactThermodynamicCandidateView stale = certified_candidate(gauge_status);
  if (rank == size - 1) ++stale.temperature.revision;
  passed &= expect(static_cast<bool>(status) && static_cast<bool>(gauge_status) &&
                       reject(stale, StatusCode::invalid_plan),
                   rank,
                   "one-rank stale exact-EOS revision reaches collective rejection without partial publication");
  if (!all_true(passed, world)) return false;

  status = rebuild_c1(8960U);
  PisoExactThermodynamicCandidateView stale_gauge =
      certified_candidate(gauge_status);
  if (rank == size - 1)
    ++stale_gauge.closed_gauge.collective_transaction;
  passed &= expect(
      static_cast<bool>(status) && static_cast<bool>(gauge_status) &&
          reject(stale_gauge, StatusCode::invalid_plan),
      rank,
      "one-rank stale closed-gauge certificate rejects collectively without partial publication");
  if (!all_true(passed, world)) return false;

  fill(correction, 0.0);
  if (rank == size - 1) {
    correction.view.unchecked({0, 0, 0}, 0U) = -2.0;
  }
  passed &= expect(static_cast<bool>(populate_exact_candidate(true)), rank,
                   "positive exact candidate closes EOS at the non-affine limiting cell");
  if (!all_true(passed, world)) return false;
  status = rebuild_c1(8961U);
  PisoStateCorrectionCertificate committed;
  candidate = certified_candidate(gauge_status);
  if (status) {
    status = gauge_status;
  }
  if (status) {
    status = fixture.coupler.correct_coupled_trial_state(
        pressure_certificate, correction.view, as_const(dh.view), state,
        candidate, fixture.trial_flux, reductions, committed);
  }
  const bool exact_positive =
      rank != size - 1 ||
      (correction.view.unchecked({0, 0, 0}, 0U) == -2.0 &&
       std::abs(fixture.density.view.unchecked({0, 0, 0}, 0U) - 0.5) <=
           64.0 * std::numeric_limits<double>::epsilon());
  passed &= expect(
      static_cast<bool>(status) && committed.valid() &&
          committed.closure == PisoStateClosure::exact_eos &&
          committed.enthalpy == enthalpy.view.revision &&
          committed.temperature == temperature.view.revision &&
          exact_positive,
      rank,
      "exact-EOS coupled correction keeps alpha one and publishes positive rho when affine rho would be negative");
  return all_true(passed, world);
}

bool test_c2_applied_candidate_audit_red(MPI_Comm world, int rank) {
  int size = 0;
  MPI_Comm_size(world, &size);
  if (size != 2) return true;

  // This fixture makes the C2 linear solution a known, rank-selective
  // correction.  The limiting cell has rho_trial=psi=1 and rho_n=rho_nm1=.1;
  // every other cell is at the accepted density, so the future applied
  // candidate has zero BDF density defect and zero corrected face flux.
  PeriodicPisoFixture fixture;
  bool passed = expect(fixture.initialize(9, world), rank,
                       "C2 applied-candidate audit fixture compiles");
  if (!all_true(passed, world)) return false;
  const Int3 cells = fixture.patch.cells;
  const Int3 global_cells = fixture.geometry.global_cells();
  constexpr double raw_pressure_correction = -2.0;
  const double expected_scale = std::nextafter(0.45, 0.0);
  const double gauge_offset =
      -(expected_scale * raw_pressure_correction) /
      static_cast<double>(global_cells.x * global_cells.y * global_cells.z);
  fill(fixture.density, 1.0);
  fill(fixture.velocity, 0.0);
  fill(fixture.momentum_diagonal, 1.0);
  fill(fixture.momentum_rhs, 0.0);

  const BdfCoefficients bdf{10.0, -15.0, 5.0, 2U};
  PisoIntermediateInput intermediate_input =
      fixture.intermediate_input(bdf, 8601U);
  PisoIntermediateCertificate intermediate_one;
  Status status = fixture.coupler.refresh(intermediate_input,
                                          intermediate_one);
  passed &= expect(static_cast<bool>(status) && intermediate_one.valid(), rank,
                   "C2 applied-candidate corrector-one refreshes");
  if (!all_true(passed, world)) return false;

  std::array<std::int32_t, 3U> bad_global_values{};
  if (rank == size - 1) {
    bad_global_values = {fixture.patch.begin.x, fixture.patch.begin.y,
                         fixture.patch.begin.z};
  }
  if (MPI_Bcast(bad_global_values.data(),
                static_cast<int>(bad_global_values.size()), MPI_INT,
                size - 1, world) != MPI_SUCCESS)
    return false;
  const Int3 bad_global{bad_global_values[0U], bad_global_values[1U],
                        bad_global_values[2U]};
  const auto same_global = [](Int3 left, Int3 right) noexcept {
    return left.x == right.x && left.y == right.y && left.z == right.z;
  };

  OwnedField accepted = make_field(50U, cells, 1U, 0U, 8602U, 9602U);
  OwnedField previous = make_field(51U, cells, 1U, 0U, 8603U, 9603U);
  OwnedField drho_dp = make_field(52U, cells, 1U, 0U, 8604U, 9604U);
  OwnedField pressure_diagonal =
      make_field(53U, cells, 1U, 0U, 8605U, 9605U);
  OwnedField pressure_rhs = make_field(54U, cells, 1U, 0U, 8606U, 9606U);
  fill(accepted, 1.0);
  fill(previous, 1.0);
  fill(drho_dp, 1.0);
  for (std::int32_t z = 0; z < cells.z; ++z) {
    for (std::int32_t y = 0; y < cells.y; ++y) {
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        const Int3 global{fixture.patch.begin.x + x,
                          fixture.patch.begin.y + y,
                          fixture.patch.begin.z + z};
        if (same_global(global, bad_global)) {
          accepted.view.unchecked(cell, 0U) = 0.1;
          previous.view.unchecked(cell, 0U) = 0.1;
        }
      }
    }
  }
  const PressureCorrectionSystemView pressure_system{
      pressure_diagonal.view, pressure_rhs.view};
  PressureCorrectionInput pressure_input;
  pressure_input.intermediate = intermediate_one;
  pressure_input.pressure_reference = intermediate_input.pressure_reference;
  pressure_input.density_trial = as_const(fixture.density.view);
  pressure_input.density_accepted = as_const(accepted.view);
  pressure_input.density_previous = as_const(previous.view);
  pressure_input.drho_dp_h_y = as_const(drho_dp.view);
  pressure_input.bdf = bdf;
  pressure_input.time = intermediate_input.momentum.time;
  pressure_input.geometry = intermediate_input.momentum.geometry;
  pressure_input.numeric_boundary = intermediate_input.numeric_boundary;
  PressureCorrectionCertificate pressure_one;
  status = fixture.coupler.assemble_pressure_system(
      pressure_input, pressure_system, pressure_one);
  passed &= expect(static_cast<bool>(status) && pressure_one.valid(), rank,
                   "C2 applied-candidate corrector-one pressure system assembles");
  if (!all_true(passed, world)) return false;

  // The first pressure solve is deliberately zero-RHS.  The pressure
  // certificate remains the public lifecycle token; the zeroing only keeps
  // the C1 state unchanged while this test later controls the C2 RHS.
  fill(pressure_rhs, 0.0);
  LinearIdentity identity_one{8607U, 8608U, 8609U, 0U, 8611U};
  const MgCoefficientIdentity coefficient_one{8612U, pressure_one.state, 0.0};
  LinearResources linear;
  passed &= expect(initialize_linear_resources(
                       world, fixture, pressure_one, pressure_system,
                       identity_one, coefficient_one, linear),
                   rank, "C2 applied-candidate linear resources compile");
  if (!all_true(passed, world)) return false;

  constexpr FieldId correction_field = 90U;
  OwnedField correction =
      make_field(correction_field, cells, 1U, 1U, 8613U, 9613U);
  fill(correction, 0.0);
  OwnedField pressure = make_field(2U, cells, 1U, 1U, 8614U, 9614U);
  fill(pressure, gauge_offset);
  PisoPressureSolveEpoch epoch;
  status = epoch.begin(fixture.piso);
  if (status) {
    status = epoch.solve(fixture.piso, 1U, pressure_one, identity_one,
                         coefficient_one, fixture.coupler,
                         linear.linear_operator, linear.mg, pressure_system,
                         correction.view, linear.workspace, linear.reductions,
                         nullptr, &linear.mg_counters);
  }
  passed &= expect(static_cast<bool>(status) && epoch.solve_calls() == 1U,
                   rank, "C2 applied-candidate C1 solve establishes lifecycle");
  if (!all_true(passed, world)) return false;

  FieldView velocity = fixture.velocity.view;
  FieldView density = fixture.density.view;
  ++velocity.revision;
  ++density.revision;
  const std::vector<double> c1_velocity_before = fixture.velocity.storage;
  const std::vector<double> c1_pressure_before = pressure.storage;
  const std::vector<double> c1_density_before = fixture.density.storage;
  const std::vector<double> c1_flux_before =
      face_values(as_const(fixture.trial_flux));
  PisoStateCorrectionCertificate corrected_one;
  status = fixture.coupler.correct_trial_state(
      pressure_one, correction.view,
      {velocity, pressure.view, density, as_const(drho_dp.view)},
      fixture.trial_flux, linear.reductions, corrected_one);
  passed &= expect(
      status.code == StatusCode::invalid_plan && !corrected_one.valid() &&
          fixture.velocity.storage == c1_velocity_before &&
          pressure.storage == c1_pressure_before &&
          fixture.density.storage == c1_density_before &&
          face_values(as_const(fixture.trial_flux)) == c1_flux_before,
      rank,
      "closed-mass affine C1 is retired before the applied-candidate lifecycle");
  return all_true(passed, world);

}

bool test_c1_normalized_committed_face_history_across_decomposition_seams(
    MPI_Comm world, int rank) {
  int size = 0;
  MPI_Comm_size(world, &size);
  constexpr std::int32_t global_cell_count = 8;
  PeriodicPisoFixture fixture;
  bool passed = expect(fixture.initialize(global_cell_count, world), rank,
                       "normalized-history periodic fixture initializes");
  if (!all_true(passed, world)) return false;

  // This is variable-step BDF2 with dt / dt_previous = 1/2 and dt = 4/3.
  // Choosing a0=1 makes the face-history oracle independent of a time-scale
  // multiplication, while all coefficients and expected values remain exact
  // binary fractions.
  constexpr BdfCoefficients bdf{1.0, -1.125, 0.125, 2U};
  constexpr std::array<double, 3U> previous_flux{1.0, 2.0, 3.0};
  constexpr std::array<double, 3U> accepted_flux{4.0, 6.0, 8.0};
  constexpr std::array<double, 3U> paired_flux{9.0, 10.0, 11.0};
  constexpr std::array<double, 3U> expected_normalized_history{
      4.375, 6.5, 8.625};
  constexpr double inverse_face_area =
      static_cast<double>(global_cell_count * global_cell_count);
  const std::array<double, 3U> previous_velocity{
      previous_flux[0U] * inverse_face_area,
      previous_flux[1U] * inverse_face_area,
      previous_flux[2U] * inverse_face_area};
  const std::array<double, 3U> accepted_velocity{
      accepted_flux[0U] * inverse_face_area,
      accepted_flux[1U] * inverse_face_area,
      accepted_flux[2U] * inverse_face_area};
  passed &= expect(
      fixture.commit_uniform_flux_history(
          world, previous_velocity, accepted_velocity, 8811U, 8821U),
      rank, "nonzero committed face-history publishes through its writer");
  if (!all_true(passed, world)) return false;

  fill(fixture.density, 1.0);
  fill(fixture.velocity, 0.0);
  fill(fixture.momentum_rhs, 0.0);
  const Int3 cells = fixture.patch.cells;
  for (std::int32_t z = 0; z < cells.z; ++z) {
    for (std::int32_t y = 0; y < cells.y; ++y) {
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        const double diagonal = bdf.a0 * cell_volume(fixture, cell);
        for (std::uint8_t component = 0U; component < 3U; ++component)
          fixture.momentum_diagonal.view.unchecked(cell, component) =
              diagonal;
      }
    }
  }

  const FaceFieldView temporal_faces[]{fixture.phi_h_by_a.x,
                                       fixture.phi_h_by_a.y,
                                       fixture.phi_h_by_a.z};
  const FaceFieldView paired_faces[]{fixture.trial_flux.x,
                                     fixture.trial_flux.y,
                                     fixture.trial_flux.z};
  for (std::size_t axis_index = 0U; axis_index < 3U; ++axis_index) {
    const FaceFieldView temporal = temporal_faces[axis_index];
    const FaceFieldView paired = paired_faces[axis_index];
    for (std::int32_t z = 0; z < paired.extents.z; ++z) {
      for (std::int32_t y = 0; y < paired.extents.y; ++y) {
        for (std::int32_t x = 0; x < paired.extents.x; ++x) {
          const Int3 face{x, y, z};
          temporal.unchecked(face) = 0.0;
          paired.unchecked(face) = paired_flux[axis_index];
        }
      }
    }
  }

  PisoIntermediateInput input = fixture.intermediate_input(bdf, 8831U);
  PisoIntermediateCertificate certificate;
  const Status status = fixture.coupler.refresh(input, certificate);
  passed &= expect(static_cast<bool>(status) && certificate.valid() &&
                       certificate.corrector == 1U,
                   rank,
                   "C1 refresh accepts certified nonzero face histories");
  if (!all_true(passed, world)) return false;

  const auto axis_value = [](Int3 value, std::size_t axis) noexcept {
    return axis == 0U ? value.x : (axis == 1U ? value.y : value.z);
  };
  const ConstFaceFieldView output_faces[]{as_const(fixture.phi_h_by_a.x),
                                          as_const(fixture.phi_h_by_a.y),
                                          as_const(fixture.phi_h_by_a.z)};
  std::uint64_t local_seam_faces = 0U;
  double local_max_error = 0.0;
  bool local_exact = true;
  for (std::size_t axis_index = 0U; axis_index < 3U; ++axis_index) {
    const ConstFaceFieldView output = output_faces[axis_index];
    const std::int32_t local_extent = axis_value(cells, axis_index);
    const std::int32_t patch_begin =
        axis_value(fixture.patch.begin, axis_index);
    for (std::int32_t z = 0; z < output.extents.z; ++z) {
      for (std::int32_t y = 0; y < output.extents.y; ++y) {
        for (std::int32_t x = 0; x < output.extents.x; ++x) {
          const Int3 face{x, y, z};
          const double value = output.unchecked(face);
          const double error =
              std::abs(value - expected_normalized_history[axis_index]);
          local_max_error = std::max(local_max_error, error);
          local_exact = local_exact &&
                        value == expected_normalized_history[axis_index];
          const std::int32_t normal = axis_value(face, axis_index);
          const bool low_decomposition_seam =
              normal == 0 && patch_begin > 0;
          const bool high_decomposition_seam =
              normal == local_extent &&
              patch_begin + local_extent < global_cell_count;
          if (low_decomposition_seam || high_decomposition_seam)
            ++local_seam_faces;
        }
      }
    }
  }

  double global_max_error = 0.0;
  std::uint64_t global_seam_faces = 0U;
  const int error_reduce = MPI_Allreduce(&local_max_error, &global_max_error,
                                         1, MPI_DOUBLE, MPI_MAX, world);
  const int seam_reduce = MPI_Allreduce(&local_seam_faces,
                                        &global_seam_faces, 1, MPI_UINT64_T,
                                        MPI_SUM, world);
  const Int3 process_grid = fixture.patch.process_grid;
  const std::uint64_t expected_global_seam_faces =
      2U * static_cast<std::uint64_t>(global_cell_count) *
      static_cast<std::uint64_t>(global_cell_count) *
      static_cast<std::uint64_t>((process_grid.x - 1) +
                                 (process_grid.y - 1) +
                                 (process_grid.z - 1));
  const bool discriminates_paired_ex2 =
      expected_normalized_history[0U] != paired_flux[0U] &&
      expected_normalized_history[1U] != paired_flux[1U] &&
      expected_normalized_history[2U] != paired_flux[2U];
  passed &= expect(
      local_exact && error_reduce == MPI_SUCCESS &&
          seam_reduce == MPI_SUCCESS && global_max_error == 0.0 &&
          global_seam_faces == expected_global_seam_faces &&
          (size == 1 || local_seam_faces > 0U) &&
          discriminates_paired_ex2,
      rank,
      "C1 uses normalized BDF history exactly on every rank and decomposition seam");
  return all_true(passed, world);
}

bool test_full_distributed_piso(MPI_Comm world, int rank) {
  int size = 0;
  MPI_Comm_size(world, &size);
  PeriodicPisoFixture fixture;
  bool passed = expect(fixture.initialize(9, world), rank,
                       "non-divisible periodic PISO fixture compiles");
  if (!all_true(passed, world)) {
    return false;
  }
  const Int3 cells = fixture.patch.cells;
  const Int3 global_cells = fixture.geometry.global_cells();
  fill(fixture.momentum_diagonal, 0.0);
  fill(fixture.momentum_rhs, 0.0);
  constexpr std::array<double, 3U> component_diagonal{2.0, 4.0, 8.0};
  for (std::int32_t z = 0; z < cells.z; ++z) {
    for (std::int32_t y = 0; y < cells.y; ++y) {
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 local{x, y, z};
        const Int3 global{fixture.patch.begin.x + x,
                          fixture.patch.begin.y + y,
                          fixture.patch.begin.z + z};
        const std::array<double, 3U> h_by_a{
            periodic_value(global.x, global_cells.x, 0.0),
            periodic_value(global.y, global_cells.y, 0.31),
            periodic_value(global.z, global_cells.z, -0.27)};
        for (std::uint8_t component = 0U; component < 3U; ++component) {
          const double diagonal =
              cell_volume(fixture, local) * component_diagonal[component];
          fixture.momentum_diagonal.view.unchecked(local, component) =
              diagonal;
          fixture.momentum_rhs.view.unchecked(local, component) =
              diagonal * h_by_a[component];
          fixture.velocity.view.unchecked(local, component) =
              h_by_a[component];
        }
      }
    }
  }

  const BdfCoefficients bdf{10.0, -15.0, 5.0, 2U};
  constexpr double pressure_reference_value = 101325.0;
  PisoIntermediateInput intermediate_input =
      fixture.intermediate_input(bdf, 7201U);
  intermediate_input.thermophysical_boundary.binding.pressure_reference =
      pressure_reference_value;
  PisoIntermediateCertificate intermediate_one;
  Status status = fixture.coupler.refresh(intermediate_input,
                                          intermediate_one);
  passed &= expect(static_cast<bool>(status) && intermediate_one.valid(), rank,
                   "distributed corrector-one intermediates refresh");

  OwnedField accepted = make_field(50U, cells, 1U, 0U, 7202U, 8202U);
  OwnedField previous = make_field(51U, cells, 1U, 0U, 7203U, 8203U);
  OwnedField drho_dp = make_field(6U, cells, 1U, 0U, 7204U, 8204U);
  OwnedField pressure_diagonal =
      make_field(53U, cells, 1U, 0U, 7205U, 8205U);
  OwnedField pressure_rhs =
      make_field(54U, cells, 1U, 0U, 7206U, 8206U);
  fill(accepted, 0.9);
  fill(previous, 0.7);
  fill(drho_dp, 0.02);
  const PressureCorrectionSystemView pressure_system{
      pressure_diagonal.view, pressure_rhs.view};
  PressureCorrectionInput pressure_input;
  pressure_input.intermediate = intermediate_one;
  pressure_input.pressure_reference = intermediate_input.pressure_reference;
  pressure_input.density_trial = as_const(fixture.density.view);
  pressure_input.density_accepted = as_const(accepted.view);
  pressure_input.density_previous = as_const(previous.view);
  pressure_input.drho_dp_h_y = as_const(drho_dp.view);
  pressure_input.bdf = bdf;
  pressure_input.time = intermediate_input.momentum.time;
  pressure_input.geometry = intermediate_input.momentum.geometry;
  pressure_input.numeric_boundary = intermediate_input.numeric_boundary;
  PressureCorrectionCertificate pressure_one;
  status = fixture.coupler.assemble_pressure_system(
      pressure_input, pressure_system, pressure_one);
  passed &= expect(static_cast<bool>(status), rank,
                   "distributed corrector-one pressure system assembles");
  if (!all_true(passed, world)) {
    return false;
  }

  LinearIdentity identity_one{7301U, 7302U, 7303U, 0U, 7305U};
  const MgCoefficientIdentity coefficient_one{7306U, pressure_one.state, 0.0};
  LinearResources linear;
  passed &= expect(initialize_linear_resources(
                       world, fixture, pressure_one, pressure_system,
                       identity_one, coefficient_one, linear),
                   rank, "distributed pressure linear resources compile cold");
  if (!all_true(passed, world)) {
    return false;
  }

  constexpr FieldId correction_field = 90U;
  OwnedField correction = make_field(correction_field, cells, 1U, 1U, 7307U,
                                     8307U);
  fill(correction, 0.0);
  const std::vector<double> failed_c1_caller_before = correction.storage;
  PisoPressureSolveEpoch failed_c1_epoch;
  Status failed_c1_status = failed_c1_epoch.begin(fixture.piso);
  NthHaloFailureOperator failed_c1_operator(
      linear.linear_operator, detail::HaloFailurePoint::completion, size - 1,
      1U);
  if (failed_c1_status) {
    failed_c1_status = failed_c1_epoch.solve(
        fixture.piso, 1U, pressure_one, identity_one, coefficient_one,
        fixture.coupler, linear.linear_operator, failed_c1_operator, linear.mg,
        pressure_system, correction.view, linear.workspace, linear.reductions,
        nullptr, &linear.mg_counters);
  }
  detail::clear_halo_failure_for_test();
  const Status failed_c1_stale_projection =
      linear.workspace.recycle_begin_projection_for_test(
          cells, identity_one.fingerprint);
  const bool failed_c1_contract =
      failed_c1_status.code == StatusCode::mpi_failure &&
          failed_c1_status.detail ==
              detail::halo_detail_completion_failure &&
          failed_c1_epoch.solve_calls() == 1U &&
          correction.storage == failed_c1_caller_before &&
          linear.workspace.recycle_correction_count_for_test() == 0U &&
          failed_c1_stale_projection.code == StatusCode::invalid_plan &&
          !linear.operator_halo.active() && linear.operator_halo.ready();
  passed &= expect(
      failed_c1_contract, rank,
      "real C1 Pressure/Halo failure records provenance without exposing a session or caller publication");
  if (!all_true(passed, world)) return false;

  PisoPressureSolveEpoch epoch;
  status = epoch.begin(fixture.piso);
  if (status) {
    status = epoch.solve(fixture.piso, 1U, pressure_one, identity_one,
                         coefficient_one, fixture.coupler,
                         linear.linear_operator, linear.mg, pressure_system,
                         correction.view, linear.workspace, linear.reductions,
                         nullptr, &linear.mg_counters);
  }
  if (!status) {
    std::cerr << "rank " << rank << " solve1 status="
              << static_cast<unsigned>(status.code)
              << " detail=" << status.detail
              << " calls=" << static_cast<unsigned>(epoch.solve_calls())
              << '\n';
  }
  passed &= expect(static_cast<bool>(status) && epoch.solve_calls() == 1U, rank,
                   "distributed corrector one retries the same resources and performs one FGMRES+MG solve");
  if (!all_true(passed, world)) {
    return false;
  }
  const std::size_t captured_cycle_count =
      linear.workspace.recycle_correction_count_for_test();
  std::array<std::vector<double>, kLinearRecycleMaximumDirections>
      captured_cycle_directions;
  const std::size_t local_cell_count =
      static_cast<std::size_t>(cells.x) * cells.y * cells.z;
  bool captured_cycle_copy =
      captured_cycle_count > 0U &&
      captured_cycle_count <= kLinearRecycleMaximumDirections;
  for (std::size_t direction = 0U;
       direction < captured_cycle_count && captured_cycle_copy; ++direction) {
    const ConstFieldView captured =
        linear.workspace.recycle_correction_for_test(direction, cells);
    captured_cycle_copy = captured.base != nullptr;
    captured_cycle_directions[direction].resize(local_cell_count);
    for (std::int32_t z = 0; z < cells.z && captured_cycle_copy; ++z) {
      for (std::int32_t y = 0; y < cells.y; ++y) {
        for (std::int32_t x = 0; x < cells.x; ++x) {
          const std::size_t local =
              (static_cast<std::size_t>(z) * cells.y + y) * cells.x + x;
          captured_cycle_directions[direction][local] =
              captured.unchecked({x, y, z}, 0U);
        }
      }
    }
  }
  passed &= expect(captured_cycle_copy && captured_cycle_count == 3U, rank,
                   "distributed C1 corrections are snapshotted for independent C2 failure/retry probes");
  if (!all_true(passed, world)) return false;

  OwnedField pressure = make_field(2U, cells, 1U, 1U, 7308U, 8308U);
  fill(pressure, 0.0);
  FaceFluxStorage trial_flux_storage;
  FaceFluxView trial_flux;
  passed &= expect(FaceFluxStorage::allocate_workspace(
                       cells, 1U, trial_flux_storage) &&
                       trial_flux_storage.workspace_view(0U, 7309U,
                                                         trial_flux),
                   rank, "distributed trial flux allocates cold");
  FieldView velocity = fixture.velocity.view;
  FieldView density = fixture.density.view;
  ++velocity.revision;
  ++density.revision;
  OwnedField enthalpy = make_field(3U, cells, 1U, 1U, 7317U, 8317U);
  OwnedField temperature = make_field(4U, cells, 1U, 1U, 7318U, 8318U);
  OwnedField enthalpy_correction =
      make_field(91U, cells, 1U, 0U, 7319U, 8319U);
  OwnedField candidate_enthalpy =
      make_field(92U, cells, 1U, 0U, 7320U, 8320U);
  OwnedField candidate_density =
      make_field(93U, cells, 1U, 0U, 7321U, 8321U);
  OwnedField candidate_temperature =
      make_field(94U, cells, 1U, 0U, 7322U, 8322U);
  OwnedField candidate_compressibility =
      make_field(99U, cells, 1U, 0U, 7323U, 8323U);
  fill(enthalpy, 0.0);
  fill(temperature, 300.0);
  fill(enthalpy_correction, 0.0);
  fill(candidate_enthalpy, 0.0);
  fill(candidate_temperature, 0.0);
  fill(candidate_compressibility, 0.0);
  double exact_reference_h = 0.0;
  double exact_cp = 0.0;
  double exact_gas_constant = 0.0;
  Status exact_profile = fixture.thermodynamics.mixture_enthalpy(
      300.0, {}, exact_reference_h, exact_cp, exact_gas_constant);
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        const double target_density =
            density.unchecked(cell, 0U) +
            drho_dp.view.unchecked(cell, 0U) *
                correction.view.unchecked(cell, 0U);
        const double absolute_pressure =
            pressure_reference_value + pressure.view.unchecked(cell, 0U) +
            correction.view.unchecked(cell, 0U);
        if (!exact_profile || !(target_density > 0.0) ||
            !(absolute_pressure > 0.0)) {
          exact_profile = {StatusCode::numerical_failure, 1U};
          continue;
        }
        const double target_temperature =
            absolute_pressure / (exact_gas_constant * target_density);
        double target_enthalpy = 0.0;
        double candidate_cp = 0.0;
        double candidate_gas_constant = 0.0;
        exact_profile = fixture.thermodynamics.mixture_enthalpy(
            target_temperature, {}, target_enthalpy, candidate_cp,
            candidate_gas_constant);
        if (!exact_profile) continue;
        ThermoState replayed;
        exact_profile = fixture.thermodynamics.evaluate_from_reference_pressure(
            pressure_reference_value,
            pressure.view.unchecked(cell, 0U) +
                correction.view.unchecked(cell, 0U),
            target_enthalpy, {},
            {velocity.unchecked(cell, 0U), velocity.unchecked(cell, 1U),
             velocity.unchecked(cell, 2U)},
            replayed, target_temperature);
        if (!exact_profile) continue;
        enthalpy.view.unchecked(cell, 0U) = target_enthalpy;
        candidate_enthalpy.view.unchecked(cell, 0U) = target_enthalpy;
        candidate_density.view.unchecked(cell, 0U) = target_density;
        candidate_temperature.view.unchecked(cell, 0U) =
            replayed.temperature;
        candidate_compressibility.view.unchecked(cell, 0U) =
            replayed.drho_dp_hY;
      }
  passed &= expect(static_cast<bool>(exact_profile), rank,
                   "distributed C1 candidate is replayed from physical p/h/Y/U");
  PisoExactEosClosureIdentity c1_closure;
  c1_closure.thermodynamics =
      pressure_input.pressure_reference.thermodynamics;
  c1_closure.pressure_reference =
      pressure_input.pressure_reference.pressure_reference;
  c1_closure.composition = exact_composition_identity_for_test(
      c1_closure.thermodynamics, {}, cells);
  c1_closure.pressure_state =
      make_piso_field_revision_identity(as_const(pressure.view));
  c1_closure.pressure_correction =
      make_piso_field_revision_identity(as_const(correction.view));
  c1_closure.enthalpy_state =
      make_piso_field_revision_identity(as_const(enthalpy.view));
  c1_closure.enthalpy_correction = make_piso_field_revision_identity(
      as_const(enthalpy_correction.view));
  c1_closure.candidate_enthalpy = make_piso_field_revision_identity(
      as_const(candidate_enthalpy.view));
  c1_closure.candidate_density = make_piso_field_revision_identity(
      as_const(candidate_density.view));
  c1_closure.candidate_temperature = make_piso_field_revision_identity(
      as_const(candidate_temperature.view));
  c1_closure.closure = 8324U;
  PisoExactThermodynamicCandidateView c1_candidate;
  c1_candidate.enthalpy = as_const(candidate_enthalpy.view);
  c1_candidate.density = as_const(candidate_density.view);
  c1_candidate.temperature = as_const(candidate_temperature.view);
  c1_candidate.closure = c1_closure;
  c1_candidate.pressure_compressibility =
      as_const(candidate_compressibility.view);
  c1_candidate.pressure_compressibility.field = drho_dp.view.field;
  status = fixture.prepare_closed_gauge(
      pressure_one, pressure_input.pressure_reference,
      pressure_reference_value, as_const(pressure.view),
      as_const(correction.view), c1_candidate.pressure_compressibility,
      c1_candidate.closure.closure, linear.reductions,
      c1_candidate.closed_gauge);
  passed &= expect(static_cast<bool>(status), rank,
                   "distributed exact C1 gauge transaction prepares");
  PisoStateCorrectionCertificate corrected_one;
  if (status) {
    status = fixture.coupler.correct_coupled_trial_state(
        pressure_one, correction.view, as_const(enthalpy_correction.view),
        {velocity, pressure.view, enthalpy.view, density, temperature.view},
        c1_candidate, trial_flux, linear.reductions, corrected_one);
  }
  passed &= expect(static_cast<bool>(status) && corrected_one.valid(), rank,
                   "distributed exact corrector one updates all five fields and flux");

  intermediate_input.corrector = 2U;
  intermediate_input.temporal_reference = {};
  intermediate_input.committed_face_history = {};
  intermediate_input.prior_corrector = corrected_one.state;
  intermediate_input.pressure_reference =
      corrected_one.output_pressure_reference;
  intermediate_input.thermophysical_boundary.binding.pressure_reference =
      c1_candidate.closed_gauge.next_pressure_reference;
  intermediate_input.trial_velocity = as_const(velocity);
  intermediate_input.density = density;
  ConstFaceFluxView corrected_trial_flux = as_const(trial_flux);
  corrected_trial_flux.revision = corrected_one.state;
  intermediate_input.trial_flux = corrected_trial_flux;
  ++intermediate_input.momentum.state;
  PisoIntermediateCertificate intermediate_two;
  status = fixture.coupler.refresh(intermediate_input, intermediate_two);
  passed &= expect(static_cast<bool>(status) &&
                       intermediate_two.corrector == 2U,
                   rank, "distributed corrector two refreshes current trial state");
  pressure_input.intermediate = intermediate_two;
  pressure_input.pressure_reference = intermediate_input.pressure_reference;
  pressure_input.density_trial = as_const(density);
  PressureCorrectionCertificate pressure_two;
  status = fixture.coupler.assemble_pressure_system(
      pressure_input, pressure_system, pressure_two);
  passed &= expect(static_cast<bool>(status), rank,
                   "distributed corrector-two pressure system assembles");
  if (!all_true(passed, world)) {
    return false;
  }

  LinearIdentity identity_two = identity_one;
  identity_two.numeric = 7312U;
  identity_two.hierarchy = 7313U;
  identity_two.fingerprint = 7315U;
  const MgCoefficientIdentity coefficient_two{7316U, pressure_two.state, 0.0};
  fill(correction, 0.0);
  status = epoch.solve(fixture.piso, 2U, pressure_two, identity_two,
                       coefficient_two, fixture.coupler,
                       linear.linear_operator, linear.mg, pressure_system,
                       correction.view, linear.workspace, linear.reductions,
                       nullptr, &linear.mg_counters);
  if (!status) {
    std::cerr << "rank " << rank << " solve2 status="
              << static_cast<unsigned>(status.code)
              << " detail=" << status.detail
              << " calls=" << static_cast<unsigned>(epoch.solve_calls())
              << '\n';
  }
  passed &= expect(static_cast<bool>(status) && epoch.solve_calls() == 2U, rank,
                   "distributed corrector two performs one refreshed solve");
  if (!all_true(passed, world)) {
    return false;
  }
  PisoAttemptReport solved_pressure_report;
  const Status solved_pressure_observed = epoch.observe(solved_pressure_report);
  passed &= expect(
      static_cast<bool>(solved_pressure_observed) &&
          solved_pressure_report.pressure_solve_calls == 2U &&
          solved_pressure_report.pressure[1U].recycle_projection_accepted,
      rank, "distributed solved pressure state is observable before downstream publication");

  // Exercise current-A2 Halo failure and retry while the exact C2 system and
  // Native-MG numerics are still the ones just solved above.  Later state and
  // final-flux work is intentionally downstream of this gate.
  passed &= test_real_pressure_halo_failure_provenance(
      world, rank, cells, as_const(pressure_system.rhs), identity_two,
      fixture.piso.pressure_solve(), linear.linear_operator, linear,
      captured_cycle_directions, captured_cycle_count, true);
  passed &= test_real_pressure_projection_failure_lifecycle(
      world, rank, cells, as_const(pressure_system.rhs), identity_two,
      fixture.piso.pressure_solve(), linear.linear_operator, linear,
      captured_cycle_directions, captured_cycle_count);
  IbmForceFixture ibm_fixture;
  const bool ibm_fixture_ready =
      ibm_fixture.initialize(world, global_cells.x) &&
      ibm_fixture.patch.cells.x == cells.x &&
      ibm_fixture.patch.cells.y == cells.y &&
      ibm_fixture.patch.cells.z == cells.z;
  passed &= expect(ibm_fixture_ready, rank,
                   "distributed IBM provenance fixture matches pressure partition");
  if (all_true(ibm_fixture_ready, world)) {
    OwnedCoefficientFace ibm_x = make_coefficient_face(
        CartesianAxis::x, cells, 8511U);
    OwnedCoefficientFace ibm_y = make_coefficient_face(
        CartesianAxis::y, cells, 8512U);
    OwnedCoefficientFace ibm_z = make_coefficient_face(
        CartesianAxis::z, cells, 8513U);
    IbmPressureOperator ibm_operator;
    const bool ibm_bound = static_cast<bool>(IbmPressureOperator::bind(
        linear.linear_operator, ibm_fixture.topology, ibm_fixture.boundary,
        as_const(ibm_x.view), as_const(ibm_y.view), as_const(ibm_z.view),
        ibm_fixture.geometry.topology_revision(), ibm_operator));
    passed &= expect(ibm_bound, rank,
                     "IBM pressure adapter binds over the real pressure operator");
    if (all_true(ibm_bound, world)) {
      passed &= test_real_pressure_halo_failure_provenance(
          world, rank, cells, as_const(pressure_system.rhs), identity_two,
          fixture.piso.pressure_solve(), ibm_operator, linear,
          captured_cycle_directions, captured_cycle_count, false);
    }
  }
  if (!all_true(passed, world)) return false;

  FieldRegistry transaction_registry;
  FieldSchema transaction_schema;
  FieldId dependency = 0U;
  const bool schema_ok =
      transaction_registry.declare_field("piso.final_state", 1U, 0U,
                                         dependency) &&
      transaction_registry.freeze(transaction_schema);
  const std::array requests{ArenaFieldRequest{
      dependency, {1, 1, 1}, {0U}, FieldLifetime::state_layer}};
  ArenaLayout layout;
  StateLayers layers;
  AttemptTransaction transaction;
  FaceFluxStorage final_flux_storage;
  FinalFaceFluxAuthority final_authority;
  FinalFaceFluxWriter final_writer;
  passed &= expect(
      schema_ok &&
          ArenaLayout::compile(transaction_schema,
                               {requests.data(), requests.size()}, layout) &&
          StateLayers::allocate(layout, layers) &&
          AttemptTransaction::create(layers.field_count(), 1U,
                                     layers.field_count(), transaction) &&
          FaceFluxStorage::allocate_final(cells, final_flux_storage) &&
          final_authority.claim(fixture.piso.pressure_stage(),
                                fixture.piso.final_flux_slot(), transaction,
                                final_writer) &&
          transaction.begin(layers) && transaction.revise_trial(dependency),
      rank, "distributed final-flux transaction prepares cold");
  const RevisionDependency final_dependency{
      AttemptTransaction::field_revision_source(dependency),
      transaction.trial_revision(dependency)};
  PendingFaceFluxView pending;
  passed &= expect(static_cast<bool>(final_writer.begin_pending(
                       transaction, final_flux_storage, pending)),
                   rank, "distributed corrector two acquires pending final flux");
  ++velocity.revision;
  ++density.revision;
  ++pressure.view.revision;
  ++enthalpy.view.revision;
  ++temperature.view.revision;
  OwnedField c2_enthalpy_correction =
      make_field(95U, cells, 1U, 0U, 8325U, 8425U);
  OwnedField c2_candidate_enthalpy =
      make_field(96U, cells, 1U, 0U, 8326U, 8426U);
  OwnedField c2_candidate_density =
      make_field(97U, cells, 1U, 0U, 8327U, 8427U);
  OwnedField c2_candidate_temperature =
      make_field(98U, cells, 1U, 0U, 8328U, 8428U);
  OwnedField c2_candidate_compressibility =
      make_field(100U, cells, 1U, 0U, 8330U, 8430U);
  fill(c2_enthalpy_correction, 0.0);
  fill(c2_candidate_enthalpy, 0.0);
  fill(c2_candidate_temperature, 0.0);
  fill(c2_candidate_compressibility, 0.0);
  exact_profile = {};
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        const double target_density =
            density.unchecked(cell, 0U) +
            drho_dp.view.unchecked(cell, 0U) *
                correction.view.unchecked(cell, 0U);
        const double absolute_pressure =
            c1_candidate.closed_gauge.next_pressure_reference +
            pressure.view.unchecked(cell, 0U) +
            correction.view.unchecked(cell, 0U);
        if (!(target_density > 0.0) || !(absolute_pressure > 0.0)) {
          exact_profile = {StatusCode::numerical_failure, 1U};
          continue;
        }
        const double target_temperature =
            absolute_pressure / (exact_gas_constant * target_density);
        double target_enthalpy = 0.0;
        double candidate_cp = 0.0;
        double candidate_gas_constant = 0.0;
        exact_profile = fixture.thermodynamics.mixture_enthalpy(
            target_temperature, {}, target_enthalpy, candidate_cp,
            candidate_gas_constant);
        if (!exact_profile) continue;
        ThermoState replayed;
        exact_profile = fixture.thermodynamics.evaluate_from_reference_pressure(
            c1_candidate.closed_gauge.next_pressure_reference,
            pressure.view.unchecked(cell, 0U) +
                correction.view.unchecked(cell, 0U),
            target_enthalpy, {},
            {velocity.unchecked(cell, 0U), velocity.unchecked(cell, 1U),
             velocity.unchecked(cell, 2U)},
            replayed, target_temperature);
        if (!exact_profile) continue;
        c2_candidate_enthalpy.view.unchecked(cell, 0U) = target_enthalpy;
        c2_enthalpy_correction.view.unchecked(cell, 0U) =
            target_enthalpy - enthalpy.view.unchecked(cell, 0U);
        c2_candidate_density.view.unchecked(cell, 0U) = target_density;
        c2_candidate_temperature.view.unchecked(cell, 0U) =
            replayed.temperature;
        c2_candidate_compressibility.view.unchecked(cell, 0U) =
            replayed.drho_dp_hY;
      }
  passed &= expect(static_cast<bool>(exact_profile), rank,
                   "distributed C2 candidate is replayed from physical p/h/Y/U");
  PisoExactEosClosureIdentity c2_closure;
  c2_closure.thermodynamics =
      pressure_input.pressure_reference.thermodynamics;
  c2_closure.pressure_reference =
      pressure_input.pressure_reference.pressure_reference;
  c2_closure.composition = c1_closure.composition;
  c2_closure.pressure_state =
      make_piso_field_revision_identity(as_const(pressure.view));
  c2_closure.pressure_correction =
      make_piso_field_revision_identity(as_const(correction.view));
  c2_closure.enthalpy_state =
      make_piso_field_revision_identity(as_const(enthalpy.view));
  c2_closure.enthalpy_correction = make_piso_field_revision_identity(
      as_const(c2_enthalpy_correction.view));
  c2_closure.candidate_enthalpy = make_piso_field_revision_identity(
      as_const(c2_candidate_enthalpy.view));
  c2_closure.candidate_density = make_piso_field_revision_identity(
      as_const(c2_candidate_density.view));
  c2_closure.candidate_temperature = make_piso_field_revision_identity(
      as_const(c2_candidate_temperature.view));
  c2_closure.closure = 8429U;
  PisoExactThermodynamicCandidateView c2_candidate;
  c2_candidate.enthalpy = as_const(c2_candidate_enthalpy.view);
  c2_candidate.density = as_const(c2_candidate_density.view);
  c2_candidate.temperature = as_const(c2_candidate_temperature.view);
  c2_candidate.closure = c2_closure;
  c2_candidate.pressure_compressibility =
      as_const(c2_candidate_compressibility.view);
  c2_candidate.pressure_compressibility.field = drho_dp.view.field;
  status = fixture.prepare_closed_gauge(
      pressure_two, pressure_input.pressure_reference,
      c1_candidate.closed_gauge.next_pressure_reference,
      as_const(pressure.view), as_const(correction.view),
      c2_candidate.pressure_compressibility, c2_candidate.closure.closure,
      linear.reductions, c2_candidate.closed_gauge);
  passed &= expect(static_cast<bool>(status), rank,
                   "distributed exact C2 gauge transaction prepares");
  PisoStateCorrectionCertificate corrected_two;
  if (status) {
    status = fixture.coupler.correct_coupled_pending_state(
        pressure_two, correction.view,
        as_const(c2_enthalpy_correction.view),
        {velocity, pressure.view, enthalpy.view, density, temperature.view},
        c2_candidate, pending, linear.reductions, corrected_two);
  }
  passed &= expect(static_cast<bool>(status) && corrected_two.valid(), rank,
                   "distributed exact corrector two writes five fields and pending pressure flux");
  if (!all_true(passed, world)) {
    return false;
  }

  PisoTerminalAuditInput audit;
  audit.correction = corrected_two;
  audit.pressure_reference = corrected_two.output_pressure_reference;
  audit.density = as_const(density);
  audit.eos_density = as_const(density);
  audit.density_accepted = as_const(accepted.view);
  audit.density_previous = as_const(previous.view);
  audit.pressure_perturbation = as_const(pressure.view);
  audit.drho_dp_h_y = c2_candidate.pressure_compressibility;
  audit.bdf = bdf;
  audit.step_dt = time_step_for_bdf(bdf);
  audit.convective_cfl_limit = 1.0e6;
  double local_mass = 0.0;
  for (std::int32_t z = 0; z < cells.z; ++z) {
    for (std::int32_t y = 0; y < cells.y; ++y) {
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        local_mass += density.unchecked(cell, 0U) * cell_volume(fixture, cell);
      }
    }
  }
  MPI_Allreduce(&local_mass, &audit.closed_mass_target, 1, MPI_DOUBLE,
                MPI_SUM, world);
  PisoAttemptReport report;
  PisoTerminalCertificate terminal;
  if (size > 1) {
    PisoTerminalAuditInput rank_local_invalid = audit;
    if (rank == size - 1)
      rank_local_invalid.boundary_closure_residual =
          std::numeric_limits<double>::quiet_NaN();
    PisoAttemptReport rejected_report;
    PisoTerminalCertificate rejected_terminal;
    const Status rejected = fixture.coupler.audit_pending_final(
        rank_local_invalid, pending, linear.reductions, rejected_report,
        rejected_terminal);
    passed &= expect(
        rejected.code == StatusCode::invalid_plan &&
            rejected.detail == 1503U && !rejected_terminal.valid(),
        rank,
        "rank-local terminal-input failure reaches collective consensus");
    if (!all_true(passed, world)) return false;
  }
  status = fixture.coupler.audit_pending_final(
      audit, pending, linear.reductions, report, terminal);
  if (!status) {
    std::cerr << "rank " << rank << " audit status="
              << static_cast<unsigned>(status.code)
              << " eos=" << report.eos_residual
              << " continuity=" << report.continuity_residual
              << " mass=" << report.closed_mass_residual
              << " gauge=" << report.gauge_residual << '\n';
  }
  passed &= expect(static_cast<bool>(status) && terminal.valid(), rank,
                   "distributed terminal EOS/continuity/mass/gauge audit passes");
  const std::array dependencies{final_dependency};
  if (status) {
    status = fixture.coupler.publish_pending_final(
        terminal, {dependencies.data(), dependencies.size()}, {},
        linear.reductions, final_writer, pending);
  }
  if (status) {
    status = transaction.collective_finish(world, Status{});
  }
  ConstFaceFluxView committed;
  passed &= expect(static_cast<bool>(status) &&
                       final_writer.committed(final_flux_storage, committed) &&
                       committed.revision == corrected_two.face_flux,
                   rank, "distributed pressure-equation flux commits collectively");
  if (status) {
    status = epoch.finalize(report);
  }
  passed &= expect(static_cast<bool>(status) &&
                       report.pressure_solve_calls == 2U &&
                       report.pressure[0U].status &&
                       report.pressure[1U].status &&
                       report.pressure[0U].recycle_offered_directions == 0U &&
                       !report.pressure[0U].recycle_projection_attempted &&
                       !report.pressure[0U].recycle_projection_accepted &&
                       report.pressure[0U].recycle_capture_cycle_attempts >=
                           report.pressure[0U].recycle_cycle_corrections &&
                       report.pressure[0U].recycle_capture_vector_passes ==
                           2U *
                               report.pressure[0U]
                                   .recycle_capture_cycle_attempts &&
                       report.pressure[0U].recycle_capture_reduction_calls ==
                           report.pressure[0U]
                               .recycle_capture_cycle_attempts &&
                       report.pressure[0U]
                               .recycle_capture_blocking_operations ==
                           2U *
                               report.pressure[0U]
                                   .recycle_capture_cycle_attempts &&
                       report.pressure[0U].recycle_cycle_corrections == 3U &&
                       report.pressure[0U].recycle_capture_cycle_attempts ==
                           3U &&
                       report.pressure[0U].recycle_capture_vector_passes ==
                           6U &&
                       report.pressure[0U].recycle_capture_reduction_calls ==
                           3U &&
                       report.pressure[0U]
                               .recycle_capture_blocking_operations == 6U &&
                       report.pressure[1U].recycle_offered_directions == 3U &&
                       report.pressure[1U].recycle_projection_attempted &&
                       report.pressure[1U].recycle_retained_directions == 3U &&
                       report.pressure[1U].recycle_operator_applies == 4U &&
                       report.pressure[1U].recycle_reduction_calls == 17U &&
                       report.pressure[1U].recycle_projection_accepted &&
                       report.pressure[1U].recycle_cycle_corrections == 0U &&
                       report.pressure[1U]
                               .recycle_capture_vector_passes == 0U &&
                       report.pressure[1U]
                               .recycle_capture_cycle_attempts == 0U &&
                       report.pressure[1U]
                               .recycle_capture_reduction_calls == 0U &&
                       report.pressure[1U]
                               .recycle_capture_blocking_operations == 0U &&
                       report.pressure[1U].recycle_projected_true_residual <
                           report.pressure[1U].initial_true_residual &&
                       report.final_flux_revision == committed.revision &&
                       linear.mg.counters().numeric_refreshes == 1U,
                   rank, "distributed report records exactly two solves and one numeric refresh");
  return all_true(passed, world);
}

bool test_candidate_boundary_finalizer_collective(MPI_Comm world, int rank) {
  int size = 0;
  MPI_Comm_size(world, &size);
  bool passed = true;

  CandidateBoundaryFixture fixture;
  CandidateBoundaryFixtureSpec velocity_spec;
  velocity_spec.inlet = CandidateBoundaryInlet::velocity;
  velocity_spec.inlet_velocity = 1.25;
  velocity_spec.multispecies = true;
  CandidateBoundaryScratch baseline;
  passed &= expect(
      fixture.initialize(world, velocity_spec) &&
          fixture.stage(0.0, 16.0, 4.0, 12000U, baseline) &&
          baseline.final_boundary.valid(),
      rank, "distributed real-chain velocity-inlet finalizer succeeds");
  if (!all_true(passed, world)) return false;

  PisoFrozenMomentumExactCandidateCertificate multispecies_exact_baseline;
  const Status multispecies_exact_status =
      fixture.coupler.certify_frozen_momentum_exact_baseline(
          fixture.authority, baseline.pressure_stage,
          baseline.velocity_stage, baseline.flux_stage,
          fixture.exact_input(baseline), fixture.reductions,
          multispecies_exact_baseline);
  passed &= expect(
      multispecies_exact_status && multispecies_exact_baseline.valid() &&
          baseline.independent_species.size() == 2U &&
          baseline.thermophysical_species_aliases.size() == 2U &&
          baseline.independent_species[0U].view.field == 230U &&
          baseline.independent_species[1U].view.field == 231U &&
          baseline.thermophysical_species_aliases[0U].field == 8U &&
          baseline.thermophysical_species_aliases[1U].field == 9U,
      rank,
      "distributed three-species/two-independent baseline binds distinct candidate and semantic scalar IDs into exact publication");
  if (!all_true(passed, world)) return false;

  const Int3 cells = fixture.patch.cells;
  const bool inlet_owner =
      fixture.local_face_owner(CartesianFace::x_min);
  const bool outlet_owner =
      fixture.local_face_owner(CartesianFace::x_max);
  const double inlet_final =
      baseline.final_flux.x.unchecked({0, 0, 0});
  const double inlet_mechanical =
      baseline.mechanical_flux.x.unchecked({0, 0, 0});
  const double outlet_final =
      baseline.final_flux.x.unchecked({cells.x, 0, 0});
  const double outlet_mechanical =
      baseline.mechanical_flux.x.unchecked({cells.x, 0, 0});
  passed &= expect(
      (inlet_owner ? (inlet_final > 0.0 && inlet_final != inlet_mechanical)
                   : inlet_final == inlet_mechanical) &&
          (!outlet_owner || outlet_final == outlet_mechanical) &&
          baseline.final_flux.x.unchecked({1, 0, 0}) ==
              baseline.mechanical_flux.x.unchecked({1, 0, 0}) &&
          baseline.final_flux.y.unchecked({0, 0, 0}) == 0.0 &&
          baseline.final_flux.z.unchecked({0, 0, 0}) == 0.0,
      rank,
      "only the physical owner rewrites inlet/outlet faces while partition/internal flux is bit-preserved");
  if (!all_true(passed, world)) return false;

  const std::vector<double> finalized_before_rejection =
      face_values(as_const(baseline.final_flux));
  PressureEnergyCandidateBoundaryFinalizeInput reference_poison =
      fixture.finalizer_input(baseline);
  if (rank == size - 1)
    ++reference_poison.pressure_reference.pressure_reference;
  FinalBoundaryFluxCertificate cleared = baseline.final_boundary;
  const Status reference_rejected = fixture.finalizer.finalize(
      reference_poison, fixture.reductions, cleared);
  passed &= expect(
      reference_rejected.code == StatusCode::invalid_plan &&
          !cleared.valid() &&
          face_values(as_const(baseline.final_flux)) ==
              finalized_before_rejection,
      rank,
      "one-rank pressure-reference poison rejects collectively with cleared certificate and zero flux write");

  HaloEngine foreign_state_halo;
  const std::array<HaloFieldSpec, 1U> foreign_state_contract{{
      {baseline.pressure.view.field, 1U, 1U}}};
  passed &= expect(
      static_cast<bool>(foreign_state_halo.reserve(
          world, fixture.patch,
          {foreign_state_contract.data(), foreign_state_contract.size()},
          fixture.boundary.halo_topology())),
      rank, "foreign candidate state halo fixture reserves collectively");
  PressureEnergyCandidateBoundaryFinalizeInput halo_poison =
      fixture.finalizer_input(baseline);
  if (rank == size - 1) halo_poison.state_halo = &foreign_state_halo;
  FinalBoundaryFluxCertificate halo_cleared = baseline.final_boundary;
  const Status halo_rejected = fixture.finalizer.finalize(
      halo_poison, fixture.reductions, halo_cleared);
  passed &= expect(
      halo_rejected.code == StatusCode::invalid_plan &&
          !halo_cleared.valid() &&
          face_values(as_const(baseline.final_flux)) ==
              finalized_before_rejection,
      rank,
      "one-rank foreign candidate-state halo contract rejects collectively without entering divergent finalizer collectives");

  PressureEnergyCandidateBoundaryFinalizeInput species_span_poison =
      fixture.finalizer_input(baseline);
  if (rank == size - 1 && species_span_poison.independent_species.size > 0U)
    --species_span_poison.independent_species.size;
  FinalBoundaryFluxCertificate species_span_cleared =
      baseline.final_boundary;
  const Status species_span_rejected = fixture.finalizer.finalize(
      species_span_poison, fixture.reductions, species_span_cleared);
  passed &= expect(
      species_span_rejected.code == StatusCode::invalid_plan &&
          !species_span_cleared.valid() &&
          face_values(as_const(baseline.final_flux)) ==
              finalized_before_rejection,
      rank,
      "one-rank short composition span rejects collectively before any species dereference or flux write");

  const Int3 poison_cell{0, 0, 0};
  const double velocity_before =
      baseline.velocity.view.unchecked(poison_cell, 0U);
  if (rank == size - 1)
    baseline.velocity.view.unchecked(poison_cell, 0U) =
        velocity_before + 0.125;
  FinalBoundaryFluxCertificate state_cleared = baseline.final_boundary;
  const Status state_rejected = fixture.finalizer.finalize(
      fixture.finalizer_input(baseline), fixture.reductions, state_cleared);
  passed &= expect(
      state_rejected.code == StatusCode::invalid_plan &&
          !state_cleared.valid() &&
          face_values(as_const(baseline.final_flux)) ==
              finalized_before_rejection,
      rank,
      "one-rank unrevisioned candidate-state poison rejects collectively with zero flux write");
  baseline.velocity.view.unchecked(poison_cell, 0U) = velocity_before;
  passed &= expect(
      static_cast<bool>(fixture.finalizer.finalize(
          fixture.finalizer_input(baseline), fixture.reductions,
          baseline.final_boundary)) &&
          baseline.final_boundary.valid(),
      rank, "failed distributed finalizer probes leave the stage replayable");
  if (!all_true(passed, world)) return false;

  const double species_before =
      baseline.independent_species[0U].view.unchecked(poison_cell, 0U);
  if (rank == size - 1)
    baseline.independent_species[0U].view.unchecked(poison_cell, 0U) =
        std::nextafter(species_before,
                       std::numeric_limits<double>::infinity());
  PisoFrozenMomentumExactCandidateCertificate poisoned_exact =
      multispecies_exact_baseline;
  const Status exact_species_rejected =
      fixture.coupler.certify_frozen_momentum_exact_baseline(
          fixture.authority, baseline.pressure_stage,
          baseline.velocity_stage, baseline.flux_stage,
          fixture.exact_input(baseline), fixture.reductions,
          poisoned_exact);
  passed &= expect(
      exact_species_rejected.code == StatusCode::invalid_plan &&
          !poisoned_exact.valid(),
      rank,
      "one-rank no-revision composition poison invalidates final-boundary replay collectively");
  baseline.independent_species[0U].view.unchecked(poison_cell, 0U) =
      species_before;
  if (!all_true(passed, world)) return false;

  CandidateBoundaryFixture mass_fixture;
  CandidateBoundaryFixtureSpec mass_spec;
  mass_spec.inlet = CandidateBoundaryInlet::mass_flow;
  mass_spec.mass_flow_rate = 0.375;
  CandidateBoundaryScratch mass_candidate;
  passed &= expect(
      mass_fixture.initialize(world, mass_spec) &&
          mass_fixture.stage(0.0, 8.0, 0.0, 13000U, mass_candidate),
      rank, "distributed real-chain mass-flow inlet finalizer succeeds");
  double local_mass = 0.0;
  if (mass_fixture.local_face_owner(CartesianFace::x_min)) {
    for (std::int32_t z = 0; z < mass_fixture.patch.cells.z; ++z)
      for (std::int32_t y = 0; y < mass_fixture.patch.cells.y; ++y)
        local_mass += mass_candidate.final_flux.x.unchecked({0, y, z});
  }
  double global_mass = 0.0;
  MPI_Allreduce(&local_mass, &global_mass, 1, MPI_DOUBLE, MPI_SUM, world);
  passed &= expect(
      std::abs(global_mass - mass_spec.mass_flow_rate) <=
          64.0 * std::numeric_limits<double>::epsilon(),
      rank, "distributed mass-flow capacity reaches one global target");
  if (!all_true(passed, world)) return false;

  CandidateBoundaryFixture backflow_fixture;
  CandidateBoundaryFixtureSpec backflow_spec;
  backflow_spec.allow_backflow = true;
  backflow_spec.multispecies = true;
  backflow_spec.backflow_velocity = -0.4;
  backflow_spec.backflow_temperature = 360.0;
  CandidateBoundaryScratch inward;
  CandidateBoundaryScratch recovered;
  passed &= expect(
      backflow_fixture.initialize(world, backflow_spec) &&
          backflow_fixture.stage(1.0, -1000.0, 0.0, 14000U, inward),
      rank, "distributed pressure outlet accepts configured backflow");
  if (backflow_fixture.local_face_owner(CartesianFace::x_max))
    passed &= expect(
        inward.mechanical_flux.x.unchecked(
            {backflow_fixture.patch.cells.x, 0, 0}) < 0.0 &&
            inward.final_flux.x.unchecked(
                {backflow_fixture.patch.cells.x, 0, 0}) < 0.0,
        rank, "outlet owner classifies inward provisional flux");
  passed &= expect(
      backflow_fixture.stage(1.0, 1000.0, 0.0, 15000U, recovered), rank,
      "distributed pressure outlet recovers to outflow");
  if (backflow_fixture.local_face_owner(CartesianFace::x_max))
    passed &= expect(
        recovered.mechanical_flux.x.unchecked(
            {backflow_fixture.patch.cells.x, 0, 0}) > 0.0 &&
            recovered.final_flux.x.unchecked(
                {backflow_fixture.patch.cells.x, 0, 0}) ==
                recovered.mechanical_flux.x.unchecked(
                    {backflow_fixture.patch.cells.x, 0, 0}),
        rank, "outlet recovery bit-preserves provisional outflow");
  if (!all_true(passed, world)) return false;

  CandidateBoundaryFixture disabled_fixture;
  CandidateBoundaryFixtureSpec disabled_spec;
  disabled_spec.allow_backflow = false;
  CandidateBoundaryScratch rejected_inward;
  CandidateBoundaryScratch disabled_replay;
  passed &= expect(disabled_fixture.initialize(world, disabled_spec), rank,
                   "disabled-backflow real compiler fixture initializes");
  const bool rejected_stage = disabled_fixture.stage(
      1.0, -1000.0, 0.0, 16000U, rejected_inward);
  const std::vector<double> rejected_values =
      face_values(as_const(rejected_inward.final_flux));
  passed &= expect(
      !rejected_stage &&
          disabled_fixture.diagnostic_status.code ==
              StatusCode::rejected_step &&
          !rejected_inward.final_boundary.valid() &&
          std::all_of(rejected_values.begin(), rejected_values.end(),
                      [](double value) { return value == -991.0; }) &&
          disabled_fixture.stage(1.0, 1000.0, 0.0, 17000U,
                                 disabled_replay),
      rank,
      "first disabled backflow rejects on every rank, zero-writes, and preserves replay authority");
  if (!all_true(passed, world)) return false;

  CandidateBoundaryFixture ibm_fixture;
  CandidateBoundaryFixtureSpec ibm_spec;
  ibm_spec.immersed = true;
  ibm_spec.cells_per_axis = 16;
  CandidateBoundaryScratch ibm_baseline;
  passed &= expect(
      ibm_fixture.initialize(world, ibm_spec) &&
          ibm_fixture.stage(0.0, 16.0, 0.0, 17500U, ibm_baseline) &&
          ibm_baseline.final_boundary.valid() &&
          static_cast<bool>(ibm_fixture.validate_zero_interface_flux(
              as_const(ibm_baseline.final_flux))),
      rank,
      "distributed real-chain IBM finalizer seals exact zero interface flux");
  std::uint64_t local_ibm_links = ibm_fixture.immersed_link_count();
  std::uint64_t global_ibm_links = 0U;
  MPI_Allreduce(&local_ibm_links, &global_ibm_links, 1, MPI_UINT64_T,
                MPI_SUM, world);
  passed &= expect(global_ibm_links > 0U, rank,
                   "distributed IBM fixture exercises control-face links");
  if (!all_true(passed, world)) return false;

  const std::vector<double> ibm_final_before =
      face_values(as_const(ibm_baseline.final_flux));
  const double ibm_density_before =
      ibm_baseline.density.view.unchecked(poison_cell, 0U);
  if (rank == size - 1)
    ibm_baseline.density.view.unchecked(poison_cell, 0U) =
        ibm_density_before * 1.01;
  FinalBoundaryFluxCertificate ibm_cleared = ibm_baseline.final_boundary;
  const Status ibm_poison_rejected = ibm_fixture.finalizer.finalize(
      ibm_fixture.finalizer_input(ibm_baseline), ibm_fixture.reductions,
      ibm_cleared);
  passed &= expect(
      ibm_poison_rejected.code == StatusCode::invalid_plan &&
          !ibm_cleared.valid() &&
          face_values(as_const(ibm_baseline.final_flux)) == ibm_final_before,
      rank,
      "one-rank IBM candidate poison rejects collectively with zero final-flux write");
  ibm_baseline.density.view.unchecked(poison_cell, 0U) = ibm_density_before;
  passed &= expect(
      static_cast<bool>(ibm_fixture.finalizer.finalize(
          ibm_fixture.finalizer_input(ibm_baseline), ibm_fixture.reductions,
          ibm_baseline.final_boundary)) &&
          ibm_baseline.final_boundary.valid() &&
          static_cast<bool>(ibm_fixture.validate_zero_interface_flux(
              as_const(ibm_baseline.final_flux))),
      rank, "failed IBM finalizer probe preserves replay authority");
  if (!all_true(passed, world)) return false;

  PisoFrozenMomentumExactCandidateCertificate exact_baseline;
  Status status = fixture.coupler.certify_frozen_momentum_exact_baseline(
      fixture.authority, baseline.pressure_stage, baseline.velocity_stage,
      baseline.flux_stage, fixture.exact_input(baseline),
      fixture.reductions, exact_baseline);
  passed &= expect(static_cast<bool>(status) && exact_baseline.valid(), rank,
                   "distributed open exact baseline certifies");
  CandidateBoundaryScratch selected;
  PisoFrozenMomentumExactCandidateCertificate exact_selected;
  if (status)
    status = fixture.stage(0.5, 16.0, 4.0, 18000U, selected)
                 ? Status{}
                 : fixture.diagnostic_status;
  if (status)
    status = fixture.coupler.certify_frozen_momentum_exact_candidate(
        fixture.authority, exact_baseline, selected.pressure_stage,
        selected.velocity_stage, selected.flux_stage,
        fixture.exact_input(selected), fixture.reductions, exact_selected);
  passed &= expect(static_cast<bool>(status) && exact_selected.valid(), rank,
                   "distributed open positive exact candidate certifies");
  if (!all_true(passed, world)) return false;

  PressureEnergyGlobalizationSample baseline_sample;
  baseline_sample.alpha = 0.0;
  baseline_sample.global_normalized_continuity = 1.0;
  baseline_sample.global_normalized_energy = 1.0;
  baseline_sample.thermodynamically_admissible = true;
  baseline_sample.state_and_flux_finite = true;
  baseline_sample.corrector = exact_baseline.corrector();
  baseline_sample.target_time = exact_baseline.target_time();
  baseline_sample.correction_direction = exact_baseline.correction_direction();
  baseline_sample.state_provenance =
      exact_baseline.candidate_state_provenance();
  baseline_sample.mass_flux_provenance =
      exact_baseline.candidate_mass_flux_provenance();
  std::array<PressureEnergyGlobalizationSample,
             kPressureEnergyGlobalizationCandidateCount>
      samples{};
  for (std::size_t index = 0U; index < samples.size(); ++index) {
    samples[index] = baseline_sample;
    samples[index].alpha = std::ldexp(1.0, -static_cast<int>(index));
    samples[index].global_normalized_continuity = 2.0;
    samples[index].global_normalized_energy = 2.0;
    samples[index].state_provenance = 19000U + index;
    samples[index].mass_flux_provenance = 20000U + index;
  }
  samples[1U].global_normalized_continuity = 0.2;
  samples[1U].global_normalized_energy = 0.2;
  samples[1U].state_provenance = exact_selected.candidate_state_provenance();
  samples[1U].mass_flux_provenance =
      exact_selected.candidate_mass_flux_provenance();
  PressureEnergyGlobalizationSelectionCertificate selection;
  passed &= expect(
      static_cast<bool>(select_pressure_energy_globalization(
          baseline_sample, {samples.data(), samples.size()}, selection)) &&
          selection.valid(),
      rank, "distributed open selector binds physical final flux");
  if (!all_true(passed, world)) return false;

  FaceFluxStorage commit_storage;
  FaceFluxView committed_flux;
  passed &= expect(
      static_cast<bool>(FaceFluxStorage::allocate_workspace(
          cells, 1U, commit_storage)) &&
          static_cast<bool>(commit_storage.workspace_view(
              0U, 21000U, committed_flux)),
      rank, "distributed open commit flux allocates independently");
  const std::vector<double> pressure_before = fixture.pressure.storage;
  const std::vector<double> enthalpy_before = fixture.enthalpy.storage;
  const std::vector<double> density_before = fixture.density.storage;
  const std::vector<double> temperature_before = fixture.temperature.storage;
  const std::vector<double> live_velocity_before = fixture.velocity.storage;
  const std::vector<double> committed_before =
      face_values(as_const(committed_flux));
  const double selected_velocity_before =
      selected.velocity.view.unchecked(poison_cell, 0U);
  if (rank == size - 1)
    selected.velocity.view.unchecked(poison_cell, 0U) =
        selected_velocity_before + 0.25;
  PisoStateCorrectionCertificate rejected_commit;
  const Status commit_rejected =
      fixture.coupler.commit_frozen_momentum_coupled_trial_state(
          fixture.authority, exact_selected, selection,
          {fixture.velocity.view, fixture.pressure.view,
           fixture.enthalpy.view, fixture.density.view,
           fixture.temperature.view},
          committed_flux, fixture.reductions, rejected_commit);
  passed &= expect(
      commit_rejected.code == StatusCode::invalid_plan &&
          !rejected_commit.valid() &&
          fixture.pressure.storage == pressure_before &&
          fixture.enthalpy.storage == enthalpy_before &&
          fixture.density.storage == density_before &&
          fixture.temperature.storage == temperature_before &&
          fixture.velocity.storage == live_velocity_before &&
          face_values(as_const(committed_flux)) == committed_before,
      rank,
      "one-rank selected-state poison rejects publication collectively with zero live write");
  selected.velocity.view.unchecked(poison_cell, 0U) =
      selected_velocity_before;
  PisoStateCorrectionCertificate committed;
  const Status commit_status =
      fixture.coupler.commit_frozen_momentum_coupled_trial_state(
          fixture.authority, exact_selected, selection,
          {fixture.velocity.view, fixture.pressure.view,
           fixture.enthalpy.view, fixture.density.view,
           fixture.temperature.view},
          committed_flux, fixture.reductions, committed);
  passed &= expect(
      static_cast<bool>(commit_status) && committed.valid() &&
          face_values(as_const(committed_flux)) ==
              face_values(as_const(selected.final_flux)),
      rank,
      "distributed open selection atomically publishes candidate state/final flux");
  return all_true(passed, world);
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    return 2;
  }
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  const bool normalized_history =
      test_c1_normalized_committed_face_history_across_decomposition_seams(
          MPI_COMM_WORLD, rank);
  const bool pressure_fraction =
      test_pressure_correction_fraction_red(MPI_COMM_WORLD, rank);
  const bool exact_eos_transaction =
      test_exact_eos_correction_collective_transaction(MPI_COMM_WORLD, rank);
  const bool c2_applied_candidate =
      test_c2_applied_candidate_audit_red(MPI_COMM_WORLD, rank);
  const bool negative_density =
      test_negative_density_refresh_collective(MPI_COMM_WORLD, rank);
  const bool post_halo =
      test_post_halo_thermophysical_revalidation_collective(MPI_COMM_WORLD,
                                                            rank);
  const bool frozen_candidate =
      test_frozen_momentum_candidate_collective(MPI_COMM_WORLD, rank);
  const bool candidate_boundary =
      test_candidate_boundary_finalizer_collective(MPI_COMM_WORLD, rank);
  const bool full_piso = test_full_distributed_piso(MPI_COMM_WORLD, rank);
  const bool passed =
      normalized_history && pressure_fraction && exact_eos_transaction &&
      negative_density && post_halo && frozen_candidate &&
      candidate_boundary && full_piso && c2_applied_candidate;
  if (!passed && rank == 0) {
    std::cerr << "solver-piso-mpi summary normalized/history="
              << normalized_history << " pressure-fraction="
              << pressure_fraction << " exact-eos=" << exact_eos_transaction
              << " negative-density=" << negative_density
              << " post-halo=" << post_halo
              << " frozen-candidate=" << frozen_candidate
              << " candidate-boundary=" << candidate_boundary
              << " full-piso=" << full_piso
              << " c2-applied=" << c2_applied_candidate << '\n';
  }
  MPI_Finalize();
  return passed ? 0 : 1;
}
