// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_product.hpp"
#include "hundun/v04_app.hpp"

#include "core_product_freeze_detail.hpp"
#include "../support/product_fixture.hpp"

#include <mpi.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <limits>
#include <string_view>
#include <vector>

namespace {

using namespace hundun::v04;

bool expect(bool condition, const char* description) {
  if (!condition) std::cerr << "FAIL: " << description << '\n';
  return condition;
}

bool test_incomplete_registration() {
  bool passed = true;
  const std::uint32_t complete =
      detail::kProductRequiredCapabilities | detail::product_ibm_donors;
  passed &= expect(static_cast<bool>(
                       detail::validate_product_capabilities_for_test(
                           complete, true)),
                   "complete capability manifest passes");
  for (std::uint32_t bit = 1U;
       bit <= static_cast<std::uint32_t>(detail::product_collective_epochs);
       bit <<= 1U) {
    passed &= expect(
        !detail::validate_product_capabilities_for_test(complete & ~bit, true),
        "each missing capability rejects before allocation");
  }
  passed &= expect(static_cast<bool>(
                       detail::validate_product_capabilities_for_test(
                           complete & ~detail::product_ibm_donors, false)),
                   "non-IBM manifest does not require donor registration");
  return passed;
}

bool test_pressure_energy_temporal_operand_scale() {
  using State = detail::ProductPressureEnergyTemporalState;
  bool passed = true;
  double scale = 0.0;
  passed &= expect(
      detail::product_pressure_energy_temporal_operand_scale(
          {15.0, -20.0, 5.0, 2U}, 2.0,
          State{2.0, 3.0, 4.0}, State{1.0, 5.0, 7.0},
          State{4.0, 2.0, 3.0}, scale) &&
          scale == 890.0,
      "BDF2 energy scale retains all rho-h and pressure history operands");

  const double nan = std::numeric_limits<double>::quiet_NaN();
  passed &= expect(
      detail::product_pressure_energy_temporal_operand_scale(
          {10.0, -10.0, 0.0, 1U}, 2.0,
          State{2.0, 3.0, 4.0}, State{1.0, 5.0, 7.0},
          State{nan, nan, nan}, scale) &&
          scale == 440.0,
      "BE recovery ignores the unbound previous energy layer");

  constexpr double inverse_dt = 1.0e12;
  const State target{1.0, 1.0, 1.0};
  const State large_history{1.0, 1.0e9, 1.0e9};
  passed &= expect(
      detail::product_pressure_energy_temporal_operand_scale(
          {1.5 * inverse_dt, -2.0 * inverse_dt, 0.5 * inverse_dt, 2U},
          1.0, target, large_history, large_history, scale) &&
          scale > 1.0e6 * 1.5 * inverse_dt * 1.0,
      "small-dt cancellation is scaled by large BDF history operands");
  passed &= expect(
      !detail::product_pressure_energy_temporal_operand_scale(
          {1.5 * inverse_dt, -2.0 * inverse_dt, 0.5 * inverse_dt, 2U},
          1.0, target, large_history, State{nan, nan, nan}, scale),
      "BDF2 rejects an invalid required history operand");
  return passed;
}

bool test_pressure_inexact_forcing_policy() {
  const LinearSolveControl base{1.0e-8, 1.0e-6, 500U, 16U, 30U};
  const auto same_fixed_control = [&](const LinearSolveControl& value) {
    return value.absolute_tolerance == base.absolute_tolerance &&
           value.maximum_iterations == base.maximum_iterations &&
           value.true_residual_interval == base.true_residual_interval &&
           value.restart == base.restart;
  };
  bool passed = true;

  const LinearSolveControl cold =
      detail::product_pressure_inexact_forcing_control(
          base, {0.0, 0.0, 0.0, 1.0e-6, false, false});
  passed &= expect(same_fixed_control(cold) &&
                       cold.relative_tolerance == 1.0e-4,
                   "cold inexact pressure solve uses the guarded ceiling");

  const LinearSolveControl intermediate =
      detail::product_pressure_inexact_forcing_control(
          base, {5.0e-5, 4.0e-5, 0.0, 1.0e-6, true, false});
  passed &= expect(same_fixed_control(intermediate) &&
                       intermediate.relative_tolerance == 5.0e-6,
                   "inexact pressure tolerance follows nonlinear residual");

  const LinearSolveControl terminal =
      detail::product_pressure_inexact_forcing_control(
          base, {9.0e-6, 8.0e-6, 0.0, 1.0e-6, true, false});
  passed &= expect(same_fixed_control(terminal) &&
                       terminal.relative_tolerance ==
                           base.relative_tolerance,
                   "terminal band restores the case pressure tolerance");

  const LinearSolveControl stagnated =
      detail::product_pressure_inexact_forcing_control(
          base, {1.0e-3, 8.0e-4, 1.05e-3, 1.0e-6, true, true});
  passed &= expect(same_fixed_control(stagnated) &&
                       stagnated.relative_tolerance ==
                           base.relative_tolerance,
                   "stagnating nonlinear residual restores full accuracy");

  const LinearSolveControl contracting =
      detail::product_pressure_inexact_forcing_control(
          base, {8.0e-4, 7.0e-4, 1.0e-3, 1.0e-6, true, true});
  passed &= expect(same_fixed_control(contracting) &&
                       contracting.relative_tolerance == 8.0e-5,
                   "contracting far-field residual retains inexact solve");

  const double nan = std::numeric_limits<double>::quiet_NaN();
  const LinearSolveControl invalid =
      detail::product_pressure_inexact_forcing_control(
          base, {nan, 1.0e-3, 0.0, 1.0e-6, true, false});
  passed &= expect(same_fixed_control(invalid) &&
                       invalid.relative_tolerance ==
                           base.relative_tolerance,
                   "non-finite forcing input fails closed");

  const LinearSolveControl already_loose{1.0e-8, 1.0e-4, 500U, 16U, 30U};
  const LinearSolveControl unchanged =
      detail::product_pressure_inexact_forcing_control(
          already_loose, {1.0, 1.0, 0.0, 1.0e-4, false, false});
  passed &= expect(unchanged.absolute_tolerance ==
                           already_loose.absolute_tolerance &&
                       unchanged.relative_tolerance ==
                           already_loose.relative_tolerance &&
                       unchanged.maximum_iterations ==
                           already_loose.maximum_iterations &&
                       unchanged.true_residual_interval ==
                           already_loose.true_residual_interval &&
                       unchanged.restart == already_loose.restart,
                   "an already-loose case authority is not changed");

  const LinearSolveControl deliberately_strict{1.0e-15, 1.0e-13, 800U, 4U,
                                                64U};
  const LinearSolveControl strict_unchanged =
      detail::product_pressure_inexact_forcing_control(
          deliberately_strict,
          {0.0, 0.0, 0.0, 1.0e-12, false, false});
  passed &= expect(
      strict_unchanged.relative_tolerance ==
              deliberately_strict.relative_tolerance &&
          strict_unchanged.absolute_tolerance ==
              deliberately_strict.absolute_tolerance &&
          strict_unchanged.maximum_iterations ==
              deliberately_strict.maximum_iterations &&
          strict_unchanged.true_residual_interval ==
              deliberately_strict.true_residual_interval &&
          strict_unchanged.restart == deliberately_strict.restart,
      "a case tolerance stricter than its terminal gate disables forcing");
  return passed;
}

bool test_pressure_aitken_initial_alpha() {
  bool passed = true;
  passed &= expect(
      std::abs(detail::product_pressure_aitken_initial_alpha(1.0, 0.2, 1.0) -
               1.25) <= 8.0 * std::numeric_limits<double>::epsilon(),
      "contracting nonlinear merit produces the scalar Aitken factor");
  passed &= expect(
      std::abs(detail::product_pressure_aitken_initial_alpha(1.0, 0.2, 1.4) -
                   1.75) <=
          8.0 * std::numeric_limits<double>::epsilon(),
      "the previous relaxation removes the observed-contraction lag");
  passed &= expect(
      detail::product_pressure_aitken_initial_alpha(1.0, 0.5, 1.0) == 2.0 &&
          detail::product_pressure_aitken_initial_alpha(1.0, 0.8, 1.5) == 2.0,
      "Aitken extrapolation is capped at two");
  passed &= expect(
      detail::product_pressure_aitken_initial_alpha(1.0, 0.95, 1.5) == 1.0 &&
          detail::product_pressure_aitken_initial_alpha(1.0, 1.0, 1.5) == 1.0,
      "stagnating merit disables extrapolation");
  const double nan = std::numeric_limits<double>::quiet_NaN();
  passed &= expect(
      detail::product_pressure_aitken_initial_alpha(0.0, 0.5, 1.0) == 1.0 &&
          detail::product_pressure_aitken_initial_alpha(1.0, nan, 1.0) == 1.0 &&
          detail::product_pressure_aitken_initial_alpha(1.0, 0.2, 0.5) == 1.0,
      "invalid Aitken evidence fails closed to the legacy full step");
  return passed;
}

bool test_freeze() {
  ValidatedModel model = test::product_model();
  model.solver.pressure = {1.0e-8, 2.0e-7, 333U, 7U, 16U};
  model.solver.terminal = {3.0e-6, 4.0e-6, 5.0e-6, 6.0e-6};
  CompiledCasePlan plan;
  Status status =
      ProductCompiler::compile(MPI_COMM_SELF, model, {}, plan);
  if (!status) {
    std::cerr << "status=" << static_cast<unsigned>(status.code)
              << " detail=" << status.detail << '\n';
  }
  bool passed = true;
  passed &= expect(static_cast<bool>(status), "complete product freezes");
  if (!status) return false;
  const PlanSummary summary = plan.summary();
  const Span<const ProductFreezePhase> order = plan.freeze_order();
  passed &= expect(plan.fingerprint() != 0U, "sealed fingerprint is valid");
  passed &= expect(order.size == test::kFreezeOrder.size(),
                   "freeze order records every phase");
  for (std::size_t index = 0U;
       index < order.size && index < test::kFreezeOrder.size(); ++index) {
    passed &= expect(order.data[index] == test::kFreezeOrder[index],
                     "freeze phase order is exact");
  }
  passed &= expect(summary.global_cells.x == 17 &&
                       summary.global_cells.y == 11 &&
                       summary.global_cells.z == 7,
                   "geometry precedes production sizing");
  passed &= expect(summary.field_count != 0U &&
                       summary.arena_doubles != 0U &&
                       summary.graph_stage_count == 15U &&
                       summary.graph_node_count != 0U,
                   "schema, arena and graph are complete");
  passed &= expect(summary.pressure_correctors == 2U,
                   "exactly two pressure correctors are sealed");
  passed &= expect(summary.pressure_absolute_tolerance == 1.0e-8 &&
                       summary.pressure_relative_tolerance == 2.0e-7 &&
                       summary.pressure_maximum_iterations == 333U &&
                       summary.pressure_true_residual_interval == 7U &&
                       summary.pressure_krylov_restart == 16U &&
                       summary.terminal_eos_tolerance == 3.0e-6 &&
                       summary.terminal_continuity_tolerance == 4.0e-6 &&
                       summary.terminal_closed_mass_tolerance == 5.0e-6 &&
                       summary.terminal_gauge_tolerance == 6.0e-6,
                   "case solver controls are frozen into the product plan");
  passed &= expect(!summary.exact_numeric_certified &&
                       !summary.preconditioner_setup_certified,
                   "uninitialized numeric capacities remain uncertified");
  passed &= expect(summary.sealed && !summary.immersed,
                   "non-IBM product is sealed");
  passed &= expect(plan.field_schema() != nullptr &&
                       plan.arena_layout() != nullptr &&
                       plan.execution_graph() != nullptr &&
                       plan.io_services() != nullptr &&
                       plan.contribution_plan() != nullptr &&
                       plan.contribution_plan()->valid(),
                   "all immutable bundles are bound");
  passed &= expect(plan.cpu_execution_plan() != nullptr &&
                       plan.cpu_execution_plan()->pure_mpi() &&
                       plan.cpu_execution_plan()->worker_count() == 1U &&
                       plan.cpu_plan_fingerprint() != 0U &&
                       plan.stl_fingerprint() == 0U,
                   "pure-MPI CPU execution policy is frozen into product identity");
  passed &= expect(plan.state_storage_address() != 0U &&
                       plan.krylov_storage_address() != 0U &&
                       plan.mg_storage_address() != 0U,
                   "persistent storage is allocated and bound");
  passed &= expect(plan.execution_graph()->stage(40U) != nullptr &&
                       plan.execution_graph()->stage(50U) != nullptr &&
                       plan.execution_graph()->stage(51U) == nullptr,
                   "graph owns exactly the two PISO stages");
  FieldId pressure_correction{};
  bool found_pressure_correction = false;
  for (const FieldDescriptor& field : *plan.field_schema()) {
    if (field.stable_name == "delta_pi") {
      pressure_correction = field.id;
      found_pressure_correction = true;
      break;
    }
  }
  const FieldAccessSpec pressure_correction_workspace{
      pressure_correction, StateVisibility::workspace};
  const auto contains_pressure_correction =
      [&](Span<const FieldAccessSpec> accesses) {
        for (std::size_t index = 0U; index < accesses.size; ++index) {
          if (accesses.data[index] == pressure_correction_workspace)
            return true;
        }
        return false;
      };
  const FrozenExecutionGraph& graph = *plan.execution_graph();
  passed &= expect(
      found_pressure_correction &&
          contains_pressure_correction(graph.reads(40U)) &&
          contains_pressure_correction(graph.writes(40U)) &&
          contains_pressure_correction(graph.invalidations(40U)) &&
          contains_pressure_correction(graph.reads(50U)) &&
          contains_pressure_correction(graph.writes(50U)) &&
          contains_pressure_correction(graph.invalidations(50U)),
      "both PISO stages declare pressure correction read-invalidate-write");
  passed &= expect(plan.execution_graph()->stage(12U) != nullptr &&
                       plan.execution_graph()->stage(45U) != nullptr,
                   "graph freezes C1 and C2 model-integration seams");
  struct PressureEnergyWorkspaceSpec {
    std::string_view stable_name;
    std::uint8_t ghost_width;
  };
  constexpr std::array<PressureEnergyWorkspaceSpec, 12U>
      pressure_energy_workspaces{{
          // Physical thermodynamic ghosts are reconstructed from the same
          // p/h/Y state; rho_h therefore carries the state ghost width.
          {"drho_dh_pY", 2U},
          {"pressure_energy_C_h", 0U},
          {"pressure_energy_C_h_row_scale", 0U},
          {"pressure_energy_E_p", 0U},
          {"pressure_energy_E_h", 0U},
          {"pressure_energy_R_C", 0U},
          {"pressure_energy_R_E", 0U},
          {"delta_h", 0U},
          // E_h differentiates temperature-space conduction as
          // lambda*grad(delta_h/cp); it cannot alias the Schur delta_h view.
          {"pressure_energy_delta_T", 1U},
          {"pressure_energy_schur_continuity_response", 0U},
          {"pressure_energy_schur_eliminated_h", 2U},
          {"pressure_energy_schur_energy_response", 0U},
      }};
  std::array<FieldId, pressure_energy_workspaces.size()> workspace_ids{};
  for (std::size_t expected = 0U;
       expected < pressure_energy_workspaces.size(); ++expected) {
    const PressureEnergyWorkspaceSpec specification =
        pressure_energy_workspaces[expected];
    const FieldDescriptor* descriptor = nullptr;
    for (const FieldDescriptor& field : *plan.field_schema()) {
      if (field.stable_name == specification.stable_name) {
        descriptor = &field;
        break;
      }
    }
    passed &= expect(descriptor != nullptr,
                     "pressure-energy workspace is cold-registered");
    if (descriptor == nullptr) continue;
    workspace_ids[expected] = descriptor->id;
    const ArenaFieldLayout* layout = plan.arena_layout()->field(descriptor->id);
    passed &= expect(descriptor->components == 1U &&
                         descriptor->ghost_width == specification.ghost_width &&
                         layout != nullptr &&
                         layout->lifetime == FieldLifetime::persistent_workspace &&
                         layout->replicas == 1U,
                     "pressure-energy storage is an independent persistent workspace");
  }
  const auto contains_workspace = [](Span<const FieldAccessSpec> accesses,
                                     FieldId field) {
    const FieldAccessSpec expected{field, StateVisibility::workspace};
    for (std::size_t index = 0U; index < accesses.size; ++index)
      if (accesses.data[index] == expected) return true;
    return false;
  };
  const auto contains_trial = [](Span<const FieldAccessSpec> accesses,
                                 FieldId field) {
    const FieldAccessSpec expected{field, StateVisibility::trial};
    for (std::size_t index = 0U; index < accesses.size; ++index)
      if (accesses.data[index] == expected) return true;
    return false;
  };
  const FrozenStage* c1 = graph.stage(40U);
  const FrozenStage* c2 = graph.stage(50U);
  passed &= expect(c1 != nullptr && c2 != nullptr &&
                       c1->resources.numeric_refills == 3U &&
                       c1->resources.linear_iterations == 400U &&
                       c2->resources.numeric_refills == 3U &&
                       c2->resources.linear_iterations == 400U,
                   "C1/C2 reserve the coupled block numeric lifecycle");
  for (FieldId field : workspace_ids) {
    passed &= expect(field != 0U &&
                         contains_workspace(graph.writes(40U), field) &&
                         contains_workspace(graph.reads(50U), field) &&
                         contains_workspace(graph.writes(50U), field),
                     "C1 publishes and C2 refreshes every pressure-energy workspace");
  }
  struct PressureEnergyCandidateWorkspaceSpec {
    std::string_view stable_name;
    std::uint8_t components;
    std::uint8_t ghost_width;
  };
  constexpr std::array<PressureEnergyCandidateWorkspaceSpec, 5U>
      candidate_state_workspaces{{
          {"pressure_energy_candidate_pi", 1U, 2U},
          {"pressure_energy_candidate_h", 1U, 2U},
          {"pressure_energy_candidate_rho", 1U, 2U},
          {"pressure_energy_candidate_T", 1U, 2U},
          {"pressure_energy_candidate_U", 3U, 2U},
      }};
  for (const PressureEnergyCandidateWorkspaceSpec& specification :
       candidate_state_workspaces) {
    const FieldDescriptor* descriptor = nullptr;
    for (const FieldDescriptor& field : *plan.field_schema()) {
      if (field.stable_name == specification.stable_name) {
        descriptor = &field;
        break;
      }
    }
    passed &= expect(descriptor != nullptr,
                     "complete pressure-energy candidate state is registered");
    if (descriptor == nullptr) continue;
    const ArenaFieldLayout* layout = plan.arena_layout()->field(descriptor->id);
    passed &= expect(
        descriptor->components == specification.components &&
            descriptor->ghost_width == specification.ghost_width &&
            layout != nullptr &&
            layout->lifetime == FieldLifetime::persistent_workspace &&
            layout->replicas == 1U &&
            contains_workspace(graph.writes(40U), descriptor->id) &&
            contains_workspace(graph.reads(50U), descriptor->id) &&
            contains_workspace(graph.writes(50U), descriptor->id),
        "candidate p/h/rho/T/U has frozen residual reach and C1/C2 lifecycle");
  }
  constexpr std::array<PressureEnergyCandidateWorkspaceSpec, 9U>
      candidate_scratch_workspaces{{
          {"pressure_energy_candidate_delta_pi", 1U, 2U},
          {"pressure_energy_candidate_mu", 1U, 2U},
          {"pressure_energy_candidate_mu_eff", 1U, 2U},
          {"pressure_energy_candidate_grad_U", 9U, 2U},
          {"pressure_energy_candidate_drho_dp_hY", 1U, 2U},
          {"pressure_energy_candidate_drho_dh_pY", 1U, 2U},
          {"pressure_energy_candidate_lambda", 1U, 2U},
          {"pressure_energy_candidate_cp", 1U, 2U},
          {"pressure_energy_candidate_lambda_over_cp", 1U, 2U},
      }};
  for (const PressureEnergyCandidateWorkspaceSpec& specification :
       candidate_scratch_workspaces) {
    const FieldDescriptor* descriptor = nullptr;
    for (const FieldDescriptor& field : *plan.field_schema())
      if (field.stable_name == specification.stable_name) {
        descriptor = &field;
        break;
      }
    passed &= expect(descriptor != nullptr,
                     "candidate correction/material scratch is registered");
    if (descriptor == nullptr) continue;
    const ArenaFieldLayout* layout = plan.arena_layout()->field(descriptor->id);
    passed &= expect(
        descriptor->components == specification.components &&
            descriptor->ghost_width == specification.ghost_width &&
            layout != nullptr &&
            layout->lifetime == FieldLifetime::persistent_workspace &&
            layout->replicas == 1U &&
            contains_workspace(graph.writes(40U), descriptor->id) &&
            contains_workspace(graph.reads(50U), descriptor->id) &&
            contains_workspace(graph.writes(50U), descriptor->id),
        "candidate correction/material scratch has a complete C1/C2 lifecycle");
  }
  for (std::size_t index : {3U, 4U}) {
    const FieldId field = workspace_ids[index];
    passed &= expect(
        field != 0U && contains_workspace(graph.reads(40U), field) &&
            contains_workspace(graph.writes(40U), field) &&
            contains_workspace(graph.invalidations(40U), field) &&
            contains_workspace(graph.reads(50U), field) &&
            contains_workspace(graph.writes(50U), field) &&
            contains_workspace(graph.invalidations(50U), field),
        "C1/C2 explicitly read-invalidate-write both energy Jacobian blocks");
  }
  constexpr std::array<std::string_view, 5U> coupled_state_names{{
      "rho", "U", "pi", "h", "T"}};
  for (std::string_view stable_name : coupled_state_names) {
    FieldId field = 0U;
    bool found = false;
    for (const FieldDescriptor& descriptor : *plan.field_schema()) {
      if (descriptor.stable_name == stable_name) {
        field = descriptor.id;
        found = true;
        break;
      }
    }
    passed &= expect(found && contains_trial(graph.reads(40U), field) &&
                         contains_trial(graph.writes(40U), field) &&
                         contains_trial(graph.invalidations(40U), field) &&
                         contains_trial(graph.reads(50U), field) &&
                         contains_trial(graph.writes(50U), field) &&
                         contains_trial(graph.invalidations(50U), field),
                     "C1/C2 declare the complete same-target primitive state");
  }
  passed &= expect(plan.io_services()->snapshot_fields().size == 3U &&
                       plan.io_services()->services().size == 5U,
                   "all cold services share the committed snapshot schema");

  const PlanFingerprint fingerprint = plan.fingerprint();
  const std::uintptr_t state = plan.state_storage_address();
  CompiledCasePlan moved = std::move(plan);
  passed &= expect(plan.fingerprint() == 0U &&
                       moved.fingerprint() == fingerprint &&
                       moved.state_storage_address() == state,
                   "move preserves the sealed product and invalidates source");
  status = ProductCompiler::compile(MPI_COMM_SELF, model, {}, moved);
  passed &= expect(status.code == StatusCode::invalid_plan &&
                       moved.fingerprint() == fingerprint,
                   "a sealed product rejects a second freeze atomically");
  ValidatedModel invalid = model;
  invalid.fingerprint = 0U;
  CompiledCasePlan empty;
  status = ProductCompiler::compile(MPI_COMM_SELF, invalid, {}, empty);
  passed &= expect(status.code == StatusCode::invalid_plan &&
                       empty.fingerprint() == 0U,
                   "invalid input rejects without publishing a product");
  return passed;
}

bool test_live_thermal_halo_resource_contract() {
  constexpr Int3 cells{17, 11, 7};
  CompiledCasePlan plan;
  const Status status =
      ProductCompiler::compile(MPI_COMM_SELF, test::product_model(cells), {},
                               plan);
  const FrozenExecutionGraph* graph = plan.execution_graph();
  const FrozenStage* stage20 = graph == nullptr ? nullptr : graph->stage(20U);
  const FrozenStage* stage40 = graph == nullptr ? nullptr : graph->stage(40U);
  const FrozenStage* stage45 = graph == nullptr ? nullptr : graph->stage(45U);
  const FrozenStage* stage50 = graph == nullptr ? nullptr : graph->stage(50U);
  const FrozenStage* stage60 = graph == nullptr ? nullptr : graph->stage(60U);

  std::size_t thermal_bytes = 0U;
  std::size_t turbulence_bytes = 0U;
  std::size_t pressure_bytes = 0U;
  std::size_t candidate_correction_bytes = 0U;
  std::size_t candidate_state_full_bytes = 0U;
  std::size_t candidate_state_face_bytes = 0U;
  std::size_t candidate_finalizer_bytes = 0U;
  std::size_t force_bytes = 0U;
  const bool sized =
      detail::product_halo_bytes(cells, 2U, 1U, thermal_bytes) &&
      detail::product_halo_bytes(cells, 6U, 2U, turbulence_bytes) &&
      detail::product_halo_bytes(cells, 8U, 2U, pressure_bytes) &&
      detail::product_halo_bytes(cells, 1U, 1U,
                                 candidate_correction_bytes) &&
      detail::product_halo_bytes(cells, 6U, 2U,
                                 candidate_state_full_bytes) &&
      detail::product_halo_bytes(cells, 1U, 1U,
                                 candidate_state_face_bytes) &&
      detail::product_halo_bytes(cells, 7U, 1U,
                                 candidate_finalizer_bytes) &&
      detail::product_halo_bytes(cells, 15U, 2U, force_bytes);
  if (!expect(status && graph != nullptr && stage20 != nullptr &&
                  stage40 != nullptr && stage45 != nullptr &&
                  stage50 != nullptr && stage60 != nullptr && sized,
              "thermal resource-contract fixture freezes"))
    return false;

  constexpr std::uint64_t candidate_evaluations =
      static_cast<std::uint64_t>(
          kPressureEnergyGlobalizationCandidateCount + 2U);
  const std::uint64_t candidate_bytes_per_evaluation =
      candidate_correction_bytes + candidate_state_full_bytes +
      candidate_state_face_bytes + thermal_bytes +
      candidate_finalizer_bytes;
  const std::uint64_t candidate_messages_per_evaluation =
      candidate_bytes_per_evaluation / sizeof(double);
  const std::uint64_t thermal_messages =
      thermal_bytes / sizeof(double);
  const std::uint64_t stage40_bytes =
      pressure_bytes +
      candidate_evaluations * candidate_bytes_per_evaluation + thermal_bytes;
  const std::uint64_t stage40_messages =
      6U + candidate_evaluations * candidate_messages_per_evaluation +
      thermal_messages;
  const std::uint64_t refinement_exchanges =
      static_cast<std::uint64_t>(kPressureEnergyRefinementCapacity);

  return expect(
      stage20->resources.merged_halo_bytes ==
              turbulence_bytes + thermal_bytes &&
          stage20->resources.merged_halo_messages ==
              6U + thermal_messages &&
          stage40->resources.merged_halo_bytes == stage40_bytes &&
          stage40->resources.merged_halo_messages == stage40_messages &&
          stage45->resources.merged_halo_bytes == 2U * thermal_bytes &&
          stage45->resources.merged_halo_messages ==
              2U * thermal_messages &&
          stage50->resources.merged_halo_bytes ==
              stage40_bytes +
                  refinement_exchanges *
                      (candidate_evaluations *
                           candidate_bytes_per_evaluation +
                       thermal_bytes) &&
          stage50->resources.merged_halo_messages ==
              stage40_messages +
                  refinement_exchanges *
                      (candidate_evaluations *
                           candidate_messages_per_evaluation +
                       thermal_messages) &&
          stage60->resources.merged_halo_bytes ==
              force_bytes + thermal_bytes &&
          stage60->resources.merged_halo_messages ==
              6U + thermal_messages,
      "every live and refinement-candidate halo is present in its stage seal");
}

bool test_pressure_energy_restart_schema() {
  ValidatedModel model = test::product_model({17, 11, 7});
  CompiledCasePlan plan;
  Status status = ProductCompiler::compile(MPI_COMM_SELF, model, {}, plan);
  std::array<FieldId, 26U> pressure_energy_workspaces{};
  if (status) {
    constexpr std::array<std::string_view, 26U> names{{
        "drho_dh_pY",
        "pressure_energy_C_h",
        "pressure_energy_C_h_row_scale",
        "pressure_energy_E_p",
        "pressure_energy_E_h",
        "pressure_energy_R_C",
        "pressure_energy_R_E",
        "delta_h",
        "pressure_energy_delta_T",
        "pressure_energy_schur_continuity_response",
        "pressure_energy_schur_eliminated_h",
        "pressure_energy_schur_energy_response",
        "pressure_energy_candidate_pi",
        "pressure_energy_candidate_h",
        "pressure_energy_candidate_rho",
        "pressure_energy_candidate_T",
        "pressure_energy_candidate_U",
        "pressure_energy_candidate_delta_pi",
        "pressure_energy_candidate_mu",
        "pressure_energy_candidate_mu_eff",
        "pressure_energy_candidate_grad_U",
        "pressure_energy_candidate_drho_dp_hY",
        "pressure_energy_candidate_drho_dh_pY",
        "pressure_energy_candidate_lambda",
        "pressure_energy_candidate_cp",
        "pressure_energy_candidate_lambda_over_cp",
    }};
    for (std::size_t expected = 0U; expected < names.size(); ++expected)
      for (const FieldDescriptor& field : *plan.field_schema())
        if (field.stable_name == names[expected])
          pressure_energy_workspaces[expected] = field.id;
  }
  ProductDriver driver;
  if (status)
    status = ProductDriver::create(MPI_COMM_SELF, std::move(plan), driver);
  RestartExpected expected;
  if (status) status = driver.restart_expected(expected);
  bool restart_schema_is_state_only = status && expected.fields.size == 3U;
  for (std::size_t index = 0U;
       index < expected.fields.size && restart_schema_is_state_only; ++index)
    for (FieldId workspace : pressure_energy_workspaces)
      restart_schema_is_state_only &=
          workspace != 0U && expected.fields.data[index].field != workspace;
  return expect(status && restart_schema_is_state_only,
                "pressure-energy workspaces never enter Restart state");
}

bool test_pressure_energy_candidate_storage_lineage() {
  ValidatedModel model = test::product_model({17, 11, 7});
  CompiledCasePlan plan;
  Status status = ProductCompiler::compile(MPI_COMM_SELF, model, {}, plan);
  detail::PressureEnergyCandidateStorageDiagnostic diagnostic;
  bool passed = status &&
                detail::pressure_energy_candidate_storage_diagnostic_for_test(
                    diagnostic);
  passed &= expect(passed && diagnostic.valid &&
                       diagnostic.plan == plan.fingerprint() &&
                       diagnostic.lineage_fingerprint != 0U,
                   "sealed plan publishes candidate storage lineage");
  const auto disjoint = [](std::uintptr_t left_begin,
                           std::uintptr_t left_end,
                           std::uintptr_t right_begin,
                           std::uintptr_t right_end) {
    return left_begin != 0U && left_begin < left_end && right_begin != 0U &&
           right_begin < right_end &&
           (left_end <= right_begin || right_end <= left_begin);
  };
  constexpr std::array<std::uint8_t,
                       detail::kPressureEnergyCandidateFieldCount>
      components{{1U, 1U, 1U, 1U, 3U}};
  constexpr std::array<std::uint8_t,
                       detail::kPressureEnergyCandidateFieldCount>
      ghosts{{2U, 2U, 2U, 2U, 2U}};
  for (std::size_t index = 0U; index < diagnostic.fields.size(); ++index) {
    const auto& field = diagnostic.fields[index];
    passed &= expect(
        field.candidate_field != field.trial_field &&
            field.candidate_revision != 0U && field.trial_revision != 0U &&
            field.candidate_revision != field.trial_revision &&
            field.candidate_storage == field.trial_storage &&
            field.candidate_revision_domain == field.trial_revision_domain &&
            field.components == components[index] &&
            field.ghost_width == ghosts[index] &&
            disjoint(field.candidate_begin, field.candidate_end,
                     field.trial_begin, field.trial_end),
        "candidate field has a disjoint base and distinct revision lineage");
    for (std::size_t other = index + 1U;
         other < diagnostic.fields.size(); ++other) {
      const auto& right = diagnostic.fields[other];
      passed &= expect(
          disjoint(field.candidate_begin, field.candidate_end,
                   right.candidate_begin, right.candidate_end),
          "candidate primitive workspaces do not alias one another");
    }
  }
  constexpr std::array<std::uint8_t,
                       detail::kPressureEnergyCandidateScratchFieldCount>
      scratch_ghosts{{2U, 2U, 2U, 2U, 2U, 2U, 2U, 2U, 2U}};
  constexpr std::array<std::uint8_t,
                       detail::kPressureEnergyCandidateScratchFieldCount>
      scratch_components{{1U, 1U, 1U, 9U, 1U, 1U, 1U, 1U, 1U}};
  for (std::size_t index = 0U; index < diagnostic.scratch_fields.size();
       ++index) {
    const auto& field = diagnostic.scratch_fields[index];
    passed &= expect(
        field.candidate_field != field.trial_field &&
            field.candidate_revision != 0U && field.trial_revision != 0U &&
            field.candidate_revision != field.trial_revision &&
            field.candidate_storage == field.trial_storage &&
            field.candidate_revision_domain == field.trial_revision_domain &&
            field.components == scratch_components[index] &&
            field.ghost_width == scratch_ghosts[index] &&
            disjoint(field.candidate_begin, field.candidate_end,
                     field.trial_begin, field.trial_end),
        "candidate correction/material scratch cannot overwrite live workspace");
    for (const auto& primitive : diagnostic.fields)
      passed &= expect(
          disjoint(field.candidate_begin, field.candidate_end,
                   primitive.candidate_begin, primitive.candidate_end),
          "candidate material scratch is disjoint from candidate primitives");
    for (std::size_t other = index + 1U;
         other < diagnostic.scratch_fields.size(); ++other) {
      const auto& right = diagnostic.scratch_fields[other];
      passed &= expect(
          disjoint(field.candidate_begin, field.candidate_end,
                   right.candidate_begin, right.candidate_end),
          "candidate material workspaces do not alias one another");
    }
  }
  passed &= expect(
      diagnostic.coupled_state_halo != 0U &&
          diagnostic.candidate_state_halo != 0U &&
          diagnostic.coupled_state_halo != diagnostic.candidate_state_halo,
      "candidate halo owns an independent engine lineage");
  passed &= expect(
      diagnostic.coupled_thermal_halo != 0U &&
          diagnostic.candidate_thermal_halo != 0U &&
          diagnostic.coupled_thermal_halo !=
              diagnostic.candidate_thermal_halo &&
          diagnostic.coupled_thermal_halo != diagnostic.coupled_state_halo &&
          diagnostic.candidate_thermal_halo !=
              diagnostic.candidate_state_halo,
      "effective thermal ghosts own live and candidate halo lineages");
  passed &= expect(
      diagnostic.candidate_finalizer_state_halo != 0U &&
          diagnostic.candidate_finalizer_state_halo !=
              diagnostic.candidate_state_halo &&
          diagnostic.candidate_finalizer_state_halo !=
              diagnostic.coupled_state_halo,
      "candidate physical-boundary finalizer owns an exact independent halo");
  passed &= expect(
      diagnostic.correction_halo != 0U &&
          diagnostic.candidate_correction_halo != 0U &&
          diagnostic.correction_halo !=
              diagnostic.candidate_correction_halo,
      "candidate delta-pi owns an independent correction halo lineage");
  const auto same_halo_capacity = [](const HaloPlanStats& left,
                                     const HaloPlanStats& right) {
    return left.transport_peer_count == right.transport_peer_count &&
           left.local_peer_count == right.local_peer_count &&
           left.persistent_request_count == right.persistent_request_count &&
           left.send_capacity_doubles == right.send_capacity_doubles &&
           left.receive_capacity_doubles == right.receive_capacity_doubles &&
           left.maximum_messages_per_exchange ==
               right.maximum_messages_per_exchange &&
           left.maximum_bytes_per_exchange ==
               right.maximum_bytes_per_exchange &&
           left.maximum_tag == right.maximum_tag;
  };
  passed &= expect(
      same_halo_capacity(diagnostic.coupled_state_halo_plan,
                         diagnostic.candidate_state_halo_plan) &&
          same_halo_capacity(diagnostic.coupled_thermal_halo_plan,
                             diagnostic.candidate_thermal_halo_plan) &&
          same_halo_capacity(diagnostic.correction_halo_plan,
                             diagnostic.candidate_correction_halo_plan),
      "candidate halos freeze the live geometric reach under new field lineage");
  passed &= expect(
      diagnostic.workspace_flux_capacity.aligned_payload_allocations == 1U &&
          diagnostic.workspace_flux_capacity.replicas == 6U &&
          diagnostic.workspace_flux_capacity.directional_blocks == 18U &&
          diagnostic.final_flux_capacity.aligned_payload_allocations == 1U &&
          diagnostic.final_flux_capacity.replicas == 3U &&
          diagnostic.final_flux_capacity.directional_blocks == 9U &&
          diagnostic.field_storage_capacity.aligned_payload_allocations ==
              1U &&
          diagnostic.arena_doubles == plan.summary().arena_doubles,
      "all candidate state and flux capacity is allocated during freeze");
  for (std::size_t left = 0U; left < diagnostic.flux_replicas.size(); ++left) {
    const auto& flux = diagnostic.flux_replicas[left];
    passed &= expect(flux.storage != 0U && flux.revision_domain != 0U &&
                         flux.revision != 0U &&
                         flux.replica_begin < flux.replica_end &&
                         flux.bases[0U] != 0U && flux.bases[1U] != 0U &&
                         flux.bases[2U] != 0U,
                     "each workspace flux replica has frozen storage");
    for (std::size_t right = left + 1U;
         right < diagnostic.flux_replicas.size(); ++right)
      passed &= expect(
          disjoint(flux.replica_begin, flux.replica_end,
                   diagnostic.flux_replicas[right].replica_begin,
                   diagnostic.flux_replicas[right].replica_end),
          "candidate flux does not alias predictor or provisional flux");
  }
  ProductDriver driver;
  const char* lifecycle_stage = "create";
  if (status)
    status = ProductDriver::create(MPI_COMM_SELF, std::move(plan), driver);
  DriverInitialState initial;
  initial.pressure_reference = 101325.0;
  initial.temperature = 300.0;
  if (status) {
    lifecycle_stage = "initialize";
    status = driver.initialize(initial);
  }
  DriverStepReport step;
  if (status) {
    lifecycle_stage = "advance";
    status = driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, step);
    passed &= expect(status && step.accepted,
                     "candidate capacity remains allocation-free in an attempt");
  }
  RestartSnapshot snapshot;
  if (status) {
    lifecycle_stage = "snapshot";
    status = driver.committed_restart_snapshot(snapshot);
  }
  const auto& candidate_flux = diagnostic.flux_replicas[2U];
  if (status) {
    const std::array<ConstFaceFieldView, 3U> final_flux{{
        snapshot.final_mass_flux.x, snapshot.final_mass_flux.y,
        snapshot.final_mass_flux.z}};
    for (std::size_t axis = 0U; axis < final_flux.size(); ++axis)
      passed &= expect(
          final_flux[axis].base != nullptr &&
              reinterpret_cast<std::uintptr_t>(final_flux[axis].base) !=
                  candidate_flux.bases[axis] &&
              final_flux[axis].storage_identity != candidate_flux.storage,
          "candidate flux storage is independent of final flux authority");
  }
  if (!status)
    std::cerr << "candidate storage lifecycle status="
              << static_cast<unsigned>(status.code) << '/' << status.detail
              << " stage=" << lifecycle_stage << '\n';
  return expect(status && passed,
                "candidate state/halo/flux storage lineage is sealed");
}

struct CommittedDriverBits {
  double time{};
  double dt{};
  double pressure_reference{};
  std::uint64_t step{};
  std::vector<std::uint64_t> field_metadata;
  std::vector<std::uint64_t> field_values;
  std::array<std::vector<std::uint64_t>, 3U> flux_values;
  RevisionToken flux_revision{};
  FaceFluxCertificate flux_certificate{};
};

std::uint64_t double_bits(double value) {
  std::uint64_t bits = 0U;
  static_assert(sizeof(bits) == sizeof(value));
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

bool capture_committed_driver_bits(ProductDriver& driver,
                                   CommittedDriverBits& out) {
  RestartSnapshot snapshot;
  const Status status = driver.committed_restart_snapshot(snapshot);
  if (!status) return false;
  CommittedDriverBits candidate;
  candidate.time = snapshot.time;
  candidate.dt = snapshot.dt;
  candidate.pressure_reference = snapshot.pressure_reference;
  candidate.step = snapshot.step;
  for (std::size_t field = 0U; field < snapshot.fields.size; ++field) {
    const RestartFieldView& view = snapshot.fields.data[field];
    candidate.field_metadata.push_back(
        static_cast<std::uint64_t>(view.role));
    candidate.field_metadata.push_back(view.values.field);
    candidate.field_metadata.push_back(view.values.revision);
    candidate.field_metadata.push_back(view.values.storage_identity);
    candidate.field_metadata.push_back(view.values.revision_domain);
    for (std::uint8_t component = 0U;
         component < view.values.components; ++component)
      for (std::int32_t z = 0; z < view.values.interior.z; ++z)
        for (std::int32_t y = 0; y < view.values.interior.y; ++y)
          for (std::int32_t x = 0; x < view.values.interior.x; ++x)
            candidate.field_values.push_back(double_bits(
                view.values.unchecked({x, y, z}, component)));
  }
  const std::array<ConstFaceFieldView, 3U> flux{{
      snapshot.final_mass_flux.x, snapshot.final_mass_flux.y,
      snapshot.final_mass_flux.z}};
  for (std::size_t axis = 0U; axis < flux.size(); ++axis)
    for (std::int32_t z = 0; z < flux[axis].extents.z; ++z)
      for (std::int32_t y = 0; y < flux[axis].extents.y; ++y)
        for (std::int32_t x = 0; x < flux[axis].extents.x; ++x)
          candidate.flux_values[axis].push_back(
              double_bits(flux[axis].unchecked({x, y, z})));
  candidate.flux_revision = snapshot.final_mass_flux.revision;
  candidate.flux_certificate = snapshot.final_mass_flux.certificate;
  out = std::move(candidate);
  return true;
}

bool same_committed_driver_bits(const CommittedDriverBits& left,
                                const CommittedDriverBits& right) {
  return double_bits(left.time) == double_bits(right.time) &&
         double_bits(left.dt) == double_bits(right.dt) &&
         double_bits(left.pressure_reference) ==
             double_bits(right.pressure_reference) &&
         left.step == right.step &&
         left.field_metadata == right.field_metadata &&
         left.field_values == right.field_values &&
         left.flux_values == right.flux_values &&
         left.flux_revision == right.flux_revision &&
         left.flux_certificate == right.flux_certificate;
}

bool test_fatal_time_finish_is_recoverable() {
  ValidatedModel model = test::product_model({17, 11, 7});
  model.time.initial_dt = 4.0;
  model.time.maximum_dt = 4.0;
  CompiledCasePlan plan;
  Status status = ProductCompiler::compile(MPI_COMM_SELF, model, {}, plan);
  ProductDriver driver;
  if (status)
    status = ProductDriver::create(MPI_COMM_SELF, std::move(plan), driver);
  DriverInitialState initial;
  initial.pressure_reference = 98000.0;
  initial.temperature = 315.0;
  initial.start_time = 1.0e16;
  if (status) status = driver.initialize(initial);

  DriverStepReport accepted;
  detail::arm_pressure_energy_candidate_globalization_once_for_test();
  if (status)
    status = driver.advance({std::numeric_limits<double>::infinity(),
                             std::numeric_limits<double>::infinity(),
                             std::numeric_limits<double>::infinity(),
                             std::numeric_limits<double>::infinity(),
                             std::numeric_limits<double>::infinity()},
                            accepted);
  detail::PressureEnergyCandidateGlobalizationDiagnostic
      stationary_diagnostic;
  const bool captured_stationary_diagnostic =
      detail::pressure_energy_candidate_globalization_diagnostic_for_test(
          stationary_diagnostic);
  detail::clear_pressure_energy_candidate_globalization_for_test();
  CommittedDriverBits before;
  if (status && !capture_committed_driver_bits(driver, before))
    status = {StatusCode::invalid_plan, 1U};

  DriverStepReport rejected;
  Status no_progress = status;
  if (status)
    no_progress = driver.advance(
        {2.0e-8, 2.0e-8, 2.0e-8, 2.0e-8, 2.0e-8}, rejected);
  CommittedDriverBits after;
  const bool captured_after = capture_committed_driver_bits(driver, after);
  const bool fatal_rolled_back =
      captured_after && same_committed_driver_bits(before, after) &&
      double_bits(driver.pressure_reference()) ==
          double_bits(before.pressure_reference);
  DriverStepReport recovery;
  Status recovered = no_progress;
  if (captured_after)
    recovered = driver.advance({std::numeric_limits<double>::infinity(),
                                std::numeric_limits<double>::infinity(),
                                std::numeric_limits<double>::infinity(),
                                std::numeric_limits<double>::infinity(),
                                std::numeric_limits<double>::infinity()},
                               recovery);
  detail::PressureCorrectionWarmStartDiagnostic recovery_warm_start;
  const bool captured_recovery_warm_start =
      detail::pressure_correction_warm_start_diagnostic_for_test(
          recovery_warm_start);
  const bool rolled_back_and_recoverable =
      status && accepted.accepted && accepted.accepted_step == 1U &&
      captured_stationary_diagnostic &&
      stationary_diagnostic.production_candidate_loop &&
      stationary_diagnostic.baseline_commit &&
      !stationary_diagnostic.selection_valid &&
      stationary_diagnostic.replay_valid && stationary_diagnostic.committed &&
      stationary_diagnostic.corrector == 2U &&
      stationary_diagnostic.maximum_absolute_pressure_correction > 0.0 &&
      stationary_diagnostic.maximum_absolute_enthalpy_correction > 0.0 &&
      no_progress.code == StatusCode::invalid_plan &&
      no_progress.detail == 453U &&
      rejected.attempts == 1U &&
      rejected.failure.code == StatusCode::ok &&
      rejected.piso.pressure_solve_calls == 2U && fatal_rolled_back &&
      recovered &&
      recovery.accepted && recovery.attempts == 1U &&
      recovery.proposal.origin == StepOrigin::retry &&
      recovery.proposal.bdf.order == 1U && captured_recovery_warm_start &&
      recovery_warm_start.origin == StepOrigin::retry &&
      recovery_warm_start.attempt == 0U &&
      !recovery_warm_start.authority_available &&
      !recovery_warm_start.used;
  if (!rolled_back_and_recoverable) {
    std::cerr << "time-finish rollback setup="
              << static_cast<unsigned>(status.code) << ':' << status.detail
              << " accepted=" << accepted.accepted
              << " accepted-step=" << accepted.accepted_step
              << " stationary-diagnostic=" << captured_stationary_diagnostic
              << '/' << stationary_diagnostic.production_candidate_loop << '/'
              << stationary_diagnostic.baseline_commit << '/'
              << stationary_diagnostic.selection_valid << '/'
              << stationary_diagnostic.replay_valid << '/'
              << stationary_diagnostic.committed
              << " direction="
              << stationary_diagnostic.maximum_absolute_pressure_correction
              << '/'
              << stationary_diagnostic.maximum_absolute_enthalpy_correction
              << " residual="
              << stationary_diagnostic.baseline_normalized_continuity << '/'
              << stationary_diagnostic.baseline_normalized_energy
              << " finish=" << static_cast<unsigned>(no_progress.code)
              << ':' << no_progress.detail
              << " attempts=" << rejected.attempts
              << " solves="
              << static_cast<unsigned>(rejected.piso.pressure_solve_calls)
              << " captured=" << captured_after
              << " fatal-rolled-back=" << fatal_rolled_back
              << " recovery=" << static_cast<unsigned>(recovered.code)
              << ':' << recovered.detail
              << " recovery-accepted=" << recovery.accepted
              << " recovery-attempts=" << recovery.attempts
              << " recovery-origin="
              << static_cast<unsigned>(recovery.proposal.origin)
              << " recovery-order="
              << static_cast<unsigned>(recovery.proposal.bdf.order)
              << " recovery-warm=" << captured_recovery_warm_start << '/'
              << recovery_warm_start.authority_available << '/'
              << recovery_warm_start.used
              << '\n';
  }
  return expect(
      rolled_back_and_recoverable,
      "fatal time finish rolls back the product transaction, retires its "
      "proposal, and permits a fresh BE recovery advance");
}

bool test_product_step_role_matrix() {
  ValidatedModel model = test::product_model({17, 11, 7});
  CompiledCasePlan plan;
  Status status = ProductCompiler::compile(MPI_COMM_SELF, model, {}, plan);
  ProductDriver driver;
  if (status)
    status = ProductDriver::create(MPI_COMM_SELF, std::move(plan), driver);
  DriverInitialState initial;
  initial.pressure_reference = 98000.0;
  initial.temperature = 315.0;
  if (status) status = driver.initialize(initial);
  DriverStepReport first;
  if (status)
    status = driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, first);
  DriverStepReport second;
  if (status)
    status = driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, second);
  RestartSnapshot restart;
  if (status) status = driver.committed_restart_snapshot(restart);
  const auto coupled_c1_role = [](const LinearSolveResult& result) {
    return result.status && result.convergence_audits == 0U &&
           result.recycle_offered_directions == 0U &&
           result.recycle_retained_directions == 0U &&
           result.recycle_operator_applies == 0U &&
           result.recycle_reduction_calls == 0U &&
           !result.recycle_projection_attempted &&
           !result.recycle_projection_accepted &&
           result.recycle_projected_true_residual == 0.0 &&
           result.recycle_capture_cycle_attempts >=
               result.recycle_cycle_corrections &&
           result.recycle_capture_vector_passes ==
               2U * result.recycle_capture_cycle_attempts &&
           result.recycle_capture_reduction_calls ==
               result.recycle_capture_cycle_attempts &&
           result.recycle_capture_blocking_operations ==
               2U * result.recycle_capture_cycle_attempts;
  };
  const auto coupled_c2_role = [](const LinearSolveResult& result) {
    const bool projected_direction =
        result.recycle_retained_directions == 1U &&
        result.recycle_operator_applies == 2U &&
        result.recycle_reduction_calls == 5U &&
        result.recycle_projection_attempted &&
        result.recycle_projection_accepted &&
        result.recycle_projected_true_residual <
            result.initial_true_residual;
    // A terminal uniform baseline still offers C1's direction to C2, but the
    // zero-residual solve has no direction to retain or project.  Its typed
    // stationary commit is checked independently by the time-finish test.
    const bool stationary_no_op =
        result.recycle_retained_directions == 0U &&
        result.recycle_operator_applies == 0U &&
        result.recycle_reduction_calls == 0U &&
        !result.recycle_projection_attempted &&
        !result.recycle_projection_accepted;
    return result.status && result.convergence_audits == 0U &&
           result.convergence_rejections == 0U &&
           result.recycle_cycle_corrections == 0U &&
           result.recycle_capture_cycle_attempts == 0U &&
           result.recycle_capture_vector_passes == 0U &&
           result.recycle_capture_reduction_calls == 0U &&
           result.recycle_capture_blocking_operations == 0U &&
           result.recycle_offered_directions == 1U &&
           (projected_direction || stationary_no_op);
  };
  bool uniform_state = restart.fields.size == 3U;
  for (std::size_t field = 0U;
       field < restart.fields.size && uniform_state; ++field) {
    const RestartFieldView& snapshot = restart.fields.data[field];
    double minimum = std::numeric_limits<double>::infinity();
    double maximum = -std::numeric_limits<double>::infinity();
    double maximum_absolute = 0.0;
    for (std::uint8_t component = 0U;
         component < snapshot.values.components; ++component)
      for (std::int32_t z = 0; z < snapshot.values.interior.z; ++z)
        for (std::int32_t y = 0; y < snapshot.values.interior.y; ++y)
          for (std::int32_t x = 0; x < snapshot.values.interior.x; ++x) {
            const double value =
                snapshot.values.unchecked({x, y, z}, component);
            minimum = std::min(minimum, value);
            maximum = std::max(maximum, value);
            maximum_absolute = std::max(maximum_absolute, std::abs(value));
          }
    if (snapshot.role == RestartFieldRole::velocity)
      uniform_state &= maximum_absolute < 1.0e-10;
    else if (snapshot.role == RestartFieldRole::pressure_perturbation)
      uniform_state &= maximum_absolute < 1.0e-8;
    else if (snapshot.role == RestartFieldRole::enthalpy)
      uniform_state &= maximum - minimum < 1.0e-8;
  }
  const bool role_matrix =
      static_cast<bool>(status) && first.accepted && second.accepted &&
          first.attempts == 1U && second.attempts == 1U &&
          first.piso.pressure_solve_calls == 2U &&
          second.piso.pressure_solve_calls == 2U &&
          first.failure.code == StatusCode::ok &&
          second.failure.code == StatusCode::ok &&
          coupled_c1_role(first.piso.pressure[0U]) &&
          coupled_c2_role(first.piso.pressure[1U]) &&
          coupled_c1_role(second.piso.pressure[0U]) &&
          coupled_c2_role(second.piso.pressure[1U]) && uniform_state &&
          second.accepted_step == 2U &&
          restart.final_mass_flux.certificate.valid() &&
          restart.final_mass_flux.revision == second.piso.final_flux_revision;
  if (!role_matrix) {
    const auto dump_recycle = [](const LinearSolveResult& result) {
      std::cerr << " {term=" << static_cast<unsigned>(result.termination)
                << " offered=" << result.recycle_offered_directions
                << " retained=" << result.recycle_retained_directions
                << " op=" << result.recycle_operator_applies
                << " red=" << result.recycle_reduction_calls
                << " proj=" << result.recycle_projection_attempted
                << '/' << result.recycle_projection_accepted
                << " cycle=" << result.recycle_cycle_corrections
                << " capture=" << result.recycle_capture_cycle_attempts
                << '/' << result.recycle_capture_vector_passes
                << '/' << result.recycle_capture_reduction_calls << '}';
    };
    std::cerr << "product-step status="
              << static_cast<unsigned>(status.code)
              << " detail=" << status.detail
              << " first.accepted=" << first.accepted
              << " first.attempts=" << first.attempts
              << " first.failed_stage=" << first.failed_stage
              << " first.failure="
              << static_cast<unsigned>(first.failure.code)
              << ':' << first.failure.detail
              << " first.solves="
              << static_cast<unsigned>(first.piso.pressure_solve_calls)
              << " second.accepted=" << second.accepted
              << " second.attempts=" << second.attempts
              << " second.failed_stage=" << second.failed_stage
              << " second.failure="
              << static_cast<unsigned>(second.failure.code)
              << ':' << second.failure.detail
              << " second.solves="
              << static_cast<unsigned>(second.piso.pressure_solve_calls)
              << " recycle="
              << coupled_c1_role(first.piso.pressure[0U]) << ','
              << coupled_c2_role(first.piso.pressure[1U]) << ','
              << coupled_c1_role(second.piso.pressure[0U]) << ','
              << coupled_c2_role(second.piso.pressure[1U])
              << " uniform=" << uniform_state
              << " accepted_step=" << second.accepted_step
              << " restart.valid="
              << restart.final_mass_flux.certificate.valid()
              << " restart.rev=" << restart.final_mass_flux.revision
              << " report.rev=" << second.piso.final_flux_revision
              << '\n';
    dump_recycle(first.piso.pressure[0U]);
    dump_recycle(first.piso.pressure[1U]);
    dump_recycle(second.piso.pressure[0U]);
    dump_recycle(second.piso.pressure[1U]);
    std::cerr << '\n';
    for (std::size_t field = 0U; field < restart.fields.size; ++field) {
      const RestartFieldView& snapshot = restart.fields.data[field];
      double minimum = std::numeric_limits<double>::infinity();
      double maximum = -std::numeric_limits<double>::infinity();
      for (std::uint8_t component = 0U;
           component < snapshot.values.components; ++component)
        for (std::int32_t z = 0; z < snapshot.values.interior.z; ++z)
          for (std::int32_t y = 0; y < snapshot.values.interior.y; ++y)
            for (std::int32_t x = 0; x < snapshot.values.interior.x; ++x) {
              const double value =
                  snapshot.values.unchecked({x, y, z}, component);
              minimum = std::min(minimum, value);
              maximum = std::max(maximum, value);
            }
      std::cerr << " restart-field role="
                << static_cast<unsigned>(snapshot.role)
                << " min=" << minimum << " max=" << maximum << '\n';
    }
  }
  return expect(
      role_matrix,
      "product step report keeps C1/C2 roles, two solves and committed flux");
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 2;
  const bool passed = test_incomplete_registration() &&
                      test_pressure_energy_temporal_operand_scale() &&
                      test_pressure_inexact_forcing_policy() &&
                      test_pressure_aitken_initial_alpha() &&
                      test_freeze() &&
                      test_live_thermal_halo_resource_contract() &&
                      test_pressure_energy_restart_schema() &&
                      test_pressure_energy_candidate_storage_lineage() &&
                      test_fatal_time_finish_is_recoverable() &&
                      test_product_step_role_matrix();
  MPI_Finalize();
  return passed ? 0 : 1;
}
