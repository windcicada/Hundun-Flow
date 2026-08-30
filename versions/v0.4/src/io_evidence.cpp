// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_io.hpp"

#include "io_output_detail.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <locale>
#include <new>
#include <sstream>
#include <string>

namespace hundun::v04 {
namespace {

std::string_view termination_name(LinearTermination termination) noexcept {
  switch (termination) {
    case LinearTermination::converged:
      return "converged";
    case LinearTermination::zero_rhs:
      return "zero_rhs";
    case LinearTermination::maximum_iterations:
      return "maximum_iterations";
    case LinearTermination::breakdown:
      return "breakdown";
    case LinearTermination::non_finite:
      return "non_finite";
    case LinearTermination::operator_failure:
      return "operator_failure";
    case LinearTermination::preconditioner_failure:
      return "preconditioner_failure";
    case LinearTermination::collective_failure:
      return "collective_failure";
    case LinearTermination::invalid_plan:
      return "invalid_plan";
    case LinearTermination::convergence_audit_failure:
      return "convergence_audit_failure";
  }
  return "invalid";
}

std::string_view pressure_solve_contract_name(
    RuntimePressureSolveContract contract) noexcept {
  switch (contract) {
    case RuntimePressureSolveContract::pressure_continuity:
      return "pressure_continuity";
    case RuntimePressureSolveContract::continuity_energy_coupled:
      return "continuity_energy_coupled";
    case RuntimePressureSolveContract::invalid:
      break;
  }
  return "invalid";
}

std::string_view pressure_energy_refinement_termination_name(
    RuntimePressureEnergyRefinementTermination termination) noexcept {
  switch (termination) {
    case RuntimePressureEnergyRefinementTermination::
        component_residuals_converged:
      return "component_residuals_converged";
    case RuntimePressureEnergyRefinementTermination::
        iteration_capacity_exhausted:
      return "iteration_capacity_exhausted";
    case RuntimePressureEnergyRefinementTermination::rejected_candidate:
      return "rejected_candidate";
    case RuntimePressureEnergyRefinementTermination::none:
      return "none";
  }
  return "invalid";
}

bool accepted_terminal_metric(double residual, double tolerance,
                              bool enabled_required = true) noexcept {
  if (!std::isfinite(residual) || residual < 0.0 ||
      !std::isfinite(tolerance) || tolerance < 0.0 ||
      (enabled_required && tolerance == 0.0))
    return false;
  // A zero tolerance is the V2 pressure-continuity convention for a disabled
  // energy gate.  Coupled V3 records never enter this branch.
  return tolerance == 0.0 || residual <= tolerance;
}

std::string_view predictor_constraint_name(
    std::uint8_t constraint) noexcept {
  switch (constraint) {
    case 0U:
      return "none";
    case 1U:
      return "density";
    case 2U:
      return "independent_species";
    case 3U:
      return "dependent_species";
    case 4U:
      return "enthalpy_lower";
    case 5U:
      return "enthalpy_upper";
  }
  return "invalid";
}

std::string_view predictor_low_state_name(std::uint8_t low_state) noexcept {
  switch (low_state) {
    case 0U:
      return "none";
    case 1U:
      return "bdf_accepted_rate";
    case 2U:
      return "scaled_euler";
    case 3U:
      return "implicit_upwind";
    case 4U:
      return "implicit_upwind_source_limited";
    case 5U:
      return "bdf_accepted_rate_homotopy";
    case 6U:
      return "bdf_local_donor_flux";
    case 7U:
      return "bdf_local_donor_flux_source_limited";
  }
  return "invalid";
}

Status validate_record(const IoServicePlan& services,
                       const RuntimeEvidenceRecord& record) noexcept {
  const bool pressure_continuity_contract =
      record.pressure_solve_contract ==
      RuntimePressureSolveContract::pressure_continuity;
  const bool continuity_energy_coupled_contract =
      record.pressure_solve_contract ==
      RuntimePressureSolveContract::continuity_energy_coupled;
  const RuntimeTerminalPhysicalAudit& terminal =
      record.terminal_physical_audit;
  const bool valid_terminal_physical_audit =
      terminal.present && terminal.final_flux_revision != 0U &&
      accepted_terminal_metric(terminal.eos_residual,
                               terminal.eos_tolerance) &&
      accepted_terminal_metric(terminal.continuity_residual,
                               terminal.continuity_tolerance) &&
      accepted_terminal_metric(terminal.energy_residual,
                               terminal.energy_tolerance,
                               continuity_energy_coupled_contract) &&
      accepted_terminal_metric(terminal.closed_mass_residual,
                               terminal.closed_mass_tolerance) &&
      accepted_terminal_metric(terminal.gauge_residual,
                               terminal.gauge_tolerance);
  const bool valid_temporal_method =
      (record.requested_bdf_order == 1U ||
       record.requested_bdf_order == 2U) &&
      (record.bdf_order == 1U || record.bdf_order == 2U) &&
      (record.temporal_method_fallback
           ? record.requested_bdf_order == 2U && record.bdf_order == 1U &&
                 record.thermophysical_predictor_calls == 2U
           : record.requested_bdf_order == record.bdf_order &&
                 record.thermophysical_predictor_calls == 1U);
  const bool valid_fast_predictor_collectives =
      record.temporal_method_fallback
          ? record.predictor_blocking_collectives > 1U
          : record.predictor_blocking_collectives == 1U;
  const bool enthalpy_endpoint_converged =
      record.predictor_enthalpy_endpoint.status &&
      (record.predictor_enthalpy_endpoint.termination ==
           LinearTermination::converged ||
       record.predictor_enthalpy_endpoint.termination ==
           LinearTermination::zero_rhs);
  const bool uses_enthalpy_endpoint = record.predictor_low_state == 3U ||
                                      record.predictor_low_state == 4U;
  const bool valid_enthalpy_endpoint =
      uses_enthalpy_endpoint
          ? record.predictor_enthalpy_solve_calls == 1U &&
                enthalpy_endpoint_converged &&
                record.predictor_enthalpy_endpoint.operator_applies > 0U &&
                record.predictor_enthalpy_endpoint.reduction_calls > 0U &&
                std::isfinite(record.predictor_enthalpy_endpoint_alpha) &&
                record.predictor_enthalpy_endpoint_alpha >= 0.0 &&
                record.predictor_enthalpy_endpoint_alpha <= 1.0 &&
                (record.predictor_low_state != 3U ||
                 record.predictor_enthalpy_endpoint_alpha == 1.0) &&
                (record.predictor_low_state != 4U ||
                 record.predictor_enthalpy_endpoint_alpha < 1.0)
          : record.predictor_enthalpy_solve_calls == 0U &&
                record.predictor_enthalpy_endpoint_alpha == 1.0;
  const bool valid_bdf_endpoint =
      std::isfinite(record.predictor_bdf_endpoint_alpha) &&
      record.predictor_bdf_endpoint_alpha >= 0.0 &&
      record.predictor_bdf_endpoint_alpha <= 1.0 &&
      (record.predictor_low_state == 5U
           ? record.predictor_bdf_endpoint_alpha < 1.0
           : record.predictor_bdf_endpoint_alpha == 1.0);
  const bool valid_source_endpoint =
      std::isfinite(record.predictor_source_endpoint_alpha) &&
      record.predictor_source_endpoint_alpha >= 0.0 &&
      record.predictor_source_endpoint_alpha <= 1.0 &&
      (record.predictor_low_state == 7U
           ? record.predictor_source_endpoint_alpha < 1.0
           : record.predictor_source_endpoint_alpha == 1.0);
  const bool valid_momentum_predictor =
      std::isfinite(record.momentum_predictor_theta) &&
      record.momentum_predictor_theta >= 0.0 &&
      record.momentum_predictor_theta <= 1.0 &&
      record.momentum_predictor_activations ==
          (record.momentum_predictor_limited ? 1U : 0U) &&
      record.momentum_predictor_limited ==
          (record.momentum_predictor_theta < 1.0);
  const bool valid_predictor =
      std::isfinite(record.predictor_theta) &&
      record.predictor_theta >= 0.0 && record.predictor_theta <= 1.0 &&
      std::isfinite(record.predictor_mass_flux_scale) &&
      record.predictor_mass_flux_scale >= 0.0 &&
      record.predictor_mass_flux_scale <= 1.0 &&
      (!record.predictor_limited
           ? record.predictor_theta == 1.0 &&
                 record.predictor_mass_flux_scale == 1.0 &&
                 record.predictor_low_margin == 0.0 &&
                 record.predictor_high_margin == 0.0 &&
                 record.predictor_low_order_substeps == 0U &&
                 record.predictor_low_order_halo_exchanges == 0U &&
                 valid_fast_predictor_collectives &&
                 record.predictor_low_order_transport_passes == 0U &&
                 record.predictor_limiting_rank == -1 &&
                 record.predictor_constraint == 0U &&
                 record.predictor_low_state == 0U
           : record.predictor_theta < 1.0 &&
                 std::isfinite(record.predictor_low_margin) &&
                 record.predictor_low_margin >= 0.0 &&
                 std::isfinite(record.predictor_high_margin) &&
                 record.predictor_high_margin < 0.0 &&
                 record.predictor_low_order_substeps > 0U &&
                 ((record.predictor_low_state == 6U ||
                   record.predictor_low_state == 7U)
                      ? record.predictor_low_order_halo_exchanges == 1U
                      : record.predictor_low_order_halo_exchanges ==
                            record.predictor_low_order_substeps - 1U) &&
                 record.predictor_blocking_collectives > 1U &&
                 record.predictor_low_order_transport_passes > 0U &&
                 record.predictor_limiting_cell_x >= 0 &&
                 record.predictor_limiting_cell_y >= 0 &&
                 record.predictor_limiting_cell_z >= 0 &&
                 record.predictor_limiting_rank >= 0 &&
                 record.predictor_constraint >= 1U &&
                 record.predictor_constraint <= 5U &&
                 record.predictor_low_state >= 1U &&
                 record.predictor_low_state <= 7U &&
                 (record.predictor_low_state != 1U ||
                  record.predictor_mass_flux_scale == 1.0) &&
                 (record.predictor_low_state != 5U ||
                  record.predictor_mass_flux_scale == 1.0) &&
                 ((record.predictor_low_state != 6U &&
                   record.predictor_low_state != 7U) ||
                  record.predictor_mass_flux_scale >= 0.0) &&
                 ((record.predictor_low_state != 3U &&
                   record.predictor_low_state != 4U) ||
                  record.predictor_mass_flux_scale == 1.0)) &&
      valid_enthalpy_endpoint && valid_bdf_endpoint && valid_source_endpoint;
  const std::size_t refinement_count =
      record.pressure_energy_refinement_solve_calls;
  bool valid_refinement =
      refinement_count <= record.pressure_energy_refinement.size() &&
      record.pressure_energy_refinement_termination ==
          RuntimePressureEnergyRefinementTermination::
              component_residuals_converged;
  RevisionToken refinement_target_generation = 0U;
  std::uint64_t refinement_linear_iterations = 0U;
  for (std::size_t index = 0U; valid_refinement && index < refinement_count;
       ++index) {
    const RuntimePressureEnergyRefinementSolve& entry =
        record.pressure_energy_refinement[index];
    const LinearSolveResult& solve = entry.solve;
    const bool accepted_termination =
        solve.termination == LinearTermination::converged ||
        solve.termination == LinearTermination::zero_rhs;
    const bool finite = std::isfinite(solve.initial_true_residual) &&
                        solve.initial_true_residual >= 0.0 &&
                        std::isfinite(solve.final_true_residual) &&
                        solve.final_true_residual >= 0.0 &&
                        std::isfinite(solve.recursive_residual) &&
                        solve.recursive_residual >= 0.0;
    const bool valid_capture_counts =
        solve.recycle_capture_cycle_attempts >=
            solve.recycle_cycle_corrections &&
        solve.recycle_capture_cycle_attempts <=
            std::numeric_limits<std::uint64_t>::max() / 2U &&
        solve.recycle_capture_vector_passes ==
            2U * solve.recycle_capture_cycle_attempts &&
        solve.recycle_capture_reduction_calls ==
            solve.recycle_capture_cycle_attempts &&
        solve.recycle_capture_blocking_operations ==
            2U * solve.recycle_capture_cycle_attempts;
    const bool valid_capture_role =
        solve.convergence_audits == 0U &&
        solve.convergence_rejections == 0U &&
        solve.final_convergence_metric == 0.0 &&
        solve.convergence_limit == 0.0 &&
        solve.recycle_offered_directions == 0U &&
        solve.recycle_retained_directions == 0U &&
        solve.recycle_operator_applies == 0U &&
        solve.recycle_reduction_calls == 0U &&
        !solve.recycle_projection_attempted &&
        !solve.recycle_projection_accepted &&
        solve.recycle_projected_true_residual == 0.0;
    const bool same_target =
        entry.target_generation != 0U &&
        (index == 0U ||
         entry.target_generation == refinement_target_generation);
    bool unique_lineage = entry.collective_lineage != 0U;
    for (std::size_t prior = 0U; unique_lineage && prior < index; ++prior)
      unique_lineage =
          record.pressure_energy_refinement[prior].collective_lineage !=
          entry.collective_lineage;
    valid_refinement =
        solve.status && accepted_termination && finite &&
        valid_capture_counts && valid_capture_role && same_target &&
        unique_lineage && entry.ordinal == index + 1U;
    if (index == 0U)
      refinement_target_generation = entry.target_generation;
    if (solve.iterations >
        std::numeric_limits<std::uint64_t>::max() -
            refinement_linear_iterations) {
      valid_refinement = false;
    } else {
      refinement_linear_iterations += solve.iterations;
    }
  }
  valid_refinement =
      valid_refinement &&
      record.linear_iterations >= refinement_linear_iterations;
  if (detail::output_service(services, RuntimeServiceKind::evidence) == nullptr ||
      record.build == 0U || record.binary == 0U || record.case_model == 0U ||
      record.product == 0U || record.step == 0U ||
      !std::isfinite(record.time) ||
      !valid_temporal_method ||
      ((record.startup || record.retry || record.restart_recovery) &&
       record.statistics_eligible) ||
      record.pressure_solve_calls != 2U ||
      (!pressure_continuity_contract &&
       !continuity_energy_coupled_contract) ||
      !valid_refinement ||
      !valid_terminal_physical_audit ||
      record.momentum_predictor_solve_calls != 3U ||
      record.blocking_collectives < record.predictor_blocking_collectives ||
      !valid_momentum_predictor || !valid_predictor ||
      (record.stages.size != 0U && record.stages.data == nullptr))
    return {StatusCode::invalid_plan, detail::kOutputInput};
  for (std::size_t corrector = 0U; corrector < record.pressure.size();
       ++corrector) {
    const LinearSolveResult& solve = record.pressure[corrector];
    const bool accepted_termination =
        solve.termination == LinearTermination::converged ||
        solve.termination == LinearTermination::zero_rhs;
    const bool finite = std::isfinite(solve.initial_true_residual) &&
                        solve.initial_true_residual >= 0.0 &&
                        std::isfinite(solve.final_true_residual) &&
                        solve.final_true_residual >= 0.0 &&
                        std::isfinite(solve.recursive_residual) &&
                        solve.recursive_residual >= 0.0;
    const bool no_linear_audit =
        solve.convergence_audits == 0U &&
        solve.convergence_rejections == 0U &&
        solve.final_convergence_metric == 0.0 &&
        solve.convergence_limit == 0.0;
    const bool accepted_linear_audit =
        solve.convergence_audits > 0U &&
        solve.convergence_rejections <= solve.convergence_audits &&
        std::isfinite(solve.final_convergence_metric) &&
        solve.final_convergence_metric >= 0.0 &&
        std::isfinite(solve.convergence_limit) &&
        solve.convergence_limit > 0.0 &&
        solve.final_convergence_metric <= solve.convergence_limit;
    const bool valid_audit =
        corrector == 0U
            ? no_linear_audit
            : (pressure_continuity_contract ? accepted_linear_audit
                                            : no_linear_audit);
    const bool valid_recycle_counts =
        solve.recycle_offered_directions <= kLinearRecycleMaximumDirections &&
        solve.recycle_retained_directions <=
            solve.recycle_offered_directions &&
        solve.recycle_capture_cycle_attempts >=
            solve.recycle_cycle_corrections &&
        solve.recycle_capture_cycle_attempts <=
            std::numeric_limits<std::uint64_t>::max() / 2U &&
        solve.recycle_capture_vector_passes ==
            2U * solve.recycle_capture_cycle_attempts &&
        solve.recycle_capture_reduction_calls ==
            solve.recycle_capture_cycle_attempts &&
        solve.recycle_capture_blocking_operations ==
            2U * solve.recycle_capture_cycle_attempts;
    const bool valid_recycle_role =
        corrector == 0U
            ? solve.recycle_offered_directions == 0U &&
                  solve.recycle_retained_directions == 0U &&
                  solve.recycle_operator_applies == 0U &&
                  solve.recycle_reduction_calls == 0U &&
                  !solve.recycle_projection_attempted &&
                  !solve.recycle_projection_accepted &&
                  solve.recycle_projected_true_residual == 0.0
            : solve.recycle_capture_vector_passes == 0U &&
                  solve.recycle_cycle_corrections == 0U &&
                  solve.recycle_capture_cycle_attempts == 0U &&
                  solve.recycle_capture_reduction_calls == 0U &&
                  solve.recycle_capture_blocking_operations == 0U &&
                  ((!solve.recycle_projection_attempted &&
                    solve.recycle_retained_directions == 0U &&
                    solve.recycle_operator_applies == 0U &&
                    solve.recycle_reduction_calls == 0U &&
                    !solve.recycle_projection_accepted &&
                    solve.recycle_projected_true_residual == 0.0) ||
                   (solve.recycle_projection_attempted &&
                    solve.recycle_offered_directions > 0U &&
                    solve.recycle_reduction_calls > 0U &&
                    solve.recycle_operator_applies ==
                        solve.recycle_offered_directions +
                            (solve.recycle_retained_directions != 0U ? 1U : 0U) &&
                    std::isfinite(solve.recycle_projected_true_residual) &&
                    solve.recycle_projected_true_residual >= 0.0 &&
                    (solve.recycle_retained_directions != 0U ||
                     (!solve.recycle_projection_accepted &&
                      solve.recycle_reduction_calls ==
                          solve.recycle_offered_directions &&
                      solve.recycle_projected_true_residual == 0.0))));
    const bool projection_improved =
        solve.recycle_projection_attempted &&
        solve.recycle_retained_directions > 0U &&
        solve.recycle_projected_true_residual <
            solve.initial_true_residual;
    const bool valid_recycle_admission =
        solve.recycle_projection_accepted == projection_improved &&
        (!solve.recycle_projection_accepted || solve.iterations != 0U ||
         solve.final_true_residual ==
             solve.recycle_projected_true_residual);
    if (!solve.status || !accepted_termination || !finite || !valid_audit ||
        !valid_recycle_counts || !valid_recycle_role ||
        !valid_recycle_admission)
      return {StatusCode::invalid_plan, detail::kOutputInput};
  }
  for (const LinearSolveResult& solve : record.momentum_predictor) {
    const bool accepted_termination =
        solve.termination == LinearTermination::converged ||
        solve.termination == LinearTermination::zero_rhs;
    const bool finite = std::isfinite(solve.initial_true_residual) &&
                        solve.initial_true_residual >= 0.0 &&
                        std::isfinite(solve.final_true_residual) &&
                        solve.final_true_residual >= 0.0 &&
                        std::isfinite(solve.recursive_residual) &&
                        solve.recursive_residual >= 0.0;
    if (!solve.status || !accepted_termination || !finite ||
        solve.final_true_residual > solve.initial_true_residual ||
        solve.convergence_audits != 0U ||
        solve.convergence_rejections != 0U ||
        solve.recycle_offered_directions != 0U ||
        solve.recycle_retained_directions != 0U ||
        solve.recycle_operator_applies != 0U ||
        solve.recycle_reduction_calls != 0U ||
        solve.recycle_projection_attempted ||
        solve.recycle_projection_accepted ||
        solve.recycle_projected_true_residual != 0.0 ||
        solve.recycle_cycle_corrections != 0U ||
        solve.recycle_capture_vector_passes != 0U ||
        solve.recycle_capture_cycle_attempts != 0U ||
        solve.recycle_capture_reduction_calls != 0U ||
        solve.recycle_capture_blocking_operations != 0U)
      return {StatusCode::invalid_plan, detail::kOutputInput};
  }
  const std::uint64_t expected_offered_directions =
      std::min<std::uint64_t>(
          record.pressure[0U].recycle_cycle_corrections,
          static_cast<std::uint64_t>(kLinearRecycleMaximumDirections));
  if (record.pressure[1U].recycle_offered_directions !=
      expected_offered_directions)
    return {StatusCode::invalid_plan, detail::kOutputInput};
  for (std::size_t index = 0U; index < record.stages.size; ++index) {
    const StageTimingRecord stage = record.stages.data[index];
    if (stage.stage == 0U ||
        stage.minimum_nanoseconds > stage.mean_nanoseconds ||
        stage.mean_nanoseconds > stage.maximum_nanoseconds)
      return {StatusCode::invalid_plan, detail::kOutputInput};
    for (std::size_t prior = 0U; prior < index; ++prior)
      if (record.stages.data[prior].stage == stage.stage)
        return {StatusCode::invalid_plan, detail::kOutputInput};
  }
  return {};
}

std::string encode_record(const RuntimeEvidenceRecord& record) {
  std::ostringstream json;
  json.imbue(std::locale::classic());
  json << std::setprecision(17)
       << "{\"schema\":\"HUNDUN_V04_EVIDENCE_V4\""
       << ",\"build\":" << record.build
       << ",\"binary\":" << record.binary
       << ",\"case\":" << record.case_model
       << ",\"stl\":" << record.stl
       << ",\"product\":" << record.product
       << ",\"cpu_plan\":" << record.cpu_plan
       << ",\"step\":" << record.step << ",\"time\":" << record.time
       << ",\"requested_bdf_order\":"
       << static_cast<unsigned>(record.requested_bdf_order)
       << ",\"bdf_order\":" << static_cast<unsigned>(record.bdf_order)
       << ",\"temporal_method_fallback\":"
       << (record.temporal_method_fallback ? "true" : "false")
       << ",\"thermophysical_predictor_calls\":"
       << static_cast<unsigned>(record.thermophysical_predictor_calls)
       << ",\"launcher_ns\":" << record.launcher_nanoseconds
       << ",\"max_rank_step_ns\":"
       << record.maximum_rank_step_nanoseconds
       << ",\"max_rank_rss_bytes\":" << record.maximum_rank_rss_bytes
       << ",\"max_node_rss_bytes\":" << record.maximum_node_rss_bytes
       << ",\"structured_messages\":" << record.structured_messages
       << ",\"structured_bytes\":" << record.structured_bytes
       << ",\"ibm_messages\":" << record.ibm_messages
       << ",\"ibm_bytes\":" << record.ibm_bytes
       << ",\"blocking_collectives\":" << record.blocking_collectives
       << ",\"nonblocking_collectives\":"
       << record.nonblocking_collectives
       << ",\"reduction_ns\":" << record.reduction_nanoseconds
       << ",\"linear_iterations\":" << record.linear_iterations
       << ",\"exact_numeric_refills\":" << record.exact_numeric_refills
       << ",\"coarse_numeric_refills\":" << record.coarse_numeric_refills
       << ",\"preconditioner_setups\":" << record.preconditioner_setups
       << ",\"preconditioner_reuses\":" << record.preconditioner_reuses
       << ",\"heap_allocations\":" << record.heap_allocations
       << ",\"pressure_solve_calls\":"
       << static_cast<unsigned>(record.pressure_solve_calls)
       << ",\"pressure_solve_contract\":\""
       << pressure_solve_contract_name(record.pressure_solve_contract)
       << "\""
       << ",\"pressure_energy_refinement_solve_calls\":"
       << static_cast<unsigned>(
              record.pressure_energy_refinement_solve_calls)
       << ",\"pressure_energy_refinement_termination\":\""
       << pressure_energy_refinement_termination_name(
              record.pressure_energy_refinement_termination)
       << "\""
       << ",\"terminal_physical_audit\":{\"present\":"
       << (record.terminal_physical_audit.present ? "true" : "false")
       << ",\"final_flux_revision\":"
       << record.terminal_physical_audit.final_flux_revision
       << ",\"eos_residual\":"
       << record.terminal_physical_audit.eos_residual
       << ",\"eos_tolerance\":"
       << record.terminal_physical_audit.eos_tolerance
       << ",\"continuity_residual\":"
       << record.terminal_physical_audit.continuity_residual
       << ",\"continuity_tolerance\":"
       << record.terminal_physical_audit.continuity_tolerance
       << ",\"energy_residual\":"
       << record.terminal_physical_audit.energy_residual
       << ",\"energy_tolerance\":"
       << record.terminal_physical_audit.energy_tolerance
       << ",\"closed_mass_residual\":"
       << record.terminal_physical_audit.closed_mass_residual
       << ",\"closed_mass_tolerance\":"
       << record.terminal_physical_audit.closed_mass_tolerance
       << ",\"gauge_residual\":"
       << record.terminal_physical_audit.gauge_residual
       << ",\"gauge_tolerance\":"
       << record.terminal_physical_audit.gauge_tolerance << '}'
       << ",\"momentum_predictor_solve_calls\":"
       << static_cast<unsigned>(record.momentum_predictor_solve_calls)
       << ",\"momentum_predictor_limiter\":{\"limited\":"
       << (record.momentum_predictor_limited ? "true" : "false")
       << ",\"theta\":" << record.momentum_predictor_theta
       << ",\"activations\":"
       << record.momentum_predictor_activations << '}'
       << ",\"thermophysical_predictor\":{\"limited\":"
       << (record.predictor_limited ? "true" : "false")
       << ",\"theta\":" << record.predictor_theta
       << ",\"low_state\":\""
       << predictor_low_state_name(record.predictor_low_state)
       << "\",\"mass_flux_scale\":"
       << record.predictor_mass_flux_scale
       << ",\"constraint\":\""
       << predictor_constraint_name(record.predictor_constraint)
       << "\",\"limiting_rank\":"
       << record.predictor_limiting_rank
       << ",\"limiting_global_cell\":["
       << record.predictor_limiting_cell_x << ','
       << record.predictor_limiting_cell_y << ','
       << record.predictor_limiting_cell_z
       << "],\"low_margin\":"
       << record.predictor_low_margin
       << ",\"high_margin\":"
       << record.predictor_high_margin
       << ",\"low_order_substeps\":"
       << record.predictor_low_order_substeps
       << ",\"low_order_transport_passes\":"
       << record.predictor_low_order_transport_passes
       << ",\"low_order_halo_exchanges\":"
       << record.predictor_low_order_halo_exchanges
       << ",\"blocking_collectives\":"
       << record.predictor_blocking_collectives
       << ",\"enthalpy_endpoint_alpha\":"
       << record.predictor_enthalpy_endpoint_alpha
       << ",\"bdf_endpoint_alpha\":"
       << record.predictor_bdf_endpoint_alpha
       << ",\"source_endpoint_alpha\":"
       << record.predictor_source_endpoint_alpha
       << ",\"enthalpy_solve_calls\":"
       << static_cast<unsigned>(record.predictor_enthalpy_solve_calls) << '}'
       << ",\"thermophysical_enthalpy_endpoint\":{\"status_code\":"
       << static_cast<unsigned>(record.predictor_enthalpy_endpoint.status.code)
       << ",\"termination\":\""
       << termination_name(record.predictor_enthalpy_endpoint.termination)
       << "\",\"iterations\":"
       << record.predictor_enthalpy_endpoint.iterations
       << ",\"initial_true_residual\":"
       << record.predictor_enthalpy_endpoint.initial_true_residual
       << ",\"final_true_residual\":"
       << record.predictor_enthalpy_endpoint.final_true_residual
       << ",\"recursive_residual\":"
       << record.predictor_enthalpy_endpoint.recursive_residual
       << ",\"reduction_calls\":"
       << record.predictor_enthalpy_endpoint.reduction_calls
       << ",\"operator_applies\":"
       << record.predictor_enthalpy_endpoint.operator_applies
       << ",\"preconditioner_applies\":"
       << record.predictor_enthalpy_endpoint.preconditioner_applies
       << ",\"norm_breakdown_restarts\":"
       << record.predictor_enthalpy_endpoint.norm_breakdown_restarts << '}'
       << ",\"pressure\":[";
  for (std::size_t corrector = 0U; corrector < record.pressure.size();
       ++corrector) {
    if (corrector != 0U) json << ',';
    const LinearSolveResult& solve = record.pressure[corrector];
    json << "{\"corrector\":" << corrector + 1U
         << ",\"status_code\":"
         << static_cast<unsigned>(solve.status.code)
         << ",\"termination\":\"" << termination_name(solve.termination)
         << "\",\"iterations\":" << solve.iterations
         << ",\"initial_true_residual\":" << solve.initial_true_residual
         << ",\"final_true_residual\":" << solve.final_true_residual
         << ",\"recursive_residual\":" << solve.recursive_residual
         << ",\"reduction_calls\":" << solve.reduction_calls
         << ",\"operator_applies\":" << solve.operator_applies
         << ",\"preconditioner_applies\":"
         << solve.preconditioner_applies
         << ",\"norm_breakdown_restarts\":"
         << solve.norm_breakdown_restarts
         << ",\"convergence_audits\":" << solve.convergence_audits
         << ",\"convergence_rejections\":"
         << solve.convergence_rejections
         << ",\"final_convergence_metric\":"
         << solve.final_convergence_metric
         << ",\"convergence_limit\":" << solve.convergence_limit
         << ",\"recycle_offered_directions\":"
         << solve.recycle_offered_directions
         << ",\"recycle_retained_directions\":"
         << solve.recycle_retained_directions
         << ",\"recycle_operator_applies\":"
         << solve.recycle_operator_applies
         << ",\"recycle_reduction_calls\":"
         << solve.recycle_reduction_calls
         << ",\"recycle_projection_attempted\":"
         << (solve.recycle_projection_attempted ? "true" : "false")
         << ",\"recycle_projection_accepted\":"
         << (solve.recycle_projection_accepted ? "true" : "false")
         << ",\"recycle_projected_true_residual\":"
         << solve.recycle_projected_true_residual
         << ",\"recycle_cycle_corrections\":"
         << solve.recycle_cycle_corrections
         << ",\"recycle_capture_vector_passes\":"
         << solve.recycle_capture_vector_passes
         << ",\"recycle_capture_cycle_attempts\":"
         << solve.recycle_capture_cycle_attempts
         << ",\"recycle_capture_reduction_calls\":"
         << solve.recycle_capture_reduction_calls
         << ",\"recycle_capture_blocking_operations\":"
         << solve.recycle_capture_blocking_operations << '}';
  }
  json << ']' << ",\"pressure_energy_refinement\":[";
  for (std::size_t index = 0U;
       index < record.pressure_energy_refinement_solve_calls; ++index) {
    if (index != 0U) json << ',';
    const RuntimePressureEnergyRefinementSolve& entry =
        record.pressure_energy_refinement[index];
    const LinearSolveResult& solve = entry.solve;
    json << "{\"ordinal\":" << static_cast<unsigned>(entry.ordinal)
         << ",\"target_generation\":" << entry.target_generation
         << ",\"collective_lineage\":" << entry.collective_lineage
         << ",\"status_code\":"
         << static_cast<unsigned>(solve.status.code)
         << ",\"termination\":\"" << termination_name(solve.termination)
         << "\",\"iterations\":" << solve.iterations
         << ",\"initial_true_residual\":" << solve.initial_true_residual
         << ",\"final_true_residual\":" << solve.final_true_residual
         << ",\"recursive_residual\":" << solve.recursive_residual
         << ",\"reduction_calls\":" << solve.reduction_calls
         << ",\"operator_applies\":" << solve.operator_applies
         << ",\"preconditioner_applies\":"
         << solve.preconditioner_applies
         << ",\"norm_breakdown_restarts\":"
         << solve.norm_breakdown_restarts
         << ",\"convergence_audits\":" << solve.convergence_audits
         << ",\"convergence_rejections\":"
         << solve.convergence_rejections
         << ",\"final_convergence_metric\":"
         << solve.final_convergence_metric
         << ",\"convergence_limit\":" << solve.convergence_limit
         << ",\"recycle_offered_directions\":"
         << solve.recycle_offered_directions
         << ",\"recycle_retained_directions\":"
         << solve.recycle_retained_directions
         << ",\"recycle_operator_applies\":"
         << solve.recycle_operator_applies
         << ",\"recycle_reduction_calls\":"
         << solve.recycle_reduction_calls
         << ",\"recycle_projection_attempted\":"
         << (solve.recycle_projection_attempted ? "true" : "false")
         << ",\"recycle_projection_accepted\":"
         << (solve.recycle_projection_accepted ? "true" : "false")
         << ",\"recycle_projected_true_residual\":"
         << solve.recycle_projected_true_residual
         << ",\"recycle_cycle_corrections\":"
         << solve.recycle_cycle_corrections
         << ",\"recycle_capture_vector_passes\":"
         << solve.recycle_capture_vector_passes
         << ",\"recycle_capture_cycle_attempts\":"
         << solve.recycle_capture_cycle_attempts
         << ",\"recycle_capture_reduction_calls\":"
         << solve.recycle_capture_reduction_calls
         << ",\"recycle_capture_blocking_operations\":"
         << solve.recycle_capture_blocking_operations << '}';
  }
  json << ']' << ",\"momentum_predictor\":[";
  for (std::size_t component = 0U;
       component < record.momentum_predictor.size(); ++component) {
    if (component != 0U) json << ',';
    const LinearSolveResult& solve = record.momentum_predictor[component];
    json << "{\"component\":" << component
         << ",\"status_code\":"
         << static_cast<unsigned>(solve.status.code)
         << ",\"termination\":\"" << termination_name(solve.termination)
         << "\",\"iterations\":" << solve.iterations
         << ",\"initial_true_residual\":" << solve.initial_true_residual
         << ",\"final_true_residual\":" << solve.final_true_residual
         << ",\"recursive_residual\":" << solve.recursive_residual
         << ",\"reduction_calls\":" << solve.reduction_calls
         << ",\"operator_applies\":" << solve.operator_applies
         << ",\"preconditioner_applies\":"
         << solve.preconditioner_applies
         << ",\"norm_breakdown_restarts\":"
         << solve.norm_breakdown_restarts << '}';
  }
  json << ']'
       << ",\"startup\":" << (record.startup ? "true" : "false")
       << ",\"retry\":" << (record.retry ? "true" : "false")
       << ",\"restart_recovery\":"
       << (record.restart_recovery ? "true" : "false")
       << ",\"statistics_eligible\":"
       << (record.statistics_eligible ? "true" : "false")
       << ",\"stages\":[";
  for (std::size_t index = 0U; index < record.stages.size; ++index) {
    if (index != 0U) json << ',';
    const StageTimingRecord stage = record.stages.data[index];
    json << "{\"id\":" << stage.stage << ",\"min_ns\":"
         << stage.minimum_nanoseconds << ",\"mean_ns\":"
         << stage.mean_nanoseconds << ",\"max_ns\":"
         << stage.maximum_nanoseconds << '}';
  }
  json << "]}\n";
  return json.str();
}

std::uint64_t hash_text(std::string_view text) noexcept {
  std::uint64_t hash = UINT64_C(1469598103934665603);
  for (const unsigned char value : text) {
    hash ^= value;
    hash *= UINT64_C(1099511628211);
  }
  return hash == 0U ? 1U : hash;
}

}  // namespace

Status EvidenceWriter::append(MPI_Comm communicator,
                              const std::filesystem::path& evidence_file,
                              const IoServicePlan& services,
                              const RuntimeEvidenceRecord& record) noexcept try {
  int rank = 0;
  if (MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS ||
      evidence_file.empty() || evidence_file.parent_path().empty())
    return {StatusCode::invalid_plan, detail::kOutputInput};
  Status status = validate_record(services, record);
  std::string text;
  if (status) text = encode_record(record);
  std::uint64_t local_hash = status ? hash_text(text) : 0U;
  std::uint64_t minimum = local_hash;
  std::uint64_t maximum = local_hash;
  if (MPI_Allreduce(MPI_IN_PLACE, &minimum, 1, MPI_UINT64_T, MPI_MIN,
                    communicator) != MPI_SUCCESS ||
      MPI_Allreduce(MPI_IN_PLACE, &maximum, 1, MPI_UINT64_T, MPI_MAX,
                    communicator) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, detail::kOutputCollective};
  if (status && (minimum == 0U || minimum != maximum))
    status = {StatusCode::invalid_plan, detail::kOutputInput};
  status = detail::output_collective_status(communicator, status);
  if (!status) return status;
  status = detail::output_create_directory(communicator, rank,
                                           evidence_file.parent_path());
  if (!status) return status;
  if (rank == 0) {
    const RuntimeServiceCapacity* capacity =
        detail::output_service(services, RuntimeServiceKind::evidence);
    status = text.size() <= capacity->maximum_staging_bytes_per_rank &&
                     detail::output_write_file(evidence_file, text, true)
                 ? Status{}
                 : Status{StatusCode::io_failure, detail::kOutputFile};
  }
  return detail::output_collective_status(communicator, status);
} catch (const std::bad_alloc&) {
  return {StatusCode::allocation_failure, detail::kOutputCapacity};
} catch (...) {
  return {StatusCode::io_failure, detail::kOutputFile};
}

}  // namespace hundun::v04
