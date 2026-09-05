// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_app.hpp"
#include "hundun/v04_mpi_runtime.hpp"

#include <mpi.h>

#include <cstdint>
#include <charconv>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>
#include <string_view>

namespace {

using hundun::v04::ApplicationService;
using hundun::v04::ApplicationRunOptions;
using hundun::v04::ApplicationRunReport;
using hundun::v04::CaseValidationReport;
using hundun::v04::Status;
using hundun::v04::StatusCode;

const char* failure_phase_name(hundun::v04::ApplicationFailurePhase phase) {
  using Phase = hundun::v04::ApplicationFailurePhase;
  switch (phase) {
    case Phase::none:
      return "setup";
    case Phase::advance:
      return "advance";
    case Phase::visit:
      return "Visit";
    case Phase::screen:
      return "Screen";
    case Phase::monitor:
      return "Monitor";
    case Phase::restart:
      return "Restart";
    case Phase::resources:
      return "resources";
    case Phase::evidence:
      return "Evidence";
  }
  return "unknown";
}

int finish(Status status, int rank) {
  if (!status && rank == 0) {
    std::cerr << hundun::v04::status_message(status)
              << " detail=" << status.detail << '\n';
  }
  return status ? 0 : 1;
}

Status broadcast_root_status(Status status, int rank) {
  std::uint64_t wire = rank == 0
                           ? (static_cast<std::uint64_t>(status.code) << 32U) |
                                 status.detail
                           : 0U;
  if (MPI_Bcast(&wire, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, 10504U};
  return {static_cast<StatusCode>(wire >> 32U),
          static_cast<std::uint32_t>(wire)};
}

void usage(int rank) {
  if (rank == 0) {
    std::cerr << "usage:\n"
              << "  hundun --version\n"
              << "  hundun validate <case-dir> [--dry-plan]\n"
              << "  hundun run <case-dir> --output <run-dir> --steps <N>"
                 " [--restart <restart-dir>] [--output-interval <N>]"
                 " [--restart-interval <N>]\n"
              << "    interval 0 disables Visit/screen/monitor or Restart;"
                 " evidence remains enabled outside the timed step\n"
              << "  hundun init-case --output <case-dir>\n";
  }
}

void print_numerical_failure(
    const hundun::v04::NumericalFailureContext& failure) {
  if (!failure.valid) return;
  std::cerr << std::setprecision(17)
            << "numerical_failure_context_v1"
            << " field=" << static_cast<unsigned>(failure.field)
            << " field_id=" << failure.field_id
            << " status=" << static_cast<unsigned>(failure.failure.code)
            << '/' << failure.failure.detail
            << " stage=" << failure.stage
            << " global_cell=" << failure.global_cell
            << " global_index=" << failure.global_index.x << ','
            << failure.global_index.y << ',' << failure.global_index.z
            << " rank=" << failure.rank
            << " attempted_step=" << failure.attempted_step
            << " generation=" << failure.generation
            << " time=" << failure.time
            << " target_time=" << failure.target_time
            << " dt_before=" << failure.dt_before
            << " method_before="
            << static_cast<unsigned>(failure.method_before.order) << ','
            << failure.method_before.a0 << ',' << failure.method_before.a1
            << ',' << failure.method_before.a2
            << " origin_before="
            << static_cast<unsigned>(failure.origin_before)
            << " retry_proposed=" << (failure.retry_proposed ? 1 : 0)
            << " dt_after=" << failure.dt_after
            << " method_after="
            << static_cast<unsigned>(failure.method_after.order) << ','
            << failure.method_after.a0 << ',' << failure.method_after.a1
            << ',' << failure.method_after.a2
            << " origin_after=" << static_cast<unsigned>(failure.origin_after)
            << " failed_value=" << failure.failed_value
            << " allowed=" << failure.allowed_minimum << ','
            << failure.allowed_maximum
            << " p_abs=" << failure.pressure_absolute
            << " T_before=" << failure.temperature_before
            << " T_estimate=" << failure.temperature_estimate
            << " rho=" << failure.density_before << ','
            << failure.density_previous << ',' << failure.density_predicted
            << " h=" << failure.enthalpy_before << ','
            << failure.enthalpy_previous
            << " cp_before=" << failure.cp_before
            << " Y=";
  for (std::size_t index = 0U; index < failure.mass_fraction_count; ++index) {
    if (index != 0U) std::cerr << ',';
    std::cerr << failure.mass_fractions[index];
  }
  std::cerr << " Y_truncated="
            << (failure.mass_fractions_truncated ? 1 : 0)
            << " revisions=" << failure.enthalpy_accepted_revision << ','
            << failure.enthalpy_previous_revision << ','
            << failure.nonadvective_accepted_revision << ','
            << failure.nonadvective_previous_revision << ','
            << failure.mass_flux_accepted_revision << ','
            << failure.mass_flux_previous_revision << ','
            << failure.temperature_accepted_revision << ','
            << failure.conductivity_revision << ','
            << failure.velocity_gradient_revision << ','
            << failure.effective_viscosity_revision
            << " boundary_stencil="
            << (failure.physical_boundary_stencil ? 1 : 0) << ','
            << (failure.mpi_boundary_stencil ? 1 : 0) << ','
            << (failure.immersed_interface_cell ? 1 : 0)
            << " mass_div=" << failure.mass_divergence_accepted << ','
            << failure.mass_divergence_previous
            << " advection=" << failure.advection_accepted << ','
            << failure.advection_previous
            << " nonadvective=" << failure.nonadvective_accepted << ','
            << failure.nonadvective_previous
            << " predictor_terms=" << failure.conservative_history_value
            << ',' << failure.accepted_advection_delta << ','
            << failure.accepted_nonadvective_delta << ','
            << failure.previous_advection_delta << ','
            << failure.previous_nonadvective_delta << ','
            << failure.reconstructed_value
            << " rate_terms=" << failure.diffusion_accepted << ','
            << failure.pressure_work_accepted << ','
            << failure.viscous_dissipation_accepted << ','
            << failure.explicit_source_accepted << ','
            << failure.implicit_sink_accepted
            << " rate_complete="
            << (failure.rate_breakdown_complete ? 1 : 0)
            << " first_bad="
            << static_cast<unsigned>(failure.first_bad_contributor)
            << " tvd_envelope="
            << (failure.face_envelope_checked ? 1 : 0) << ','
            << (failure.face_envelope_valid ? 1 : 0) << ','
            << failure.maximum_face_envelope_violation << ','
            << failure.selected_face_value << ','
            << failure.selected_donor_minimum << ','
            << failure.selected_donor_maximum
            << " be_counterfactual=" << failure.counterfactual_be_density
            << ',' << failure.counterfactual_be_enthalpy << ','
            << (failure.counterfactual_be_admissible ? 1 : 0) << '\n';
}

void print_predictor_failure(
    const hundun::v04::ThermophysicalPredictorFailure& failure) {
  if (!failure.valid) return;
  std::cerr << std::setprecision(17)
            << "predictor_failure_context_v1"
            << " reason=" << static_cast<unsigned>(failure.reason)
            << " field=" << static_cast<unsigned>(failure.field)
            << " field_index=" << failure.field_index
            << " constraint=" << static_cast<unsigned>(failure.constraint)
            << " rank=" << failure.rank
            << " has_cell=" << (failure.has_cell ? 1 : 0)
            << " global_index=" << failure.global_index.x << ','
            << failure.global_index.y << ',' << failure.global_index.z
            << " has_substep=" << (failure.has_substep ? 1 : 0)
            << " substep=" << failure.substep
            << " has_margins=" << (failure.has_margins ? 1 : 0)
            << " margins=" << failure.low_margin << ','
            << failure.high_margin << " scalar_mask=" << failure.scalar_mask
            << " dt_substep=" << failure.dt_substep
            << " maximum_cfl=" << failure.maximum_cfl
            << " local_cfl=" << failure.local_cfl
            << " density_current=" << failure.density_current
            << " density_next=" << failure.density_next
            << " quantity_current=" << failure.quantity_current
            << " quantity_previous=" << failure.quantity_previous
            << " divergence=" << failure.divergence
            << " nonadvective_rhs=" << failure.nonadvective_rhs
            << " conserved_next=" << failure.conserved_next
            << " observed_value=" << failure.observed_value
            << " allowed=" << failure.allowed_lower << ','
            << failure.allowed_upper << '\n';
}

bool positive_integer(std::string_view text, std::uint64_t& out) noexcept {
  if (text.empty()) return false;
  std::uint64_t candidate = 0U;
  const char* begin = text.data();
  const char* end = begin + text.size();
  const auto parsed = std::from_chars(begin, end, candidate);
  if (parsed.ec != std::errc{} || parsed.ptr != end || candidate == 0U)
    return false;
  out = candidate;
  return true;
}

bool nonnegative_integer(std::string_view text, std::uint64_t& out) noexcept {
  if (text.empty()) return false;
  std::uint64_t candidate = 0U;
  const char* begin = text.data();
  const char* end = begin + text.size();
  const auto parsed = std::from_chars(begin, end, candidate);
  if (parsed.ec != std::errc{} || parsed.ptr != end) return false;
  out = candidate;
  return true;
}

}  // namespace

int main(int argc, char* argv[]) {
  if (!hundun::v04::prepare_mpi_runtime_environment()) return 2;
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 2;
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int result = 2;
  if (argc == 2 && std::string_view{argv[1]} == "--version") {
    if (rank == 0) std::cout << "HUNDUN-FLOW 1.0.0 source=v0.4\n";
    result = 0;
  } else if ((argc == 3 || argc == 4) &&
             std::string_view{argv[1]} == "validate" &&
             (argc == 3 || std::string_view{argv[3]} == "--dry-plan")) {
    CaseValidationReport report;
    const Status status =
        ApplicationService::validate(MPI_COMM_WORLD, argv[2], report);
    if (status && rank == 0) {
      std::cout << "VALID case=" << report.case_model
                << " product=" << report.product
                << " global=" << report.summary.global_cells.x << 'x'
                << report.summary.global_cells.y << 'x'
                << report.summary.global_cells.z
                << " local=" << report.summary.local_cells.x << 'x'
                << report.summary.local_cells.y << 'x'
                << report.summary.local_cells.z
                << " fields=" << report.summary.field_count
                << " arena_doubles=" << report.summary.arena_doubles
                << " stages=" << report.summary.graph_stage_count
                << " correctors="
                << static_cast<unsigned>(report.summary.pressure_correctors)
                << " pressure_rtol="
                << report.summary.pressure_relative_tolerance
                << " pressure_atol="
                << report.summary.pressure_absolute_tolerance
                << " pressure_maxit="
                << report.summary.pressure_maximum_iterations
                << " pressure_restart="
                << report.summary.pressure_krylov_restart
                << " sealed=" << (report.summary.sealed ? 1 : 0) << '\n';
    }
    result = finish(status, rank);
  } else if (argc == 4 && std::string_view{argv[1]} == "init-case" &&
             std::string_view{argv[2]} == "--output") {
    Status status;
    if (rank == 0) {
      const std::filesystem::path output = argv[3];
      status = ApplicationService::validate_run_directories(
          output.parent_path() / "__hundun_input_guard__", output,
          HUNDUN_V04_SOURCE_ROOT);
      if (status)
        status = ApplicationService::initialize_case_directory(output);
    }
    status = broadcast_root_status(status, rank);
    if (status && rank == 0) std::cout << "INITIALIZED " << argv[3] << '\n';
    result = finish(status, rank);
  } else if (argc >= 7 && std::string_view{argv[1]} == "run") {
    ApplicationRunOptions options;
    options.case_root = argv[2];
    options.source_root = HUNDUN_V04_SOURCE_ROOT;
    bool parsed = true;
    bool saw_output = false;
    bool saw_steps = false;
    for (int index = 3; index < argc && parsed; index += 2) {
      if (index + 1 >= argc) {
        parsed = false;
        break;
      }
      const std::string_view flag{argv[index]};
      const std::string_view value{argv[index + 1]};
      if (flag == "--output" && !saw_output) {
        options.run_directory = std::string{value};
        saw_output = true;
      } else if (flag == "--steps" && !saw_steps) {
        parsed = positive_integer(value, options.steps);
        saw_steps = true;
      } else if (flag == "--output-interval") {
        parsed = nonnegative_integer(value, options.output_interval);
      } else if (flag == "--restart-interval") {
        parsed = nonnegative_integer(value, options.restart_interval);
      } else if (flag == "--restart" &&
                 options.restart_directory.empty()) {
        options.restart_directory = std::string{value};
      } else {
        parsed = false;
      }
    }
    if (!parsed || !saw_output || !saw_steps) {
      usage(rank);
      result = 2;
    } else {
      ApplicationRunReport report;
      const Status status =
          ApplicationService::run(MPI_COMM_WORLD, options, report);
      if (status && rank == 0) {
        std::cout << "COMPLETED steps=" << report.accepted_steps
                  << " time=" << report.final_time
                  << " product=" << report.product
                  << " predictor_limiter_activations="
                  << report.predictor_limiter_activations
                  << " predictor_min_theta="
                  << report.minimum_predictor_theta
                  << " predictor_low_order_transport_passes="
                  << report.predictor_low_order_transport_passes
                  << " predictor_low_order_halo_exchanges="
                  << report.predictor_low_order_halo_exchanges
                  << " advective_convective_cfl_out_max="
                  << report.maximum_advective_convective_cfl_out
                  << " advective_convective_cfl_abs_max="
                  << report.maximum_advective_convective_cfl_abs
                  << " advective_convective_cfl_limit="
                  << report.advective_convective_cfl_limit
                  << " committed_convective_cfl_out_max="
                  << report.maximum_committed_convective_cfl_out
                  << " committed_convective_cfl_abs_max="
                  << report.maximum_committed_convective_cfl_abs
                  << " committed_convective_cfl_limit="
                  << report.committed_convective_cfl_limit
                  << " momentum_afc_applicability="
                  << (report.momentum_predictor_limiter
                              .correction_metrics_applicable
                          ? "applicable"
                          : "not_applicable")
                  << " momentum_afc_active_faces="
                  << report.momentum_predictor_limiter
                         .active_correction_faces
                  << " momentum_afc_limited_faces="
                  << report.momentum_predictor_limiter.activations;
        if (report.momentum_predictor_limiter
                .correction_metrics_applicable) {
          std::cout << " momentum_afc_retained_l1="
                    << report.momentum_predictor_limiter.theta
                    << " momentum_afc_min_alpha="
                    << report.momentum_predictor_limiter.minimum_face_alpha
                    << " momentum_afc_limited_fraction="
                    << report.momentum_predictor_limiter
                           .limited_face_fraction;
        }
        std::cout << '\n';
      }
      if (!status && rank == 0) {
        std::cerr << "termination_phase="
                  << failure_phase_name(report.failure_phase)
                  << " committed_step=" << report.accepted_steps
                  << " committed_time=" << std::setprecision(17)
                  << report.final_time
                  << " output_root=" << options.run_directory << '\n';
      }
      if (!status && rank == 0 && report.failed_stage != 0U) {
        std::cerr << "failed_stage=" << report.failed_stage
                  << " attempts=" << report.attempts;
        const bool terminal_physical_audit_present =
            report.piso.final_flux_revision != 0U;
        if (terminal_physical_audit_present) {
          std::cerr << " terminal_audit=available"
                    << " eos=" << report.piso.eos_residual
                    << " continuity=" << report.piso.continuity_residual;
          if (report.piso.continuity_witness.valid)
            std::cerr << " continuity_witness="
                      << report.piso.continuity_witness.global_index.x << ','
                      << report.piso.continuity_witness.global_index.y << ','
                      << report.piso.continuity_witness.global_index.z << '/'
                      << report.piso.continuity_witness.rank << '/'
                      << report.piso.continuity_witness.raw_balance << '/'
                      << report.piso.continuity_witness.unsteady << '/'
                      << report.piso.continuity_witness.flux_divergence << '/'
                      << report.piso.continuity_witness.scale << '/'
                      << report.piso.continuity_witness.global_domain_edge
                      << '/' << report.piso.continuity_witness.density << '/'
                      << report.piso.continuity_witness.accepted_density << '/'
                      << report.piso.continuity_witness.face_fluxes[0U] << ','
                      << report.piso.continuity_witness.face_fluxes[1U] << ','
                      << report.piso.continuity_witness.face_fluxes[2U] << ','
                      << report.piso.continuity_witness.face_fluxes[3U] << ','
                      << report.piso.continuity_witness.face_fluxes[4U] << ','
                      << report.piso.continuity_witness.face_fluxes[5U];
          std::cerr << " energy=" << report.piso.energy_residual
                    << " closed_mass=" << report.piso.closed_mass_residual
                    << " gauge=" << report.piso.gauge_residual
                    << " committed_cfl_out="
                    << report.piso.committed_convective_cfl_out_max
                    << " committed_cfl_abs="
                    << report.piso.committed_convective_cfl_abs_max
                    << " committed_cfl_limit="
                    << report.piso.committed_convective_cfl_limit;
        } else {
          std::cerr << " terminal_audit=unavailable";
        }
        const hundun::v04::CommittedConvectiveCflCertificate& committed_cfl =
            report.piso.committed_convective_cfl;
        if (committed_cfl.failure_witness.valid) {
          const auto emit_committed_winner = [&](
              const char* label,
              const hundun::v04::ConvectiveCflFailureWitness& winner) {
            std::cerr << ' ' << label << '=' << winner.global_cell.x << ','
                      << winner.global_cell.y << ',' << winner.global_cell.z
                      << '/' << winner.rank << '/' << winner.out << '/'
                      << winner.absolute << '/' << winner.density_volume
                      << '/' << winner.outgoing_mass_flow << '/'
                      << winner.absolute_mass_flow;
          };
          emit_committed_winner(
              "committed_cfl_out_winner",
              committed_cfl.failure_witness.out_winner);
          emit_committed_winner(
              "committed_cfl_abs_winner",
              committed_cfl.failure_witness.absolute_winner);
        }
        if (report.momentum_predictor_limiter.advective_cfl.valid()) {
          const hundun::v04::MomentumAdvectiveCflCertificate& cfl =
              report.momentum_predictor_limiter.advective_cfl;
          std::cerr << " advective_cfl_out=" << cfl.out_max
                    << " advective_cfl_abs=" << cfl.absolute_max
                    << " advective_cfl_limit=" << cfl.limit
                    << " advective_cfl_flux_revision=" << cfl.face_flux;
          if (cfl.failure_witness.valid)
            std::cerr << " advective_cfl_witness="
                      << cfl.failure_witness.global_cell.x << ','
                      << cfl.failure_witness.global_cell.y << ','
                      << cfl.failure_witness.global_cell.z << '/'
                      << cfl.failure_witness.rank << '/'
                      << cfl.failure_witness.out << '/'
                      << cfl.failure_witness.absolute << '/'
                      << cfl.failure_witness.density_volume << '/'
                      << cfl.failure_witness.outgoing_mass_flow << '/'
                      << cfl.failure_witness.absolute_mass_flow;
        }
        std::cerr << " pressure_calls="
                  << static_cast<unsigned>(report.piso.pressure_solve_calls)
                  << " pe_refinement_calls="
                  << static_cast<unsigned>(
                         report.piso.pressure_energy_refinement_solve_calls)
                  << " pe_refinement_termination="
                  << static_cast<unsigned>(
                         report.piso.pressure_energy_refinement_termination)
                  << " p1_iterations=" << report.piso.pressure[0U].iterations
                  << " p1_true="
                  << report.piso.pressure[0U].final_true_residual
                  << " p2_iterations=" << report.piso.pressure[1U].iterations
                  << " p2_termination="
                  << static_cast<unsigned>(
                         report.piso.pressure[1U].termination)
                  << " p2_initial="
                  << report.piso.pressure[1U].initial_true_residual
                  << " p2_true="
                  << report.piso.pressure[1U].final_true_residual
                  << " p2_recursive="
                  << report.piso.pressure[1U].recursive_residual
                  << " p2_breakdowns="
                  << report.piso.pressure[1U].norm_breakdown_restarts
                  << " p2_audits="
                  << report.piso.pressure[1U].convergence_audits
                  << " p2_continuity="
                  << report.piso.pressure[1U].final_convergence_metric
                  << " p2_continuity_limit="
                  << report.piso.pressure[1U].convergence_limit
                  << " p2_continuity_rejections="
                  << report.piso.pressure[1U].convergence_rejections
                  << " p2_application_alpha="
                  << report.piso.pressure[1U].convergence_application_scale
                  << " p2_unscaled_continuity="
                  << report.piso.pressure[1U].convergence_unscaled_metric
                  << " p2_maximum_depletion="
                  << report.piso.pressure[1U].convergence_maximum_depletion
                  << " p2_operator_parity_error="
                  << report.piso.pressure[1U]
                         .convergence_operator_parity_error
                  << " p2_failure_provenance_valid="
                  << report.piso.pressure[1U]
                         .convergence_failure_provenance.valid
                  << " p2_predecessor_c1_application_alpha="
                  << report.piso.pressure[1U]
                         .convergence_failure_provenance
                         .predecessor_application_scale
                  << " p2_predecessor_c1_maximum_depletion="
                  << report.piso.pressure[1U]
                         .convergence_failure_provenance
                         .predecessor_maximum_depletion
                  << " p2_limiting_global_cell="
                  << report.piso.pressure[1U]
                         .convergence_failure_provenance.global_cell
                  << " p2_limiting_global_x="
                  << report.piso.pressure[1U]
                         .convergence_failure_provenance.global_index.x
                  << " p2_limiting_global_y="
                  << report.piso.pressure[1U]
                         .convergence_failure_provenance.global_index.y
                  << " p2_limiting_global_z="
                  << report.piso.pressure[1U]
                         .convergence_failure_provenance.global_index.z
                  << " p2_limiting_owner_rank="
                  << report.piso.pressure[1U]
                         .convergence_failure_provenance.owner_rank
                  << " p2_limiting_pressure_activity="
                  << static_cast<unsigned>(
                         report.piso.pressure[1U]
                             .convergence_failure_provenance
                             .pressure_activity)
                  << " p2_limiting_rho_c1="
                  << report.piso.pressure[1U]
                         .convergence_failure_provenance.density
                  << " p2_limiting_psi="
                  << report.piso.pressure[1U]
                         .convergence_failure_provenance.psi
                  << " p2_limiting_raw_delta_p="
                  << report.piso.pressure[1U]
                         .convergence_failure_provenance.raw_correction
                  << " p2_limiting_depletion="
                  << report.piso.pressure[1U]
                         .convergence_failure_provenance.depletion
                  << " p2_limiting_bdf_storage="
                  << report.piso.pressure[1U]
                         .convergence_failure_provenance.bdf_storage
                  << " p2_limiting_c1_flux_divergence="
                  << report.piso.pressure[1U]
                         .convergence_failure_provenance.flux_divergence
                  << " p2_limiting_reconstructed_rhs="
                  << report.piso.pressure[1U]
                         .convergence_failure_provenance.reconstructed_rhs
                  << " p2_limiting_sealed_rhs="
                  << report.piso.pressure[1U]
                         .convergence_failure_provenance.sealed_rhs
                  << " p2_limiting_rhs_absolute_mismatch="
                  << report.piso.pressure[1U]
                         .convergence_failure_provenance
                         .rhs_absolute_mismatch
                  << " p2_limiting_rhs_relative_mismatch="
                  << report.piso.pressure[1U]
                         .convergence_failure_provenance
                         .rhs_relative_mismatch
                  << " momentum_calls="
                  << static_cast<unsigned>(
                         report.momentum_predictor_solve.solve_calls)
                  << " u1_iterations="
                  << report.momentum_predictor_solve.components[0U].iterations
                  << " u1_true="
                  << report.momentum_predictor_solve.components[0U]
                         .final_true_residual
                  << " u2_iterations="
                  << report.momentum_predictor_solve.components[1U].iterations
                  << " u2_true="
                  << report.momentum_predictor_solve.components[1U]
                         .final_true_residual
                  << " u3_iterations="
                  << report.momentum_predictor_solve.components[2U].iterations
                  << " u3_true="
                  << report.momentum_predictor_solve.components[2U]
                         .final_true_residual
                  << " requested_bdf_order="
                  << static_cast<unsigned>(report.requested_bdf.order)
                  << " effective_bdf_order="
                  << static_cast<unsigned>(report.effective_bdf.order)
                  << " temporal_method_fallback="
                  << (report.temporal_method_fallback ? 1 : 0)
                  << " predictor_calls="
                  << static_cast<unsigned>(
                         report.thermophysical_predictor_calls)
                  << " predictor_theta="
                  << report.thermophysical_predictor.theta
                  << " predictor_constraint="
                  << static_cast<unsigned>(
                         report.thermophysical_predictor.constraint)
                  << " predictor_bdf_alpha="
                  << report.thermophysical_predictor.bdf_endpoint_alpha
                  << " predictor_source_alpha="
                  << report.thermophysical_predictor.source_endpoint_alpha
                  << " predictor_high_margin="
                  << report.thermophysical_predictor.high_margin
                  << " predictor_rank="
                  << report.thermophysical_predictor.limiting_rank
                  << " predictor_cell="
                  << report.thermophysical_predictor.limiting_cell.x << ','
                  << report.thermophysical_predictor.limiting_cell.y << ','
                  << report.thermophysical_predictor.limiting_cell.z;
      }
      if (!status && rank == 0 && report.failed_stage != 0U &&
          report.pressure_energy_globalization.valid) {
        const hundun::v04::PressureEnergyGlobalizationAttemptReport&
            globalization =
            report.pressure_energy_globalization;
        std::cerr << " pe_corrector="
                  << static_cast<unsigned>(globalization.corrector)
                  << " pe_trajectory_count="
                  << static_cast<unsigned>(globalization.trajectory_count)
                  << " pe_max_dp="
                  << globalization.maximum_absolute_pressure_correction
                  << " pe_max_dh="
                  << globalization.maximum_absolute_enthalpy_correction
                  << " pe_baseline="
                  << globalization.baseline.global_normalized_continuity
                  << '/'
                  << globalization.baseline.global_normalized_energy;
        for (std::size_t index = 0U;
             index < globalization.sample_count; ++index) {
          const hundun::v04::PressureEnergyGlobalizationSample& sample =
              globalization.candidates[index];
          std::cerr << " pe" << index << '=' << sample.alpha << '/'
                    << sample.global_normalized_continuity << '/'
                    << sample.global_normalized_energy << '/'
                    << sample.thermodynamically_admissible << '/'
                    << sample.state_and_flux_finite;
        }
        for (std::size_t index = 0U;
             index < globalization.trajectory_count &&
             index < globalization.trajectory.size(); ++index) {
          const auto& iteration = globalization.trajectory[index];
          std::cerr
              << " pe_path" << index << '='
              << static_cast<unsigned>(iteration.corrector) << '/'
              << static_cast<unsigned>(iteration.refinement_iteration) << '/'
              << iteration.baseline.global_normalized_continuity << '/'
              << iteration.baseline.global_normalized_energy << '/'
              << iteration.selected.alpha << '/'
              << iteration.selected.global_normalized_continuity << '/'
              << iteration.selected.global_normalized_energy << '/'
              << iteration.maximum_absolute_pressure_correction << '/'
              << iteration.maximum_absolute_enthalpy_correction;
        }
      }
      if (!status && rank == 0 && report.failed_stage != 0U)
        std::cerr << '\n';
      if (!status && rank == 0)
        print_numerical_failure(report.numerical_failure);
      if (!status && rank == 0)
        print_predictor_failure(report.thermophysical_predictor.failure);
      result = finish(status, rank);
    }
  } else {
    usage(rank);
  }
  MPI_Finalize();
  return result;
}
