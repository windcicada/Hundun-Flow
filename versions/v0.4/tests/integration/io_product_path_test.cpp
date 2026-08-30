// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_io.hpp"

#include "../support/piso_fixture.hpp"

#include <mpi.h>
#include <unistd.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <string>
#include <string_view>

namespace {

namespace fs = std::filesystem;
using namespace hundun::v04;

bool expect(bool condition, std::string_view description) {
  if (!condition) std::cerr << "FAIL: " << description << '\n';
  return condition;
}

std::string read(const fs::path& path) {
  std::ifstream input(path, std::ios::binary);
  return {std::istreambuf_iterator<char>(input),
          std::istreambuf_iterator<char>()};
}

CartesianMeshSpec mesh(bool stretched) {
  CartesianMeshSpec spec;
  spec.kind = stretched ? GeometryKind::tensor_stretched
                        : GeometryKind::uniform;
  spec.lower = {0.0, 0.0, 0.0};
  spec.upper = {1.0, 1.0, 1.0};
  spec.has_exact_cells = true;
  spec.exact_cells = stretched ? Int3{17, 11, 7} : Int3{6, 5, 4};
  spec.has_base_spacing = stretched;
  spec.base_spacing = {1.04 / 17.0, 1.04 / 11.0, 1.04 / 7.0};
  spec.minimum_spacing = stretched
                             ? Real3{0.96 / 17.0, 0.96 / 11.0, 0.96 / 7.0}
                             : Real3{1.0 / 6.0, 1.0 / 5.0, 1.0 / 4.0};
  spec.max_growth_ratio = stretched ? 1.0 + 0.8 / 17.0 : 1.0;
  if (stretched)
    spec.focus_regions.push_back(
        {{0.35, 0.35, 0.35},
         {0.65, 0.65, 0.65},
         {0.98 / 17.0, 0.98 / 11.0, 0.98 / 7.0}});
  spec.limits.max_global_cells = stretched ? 17U * 11U * 7U : 120U;
  spec.limits.max_memory_bytes_per_rank = UINT64_C(134217728);
  return spec;
}

bool run_case(const fs::path& root, bool stretched) {
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  if (!CartesianGeometryCompiler::compile(MPI_COMM_SELF, mesh(stretched), {},
                                           geometry, patch)) {
    std::cerr << "geometry failed stretched=" << stretched << '\n';
    return false;
  }
  test::OwnedField velocity =
      test::make_field(0U, patch.cells, 3U, 0U, 11U, 21U);
  test::OwnedField pressure =
      test::make_field(1U, patch.cells, 1U, 0U, 12U, 22U);
  for (double& value : velocity.storage) value = 1.25;
  for (double& value : pressure.storage) value = 2.5;
  const std::array<SnapshotFieldSpec, 2U> specs{{{0U, 3U}, {1U, 1U}}};
  const std::size_t local_cells =
      static_cast<std::size_t>(patch.cells.x) * patch.cells.y * patch.cells.z;
  std::array<RuntimeServiceCapacity, 5U> capacities{};
  for (std::size_t index = 0U; index < capacities.size(); ++index)
    capacities[index] = {static_cast<RuntimeServiceKind>(index),
                         static_cast<StageId>(200U + index),
                         local_cells * 4U * sizeof(double),
                         UINT64_C(1048576), 8U};
  IoServicePlan services;
  if (!IoServicePlan::compile({specs.data(), specs.size()},
                              {capacities.data(), capacities.size()},
                              local_cells,
                              services))
    return false;
  const std::array<SnapshotFieldView, 2U> fields{{
      {"U", as_const(velocity.view), velocity.view.revision},
      {"pi", as_const(pressure.view), pressure.view.revision}}};
  CommittedOutputSnapshot snapshot{&geometry,
                                   patch,
                                   UINT64_C(0x19020001),
                                   UINT64_C(0x19020002),
                                   0.25,
                                   25U,
                                   {fields.data(), fields.size()},
                                   true};
  const fs::path visit = root / (stretched ? "Visit-stretched" : "Visit-uniform");
  const Status visit_status =
      VisitWriter::write(MPI_COMM_SELF, visit, services, snapshot);
  if (!visit_status)
    std::cerr << "visit failed stretched=" << stretched << " code="
              << static_cast<unsigned>(visit_status.code)
              << " detail=" << visit_status.detail << '\n';
  bool passed = static_cast<bool>(visit_status);
  const std::string extension = stretched ? ".vtr" : ".vti";
  const fs::path data = visit /
      ("step-00000000000000000025-rank-00000000" + extension);
  const std::string bytes = read(data);
  passed &= !bytes.empty() &&
            bytes.find(stretched ? "RectilinearGrid" : "ImageData") !=
                std::string::npos &&
            bytes.find("AppendedData encoding=\"raw\"") != std::string::npos;
  passed &= read(visit / "step-00000000000000000025.visit").find(
                "!NBLOCKS 1") != std::string::npos;
  passed &= static_cast<bool>(ScreenWriter::append(
      MPI_COMM_SELF, root / "screen.log", services, snapshot,
      "continuity=0 eos=0"));
  passed &= static_cast<bool>(MonitorWriter::append(
      MPI_COMM_SELF, root / "monitor" / "flow.jsonl", services, snapshot,
      "{\"mass\":1.0}"));
  passed &= read(root / "screen.log").find("step=25") != std::string::npos;
  passed &= read(root / "monitor" / "flow.jsonl").find("\"payload\"") !=
            std::string::npos;

  const std::array<StageTimingRecord, 2U> timings{{
      {30U, 10U, 20U, 30U}, {50U, 40U, 50U, 60U}}};
  RuntimeEvidenceRecord evidence;
  evidence.build = 1U;
  evidence.binary = 2U;
  evidence.case_model = 3U;
  evidence.product = snapshot.plan;
  evidence.step = snapshot.step;
  evidence.time = snapshot.time;
  evidence.requested_bdf_order = 2U;
  evidence.bdf_order = 2U;
  evidence.thermophysical_predictor_calls = 1U;
  evidence.maximum_rank_step_nanoseconds = 100U;
  evidence.blocking_collectives = 1U;
  evidence.pressure_solve_calls = 2U;
  evidence.pressure_solve_contract =
      RuntimePressureSolveContract::pressure_continuity;
  evidence.pressure_energy_refinement_termination =
      RuntimePressureEnergyRefinementTermination::
          component_residuals_converged;
  evidence.terminal_physical_audit.present = true;
  evidence.terminal_physical_audit.final_flux_revision = 17U;
  evidence.terminal_physical_audit.eos_tolerance = 1.0e-6;
  evidence.terminal_physical_audit.continuity_tolerance = 1.0e-6;
  evidence.terminal_physical_audit.energy_tolerance = 1.0e-6;
  evidence.terminal_physical_audit.closed_mass_tolerance = 1.0e-6;
  evidence.terminal_physical_audit.gauge_tolerance = 1.0e-6;
  evidence.terminal_physical_audit.committed_convective_cfl_out_max = 0.25;
  evidence.terminal_physical_audit.committed_convective_cfl_abs_max = 0.125;
  evidence.terminal_physical_audit.committed_convective_cfl_limit = 0.5;
  evidence.momentum_advective_cfl = {
      true, 41U, 14U, 15U, 16U, 0U, 0.25, 0.125, 0.25, 0.5};
  evidence.momentum_predictor_solve_calls = 3U;
  evidence.predictor_blocking_collectives = 1U;
  evidence.pressure[0U].termination = LinearTermination::converged;
  evidence.pressure[1U].termination = LinearTermination::converged;
  evidence.pressure[1U].convergence_audits = 1U;
  evidence.pressure[1U].convergence_limit = 1.0e-6;
  for (LinearSolveResult& solve : evidence.momentum_predictor) {
    solve.termination = LinearTermination::converged;
  }
  evidence.stages = {timings.data(), timings.size()};
  evidence.startup = true;
  evidence.statistics_eligible = false;
  passed &= static_cast<bool>(EvidenceWriter::append(
      MPI_COMM_SELF, root / "evidence" / "runtime.jsonl", services,
      evidence));
  const std::string evidence_text =
      read(root / "evidence" / "runtime.jsonl");
  passed &= evidence_text.find("HUNDUN_V04_EVIDENCE_V5") !=
                std::string::npos &&
            evidence_text.find("\"statistics_eligible\":false") !=
                std::string::npos &&
            evidence_text.find("\"pressure_solve_calls\":2") !=
                std::string::npos &&
            evidence_text.find("\"momentum_predictor_solve_calls\":3") !=
                std::string::npos &&
            evidence_text.find("\"requested_bdf_order\":2") !=
                std::string::npos &&
            evidence_text.find("\"bdf_order\":2") != std::string::npos &&
            evidence_text.find("\"temporal_method_fallback\":false") !=
                std::string::npos &&
            evidence_text.find("\"thermophysical_predictor_calls\":1") !=
                std::string::npos &&
            evidence_text.find("\"momentum_predictor\":[") !=
                std::string::npos &&
            evidence_text.find(
                "\"momentum_predictor_limiter\":{\"scheme\":\"common_face_afc_v2\",\"limited\":false,\"retained_correction_l1_ratio\":1,\"limited_faces\":0,\"advective_cfl\":{\"present\":true,\"plan\":41,\"time_revision\":14,\"density_revision\":15,\"face_flux_revision\":16,\"activity_collective\":0,\"dt\":0.25,\"out_max\":0.125,\"abs_max\":0.25,\"limit\":0.5}}") !=
                std::string::npos &&
            evidence_text.find("\"low_state\":\"none\"") !=
                std::string::npos &&
            evidence_text.find("\"enthalpy_solve_calls\":0") !=
                std::string::npos &&
            evidence_text.find("\"thermophysical_enthalpy_endpoint\":{") !=
                std::string::npos &&
            evidence_text.find("\"mass_flux_scale\":1") !=
                std::string::npos &&
            evidence_text.find("\"bdf_endpoint_alpha\":1") !=
                std::string::npos &&
            evidence_text.find("\"source_endpoint_alpha\":1") !=
                std::string::npos &&
            evidence_text.find(
                "\"pressure_solve_contract\":\"pressure_continuity\"") !=
                std::string::npos &&
            evidence_text.find(
                "\"pressure_energy_refinement_solve_calls\":0") !=
                std::string::npos &&
            evidence_text.find(
                "\"pressure_energy_refinement_termination\":\"component_residuals_converged\"") !=
                std::string::npos &&
            evidence_text.find("\"pressure_energy_refinement\":[]") !=
                std::string::npos &&
            evidence_text.find("\"terminal_physical_audit\":{") !=
                std::string::npos &&
            evidence_text.find("\"final_flux_revision\":17") !=
                std::string::npos &&
            evidence_text.find(
                "\"committed_convective_cfl\":{\"out_max\":0.25,\"abs_max\":0.125,\"limit\":0.5}") !=
                std::string::npos &&
            evidence_text.find("\"convergence_audits\":1") !=
                std::string::npos;

  // A zero right-hand side terminates before the optional terminal
  // convergence audit is invoked.  This is the legitimate C2 result emitted
  // by the continuity-energy coupled pressure solve during rank-change
  // restart recovery.
  RuntimeEvidenceRecord coupled_c2_zero_rhs = evidence;
  coupled_c2_zero_rhs.pressure_solve_contract =
      RuntimePressureSolveContract::continuity_energy_coupled;
  coupled_c2_zero_rhs.pressure[1U].termination =
      LinearTermination::zero_rhs;
  coupled_c2_zero_rhs.pressure[1U].convergence_audits = 0U;
  coupled_c2_zero_rhs.pressure[1U].convergence_rejections = 0U;
  coupled_c2_zero_rhs.pressure[1U].final_convergence_metric = 0.0;
  coupled_c2_zero_rhs.pressure[1U].convergence_limit = 0.0;
  const Status coupled_c2_zero_rhs_status = EvidenceWriter::append(
      MPI_COMM_SELF, root / "evidence" / "coupled-c2-zero-rhs.jsonl",
      services, coupled_c2_zero_rhs);
  passed &= expect(
      static_cast<bool>(coupled_c2_zero_rhs_status),
      "C2 zero-RHS evidence does not require a convergence audit");
  if (!coupled_c2_zero_rhs_status)
    std::cerr << "C2 zero-RHS evidence code="
              << static_cast<unsigned>(coupled_c2_zero_rhs_status.code)
              << " detail=" << coupled_c2_zero_rhs_status.detail << '\n';
  const std::string coupled_c2_zero_rhs_text =
      read(root / "evidence" / "coupled-c2-zero-rhs.jsonl");
  passed &= coupled_c2_zero_rhs_text.find(
                "\"pressure_solve_contract\":\"continuity_energy_coupled\"") !=
                std::string::npos &&
            coupled_c2_zero_rhs_text.find("\"termination\":\"zero_rhs\"") !=
                std::string::npos &&
            coupled_c2_zero_rhs_text.find("\"convergence_audits\":0") !=
                std::string::npos;

  RuntimeEvidenceRecord one_refinement = coupled_c2_zero_rhs;
  one_refinement.linear_iterations = 2U;
  one_refinement.pressure_energy_refinement_solve_calls = 1U;
  one_refinement.pressure_energy_refinement[0U].solve.termination =
      LinearTermination::converged;
  one_refinement.pressure_energy_refinement[0U].solve.iterations = 2U;
  one_refinement.pressure_energy_refinement[0U].solve.initial_true_residual =
      1.0;
  one_refinement.pressure_energy_refinement[0U].solve.final_true_residual =
      1.0e-8;
  one_refinement.pressure_energy_refinement[0U].solve.recursive_residual =
      1.0e-8;
  one_refinement.pressure_energy_refinement[0U].target_generation = 71U;
  one_refinement.pressure_energy_refinement[0U].collective_lineage = 81U;
  one_refinement.pressure_energy_refinement[0U].ordinal = 1U;
  const fs::path one_refinement_path =
      root / "evidence" / "one-refinement.jsonl";
  passed &= expect(
      static_cast<bool>(EvidenceWriter::append(
          MPI_COMM_SELF, one_refinement_path, services, one_refinement)),
      "V4 evidence accepts one typed same-target refinement solve");
  const std::string one_refinement_text = read(one_refinement_path);
  passed &= one_refinement_text.find(
                "\"pressure_energy_refinement_solve_calls\":1") !=
                std::string::npos &&
            one_refinement_text.find(
                "\"pressure_energy_refinement\":[{\"ordinal\":1,\"target_generation\":71,\"collective_lineage\":81") !=
                std::string::npos;

  RuntimeEvidenceRecord two_refinements = one_refinement;
  two_refinements.linear_iterations = 4U;
  two_refinements.pressure_energy_refinement_solve_calls = 2U;
  two_refinements.pressure_energy_refinement[1U] =
      two_refinements.pressure_energy_refinement[0U];
  two_refinements.pressure_energy_refinement[1U].collective_lineage = 82U;
  two_refinements.pressure_energy_refinement[1U].ordinal = 2U;

  RuntimeEvidenceRecord wrong_refinement_order = two_refinements;
  wrong_refinement_order.pressure_energy_refinement[0U].ordinal = 2U;
  wrong_refinement_order.pressure_energy_refinement[1U].ordinal = 1U;
  passed &= expect(
      EvidenceWriter::append(
          MPI_COMM_SELF, root / "evidence" / "wrong-refinement-order.jsonl",
          services, wrong_refinement_order)
              .code == StatusCode::invalid_plan,
      "V4 evidence rejects out-of-order refinement ordinals");
  RuntimeEvidenceRecord wrong_refinement_target = two_refinements;
  wrong_refinement_target.pressure_energy_refinement[1U].target_generation =
      72U;
  passed &= expect(
      EvidenceWriter::append(
          MPI_COMM_SELF, root / "evidence" / "wrong-refinement-target.jsonl",
          services, wrong_refinement_target)
              .code == StatusCode::invalid_plan,
      "V4 evidence rejects refinement solves from different targets");
  RuntimeEvidenceRecord wrong_refinement_lineage = two_refinements;
  wrong_refinement_lineage.pressure_energy_refinement[1U]
      .collective_lineage = 81U;
  passed &= expect(
      EvidenceWriter::append(
          MPI_COMM_SELF,
          root / "evidence" / "wrong-refinement-lineage.jsonl", services,
          wrong_refinement_lineage)
              .code == StatusCode::invalid_plan,
      "V4 evidence rejects duplicate refinement lineages");
  RuntimeEvidenceRecord wrong_refinement_count = one_refinement;
  wrong_refinement_count.pressure_energy_refinement_solve_calls = 2U;
  passed &= expect(
      EvidenceWriter::append(
          MPI_COMM_SELF, root / "evidence" / "wrong-refinement-count.jsonl",
          services, wrong_refinement_count)
              .code == StatusCode::invalid_plan,
      "V4 evidence rejects a count above its valid typed prefix");
  RuntimeEvidenceRecord wrong_refinement_termination = one_refinement;
  wrong_refinement_termination.pressure_energy_refinement_termination =
      RuntimePressureEnergyRefinementTermination::
          iteration_capacity_exhausted;
  passed &= expect(
      EvidenceWriter::append(
          MPI_COMM_SELF,
          root / "evidence" / "wrong-refinement-termination.jsonl", services,
          wrong_refinement_termination)
              .code == StatusCode::invalid_plan,
      "accepted V4 evidence rejects non-converged refinement termination");
  RuntimeEvidenceRecord nonfinite_refinement = one_refinement;
  nonfinite_refinement.pressure_energy_refinement[0U]
      .solve.final_true_residual = std::numeric_limits<double>::quiet_NaN();
  passed &= expect(
      EvidenceWriter::append(
          MPI_COMM_SELF, root / "evidence" / "nonfinite-refinement.jsonl",
          services, nonfinite_refinement)
              .code == StatusCode::invalid_plan,
      "V4 evidence rejects non-finite refinement residuals");
  RuntimeEvidenceRecord missing_refinement_resource = one_refinement;
  missing_refinement_resource.linear_iterations = 1U;
  passed &= expect(
      EvidenceWriter::append(
          MPI_COMM_SELF,
          root / "evidence" / "missing-refinement-resource.jsonl", services,
          missing_refinement_resource)
              .code == StatusCode::invalid_plan,
      "V4 linear-iteration resources include refinement solve work");

  RuntimeEvidenceRecord missing_terminal_audit = coupled_c2_zero_rhs;
  missing_terminal_audit.terminal_physical_audit.present = false;
  passed &= expect(
      EvidenceWriter::append(
          MPI_COMM_SELF,
          root / "evidence" / "missing-terminal-physical-audit.jsonl",
          services, missing_terminal_audit)
              .code == StatusCode::invalid_plan,
      "coupled evidence requires a terminal physical audit");
  RuntimeEvidenceRecord missing_terminal_flux = coupled_c2_zero_rhs;
  missing_terminal_flux.terminal_physical_audit.final_flux_revision = 0U;
  passed &= expect(
      EvidenceWriter::append(
          MPI_COMM_SELF,
          root / "evidence" / "missing-terminal-final-flux.jsonl", services,
          missing_terminal_flux)
              .code == StatusCode::invalid_plan,
      "terminal physical audit binds the final flux revision");
  RuntimeEvidenceRecord nonfinite_terminal_audit = coupled_c2_zero_rhs;
  nonfinite_terminal_audit.terminal_physical_audit.energy_residual =
      std::numeric_limits<double>::quiet_NaN();
  passed &= expect(
      EvidenceWriter::append(
          MPI_COMM_SELF,
          root / "evidence" / "nonfinite-terminal-physical-audit.jsonl",
          services, nonfinite_terminal_audit)
              .code == StatusCode::invalid_plan,
      "terminal physical audit rejects non-finite metrics");
  RuntimeEvidenceRecord over_limit_terminal_audit = coupled_c2_zero_rhs;
  over_limit_terminal_audit.terminal_physical_audit.continuity_residual =
      2.0e-6;
  passed &= expect(
      EvidenceWriter::append(
          MPI_COMM_SELF,
          root / "evidence" / "over-limit-terminal-physical-audit.jsonl",
          services, over_limit_terminal_audit)
              .code == StatusCode::invalid_plan,
      "terminal physical audit rejects a metric above its limit");
  RuntimeEvidenceRecord missing_committed_cfl_limit = coupled_c2_zero_rhs;
  missing_committed_cfl_limit.terminal_physical_audit
      .committed_convective_cfl_limit = 0.0;
  passed &= expect(
      EvidenceWriter::append(
          MPI_COMM_SELF,
          root / "evidence" / "missing-committed-convective-cfl-limit.jsonl",
          services, missing_committed_cfl_limit)
              .code == StatusCode::invalid_plan,
      "terminal physical audit requires a positive committed CFL limit");
  RuntimeEvidenceRecord nonfinite_committed_cfl = coupled_c2_zero_rhs;
  nonfinite_committed_cfl.terminal_physical_audit
      .committed_convective_cfl_abs_max =
      std::numeric_limits<double>::infinity();
  passed &= expect(
      EvidenceWriter::append(
          MPI_COMM_SELF,
          root / "evidence" / "nonfinite-committed-convective-cfl.jsonl",
          services, nonfinite_committed_cfl)
              .code == StatusCode::invalid_plan,
      "terminal physical audit rejects non-finite committed CFL");
  RuntimeEvidenceRecord over_limit_committed_cfl = coupled_c2_zero_rhs;
  over_limit_committed_cfl.terminal_physical_audit
      .committed_convective_cfl_out_max = 0.5000000001;
  passed &= expect(
      EvidenceWriter::append(
          MPI_COMM_SELF,
          root / "evidence" / "over-limit-committed-convective-cfl.jsonl",
          services, over_limit_committed_cfl)
              .code == StatusCode::invalid_plan,
      "terminal physical audit rejects committed outward CFL above its limit");
  constexpr double kCflSlack =
      1.0 + 64.0 * std::numeric_limits<double>::epsilon();
  RuntimeEvidenceRecord cfl_roundoff_edge = coupled_c2_zero_rhs;
  cfl_roundoff_edge.terminal_physical_audit
      .committed_convective_cfl_out_max =
      cfl_roundoff_edge.terminal_physical_audit
          .committed_convective_cfl_limit *
      kCflSlack;
  passed &= expect(
      static_cast<bool>(EvidenceWriter::append(
          MPI_COMM_SELF,
          root / "evidence" / "committed-convective-cfl-roundoff-edge.jsonl",
          services, cfl_roundoff_edge)),
      "terminal physical audit accepts the 64-epsilon CFL comparison edge");
  RuntimeEvidenceRecord cfl_beyond_roundoff = cfl_roundoff_edge;
  cfl_beyond_roundoff.terminal_physical_audit
      .committed_convective_cfl_out_max = std::nextafter(
      cfl_roundoff_edge.terminal_physical_audit
          .committed_convective_cfl_out_max,
      std::numeric_limits<double>::infinity());
  passed &= expect(
      EvidenceWriter::append(
          MPI_COMM_SELF,
          root / "evidence" / "committed-convective-cfl-beyond-roundoff.jsonl",
          services, cfl_beyond_roundoff)
              .code == StatusCode::invalid_plan,
      "terminal physical audit rejects CFL beyond the 64-epsilon edge");
  RuntimeEvidenceRecord missing_advective_cfl = coupled_c2_zero_rhs;
  missing_advective_cfl.momentum_advective_cfl.present = false;
  passed &= expect(
      EvidenceWriter::append(
          MPI_COMM_SELF,
          root / "evidence" / "missing-advective-convective-cfl.jsonl",
          services, missing_advective_cfl)
              .code == StatusCode::invalid_plan,
      "momentum limiter requires its advective CFL certificate");
  RuntimeEvidenceRecord invalid_advective_cfl = coupled_c2_zero_rhs;
  invalid_advective_cfl.momentum_advective_cfl.out_max = 0.5000000001;
  passed &= expect(
      EvidenceWriter::append(
          MPI_COMM_SELF,
          root / "evidence" / "invalid-advective-convective-cfl.jsonl",
          services, invalid_advective_cfl)
              .code == StatusCode::invalid_plan,
      "momentum limiter rejects an advective CFL above its limit");
  RuntimeEvidenceRecord same_revision_advective_cfl = coupled_c2_zero_rhs;
  same_revision_advective_cfl.momentum_advective_cfl.face_flux_revision =
      same_revision_advective_cfl.terminal_physical_audit.final_flux_revision;
  passed &= expect(
      EvidenceWriter::append(
          MPI_COMM_SELF,
          root / "evidence" / "same-revision-advective-convective-cfl.jsonl",
          services, same_revision_advective_cfl)
              .code == StatusCode::invalid_plan,
      "advective and committed CFL certificates require distinct flux revisions");
  RuntimeEvidenceRecord mismatched_limit_advective_cfl = coupled_c2_zero_rhs;
  mismatched_limit_advective_cfl.momentum_advective_cfl.limit = 0.75;
  passed &= expect(
      EvidenceWriter::append(
          MPI_COMM_SELF,
          root / "evidence" / "mismatched-advective-convective-cfl-limit.jsonl",
          services, mismatched_limit_advective_cfl)
              .code == StatusCode::invalid_plan,
      "advective and committed CFL certificates require one configured limit");
  RuntimeEvidenceRecord missing_ibm_activity = coupled_c2_zero_rhs;
  missing_ibm_activity.stl = 18U;
  passed &= expect(
      EvidenceWriter::append(
          MPI_COMM_SELF, root / "evidence" / "missing-ibm-activity.jsonl",
          services, missing_ibm_activity)
              .code == StatusCode::invalid_plan,
      "an IBM evidence row requires a nonzero activity collective");
  RuntimeEvidenceRecord fabricated_ibm_activity = coupled_c2_zero_rhs;
  fabricated_ibm_activity.momentum_advective_cfl.activity_collective = 18U;
  passed &= expect(
      EvidenceWriter::append(
          MPI_COMM_SELF, root / "evidence" / "fabricated-ibm-activity.jsonl",
          services, fabricated_ibm_activity)
              .code == StatusCode::invalid_plan,
      "a no-IBM evidence row rejects a fabricated activity collective");
  RuntimeEvidenceRecord wrong_pressure_contract = coupled_c2_zero_rhs;
  wrong_pressure_contract.pressure_solve_contract =
      RuntimePressureSolveContract::pressure_continuity;
  passed &= expect(
      EvidenceWriter::append(
          MPI_COMM_SELF, root / "evidence" / "wrong-pressure-contract.jsonl",
          services, wrong_pressure_contract)
              .code == StatusCode::invalid_plan,
      "pressure-continuity C2 still requires its linear convergence audit");
  RuntimeEvidenceRecord fabricated_coupled_audit = coupled_c2_zero_rhs;
  fabricated_coupled_audit.pressure[1U].convergence_audits = 1U;
  fabricated_coupled_audit.pressure[1U].convergence_limit = 1.0e-6;
  passed &= expect(
      EvidenceWriter::append(
          MPI_COMM_SELF, root / "evidence" / "fabricated-coupled-audit.jsonl",
          services, fabricated_coupled_audit)
              .code == StatusCode::invalid_plan,
      "coupled pressure evidence cannot fabricate a linear audit");

  RuntimeEvidenceRecord legacy_energy_gate_disabled = evidence;
  legacy_energy_gate_disabled.terminal_physical_audit.energy_residual = 2.0;
  legacy_energy_gate_disabled.terminal_physical_audit.energy_tolerance = 0.0;
  passed &= expect(
      static_cast<bool>(EvidenceWriter::append(
          MPI_COMM_SELF,
          root / "evidence" / "legacy-energy-gate-disabled.jsonl", services,
          legacy_energy_gate_disabled)),
      "pressure-continuity evidence preserves the disabled legacy energy gate");
  RuntimeEvidenceRecord coupled_energy_gate_disabled = coupled_c2_zero_rhs;
  coupled_energy_gate_disabled.terminal_physical_audit.energy_tolerance = 0.0;
  passed &= expect(
      EvidenceWriter::append(
          MPI_COMM_SELF,
          root / "evidence" / "coupled-energy-gate-disabled.jsonl", services,
          coupled_energy_gate_disabled)
              .code == StatusCode::invalid_plan,
      "continuity-energy coupled evidence requires an enabled energy gate");

  RuntimeEvidenceRecord temporal_fallback = evidence;
  temporal_fallback.bdf_order = 1U;
  temporal_fallback.thermophysical_predictor_calls = 2U;
  temporal_fallback.temporal_method_fallback = true;
  temporal_fallback.blocking_collectives = 2U;
  temporal_fallback.predictor_blocking_collectives = 2U;
  const fs::path temporal_fallback_path =
      root / "evidence" / "runtime-temporal-fallback.jsonl";
  passed &= static_cast<bool>(EvidenceWriter::append(
      MPI_COMM_SELF, temporal_fallback_path, services, temporal_fallback));
  const std::string temporal_fallback_text = read(temporal_fallback_path);
  passed &= temporal_fallback_text.find("\"requested_bdf_order\":2") !=
                std::string::npos &&
            temporal_fallback_text.find("\"bdf_order\":1") !=
                std::string::npos &&
            temporal_fallback_text.find(
                "\"temporal_method_fallback\":true") !=
                std::string::npos &&
            temporal_fallback_text.find(
                "\"thermophysical_predictor_calls\":2") !=
                std::string::npos;

  RuntimeEvidenceRecord scaled_euler = evidence;
  scaled_euler.blocking_collectives = 3U;
  scaled_euler.predictor_theta = 0.5;
  scaled_euler.predictor_mass_flux_scale = 0.25;
  scaled_euler.predictor_low_margin = 1.0;
  scaled_euler.predictor_high_margin = -1.0;
  scaled_euler.predictor_low_order_transport_passes = 1U;
  scaled_euler.predictor_low_order_substeps = 1U;
  scaled_euler.predictor_blocking_collectives = 3U;
  scaled_euler.predictor_limiting_cell_x = 1;
  scaled_euler.predictor_limiting_cell_y = 2;
  scaled_euler.predictor_limiting_cell_z = 3;
  scaled_euler.predictor_limiting_rank = 0;
  scaled_euler.predictor_constraint = 1U;
  scaled_euler.predictor_low_state = 2U;
  scaled_euler.predictor_limited = true;
  const fs::path scaled_euler_path =
      root / "evidence" / "runtime-scaled-euler.jsonl";
  passed &= static_cast<bool>(EvidenceWriter::append(
      MPI_COMM_SELF, scaled_euler_path, services, scaled_euler));
  const std::string scaled_euler_text = read(scaled_euler_path);
  passed &= scaled_euler_text.find("\"low_state\":\"scaled_euler\"") !=
                std::string::npos &&
            scaled_euler_text.find("\"mass_flux_scale\":0.25") !=
                std::string::npos;
  RuntimeEvidenceRecord bdf_homotopy = scaled_euler;
  bdf_homotopy.predictor_low_state = 5U;
  bdf_homotopy.predictor_mass_flux_scale = 1.0;
  bdf_homotopy.predictor_bdf_endpoint_alpha = 0.375;
  const fs::path bdf_homotopy_path =
      root / "evidence" / "runtime-bdf-homotopy.jsonl";
  passed &= static_cast<bool>(EvidenceWriter::append(
      MPI_COMM_SELF, bdf_homotopy_path, services, bdf_homotopy));
  const std::string bdf_homotopy_text = read(bdf_homotopy_path);
  passed &= bdf_homotopy_text.find(
                "\"low_state\":\"bdf_accepted_rate_homotopy\"") !=
                std::string::npos &&
            bdf_homotopy_text.find("\"mass_flux_scale\":1") !=
                std::string::npos &&
            bdf_homotopy_text.find("\"bdf_endpoint_alpha\":0.375") !=
                std::string::npos;
  RuntimeEvidenceRecord implicit_upwind = scaled_euler;
  implicit_upwind.predictor_mass_flux_scale = 1.0;
  implicit_upwind.predictor_low_state = 3U;
  implicit_upwind.predictor_enthalpy_solve_calls = 1U;
  implicit_upwind.predictor_enthalpy_endpoint.status = {};
  implicit_upwind.predictor_enthalpy_endpoint.termination =
      LinearTermination::converged;
  implicit_upwind.predictor_enthalpy_endpoint.iterations = 2U;
  implicit_upwind.predictor_enthalpy_endpoint.reduction_calls = 7U;
  implicit_upwind.predictor_enthalpy_endpoint.operator_applies = 4U;
  implicit_upwind.predictor_enthalpy_endpoint.preconditioner_applies = 2U;
  const fs::path implicit_upwind_path =
      root / "evidence" / "runtime-implicit-upwind.jsonl";
  passed &= static_cast<bool>(EvidenceWriter::append(
      MPI_COMM_SELF, implicit_upwind_path, services, implicit_upwind));
  const std::string implicit_upwind_text = read(implicit_upwind_path);
  passed &= implicit_upwind_text.find("\"low_state\":\"implicit_upwind\"") !=
                std::string::npos &&
            implicit_upwind_text.find("\"enthalpy_endpoint_alpha\":1") !=
                std::string::npos &&
            implicit_upwind_text.find("\"enthalpy_solve_calls\":1") !=
                std::string::npos &&
            implicit_upwind_text.find("\"operator_applies\":4") !=
                std::string::npos;
  RuntimeEvidenceRecord source_limited = implicit_upwind;
  source_limited.predictor_low_state = 4U;
  source_limited.predictor_enthalpy_endpoint_alpha = 0.625;
  const fs::path source_limited_path =
      root / "evidence" / "runtime-implicit-source-limited.jsonl";
  passed &= static_cast<bool>(EvidenceWriter::append(
      MPI_COMM_SELF, source_limited_path, services, source_limited));
  const std::string source_limited_text = read(source_limited_path);
  passed &= source_limited_text.find(
                "\"low_state\":\"implicit_upwind_source_limited\"") !=
                std::string::npos &&
            source_limited_text.find(
                "\"enthalpy_endpoint_alpha\":0.625") !=
                std::string::npos;
  RuntimeEvidenceRecord source_bundle = scaled_euler;
  source_bundle.predictor_low_state = 7U;
  source_bundle.predictor_low_order_halo_exchanges = 1U;
  source_bundle.predictor_source_endpoint_alpha = 0.375;
  const fs::path source_bundle_path =
      root / "evidence" / "runtime-source-bundle-limited.jsonl";
  passed &= static_cast<bool>(EvidenceWriter::append(
      MPI_COMM_SELF, source_bundle_path, services, source_bundle));
  const std::string source_bundle_text = read(source_bundle_path);
  passed &= source_bundle_text.find(
                "\"low_state\":\"bdf_local_donor_flux_source_limited\"") !=
                std::string::npos &&
            source_bundle_text.find("\"source_endpoint_alpha\":0.375") !=
                std::string::npos &&
            source_bundle_text.find("\"bdf_endpoint_alpha\":1") !=
                std::string::npos &&
            source_bundle_text.find("\"enthalpy_endpoint_alpha\":1") !=
                std::string::npos &&
            source_bundle_text.find("\"enthalpy_solve_calls\":0") !=
                std::string::npos &&
            source_bundle_text.find("\"low_order_halo_exchanges\":1") !=
                std::string::npos;
  RuntimeEvidenceRecord invalid_bdf_scale = scaled_euler;
  invalid_bdf_scale.predictor_low_state = 1U;
  passed &= !EvidenceWriter::append(
      MPI_COMM_SELF, root / "evidence" / "invalid-bdf-scale.jsonl",
      services, invalid_bdf_scale);
  RuntimeEvidenceRecord invalid_homotopy_alpha = bdf_homotopy;
  invalid_homotopy_alpha.predictor_bdf_endpoint_alpha = 1.0;
  passed &= !EvidenceWriter::append(
      MPI_COMM_SELF, root / "evidence" / "invalid-homotopy-alpha.jsonl",
      services, invalid_homotopy_alpha);
  RuntimeEvidenceRecord invalid_source_bundle_alpha = source_bundle;
  invalid_source_bundle_alpha.predictor_source_endpoint_alpha = 1.0;
  passed &= !EvidenceWriter::append(
      MPI_COMM_SELF, root / "evidence" / "invalid-source-bundle-alpha.jsonl",
      services, invalid_source_bundle_alpha);
  RuntimeEvidenceRecord invalid_foreign_source_alpha = scaled_euler;
  invalid_foreign_source_alpha.predictor_source_endpoint_alpha = 0.5;
  passed &= !EvidenceWriter::append(
      MPI_COMM_SELF, root / "evidence" / "invalid-foreign-source-alpha.jsonl",
      services, invalid_foreign_source_alpha);
  RuntimeEvidenceRecord invalid_nonunit_bdf_alpha = source_bundle;
  invalid_nonunit_bdf_alpha.predictor_bdf_endpoint_alpha = 0.5;
  passed &= !EvidenceWriter::append(
      MPI_COMM_SELF, root / "evidence" / "invalid-nonunit-bdf-alpha.jsonl",
      services, invalid_nonunit_bdf_alpha);

  const std::array<std::string_view, 16U> serialized_pressure_keys{{
      "\"reduction_calls\":0", "\"operator_applies\":0",
      "\"preconditioner_applies\":0", "\"norm_breakdown_restarts\":0",
      "\"recycle_offered_directions\":0",
      "\"recycle_retained_directions\":0",
      "\"recycle_operator_applies\":0",
      "\"recycle_reduction_calls\":0",
      "\"recycle_projection_attempted\":false",
      "\"recycle_projection_accepted\":false",
      "\"recycle_projected_true_residual\":0",
      "\"recycle_cycle_corrections\":0",
      "\"recycle_capture_vector_passes\":0",
      "\"recycle_capture_cycle_attempts\":0",
      "\"recycle_capture_reduction_calls\":0",
      "\"recycle_capture_blocking_operations\":0"}};
  for (const std::string_view key : serialized_pressure_keys)
    passed &= evidence_text.find(key) != std::string::npos;

  RuntimeEvidenceRecord legal_capture = evidence;
  legal_capture.pressure[0U].recycle_cycle_corrections = 1U;
  legal_capture.pressure[0U].recycle_capture_vector_passes = 4U;
  legal_capture.pressure[0U].recycle_capture_cycle_attempts = 2U;
  legal_capture.pressure[0U].recycle_capture_reduction_calls = 2U;
  legal_capture.pressure[0U].recycle_capture_blocking_operations = 4U;
  legal_capture.pressure[1U].recycle_offered_directions = 1U;
  passed &= static_cast<bool>(EvidenceWriter::append(
      MPI_COMM_SELF, root / "evidence" / "legal-capture.jsonl", services,
      legal_capture));

  RuntimeEvidenceRecord legal_accepted = evidence;
  legal_accepted.pressure[0U].recycle_cycle_corrections = 2U;
  legal_accepted.pressure[0U].recycle_capture_vector_passes = 4U;
  legal_accepted.pressure[0U].recycle_capture_cycle_attempts = 2U;
  legal_accepted.pressure[0U].recycle_capture_reduction_calls = 2U;
  legal_accepted.pressure[0U].recycle_capture_blocking_operations = 4U;
  legal_accepted.pressure[1U].initial_true_residual = 1.0;
  legal_accepted.pressure[1U].final_true_residual = 0.25;
  legal_accepted.pressure[1U].recursive_residual = 0.25;
  legal_accepted.pressure[1U].recycle_offered_directions = 2U;
  legal_accepted.pressure[1U].recycle_retained_directions = 1U;
  legal_accepted.pressure[1U].recycle_operator_applies = 3U;
  legal_accepted.pressure[1U].recycle_reduction_calls = 7U;
  legal_accepted.pressure[1U].recycle_projection_attempted = true;
  legal_accepted.pressure[1U].recycle_projection_accepted = true;
  legal_accepted.pressure[1U].recycle_projected_true_residual = 0.25;
  legal_accepted.pressure[1U].iterations = 0U;
  passed &= static_cast<bool>(EvidenceWriter::append(
      MPI_COMM_SELF, root / "evidence" / "legal-accepted.jsonl", services,
      legal_accepted));

  // Exact values from the frozen identity-operator projection fixture:
  // r0=e0+e1, U=span(e0), projected residual=e1, followed by one ordinary
  // FGMRES iteration and the canonical supplemental audit.
  RuntimeEvidenceRecord legal_accepted_ordinary = evidence;
  legal_accepted_ordinary.pressure[0U].recycle_cycle_corrections = 1U;
  legal_accepted_ordinary.pressure[0U].recycle_capture_vector_passes = 2U;
  legal_accepted_ordinary.pressure[0U].recycle_capture_cycle_attempts = 1U;
  legal_accepted_ordinary.pressure[0U].recycle_capture_reduction_calls = 1U;
  legal_accepted_ordinary.pressure[0U]
      .recycle_capture_blocking_operations = 2U;
  legal_accepted_ordinary.pressure[1U].iterations = 1U;
  legal_accepted_ordinary.pressure[1U].initial_true_residual =
      1.4142135623730951;
  legal_accepted_ordinary.pressure[1U].final_true_residual = 0.0;
  legal_accepted_ordinary.pressure[1U].recursive_residual = 0.0;
  legal_accepted_ordinary.pressure[1U].operator_applies = 3U;
  legal_accepted_ordinary.pressure[1U].preconditioner_applies = 1U;
  legal_accepted_ordinary.pressure[1U].reduction_calls = 7U;
  legal_accepted_ordinary.pressure[1U].recycle_offered_directions = 1U;
  legal_accepted_ordinary.pressure[1U].recycle_retained_directions = 1U;
  legal_accepted_ordinary.pressure[1U].recycle_operator_applies = 2U;
  legal_accepted_ordinary.pressure[1U].recycle_reduction_calls = 5U;
  legal_accepted_ordinary.pressure[1U].recycle_projection_attempted = true;
  legal_accepted_ordinary.pressure[1U].recycle_projection_accepted = true;
  legal_accepted_ordinary.pressure[1U].recycle_projected_true_residual = 1.0;
  passed &= expect(
      static_cast<bool>(EvidenceWriter::append(
          MPI_COMM_SELF,
          root / "evidence" / "legal-accepted-ordinary.jsonl", services,
          legal_accepted_ordinary)),
      "accepted projection may continue through ordinary FGMRES iterations");

  RuntimeEvidenceRecord legal_all_deflated = evidence;
  legal_all_deflated.pressure[0U].recycle_cycle_corrections = 2U;
  legal_all_deflated.pressure[0U].recycle_capture_vector_passes = 4U;
  legal_all_deflated.pressure[0U].recycle_capture_cycle_attempts = 2U;
  legal_all_deflated.pressure[0U].recycle_capture_reduction_calls = 2U;
  legal_all_deflated.pressure[0U].recycle_capture_blocking_operations = 4U;
  legal_all_deflated.pressure[1U].initial_true_residual = 1.0;
  legal_all_deflated.pressure[1U].final_true_residual = 1.0;
  legal_all_deflated.pressure[1U].recursive_residual = 1.0;
  legal_all_deflated.pressure[1U].recycle_offered_directions = 2U;
  legal_all_deflated.pressure[1U].recycle_operator_applies = 2U;
  legal_all_deflated.pressure[1U].recycle_reduction_calls = 2U;
  legal_all_deflated.pressure[1U].recycle_projection_attempted = true;
  legal_all_deflated.pressure[1U].recycle_projected_true_residual = 0.0;
  passed &= expect(
      static_cast<bool>(EvidenceWriter::append(
          MPI_COMM_SELF, root / "evidence" / "legal-all-deflated.jsonl",
          services, legal_all_deflated)),
      "all-deflated exact-zero projection has offered reductions only");

  RuntimeEvidenceRecord legal_non_improving = evidence;
  legal_non_improving.pressure[0U].recycle_cycle_corrections = 1U;
  legal_non_improving.pressure[0U].recycle_capture_vector_passes = 2U;
  legal_non_improving.pressure[0U].recycle_capture_cycle_attempts = 1U;
  legal_non_improving.pressure[0U].recycle_capture_reduction_calls = 1U;
  legal_non_improving.pressure[0U].recycle_capture_blocking_operations = 2U;
  legal_non_improving.pressure[1U].initial_true_residual = 1.0;
  legal_non_improving.pressure[1U].final_true_residual = 1.0;
  legal_non_improving.pressure[1U].recursive_residual = 1.0;
  legal_non_improving.pressure[1U].recycle_offered_directions = 1U;
  legal_non_improving.pressure[1U].recycle_retained_directions = 1U;
  legal_non_improving.pressure[1U].recycle_operator_applies = 2U;
  legal_non_improving.pressure[1U].recycle_reduction_calls = 5U;
  legal_non_improving.pressure[1U].recycle_projection_attempted = true;
  legal_non_improving.pressure[1U].recycle_projected_true_residual = 1.0;
  passed &= static_cast<bool>(EvidenceWriter::append(
      MPI_COMM_SELF, root / "evidence" / "legal-non-improving.jsonl",
      services, legal_non_improving));
  const std::string accepted_text =
      read(root / "evidence" / "legal-accepted.jsonl");
  const std::string capture_text =
      read(root / "evidence" / "legal-capture.jsonl");
  const std::string deflated_text =
      read(root / "evidence" / "legal-all-deflated.jsonl");
  const std::string non_improving_text =
      read(root / "evidence" / "legal-non-improving.jsonl");
  passed &= accepted_text.find("\"recycle_projection_accepted\":true") !=
                std::string::npos &&
            accepted_text.find("\"recycle_operator_applies\":3") !=
                std::string::npos &&
            capture_text.find("\"recycle_capture_cycle_attempts\":2") !=
                std::string::npos &&
            deflated_text.find("\"recycle_retained_directions\":0") !=
                std::string::npos &&
            non_improving_text.find(
                "\"recycle_projection_accepted\":false") !=
                std::string::npos;

  RuntimeEvidenceRecord skipped_projection = evidence;
  skipped_projection.pressure[0U].recycle_cycle_corrections = 2U;
  skipped_projection.pressure[0U].recycle_capture_vector_passes = 4U;
  skipped_projection.pressure[0U].recycle_capture_cycle_attempts = 2U;
  skipped_projection.pressure[0U].recycle_capture_reduction_calls = 2U;
  skipped_projection.pressure[0U].recycle_capture_blocking_operations = 4U;
  skipped_projection.pressure[1U].recycle_offered_directions = 2U;
  passed &= static_cast<bool>(EvidenceWriter::append(
      MPI_COMM_SELF, root / "evidence" / "skipped-projection.jsonl",
      services, skipped_projection));
  RuntimeEvidenceRecord legal_capped = evidence;
  legal_capped.pressure[0U].recycle_cycle_corrections = 5U;
  legal_capped.pressure[0U].recycle_capture_vector_passes = 10U;
  legal_capped.pressure[0U].recycle_capture_cycle_attempts = 5U;
  legal_capped.pressure[0U].recycle_capture_reduction_calls = 5U;
  legal_capped.pressure[0U].recycle_capture_blocking_operations = 10U;
  legal_capped.pressure[1U].recycle_offered_directions = 4U;
  passed &= static_cast<bool>(EvidenceWriter::append(
      MPI_COMM_SELF, root / "evidence" / "legal-capped.jsonl", services,
      legal_capped));
  RuntimeEvidenceRecord invalid_projection = skipped_projection;
  invalid_projection.pressure[1U].recycle_operator_applies = 1U;
  passed &= EvidenceWriter::append(
                MPI_COMM_SELF, root / "evidence" / "invalid-projection.jsonl",
                services, invalid_projection)
                .code == StatusCode::invalid_plan;
  RuntimeEvidenceRecord invalid_operator_formula = legal_accepted;
  invalid_operator_formula.pressure[1U].recycle_operator_applies = 2U;
  passed &= EvidenceWriter::append(
                MPI_COMM_SELF,
                root / "evidence" / "invalid-operator-formula.jsonl", services,
                invalid_operator_formula)
                .code == StatusCode::invalid_plan;
  RuntimeEvidenceRecord improved_but_rejected = legal_accepted;
  improved_but_rejected.pressure[1U].recycle_projection_accepted = false;
  passed &= expect(
      EvidenceWriter::append(
          MPI_COMM_SELF,
          root / "evidence" / "improved-but-rejected.jsonl", services,
          improved_but_rejected)
              .code == StatusCode::invalid_plan,
      "improved projection cannot be marked rejected");
  RuntimeEvidenceRecord non_improving_but_accepted = legal_non_improving;
  non_improving_but_accepted.pressure[1U].recycle_projection_accepted = true;
  passed &= expect(
      EvidenceWriter::append(
          MPI_COMM_SELF,
          root / "evidence" / "non-improving-but-accepted.jsonl", services,
          non_improving_but_accepted)
              .code == StatusCode::invalid_plan,
      "non-improving projection cannot be marked accepted");
  RuntimeEvidenceRecord zero_iteration_mismatch = legal_accepted;
  zero_iteration_mismatch.pressure[1U].final_true_residual = 0.1;
  passed &= expect(
      EvidenceWriter::append(
          MPI_COMM_SELF,
          root / "evidence" / "zero-iteration-mismatch.jsonl", services,
          zero_iteration_mismatch)
              .code == StatusCode::invalid_plan,
      "accepted zero-iteration projection requires final projected residual");
  RuntimeEvidenceRecord cap_exceeded = legal_capped;
  cap_exceeded.pressure[1U].recycle_offered_directions = 5U;
  passed &= EvidenceWriter::append(
                MPI_COMM_SELF, root / "evidence" / "cap-exceeded.jsonl",
                services, cap_exceeded)
                .code == StatusCode::invalid_plan;
  RuntimeEvidenceRecord offered_too_few = skipped_projection;
  offered_too_few.pressure[1U].recycle_offered_directions = 1U;
  passed &= EvidenceWriter::append(
                MPI_COMM_SELF, root / "evidence" / "offered-too-few.jsonl",
                services, offered_too_few)
                .code == StatusCode::invalid_plan;
  RuntimeEvidenceRecord offered_too_many = skipped_projection;
  offered_too_many.pressure[1U].recycle_offered_directions = 3U;
  passed &= EvidenceWriter::append(
                MPI_COMM_SELF, root / "evidence" / "offered-too-many.jsonl",
                services, offered_too_many)
                .code == StatusCode::invalid_plan;
  RuntimeEvidenceRecord invalid_capture = evidence;
  invalid_capture.pressure[0U].recycle_capture_cycle_attempts = 1U;
  invalid_capture.pressure[0U].recycle_capture_vector_passes = 0U;
  passed &= EvidenceWriter::append(
                MPI_COMM_SELF, root / "evidence" / "invalid-capture.jsonl",
                services, invalid_capture)
                .code == StatusCode::invalid_plan;
  RuntimeEvidenceRecord invalid_momentum_limiter = evidence;
  invalid_momentum_limiter.momentum_predictor_theta = 0.5;
  passed &= EvidenceWriter::append(
                MPI_COMM_SELF,
                root / "evidence" / "invalid-momentum-limiter.jsonl",
                services, invalid_momentum_limiter)
                .code == StatusCode::invalid_plan;
  RuntimeEvidenceRecord limited_momentum_faces = evidence;
  limited_momentum_faces.momentum_predictor_limited = true;
  limited_momentum_faces.momentum_predictor_theta = 0.5;
  limited_momentum_faces.momentum_predictor_activations = 7U;
  const fs::path limited_momentum_faces_path =
      root / "evidence" / "limited-momentum-faces.jsonl";
  passed &= static_cast<bool>(EvidenceWriter::append(
      MPI_COMM_SELF, limited_momentum_faces_path, services,
      limited_momentum_faces));
  const std::string limited_momentum_faces_text =
      read(limited_momentum_faces_path);
  passed &= limited_momentum_faces_text.find(
                "\"retained_correction_l1_ratio\":0.5,\"limited_faces\":7") !=
            std::string::npos;
  RuntimeEvidenceRecord projection_without_capture = evidence;
  projection_without_capture.pressure[1U].recycle_offered_directions = 1U;
  passed &= EvidenceWriter::append(
                MPI_COMM_SELF, root / "evidence" /
                                    "projection-without-capture.jsonl",
                services, projection_without_capture)
                .code == StatusCode::invalid_plan;

  CommittedOutputSnapshot trial = snapshot;
  trial.committed = false;
  passed &= VisitWriter::write(MPI_COMM_SELF, root / "trial", services,
                               trial)
                .code == StatusCode::invalid_plan;
  auto stale_fields = fields;
  stale_fields[0U].accepted_revision += 1U;
  CommittedOutputSnapshot stale = snapshot;
  stale.fields = {stale_fields.data(), stale_fields.size()};
  passed &= ScreenWriter::append(MPI_COMM_SELF, root / "stale.log", services,
                                 stale, "must reject")
                .code == StatusCode::invalid_plan;
  evidence.statistics_eligible = true;
  passed &= EvidenceWriter::append(MPI_COMM_SELF, root / "invalid.jsonl",
                                   services, evidence)
                .code == StatusCode::invalid_plan;
  if (!passed) std::cerr << "run_case failed stretched=" << stretched << '\n';
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 2;
  const fs::path root = fs::temp_directory_path() /
                        ("hundun-v04-io-product-" +
                         std::to_string(::getpid()));
  std::error_code error;
  fs::remove_all(root, error);
  fs::create_directories(root);
  const bool passed = run_case(root, false) && run_case(root, true);
  if (passed) fs::remove_all(root, error);
  if (!passed) std::cerr << "committed output path failure\n";
  MPI_Finalize();
  return passed ? 0 : 1;
}
