// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "tests/support/stage3_performance_evidence.hpp"

#include "hundun/bc_basic_boundary.hpp"
#include "hundun/exec_execution.hpp"
#include "hundun/flow_checkpoint_v3.hpp"
#include "hundun/flow_immersed.hpp"
#include "hundun/fvm_cell_centered.hpp"
#include "hundun/ib_domain.hpp"
#include "hundun/ib_ghost_stencil_plan.hpp"
#include "hundun/ib_local_flow_pattern.hpp"
#include "hundun/ib_surface.hpp"
#include "hundun/ib_surface_query.hpp"
#include "hundun/lin_bicgstab.hpp"
#include "hundun/lin_conjugate_gradient.hpp"
#include "hundun/lin_preconditioners.hpp"
#include "hundun/mesh_geometry.hpp"
#include "hundun/mesh_topology.hpp"
#include "hundun/rt_exchange_plan.hpp"
#include "hundun/rt_field_registry.hpp"
#include "hundun/rt_halo_exchange.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_structured_decomposition.hpp"
#include "tests/support/flow_immersed_test_access.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <optional>
#include <string>

namespace hundun::test {
namespace {

runtime::FieldDescriptor cell_field(const char *name,
                                    std::uint32_t components,
                                    int ghost_width) {
  return {name,
          "1",
          "stage3_performance",
          runtime::FunctionSpace::cell_average,
          runtime::ScalarType::float64,
          components,
          ghost_width,
          false,
          runtime::RestartPolicy::persistent,
          runtime::OutputPolicy::never};
}

runtime::FieldDescriptor face_field(const char *name,
                                    std::uint32_t components) {
  return {name,
          "1",
          "stage3_performance",
          runtime::FunctionSpace::face_value,
          runtime::ScalarType::float64,
          components,
          0,
          false,
          runtime::RestartPolicy::persistent,
          runtime::OutputPolicy::never};
}

diagnostics::Stage3PerformanceCounters snapshot(
    const runtime::MpiContext &mpi,
    const immersed::ImmersedSurface &surface,
    const immersed::SurfaceQuery &query, const immersed::ImmersedDomain &domain,
    const immersed::GhostStencilPlan &ghost,
    const immersed::WallQuadraturePlan &wall,
    const flow::FixedStepImmersedFlow &flow, const les::WaleModel &wale) {
  diagnostics::Stage3PerformanceCounters result;
  result.init_surface_triangles = mpi.rank() == 0 ? surface.triangle_count() : 0U;
  for (const auto &value : {query.performance_counters(),
                            domain.performance_counters(),
                            ghost.performance_counters(),
                            wall.performance_counters(),
                            flow.performance_counters(),
                            wale.performance_counters(),
                            flow::checkpoint_v3_performance_counters()})
    result = diagnostics::add_stage3_performance_counters(result, value);
  std::array<std::uint64_t, 17> values{
      result.init_surface_triangles,
      result.init_query_closest_calls,
      result.init_query_segment_calls,
      result.init_classification_cells,
      result.init_ghost_qr_plans,
      result.init_ghost_rejected_plans,
      result.init_ghost_donor_references,
      result.init_wall_points,
      result.step_ghost_constraints,
      result.step_lfp_transforms,
      result.step_immersed_rows,
      result.step_pressure_wall_constraints,
      result.step_wall_quadrature_evaluations,
      result.step_force_reductions,
      result.step_wale_gradient_cells,
      result.step_wale_evaluations,
      result.checkpoint_logical_io_bytes};
  std::uint64_t logical_operations[2]{result.step_force_reductions,
                                      result.step_wale_evaluations};
  if (MPI_Allreduce(MPI_IN_PLACE, values.data(),
                    static_cast<int>(values.size()), MPI_UINT64_T, MPI_SUM,
                    mpi.comm()) != MPI_SUCCESS)
    throw std::runtime_error("unable to aggregate Stage 3 work counters");
  if (MPI_Allreduce(MPI_IN_PLACE, logical_operations, 2, MPI_UINT64_T, MPI_MAX,
                    mpi.comm()) != MPI_SUCCESS)
    throw std::runtime_error(
        "unable to aggregate Stage 3 logical work counters");
  values[13] = logical_operations[0];
  values[15] = logical_operations[1];
  return {values[0], values[1], values[2], values[3], values[4], values[5],
          values[6], values[7], values[8], values[9], values[10], values[11],
          values[12], values[13], values[14], values[15], values[16]};
}

} // namespace

Stage3PerformanceEvidence
run_stage3_performance_evidence(const runtime::MpiContext &mpi,
                                const Stage3PerformanceRun &run) {
  if (run.cells <= 0 || run.measured_steps == 0U || run.repetitions != 1)
    throw std::invalid_argument("invalid Stage 3 performance run contract");
  const int cells = run.cells;
  Stage3TemporaryDirectory directory("performance-evidence");
  auto decomposition = runtime::StructuredDecomposition::create(
      mpi, {cells, cells, cells}, {true, true, true},
      runtime::DecompositionOptions{
          stage3_performance_process_grid(mpi.size())});
  mesh::MeshTopology topology(decomposition);
  mesh::MeshGeometry geometry(
      topology,
      mesh::UniformBoxMapping({-0.5, -0.5, -0.5}, {2.0, 2.0, 2.0}));
  auto boundaries = boundary::BoundaryRegistry::create(
      stage3_performance_case(mpi.size(), cells), topology);

  std::string path_text;
  if (mpi.rank() == 0) {
    const auto path = directory.path() / "body.stl";
    write_text(path, ascii_stl(stage3_performance_body(), "body"));
    path_text = path.string();
  }
  std::uint64_t path_size = path_text.size();
  if (MPI_Bcast(&path_size, 1, MPI_UINT64_T, 0, mpi.comm()) != MPI_SUCCESS)
    throw std::runtime_error("unable to broadcast performance STL path size");
  path_text.resize(static_cast<std::size_t>(path_size));
  if (MPI_Bcast(path_text.data(), static_cast<int>(path_text.size()), MPI_BYTE,
                0, mpi.comm()) != MPI_SUCCESS)
    throw std::runtime_error("unable to broadcast performance STL path");

  auto surface = immersed::ImmersedSurface::load_collective(
      std::filesystem::path(path_text), 1.0, mpi, 0);
  auto query = immersed::SurfaceQuery::create(surface);
  auto domain = immersed::ImmersedDomain::create(
      surface, query, config::ImmersedFluidSide::outside, topology, geometry,
      boundaries, mpi);
  auto ghost = immersed::GhostStencilPlan::create(
      surface, query, domain, topology, geometry, decomposition, mpi);
  auto wall = immersed::WallQuadraturePlan::create(
      surface, query, domain, topology, geometry, mpi);
  immersed::LocalFlowPatternTransform transform;
  const int halo_width = static_cast<int>(std::max(
      ghost.maximum_halo_reach(), wall.maximum_halo_reach()));

  runtime::FieldRegistry registry;
  flow::FlowFieldIds fields;
  fields.density = registry.declare_field(cell_field("rho", 1U, halo_width));
  fields.velocity =
      registry.declare_field(cell_field("velocity", 3U, halo_width));
  fields.mechanical_pressure =
      registry.declare_field(cell_field("pi", 1U, halo_width));
  fields.face_velocity =
      registry.declare_field(face_field("face_velocity", 3U));
  fields.face_mass_flux = finite_volume::declare_face_mass_flux(registry);
  registry.freeze();
  auto state = flow::FlowState::create(
      registry, {decomposition.local_extent(), topology.local_face_count()},
      fields,
      {0U, 0.0, 0.01, 0.0, flow::MomentumTimeOrder::backward_euler});
  flow::FlowLayerValues values;
  values.density.resize(topology.owned_cell_count(), 0.0);
  values.velocity.resize(topology.owned_cell_count() * 3U, 0.0);
  values.mechanical_pressure.resize(topology.owned_cell_count(), 0.0);
  values.face_velocity.resize(topology.local_face_count() * 3U, 0.0);
  values.face_mass_flux.resize(topology.local_face_count(), 0.0);
  for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count(); ++cell)
    if (domain.region(cell) == immersed::CellRegion::fluid)
      values.density[cell] = 1.0;
  state.seed_accepted_layers(values, values);

  auto exchange_plan = runtime::ExchangePlan::create(
      decomposition, decomposition.local_extent(), halo_width);
  auto halo = runtime::HaloExchange::create(decomposition, exchange_plan);
  execution::CpuReferenceContext execution;
  linear::ConjugateGradientSolver momentum_solver(execution, mpi);
  linear::BiCGStabSolver pressure_solver(execution, mpi);
  linear::JacobiPreconditioner mx(execution), my(execution), mz(execution),
      pressure_pc(execution);
  auto wale = les::WaleModel::create(
      {0.5, 0.9, 0.7}, topology, geometry,
      domain.active_cells().owned_active_count(),
      domain.active_cells().ordered_global_ids(), execution);
  auto facade = flow::FixedStepImmersedFlow::create(
      decomposition, topology, geometry, boundaries, &domain, &ghost, &wall,
      &transform, &wale, mpi, execution, halo, momentum_solver,
      {&mx, &my, &mz}, pressure_solver, pressure_pc);

  Stage3PerformanceEvidence evidence;
  evidence.warmup_steps = run.warmup_steps;
  evidence.measured_steps = run.measured_steps;
  evidence.surface_fingerprint = surface.fingerprint();
  evidence.classification_fingerprint = domain.classification_fingerprint();
  evidence.after_initialization =
      snapshot(mpi, surface, query, domain, ghost, wall, facade, wale);
  const flow::ImmersedFlowPhysics physics{config::DensityModel::constant,
                                          1.0,
                                          0.01,
                                          std::nullopt,
                                          std::nullopt,
                                          std::nullopt};
  const auto next_stencil = [&] {
    const auto metadata = state.metadata();
    return flow::make_momentum_time_stencil(
        metadata.step == 0U ? flow::MomentumTimeOrder::backward_euler
                            : flow::MomentumTimeOrder::bdf2,
        0.01, metadata.step == 0U ? 0.0 : metadata.dt_s);
  };
  const auto accepted_step = [&] {
    const auto report = facade.attempt(state, physics, next_stencil(), {}, {});
    const auto &base = std::get<flow::StepAttemptReport>(report.base);
    if (base.disposition != flow::StepAttemptDisposition::committed)
      throw std::runtime_error("Stage 3 performance step did not commit");
    if (base.pressure_corrector_count != 2U)
      throw std::runtime_error(
          "Stage 3 performance step did not execute two PISO correctors");
    evidence.committed = true;
    evidence.pressure_corrector_count = base.pressure_corrector_count;
  };

  auto after_failed_attempt = evidence.after_initialization;
  if (run.inject_failed_attempt) {
    flow::test::ImmersedFlowTestAccess::set_failure_stage(
        flow::test::ImmersedFlowAttemptFailureStage::before_commit,
        mpi.size() > 1 ? 1 : 0);
    flow::ImmersedFlowStepAttemptReport failed;
    try {
      failed = facade.attempt(state, physics, next_stencil(), {}, {});
    } catch (...) {
      flow::test::ImmersedFlowTestAccess::clear_failure_stage();
      throw;
    }
    flow::test::ImmersedFlowTestAccess::clear_failure_stage();
    const auto &failed_base = std::get<flow::StepAttemptReport>(failed.base);
    evidence.failed_attempt_rolled_back =
        failed_base.disposition != flow::StepAttemptDisposition::committed &&
        state.metadata().step == 0U && state.metadata().time_s == 0.0;
    after_failed_attempt =
        snapshot(mpi, surface, query, domain, ghost, wall, facade, wale);
    evidence.failed_attempt_delta =
        diagnostics::stage3_performance_counter_delta(
            evidence.after_initialization, after_failed_attempt);
  }

  for (std::uint64_t step = 0U; step < run.warmup_steps; ++step)
    accepted_step();
  evidence.before_measured_step =
      snapshot(mpi, surface, query, domain, ghost, wall, facade, wale);
  mpi.barrier();
  const double begin = MPI_Wtime();
  for (std::uint64_t step = 0U; step < run.measured_steps; ++step)
    accepted_step();
  evidence.elapsed_seconds = MPI_Wtime() - begin;
  if (!(evidence.elapsed_seconds > 0.0) ||
      !std::isfinite(evidence.elapsed_seconds))
    throw std::runtime_error("invalid Stage 3 performance elapsed time");
  evidence.after_measured_step =
      snapshot(mpi, surface, query, domain, ghost, wall, facade, wale);
  evidence.measured_delta = diagnostics::stage3_performance_counter_delta(
      evidence.before_measured_step, evidence.after_measured_step);
  config::ImmersedFlowCaseConfig checkpoint_config;
  checkpoint_config.common_flow = stage3_performance_case(mpi.size(), cells);
  checkpoint_config.immersed_boundary.model =
      config::ImmersedBoundaryModel::local_flow_pattern_ghost_cell;
  checkpoint_config.immersed_boundary.geometry = config::StlGeometryConfig{
      "body.stl", 1.0, config::ImmersedFluidSide::outside};
  checkpoint_config.immersed_boundary.wall =
      config::StaticImmersedWallConfig{};
  checkpoint_config.les.model = config::LesModel::wale;
  checkpoint_config.les.wale = config::WaleConfig{0.5, 0.9, 0.7};
  flow::CheckpointV3WriteModules modules;
  modules.presence = flow::CheckpointV3Presence::constant_static_ibm_wale;
  modules.surface = &surface;
  modules.query = &query;
  modules.domain = &domain;
  modules.ghost_plan = &ghost;
  modules.wall_plan = &wall;
  modules.transform = &transform;
  modules.wale = &wale;
  modules.flow = &facade;
  const auto before_checkpoint = evidence.after_measured_step;
  const auto checkpoint = flow::write_checkpoint_v3(
      mpi, decomposition, topology, geometry, boundaries, checkpoint_config,
      modules, state, {0.01, 0U},
      std::filesystem::path(path_text).parent_path() / "checkpoint");
  if (checkpoint.disposition() != flow::CheckpointV3Disposition::completed)
    throw std::runtime_error(
        "Stage 3 performance checkpoint failed: reason=" +
        std::to_string(static_cast<int>(checkpoint.reason())) +
        " phase=" + std::to_string(static_cast<int>(checkpoint.phase())) +
        " rank=" + std::to_string(checkpoint.lowest_failing_rank()));
  const auto after_checkpoint =
      snapshot(mpi, surface, query, domain, ghost, wall, facade, wale);
  evidence.checkpoint_delta = diagnostics::stage3_performance_counter_delta(
      before_checkpoint, after_checkpoint);
  std::uint64_t counts[3]{domain.links().size(),
                          domain.active_cells().owned_active_count(),
                          wall.local_points().size()};
  if (MPI_Allreduce(MPI_IN_PLACE, counts, 3, MPI_UINT64_T, MPI_SUM,
                    mpi.comm()) != MPI_SUCCESS)
    throw std::runtime_error("unable to aggregate Stage 3 evidence counts");
  evidence.immersed_link_count = counts[0];
  evidence.owned_active_cell_count = counts[1];
  evidence.local_wall_point_count = counts[2];
  return evidence;
}

} // namespace hundun::test
