// SPDX-License-Identifier: Apache-2.0

#include "tests/support/flow_immersed_test_access.hpp"
#include "hundun/bc_basic_boundary.hpp"
#include "hundun/cfg_resolved_case.hpp"
#include "hundun/cfg_resolved_case_v3.hpp"
#include "hundun/diag_checkpoint_v3.hpp"
#include "hundun/diag_immersed_module.hpp"
#include "hundun/exec_execution.hpp"
#include "hundun/fvm_cell_centered.hpp"
#include "hundun/flow_checkpoint_v3.hpp"
#include "hundun/flow_immersed.hpp"
#include "hundun/ib_ghost_stencil_plan.hpp"
#include "hundun/ib_domain.hpp"
#include "hundun/ib_surface.hpp"
#include "hundun/ib_local_flow_pattern.hpp"
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
#include "hundun/rt_mpi_environment.hpp"
#include "hundun/rt_structured_decomposition.hpp"
#include "tests/support/stage3_state_equality.hpp"
#include "tests/support/stage3_stl_fixture.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace {

using namespace hundun;

class DiagnosticSink final : public diagnostics::DiagnosticSink {
public:
  void submit(const diagnostics::DiagnosticRecord &record) override {
    records.push_back(record);
  }
  std::vector<diagnostics::DiagnosticRecord> records;
};

bool has_metric(const diagnostics::DiagnosticRecord &record,
                std::string_view id) {
  return std::any_of(record.metrics.begin(), record.metrics.end(),
                     [&](const auto &metric) { return metric.id == id; });
}

runtime::Int3 process_grid(int ranks) {
  if (ranks == 1)
    return {1, 1, 1};
  if (ranks == 2)
    return {2, 1, 1};
  if (ranks == 4)
    return {2, 2, 1};
  throw runtime::Error("unsupported Task 11 transaction rank count");
}

config::FlowCaseConfig periodic_case(int ranks) {
  config::FlowCaseConfig config{};
  config.schema_version = 2;
  config.simulation_type = config::SimulationType::variable_density_flow;
  config.density_model = config::DensityModel::constant;
  config.resources.expected_ranks = ranks;
  config.resources.process_grid = process_grid(ranks);
  config.mesh.cells = {12, 12, 12};
  config.mesh.origin_m = {0.0, 0.0, 0.0};
  config.mesh.length_m = {1.0, 1.0, 1.0};
  config.physics.rho_ref_kg_per_m3 = 1.0;
  config.physics.dynamic_viscosity_pa_s = 0.01;
  config.physics.inlet_consistency_rtol = 1.0e-12;
  constexpr std::array<config::PatchName, 6> names{
      config::PatchName::x_min, config::PatchName::x_max,
      config::PatchName::y_min, config::PatchName::y_max,
      config::PatchName::z_min, config::PatchName::z_max};
  for (std::size_t index = 0U; index < names.size(); ++index) {
    config.boundaries[index].patch = names[index];
    config.boundaries[index].type = config::BoundaryType::periodic;
  }
  return config;
}

config::ImmersedFlowCaseConfig checkpoint_v3_case(int ranks) {
  config::ImmersedFlowCaseConfig result;
  result.common_flow = periodic_case(ranks);
  result.common_flow.case_name = "task17a-checkpoint-v3";
  result.common_flow.time.mode = config::TimeMode::fixed;
  result.common_flow.time.steps = 2U;
  result.common_flow.time.initial_dt_s = 0.01;
  result.common_flow.time.min_dt_s = 0.01;
  result.common_flow.time.max_dt_s = 0.01;
  result.common_flow.time.cfl_target = 0.5;
  result.common_flow.time.diffusion_number_target = 0.25;
  result.common_flow.time.growth_factor = 1.25;
  result.common_flow.time.retry_factor = 0.5;
  result.common_flow.time.max_retries = 8U;
  result.common_flow.restart.read = false;
  result.common_flow.restart.write_directory = "checkpoints";
  result.common_flow.restart.write_interval = 1;
  result.common_flow.diagnostics.directory = "diagnostics";
  result.common_flow.diagnostics.write_interval = 1;
  result.common_flow.diagnostics.write_mesh = false;
  result.common_flow.performance.enabled = false;
  result.common_flow.performance.directory = "performance";
  result.common_flow.performance.warmup_steps = 5;
  result.common_flow.performance.measured_steps = 20;
  result.common_flow.performance.repetitions = 5;
  result.immersed_boundary.model =
      config::ImmersedBoundaryModel::local_flow_pattern_ghost_cell;
  result.immersed_boundary.geometry = config::StlGeometryConfig{
      "body.stl", 1.0, config::ImmersedFluidSide::outside};
  result.immersed_boundary.wall = config::StaticImmersedWallConfig{};
  result.les.model = config::LesModel::none;
  return result;
}

runtime::FieldDescriptor cell_field(const char *name,
                                    std::uint32_t components,
                                    int ghost_width) {
  return {name,
          "1",
          "stage3_task11",
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
          "stage3_task11",
          runtime::FunctionSpace::face_value,
          runtime::ScalarType::float64,
          components,
          0,
          false,
          runtime::RestartPolicy::persistent,
          runtime::OutputPolicy::never};
}

std::vector<test::StlFixtureTriangle> internal_cube() {
  auto coarse = test::outward_cube();
  for (auto &triangle : coarse)
    for (auto &vertex : triangle.vertices) {
      vertex.x = 0.30 + 0.40 * vertex.x;
      vertex.y = 0.30 + 0.40 * vertex.y;
      vertex.z = 0.30 + 0.40 * vertex.z;
    }
  const auto midpoint = [](runtime::Real3 left, runtime::Real3 right) {
    return runtime::Real3{0.5 * (left.x + right.x), 0.5 * (left.y + right.y),
                          0.5 * (left.z + right.z)};
  };
  std::vector<test::StlFixtureTriangle> refined;
  refined.reserve(4U * coarse.size());
  for (const auto &triangle : coarse) {
    const auto ab = midpoint(triangle.vertices[0], triangle.vertices[1]);
    const auto bc = midpoint(triangle.vertices[1], triangle.vertices[2]);
    const auto ca = midpoint(triangle.vertices[2], triangle.vertices[0]);
    refined.push_back({triangle.file_normal, {triangle.vertices[0], ab, ca}});
    refined.push_back({triangle.file_normal, {ab, triangle.vertices[1], bc}});
    refined.push_back({triangle.file_normal, {ca, bc, triangle.vertices[2]}});
    refined.push_back({triangle.file_normal, {ab, bc, ca}});
  }
  return refined;
}

struct Fixture final {
  explicit Fixture(const runtime::MpiContext &mpi)
      : directory("task11-immersed-transaction"),
        decomposition(runtime::StructuredDecomposition::create(
            mpi, {12, 12, 12}, {true, true, true},
            runtime::DecompositionOptions{process_grid(mpi.size())})),
        topology(decomposition),
        geometry(topology,
                 mesh::UniformBoxMapping({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0})),
        boundaries(boundary::BoundaryRegistry::create(periodic_case(mpi.size()),
                                                      topology)) {
    std::string path_text;
    if (mpi.rank() == 0) {
      const auto path = directory.path() / "body.stl";
      test::write_text(path, test::ascii_stl(internal_cube(), "body"));
      path_text = path.string();
    }
    std::uint64_t path_size = path_text.size();
    HUNDUN_CHECK(MPI_Bcast(&path_size, 1, MPI_UINT64_T, 0, mpi.comm()) ==
                 MPI_SUCCESS);
    path_text.resize(static_cast<std::size_t>(path_size));
    HUNDUN_CHECK(MPI_Bcast(path_text.data(), static_cast<int>(path_text.size()),
                           MPI_BYTE, 0, mpi.comm()) == MPI_SUCCESS);
    const std::filesystem::path path(path_text);
    surface.emplace(
        immersed::ImmersedSurface::load_collective(path, 1.0, mpi, 0));
    query.emplace(immersed::SurfaceQuery::create(*surface));
    domain.emplace(immersed::ImmersedDomain::create(
        *surface, *query, config::ImmersedFluidSide::outside, topology,
        geometry, boundaries, mpi));
    ghost_plan.emplace(immersed::GhostStencilPlan::create(
        *surface, *query, *domain, topology, geometry, decomposition, mpi));
    wall_plan.emplace(immersed::WallQuadraturePlan::create(
        *surface, *query, *domain, topology, geometry, mpi));
    ghost_reach = static_cast<int>(ghost_plan->maximum_halo_reach());
    wall_reach = static_cast<int>(wall_plan->maximum_halo_reach());
    halo_width = std::max(ghost_reach, wall_reach);

    fields.density = registry.declare_field(cell_field("rho", 1U, halo_width));
    fields.velocity =
        registry.declare_field(cell_field("velocity", 3U, halo_width));
    fields.mechanical_pressure =
        registry.declare_field(cell_field("pi", 1U, halo_width));
    fields.face_velocity =
        registry.declare_field(face_field("face_velocity", 3U));
    fields.face_mass_flux = finite_volume::declare_face_mass_flux(registry);
    registry.freeze();
    HUNDUN_CHECK(registry.descriptor(fields.density).ghost_width == halo_width);
    HUNDUN_CHECK(registry.descriptor(fields.velocity).ghost_width ==
                 halo_width);
    HUNDUN_CHECK(registry.descriptor(fields.mechanical_pressure).ghost_width ==
                 halo_width);
  }

  flow::FlowState make_state() const {
    const auto local = decomposition.local_extent();
    auto state = flow::FlowState::create(
        registry, {local, topology.local_face_count()}, fields,
        {0U, 0.0, 0.01, 0.0, flow::MomentumTimeOrder::backward_euler});
    flow::FlowLayerValues values;
    values.density.resize(topology.owned_cell_count(), 0.0);
    values.velocity.resize(topology.owned_cell_count() * 3U, 0.0);
    values.mechanical_pressure.resize(topology.owned_cell_count(), 0.0);
    values.face_velocity.resize(topology.local_face_count() * 3U, 0.0);
    values.face_mass_flux.resize(topology.local_face_count(), 0.0);
    for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count();
         ++cell)
      if (domain->region(cell) == immersed::CellRegion::fluid)
        values.density[cell] = 1.0;
    state.seed_accepted_layers(values, values);
    return state;
  }

  test::Stage3TemporaryDirectory directory;
  runtime::StructuredDecomposition decomposition;
  mesh::MeshTopology topology;
  mesh::MeshGeometry geometry;
  boundary::BoundaryRegistry boundaries;
  std::optional<immersed::ImmersedSurface> surface;
  std::optional<immersed::SurfaceQuery> query;
  std::optional<immersed::ImmersedDomain> domain;
  std::optional<immersed::GhostStencilPlan> ghost_plan;
  std::optional<immersed::WallQuadraturePlan> wall_plan;
  int ghost_reach{};
  int wall_reach{};
  int halo_width{};
  immersed::LocalFlowPatternTransform transform;
  runtime::FieldRegistry registry;
  flow::FlowFieldIds fields;
};

bool finite(runtime::Real3 value) {
  return std::isfinite(value.x) && std::isfinite(value.y) &&
         std::isfinite(value.z);
}

void check_force(const flow::ForceAttemptReport &force) {
  for (const auto value :
       {force.operator_force.pressure_N, force.operator_force.viscous_N,
        force.operator_force.total_N, force.budget_reaction.pressure_N,
        force.budget_reaction.viscous_N, force.budget_reaction.total_N,
        force.surface_traction.pressure_N,
        force.surface_traction.viscous_N, force.surface_traction.total_N,
        force.consistency.pressure_N, force.consistency.viscous_N,
        force.consistency.total_N})
    HUNDUN_CHECK(finite(value));
}

bool real3_bitwise_equal(runtime::Real3 left, runtime::Real3 right) {
  return test::fp64_bits(left.x) == test::fp64_bits(right.x) &&
         test::fp64_bits(left.y) == test::fp64_bits(right.y) &&
         test::fp64_bits(left.z) == test::fp64_bits(right.z);
}

bool force_bitwise_equal(const flow::ForceAttemptReport &left,
                         const flow::ForceAttemptReport &right) {
  return real3_bitwise_equal(left.operator_force.pressure_N,
                             right.operator_force.pressure_N) &&
         real3_bitwise_equal(left.operator_force.viscous_N,
                             right.operator_force.viscous_N) &&
         real3_bitwise_equal(left.operator_force.total_N,
                             right.operator_force.total_N) &&
         real3_bitwise_equal(left.budget_reaction.pressure_N,
                             right.budget_reaction.pressure_N) &&
         real3_bitwise_equal(left.budget_reaction.viscous_N,
                             right.budget_reaction.viscous_N) &&
         real3_bitwise_equal(left.budget_reaction.total_N,
                             right.budget_reaction.total_N) &&
         real3_bitwise_equal(left.surface_traction.pressure_N,
                             right.surface_traction.pressure_N) &&
         real3_bitwise_equal(left.surface_traction.viscous_N,
                             right.surface_traction.viscous_N) &&
         real3_bitwise_equal(left.surface_traction.total_N,
                             right.surface_traction.total_N) &&
         real3_bitwise_equal(left.consistency.pressure_N,
                             right.consistency.pressure_N) &&
         real3_bitwise_equal(left.consistency.viscous_N,
                             right.consistency.viscous_N) &&
         real3_bitwise_equal(left.consistency.total_N,
                             right.consistency.total_N);
}

bool accepted_state_bitwise_equal(const test::Stage3StateSnapshot &left,
                                  const test::Stage3StateSnapshot &right) {
  return test::flow_layer_values_bitwise_equal(left.history, right.history) &&
         test::flow_layer_values_bitwise_equal(left.committed,
                                               right.committed) &&
         test::accepted_step_metadata_bitwise_equal(left.metadata,
                                                    right.metadata);
}

bool wall_gradients_bitwise_equal(
    const std::vector<flow::test::ImmersedFlowWallGradientSnapshot> &left,
    const std::vector<flow::test::ImmersedFlowWallGradientSnapshot> &right) {
  if (left.size() != right.size())
    return false;
  for (std::size_t index = 0U; index < left.size(); ++index) {
    if (left[index].link != right[index].link ||
        test::fp64_bits(left[index].normal_gradient_pa_per_m) !=
            test::fp64_bits(right[index].normal_gradient_pa_per_m))
      return false;
  }
  return true;
}

bool pressure_authority_state_bitwise_equal(
    const flow::test::ImmersedFlowPressureAuthorityStateSnapshot &left,
    const flow::test::ImmersedFlowPressureAuthorityStateSnapshot &right) {
  return left.history_available == right.history_available &&
         left.committed_available == right.committed_available &&
         wall_gradients_bitwise_equal(left.history, right.history) &&
         wall_gradients_bitwise_equal(left.committed, right.committed) &&
         wall_gradients_bitwise_equal(left.pending, right.pending);
}

bool pressure_authority_equality_oracle_is_mutation_sensitive() {
  flow::test::ImmersedFlowPressureAuthorityStateSnapshot baseline;
  baseline.history = {{11U, 1.0}, {12U, -0.0}};
  baseline.committed = {{11U, 2.0}, {12U, 3.0}};
  baseline.pending = {{11U, 4.0}, {12U, 5.0}};
  baseline.history_available = true;
  baseline.committed_available = true;
  const auto exact = baseline;
  const auto changed = [&](auto mutate) {
    auto candidate = baseline;
    mutate(candidate);
    return !pressure_authority_state_bitwise_equal(baseline, candidate);
  };
  return pressure_authority_state_bitwise_equal(baseline, exact) &&
         changed([](auto &state) { state.committed_available = false; }) &&
         changed([](auto &state) {
           state.history[1].normal_gradient_pa_per_m = 0.0;
         }) &&
         changed([](auto &state) {
           state.committed[0].normal_gradient_pa_per_m = std::nextafter(
               state.committed[0].normal_gradient_pa_per_m,
               std::numeric_limits<double>::infinity());
         }) &&
         changed([](auto &state) { state.pending[0].link = 13U; }) &&
         changed([](auto &state) { state.pending.pop_back(); });
}

class RankFailingSolver final : public linear::LinearSolver {
public:
  RankFailingSolver(const linear::LinearSolver &base, int rank,
                    int failing_rank)
      : base_(&base), rank_(rank), failing_rank_(failing_rank) {}

  linear::SolveReport
  solve(const linear::LinearOperator &linear_operator,
        linear::Preconditioner &preconditioner,
        execution::VectorView<const double> b, execution::VectorView<double> x,
        const linear::SolveControl &control) const override {
    auto report = base_->solve(linear_operator, preconditioner, b, x, control);
    if (rank_ == failing_rank_) {
      report.reason = linear::SolveTerminationReason::maximum_iterations;
      report.lowest_failing_rank = failing_rank_;
    }
    return report;
  }

private:
  const linear::LinearSolver *base_;
  int rank_;
  int failing_rank_;
};

class NthRankFailingSolver final : public linear::LinearSolver {
public:
  NthRankFailingSolver(const linear::LinearSolver &base, int rank,
                       int failing_rank, std::uint64_t failing_call)
      : base_(&base), rank_(rank), failing_rank_(failing_rank),
        failing_call_(failing_call) {
    HUNDUN_CHECK(failing_call_ > 0U);
  }

  linear::SolveReport
  solve(const linear::LinearOperator &linear_operator,
        linear::Preconditioner &preconditioner,
        execution::VectorView<const double> b, execution::VectorView<double> x,
        const linear::SolveControl &control) const override {
    auto report = base_->solve(linear_operator, preconditioner, b, x, control);
    ++call_count_;
    if (rank_ == failing_rank_ && call_count_ == failing_call_) {
      report.reason = linear::SolveTerminationReason::maximum_iterations;
      report.lowest_failing_rank = failing_rank_;
    }
    return report;
  }

private:
  const linear::LinearSolver *base_;
  int rank_;
  int failing_rank_;
  std::uint64_t failing_call_{};
  mutable std::uint64_t call_count_{};
};

void run(const std::string &mode) {
  if (mode != "transaction" && mode != "halo_contract" &&
      mode != "checkpoint_v3" &&
      mode != "checkpoint_v3_collective_preflight" &&
      mode != "checkpoint_v3_invalid_state" &&
      mode != "checkpoint_v3_path_agreement" &&
      mode != "checkpoint_v3_active_attempt" && mode != "diagnostics")
    throw runtime::Error("unknown Task 11 transaction test mode");
  HUNDUN_CHECK(test::stage3_state_equality_oracle_is_mutation_sensitive());
  HUNDUN_CHECK(pressure_authority_equality_oracle_is_mutation_sensitive());

  auto mpi = runtime::MpiContext::duplicate(MPI_COMM_WORLD);
  if (mode == "halo_contract" && mpi.size() != 1)
    throw runtime::Error(
        "Task 11 transaction Halo contract requires one MPI rank");
  Fixture fixture(mpi);
  HUNDUN_CHECK(fixture.ghost_reach == 4);
  HUNDUN_CHECK(fixture.wall_reach == 4);
  HUNDUN_CHECK(fixture.halo_width ==
               std::max(fixture.ghost_reach, fixture.wall_reach));
  auto halo = runtime::HaloExchange::create(
      fixture.decomposition,
      runtime::ExchangePlan::create(fixture.decomposition,
                                    fixture.decomposition.local_extent(),
                                    fixture.halo_width));
  HUNDUN_CHECK(halo.ghost_width() == fixture.halo_width);
  execution::CpuReferenceContext execution;
  linear::ConjugateGradientSolver momentum_solver(execution, mpi);
  linear::BiCGStabSolver pressure_solver(execution, mpi);
  linear::JacobiPreconditioner mx(execution);
  linear::JacobiPreconditioner my(execution);
  linear::JacobiPreconditioner mz(execution);
  linear::JacobiPreconditioner pressure_pc(execution);
  auto immersed_flow = flow::FixedStepImmersedFlow::create(
      fixture.decomposition, fixture.topology, fixture.geometry,
      fixture.boundaries, &*fixture.domain, &*fixture.ghost_plan,
      &*fixture.wall_plan, &fixture.transform, nullptr, mpi, execution, halo,
      momentum_solver, {&mx, &my, &mz}, pressure_solver, pressure_pc);

  const auto shared_checkpoint_path = [&](std::string name) {
    std::string text;
    if (mpi.rank() == 0)
      text = (fixture.directory.path() / std::move(name)).string();
    std::uint64_t size = text.size();
    HUNDUN_CHECK(MPI_Bcast(&size, 1, MPI_UINT64_T, 0, mpi.comm()) ==
                 MPI_SUCCESS);
    text.resize(static_cast<std::size_t>(size));
    HUNDUN_CHECK(MPI_Bcast(text.data(), static_cast<int>(text.size()),
                           MPI_BYTE, 0, mpi.comm()) == MPI_SUCCESS);
    return std::filesystem::path(text);
  };

  if (mode == "checkpoint_v3_collective_preflight") {
    HUNDUN_CHECK(mpi.size() == 2);
    auto state = fixture.make_state();
    auto config = checkpoint_v3_case(mpi.size());
    if (mpi.rank() == 1)
      config.common_flow.scalars.push_back({"unexpected", 0.0});
    const auto report = flow::write_checkpoint_v3(
        mpi, fixture.decomposition, fixture.topology, fixture.geometry,
        fixture.boundaries, config, *fixture.surface, *fixture.query,
        *fixture.domain, *fixture.ghost_plan, *fixture.wall_plan,
        fixture.transform, state, immersed_flow, {0.01, 0U},
        shared_checkpoint_path("checkpoint-v3-unsupported"));
    HUNDUN_CHECK(report.disposition() ==
                 flow::CheckpointV3Disposition::failed);
    HUNDUN_CHECK(report.reason() == flow::CheckpointV3FailureReason::presence);
    HUNDUN_CHECK(report.lowest_failing_rank() == 1);
    return;
  }

  if (mode == "checkpoint_v3_invalid_state") {
    HUNDUN_CHECK(mpi.size() == 1);
    auto invalid_physics = fixture.make_state();
    auto invalid_history = invalid_physics.snapshot(flow::FlowLayer::history);
    auto invalid_committed =
        invalid_physics.snapshot(flow::FlowLayer::committed);
    bool mutated_solid = false;
    for (mesh::LocalCellId cell = 0U;
         cell < fixture.topology.owned_cell_count(); ++cell) {
      if (fixture.domain->region(cell) != immersed::CellRegion::fluid) {
        invalid_history.density[cell] = 1.0;
        invalid_committed.density[cell] = 1.0;
        mutated_solid = true;
        break;
      }
    }
    HUNDUN_CHECK(mutated_solid);
    invalid_physics.seed_accepted_layers(invalid_history, invalid_committed);
    const auto invalid_physics_report = flow::write_checkpoint_v3(
        mpi, fixture.decomposition, fixture.topology, fixture.geometry,
        fixture.boundaries, checkpoint_v3_case(mpi.size()), *fixture.surface,
        *fixture.query, *fixture.domain, *fixture.ghost_plan,
        *fixture.wall_plan, fixture.transform, invalid_physics, immersed_flow,
        {0.01, 0U}, shared_checkpoint_path("checkpoint-v3-invalid-physics"));
    HUNDUN_CHECK(invalid_physics_report.disposition() ==
                 flow::CheckpointV3Disposition::failed);
    HUNDUN_CHECK(invalid_physics_report.reason() ==
                 flow::CheckpointV3FailureReason::state);
    auto moved_from = fixture.make_state();
    auto live = std::move(moved_from);
    static_cast<void>(live);
    const auto report = flow::write_checkpoint_v3(
        mpi, fixture.decomposition, fixture.topology, fixture.geometry,
        fixture.boundaries, checkpoint_v3_case(mpi.size()), *fixture.surface,
        *fixture.query, *fixture.domain, *fixture.ghost_plan,
        *fixture.wall_plan, fixture.transform, moved_from, immersed_flow,
        {0.01, 0U}, shared_checkpoint_path("checkpoint-v3-invalid-state"));
    HUNDUN_CHECK(report.disposition() ==
                 flow::CheckpointV3Disposition::failed);
    HUNDUN_CHECK(report.reason() == flow::CheckpointV3FailureReason::state);
    return;
  }

  if (mode == "checkpoint_v3_path_agreement") {
    HUNDUN_CHECK(mpi.size() == 2);
    auto state = fixture.make_state();
    const auto config = checkpoint_v3_case(mpi.size());
    const auto directory_a = shared_checkpoint_path("checkpoint-v3-path-a");
    const auto directory_b = shared_checkpoint_path("checkpoint-v3-path-b");
    const auto written = flow::write_checkpoint_v3(
        mpi, fixture.decomposition, fixture.topology, fixture.geometry,
        fixture.boundaries, config, *fixture.surface, *fixture.query,
        *fixture.domain, *fixture.ghost_plan, *fixture.wall_plan,
        fixture.transform, state, immersed_flow, {0.01, 0U}, directory_a);
    HUNDUN_CHECK(written.disposition() ==
                 flow::CheckpointV3Disposition::completed);
    if (mpi.rank() == 0) {
      std::filesystem::copy(directory_a, directory_b,
                            std::filesystem::copy_options::recursive);
    }
    mpi.barrier();
    const auto before = test::snapshot_stage3_state(state);
    const auto authority_before =
        flow::test::ImmersedFlowTestAccess::pressure_authority_state(
            immersed_flow);
    const auto &selected = mpi.rank() == 0 ? directory_a : directory_b;
    const auto restored = flow::read_checkpoint_v3(
        mpi, fixture.decomposition, fixture.topology, fixture.geometry,
        fixture.boundaries, config, *fixture.surface, *fixture.query,
        *fixture.domain, *fixture.ghost_plan, *fixture.wall_plan,
        fixture.transform, state, immersed_flow, selected);
    HUNDUN_CHECK(!restored.restored());
    HUNDUN_CHECK(restored.report().reason() ==
                 flow::CheckpointV3FailureReason::invalid_input);
    HUNDUN_CHECK(test::stage3_state_bitwise_equal(
        before, test::snapshot_stage3_state(state)));
    HUNDUN_CHECK(pressure_authority_state_bitwise_equal(
        authority_before,
        flow::test::ImmersedFlowTestAccess::pressure_authority_state(
            immersed_flow)));
    return;
  }

  if (mode == "checkpoint_v3_active_attempt") {
    HUNDUN_CHECK(mpi.size() == 1);
    auto state = fixture.make_state();
    const auto config = checkpoint_v3_case(mpi.size());
    const auto directory =
        shared_checkpoint_path("checkpoint-v3-active-attempt");
    const auto written = flow::write_checkpoint_v3(
        mpi, fixture.decomposition, fixture.topology, fixture.geometry,
        fixture.boundaries, config, *fixture.surface, *fixture.query,
        *fixture.domain, *fixture.ghost_plan, *fixture.wall_plan,
        fixture.transform, state, immersed_flow, {0.01, 0U}, directory);
    HUNDUN_CHECK(written.disposition() ==
                 flow::CheckpointV3Disposition::completed);
    state.begin_attempt();
    const auto before = test::snapshot_stage3_state(state);
    const auto authority_before =
        flow::test::ImmersedFlowTestAccess::pressure_authority_state(
            immersed_flow);
    const auto active_write = flow::write_checkpoint_v3(
        mpi, fixture.decomposition, fixture.topology, fixture.geometry,
        fixture.boundaries, config, *fixture.surface, *fixture.query,
        *fixture.domain, *fixture.ghost_plan, *fixture.wall_plan,
        fixture.transform, state, immersed_flow, {0.01, 0U},
        shared_checkpoint_path("checkpoint-v3-active-write"));
    HUNDUN_CHECK(active_write.disposition() ==
                 flow::CheckpointV3Disposition::failed);
    HUNDUN_CHECK(active_write.reason() ==
                 flow::CheckpointV3FailureReason::state);
    const auto restored = flow::read_checkpoint_v3(
        mpi, fixture.decomposition, fixture.topology, fixture.geometry,
        fixture.boundaries, config, *fixture.surface, *fixture.query,
        *fixture.domain, *fixture.ghost_plan, *fixture.wall_plan,
        fixture.transform, state, immersed_flow, directory);
    HUNDUN_CHECK(!restored.restored());
    HUNDUN_CHECK(restored.report().reason() ==
                 flow::CheckpointV3FailureReason::state);
    HUNDUN_CHECK(test::stage3_state_bitwise_equal(
        before, test::snapshot_stage3_state(state)));
    HUNDUN_CHECK(pressure_authority_state_bitwise_equal(
        authority_before,
        flow::test::ImmersedFlowTestAccess::pressure_authority_state(
            immersed_flow)));
    state.rollback_attempt();
    return;
  }

  if (mode == "checkpoint_v3") {
    if (mpi.size() != 1 && mpi.size() != 2 && mpi.size() != 4)
      throw runtime::Error(
          "Checkpoint v3 RED requires one, two or four MPI ranks");
    const flow::ImmersedFlowPhysics physics{
        config::DensityModel::constant, 1.0, 0.01, std::nullopt,
        std::nullopt, std::nullopt};
    const auto first_stencil = flow::make_momentum_time_stencil(
        flow::MomentumTimeOrder::backward_euler, 0.01, 0.0);
    const auto second_stencil = flow::make_momentum_time_stencil(
        flow::MomentumTimeOrder::bdf2, 0.01, 0.01);
    const std::size_t active_count =
        fixture.domain->active_cells().owned_active_count();
    std::vector<double> body_source(active_count * 3U, 0.0);
    for (std::size_t row = 0U; row < active_count; ++row)
      body_source[row * 3U] = 0.125;

    auto continuous_state = fixture.make_state();
    flow::test::ImmersedFlowTestAccess::set_manufactured_body_source(
        immersed_flow, body_source);
    const auto first = immersed_flow.attempt(
        continuous_state, physics, first_stencil, {}, {});
    flow::test::ImmersedFlowTestAccess::clear_manufactured_body_source(
        immersed_flow);
    const auto &first_base = std::get<flow::StepAttemptReport>(first.base);
    HUNDUN_CHECK(first_base.disposition ==
                 flow::StepAttemptDisposition::committed);
    HUNDUN_CHECK(first_base.pressure_corrector_count == 2U);
    HUNDUN_CHECK(first.force.has_value());

    const auto checkpoint_directory =
        shared_checkpoint_path("checkpoint-v3");
    const auto case_config = checkpoint_v3_case(mpi.size());
    const flow::CheckpointV3ControlState control{0.01, 0U};
    const auto written = flow::write_checkpoint_v3(
        mpi, fixture.decomposition, fixture.topology, fixture.geometry,
        fixture.boundaries, case_config, *fixture.surface, *fixture.query,
        *fixture.domain, *fixture.ghost_plan, *fixture.wall_plan,
        fixture.transform, continuous_state, immersed_flow, control,
        checkpoint_directory);
    HUNDUN_CHECK(written.disposition() ==
                 flow::CheckpointV3Disposition::completed);
    HUNDUN_CHECK(std::filesystem::exists(checkpoint_directory / "COMPLETED"));
    DiagnosticSink written_diagnostics;
    const diagnostics::DiagnosticRequest written_request{
        diagnostics::DiagnosticLevel::summary,
        diagnostics::DiagnosticScope::collective,
        {mpi.rank(), written.step(), written.time_s(), "completed-marker"},
        {}, 0U};
    diagnostics::collect_diagnostics(written, mpi, written_request,
                                     written_diagnostics);
    HUNDUN_CHECK(written_diagnostics.records.size() == 1U);
    HUNDUN_CHECK(written_diagnostics.records[0].status ==
                 diagnostics::DiagnosticStatus::ok);

    auto restored_state = fixture.make_state();
    auto restored_halo = runtime::HaloExchange::create(
        fixture.decomposition,
        runtime::ExchangePlan::create(fixture.decomposition,
                                      fixture.decomposition.local_extent(),
                                      fixture.halo_width));
    execution::CpuReferenceContext restored_execution;
    linear::ConjugateGradientSolver restored_momentum(restored_execution, mpi);
    linear::BiCGStabSolver restored_pressure(restored_execution, mpi);
    linear::JacobiPreconditioner restored_mx(restored_execution);
    linear::JacobiPreconditioner restored_my(restored_execution);
    linear::JacobiPreconditioner restored_mz(restored_execution);
    linear::JacobiPreconditioner restored_pressure_pc(restored_execution);
    auto restored_flow = flow::FixedStepImmersedFlow::create(
        fixture.decomposition, fixture.topology, fixture.geometry,
        fixture.boundaries, &*fixture.domain, &*fixture.ghost_plan,
        &*fixture.wall_plan, &fixture.transform, nullptr, mpi,
        restored_execution, restored_halo, restored_momentum,
        {&restored_mx, &restored_my, &restored_mz}, restored_pressure,
        restored_pressure_pc);
    const auto restored = flow::read_checkpoint_v3(
        mpi, fixture.decomposition, fixture.topology, fixture.geometry,
        fixture.boundaries, case_config, *fixture.surface, *fixture.query,
        *fixture.domain, *fixture.ghost_plan, *fixture.wall_plan,
        fixture.transform, restored_state, restored_flow,
        checkpoint_directory);
    HUNDUN_CHECK(restored.restored());
    HUNDUN_CHECK(restored.control_state().proposed_next_dt_s == 0.01);
    DiagnosticSink restored_diagnostics;
    const diagnostics::DiagnosticRequest restored_request{
        diagnostics::DiagnosticLevel::invariants,
        diagnostics::DiagnosticScope::local,
        {mpi.rank(), restored.report().step(), restored.report().time_s(),
         "restore-publish"},
        {}, 0U};
    diagnostics::collect_diagnostics(restored.report(), restored_request,
                                     restored_diagnostics);
    HUNDUN_CHECK(restored_diagnostics.records.size() == 1U);
    HUNDUN_CHECK(restored_diagnostics.records[0].status ==
                 diagnostics::DiagnosticStatus::ok);
    HUNDUN_CHECK(accepted_state_bitwise_equal(
        test::snapshot_stage3_state(continuous_state),
        test::snapshot_stage3_state(restored_state)));
    HUNDUN_CHECK(pressure_authority_state_bitwise_equal(
        flow::test::ImmersedFlowTestAccess::pressure_authority_state(
            immersed_flow),
        flow::test::ImmersedFlowTestAccess::pressure_authority_state(
            restored_flow)));

    flow::test::ImmersedFlowTestAccess::set_manufactured_body_source(
        immersed_flow, body_source);
    flow::test::ImmersedFlowTestAccess::set_manufactured_body_source(
        restored_flow, body_source);
    const auto continuous_second = immersed_flow.attempt(
        continuous_state, physics, second_stencil, {}, {});
    const auto restored_second = restored_flow.attempt(
        restored_state, physics, second_stencil, {}, {});
    flow::test::ImmersedFlowTestAccess::clear_manufactured_body_source(
        immersed_flow);
    flow::test::ImmersedFlowTestAccess::clear_manufactured_body_source(
        restored_flow);
    const auto &continuous_base =
        std::get<flow::StepAttemptReport>(continuous_second.base);
    const auto &restored_base =
        std::get<flow::StepAttemptReport>(restored_second.base);
    HUNDUN_CHECK(continuous_base.disposition ==
                 flow::StepAttemptDisposition::committed);
    HUNDUN_CHECK(restored_base.disposition ==
                 flow::StepAttemptDisposition::committed);
    HUNDUN_CHECK(continuous_base.pressure_corrector_count == 2U);
    HUNDUN_CHECK(restored_base.pressure_corrector_count == 2U);
    HUNDUN_CHECK(continuous_second.force.has_value());
    HUNDUN_CHECK(restored_second.force.has_value());
    HUNDUN_CHECK(force_bitwise_equal(*continuous_second.force,
                                     *restored_second.force));
    HUNDUN_CHECK(test::stage3_state_bitwise_equal(
        test::snapshot_stage3_state(continuous_state),
        test::snapshot_stage3_state(restored_state)));

    const auto state_before_identity_mutation =
        test::snapshot_stage3_state(restored_state);
    const auto authority_before_identity_mutation =
        flow::test::ImmersedFlowTestAccess::pressure_authority_state(
            restored_flow);
    auto incompatible_config = case_config;
    incompatible_config.common_flow.physics.dynamic_viscosity_pa_s *= 2.0;
    const auto incompatible = flow::read_checkpoint_v3(
        mpi, fixture.decomposition, fixture.topology, fixture.geometry,
        fixture.boundaries, incompatible_config, *fixture.surface,
        *fixture.query, *fixture.domain, *fixture.ghost_plan,
        *fixture.wall_plan, fixture.transform, restored_state, restored_flow,
        checkpoint_directory);
    HUNDUN_CHECK(!incompatible.restored());
    HUNDUN_CHECK(incompatible.report().reason() ==
                 flow::CheckpointV3FailureReason::fingerprint);
    HUNDUN_CHECK(test::stage3_state_bitwise_equal(
        state_before_identity_mutation,
        test::snapshot_stage3_state(restored_state)));
    HUNDUN_CHECK(pressure_authority_state_bitwise_equal(
        authority_before_identity_mutation,
        flow::test::ImmersedFlowTestAccess::pressure_authority_state(
            restored_flow)));

    const auto state_before_corruption =
        test::snapshot_stage3_state(restored_state);
    const auto authority_before_corruption =
        flow::test::ImmersedFlowTestAccess::pressure_authority_state(
            restored_flow);
    const auto manifest_path = checkpoint_directory / "manifest.v3.bin";
    if (mpi.rank() == 0) {
      std::ifstream input(manifest_path, std::ios::binary);
      std::vector<std::uint8_t> manifest_bytes{
          std::istreambuf_iterator<char>(input),
          std::istreambuf_iterator<char>()};
      HUNDUN_CHECK(!manifest_bytes.empty());
      manifest_bytes.back() ^= 0x1U;
      test::write_bytes(manifest_path, manifest_bytes);
    }
    mpi.barrier();
    const auto rejected = flow::read_checkpoint_v3(
        mpi, fixture.decomposition, fixture.topology, fixture.geometry,
        fixture.boundaries, case_config, *fixture.surface, *fixture.query,
        *fixture.domain, *fixture.ghost_plan, *fixture.wall_plan,
        fixture.transform, restored_state, restored_flow,
        checkpoint_directory);
    HUNDUN_CHECK(!rejected.restored());
    HUNDUN_CHECK(rejected.report().rollback_status() ==
                 flow::CheckpointV3CheckStatus::passed);
    HUNDUN_CHECK(test::stage3_state_bitwise_equal(
        state_before_corruption, test::snapshot_stage3_state(restored_state)));
    HUNDUN_CHECK(pressure_authority_state_bitwise_equal(
        authority_before_corruption,
        flow::test::ImmersedFlowTestAccess::pressure_authority_state(
            restored_flow)));
    return;
  }

  const auto check_halo_contract = [&] {
    constexpr int three_layer_width = 3;
    auto three_layer_halo = runtime::HaloExchange::create(
        fixture.decomposition,
        runtime::ExchangePlan::create(fixture.decomposition,
                                      fixture.decomposition.local_extent(),
                                      three_layer_width));
    linear::JacobiPreconditioner three_layer_mx(execution);
    linear::JacobiPreconditioner three_layer_my(execution);
    linear::JacobiPreconditioner three_layer_mz(execution);
    linear::JacobiPreconditioner three_layer_pressure_pc(execution);
    bool rejected = false;
    std::string observed_message;
    try {
      static_cast<void>(flow::FixedStepImmersedFlow::create(
          fixture.decomposition, fixture.topology, fixture.geometry,
          fixture.boundaries, &*fixture.domain, &*fixture.ghost_plan,
          &*fixture.wall_plan, &fixture.transform, nullptr, mpi, execution,
          three_layer_halo, momentum_solver,
          {&three_layer_mx, &three_layer_my, &three_layer_mz}, pressure_solver,
          three_layer_pressure_pc));
    } catch (const runtime::Error &error) {
      rejected = true;
      observed_message = error.what();
    }
    const std::string legacy_message =
        "immersed reconstruction Halo width is insufficient (lowest failing "
        "rank 0)";
    const std::string expected_message =
        "immersed-flow wall quadrature Halo width is insufficient (lowest failing "
        "rank 0)";
    const bool legacy_generic = observed_message == legacy_message;
    const bool wall_specific = observed_message == expected_message;
    const bool max_width_construct_success =
        halo.ghost_width() == fixture.halo_width;
    std::fprintf(
        stdout,
        "TASK11_A5_TRANSACTION ghost_reach=%d wall_reach=%d "
        "halo_width=%d max_width_construct_success=%d "
        "under_width=%d under_width_exact=%d rejected=%d observed=\"%s\" "
        "legacy_generic=%d future_wall_specific=%d\n",
        fixture.ghost_reach, fixture.wall_reach, fixture.halo_width,
        max_width_construct_success ? 1 : 0, three_layer_width,
        three_layer_halo.ghost_width() == three_layer_width ? 1 : 0,
        rejected ? 1 : 0, observed_message.c_str(), legacy_generic ? 1 : 0,
        wall_specific ? 1 : 0);
    std::fflush(stdout);
    HUNDUN_CHECK(max_width_construct_success);
    HUNDUN_CHECK(three_layer_halo.ghost_width() == three_layer_width);
    HUNDUN_CHECK(rejected);
    HUNDUN_CHECK(legacy_generic || wall_specific);
    HUNDUN_CHECK(wall_specific);
  };
  if (mpi.size() == 1)
    check_halo_contract();
  if (mode == "halo_contract")
    return;

  const flow::ImmersedFlowPhysics physics{config::DensityModel::constant,
                                    1.0,
                                    0.01,
                                    std::nullopt,
                                    std::nullopt,
                                    std::nullopt};
  const auto stencil = flow::make_momentum_time_stencil(
      flow::MomentumTimeOrder::backward_euler, 0.01, 0.0);

  auto no_source_state = fixture.make_state();
  const auto no_source =
      immersed_flow.attempt(no_source_state, physics, stencil, {}, {});
  const auto &no_source_base =
      std::get<flow::StepAttemptReport>(no_source.base);
  HUNDUN_CHECK(no_source_base.disposition ==
               flow::StepAttemptDisposition::committed);
  HUNDUN_CHECK(no_source.force.has_value());
  check_force(*no_source.force);

  const std::size_t active_count =
      fixture.domain->active_cells().owned_active_count();
  std::vector<double> body_source(active_count * 3U, 0.0);
  for (std::size_t row = 0U; row < active_count; ++row)
    body_source[row * 3U] = 0.125;
  flow::test::ImmersedFlowTestAccess::set_manufactured_body_source(immersed_flow,
                                                                 body_source);
  auto source_state = fixture.make_state();
  const auto with_source =
      immersed_flow.attempt(source_state, physics, stencil, {}, {});
  flow::test::ImmersedFlowTestAccess::clear_manufactured_body_source(immersed_flow);
  const auto &source_base = std::get<flow::StepAttemptReport>(with_source.base);
  HUNDUN_CHECK(source_base.disposition ==
               flow::StepAttemptDisposition::committed);
  HUNDUN_CHECK(with_source.force.has_value());
  check_force(*with_source.force);
  HUNDUN_CHECK(flow::test::ImmersedFlowTestAccess::
                   last_wall_pressure_gradient_application_norm(immersed_flow) > 0.0);
  HUNDUN_CHECK(!test::flow_layer_values_bitwise_equal(
      no_source_state.snapshot(flow::FlowLayer::committed),
      source_state.snapshot(flow::FlowLayer::committed)));
  const auto accepted_authority =
      flow::test::ImmersedFlowTestAccess::pressure_authority_state(immersed_flow);
  HUNDUN_CHECK(accepted_authority.committed_available);
  HUNDUN_CHECK(!accepted_authority.committed.empty());

  if (mode == "diagnostics") {
    auto source = immersed_flow.diagnostic_source(source_state, with_source);
    auto mutated_report = with_source;
    std::get<flow::StepAttemptReport>(mutated_report.base)
        .pressure_corrector_count = 1U;
    bool mutated_report_rejected = false;
    try {
      static_cast<void>(
          immersed_flow.diagnostic_source(source_state, mutated_report));
    } catch (const runtime::Error &) {
      mutated_report_rejected = true;
    }
    HUNDUN_CHECK(mutated_report_rejected);
    const auto metadata = source_state.metadata();
    const auto request = [&](diagnostics::DiagnosticLevel level,
                             diagnostics::DiagnosticScope scope) {
      return diagnostics::DiagnosticRequest{
          level, scope,
          {mpi.rank(), metadata.step, metadata.time_s,
           "immersed-flow.attempt-result"},
          {}, 0U};
    };
    DiagnosticSink local_sink;
    diagnostics::collect_diagnostics(
        source, request(diagnostics::DiagnosticLevel::summary,
                        diagnostics::DiagnosticScope::local),
        local_sink);
    HUNDUN_CHECK(local_sink.records.size() == 1U);
    HUNDUN_CHECK(local_sink.records[0].status ==
                 diagnostics::DiagnosticStatus::ok);
    HUNDUN_CHECK(has_metric(local_sink.records[0],
                            "continuity.final-normalized-l2"));
    HUNDUN_CHECK(has_metric(local_sink.records[0],
                            "wall-penetration.maximum"));
    HUNDUN_CHECK(
        has_metric(local_sink.records[0], "force.operator.total.x"));
    HUNDUN_CHECK(has_metric(local_sink.records[0],
                            "force.budget-reaction.total.x"));
    HUNDUN_CHECK(has_metric(local_sink.records[0],
                            "force.surface-traction.total.x"));
    HUNDUN_CHECK(has_metric(local_sink.records[0],
                            "force.consistency.total.x"));

    DiagnosticSink counter_sink;
    diagnostics::collect_diagnostics(
        source, request(diagnostics::DiagnosticLevel::counters,
                        diagnostics::DiagnosticScope::local),
        counter_sink);
    HUNDUN_CHECK(counter_sink.records.size() == 1U);
    HUNDUN_CHECK(counter_sink.records[0].counters.size() >= 16U);

    DiagnosticSink collective_sink;
    diagnostics::collect_diagnostics(
        source, mpi,
        request(diagnostics::DiagnosticLevel::summary,
                diagnostics::DiagnosticScope::collective),
        collective_sink);
    HUNDUN_CHECK(collective_sink.records.size() == 1U);

    auto unsupported_physics = physics;
    unsupported_physics.density_model = config::DensityModel::material;
    const auto unsupported = immersed_flow.attempt(
        source_state, unsupported_physics, stencil, {}, {});
    HUNDUN_CHECK(
        std::get<flow::StepAttemptReport>(unsupported.base).disposition ==
        flow::StepAttemptDisposition::non_retryable_failure);
    bool preflight_stale_rejected = false;
    try {
      static_cast<void>(source.report());
    } catch (const runtime::Error &) {
      preflight_stale_rejected = true;
    }
    HUNDUN_CHECK(preflight_stale_rejected);
    auto preflight_source =
        immersed_flow.diagnostic_source(source_state, unsupported);
    DiagnosticSink preflight_sink;
    diagnostics::collect_diagnostics(
        preflight_source,
        request(diagnostics::DiagnosticLevel::summary,
                diagnostics::DiagnosticScope::local),
        preflight_sink);
    HUNDUN_CHECK(preflight_sink.records.size() == 1U);
    HUNDUN_CHECK(preflight_sink.records[0].status ==
                 diagnostics::DiagnosticStatus::failed);

    flow::test::ImmersedFlowTestAccess::set_failure_stage(
        flow::test::ImmersedFlowAttemptFailureStage::after_corrector_1);
    const auto failed =
        immersed_flow.attempt(source_state, physics, stencil, {}, {});
    flow::test::ImmersedFlowTestAccess::clear_failure_stage();
    HUNDUN_CHECK(std::get<flow::StepAttemptReport>(failed.base).disposition ==
                 flow::StepAttemptDisposition::recoverable_failure);
    bool stale_rejected = false;
    try {
      static_cast<void>(preflight_source.report());
    } catch (const runtime::Error &) {
      stale_rejected = true;
    }
    HUNDUN_CHECK(stale_rejected);

    auto failed_source = immersed_flow.diagnostic_source(source_state, failed);
    DiagnosticSink failed_sink;
    diagnostics::collect_diagnostics(
        failed_source,
        request(diagnostics::DiagnosticLevel::summary,
                diagnostics::DiagnosticScope::local),
        failed_sink);
    HUNDUN_CHECK(failed_sink.records.size() == 1U);
    HUNDUN_CHECK(failed_sink.records[0].status ==
                 diagnostics::DiagnosticStatus::failed);
    HUNDUN_CHECK(!has_metric(failed_sink.records[0],
                             "force.operator.total.x"));
    return;
  }

  const int failing_rank = mpi.size() > 1 ? 1 : 0;
  const auto check_failure =
      [&](flow::test::ImmersedFlowAttemptFailureStage failure_stage,
          flow::StepFailureReason reason, std::uint32_t correctors) {
        const auto before = test::snapshot_stage3_state(source_state);
        const auto authority_before =
            flow::test::ImmersedFlowTestAccess::pressure_authority_state(immersed_flow);
        flow::test::ImmersedFlowTestAccess::set_failure_stage(failure_stage);
        const auto failed =
            immersed_flow.attempt(source_state, physics, stencil, {}, {});
        flow::test::ImmersedFlowTestAccess::clear_failure_stage();
        const auto &failed_base =
            std::get<flow::StepAttemptReport>(failed.base);
        HUNDUN_CHECK(failed_base.disposition ==
                     flow::StepAttemptDisposition::recoverable_failure);
        HUNDUN_CHECK(failed_base.reason == reason);
        HUNDUN_CHECK(failed_base.lowest_failing_rank == failing_rank);
        HUNDUN_CHECK(failed_base.pressure_corrector_count == correctors);
        HUNDUN_CHECK(!failed.force.has_value());
        HUNDUN_CHECK(test::stage3_state_bitwise_equal(
            before, test::snapshot_stage3_state(source_state)));
        HUNDUN_CHECK(pressure_authority_state_bitwise_equal(
            authority_before,
            flow::test::ImmersedFlowTestAccess::pressure_authority_state(immersed_flow)));
      };
  check_failure(flow::test::ImmersedFlowAttemptFailureStage::after_corrector_1,
                flow::StepFailureReason::non_finite_trial, 1U);
  check_failure(flow::test::ImmersedFlowAttemptFailureStage::after_corrector_2,
                flow::StepFailureReason::non_finite_trial, 2U);
  check_failure(flow::test::ImmersedFlowAttemptFailureStage::final_wall_penetration,
                flow::StepFailureReason::final_continuity_residual, 2U);
  check_failure(
      flow::test::ImmersedFlowAttemptFailureStage::final_force_reconstruction,
      flow::StepFailureReason::final_conservation_defect, 2U);

  const auto before_commit_state = test::snapshot_stage3_state(source_state);
  const auto before_commit_authority =
      flow::test::ImmersedFlowTestAccess::pressure_authority_state(immersed_flow);
  flow::test::ImmersedFlowTestAccess::set_failure_stage(
      flow::test::ImmersedFlowAttemptFailureStage::before_commit, failing_rank);
  const auto commit_preparation_failed =
      immersed_flow.attempt(source_state, physics, stencil, {}, {});
  flow::test::ImmersedFlowTestAccess::clear_failure_stage();
  const auto &commit_preparation_base =
      std::get<flow::StepAttemptReport>(commit_preparation_failed.base);
  const bool local_commit_contract =
      commit_preparation_base.disposition ==
          flow::StepAttemptDisposition::non_retryable_failure &&
      commit_preparation_base.reason == flow::StepFailureReason::invalid_input &&
      commit_preparation_base.lowest_failing_rank == failing_rank &&
      commit_preparation_base.pressure_corrector_count == 2U &&
      !commit_preparation_failed.force.has_value() &&
      test::stage3_state_bitwise_equal(
          before_commit_state, test::snapshot_stage3_state(source_state)) &&
      pressure_authority_state_bitwise_equal(
          before_commit_authority,
          flow::test::ImmersedFlowTestAccess::pressure_authority_state(immersed_flow));
  int local_commit_contract_value = local_commit_contract ? 1 : 0;
  int collective_commit_contract_value = 0;
  HUNDUN_CHECK(MPI_Allreduce(&local_commit_contract_value,
                             &collective_commit_contract_value, 1, MPI_INT,
                             MPI_MIN, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(collective_commit_contract_value == 1);

  const auto check_wall_input_failure =
      [&](flow::test::ImmersedFlowWallInputFailure failure) {
        const auto before = test::snapshot_stage3_state(source_state);
        const auto authority_before =
            flow::test::ImmersedFlowTestAccess::pressure_authority_state(immersed_flow);
        flow::test::ImmersedFlowTestAccess::set_wall_input_failure(failure,
                                                                 failing_rank);
        const auto failed =
            immersed_flow.attempt(source_state, physics, stencil, {}, {});
        flow::test::ImmersedFlowTestAccess::clear_wall_input_failure();
        const auto &failed_base =
            std::get<flow::StepAttemptReport>(failed.base);
        HUNDUN_CHECK(failed_base.disposition ==
                     flow::StepAttemptDisposition::recoverable_failure);
        HUNDUN_CHECK(failed_base.reason ==
                     flow::StepFailureReason::non_finite_trial);
        HUNDUN_CHECK(failed_base.lowest_failing_rank == failing_rank);
        HUNDUN_CHECK(failed_base.pressure_corrector_count == 0U);
        HUNDUN_CHECK(!failed.force.has_value());
        HUNDUN_CHECK(test::stage3_state_bitwise_equal(
            before, test::snapshot_stage3_state(source_state)));
        HUNDUN_CHECK(pressure_authority_state_bitwise_equal(
            authority_before,
            flow::test::ImmersedFlowTestAccess::pressure_authority_state(immersed_flow)));
      };
  check_wall_input_failure(
      flow::test::ImmersedFlowWallInputFailure::non_positive_density);
  check_wall_input_failure(
      flow::test::ImmersedFlowWallInputFailure::non_positive_coefficient);

  const auto authority_before_rotation =
      flow::test::ImmersedFlowTestAccess::pressure_authority_state(immersed_flow);
  const auto accepted_rotation =
      immersed_flow.attempt(source_state, physics, stencil, {}, {});
  const auto &accepted_rotation_base =
      std::get<flow::StepAttemptReport>(accepted_rotation.base);
  HUNDUN_CHECK(accepted_rotation_base.disposition ==
               flow::StepAttemptDisposition::committed);
  const auto authority_after_rotation =
      flow::test::ImmersedFlowTestAccess::pressure_authority_state(immersed_flow);
  HUNDUN_CHECK(authority_after_rotation.history_available ==
               authority_before_rotation.committed_available);
  HUNDUN_CHECK(authority_after_rotation.committed_available);
  HUNDUN_CHECK(wall_gradients_bitwise_equal(
      authority_after_rotation.history,
      authority_before_rotation.committed));
  HUNDUN_CHECK(wall_gradients_bitwise_equal(
      authority_after_rotation.pending, authority_before_rotation.history));

  RankFailingSolver failing_momentum(momentum_solver, mpi.rank(), failing_rank);
  linear::JacobiPreconditioner fmx(execution);
  linear::JacobiPreconditioner fmy(execution);
  linear::JacobiPreconditioner fmz(execution);
  linear::JacobiPreconditioner fpressure_pc(execution);
  auto failing_immersed_flow = flow::FixedStepImmersedFlow::create(
      fixture.decomposition, fixture.topology, fixture.geometry,
      fixture.boundaries, &*fixture.domain, &*fixture.ghost_plan,
      &*fixture.wall_plan, &fixture.transform, nullptr, mpi, execution, halo,
      failing_momentum, {&fmx, &fmy, &fmz}, pressure_solver, fpressure_pc);
  auto momentum_state = fixture.make_state();
  const auto momentum_before = test::snapshot_stage3_state(momentum_state);
  const auto momentum_authority_before =
      flow::test::ImmersedFlowTestAccess::pressure_authority_state(failing_immersed_flow);
  const auto momentum_failed =
      failing_immersed_flow.attempt(momentum_state, physics, stencil, {}, {});
  const auto &momentum_base =
      std::get<flow::StepAttemptReport>(momentum_failed.base);
  HUNDUN_CHECK(momentum_base.disposition ==
               flow::StepAttemptDisposition::recoverable_failure);
  HUNDUN_CHECK(momentum_base.reason ==
               flow::StepFailureReason::momentum_linear_solve);
  HUNDUN_CHECK(momentum_base.lowest_failing_rank == failing_rank);
  HUNDUN_CHECK(momentum_base.pressure_corrector_count == 0U);
  HUNDUN_CHECK(!momentum_failed.force.has_value());
  HUNDUN_CHECK(test::stage3_state_bitwise_equal(
      momentum_before, test::snapshot_stage3_state(momentum_state)));
  HUNDUN_CHECK(pressure_authority_state_bitwise_equal(
      momentum_authority_before,
      flow::test::ImmersedFlowTestAccess::pressure_authority_state(
          failing_immersed_flow)));

  NthRankFailingSolver failing_inner_momentum(momentum_solver, mpi.rank(),
                                               failing_rank, 4U);
  linear::JacobiPreconditioner imx(execution);
  linear::JacobiPreconditioner imy(execution);
  linear::JacobiPreconditioner imz(execution);
  linear::JacobiPreconditioner inner_pressure_pc(execution);
  auto inner_failing_immersed_flow = flow::FixedStepImmersedFlow::create(
      fixture.decomposition, fixture.topology, fixture.geometry,
      fixture.boundaries, &*fixture.domain, &*fixture.ghost_plan,
      &*fixture.wall_plan, &fixture.transform, nullptr, mpi, execution, halo,
      failing_inner_momentum, {&imx, &imy, &imz}, pressure_solver,
      inner_pressure_pc);
  auto inner_state = fixture.make_state();
  const auto inner_before = test::snapshot_stage3_state(inner_state);
  const auto inner_authority_before =
      flow::test::ImmersedFlowTestAccess::pressure_authority_state(
          inner_failing_immersed_flow);
  const auto inner_failed =
      inner_failing_immersed_flow.attempt(inner_state, physics, stencil, {}, {});
  const auto &inner_base =
      std::get<flow::StepAttemptReport>(inner_failed.base);
  HUNDUN_CHECK(inner_base.disposition ==
               flow::StepAttemptDisposition::recoverable_failure);
  HUNDUN_CHECK(inner_base.reason ==
               flow::StepFailureReason::pressure_linear_solve);
  HUNDUN_CHECK(inner_base.lowest_failing_rank == failing_rank);
  HUNDUN_CHECK(inner_base.pressure_corrector_count == 0U);
  HUNDUN_CHECK(!inner_failed.force.has_value());
  HUNDUN_CHECK(test::stage3_state_bitwise_equal(
      inner_before, test::snapshot_stage3_state(inner_state)));
  HUNDUN_CHECK(pressure_authority_state_bitwise_equal(
      inner_authority_before,
      flow::test::ImmersedFlowTestAccess::pressure_authority_state(
          inner_failing_immersed_flow)));

  const auto check_pressure_solve_failure =
      [&](std::uint64_t failing_call, std::uint32_t expected_correctors) {
        NthRankFailingSolver failing_pressure(pressure_solver, mpi.rank(),
                                              failing_rank, failing_call);
        linear::JacobiPreconditioner pmx(execution);
        linear::JacobiPreconditioner pmy(execution);
        linear::JacobiPreconditioner pmz(execution);
        linear::JacobiPreconditioner ppressure_pc(execution);
        auto pressure_failing_immersed_flow = flow::FixedStepImmersedFlow::create(
            fixture.decomposition, fixture.topology, fixture.geometry,
            fixture.boundaries, &*fixture.domain, &*fixture.ghost_plan,
            &*fixture.wall_plan, &fixture.transform, nullptr, mpi, execution,
            halo, momentum_solver, {&pmx, &pmy, &pmz}, failing_pressure,
            ppressure_pc);
        auto rollback_state = fixture.make_state();
        const auto before = test::snapshot_stage3_state(rollback_state);
        const auto authority_before =
            flow::test::ImmersedFlowTestAccess::pressure_authority_state(
                pressure_failing_immersed_flow);
        const auto failed = pressure_failing_immersed_flow.attempt(
            rollback_state, physics, stencil, {}, {});
        const auto &failed_base =
            std::get<flow::StepAttemptReport>(failed.base);
        HUNDUN_CHECK(failed_base.disposition ==
                     flow::StepAttemptDisposition::recoverable_failure);
        HUNDUN_CHECK(failed_base.reason ==
                     flow::StepFailureReason::pressure_linear_solve);
        HUNDUN_CHECK(failed_base.lowest_failing_rank == failing_rank);
        HUNDUN_CHECK(failed_base.pressure_corrector_count ==
                     expected_correctors);
        HUNDUN_CHECK(!failed.force.has_value());
        HUNDUN_CHECK(test::stage3_state_bitwise_equal(
            before, test::snapshot_stage3_state(rollback_state)));
        HUNDUN_CHECK(pressure_authority_state_bitwise_equal(
            authority_before,
            flow::test::ImmersedFlowTestAccess::pressure_authority_state(
                pressure_failing_immersed_flow)));
      };
  check_pressure_solve_failure(1U, 0U);
  check_pressure_solve_failure(2U, 1U);
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  return hundun::test::run([&] {
    const std::string mode = argc > 1 ? argv[1] : "transaction";
    run(mode);
  });
}
