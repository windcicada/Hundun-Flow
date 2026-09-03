// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/bc_basic_boundary.hpp"
#include "hundun/cfg_resolved_case_v3.hpp"
#include "hundun/diag_immersed_module.hpp"
#include "hundun/diag_immersed_static.hpp"
#include "hundun/exec_execution.hpp"
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
#include "hundun/rt_mpi_environment.hpp"
#include "hundun/rt_structured_decomposition.hpp"
#include "tests/support/flow_immersed_test_access.hpp"
#include "tests/support/stage3_stl_fixture.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace {

using namespace hundun;

class RecordingSink final : public diagnostics::DiagnosticSink {
public:
  void submit(const diagnostics::DiagnosticRecord &record) override {
    records.push_back(record);
  }
  std::vector<diagnostics::DiagnosticRecord> records;
};

std::uint64_t fp64_bits(double value) noexcept {
  std::uint64_t result{};
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

runtime::Int3 process_grid(int ranks) {
  if (ranks == 1)
    return {1, 1, 1};
  if (ranks == 2)
    return {2, 1, 1};
  throw runtime::Error("unsupported Stage 3 diagnostics rank count");
}

config::FlowCaseConfig periodic_case(int ranks) {
  config::FlowCaseConfig result{};
  result.schema_version = 2;
  result.simulation_type = config::SimulationType::variable_density_flow;
  result.density_model = config::DensityModel::constant;
  result.resources.expected_ranks = ranks;
  result.resources.process_grid = process_grid(ranks);
  result.mesh.cells = {12, 12, 12};
  result.mesh.origin_m = {0.0, 0.0, 0.0};
  result.mesh.length_m = {1.0, 1.0, 1.0};
  result.physics.rho_ref_kg_per_m3 = 1.0;
  result.physics.dynamic_viscosity_pa_s = 0.01;
  result.physics.inlet_consistency_rtol = 1.0e-12;
  constexpr std::array<config::PatchName, 6> names{
      config::PatchName::x_min, config::PatchName::x_max,
      config::PatchName::y_min, config::PatchName::y_max,
      config::PatchName::z_min, config::PatchName::z_max};
  for (std::size_t index = 0U; index < names.size(); ++index) {
    result.boundaries[index].patch = names[index];
    result.boundaries[index].type = config::BoundaryType::periodic;
  }
  return result;
}

runtime::FieldDescriptor cell_field(const char *name,
                                    std::uint32_t components,
                                    int ghost_width) {
  return {name,
          "1",
          "stage3_diagnostics",
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
          "stage3_diagnostics",
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
    return runtime::Real3{0.5 * (left.x + right.x),
                          0.5 * (left.y + right.y),
                          0.5 * (left.z + right.z)};
  };
  std::vector<test::StlFixtureTriangle> result;
  result.reserve(4U * coarse.size());
  for (const auto &triangle : coarse) {
    const auto ab = midpoint(triangle.vertices[0], triangle.vertices[1]);
    const auto bc = midpoint(triangle.vertices[1], triangle.vertices[2]);
    const auto ca = midpoint(triangle.vertices[2], triangle.vertices[0]);
    result.push_back({triangle.file_normal, {triangle.vertices[0], ab, ca}});
    result.push_back({triangle.file_normal, {ab, triangle.vertices[1], bc}});
    result.push_back({triangle.file_normal, {ca, bc, triangle.vertices[2]}});
    result.push_back({triangle.file_normal, {ab, bc, ca}});
  }
  return result;
}

struct Fixture final {
  explicit Fixture(const runtime::MpiContext &mpi)
      : directory("stage3-diagnostics"),
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
    surface.emplace(immersed::ImmersedSurface::load_collective(
        std::filesystem::path(path_text), 1.0, mpi, 0));
    query.emplace(immersed::SurfaceQuery::create(*surface));
    domain.emplace(immersed::ImmersedDomain::create(
        *surface, *query, config::ImmersedFluidSide::outside, topology,
        geometry, boundaries, mpi));
    ghost.emplace(immersed::GhostStencilPlan::create(
        *surface, *query, *domain, topology, geometry, decomposition, mpi));
    wall.emplace(immersed::WallQuadraturePlan::create(
        *surface, *query, *domain, topology, geometry, mpi));
    halo_width = std::max(static_cast<int>(ghost->maximum_halo_reach()),
                          static_cast<int>(wall->maximum_halo_reach()));
    fields.density = registry.declare_field(cell_field("rho", 1U, halo_width));
    fields.velocity =
        registry.declare_field(cell_field("velocity", 3U, halo_width));
    fields.mechanical_pressure =
        registry.declare_field(cell_field("pi", 1U, halo_width));
    fields.face_velocity =
        registry.declare_field(face_field("face_velocity", 3U));
    fields.face_mass_flux = finite_volume::declare_face_mass_flux(registry);
    registry.freeze();
  }

  flow::FlowState make_state() const {
    auto state = flow::FlowState::create(
        registry,
        {decomposition.local_extent(), topology.local_face_count()}, fields,
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
  std::optional<immersed::GhostStencilPlan> ghost;
  std::optional<immersed::WallQuadraturePlan> wall;
  immersed::LocalFlowPatternTransform transform;
  int halo_width{};
  runtime::FieldRegistry registry;
  flow::FlowFieldIds fields;
};

diagnostics::DiagnosticRequest request(int rank,
                                       diagnostics::DiagnosticScope scope) {
  diagnostics::DiagnosticRequest result;
  result.level = diagnostics::DiagnosticLevel::summary;
  result.scope = scope;
  result.frame = {rank, 1U, 0.01, "immersed-flow.wall-force"};
  return result;
}

void check_force_fields(const diagnostics::DiagnosticRecord &record) {
  HUNDUN_CHECK(record.module_kind ==
               diagnostics::DiagnosticModuleKind::wall_force);
  HUNDUN_CHECK(record.module_id == "hundun.immersed.wall-force");
  HUNDUN_CHECK(record.metrics.size() == 48U);
  constexpr std::array<std::string_view, 12> force_prefixes{
      "force.budget-reaction.pressure", "force.budget-reaction.total",
      "force.budget-reaction.viscous", "force.consistency.pressure",
      "force.consistency.total", "force.consistency.viscous",
      "force.operator.pressure", "force.operator.total",
      "force.operator.viscous", "force.surface-traction.pressure",
      "force.surface-traction.total", "force.surface-traction.viscous"};
  std::size_t position = 0U;
  for (const auto prefix : force_prefixes)
    for (const auto axis : {"x", "y", "z"}) {
      HUNDUN_CHECK(record.metrics[position].id ==
                   std::string(prefix) + "." + axis);
      HUNDUN_CHECK(record.metrics[position].unit == "N");
      ++position;
    }
  for (const auto part : {"pressure", "total", "viscous"})
    for (const auto axis : {"x", "y", "z"}) {
      HUNDUN_CHECK(record.metrics[position].id ==
                   std::string("moment.") + part + "." + axis);
      HUNDUN_CHECK(record.metrics[position].unit == "N*m");
      ++position;
    }
  for (const auto axis : {"x", "y", "z"}) {
    HUNDUN_CHECK(record.metrics[position].id ==
                 std::string("surface-area-vector-closure.") + axis);
    HUNDUN_CHECK(record.metrics[position].unit == "m2");
    ++position;
  }
  HUNDUN_CHECK(position == record.metrics.size());
}

void check_string_agreement(const runtime::MpiContext &mpi,
                            const std::string &value) {
  std::uint64_t size = value.size();
  std::uint64_t minimum = size;
  std::uint64_t maximum = size;
  HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, &minimum, 1, MPI_UINT64_T, MPI_MIN,
                             mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, &maximum, 1, MPI_UINT64_T, MPI_MAX,
                             mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(minimum == maximum);
  std::string root = value;
  HUNDUN_CHECK(MPI_Bcast(root.data(), static_cast<int>(root.size()), MPI_BYTE,
                         0, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(root == value);
}

void run() {
  auto mpi = runtime::MpiContext::duplicate(MPI_COMM_WORLD);
  HUNDUN_CHECK(mpi.size() == 1 || mpi.size() == 2);
  Fixture fixture(mpi);
  const auto static_summary = diagnostics::summarize_immersed_static(
      *fixture.surface, *fixture.query, *fixture.domain, *fixture.ghost,
      *fixture.wall, fixture.transform);
  HUNDUN_CHECK(static_summary.vertex_count == fixture.surface->vertex_count());
  HUNDUN_CHECK(static_summary.triangle_count ==
               fixture.surface->triangle_count());
  HUNDUN_CHECK(static_summary.immersed_link_count ==
               fixture.domain->links().size());

  auto halo = runtime::HaloExchange::create(
      fixture.decomposition,
      runtime::ExchangePlan::create(fixture.decomposition,
                                    fixture.decomposition.local_extent(),
                                    fixture.halo_width));
  execution::CpuReferenceContext execution;
  linear::ConjugateGradientSolver momentum_solver(execution, mpi);
  linear::BiCGStabSolver pressure_solver(execution, mpi);
  linear::JacobiPreconditioner mx(execution);
  linear::JacobiPreconditioner my(execution);
  linear::JacobiPreconditioner mz(execution);
  linear::JacobiPreconditioner pressure_pc(execution);
  auto facade = flow::FixedStepImmersedFlow::create(
      fixture.decomposition, fixture.topology, fixture.geometry,
      fixture.boundaries, &*fixture.domain, &*fixture.ghost, &*fixture.wall,
      &fixture.transform, nullptr, mpi, execution, halo, momentum_solver,
      {&mx, &my, &mz}, pressure_solver, pressure_pc);
  auto state = fixture.make_state();
  const flow::ImmersedFlowPhysics physics{config::DensityModel::constant,
                                          1.0,
                                          0.01,
                                          std::nullopt,
                                          std::nullopt,
                                          std::nullopt};
  const auto stencil = flow::make_momentum_time_stencil(
      flow::MomentumTimeOrder::backward_euler, 0.01, 0.0);
  const auto report = facade.attempt(state, physics, stencil, {}, {});
  const auto &base = std::get<flow::StepAttemptReport>(report.base);
  HUNDUN_CHECK(base.disposition == flow::StepAttemptDisposition::committed);
  HUNDUN_CHECK(base.pressure_corrector_count == 2U);
  HUNDUN_CHECK(report.force.has_value());
  auto source = facade.diagnostic_source(state, report);
  HUNDUN_CHECK(source.local_flow_pattern_available());
  const auto attempt_summary = diagnostics::with_local_flow_pattern_snapshot(
      static_summary, source);
  HUNDUN_CHECK(attempt_summary.local_flow_pattern_snapshot_available);
  HUNDUN_CHECK(attempt_summary.local_flow_pattern_algorithm_fingerprint ==
               source.local_flow_pattern_algorithm_fingerprint());
  HUNDUN_CHECK(attempt_summary.local_flow_pattern_row_fingerprint ==
               source.local_flow_pattern_row_fingerprint());
  HUNDUN_CHECK(attempt_summary.replacement_group_count ==
               source.local_flow_pattern_replacement_group_count());
  HUNDUN_CHECK(attempt_summary.algebraic_occurrence_count ==
               source.local_flow_pattern_algebraic_occurrence_count());
  auto lfp_request = diagnostics::DiagnosticRequest{};
  lfp_request.level = diagnostics::DiagnosticLevel::summary;
  lfp_request.scope = diagnostics::DiagnosticScope::local;
  lfp_request.frame = {mpi.rank(), 0U, 0.0,
                       "immersed-static.local-flow-pattern"};
  RecordingSink lfp_sink;
  diagnostics::collect_diagnostics(
      attempt_summary, diagnostics::DiagnosticModuleKind::local_flow_pattern,
      lfp_request, lfp_sink);
  HUNDUN_CHECK(lfp_sink.records.size() == 1U);
  HUNDUN_CHECK(lfp_sink.records.front().metrics.size() == 1U);
  HUNDUN_CHECK(lfp_sink.records.front().metrics.front().value.bits ==
               fp64_bits(
                   source.local_flow_pattern_replacement_coefficient_l2()));
  for (const auto &entry : {
           std::pair{diagnostics::DiagnosticModuleKind::immersed_surface,
                     std::string_view("immersed-static.surface")},
           std::pair{diagnostics::DiagnosticModuleKind::ghost_stencil,
                     std::string_view("immersed-static.ghost-stencil")},
           std::pair{diagnostics::DiagnosticModuleKind::local_flow_pattern,
                     std::string_view("immersed-static.local-flow-pattern")},
       }) {
    auto collective_request = diagnostics::DiagnosticRequest{};
    collective_request.level = diagnostics::DiagnosticLevel::summary;
    collective_request.scope = diagnostics::DiagnosticScope::collective;
    collective_request.frame = {mpi.rank(), 0U, 0.0, entry.second};
    RecordingSink static_collective;
    diagnostics::collect_diagnostics(
        entry.first == diagnostics::DiagnosticModuleKind::local_flow_pattern
            ? attempt_summary
            : static_summary,
        entry.first, mpi, collective_request, static_collective);
    HUNDUN_CHECK(static_collective.records.size() == 1U);
    check_string_agreement(
        mpi, static_collective.records.front().state_fingerprint.hex);
    auto comparable = static_collective.records.front();
    comparable.rank = 0;
    check_string_agreement(mpi, diagnostics::to_canonical_json(comparable));
  }
  const auto descriptor = diagnostics::describe_diagnostics(
      source, diagnostics::DiagnosticModuleKind::wall_force);
  HUNDUN_CHECK(descriptor.module_kind ==
               diagnostics::DiagnosticModuleKind::wall_force);
  const std::vector<std::string_view> expected_fields{
      "area-closure", "density-variant", "force.budget-reaction",
      "force.consistency", "force.operator", "force.surface-traction",
      "lowest-rank", "moment", "point-count", "wale"};
  HUNDUN_CHECK(diagnostics::diagnostic_fingerprint_field_ids(
                   source, diagnostics::DiagnosticModuleKind::wall_force) ==
               expected_fields);

  RecordingSink first;
  diagnostics::collect_diagnostics(
      source, diagnostics::DiagnosticModuleKind::wall_force,
      request(mpi.rank(), diagnostics::DiagnosticScope::local), first);
  HUNDUN_CHECK(first.records.size() == 1U);
  check_force_fields(first.records.front());
  RecordingSink repeated;
  diagnostics::collect_diagnostics(
      source, diagnostics::DiagnosticModuleKind::wall_force,
      request(mpi.rank(), diagnostics::DiagnosticScope::local), repeated);
  HUNDUN_CHECK(diagnostics::to_canonical_json(first.records.front()) ==
               diagnostics::to_canonical_json(repeated.records.front()));

  RecordingSink collective;
  diagnostics::collect_diagnostics(
      source, diagnostics::DiagnosticModuleKind::wall_force, mpi,
      request(mpi.rank(), diagnostics::DiagnosticScope::collective),
      collective);
  HUNDUN_CHECK(collective.records.size() == 1U);
  check_force_fields(collective.records.front());
  auto comparable_force = collective.records.front();
  comparable_force.rank = 0;
  check_string_agreement(
      mpi, diagnostics::to_canonical_json(comparable_force));

  flow::test::ImmersedFlowTestAccess::set_failure_stage(
      flow::test::ImmersedFlowAttemptFailureStage::before_commit,
      mpi.size() > 1 ? 1 : 0);
  auto failed_state = fixture.make_state();
  const auto failed = facade.attempt(failed_state, physics, stencil, {}, {});
  flow::test::ImmersedFlowTestAccess::clear_failure_stage();
  HUNDUN_CHECK(std::get<flow::StepAttemptReport>(failed.base).disposition !=
               flow::StepAttemptDisposition::committed);
  HUNDUN_CHECK(!failed.force.has_value());
  auto failed_source = facade.diagnostic_source(failed_state, failed);
  bool unavailable = false;
  try {
    RecordingSink sink;
    diagnostics::collect_diagnostics(
        failed_source, diagnostics::DiagnosticModuleKind::wall_force,
        request(mpi.rank(), diagnostics::DiagnosticScope::local), sink);
  } catch (const diagnostics::DiagnosticCollectionError &error) {
    unavailable = error.classification() ==
                      diagnostics::DiagnosticFailureClass::unavailable &&
                  error.code() == "stage3.wall-force.diagnostics.unavailable";
  }
  HUNDUN_CHECK(unavailable);
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  return hundun::test::run(run);
}
