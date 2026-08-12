// SPDX-License-Identifier: Apache-2.0

#include "tests/support/flow_ideal_gas_closure_test_access.hpp"
#include "tests/support/flow_immersed_test_access.hpp"
#include "tests/support/stage3_state_equality.hpp"
#include "tests/support/stage3_stl_fixture.hpp"
#include "tests/support/test_main.hpp"

#include "hundun/bc_basic_boundary.hpp"
#include "hundun/cfg_resolved_case_v3.hpp"
#include "hundun/exec_execution.hpp"
#include "hundun/flow_checkpoint_v3.hpp"
#include "hundun/flow_ideal_gas_piso.hpp"
#include "hundun/flow_immersed.hpp"
#include "hundun/fvm_cell_centered.hpp"
#include "hundun/ib_domain.hpp"
#include "hundun/ib_ghost_stencil_plan.hpp"
#include "hundun/ib_local_flow_pattern.hpp"
#include "hundun/ib_surface.hpp"
#include "hundun/ib_surface_query.hpp"
#include "hundun/les_wale.hpp"
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

#include "src/flow_checkpoint_v3_detail.hpp"
#include "src/rt_checkpoint_v2_protocol_detail.hpp"

#include <mpi.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace {

using namespace hundun;

int selected_cells = 8;
bool formal_run{};
constexpr double kDt = 1.0e-4;
constexpr double kMu = 0.01;
constexpr double kCp = 1000.0;
constexpr double kGasConstant = 287.05;
constexpr double kPressure = 101325.0;

runtime::Int3 process_grid(int ranks) {
  if (ranks == 1)
    return {1, 1, 1};
  if (ranks == 2)
    return {2, 1, 1};
  if (ranks == 4)
    return {2, 2, 1};
  throw runtime::Error("unsupported Checkpoint v3 density-profile ranks");
}

bool has_ibm(flow::CheckpointV3Presence profile) noexcept {
  const auto value = static_cast<std::uint8_t>(profile);
  return value == 1U || value == 3U || value == 4U || value == 6U ||
         value == 7U || value == 9U;
}

bool has_wale(flow::CheckpointV3Presence profile) noexcept {
  const auto value = static_cast<std::uint8_t>(profile);
  return value == 2U || value == 3U || value == 5U || value == 6U ||
         value == 8U || value == 9U;
}

bool has_ideal(flow::CheckpointV3Presence profile) noexcept {
  return static_cast<std::uint8_t>(profile) >= 7U;
}

config::DensityModel density_model(flow::CheckpointV3Presence profile) {
  const auto value = static_cast<std::uint8_t>(profile);
  if (value <= 3U)
    return config::DensityModel::constant;
  if (value <= 6U)
    return config::DensityModel::material;
  return config::DensityModel::ideal_gas;
}

bool open_domain(flow::CheckpointV3Presence profile) noexcept {
  return profile == flow::CheckpointV3Presence::ideal_gas_static_ibm;
}

config::FlowCaseConfig flow_case(flow::CheckpointV3Presence profile,
                                 int ranks) {
  config::FlowCaseConfig result{};
  result.schema_version = 2;
  result.case_name = "stage3-r2-checkpoint-profile-" +
                     std::to_string(static_cast<unsigned>(profile));
  result.simulation_type = config::SimulationType::variable_density_flow;
  result.density_model = density_model(profile);
  result.resources.expected_ranks = ranks;
  result.resources.process_grid = process_grid(ranks);
  result.mesh.cells = {selected_cells, selected_cells, selected_cells};
  result.mesh.origin_m = {0.0, 0.0, 0.0};
  result.mesh.length_m = {1.0, 1.0, 1.0};
  result.physics.rho_ref_kg_per_m3 = 1.0;
  result.physics.dynamic_viscosity_pa_s = kMu;
  result.physics.inlet_consistency_rtol = 1.0e-12;
  if (density_model(profile) == config::DensityModel::material)
    result.scalars.push_back({"alpha", 0.0});
  if (has_ideal(profile)) {
    result.physics.cp_J_per_kg_K = kCp;
    result.physics.gas_constant_J_per_kg_K = kGasConstant;
    result.physics.thermodynamic_pressure_pa = kPressure;
  }
  result.time.mode = config::TimeMode::fixed;
  result.time.steps = formal_run ? 4U : 2U;
  result.time.initial_dt_s = kDt;
  result.time.min_dt_s = kDt;
  result.time.max_dt_s = kDt;
  result.time.cfl_target = 0.5;
  result.time.diffusion_number_target = 0.25;
  result.time.growth_factor = 1.25;
  result.time.retry_factor = 0.5;
  result.time.max_retries = 8U;
  result.restart.read = false;
  result.restart.write_directory = "checkpoint";
  result.restart.write_interval = 1;
  result.diagnostics.directory = "diagnostics";
  result.diagnostics.write_interval = 1;
  result.diagnostics.write_mesh = false;
  result.performance.enabled = false;
  result.performance.directory = "performance";
  result.performance.warmup_steps = 1;
  result.performance.measured_steps = 1;
  result.performance.repetitions = 1;
  constexpr std::array<config::PatchName, 6> names{
      config::PatchName::x_min, config::PatchName::x_max,
      config::PatchName::y_min, config::PatchName::y_max,
      config::PatchName::z_min, config::PatchName::z_max};
  for (std::size_t patch = 0U; patch < names.size(); ++patch) {
    result.boundaries[patch].patch = names[patch];
    result.boundaries[patch].type = config::BoundaryType::periodic;
  }
  if (open_domain(profile)) {
    result.boundaries[0].type = config::BoundaryType::velocity_inlet;
    result.boundaries[0].velocity_m_per_s = {0.0, 0.0, 0.0};
    result.boundaries[0].thermal_authority =
        config::InletThermalAuthority::temperature;
    result.boundaries[0].temperature_K = 300.0;
    result.boundaries[0].enthalpy_J_per_kg = kCp * 300.0;
    result.boundaries[0].density_kg_per_m3 =
        kPressure / (kGasConstant * 300.0);
    result.boundaries[0].scalar_values =
        std::vector<config::InletScalarValue>{};
    result.boundaries[1].type = config::BoundaryType::pressure_outlet;
    result.boundaries[1].pressure_perturbation_pa = 0.0;
    for (std::size_t patch = 2U; patch < names.size(); ++patch)
      result.boundaries[patch].type = config::BoundaryType::symmetry;
  }
  return result;
}

config::ImmersedFlowCaseConfig
checkpoint_case(flow::CheckpointV3Presence profile, int ranks) {
  config::ImmersedFlowCaseConfig result;
  result.common_flow = flow_case(profile, ranks);
  result.immersed_boundary.model =
      has_ibm(profile)
          ? config::ImmersedBoundaryModel::local_flow_pattern_ghost_cell
          : config::ImmersedBoundaryModel::none;
  if (has_ibm(profile)) {
    result.immersed_boundary.geometry = config::StlGeometryConfig{
        "body.stl", 1.0, config::ImmersedFluidSide::outside};
    result.immersed_boundary.wall = config::StaticImmersedWallConfig{};
  }
  result.les.model =
      has_wale(profile) ? config::LesModel::wale : config::LesModel::none;
  if (has_wale(profile))
    result.les.wale = config::WaleConfig{0.5, 0.9, 0.7};
  return result;
}

runtime::FieldDescriptor cell_field(const char *name, const char *unit,
                                    std::uint32_t components, int ghost_width,
                                    bool conservative) {
  return {name,
          unit,
          "stage3_r2_checkpoint_density_profiles",
          runtime::FunctionSpace::cell_average,
          runtime::ScalarType::float64,
          components,
          ghost_width,
          conservative,
          runtime::RestartPolicy::persistent,
          runtime::OutputPolicy::never};
}

runtime::FieldDescriptor face_field(const char *name,
                                    std::uint32_t components) {
  return {name,
          "1",
          "stage3_r2_checkpoint_density_profiles",
          runtime::FunctionSpace::face_value,
          runtime::ScalarType::float64,
          components,
          0,
          false,
          runtime::RestartPolicy::persistent,
          runtime::OutputPolicy::never};
}

struct ProfileContext final {
  ProfileContext(const runtime::MpiContext &mpi,
                 flow::CheckpointV3Presence selected)
      : profile(selected), config(checkpoint_case(profile, mpi.size())),
        decomposition(runtime::StructuredDecomposition::create(
            mpi, {selected_cells, selected_cells, selected_cells},
            {!open_domain(profile), !open_domain(profile),
             !open_domain(profile)},
            runtime::DecompositionOptions{process_grid(mpi.size())})),
        topology(decomposition),
        geometry(topology,
                 mesh::UniformBoxMapping({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0})),
        boundaries(boundary::BoundaryRegistry::create(config.common_flow,
                                                       topology)),
        directory("stage3-r2-checkpoint-profile-" +
                  std::to_string(static_cast<unsigned>(profile))) {
    if (has_ibm(profile)) {
      std::string path_text;
      if (mpi.rank() == 0) {
        const auto triangles = test::projected_octahedral_sphere(
            {0.5, 0.5, 0.5}, 0.2, 2U);
        const auto path = directory.path() / "body.stl";
        test::write_text(path, test::ascii_stl(triangles, "body"));
        path_text = path.string();
      }
      std::uint64_t size = path_text.size();
      HUNDUN_CHECK(MPI_Bcast(&size, 1, MPI_UINT64_T, 0, mpi.comm()) ==
                   MPI_SUCCESS);
      path_text.resize(static_cast<std::size_t>(size));
      HUNDUN_CHECK(MPI_Bcast(path_text.data(), static_cast<int>(size), MPI_BYTE,
                             0, mpi.comm()) == MPI_SUCCESS);
      surface.emplace(immersed::ImmersedSurface::load_collective(
          std::filesystem::path(path_text), 1.0, mpi, 0));
      query.emplace(immersed::SurfaceQuery::create(*surface));
      domain.emplace(immersed::ImmersedDomain::create(
          *surface, *query, config::ImmersedFluidSide::outside, topology,
          geometry, boundaries, mpi));
      ghost_plan.emplace(immersed::GhostStencilPlan::create(
          *surface, *query, *domain, topology, geometry, decomposition, mpi));
      wall_plan.emplace(immersed::WallQuadraturePlan::create(
          *surface, *query, *domain, topology, geometry, mpi));
      HUNDUN_CHECK(ghost_plan->maximum_halo_reach() <= 4U);
    }

    const int ghost_width = has_ibm(profile) ? 4 : 2;
    fields.density = registry.declare_field(
        cell_field("rho", "kg/m3", 1U, ghost_width, true));
    fields.velocity = registry.declare_field(
        cell_field("velocity", "m/s", 3U, ghost_width, false));
    fields.mechanical_pressure = registry.declare_field(
        cell_field("pi", "Pa", 1U, ghost_width, false));
    fields.face_velocity =
        registry.declare_field(face_field("face_velocity", 3U));
    fields.face_mass_flux = finite_volume::declare_face_mass_flux(registry);
    if (density_model(profile) != config::DensityModel::constant) {
      rho_h = registry.declare_field(
          cell_field("rho_h", "J/m3", 1U, ghost_width, true));
      fields.transported_cell_fields.push_back(*rho_h);
    }
    if (density_model(profile) == config::DensityModel::material) {
      rho_alpha = registry.declare_field(
          cell_field("rho_alpha", "kg/m3", 1U, ghost_width, true));
      fields.transported_cell_fields.push_back(*rho_alpha);
    }
    registry.freeze();
  }

  flow::FlowState make_state() const {
    auto state = flow::FlowState::create(
        registry,
        {decomposition.local_extent(), topology.local_face_count()}, fields,
        {0U, 0.0, kDt, 0.0, flow::MomentumTimeOrder::backward_euler});
    flow::FlowLayerValues values;
    values.density.assign(topology.owned_cell_count(), 0.0);
    values.velocity.assign(topology.owned_cell_count() * 3U, 0.0);
    values.mechanical_pressure.assign(topology.owned_cell_count(), 0.0);
    values.face_velocity.assign(topology.local_face_count() * 3U, 0.0);
    values.face_mass_flux.assign(topology.local_face_count(), 0.0);
    values.transported_cell_fields.resize(fields.transported_cell_fields.size());
    for (auto &transported : values.transported_cell_fields)
      transported.assign(topology.owned_cell_count(), 0.0);
    for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count(); ++cell) {
      if (domain.has_value() &&
          domain->region(cell) != immersed::CellRegion::fluid)
        continue;
      double rho = 1.0;
      if (has_ideal(profile))
        rho = kPressure / (kGasConstant * 300.0);
      values.density[cell] = rho;
      if (rho_h.has_value())
        values.transported_cell_fields[0][cell] = rho * kCp * 300.0;
      if (rho_alpha.has_value())
        values.transported_cell_fields[1][cell] = rho * 0.2;
    }
    state.seed_accepted_layers(values, values);
    return state;
  }

  std::vector<mesh::GlobalCellId> active_global_cells() const {
    if (domain.has_value())
      return domain->active_cells().ordered_global_ids();
    std::vector<mesh::GlobalCellId> result;
    result.reserve(topology.local_cell_count());
    for (mesh::LocalCellId cell = 0U; cell < topology.local_cell_count(); ++cell)
      result.push_back(topology.global_cell_id(cell));
    return result;
  }

  std::size_t owned_active_count() const noexcept {
    return domain.has_value() ? domain->active_cells().owned_active_count()
                              : topology.owned_cell_count();
  }

  std::filesystem::path shared_path(const runtime::MpiContext &mpi,
                                    std::string name) const {
    std::string text;
    if (mpi.rank() == 0)
      text = (directory.path() / std::move(name)).string();
    std::uint64_t size = text.size();
    HUNDUN_CHECK(MPI_Bcast(&size, 1, MPI_UINT64_T, 0, mpi.comm()) ==
                 MPI_SUCCESS);
    text.resize(static_cast<std::size_t>(size));
    HUNDUN_CHECK(MPI_Bcast(text.data(), static_cast<int>(size), MPI_BYTE, 0,
                           mpi.comm()) == MPI_SUCCESS);
    return text;
  }

  flow::CheckpointV3Presence profile;
  config::ImmersedFlowCaseConfig config;
  runtime::StructuredDecomposition decomposition;
  mesh::MeshTopology topology;
  mesh::MeshGeometry geometry;
  boundary::BoundaryRegistry boundaries;
  test::Stage3TemporaryDirectory directory;
  std::optional<immersed::ImmersedSurface> surface;
  std::optional<immersed::SurfaceQuery> query;
  std::optional<immersed::ImmersedDomain> domain;
  std::optional<immersed::GhostStencilPlan> ghost_plan;
  std::optional<immersed::WallQuadraturePlan> wall_plan;
  immersed::LocalFlowPatternTransform transform;
  runtime::FieldRegistry registry;
  flow::FlowFieldIds fields;
  std::optional<runtime::FieldId> rho_h;
  std::optional<runtime::FieldId> rho_alpha;
};

struct FlowBundle final {
  FlowBundle(const runtime::MpiContext &mpi, const ProfileContext &context,
             const flow::FlowState &state)
      : halo(runtime::HaloExchange::create(
            context.decomposition,
            runtime::ExchangePlan::create(
                context.decomposition, context.decomposition.local_extent(),
                has_ibm(context.profile) ? 4 : 2))),
        momentum_solver(execution, mpi), pressure_solver(execution, mpi),
        mx(execution), my(execution), mz(execution), pressure_pc(execution) {
    if (has_ideal(context.profile)) {
      const flow::IdealGasClosureSpec spec{*context.rho_h, kCp, kGasConstant,
                                           kPressure};
      if (has_ibm(context.profile)) {
        closure.emplace(flow::test::IdealGasClosureTestAccess::create_immersed(
            context.topology, context.geometry, context.boundaries,
            *context.domain, mpi, context.registry, context.fields, state,
            spec));
      } else {
        closure.emplace(flow::IdealGasClosure::create(
            context.topology, context.geometry, context.boundaries, mpi,
            context.registry, context.fields, state, spec));
      }
    }
    if (has_wale(context.profile)) {
      wale.emplace(les::WaleModel::create(
          {0.5, 0.9, 0.7}, context.topology, context.geometry,
          context.owned_active_count(), context.active_global_cells(),
          execution));
    }
    const auto *domain =
        context.domain.has_value() ? &*context.domain : nullptr;
    const auto *ghost =
        context.ghost_plan.has_value() ? &*context.ghost_plan : nullptr;
    const auto *wall =
        context.wall_plan.has_value() ? &*context.wall_plan : nullptr;
    const auto *transform = has_ibm(context.profile) ? &context.transform : nullptr;
    const auto *wale_ptr = wale.has_value() ? &*wale : nullptr;
    if (density_model(context.profile) == config::DensityModel::constant) {
      facade.emplace(flow::FixedStepImmersedFlow::create(
          context.decomposition, context.topology, context.geometry,
          context.boundaries, domain, ghost, wall, transform, wale_ptr, mpi,
          execution, halo, momentum_solver, {&mx, &my, &mz}, pressure_solver,
          pressure_pc));
    } else {
      const flow::MaterialDensityTransportSpec material_spec{
          *context.rho_h, 0.0,
          context.rho_alpha.has_value()
              ? std::vector<runtime::FieldId>{*context.rho_alpha}
              : std::vector<runtime::FieldId>{},
          context.rho_alpha.has_value() ? std::vector<double>{0.0}
                                        : std::vector<double>{}};
      const flow::ImmersedFlowDensitySetup setup{
          density_model(context.profile), &context.registry, context.fields,
          material_spec, closure.has_value() ? &*closure : nullptr};
      facade.emplace(flow::FixedStepImmersedFlow::create(
          context.decomposition, context.topology, context.geometry,
          context.boundaries, domain, ghost, wall, transform, wale_ptr, setup,
          mpi, execution, halo, momentum_solver, {&mx, &my, &mz},
          pressure_solver, pressure_pc));
    }
  }

  execution::CpuReferenceContext execution;
  runtime::HaloExchange halo;
  linear::ConjugateGradientSolver momentum_solver;
  linear::BiCGStabSolver pressure_solver;
  linear::JacobiPreconditioner mx;
  linear::JacobiPreconditioner my;
  linear::JacobiPreconditioner mz;
  linear::JacobiPreconditioner pressure_pc;
  std::optional<flow::IdealGasClosure> closure;
  std::optional<les::WaleModel> wale;
  std::optional<flow::FixedStepImmersedFlow> facade;
};

flow::CheckpointV3WriteModules write_modules(const ProfileContext &context,
                                             const FlowBundle &bundle) {
  flow::CheckpointV3WriteModules result;
  result.presence = context.profile;
  result.wale = bundle.wale.has_value() ? &*bundle.wale : nullptr;
  result.ideal_gas = bundle.closure.has_value() ? &*bundle.closure : nullptr;
  if (has_ibm(context.profile)) {
    result.surface = &*context.surface;
    result.query = &*context.query;
    result.domain = &*context.domain;
    result.ghost_plan = &*context.ghost_plan;
    result.wall_plan = &*context.wall_plan;
    result.transform = &context.transform;
    result.flow = &*bundle.facade;
  }
  return result;
}

flow::CheckpointV3ReadModules read_modules(const ProfileContext &context,
                                           FlowBundle &bundle) {
  flow::CheckpointV3ReadModules result;
  result.presence = context.profile;
  result.wale = bundle.wale.has_value() ? &*bundle.wale : nullptr;
  result.ideal_gas = bundle.closure.has_value() ? &*bundle.closure : nullptr;
  if (has_ibm(context.profile)) {
    result.surface = &*context.surface;
    result.query = &*context.query;
    result.domain = &*context.domain;
    result.ghost_plan = &*context.ghost_plan;
    result.wall_plan = &*context.wall_plan;
    result.transform = &context.transform;
    result.flow = &*bundle.facade;
  }
  return result;
}

const flow::StepAttemptReport &
base_report(const flow::DensityStepAttemptReport &report) {
  if (const auto *constant = std::get_if<flow::StepAttemptReport>(&report))
    return *constant;
  if (const auto *material =
          std::get_if<flow::MaterialDensityStepAttemptReport>(&report))
    return material->flow();
  return std::get<flow::IdealGasStepAttemptReport>(report).flow().flow();
}

flow::ImmersedFlowStepAttemptReport
advance(const ProfileContext &context, FlowBundle &bundle,
        flow::FlowState &state) {
  const bool ideal = has_ideal(context.profile);
  const flow::ImmersedFlowPhysics physics{
      density_model(context.profile), 1.0,
      kMu,
      ideal ? std::optional<double>{kCp} : std::nullopt,
      ideal ? std::optional<double>{kGasConstant} : std::nullopt,
      ideal ? std::optional<double>{kPressure} : std::nullopt};
  const auto metadata = state.metadata();
  const auto stencil = flow::make_momentum_time_stencil(
      metadata.step == 0U ? flow::MomentumTimeOrder::backward_euler
                          : flow::MomentumTimeOrder::bdf2,
      kDt, metadata.step == 0U ? 0.0 : kDt);
  linear::SolveControl control;
  control.max_iterations = 5000U;
  control.residual_recompute_interval = 20U;
  auto result = bundle.facade->attempt(state, physics, stencil, control, control);
  const auto &base = base_report(result.base);
  HUNDUN_CHECK(base.disposition == flow::StepAttemptDisposition::committed);
  HUNDUN_CHECK(base.pressure_corrector_count == 2U);
  HUNDUN_CHECK(result.force.has_value() == has_ibm(context.profile));
  HUNDUN_CHECK(result.wale.has_value() == has_wale(context.profile));
  if (result.wale.has_value())
    HUNDUN_CHECK(result.wale->identity.value != 0U);
  return result;
}

bool same_wall_gradients(
    const std::vector<flow::detail::ImmersedFlowWallGradient> &left,
    const std::vector<flow::detail::ImmersedFlowWallGradient> &right) noexcept {
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

bool same_pressure_authority(
    const flow::test::ImmersedFlowPressureAuthorityStateSnapshot &left,
    const flow::test::ImmersedFlowPressureAuthorityStateSnapshot &right) {
  return left.history_available == right.history_available &&
         left.committed_available == right.committed_available &&
         same_wall_gradients(left.history, right.history) &&
         same_wall_gradients(left.committed, right.committed) &&
         same_wall_gradients(left.pending, right.pending);
}

bool same_accepted_state(const test::Stage3StateSnapshot &left,
                         const test::Stage3StateSnapshot &right) noexcept {
  return test::flow_layer_values_bitwise_equal(left.history, right.history) &&
         test::flow_layer_values_bitwise_equal(left.committed,
                                               right.committed) &&
         test::accepted_step_metadata_bitwise_equal(left.metadata,
                                                    right.metadata);
}

void mix(std::uint64_t &hash, std::uint64_t value) noexcept {
  hash ^= value;
  hash *= UINT64_C(1099511628211);
}

void mix_values(std::uint64_t &hash, const std::vector<double> &values) {
  mix(hash, values.size());
  for (const double value : values)
    mix(hash, test::fp64_bits(value));
}

std::uint64_t state_bits(const test::Stage3StateSnapshot &snapshot) {
  std::uint64_t hash = UINT64_C(1469598103934665603);
  const auto layer = [&](const flow::FlowLayerValues &values) {
    mix_values(hash, values.density);
    mix_values(hash, values.velocity);
    mix_values(hash, values.mechanical_pressure);
    mix_values(hash, values.face_velocity);
    mix_values(hash, values.face_mass_flux);
    mix(hash, values.transported_cell_fields.size());
    for (const auto &field : values.transported_cell_fields)
      mix_values(hash, field);
  };
  layer(snapshot.history);
  layer(snapshot.committed);
  layer(snapshot.trial);
  return hash;
}

std::uint64_t metadata_bits(const flow::AcceptedStepMetadata &metadata) {
  std::uint64_t hash = UINT64_C(1469598103934665603);
  mix(hash, metadata.step);
  mix(hash, test::fp64_bits(metadata.time_s));
  mix(hash, test::fp64_bits(metadata.dt_s));
  mix(hash, test::fp64_bits(metadata.previous_dt_s));
  mix(hash, static_cast<std::uint8_t>(metadata.order));
  return hash;
}

std::uint64_t pressure_bits(
    const flow::test::ImmersedFlowPressureAuthorityStateSnapshot &authority) {
  std::uint64_t hash = UINT64_C(1469598103934665603);
  mix(hash, authority.history_available ? 1U : 0U);
  mix(hash, authority.committed_available ? 1U : 0U);
  for (const auto *rows : {&authority.history, &authority.committed,
                           &authority.pending}) {
    mix(hash, rows->size());
    for (const auto &row : *rows) {
      mix(hash, row.link);
      mix(hash, test::fp64_bits(row.normal_gradient_pa_per_m));
    }
  }
  return hash;
}

std::uint64_t closure_bits(const flow::IdealGasClosureState &state) {
  std::uint64_t hash = UINT64_C(1469598103934665603);
  mix(hash, static_cast<std::uint8_t>(state.mode));
  mix(hash, test::fp64_bits(state.thermodynamic_pressure_pa));
  mix(hash, state.target_mass_kg.has_value() ? 1U : 0U);
  if (state.target_mass_kg.has_value())
    mix(hash, test::fp64_bits(*state.target_mass_kg));
  mix(hash, state.revision);
  return hash;
}

std::vector<std::uint8_t> read_bytes(const std::filesystem::path &path) {
  std::ifstream input(path, std::ios::binary);
  if (!input)
    throw runtime::Error("unable to read Checkpoint v3 density manifest");
  return {std::istreambuf_iterator<char>(input),
          std::istreambuf_iterator<char>()};
}

void corrupt_manifest(const runtime::MpiContext &mpi,
                      const std::filesystem::path &directory) {
  if (mpi.rank() == 0) {
    auto bytes = read_bytes(directory / "manifest.v3.bin");
    HUNDUN_CHECK(!bytes.empty());
    bytes.back() ^= 1U;
    test::write_bytes(directory / "manifest.v3.bin", bytes);
  }
  mpi.barrier();
}

void rewrite_ideal_authority(const runtime::MpiContext &mpi,
                             const std::filesystem::path &directory) {
  if (mpi.rank() == 0) {
    const auto manifest_path = directory / "manifest.v3.bin";
    auto manifest = flow::detail::decode_checkpoint_v3_manifest(
        read_bytes(manifest_path));
    HUNDUN_CHECK(manifest.ideal_gas.has_value());
    manifest.ideal_gas->thermodynamic_pressure_pa *= 1.001;
    const auto manifest_bytes =
        flow::detail::encode_checkpoint_v3_manifest(manifest);
    test::write_bytes(manifest_path, manifest_bytes);
    auto marker = runtime::checkpoint_v2::decode_completed_marker(
        read_bytes(directory / "COMPLETED"));
    marker.manifest_actual_size = manifest_bytes.size();
    marker.manifest_crc64 =
        flow::detail::checkpoint_v3_manifest_crc64(manifest_bytes);
    test::write_bytes(directory / "COMPLETED",
                      runtime::checkpoint_v2::encode_completed_marker(marker));
  }
  mpi.barrier();
}

void check_missing_module(const runtime::MpiContext &mpi,
                          const ProfileContext &context,
                          const FlowBundle &bundle,
                          const flow::FlowState &state) {
  auto modules = write_modules(context, bundle);
  if (has_ibm(context.profile))
    modules.flow = nullptr;
  else
    modules.wale = nullptr;
  const auto path = context.shared_path(mpi, "missing-module");
  const auto report = flow::write_checkpoint_v3(
      mpi, context.decomposition, context.topology, context.geometry,
      context.boundaries, context.config, modules, state, {kDt, 0U}, path);
  HUNDUN_CHECK(report.disposition() == flow::CheckpointV3Disposition::failed);
  HUNDUN_CHECK(report.reason() == flow::CheckpointV3FailureReason::presence);
  HUNDUN_CHECK(!std::filesystem::exists(path));
}

void check_inactive_negative_zero(const runtime::MpiContext &mpi,
                                  const ProfileContext &context,
                                  const FlowBundle &bundle) {
  HUNDUN_CHECK(context.profile ==
               flow::CheckpointV3Presence::material_static_ibm);
  auto invalid = context.make_state();
  auto values = invalid.snapshot(flow::FlowLayer::committed);
  bool mutated = false;
  for (mesh::LocalCellId cell = 0U; cell < context.topology.owned_cell_count();
       ++cell) {
    if (context.domain->region(cell) != immersed::CellRegion::fluid) {
      values.density[cell] = -0.0;
      mutated = true;
      break;
    }
  }
  int local_mutated = mutated ? 1 : 0;
  int global_mutated = 0;
  HUNDUN_CHECK(MPI_Allreduce(&local_mutated, &global_mutated, 1, MPI_INT,
                             MPI_MAX, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(global_mutated == 1);
  invalid.seed_accepted_layers(values, values);
  const auto path = context.shared_path(mpi, "inactive-negative-zero");
  const auto report = flow::write_checkpoint_v3(
      mpi, context.decomposition, context.topology, context.geometry,
      context.boundaries, context.config, write_modules(context, bundle),
      invalid, {kDt, 0U}, path);
  HUNDUN_CHECK(report.disposition() == flow::CheckpointV3Disposition::failed);
  HUNDUN_CHECK(report.reason() == flow::CheckpointV3FailureReason::state);
  HUNDUN_CHECK(!std::filesystem::exists(path));
}

void run_profile(const runtime::MpiContext &mpi,
                 flow::CheckpointV3Presence profile) {
  ProfileContext context(mpi, profile);
  auto continuous_state = context.make_state();
  FlowBundle continuous(mpi, context, continuous_state);
  check_missing_module(mpi, context, continuous, continuous_state);
  if (profile == flow::CheckpointV3Presence::material_static_ibm)
    check_inactive_negative_zero(mpi, context, continuous);
  const std::uint64_t total_steps = formal_run ? 4U : 2U;
  const std::uint64_t restart_step = formal_run ? 2U : 1U;
  std::optional<flow::ImmersedFlowStepAttemptReport> continuous_report;
  for (std::uint64_t step = 0U; step < restart_step; ++step)
    continuous_report = advance(context, continuous, continuous_state);

  const auto checkpoint_directory = context.shared_path(mpi, "continuation");
  const auto written = flow::write_checkpoint_v3(
      mpi, context.decomposition, context.topology, context.geometry,
      context.boundaries, context.config, write_modules(context, continuous),
      continuous_state, {kDt, 0U}, checkpoint_directory);
  if (written.disposition() != flow::CheckpointV3Disposition::completed &&
      mpi.rank() == 0) {
    std::cerr << "CHECKPOINT_PROFILE_WRITE_FAILED profile="
              << static_cast<unsigned>(profile)
              << " reason=" << static_cast<unsigned>(written.reason())
              << " phase=" << static_cast<unsigned>(written.phase()) << '\n';
  }
  HUNDUN_CHECK(written.disposition() ==
               flow::CheckpointV3Disposition::completed);
  HUNDUN_CHECK(written.presence() == profile);

  auto restarted_state = context.make_state();
  FlowBundle restarted(mpi, context, restarted_state);
  const auto before_wrong_profile = test::snapshot_stage3_state(restarted_state);
  auto wrong = read_modules(context, restarted);
  wrong.presence = static_cast<flow::CheckpointV3Presence>(255U);
  const auto wrong_result = flow::read_checkpoint_v3(
      mpi, context.decomposition, context.topology, context.geometry,
      context.boundaries, context.config, wrong, restarted_state,
      checkpoint_directory);
  HUNDUN_CHECK(!wrong_result.restored());
  HUNDUN_CHECK(wrong_result.report().reason() ==
               flow::CheckpointV3FailureReason::presence);
  HUNDUN_CHECK(test::stage3_state_bitwise_equal(
      before_wrong_profile, test::snapshot_stage3_state(restarted_state)));

  const auto restored = flow::read_checkpoint_v3(
      mpi, context.decomposition, context.topology, context.geometry,
      context.boundaries, context.config, read_modules(context, restarted),
      restarted_state, checkpoint_directory);
  HUNDUN_CHECK(restored.restored());
  HUNDUN_CHECK(restored.report().presence() == profile);
  HUNDUN_CHECK(same_accepted_state(
      test::snapshot_stage3_state(continuous_state),
      test::snapshot_stage3_state(restarted_state)));
  if (has_ibm(profile)) {
    HUNDUN_CHECK(same_pressure_authority(
        flow::test::ImmersedFlowTestAccess::pressure_authority_state(
            *continuous.facade),
        flow::test::ImmersedFlowTestAccess::pressure_authority_state(
            *restarted.facade)));
  }
  if (has_ideal(profile)) {
    HUNDUN_CHECK(test::ideal_gas_closure_state_bitwise_equal(
        continuous.closure->state(), restarted.closure->state()));
  }

  for (std::uint64_t step = restart_step; step < total_steps; ++step) {
    continuous_report = advance(context, continuous, continuous_state);
    const auto restarted_report = advance(context, restarted, restarted_state);
    if (has_wale(profile)) {
      HUNDUN_CHECK(continuous_report->wale->identity ==
                   restarted_report.wale->identity);
      HUNDUN_CHECK(test::fp64_bits(
                       continuous_report->wale->l2_nu_t_m2_per_s) ==
                   test::fp64_bits(
                       restarted_report.wale->l2_nu_t_m2_per_s));
    }
  }
  const auto final_continuous = test::snapshot_stage3_state(continuous_state);
  const auto final_restarted = test::snapshot_stage3_state(restarted_state);
  HUNDUN_CHECK(test::stage3_state_bitwise_equal(final_continuous,
                                                final_restarted));

  std::optional<flow::test::ImmersedFlowPressureAuthorityStateSnapshot>
      final_authority;
  if (has_ibm(profile)) {
    final_authority =
        flow::test::ImmersedFlowTestAccess::pressure_authority_state(
            *continuous.facade);
    HUNDUN_CHECK(same_pressure_authority(
        *final_authority,
        flow::test::ImmersedFlowTestAccess::pressure_authority_state(
            *restarted.facade)));
  }
  std::optional<flow::IdealGasClosureState> final_closure;
  if (has_ideal(profile)) {
    final_closure = continuous.closure->state();
    HUNDUN_CHECK(test::ideal_gas_closure_state_bitwise_equal(
        *final_closure, restarted.closure->state()));
    HUNDUN_CHECK(final_closure->mode ==
                 (open_domain(profile) ? flow::IdealGasPressureMode::open_fixed
                                       : flow::IdealGasPressureMode::closed_dynamic));
    HUNDUN_CHECK(final_closure->target_mass_kg.has_value() ==
                 !open_domain(profile));
  }

  if (profile == flow::CheckpointV3Presence::ideal_gas_static_ibm_wale) {
    const auto prepare_failure_directory =
        context.shared_path(mpi, "restore-prepare-failure");
    const auto prepare_failure_written = flow::write_checkpoint_v3(
        mpi, context.decomposition, context.topology, context.geometry,
        context.boundaries, context.config, write_modules(context, continuous),
        continuous_state, {kDt, 0U}, prepare_failure_directory);
    HUNDUN_CHECK(prepare_failure_written.disposition() ==
                 flow::CheckpointV3Disposition::completed);
    const auto before_snapshot_failure =
        test::snapshot_stage3_state(restarted_state);
    const auto authority_before_snapshot_failure =
        flow::test::ImmersedFlowTestAccess::pressure_authority_state(
            *restarted.facade);
    const auto closure_before_snapshot_failure = restarted.closure->state();
    flow::test::IdealGasClosureTestAccess::set_restore_preparation_fault(0);
    const auto snapshot_failure = flow::read_checkpoint_v3(
        mpi, context.decomposition, context.topology, context.geometry,
        context.boundaries, context.config, read_modules(context, restarted),
        restarted_state, prepare_failure_directory);
    HUNDUN_CHECK(!snapshot_failure.restored());
    HUNDUN_CHECK(snapshot_failure.report().reason() ==
                 flow::CheckpointV3FailureReason::state);
    HUNDUN_CHECK(snapshot_failure.report().phase() ==
                 flow::CheckpointV3Phase::restore_prepare);
    HUNDUN_CHECK(test::stage3_state_bitwise_equal(
        before_snapshot_failure, test::snapshot_stage3_state(restarted_state)));
    HUNDUN_CHECK(same_pressure_authority(
        authority_before_snapshot_failure,
        flow::test::ImmersedFlowTestAccess::pressure_authority_state(
            *restarted.facade)));
    HUNDUN_CHECK(test::ideal_gas_closure_state_bitwise_equal(
        closure_before_snapshot_failure, restarted.closure->state()));

    rewrite_ideal_authority(mpi, prepare_failure_directory);
    const auto before_prepare_failure =
        test::snapshot_stage3_state(restarted_state);
    const auto authority_before_prepare_failure =
        flow::test::ImmersedFlowTestAccess::pressure_authority_state(
            *restarted.facade);
    const auto closure_before_prepare_failure = restarted.closure->state();
    const auto prepare_failure = flow::read_checkpoint_v3(
        mpi, context.decomposition, context.topology, context.geometry,
        context.boundaries, context.config, read_modules(context, restarted),
        restarted_state, prepare_failure_directory);
    HUNDUN_CHECK(!prepare_failure.restored());
    HUNDUN_CHECK(prepare_failure.report().reason() ==
                 flow::CheckpointV3FailureReason::state);
    HUNDUN_CHECK(prepare_failure.report().phase() ==
                 flow::CheckpointV3Phase::restore_prepare);
    HUNDUN_CHECK(test::stage3_state_bitwise_equal(
        before_prepare_failure, test::snapshot_stage3_state(restarted_state)));
    HUNDUN_CHECK(same_pressure_authority(
        authority_before_prepare_failure,
        flow::test::ImmersedFlowTestAccess::pressure_authority_state(
            *restarted.facade)));
    HUNDUN_CHECK(test::ideal_gas_closure_state_bitwise_equal(
        closure_before_prepare_failure, restarted.closure->state()));
  }

  const auto before_corruption = test::snapshot_stage3_state(restarted_state);
  const auto closure_before_corruption =
      restarted.closure.has_value()
          ? std::optional<flow::IdealGasClosureState>(restarted.closure->state())
          : std::nullopt;
  const auto authority_before_corruption =
      has_ibm(profile)
          ? std::optional<flow::test::ImmersedFlowPressureAuthorityStateSnapshot>(
                flow::test::ImmersedFlowTestAccess::pressure_authority_state(
                    *restarted.facade))
          : std::nullopt;
  corrupt_manifest(mpi, checkpoint_directory);
  const auto corrupted = flow::read_checkpoint_v3(
      mpi, context.decomposition, context.topology, context.geometry,
      context.boundaries, context.config, read_modules(context, restarted),
      restarted_state, checkpoint_directory);
  HUNDUN_CHECK(!corrupted.restored());
  HUNDUN_CHECK(corrupted.report().reason() ==
               flow::CheckpointV3FailureReason::file_integrity);
  HUNDUN_CHECK(corrupted.report().rollback_status() ==
               flow::CheckpointV3CheckStatus::passed);
  HUNDUN_CHECK(test::stage3_state_bitwise_equal(
      before_corruption, test::snapshot_stage3_state(restarted_state)));
  if (closure_before_corruption.has_value())
    HUNDUN_CHECK(test::ideal_gas_closure_state_bitwise_equal(
        *closure_before_corruption, restarted.closure->state()));
  if (authority_before_corruption.has_value())
    HUNDUN_CHECK(same_pressure_authority(
        *authority_before_corruption,
        flow::test::ImmersedFlowTestAccess::pressure_authority_state(
            *restarted.facade)));

  if (formal_run && mpi.rank() == 0) {
    const auto grid = process_grid(mpi.size());
    std::cout << "profile=" << static_cast<unsigned>(profile)
              << " cells=" << selected_cells << " ranks=" << mpi.size()
              << " process_grid=" << grid.x << ',' << grid.y << ',' << grid.z
              << " continuous_steps=" << total_steps
              << " restart_step=" << restart_step << " state_bits="
              << std::hex << state_bits(final_continuous)
              << " metadata_bits=" << metadata_bits(final_continuous.metadata)
              << " pressure_authority_bits="
              << (final_authority.has_value() ? "1:" : "0:")
              << (final_authority.has_value() ? pressure_bits(*final_authority)
                                              : 0U)
              << " optional_closure_bits="
              << (final_closure.has_value() ? "1:" : "0:")
              << (final_closure.has_value() ? closure_bits(*final_closure) : 0U)
              << std::dec << " status=pass\n";
  }
}

void run() {
  auto mpi = runtime::MpiContext::duplicate(MPI_COMM_WORLD);
  HUNDUN_CHECK(mpi.size() == 1 || mpi.size() == 2 || mpi.size() == 4);
  HUNDUN_CHECK(test::stage3_state_equality_oracle_is_mutation_sensitive());
  for (std::uint8_t value = 1U; value <= 9U; ++value)
    run_profile(mpi, static_cast<flow::CheckpointV3Presence>(value));
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  if (argc == 2 && std::string(argv[1]) == "fast") {
    selected_cells = 8;
    formal_run = false;
  } else if (argc == 3 && std::string(argv[1]) == "formal" &&
             std::string(argv[2]) == "12") {
    selected_cells = 12;
    formal_run = true;
  } else {
    return 2;
  }
  return hundun::test::run(run);
}
