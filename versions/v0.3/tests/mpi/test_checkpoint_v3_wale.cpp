// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "tests/support/flow_immersed_test_access.hpp"
#include "tests/support/stage3_mms.hpp"
#include "tests/support/stage3_state_equality.hpp"
#include "tests/support/stage3_stl_fixture.hpp"
#include "tests/support/test_main.hpp"

#include "hundun/bc_basic_boundary.hpp"
#include "hundun/cfg_resolved_case_v3.hpp"
#include "hundun/exec_execution.hpp"
#include "hundun/flow_checkpoint_v3.hpp"
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

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace {

using namespace hundun;

constexpr int kCells = 8;
constexpr int kBodyFittedHaloWidth = 2;
constexpr int kIbmHaloWidth = 4;
constexpr double kRho = 1.0;
constexpr double kMu = 0.01;
constexpr double kDt = 1.0e-4;
constexpr double kPi = 3.141592653589793238462643383279502884;

runtime::Int3 process_grid(int ranks) {
  if (ranks == 1) {
    return {1, 1, 1};
  }
  if (ranks == 2) {
    return {2, 1, 1};
  }
  throw runtime::Error("WALE Checkpoint v3 requires one or two ranks");
}

config::FlowCaseConfig periodic_case(int ranks) {
  config::FlowCaseConfig result{};
  result.schema_version = 2;
  result.case_name = "stage3-r1-checkpoint-wale";
  result.simulation_type = config::SimulationType::variable_density_flow;
  result.density_model = config::DensityModel::constant;
  result.resources.expected_ranks = ranks;
  result.resources.process_grid = process_grid(ranks);
  result.mesh.cells = {kCells, kCells, kCells};
  result.mesh.origin_m = {0.0, 0.0, 0.0};
  result.mesh.length_m = {1.0, 1.0, 1.0};
  result.physics.rho_ref_kg_per_m3 = kRho;
  result.physics.dynamic_viscosity_pa_s = kMu;
  result.physics.inlet_consistency_rtol = 1.0e-12;
  result.time.mode = config::TimeMode::fixed;
  result.time.steps = 2U;
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
  for (std::size_t index = 0U; index < names.size(); ++index) {
    result.boundaries[index].patch = names[index];
    result.boundaries[index].type = config::BoundaryType::periodic;
  }
  return result;
}

config::ImmersedFlowCaseConfig checkpoint_case(int ranks, bool immersed) {
  config::ImmersedFlowCaseConfig result;
  result.common_flow = periodic_case(ranks);
  result.immersed_boundary.model =
      immersed ? config::ImmersedBoundaryModel::local_flow_pattern_ghost_cell
               : config::ImmersedBoundaryModel::none;
  if (immersed) {
    result.immersed_boundary.geometry = config::StlGeometryConfig{
        "body.stl", 1.0, config::ImmersedFluidSide::outside};
    result.immersed_boundary.wall = config::StaticImmersedWallConfig{};
  }
  result.les.model = config::LesModel::wale;
  result.les.wale = config::WaleConfig{0.5, 0.9, 0.7};
  return result;
}

runtime::FieldDescriptor cell_field(const char *name,
                                    std::uint32_t components,
                                    int ghost_width) {
  return {name,
          "1",
          "stage3_r1_checkpoint_wale",
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
          "stage3_r1_checkpoint_wale",
          runtime::FunctionSpace::face_value,
          runtime::ScalarType::float64,
          components,
          0,
          false,
          runtime::RestartPolicy::persistent,
          runtime::OutputPolicy::never};
}

runtime::Real3 velocity(runtime::Real3 point) {
  constexpr double amplitude = 0.05;
  const double x = 2.0 * kPi * point.x;
  const double y = 2.0 * kPi * point.y;
  const double z = 2.0 * kPi * point.z;
  return {amplitude * std::sin(x) * std::cos(y) * std::cos(z),
          -amplitude * std::cos(x) * std::sin(y) * std::cos(z), 0.0};
}

double dot(runtime::Real3 left, runtime::Real3 right) noexcept {
  return left.x * right.x + left.y * right.y + left.z * right.z;
}

struct CaseContext final {
  CaseContext(const runtime::MpiContext &mpi, bool use_immersed)
      : has_immersed(use_immersed),
        halo_width(use_immersed ? kIbmHaloWidth : kBodyFittedHaloWidth),
        directory(use_immersed ? "s3-r1-ibm-wale-checkpoint"
                               : "s3-r1-wale-checkpoint"),
        decomposition(runtime::StructuredDecomposition::create(
            mpi, {kCells, kCells, kCells}, {true, true, true},
            runtime::DecompositionOptions{process_grid(mpi.size())})),
        topology(decomposition),
        geometry(topology,
                 mesh::UniformBoxMapping({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0})),
        boundaries(boundary::BoundaryRegistry::create(
            periodic_case(mpi.size()), topology)) {
    if (has_immersed) {
      std::string path_text;
      if (mpi.rank() == 0) {
        const auto path = directory.path() / "body.stl";
        const auto body =
            test::stage3::approved_body(test::stage3::BodyKind::sphere);
        const auto surface_fixture = test::stage3::make_manufactured_surface(
            body, 1.0 / static_cast<double>(kCells));
        test::write_text(
            path, test::ascii_stl(surface_fixture.triangles, "body"));
        path_text = path.string();
      }
      std::uint64_t path_size = path_text.size();
      HUNDUN_CHECK(MPI_Bcast(&path_size, 1, MPI_UINT64_T, 0, mpi.comm()) ==
                   MPI_SUCCESS);
      HUNDUN_CHECK(path_size <=
                   static_cast<std::uint64_t>(std::numeric_limits<int>::max()));
      path_text.resize(static_cast<std::size_t>(path_size));
      HUNDUN_CHECK(MPI_Bcast(path_text.data(), static_cast<int>(path_size),
                             MPI_BYTE, 0, mpi.comm()) == MPI_SUCCESS);
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
      HUNDUN_CHECK(ghost_plan->maximum_halo_reach() <=
                   static_cast<std::uint32_t>(halo_width));
      HUNDUN_CHECK(wall_plan->maximum_halo_reach() <=
                   static_cast<std::uint32_t>(halo_width));
    }

    fields.density =
        registry.declare_field(cell_field("rho", 1U, halo_width));
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
        {0U, 0.0, kDt, 0.0, flow::MomentumTimeOrder::backward_euler});
    flow::FlowLayerValues values;
    values.density.assign(topology.owned_cell_count(), 0.0);
    values.velocity.assign(topology.owned_cell_count() * 3U, 0.0);
    values.mechanical_pressure.assign(topology.owned_cell_count(), 0.0);
    values.face_velocity.assign(topology.local_face_count() * 3U, 0.0);
    values.face_mass_flux.assign(topology.local_face_count(), 0.0);
    for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count();
         ++cell) {
      if (domain.has_value() &&
          domain->region(cell) != immersed::CellRegion::fluid) {
        continue;
      }
      const auto value = velocity(geometry.cell_center_m(cell));
      values.density[cell] = kRho;
      values.velocity[cell * 3U] = value.x;
      values.velocity[cell * 3U + 1U] = value.y;
      values.velocity[cell * 3U + 2U] = value.z;
    }
    for (mesh::LocalFaceId face = 0U; face < topology.local_face_count();
         ++face) {
      const bool owner_active =
          !domain.has_value() ||
          domain->region(topology.owner(face)) == immersed::CellRegion::fluid;
      const auto neighbour = topology.neighbour(face);
      const bool neighbour_active =
          !domain.has_value() || !neighbour.has_value() ||
          domain->region(*neighbour) == immersed::CellRegion::fluid;
      if (!owner_active || !neighbour_active) {
        continue;
      }
      const auto value = velocity(geometry.face_center_m(face));
      values.face_velocity[face * 3U] = value.x;
      values.face_velocity[face * 3U + 1U] = value.y;
      values.face_velocity[face * 3U + 2U] = value.z;
      values.face_mass_flux[face] = dot(
          value, geometry.face_area_vector_m2(face, mesh::FaceSide::owner));
    }
    state.seed_accepted_layers(values, values);
    return state;
  }

  std::vector<mesh::GlobalCellId> active_global_cells() const {
    if (domain.has_value()) {
      return domain->active_cells().ordered_global_ids();
    }
    std::vector<mesh::GlobalCellId> result;
    result.reserve(topology.local_cell_count());
    for (mesh::LocalCellId cell = 0U; cell < topology.local_cell_count();
         ++cell) {
      result.push_back(topology.global_cell_id(cell));
    }
    return result;
  }

  std::size_t owned_active_count() const noexcept {
    return domain.has_value() ? domain->active_cells().owned_active_count()
                              : topology.owned_cell_count();
  }

  std::filesystem::path shared_path(const runtime::MpiContext &mpi,
                                    const std::string &name) const {
    std::string text;
    if (mpi.rank() == 0) {
      text = (directory.path() / name).string();
    }
    std::uint64_t size = text.size();
    HUNDUN_CHECK(MPI_Bcast(&size, 1, MPI_UINT64_T, 0, mpi.comm()) ==
                 MPI_SUCCESS);
    HUNDUN_CHECK(size <=
                 static_cast<std::uint64_t>(std::numeric_limits<int>::max()));
    text.resize(static_cast<std::size_t>(size));
    HUNDUN_CHECK(MPI_Bcast(text.data(), static_cast<int>(size), MPI_BYTE, 0,
                           mpi.comm()) == MPI_SUCCESS);
    return std::filesystem::path(text);
  }

  bool has_immersed{};
  int halo_width{};
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
  immersed::LocalFlowPatternTransform transform;
  runtime::FieldRegistry registry;
  flow::FlowFieldIds fields;
};

struct FlowBundle final {
  FlowBundle(const CaseContext &context, const runtime::MpiContext &mpi)
      : halo(runtime::HaloExchange::create(
            context.decomposition,
            runtime::ExchangePlan::create(context.decomposition,
                                          context.decomposition.local_extent(),
                                          context.halo_width))),
        momentum_solver(execution, mpi), pressure_solver(execution, mpi),
        mx(execution), my(execution), mz(execution), pressure_pc(execution) {
    wale.emplace(les::WaleModel::create(
        {0.5, 0.9, 0.7}, context.topology, context.geometry,
        context.owned_active_count(), context.active_global_cells(),
        execution));
    facade.emplace(flow::FixedStepImmersedFlow::create(
        context.decomposition, context.topology, context.geometry,
        context.boundaries,
        context.domain.has_value() ? &*context.domain : nullptr,
        context.ghost_plan.has_value() ? &*context.ghost_plan : nullptr,
        context.wall_plan.has_value() ? &*context.wall_plan : nullptr,
        context.domain.has_value() ? &context.transform : nullptr, &*wale, mpi,
        execution, halo, momentum_solver, {&mx, &my, &mz}, pressure_solver,
        pressure_pc));
  }

  execution::CpuReferenceContext execution;
  runtime::HaloExchange halo;
  linear::ConjugateGradientSolver momentum_solver;
  linear::BiCGStabSolver pressure_solver;
  linear::JacobiPreconditioner mx;
  linear::JacobiPreconditioner my;
  linear::JacobiPreconditioner mz;
  linear::JacobiPreconditioner pressure_pc;
  std::optional<les::WaleModel> wale;
  std::optional<flow::FixedStepImmersedFlow> facade;
};

flow::CheckpointV3WriteModules write_modules(const CaseContext &context,
                                             const FlowBundle &bundle) {
  flow::CheckpointV3WriteModules result;
  result.presence =
      context.has_immersed
          ? flow::CheckpointV3Presence::constant_static_ibm_wale
          : flow::CheckpointV3Presence::constant_body_fitted_wale;
  result.wale = &*bundle.wale;
  if (context.has_immersed) {
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

flow::CheckpointV3ReadModules read_modules(const CaseContext &context,
                                           FlowBundle &bundle) {
  flow::CheckpointV3ReadModules result;
  result.presence =
      context.has_immersed
          ? flow::CheckpointV3Presence::constant_static_ibm_wale
          : flow::CheckpointV3Presence::constant_body_fitted_wale;
  result.wale = &*bundle.wale;
  if (context.has_immersed) {
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

const flow::StepAttemptReport &require_committed(
    const flow::ImmersedFlowStepAttemptReport &report, bool immersed) {
  const auto &base = std::get<flow::StepAttemptReport>(report.base);
  HUNDUN_CHECK(base.disposition == flow::StepAttemptDisposition::committed);
  HUNDUN_CHECK(base.pressure_corrector_count == 2U);
  HUNDUN_CHECK(report.wale.has_value());
  HUNDUN_CHECK(report.wale->identity.value != 0U);
  HUNDUN_CHECK(report.force.has_value() == immersed);
  return base;
}

bool same_wale_summary(const les::WaleSummary &left,
                       const les::WaleSummary &right) noexcept {
  return left.identity == right.identity &&
         test::fp64_bits(left.minimum_nu_t_m2_per_s) ==
             test::fp64_bits(right.minimum_nu_t_m2_per_s) &&
         test::fp64_bits(left.maximum_nu_t_m2_per_s) ==
             test::fp64_bits(right.maximum_nu_t_m2_per_s) &&
         test::fp64_bits(left.l2_nu_t_m2_per_s) ==
             test::fp64_bits(right.l2_nu_t_m2_per_s) &&
         left.exact_zero_count == right.exact_zero_count &&
         left.owned_active_count == right.owned_active_count;
}

bool same_accepted_state(const test::Stage3StateSnapshot &left,
                         const test::Stage3StateSnapshot &right) noexcept {
  return test::flow_layer_values_bitwise_equal(left.history, right.history) &&
         test::flow_layer_values_bitwise_equal(left.committed,
                                               right.committed) &&
         test::accepted_step_metadata_bitwise_equal(left.metadata,
                                                    right.metadata);
}

std::vector<std::uint8_t> read_bytes(const std::filesystem::path &path) {
  std::ifstream input(path, std::ios::binary);
  if (!input) {
    throw runtime::Error("unable to read WALE checkpoint fixture");
  }
  return {std::istreambuf_iterator<char>(input),
          std::istreambuf_iterator<char>()};
}

bool contains_u64_le(const std::vector<std::uint8_t> &bytes,
                     std::uint64_t value) {
  std::array<std::uint8_t, 8> needle{};
  for (std::size_t byte = 0U; byte < needle.size(); ++byte) {
    needle[byte] = static_cast<std::uint8_t>(value >> (8U * byte));
  }
  return std::search(bytes.begin(), bytes.end(), needle.begin(), needle.end()) !=
         bytes.end();
}

void check_module_rejection_before_io(
    const runtime::MpiContext &mpi, const CaseContext &context,
    FlowBundle &bundle, const config::ImmersedFlowCaseConfig &config,
    flow::FlowState &state) {
  const auto missing_path = context.shared_path(mpi, "missing-module");
  auto missing = write_modules(context, bundle);
  if (context.has_immersed) {
    missing.flow = nullptr;
  } else {
    missing.wale = nullptr;
  }
  const auto missing_report = flow::write_checkpoint_v3(
      mpi, context.decomposition, context.topology, context.geometry,
      context.boundaries, config, missing, state, {kDt, 0U}, missing_path);
  HUNDUN_CHECK(missing_report.disposition() ==
               flow::CheckpointV3Disposition::failed);
  HUNDUN_CHECK(missing_report.reason() ==
               flow::CheckpointV3FailureReason::presence);
  HUNDUN_CHECK(!std::filesystem::exists(missing_path));

  const auto extra_path = context.shared_path(mpi, "extra-module");
  auto extra = write_modules(context, bundle);
  if (!context.has_immersed) {
    extra.flow = &*bundle.facade;
  } else {
    extra.presence =
        flow::CheckpointV3Presence::constant_body_fitted_wale;
  }
  const auto extra_report = flow::write_checkpoint_v3(
      mpi, context.decomposition, context.topology, context.geometry,
      context.boundaries, config, extra, state, {kDt, 0U}, extra_path);
  HUNDUN_CHECK(extra_report.disposition() ==
               flow::CheckpointV3Disposition::failed);
  HUNDUN_CHECK(extra_report.reason() ==
               flow::CheckpointV3FailureReason::presence);
  HUNDUN_CHECK(!std::filesystem::exists(extra_path));

  const auto before = test::snapshot_stage3_state(state);
  auto unknown = read_modules(context, bundle);
  unknown.presence = static_cast<flow::CheckpointV3Presence>(255U);
  const auto unknown_result = flow::read_checkpoint_v3(
      mpi, context.decomposition, context.topology, context.geometry,
      context.boundaries, config, unknown, state,
      context.shared_path(mpi, "unknown-profile"));
  HUNDUN_CHECK(!unknown_result.restored());
  HUNDUN_CHECK(unknown_result.report().reason() ==
               flow::CheckpointV3FailureReason::presence);
  HUNDUN_CHECK(test::stage3_state_bitwise_equal(
      before, test::snapshot_stage3_state(state)));
}

void check_manifest_omits_attempt_identity(
    const runtime::MpiContext &mpi, const std::filesystem::path &directory,
    std::uint64_t identity) {
  int absent = 1;
  if (mpi.rank() == 0) {
    const auto bytes = read_bytes(directory / "manifest.v3.bin");
    absent = contains_u64_le(bytes, identity) ? 0 : 1;
  }
  HUNDUN_CHECK(MPI_Bcast(&absent, 1, MPI_INT, 0, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(absent == 1);
}

void corrupt_manifest(const runtime::MpiContext &mpi,
                      const std::filesystem::path &directory) {
  if (mpi.rank() == 0) {
    auto bytes = read_bytes(directory / "manifest.v3.bin");
    HUNDUN_CHECK(!bytes.empty());
    bytes.back() ^= 0x1U;
    test::write_bytes(directory / "manifest.v3.bin", bytes);
  }
  mpi.barrier();
}

void run_profile(const runtime::MpiContext &mpi, bool immersed) {
  CaseContext context(mpi, immersed);
  const auto config = checkpoint_case(mpi.size(), immersed);
  FlowBundle continuous(context, mpi);
  auto continuous_state = context.make_state();
  check_module_rejection_before_io(mpi, context, continuous, config,
                                   continuous_state);

  const flow::ImmersedFlowPhysics physics{
      config::DensityModel::constant, kRho, kMu, std::nullopt, std::nullopt,
      std::nullopt};
  const auto first_stencil = flow::make_momentum_time_stencil(
      flow::MomentumTimeOrder::backward_euler, kDt, 0.0);
  const auto second_stencil = flow::make_momentum_time_stencil(
      flow::MomentumTimeOrder::bdf2, kDt, kDt);
  linear::SolveControl solve;
  solve.max_iterations = 5000U;
  solve.residual_recompute_interval = 20U;

  const auto first = continuous.facade->attempt(
      continuous_state, physics, first_stencil, solve, solve);
  static_cast<void>(require_committed(first, immersed));
  const auto checkpoint_directory = context.shared_path(
      mpi, immersed ? "profile-3" : "profile-2");
  const auto written = flow::write_checkpoint_v3(
      mpi, context.decomposition, context.topology, context.geometry,
      context.boundaries, config, write_modules(context, continuous),
      continuous_state, {kDt, 0U}, checkpoint_directory);
  HUNDUN_CHECK(written.disposition() ==
               flow::CheckpointV3Disposition::completed);
  HUNDUN_CHECK(written.presence() ==
               write_modules(context, continuous).presence);
  HUNDUN_CHECK(std::filesystem::exists(checkpoint_directory / "COMPLETED"));
  check_manifest_omits_attempt_identity(
      mpi, checkpoint_directory, first.wale->identity.value);

  FlowBundle restored(context, mpi);
  auto restored_state = context.make_state();
  HUNDUN_CHECK(flow::test::ImmersedFlowTestAccess::
                   wale_coefficient_identity(*restored.facade)
                   .value == 0U);
  auto incompatible_config = config;
  incompatible_config.common_flow.physics.dynamic_viscosity_pa_s *= 2.0;
  const auto before_incompatible_read =
      test::snapshot_stage3_state(restored_state);
  const auto wale_before_incompatible_read =
      flow::test::ImmersedFlowTestAccess::wale_coefficient_identity(
          *restored.facade);
  const auto incompatible = flow::read_checkpoint_v3(
      mpi, context.decomposition, context.topology, context.geometry,
      context.boundaries, incompatible_config,
      read_modules(context, restored), restored_state, checkpoint_directory);
  HUNDUN_CHECK(!incompatible.restored());
  HUNDUN_CHECK(incompatible.report().reason() ==
               flow::CheckpointV3FailureReason::fingerprint);
  HUNDUN_CHECK(incompatible.report().fingerprint_status() ==
               flow::CheckpointV3CheckStatus::failed);
  HUNDUN_CHECK(incompatible.report().rollback_status() ==
               flow::CheckpointV3CheckStatus::passed);
  HUNDUN_CHECK(test::stage3_state_bitwise_equal(
      before_incompatible_read,
      test::snapshot_stage3_state(restored_state)));
  HUNDUN_CHECK(
      flow::test::ImmersedFlowTestAccess::wale_coefficient_identity(
          *restored.facade) == wale_before_incompatible_read);

  const auto restored_result = flow::read_checkpoint_v3(
      mpi, context.decomposition, context.topology, context.geometry,
      context.boundaries, config, read_modules(context, restored),
      restored_state, checkpoint_directory);
  HUNDUN_CHECK(restored_result.restored());
  HUNDUN_CHECK(restored_result.report().presence() == written.presence());
  HUNDUN_CHECK(restored_result.control_state().proposed_next_dt_s == kDt);
  HUNDUN_CHECK(same_accepted_state(
      test::snapshot_stage3_state(continuous_state),
      test::snapshot_stage3_state(restored_state)));
  HUNDUN_CHECK(flow::test::ImmersedFlowTestAccess::
                   wale_coefficient_identity(*restored.facade)
                   .value == 0U);

  const auto continuous_second = continuous.facade->attempt(
      continuous_state, physics, second_stencil, solve, solve);
  const auto restored_second = restored.facade->attempt(
      restored_state, physics, second_stencil, solve, solve);
  static_cast<void>(require_committed(continuous_second, immersed));
  static_cast<void>(require_committed(restored_second, immersed));
  HUNDUN_CHECK(same_wale_summary(*continuous_second.wale,
                                 *restored_second.wale));
  HUNDUN_CHECK(test::stage3_state_bitwise_equal(
      test::snapshot_stage3_state(continuous_state),
      test::snapshot_stage3_state(restored_state)));

  const auto before_failed_read = test::snapshot_stage3_state(restored_state);
  const auto wale_before_failed_read =
      flow::test::ImmersedFlowTestAccess::wale_coefficient_identity(
          *restored.facade);
  corrupt_manifest(mpi, checkpoint_directory);
  const auto failed = flow::read_checkpoint_v3(
      mpi, context.decomposition, context.topology, context.geometry,
      context.boundaries, config, read_modules(context, restored),
      restored_state, checkpoint_directory);
  HUNDUN_CHECK(!failed.restored());
  HUNDUN_CHECK(failed.report().reason() ==
               flow::CheckpointV3FailureReason::file_integrity);
  HUNDUN_CHECK(failed.report().rollback_status() ==
               flow::CheckpointV3CheckStatus::passed);
  HUNDUN_CHECK(test::stage3_state_bitwise_equal(
      before_failed_read, test::snapshot_stage3_state(restored_state)));
  HUNDUN_CHECK(
      flow::test::ImmersedFlowTestAccess::wale_coefficient_identity(
          *restored.facade) == wale_before_failed_read);
}

void run() {
  auto mpi = runtime::MpiContext::duplicate(MPI_COMM_WORLD);
  HUNDUN_CHECK(mpi.size() == 1 || mpi.size() == 2);
  HUNDUN_CHECK(test::stage3_state_equality_oracle_is_mutation_sensitive());
  run_profile(mpi, false);
  run_profile(mpi, true);
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  if (argc != 1) {
    return 2;
  }
  return hundun::test::run(run);
}
