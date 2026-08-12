// SPDX-License-Identifier: Apache-2.0

#include "tests/support/flow_ideal_gas_closure_test_access.hpp"
#include "tests/support/flow_immersed_test_access.hpp"
#include "tests/support/flow_material_density_piso_test_access.hpp"
#include "tests/support/flow_material_density_transport_test_access.hpp"
#include "tests/support/flow_state_equality.hpp"
#include "tests/support/stage3_mms.hpp"
#include "tests/support/stage3_stl_fixture.hpp"
#include "tests/support/test_main.hpp"

#include "hundun/bc_basic_boundary.hpp"
#include "hundun/diag_immersed_module.hpp"
#include "hundun/exec_execution.hpp"
#include "hundun/flow_adaptive_time_control.hpp"
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
#include "hundun/rt_halo_exchange.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "hundun/rt_structured_decomposition.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <limits>
#include <optional>
#include <string>
#include <utility>

namespace {

using namespace hundun;

constexpr int kCells = 8;
constexpr int kOpenCells = 8;
constexpr double kDt = 1.0e-4;
constexpr double kMu = 0.01;
constexpr double kCp = 1000.0;
constexpr double kGasConstant = 287.05;
constexpr double kConfiguredPressure = 101325.0;

runtime::Int3 process_grid(int ranks) {
  if (ranks == 1)
    return {1, 1, 1};
  if (ranks == 2)
    return {2, 1, 1};
  throw runtime::Error("unsupported immersed ideal-gas rank count");
}

config::FlowCaseConfig ideal_case(int ranks, bool open) {
  config::FlowCaseConfig result{};
  result.schema_version = 2;
  result.case_name = open ? "stage3-d2-immersed-ideal-open"
                          : "stage3-d2-immersed-ideal-closed";
  result.simulation_type = config::SimulationType::variable_density_flow;
  result.density_model = config::DensityModel::ideal_gas;
  result.resources.expected_ranks = ranks;
  result.resources.process_grid = process_grid(ranks);
  const int cells = open ? kOpenCells : kCells;
  result.mesh.cells = {cells, cells, cells};
  result.mesh.origin_m = {0.0, 0.0, 0.0};
  result.mesh.length_m = {1.0, 1.0, 1.0};
  result.physics.rho_ref_kg_per_m3 = 1.0;
  result.physics.dynamic_viscosity_pa_s = open ? 0.0 : kMu;
  result.physics.inlet_consistency_rtol = 1.0e-12;
  result.physics.cp_J_per_kg_K = kCp;
  result.physics.gas_constant_J_per_kg_K = kGasConstant;
  result.physics.thermodynamic_pressure_pa = kConfiguredPressure;
  constexpr std::array<config::PatchName, 6> names{
      config::PatchName::x_min, config::PatchName::x_max,
      config::PatchName::y_min, config::PatchName::y_max,
      config::PatchName::z_min, config::PatchName::z_max};
  for (std::size_t patch = 0U; patch < names.size(); ++patch) {
    result.boundaries[patch].patch = names[patch];
    result.boundaries[patch].type = config::BoundaryType::periodic;
  }
  if (open) {
    result.boundaries[0].type = config::BoundaryType::velocity_inlet;
    result.boundaries[0].velocity_m_per_s = {0.0, 0.0, 0.0};
    result.boundaries[0].thermal_authority =
        config::InletThermalAuthority::temperature;
    result.boundaries[0].temperature_K = 300.0;
    result.boundaries[0].enthalpy_J_per_kg = kCp * 300.0;
    result.boundaries[0].density_kg_per_m3 =
        kConfiguredPressure / (kGasConstant * 300.0);
    result.boundaries[0].scalar_values =
        std::vector<config::InletScalarValue>{};
    result.boundaries[1].type = config::BoundaryType::pressure_outlet;
    result.boundaries[1].pressure_perturbation_pa = 0.0;
    for (std::size_t patch = 2U; patch < names.size(); ++patch)
      result.boundaries[patch].type = config::BoundaryType::symmetry;
  }
  return result;
}

runtime::FieldDescriptor cell_field(const char *name, const char *unit,
                                    std::uint32_t components = 1U,
                                    bool conservative = true) {
  return {name,
          unit,
          "stage3_d2_immersed_ideal",
          runtime::FunctionSpace::cell_average,
          runtime::ScalarType::float64,
          components,
          4,
          conservative,
          runtime::RestartPolicy::persistent,
          runtime::OutputPolicy::never};
}

runtime::FieldDescriptor face_field(const char *name,
                                    std::uint32_t components) {
  return {name,
          "1",
          "stage3_d2_immersed_ideal",
          runtime::FunctionSpace::face_value,
          runtime::ScalarType::float64,
          components,
          0,
          false,
          runtime::RestartPolicy::persistent,
          runtime::OutputPolicy::never};
}

void check_mpi(int code, const char *operation) {
  if (code != MPI_SUCCESS)
    throw runtime::Error(std::string(operation) + " failed");
}

std::string collective_surface_path(
    const runtime::MpiContext &mpi, bool open,
    std::optional<test::Stage3TemporaryDirectory> &directory) {
  std::string path_text;
  if (mpi.rank() == 0) {
    directory.emplace(open ? "stage3-d2-ideal-open"
                           : "stage3-d2-ideal-closed");
    const auto body = test::stage3::approved_body(test::stage3::BodyKind::sphere);
    const int cells = open ? kOpenCells : kCells;
    const auto surface = test::stage3::make_manufactured_surface(
        body, 1.0 / static_cast<double>(cells));
    const auto path = directory->path() / "sphere.stl";
    test::write_text(path, test::ascii_stl(surface.triangles, "sphere"));
    path_text = path.string();
  }
  std::uint64_t size = path_text.size();
  check_mpi(MPI_Bcast(&size, 1, MPI_UINT64_T, 0, mpi.comm()),
            "MPI_Bcast(immersed ideal-gas path size)");
  HUNDUN_CHECK(size <=
               static_cast<std::uint64_t>(std::numeric_limits<int>::max()));
  path_text.resize(static_cast<std::size_t>(size));
  check_mpi(MPI_Bcast(path_text.data(), static_cast<int>(path_text.size()),
                      MPI_BYTE, 0, mpi.comm()),
            "MPI_Bcast(immersed ideal-gas path)");
  return path_text;
}

double relative_error(double left, double right) noexcept {
  return std::abs(left - right) /
         std::max({std::abs(left), std::abs(right),
                   std::numeric_limits<double>::min()});
}

double temperature(runtime::Real3 point) noexcept {
  const auto body = test::stage3::approved_body(test::stage3::BodyKind::sphere);
  const double base = point.x * (1.0 - point.x) * point.y * (1.0 - point.y) *
                      point.z * (1.0 - point.z);
  const double envelope = base * base;
  const double dx = point.x - body.centre_m.x;
  const double dy = point.y - body.centre_m.y;
  const double dz = point.z - body.centre_m.z;
  const double body_factor =
      dx * dx + dy * dy + dz * dz - body.radius_m * body.radius_m;
  return 300.0 + 20.0 * envelope * body_factor * body_factor;
}

flow::FlowLayerValues initial_values(
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    const immersed::ImmersedDomain &domain, bool open) {
  flow::FlowLayerValues result;
  result.density.assign(topology.owned_cell_count(), 0.0);
  result.velocity.assign(topology.owned_cell_count() * 3U, 0.0);
  result.mechanical_pressure.assign(topology.owned_cell_count(), 0.0);
  result.face_velocity.assign(topology.local_face_count() * 3U, 0.0);
  result.face_mass_flux.assign(topology.local_face_count(), 0.0);
  result.transported_cell_fields = {
      std::vector<double>(topology.owned_cell_count(), 0.0)};
  for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count(); ++cell) {
    if (domain.region(cell) != immersed::CellRegion::fluid)
      continue;
    const double t = open ? 300.0 : temperature(geometry.cell_center_m(cell));
    const double rho = kConfiguredPressure / (kGasConstant * t);
    result.density[cell] = rho;
    result.transported_cell_fields.front()[cell] = rho * kCp * t;
  }
  return result;
}

struct ExactState final {
  flow::FlowLayerValues history;
  flow::FlowLayerValues committed;
  flow::FlowLayerValues trial;
  flow::AcceptedStepMetadata metadata;
};

ExactState capture(const flow::FlowState &state) {
  return {state.snapshot(flow::FlowLayer::history),
          state.snapshot(flow::FlowLayer::committed),
          state.snapshot(flow::FlowLayer::trial), state.metadata()};
}

void require_exact(const ExactState &expected, const flow::FlowState &state) {
  HUNDUN_CHECK(test::flow_layer_values_bitwise_equal(
      expected.history, state.snapshot(flow::FlowLayer::history)));
  HUNDUN_CHECK(test::flow_layer_values_bitwise_equal(
      expected.committed, state.snapshot(flow::FlowLayer::committed)));
  HUNDUN_CHECK(test::flow_layer_values_bitwise_equal(
      expected.trial, state.snapshot(flow::FlowLayer::trial)));
  HUNDUN_CHECK(test::accepted_step_metadata_bitwise_equal(expected.metadata,
                                                          state.metadata()));
}

std::optional<std::size_t>
first_active_cell(const immersed::ImmersedDomain &domain,
                  const mesh::MeshTopology &topology) {
  for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count(); ++cell)
    if (domain.region(cell) == immersed::CellRegion::fluid)
      return static_cast<std::size_t>(cell);
  return std::nullopt;
}

void run_domain(bool open, bool transaction_checks) {
  auto mpi = runtime::MpiContext::duplicate(MPI_COMM_WORLD);
  HUNDUN_CHECK(mpi.size() == 1 || mpi.size() == 2);
  const int cells = open ? kOpenCells : kCells;
  const double dt = open ? 1.0e-6 : kDt;
  const auto decomposition = runtime::StructuredDecomposition::create(
      mpi, {cells, cells, cells},
      open ? std::array<bool, 3>{false, false, false}
           : std::array<bool, 3>{true, true, true},
      runtime::DecompositionOptions{process_grid(mpi.size())});
  const mesh::MeshTopology topology(decomposition);
  const mesh::MeshGeometry geometry(
      topology, mesh::UniformBoxMapping({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}));
  const auto boundaries =
      boundary::BoundaryRegistry::create(ideal_case(mpi.size(), open), topology);

  std::optional<test::Stage3TemporaryDirectory> surface_directory;
  const auto surface_path =
      collective_surface_path(mpi, open, surface_directory);
  const auto surface = immersed::ImmersedSurface::load_collective(
      std::filesystem::path(surface_path), 1.0, mpi, 0);
  const auto query = immersed::SurfaceQuery::create(surface);
  const auto domain = immersed::ImmersedDomain::create(
      surface, query, config::ImmersedFluidSide::outside, topology, geometry,
      boundaries, mpi);
  const auto ghost_plan = immersed::GhostStencilPlan::create(
      surface, query, domain, topology, geometry, decomposition, mpi);
  const auto wall_plan = immersed::WallQuadraturePlan::create(
      surface, query, domain, topology, geometry, mpi);
  HUNDUN_CHECK(ghost_plan.maximum_halo_reach() <= 4U);

  runtime::FieldRegistry registry;
  flow::FlowFieldIds fields;
  fields.density = registry.declare_field(cell_field("rho", "kg/m3"));
  fields.velocity =
      registry.declare_field(cell_field("velocity", "m/s", 3U, false));
  fields.mechanical_pressure =
      registry.declare_field(cell_field("pi", "Pa", 1U, false));
  fields.face_velocity =
      registry.declare_field(face_field("face_velocity", 3U));
  fields.face_mass_flux = finite_volume::declare_face_mass_flux(registry);
  const auto rho_h =
      registry.declare_field(cell_field("rho_h", "J/m3"));
  fields.transported_cell_fields = {rho_h};
  registry.freeze();
  const flow::MaterialDensityTransportSpec material_spec{rho_h, 0.0, {}, {}};
  const auto initial = initial_values(topology, geometry, domain, open);
  const auto make_state = [&] {
    auto state = flow::FlowState::create(
        registry,
        {decomposition.local_extent(), topology.local_face_count()}, fields,
        {0U, 0.0, dt, 0.0, flow::MomentumTimeOrder::backward_euler});
    state.seed_accepted_layers(initial, initial);
    return state;
  };
  auto state = make_state();
  auto closure = flow::test::IdealGasClosureTestAccess::create_immersed(
      topology, geometry, boundaries, domain, mpi, registry, fields, state,
      {rho_h, kCp, kGasConstant, kConfiguredPressure});

  execution::CpuReferenceContext execution;
  auto halo = runtime::HaloExchange::create(
      decomposition, runtime::ExchangePlan::create(
                         decomposition, decomposition.local_extent(), 4));
  linear::ConjugateGradientSolver momentum_solver(execution, mpi);
  linear::BiCGStabSolver pressure_solver(execution, mpi);
  linear::JacobiPreconditioner mx(execution), my(execution), mz(execution),
      pressure_preconditioner(execution);
  immersed::LocalFlowPatternTransform transform;
  const flow::ImmersedFlowDensitySetup density_setup{
      config::DensityModel::ideal_gas, &registry, fields, material_spec,
      &closure};
  auto facade = flow::FixedStepImmersedFlow::create(
      decomposition, topology, geometry, boundaries, &domain, &ghost_plan,
      &wall_plan, &transform, nullptr, density_setup, mpi, execution, halo,
      momentum_solver, {&mx, &my, &mz}, pressure_solver,
      pressure_preconditioner);
  const flow::ImmersedFlowPhysics physics{
      config::DensityModel::ideal_gas, 1.0, open ? 0.0 : kMu, kCp,
      kGasConstant,
      kConfiguredPressure};
  const auto stencil = flow::make_momentum_time_stencil(
      flow::MomentumTimeOrder::backward_euler, dt, 0.0);
  linear::SolveControl control;
  control.residual_recompute_interval = 20U;

  if (open) {
    state.begin_attempt();
    flow::test::IdealGasClosureTestAccess::begin_attempt(closure, state,
                                                         0x53334432U);
    for (const auto stage : {flow::IdealGasClosureStage::predictor,
                             flow::IdealGasClosureStage::provisional,
                             flow::IdealGasClosureStage::final}) {
      const auto report =
          flow::test::IdealGasClosureTestAccess::evaluate(closure, state,
                                                          stage);
      HUNDUN_CHECK(report.disposition() ==
                   flow::IdealGasClosureDisposition::closed);
      HUNDUN_CHECK(test::fp64_bits(report.candidate_pressure_pa()) ==
                   test::fp64_bits(kConfiguredPressure));
    }
    const auto old = state.metadata();
    flow::test::MaterialDensityPisoTestAccess::prepare_state_commit(
        state,
        {old.step + 1U, old.time_s + dt, dt, old.dt_s,
         flow::MomentumTimeOrder::backward_euler});
    HUNDUN_CHECK(
        flow::test::IdealGasClosureTestAccess::prepare_commit(closure) < 0);
    flow::test::MaterialDensityPisoTestAccess::publish_state_commit(state);
    flow::test::IdealGasClosureTestAccess::publish_commit(closure);
    HUNDUN_CHECK(closure.state().revision == 1U);
    HUNDUN_CHECK(test::fp64_bits(closure.state().thermodynamic_pressure_pa) ==
                 test::fp64_bits(kConfiguredPressure));
  }
  if (transaction_checks) {
    const auto state_before = capture(state);
    const auto closure_before = closure.state();
    auto mismatched_physics = physics;
    mismatched_physics.thermodynamic_pressure_pa =
        std::nextafter(kConfiguredPressure,
                       std::numeric_limits<double>::infinity());
    const auto rejected =
        facade.attempt(state, mismatched_physics, stencil, control, control);
    HUNDUN_CHECK(
        std::holds_alternative<flow::IdealGasStepAttemptReport>(rejected.base));
    const auto &rejected_ideal =
        std::get<flow::IdealGasStepAttemptReport>(rejected.base);
    HUNDUN_CHECK(flow::test::IdealGasClosureTestAccess::report_authenticated(
        rejected_ideal));
    HUNDUN_CHECK(rejected_ideal.flow().flow().disposition ==
                 flow::StepAttemptDisposition::non_retryable_failure);
    HUNDUN_CHECK(rejected_ideal.flow().flow().reason ==
                 flow::StepFailureReason::invalid_input);
    HUNDUN_CHECK(!rejected.force.has_value());
    HUNDUN_CHECK(!rejected.wale.has_value());
    require_exact(state_before, state);
    HUNDUN_CHECK(test::ideal_gas_closure_state_bitwise_equal(
        closure_before, closure.state()));
  }
  const auto before_facade_attempt = capture(state);
  const auto accepted = facade.attempt(state, physics, stencil, control, control);
  if (!std::holds_alternative<flow::IdealGasStepAttemptReport>(accepted.base)) {
    if (mpi.rank() == 0)
      std::cerr << "IMMERSED_IDEAL_GAS_UNSUPPORTED\n";
    HUNDUN_CHECK(false);
  }
  const auto &ideal = std::get<flow::IdealGasStepAttemptReport>(accepted.base);
  HUNDUN_CHECK(flow::test::IdealGasClosureTestAccess::report_authenticated(ideal));
  const bool committed_attempt =
      ideal.flow().flow().disposition ==
      flow::StepAttemptDisposition::committed;
  const bool expected_open_backflow =
      open && ideal.flow().flow().disposition ==
                  flow::StepAttemptDisposition::recoverable_failure &&
      ideal.flow().flow().reason == flow::StepFailureReason::boundary_backflow;
  if (!committed_attempt && !expected_open_backflow &&
      mpi.rank() == 0) {
    const auto &failed = ideal.flow().flow();
    std::cerr << "IMMERSED_IDEAL_GAS_FAILURE reason="
              << static_cast<int>(failed.reason)
              << " rank=" << failed.lowest_failing_rank
              << " correctors=" << failed.pressure_corrector_count
              << " momentum=" << failed.final_momentum_normalized_l2[0]
              << ',' << failed.final_momentum_normalized_l2[1] << ','
              << failed.final_momentum_normalized_l2[2]
              << " conservation="
              << failed.final_momentum_relative_conservation_defect[0] << ','
              << failed.final_momentum_relative_conservation_defect[1] << ','
              << failed.final_momentum_relative_conservation_defect[2]
              << '\n';
  }
  HUNDUN_CHECK(committed_attempt || expected_open_backflow);
  if (expected_open_backflow)
    require_exact(before_facade_attempt, state);
  HUNDUN_CHECK(ideal.flow().flow().pressure_corrector_count == 2U);
  HUNDUN_CHECK(ideal.closure_report_available());
  HUNDUN_CHECK(ideal.closure_report().stage() == flow::IdealGasClosureStage::final);
  HUNDUN_CHECK(ideal.closure_report().final_metrics_available());
  HUNDUN_CHECK(ideal.closure_report().enthalpy_temperature_max_relative_error() <=
               1.0e-12);
  HUNDUN_CHECK(ideal.closure_report().eos_max_relative_error() <= 1.0e-12);
  HUNDUN_CHECK(flow::test::IdealGasClosureTestAccess::
                   post_eos_evidence_authenticated(ideal.flow()));
  HUNDUN_CHECK(accepted.force.has_value() == committed_attempt);
  HUNDUN_CHECK(!accepted.wale.has_value());

  const auto committed = state.snapshot(flow::FlowLayer::committed);
  double local_mass = 0.0;
  double local_denominator = 0.0;
  double local_wrong_denominator = 0.0;
  double maximum_h_t_error = 0.0;
  double maximum_eos_error = 0.0;
  double maximum_pi = 0.0;
  for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count(); ++cell) {
    const double volume = geometry.cell_volume_m3(cell);
    if (domain.region(cell) != immersed::CellRegion::fluid) {
      HUNDUN_CHECK(test::fp64_bits(committed.density[cell]) ==
                   test::fp64_bits(0.0));
      HUNDUN_CHECK(test::fp64_bits(committed.transported_cell_fields[0][cell]) ==
                   test::fp64_bits(0.0));
      local_wrong_denominator += volume / 300.0;
      continue;
    }
    const double rho = committed.density[cell];
    const double h = committed.transported_cell_fields[0][cell] / rho;
    const double t = h / kCp;
    HUNDUN_CHECK(rho > 0.0 && t > 0.0 && std::isfinite(rho) &&
                 std::isfinite(t));
    local_mass += volume * rho;
    local_denominator += volume / t;
    local_wrong_denominator += volume / t;
    maximum_h_t_error =
        std::max(maximum_h_t_error, relative_error(h, kCp * t));
    maximum_eos_error = std::max(
        maximum_eos_error,
        relative_error(rho * kGasConstant * t,
                       closure.state().thermodynamic_pressure_pa));
    maximum_pi = std::max(maximum_pi,
                          std::abs(committed.mechanical_pressure[cell]));
  }
  double sums[3]{local_mass, local_denominator, local_wrong_denominator};
  mpi.allreduce_fp64_in_place(sums, 3U,
                              runtime::Fp64ReductionOperation::sum);
  double maxima[3]{maximum_h_t_error, maximum_eos_error, maximum_pi};
  mpi.allreduce_fp64_in_place(maxima, 3U,
                              runtime::Fp64ReductionOperation::maximum);
  const auto closure_state = closure.state();
  if (open) {
    HUNDUN_CHECK(closure_state.mode == flow::IdealGasPressureMode::open_fixed);
    HUNDUN_CHECK(!closure_state.target_mass_kg.has_value());
    HUNDUN_CHECK(test::fp64_bits(closure_state.thermodynamic_pressure_pa) ==
                 test::fp64_bits(kConfiguredPressure));
    HUNDUN_CHECK(maxima[2] <= 1.0e-8 * kConfiguredPressure);
  } else {
    HUNDUN_CHECK(closure_state.mode ==
                 flow::IdealGasPressureMode::closed_dynamic);
    HUNDUN_CHECK(closure_state.target_mass_kg.has_value());
    HUNDUN_CHECK(relative_error(sums[0], *closure_state.target_mass_kg) <=
                 5.0e-12);
    const double expected_pressure =
        *closure_state.target_mass_kg * kGasConstant / sums[1];
    HUNDUN_CHECK(relative_error(closure_state.thermodynamic_pressure_pa,
                                expected_pressure) <= 5.0e-12);
    const double inactive_volume_pressure =
        *closure_state.target_mass_kg * kGasConstant / sums[2];
    HUNDUN_CHECK(relative_error(closure_state.thermodynamic_pressure_pa,
                                inactive_volume_pressure) > 1.0e-6);
  }
  HUNDUN_CHECK(maxima[0] <= 1.0e-12);
  HUNDUN_CHECK(maxima[1] <= 1.0e-12);

  if (!transaction_checks)
    return;

  static_cast<void>(facade.diagnostic_source(state, accepted));
  for (std::uint16_t mutation = 0U;
       mutation < static_cast<std::uint16_t>(
                      flow::test::IdealGasPostEvidenceMutation::count);
       ++mutation) {
    HUNDUN_CHECK(flow::test::IdealGasClosureTestAccess::
                     post_evidence_mutation_rejected(
                         ideal,
                         static_cast<flow::test::IdealGasPostEvidenceMutation>(
                             mutation)));
  }
  for (std::uint8_t field = 0U; field < 28U; ++field)
    HUNDUN_CHECK(flow::test::IdealGasClosureTestAccess::
                     closure_field_mutation_rejected(ideal, field));
  for (const auto change :
       {flow::test::IdealGasStepReportCorruption::nested_material,
        flow::test::IdealGasStepReportCorruption::closure_stage,
        flow::test::IdealGasStepReportCorruption::closure_candidate_pressure,
        flow::test::IdealGasStepReportCorruption::closure_seal,
        flow::test::IdealGasStepReportCorruption::outer_attempt_identity,
        flow::test::IdealGasStepReportCorruption::outer_closure_presence,
        flow::test::IdealGasStepReportCorruption::outer_seal}) {
    auto changed = accepted;
    auto &changed_ideal =
        std::get<flow::IdealGasStepAttemptReport>(changed.base);
    flow::test::IdealGasClosureTestAccess::corrupt_report(changed_ideal,
                                                          change);
    HUNDUN_CHECK(
        !flow::test::IdealGasClosureTestAccess::report_authenticated(
            changed_ideal));
    bool rejected = false;
    try {
      static_cast<void>(facade.diagnostic_source(state, changed));
    } catch (const runtime::Error &) {
      rejected = true;
    }
    HUNDUN_CHECK(rejected);
  }

  const config::FlowTimeConfig time{
      config::TimeMode::adaptive, 2, kDt, kDt / 8.0, kDt, 0.5, 0.25,
      1.25, 0.5, 8};
  const int injected_rank = mpi.size() > 1 ? 1 : 0;
  const auto require_failed_attempt_neutral = [&](flow::FlowState &candidate,
                                                   flow::Bdf2RetryController &
                                                       controller) {
    const auto before = capture(candidate);
    const auto closure_before = closure.state();
    const auto controller_before = controller.state();
    const auto failed =
        facade.attempt(candidate, physics, stencil, control, control);
    HUNDUN_CHECK(
        std::holds_alternative<flow::IdealGasStepAttemptReport>(failed.base));
    const auto &failed_ideal =
        std::get<flow::IdealGasStepAttemptReport>(failed.base);
    HUNDUN_CHECK(flow::test::IdealGasClosureTestAccess::report_authenticated(
        failed_ideal));
    HUNDUN_CHECK(failed_ideal.flow().flow().disposition !=
                 flow::StepAttemptDisposition::committed);
    HUNDUN_CHECK(!failed.force.has_value());
    HUNDUN_CHECK(!failed.wale.has_value());
    require_exact(before, candidate);
    HUNDUN_CHECK(test::ideal_gas_closure_state_bitwise_equal(
        closure_before, closure.state()));
    HUNDUN_CHECK(test::time_control_state_bitwise_equal(controller_before,
                                                        controller.state()));
    return failed_ideal.flow().flow().lowest_failing_rank;
  };

  {
    auto invalid_temperature = make_state();
    const auto cell = first_active_cell(domain, topology);
    HUNDUN_CHECK(cell.has_value());
    if (mpi.rank() == injected_rank)
      flow::test::MaterialDensityTransportTestAccess::
          set_accepted_transport_value(invalid_temperature, false, 0U, *cell,
                                       0.0);
    auto controller = flow::Bdf2RetryController::create(
        time, config::DensityModel::ideal_gas, topology, geometry, mpi,
        invalid_temperature);
    HUNDUN_CHECK(require_failed_attempt_neutral(invalid_temperature,
                                                controller) == injected_rank);
  }
  {
    auto invalid_density = make_state();
    const auto cell = first_active_cell(domain, topology);
    HUNDUN_CHECK(cell.has_value());
    if (mpi.rank() == injected_rank)
      flow::test::MaterialDensityTransportTestAccess::
          set_accepted_density_value(invalid_density, false, *cell, 0.0);
    auto controller = flow::Bdf2RetryController::create(
        time, config::DensityModel::ideal_gas, topology, geometry, mpi,
        invalid_density);
    HUNDUN_CHECK(require_failed_attempt_neutral(invalid_density, controller) ==
                 injected_rank);
  }
  {
    auto mismatch = make_state();
    auto controller = flow::Bdf2RetryController::create(
        time, config::DensityModel::ideal_gas, topology, geometry, mpi,
        mismatch);
    flow::test::IdealGasClosureTestAccess::set_attempt_layout_fault(
        closure, injected_rank);
    HUNDUN_CHECK(require_failed_attempt_neutral(mismatch, controller) ==
                 injected_rank);
  }
  {
    auto post_pressure_failure = make_state();
    auto controller = flow::Bdf2RetryController::create(
        time, config::DensityModel::ideal_gas, topology, geometry, mpi,
        post_pressure_failure);
    flow::test::ImmersedFlowTestAccess::set_failure_stage(
        flow::test::ImmersedFlowAttemptFailureStage::after_final_transport,
        injected_rank);
    const int failed_rank =
        require_failed_attempt_neutral(post_pressure_failure, controller);
    flow::test::ImmersedFlowTestAccess::clear_failure_stage();
    HUNDUN_CHECK(failed_rank == injected_rank);
  }
  {
    auto before_commit_failure = make_state();
    auto controller = flow::Bdf2RetryController::create(
        time, config::DensityModel::ideal_gas, topology, geometry, mpi,
        before_commit_failure);
    flow::test::ImmersedFlowTestAccess::set_failure_stage(
        flow::test::ImmersedFlowAttemptFailureStage::before_commit,
        injected_rank);
    const int failed_rank =
        require_failed_attempt_neutral(before_commit_failure, controller);
    flow::test::ImmersedFlowTestAccess::clear_failure_stage();
    HUNDUN_CHECK(failed_rank == injected_rank);
  }
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  try {
#ifdef HUNDUN_STAGE3_IMMERSED_IDEAL_GAS_TRANSACTION
    run_domain(false, true);
    run_domain(true, false);
#else
    run_domain(false, false);
#endif
    return 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 8;
  }
}
