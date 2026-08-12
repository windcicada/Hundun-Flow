// SPDX-License-Identifier: Apache-2.0

#include "tests/support/flow_constant_density_piso_test_access.hpp"
#include "tests/support/flow_state_equality.hpp"
#include "tests/support/stage3_scientific_row.hpp"
#include "tests/support/test_main.hpp"

#include "hundun/bc_basic_boundary.hpp"
#include "hundun/exec_execution.hpp"
#include "hundun/flow_immersed.hpp"
#include "hundun/fvm_cell_centered.hpp"
#include "hundun/les_wale.hpp"
#include "hundun/lin_conjugate_gradient.hpp"
#include "hundun/lin_preconditioners.hpp"
#include "hundun/mesh_geometry.hpp"
#include "hundun/mesh_topology.hpp"
#include "hundun/rt_exchange_plan.hpp"
#include "hundun/rt_halo_exchange.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "hundun/rt_structured_decomposition.hpp"

#include <array>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

namespace {

using namespace hundun;
runtime::Int3 selected_extent{8, 8, 4};
bool screen_only{};
bool formal_scientific{};
constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kDt = 1.0e-4;

std::uint64_t fp64_bits(double value) noexcept {
  std::uint64_t result{};
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

std::uint64_t mix(std::uint64_t hash, std::uint64_t value) noexcept {
  constexpr std::uint64_t prime = UINT64_C(1099511628211);
  for (std::size_t byte = 0U; byte < sizeof(value); ++byte) {
    hash ^= (value >> (byte * 8U)) & UINT64_C(0xff);
    hash *= prime;
  }
  return hash;
}

std::uint64_t field_fingerprint(const flow::FlowLayerValues &values) {
  std::uint64_t hash = UINT64_C(1469598103934665603);
  const auto append = [&](const std::vector<double> &field) {
    hash = mix(hash, field.size());
    for (const double value : field)
      hash = mix(hash, fp64_bits(value));
  };
  append(values.density);
  append(values.velocity);
  append(values.mechanical_pressure);
  append(values.face_velocity);
  append(values.face_mass_flux);
  hash = mix(hash, values.transported_cell_fields.size());
  for (const auto &field : values.transported_cell_fields)
    append(field);
  return hash;
}

std::uint64_t report_fingerprint(
    const flow::StepAttemptReport &report,
    const les::WaleSummary &wale) noexcept {
  std::uint64_t hash = UINT64_C(1469598103934665603);
  const auto append_solve = [&](const linear::SolveReport &solve) {
    hash = mix(hash, static_cast<std::uint64_t>(solve.reason));
    hash = mix(hash, solve.iterations);
    hash = mix(hash, fp64_bits(solve.initial_residual));
    hash = mix(hash, fp64_bits(solve.recursive_residual));
    hash = mix(hash, fp64_bits(solve.final_residual));
    hash = mix(hash, solve.matvec_count);
    hash = mix(hash, solve.preconditioner_apply_count);
    hash = mix(hash, solve.global_reduction_count);
    hash = mix(hash,
               static_cast<std::uint64_t>(solve.lowest_failing_rank + 1));
  };
  hash = mix(hash, static_cast<std::uint64_t>(report.disposition));
  hash = mix(hash, static_cast<std::uint64_t>(report.reason));
  hash = mix(hash,
             static_cast<std::uint64_t>(report.lowest_failing_rank + 1));
  hash = mix(hash, report.pressure_corrector_count);
  hash = mix(hash, fp64_bits(report.attempted_dt_s));
  hash = mix(hash, fp64_bits(report.suggested_dt_s));
  for (const auto &solve : report.momentum.components)
    append_solve(solve);
  for (const auto &solve : report.pressure)
    append_solve(solve);
  hash = mix(hash, fp64_bits(report.final_continuity_normalized_l2));
  hash = mix(hash, fp64_bits(report.final_pressure_residual_l2));
  for (const double value : report.final_momentum_normalized_l2)
    hash = mix(hash, fp64_bits(value));
  hash = mix(hash, report.final_transport_normalized_l2.size());
  for (const double value : report.final_transport_normalized_l2)
    hash = mix(hash, fp64_bits(value));
  hash = mix(hash, fp64_bits(report.final_mass_relative_conservation_defect));
  for (const double value : report.final_momentum_relative_conservation_defect)
    hash = mix(hash, fp64_bits(value));
  hash = mix(hash, report.final_transport_relative_conservation_defect.size());
  for (const double value : report.final_transport_relative_conservation_defect)
    hash = mix(hash, fp64_bits(value));
  hash = mix(hash, wale.identity.value);
  hash = mix(hash, fp64_bits(wale.minimum_nu_t_m2_per_s));
  hash = mix(hash, fp64_bits(wale.maximum_nu_t_m2_per_s));
  hash = mix(hash, fp64_bits(wale.l2_nu_t_m2_per_s));
  hash = mix(hash, wale.exact_zero_count);
  hash = mix(hash, wale.owned_active_count);
  return hash;
}

runtime::Int3 process_grid(int ranks) {
  return ranks == 1 ? runtime::Int3{1, 1, 1} : runtime::Int3{2, 1, 1};
}

config::FlowCaseConfig periodic_case() {
  config::FlowCaseConfig result{};
  result.schema_version = 2;
  result.simulation_type = config::SimulationType::variable_density_flow;
  result.density_model = config::DensityModel::constant;
  result.physics.rho_ref_kg_per_m3 = 1.0;
  result.physics.dynamic_viscosity_pa_s = 0.1;
  result.physics.inlet_consistency_rtol = 1.0e-12;
  result.scalars.push_back({"alpha", 0.0});
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
                                    std::uint32_t components) {
  return {name,
          "1",
          "task13-wale-body-fitted",
          runtime::FunctionSpace::cell_average,
          runtime::ScalarType::float64,
          components,
          2,
          false,
          runtime::RestartPolicy::persistent,
          runtime::OutputPolicy::never};
}

runtime::FieldDescriptor face_field(const char *name,
                                    std::uint32_t components) {
  return {name,
          "1",
          "task13-wale-body-fitted",
          runtime::FunctionSpace::face_value,
          runtime::ScalarType::float64,
          components,
          0,
          false,
          runtime::RestartPolicy::persistent,
          runtime::OutputPolicy::never};
}

flow::FlowLayerValues initial_values(const mesh::MeshTopology &topology,
                                     const mesh::MeshGeometry &geometry) {
  flow::FlowLayerValues result;
  result.density.assign(topology.owned_cell_count(), 1.0);
  result.velocity.resize(topology.owned_cell_count() * 3U);
  result.mechanical_pressure.resize(topology.owned_cell_count());
  result.face_velocity.resize(topology.local_face_count() * 3U);
  result.face_mass_flux.resize(topology.local_face_count());
  const auto velocity = [](runtime::Real3 point) {
    constexpr double perturbation = 0.1;
    return runtime::Real3{
        0.3 + std::sin(point.x) * std::cos(point.y) +
            perturbation * (std::sin(point.z) + std::cos(point.y)),
        -0.2 - std::cos(point.x) * std::sin(point.y) +
            perturbation * (std::sin(point.x) + std::cos(point.z)),
        0.25 + perturbation * (std::sin(point.y) + std::cos(point.x))};
  };
  for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count(); ++cell) {
    const auto value = velocity(geometry.cell_center_m(cell));
    result.velocity[cell * 3U] = value.x;
    result.velocity[cell * 3U + 1U] = value.y;
    result.velocity[cell * 3U + 2U] = value.z;
    const auto point = geometry.cell_center_m(cell);
    result.mechanical_pressure[cell] =
        0.25 * (std::cos(2.0 * point.x) + std::cos(2.0 * point.y));
  }
  for (mesh::LocalFaceId face = 0U; face < topology.local_face_count(); ++face) {
    const auto value = velocity(geometry.face_center_m(face));
    result.face_velocity[face * 3U] = value.x;
    result.face_velocity[face * 3U + 1U] = value.y;
    result.face_velocity[face * 3U + 2U] = value.z;
    const auto area =
        geometry.face_area_vector_m2(face, mesh::FaceSide::owner);
    result.face_mass_flux[face] =
        value.x * area.x + value.y * area.y + value.z * area.z;
  }
  return result;
}

void run() {
  auto mpi = runtime::MpiContext::duplicate(MPI_COMM_WORLD);
  HUNDUN_CHECK(mpi.size() == 1 || mpi.size() == 2);
  auto decomposition = runtime::StructuredDecomposition::create(
      mpi, selected_extent, {true, true, true},
      runtime::DecompositionOptions{process_grid(mpi.size())});
  mesh::MeshTopology topology(decomposition);
  mesh::MeshGeometry geometry(
      topology, mesh::UniformBoxMapping({0.0, 0.0, 0.0},
                                        {2.0 * kPi, 2.0 * kPi,
                                         2.0 * kPi}));
  auto boundaries =
      boundary::BoundaryRegistry::create(periodic_case(), topology);

  runtime::FieldRegistry registry;
  flow::FlowFieldIds fields;
  fields.density = registry.declare_field(cell_field("rho", 1U));
  fields.velocity = registry.declare_field(cell_field("u", 3U));
  fields.mechanical_pressure = registry.declare_field(cell_field("pi", 1U));
  fields.face_velocity = registry.declare_field(face_field("uf", 3U));
  fields.face_mass_flux = finite_volume::declare_face_mass_flux(registry);
  registry.freeze();
  const auto make_state = [&] {
    auto result = flow::FlowState::create(
        registry,
        {decomposition.local_extent(), topology.local_face_count()}, fields,
        {0U, 0.0, kDt, 0.0, flow::MomentumTimeOrder::backward_euler});
    const auto initial = initial_values(topology, geometry);
    result.seed_accepted_layers(initial, initial);
    return result;
  };

  std::vector<mesh::GlobalCellId> active;
  active.reserve(topology.local_cell_count());
  for (mesh::LocalCellId cell = 0U; cell < topology.local_cell_count(); ++cell)
    active.push_back(topology.global_cell_id(cell));
  execution::CpuReferenceContext execution;
  auto wale = les::WaleModel::create(
      {0.5, 0.9, 0.7}, topology, geometry, topology.owned_cell_count(), active,
      execution);
  auto halo = runtime::HaloExchange::create(
      decomposition,
      runtime::ExchangePlan::create(decomposition,
                                    decomposition.local_extent(), 2));
  linear::ConjugateGradientSolver momentum_solver(execution, mpi);
  linear::ConjugateGradientSolver pressure_solver(execution, mpi);
  linear::JacobiPreconditioner mx(execution), my(execution), mz(execution),
      pressure_preconditioner(execution);
  auto facade = flow::FixedStepImmersedFlow::create(
      decomposition, topology, geometry, boundaries, nullptr, nullptr, nullptr,
      nullptr, &wale, mpi, execution, halo, momentum_solver,
      {&mx, &my, &mz}, pressure_solver, pressure_preconditioner);
  const flow::ImmersedFlowPhysics physics{
      config::DensityModel::constant, 1.0, 0.1, std::nullopt, std::nullopt,
      std::nullopt};
  const auto stencil = flow::make_momentum_time_stencil(
      flow::MomentumTimeOrder::backward_euler, kDt, 0.0);

  auto molecular_facade = flow::FixedStepImmersedFlow::create(
      decomposition, topology, geometry, boundaries, nullptr, nullptr, nullptr,
      nullptr, nullptr, mpi, execution, halo, momentum_solver,
      {&mx, &my, &mz}, pressure_solver, pressure_preconditioner);
  auto molecular_state = make_state();
  const auto molecular_report =
      molecular_facade.attempt(molecular_state, physics, stencil, {}, {});
  const auto &molecular_base =
      std::get<flow::StepAttemptReport>(molecular_report.base);
  if (molecular_base.disposition !=
          flow::StepAttemptDisposition::committed &&
      mpi.rank() == 0)
    std::cerr << "MOLECULAR_BODY_FAILURE reason="
              << static_cast<int>(molecular_base.reason)
              << " momentum_defect="
              << molecular_base.final_momentum_relative_conservation_defect[0]
              << ','
              << molecular_base.final_momentum_relative_conservation_defect[1]
              << ','
              << molecular_base.final_momentum_relative_conservation_defect[2]
              << '\n';
  HUNDUN_CHECK(molecular_base.disposition ==
               flow::StepAttemptDisposition::committed);

  auto state = make_state();
  const auto before_history = state.snapshot(flow::FlowLayer::history);
  const auto before_committed = state.snapshot(flow::FlowLayer::committed);
  const auto before_metadata = state.metadata();
  flow::test::ConstantDensityPisoTestAccess::reset();
  flow::test::ConstantDensityPisoTestAccess::set_attempt_failure_stage(
      flow::test::AttemptFailureStage::after_momentum);
  const auto failed = facade.attempt(state, physics, stencil, {}, {});
  const auto &failed_base = std::get<flow::StepAttemptReport>(failed.base);
  HUNDUN_CHECK(failed_base.disposition ==
               flow::StepAttemptDisposition::recoverable_failure);
  HUNDUN_CHECK(!failed.wale.has_value());
  HUNDUN_CHECK(test::flow_layer_values_bitwise_equal(
      before_history, state.snapshot(flow::FlowLayer::history)));
  HUNDUN_CHECK(test::flow_layer_values_bitwise_equal(
      before_committed, state.snapshot(flow::FlowLayer::committed)));
  HUNDUN_CHECK(test::accepted_step_metadata_bitwise_equal(before_metadata,
                                                          state.metadata()));

  flow::test::ConstantDensityPisoTestAccess::set_attempt_failure_stage(
      flow::test::AttemptFailureStage::none);
  const auto accepted = facade.attempt(state, physics, stencil, {}, {});
  const auto &accepted_base = std::get<flow::StepAttemptReport>(accepted.base);
  if (accepted_base.disposition !=
          flow::StepAttemptDisposition::committed &&
      mpi.rank() == 0)
    std::cerr << "WALE_BODY_FAILURE reason="
              << static_cast<int>(accepted_base.reason)
              << " rank=" << accepted_base.lowest_failing_rank
              << " correctors=" << accepted_base.pressure_corrector_count
              << " continuity="
              << accepted_base.final_continuity_normalized_l2
              << " pressure=" << accepted_base.final_pressure_residual_l2
              << " momentum="
              << accepted_base.final_momentum_normalized_l2[0] << ','
              << accepted_base.final_momentum_normalized_l2[1] << ','
              << accepted_base.final_momentum_normalized_l2[2]
              << " mass_defect="
              << accepted_base.final_mass_relative_conservation_defect
              << " momentum_defect="
              << accepted_base.final_momentum_relative_conservation_defect[0]
              << ','
              << accepted_base.final_momentum_relative_conservation_defect[1]
              << ','
              << accepted_base.final_momentum_relative_conservation_defect[2]
              << '\n';
  HUNDUN_CHECK(accepted_base.disposition ==
               flow::StepAttemptDisposition::committed);
  HUNDUN_CHECK(accepted_base.pressure_corrector_count == 2U);
  HUNDUN_CHECK(accepted.wale.has_value());
  HUNDUN_CHECK(accepted.wale->identity.value != 0U);
  HUNDUN_CHECK(accepted.wale->owned_active_count ==
               topology.owned_cell_count());
  double maximum = accepted.wale->maximum_nu_t_m2_per_s;
  mpi.allreduce_fp64_in_place(&maximum, 1U,
                              runtime::Fp64ReductionOperation::maximum);
  HUNDUN_CHECK(maximum > 0.0);
  if (formal_scientific) {
    HUNDUN_CHECK(mpi.size() == 1);
    double nu_t_square = accepted.wale->l2_nu_t_m2_per_s *
                         accepted.wale->l2_nu_t_m2_per_s *
                         geometry.cell_volume_m3(0U);
    mpi.allreduce_fp64_in_place(&nu_t_square, 1U,
                                runtime::Fp64ReductionOperation::sum);
    const double nu_t_l2 = std::sqrt(nu_t_square);
    double conservation = accepted_base.final_mass_relative_conservation_defect;
    for (const double value :
         accepted_base.final_momentum_relative_conservation_defect)
      conservation = std::max(conservation, value);
    const test::stage3::ScientificRow row{
        "wale-channel-n48-r1",
        selected_extent,
        mpi.size(),
        process_grid(mpi.size()),
        1U,
        kDt,
        kDt,
        {false, 0.0},
        {false, 0.0},
        {true, nu_t_l2},
        {true, accepted_base.final_continuity_normalized_l2},
        {true, conservation},
        {false, 0.0},
        {false, 0.0},
        {true, accepted.wale->identity.value},
        "pass"};
    HUNDUN_CHECK(test::stage3::validate_scientific_row(row));
    if (mpi.rank() == 0)
      std::cout << test::stage3::serialize_scientific_row(row) << '\n';
  }
  if (!screen_only && mpi.size() == 1) {
    HUNDUN_CHECK(
        field_fingerprint(state.snapshot(flow::FlowLayer::committed)) ==
        UINT64_C(15824403037730263857));
    HUNDUN_CHECK(report_fingerprint(accepted_base, *accepted.wale) ==
                 UINT64_C(18251920989929385807));
  }
  if (screen_only)
    return;

  auto repeat_state = make_state();
  const auto repeated = facade.attempt(repeat_state, physics, stencil, {}, {});
  const auto &repeated_base = std::get<flow::StepAttemptReport>(repeated.base);
  HUNDUN_CHECK(repeated_base.disposition ==
               flow::StepAttemptDisposition::committed);
  HUNDUN_CHECK(repeated.wale.has_value());
  HUNDUN_CHECK(repeated.wale->identity == accepted.wale->identity);

  runtime::FieldRegistry transport_registry;
  flow::FlowFieldIds transport_fields;
  transport_fields.density =
      transport_registry.declare_field(cell_field("transport_rho", 1U));
  transport_fields.velocity =
      transport_registry.declare_field(cell_field("transport_u", 3U));
  transport_fields.mechanical_pressure =
      transport_registry.declare_field(cell_field("transport_pi", 1U));
  transport_fields.face_velocity =
      transport_registry.declare_field(face_field("transport_uf", 3U));
  transport_fields.face_mass_flux =
      finite_volume::declare_face_mass_flux(transport_registry);
  const auto enthalpy =
      transport_registry.declare_field(cell_field("h", 1U));
  const auto scalar =
      transport_registry.declare_field(cell_field("alpha", 1U));
  transport_fields.transported_cell_fields = {enthalpy, scalar};
  transport_registry.freeze();
  const auto make_transport_state = [&] {
    auto result = flow::FlowState::create(
        transport_registry,
        {decomposition.local_extent(), topology.local_face_count()},
        transport_fields,
        {0U, 0.0, kDt, 0.0, flow::MomentumTimeOrder::backward_euler});
    auto initial = initial_values(topology, geometry);
    initial.transported_cell_fields.resize(2U);
    initial.transported_cell_fields[0].resize(topology.owned_cell_count());
    initial.transported_cell_fields[1].resize(topology.owned_cell_count());
    for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count();
         ++cell) {
      const auto point = geometry.cell_center_m(cell);
      initial.transported_cell_fields[0][cell] =
          1.0 + 0.2 * std::sin(point.x) * std::sin(point.y);
      initial.transported_cell_fields[1][cell] =
          0.5 + 0.1 * std::cos(point.x) * std::cos(point.y);
    }
    result.seed_accepted_layers(initial, initial);
    return result;
  };
  auto transport_flow = flow::FixedStepConstantDensityFlow::create(
      decomposition, topology, geometry, boundaries, mpi, execution, halo,
      momentum_solver, {&mx, &my, &mz}, pressure_solver,
      pressure_preconditioner,
      {{enthalpy, finite_volume::FiniteVolumeQuantity::enthalpy(), 0.0},
       {scalar, finite_volume::FiniteVolumeQuantity::scalar(0U), 0.0}});
  auto low_divisor_wale = les::WaleModel::create(
      {0.5, 0.5, 0.5}, topology, geometry, topology.owned_cell_count(), active,
      execution);
  auto high_divisor_wale = les::WaleModel::create(
      {0.5, 2.0, 2.0}, topology, geometry, topology.owned_cell_count(), active,
      execution);
  auto low_state = make_transport_state();
  auto high_state = make_transport_state();
  les::WaleSummary low_summary, high_summary;
  const auto low_report = transport_flow.attempt_with_wale(
      low_state, 1.0, 0.1, stencil, {}, {}, low_divisor_wale, low_summary);
  const auto high_report = transport_flow.attempt_with_wale(
      high_state, 1.0, 0.1, stencil, {}, {}, high_divisor_wale, high_summary);
  HUNDUN_CHECK(low_report.disposition ==
               flow::StepAttemptDisposition::committed);
  HUNDUN_CHECK(high_report.disposition ==
               flow::StepAttemptDisposition::committed);
  const auto low_values = low_state.snapshot(flow::FlowLayer::committed);
  const auto high_values = high_state.snapshot(flow::FlowLayer::committed);
  double maximum_transport_difference = 0.0;
  for (std::size_t field = 0U;
       field < low_values.transported_cell_fields.size(); ++field)
    for (std::size_t cell = 0U;
         cell < low_values.transported_cell_fields[field].size(); ++cell)
      maximum_transport_difference = std::max(
          maximum_transport_difference,
          std::abs(low_values.transported_cell_fields[field][cell] -
                   high_values.transported_cell_fields[field][cell]));
  mpi.allreduce_fp64_in_place(&maximum_transport_difference, 1U,
                              runtime::Fp64ReductionOperation::maximum);
  HUNDUN_CHECK(maximum_transport_difference > 1.0e-14);
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  if (argc == 3 && std::string(argv[1]) == "formal" &&
      std::string(argv[2]) == "channel") {
    selected_extent = {48, 48, 48};
    screen_only = true;
    formal_scientific = true;
  } else if (argc == 2 &&
             (std::string(argv[1]) == "12" || std::string(argv[1]) == "24")) {
    const int cells = std::string(argv[1]) == "12" ? 12 : 24;
    selected_extent = {cells, cells, cells};
    screen_only = true;
  } else if (argc != 1) {
    return 2;
  }
  return hundun::test::run(run);
}
