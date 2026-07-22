// SPDX-License-Identifier: Apache-2.0

#include "flow/src/ideal_gas_closure_test_access.hpp"
#include "hundun/boundary/basic_boundary.hpp"
#include "hundun/config/resolved_case.hpp"
#include "hundun/execution/execution.hpp"
#include "hundun/finite_volume/cell_centered_fvm.hpp"
#include "hundun/flow/ideal_gas_piso.hpp"
#include "hundun/linear/conjugate_gradient.hpp"
#include "hundun/linear/preconditioners.hpp"
#include "hundun/mesh/mesh_geometry.hpp"
#include "hundun/mesh/mesh_topology.hpp"
#include "hundun/runtime/exchange_plan.hpp"
#include "hundun/runtime/halo_exchange.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/mpi_environment.hpp"
#include "hundun/runtime/structured_decomposition.hpp"
#include "tests/support/flow_state_equality.hpp"
#include "tests/support/test_main.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <string_view>
#include <vector>

namespace {

#ifndef HUNDUN_TASK21_DIAGNOSTIC_MATRIX
#define HUNDUN_TASK21_DIAGNOSTIC_MATRIX(source, mpi, state, report)           \
  static_cast<void>(source)
#endif

hundun::runtime::Int3 grid(int ranks) { return {ranks, 1, 1}; }

hundun::config::FlowCaseConfig periodic_case() {
  hundun::config::FlowCaseConfig config{};
  config.schema_version = 2;
  config.simulation_type =
      hundun::config::SimulationType::variable_density_flow;
  config.density_model = hundun::config::DensityModel::ideal_gas;
  config.physics.rho_ref_kg_per_m3 = 1.0;
  config.physics.inlet_consistency_rtol = 1.0e-12;
  config.physics.cp_J_per_kg_K = 1000.0;
  config.physics.gas_constant_J_per_kg_K = 287.05;
  config.physics.thermodynamic_pressure_pa = 101325.0;
  constexpr std::array names{
      hundun::config::PatchName::x_min, hundun::config::PatchName::x_max,
      hundun::config::PatchName::y_min, hundun::config::PatchName::y_max,
      hundun::config::PatchName::z_min, hundun::config::PatchName::z_max};
  for (std::size_t patch = 0; patch < names.size(); ++patch) {
    config.boundaries[patch].patch = names[patch];
    config.boundaries[patch].type = hundun::config::BoundaryType::periodic;
  }
  return config;
}

hundun::config::FlowCaseConfig open_case() {
  auto config = periodic_case();
  config.boundaries[0].type = hundun::config::BoundaryType::velocity_inlet;
  config.boundaries[0].velocity_m_per_s = {1.0, 0.0, 0.0};
  config.boundaries[0].thermal_authority =
      hundun::config::InletThermalAuthority::temperature;
  config.boundaries[0].temperature_K = 300.0;
  config.boundaries[0].enthalpy_J_per_kg = 300000.0;
  config.boundaries[0].density_kg_per_m3 = 101325.0 / (287.05 * 300.0);
  config.boundaries[0].scalar_values =
      std::vector<hundun::config::InletScalarValue>{};
  config.boundaries[1].type = hundun::config::BoundaryType::pressure_outlet;
  config.boundaries[1].pressure_perturbation_pa = 0.0;
  for (std::size_t patch = 2U; patch < config.boundaries.size(); ++patch)
    config.boundaries[patch].type = hundun::config::BoundaryType::symmetry;
  return config;
}

hundun::runtime::FieldDescriptor cell(const char *name, const char *unit,
                                      std::uint32_t components,
                                      bool conservative) {
  return {name,
          unit,
          "task21",
          hundun::runtime::FunctionSpace::cell_average,
          hundun::runtime::ScalarType::float64,
          components,
          2,
          conservative,
          hundun::runtime::RestartPolicy::persistent,
          hundun::runtime::OutputPolicy::never};
}

hundun::runtime::FieldDescriptor face(const char *name, const char *unit,
                                      std::uint32_t components) {
  return {name,
          unit,
          "task21",
          hundun::runtime::FunctionSpace::face_value,
          hundun::runtime::ScalarType::float64,
          components,
          0,
          false,
          hundun::runtime::RestartPolicy::persistent,
          hundun::runtime::OutputPolicy::never};
}

void run(const hundun::runtime::MpiContext &mpi, bool open_domain,
         bool exhaust_generation) {
  constexpr hundun::runtime::Int3 extent{8, 4, 4};
  auto decomposition = hundun::runtime::StructuredDecomposition::create(
      mpi, extent, {!open_domain, !open_domain, !open_domain},
      hundun::runtime::DecompositionOptions{grid(mpi.size())});
  hundun::mesh::MeshTopology topology(decomposition);
  hundun::mesh::MeshGeometry geometry(
      topology,
      hundun::mesh::UniformBoxMapping({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}));
  auto boundaries = hundun::boundary::BoundaryRegistry::create(
      open_domain ? open_case() : periodic_case(), topology);

  hundun::runtime::FieldRegistry registry;
  hundun::flow::FlowFieldIds fields;
  fields.density = registry.declare_field(cell("rho", "kg/m3", 1U, true));
  fields.velocity = registry.declare_field(cell("velocity", "m/s", 3U, false));
  fields.mechanical_pressure =
      registry.declare_field(cell("pi", "Pa", 1U, false));
  fields.face_velocity =
      registry.declare_field(face("face_velocity", "m/s", 3U));
  fields.face_mass_flux =
      hundun::finite_volume::declare_face_mass_flux(registry);
  const auto rho_h = registry.declare_field(cell("rho_h", "J/m3", 1U, true));
  fields.transported_cell_fields = {rho_h};
  registry.freeze();

  constexpr double dt = 1.0e-3;
  constexpr double cp = 1000.0;
  constexpr double gas_constant = 287.05;
  constexpr double initial_temperature = 300.0;
  constexpr double pressure = 101325.0;
  const double density = pressure / (gas_constant * initial_temperature);
  auto state = hundun::flow::FlowState::create(
      registry, {decomposition.local_extent(), topology.local_face_count()},
      fields,
      {0U, 0.0, dt, 0.0, hundun::flow::MomentumTimeOrder::backward_euler});
  hundun::flow::FlowLayerValues initial;
  initial.density.assign(topology.owned_cell_count(), density);
  initial.velocity.assign(topology.owned_cell_count() * 3U, 0.0);
  if (open_domain)
    for (std::size_t cell_index = 0; cell_index < topology.owned_cell_count();
         ++cell_index)
      initial.velocity[cell_index * 3U] = 1.0;
  initial.mechanical_pressure.assign(topology.owned_cell_count(), 0.0);
  initial.face_velocity.assign(topology.local_face_count() * 3U, 0.0);
  initial.face_mass_flux.assign(topology.local_face_count(), 0.0);
  if (open_domain) {
    for (hundun::mesh::LocalFaceId face = 0;
         face < topology.local_face_count(); ++face) {
      const auto area = geometry.face_area_vector_m2(
          face, hundun::mesh::FaceSide::owner);
      initial.face_velocity[face * 3U] = 1.0;
      initial.face_mass_flux[face] = density * area.x;
    }
  }
  initial.transported_cell_fields = {std::vector<double>(
      topology.owned_cell_count(), density * cp * initial_temperature)};
  state.seed_accepted_layers(initial, initial);

  auto closure = hundun::flow::IdealGasClosure::create(
      topology, geometry, boundaries, mpi, registry, fields, state,
      {rho_h, cp, gas_constant, pressure});
  auto halo = hundun::runtime::HaloExchange::create(
      decomposition, hundun::runtime::ExchangePlan::create(
                         decomposition, decomposition.local_extent(), 2));
  hundun::execution::CpuReferenceContext execution;
  hundun::linear::ConjugateGradientSolver momentum_solver(execution, mpi);
  hundun::linear::ConjugateGradientSolver pressure_solver(execution, mpi);
  hundun::linear::JacobiPreconditioner mx(execution), my(execution),
      mz(execution), pressure_preconditioner(execution);
  hundun::flow::MaterialDensityTransportSpec material_spec;
  material_spec.enthalpy_density = rho_h;
  material_spec.enthalpy_diffusivity_kg_per_m_s = 0.0;
  auto flow = hundun::flow::FixedStepIdealGasFlow::create(
      decomposition, topology, geometry, boundaries, mpi, execution, halo,
      momentum_solver, {&mx, &my, &mz}, pressure_solver,
      pressure_preconditioner, registry, fields, material_spec,
      std::move(closure));
  hundun::flow::test::IdealGasClosureTestAccess::set_uniform_enthalpy_rate(
      flow, open_domain ? 0.0 : cp * 5.0 / dt);

  HUNDUN_CHECK(
      hundun::test::ideal_gas_state_equality_oracle_is_mutation_sensitive());
  const auto stencil = hundun::flow::make_momentum_time_stencil(
      hundun::flow::MomentumTimeOrder::backward_euler, dt, 0.0);
  const auto report = flow.attempt(state, 0.0, stencil, {}, {});
  HUNDUN_CHECK(report.flow().flow().disposition ==
               hundun::flow::StepAttemptDisposition::committed);
  HUNDUN_CHECK(report.flow().flow().reason ==
               hundun::flow::StepFailureReason::none);
  HUNDUN_CHECK(report.flow().flow().pressure_corrector_count == 2U);
  HUNDUN_CHECK(report.closure_report_available());
  const auto &closure_report = report.closure_report();
  HUNDUN_CHECK(closure_report.disposition() ==
               hundun::flow::IdealGasClosureDisposition::closed);
  HUNDUN_CHECK(closure_report.stage() ==
               hundun::flow::IdealGasClosureStage::final);
  HUNDUN_CHECK(closure_report.evaluation_count() == 3U);
  HUNDUN_CHECK(closure_report.collective_count() == 14U);
  HUNDUN_CHECK(closure_report.final_metrics_available());
  HUNDUN_CHECK(closure_report.rho_remap_normalized_l2() <= 1.0e-10);
  HUNDUN_CHECK(closure_report.rho_h_remap_normalized_l2() <= 1.0e-9);
  HUNDUN_CHECK(closure_report.eos_max_relative_error() <= 1.0e-12);
  HUNDUN_CHECK(
      hundun::flow::test::IdealGasClosureTestAccess::report_authenticated(
          report));
  HUNDUN_CHECK(hundun::flow::test::IdealGasClosureTestAccess::
                   post_eos_evidence_authenticated(report.flow()));
  HUNDUN_CHECK(state.metadata().step == 1U);
  auto closure_source = flow.closure_diagnostic_source(state, report);
  HUNDUN_CHECK(closure_source.fingerprint_field_count() == 3U);
  HUNDUN_CHECK(closure_source.fingerprint_field_id(0U) ==
               std::string_view("p0"));
  HUNDUN_CHECK(closure_source.fingerprint_field_id(1U) ==
               std::string_view("rho"));
  HUNDUN_CHECK(closure_source.fingerprint_field_id(2U) ==
               std::string_view("rho_h"));
  HUNDUN_CHECK(closure_source.fingerprint_field_item_count(0U) ==
               (mpi.rank() == 0 ? 1U : 0U));
  HUNDUN_CHECK(closure_source.sample_field_count() == 5U);
  HUNDUN_CHECK(closure_source.sample_field_id(4U) ==
               std::string_view("temperature"));
  const double first_temperature = open_domain ? 300.0 : 305.0;
  HUNDUN_CHECK_NEAR(closure_source.sample_field_value(4U, 0U),
                    first_temperature, 5.0e-12 * first_temperature);
  HUNDUN_TASK21_DIAGNOSTIC_MATRIX(closure_source, mpi, state, report);
  auto material_source = flow.flow_diagnostic_source(state, report);
  HUNDUN_CHECK(material_source.report().attempt_identity() ==
               report.attempt_identity());
  const auto closure_state = flow.closure_state();
  HUNDUN_CHECK(closure_state.revision == 1U);
  HUNDUN_CHECK_NEAR(closure_state.thermodynamic_pressure_pa,
                    open_domain ? pressure : pressure * 305.0 / 300.0,
                    1.0e-12 * pressure);
  HUNDUN_CHECK(closure_state.target_mass_kg.has_value() != open_domain);
  const auto committed = state.snapshot(hundun::flow::FlowLayer::committed);
  HUNDUN_CHECK(std::all_of(
      committed.density.begin(), committed.density.end(), [&](double value) {
        return std::abs(value - density) <= 5.0e-12 * density;
      }));
  for (std::size_t cell_index = 0; cell_index < committed.density.size();
       ++cell_index) {
    const double h = committed.transported_cell_fields.front()[cell_index] /
                     committed.density[cell_index];
    HUNDUN_CHECK_NEAR(h / cp, first_temperature, 5.0e-12 * first_temperature);
  }

  const auto bdf2 = hundun::flow::make_momentum_time_stencil(
      hundun::flow::MomentumTimeOrder::bdf2, dt, dt);
  const auto second = flow.attempt(state, 0.0, bdf2, {}, {});
  HUNDUN_CHECK(second.flow().flow().disposition ==
               hundun::flow::StepAttemptDisposition::committed);
  HUNDUN_CHECK(second.closure_report().evaluation_count() == 3U);
  HUNDUN_CHECK(second.closure_report().collective_count() == 14U);
  HUNDUN_CHECK(flow.closure_state().revision == 2U);
  HUNDUN_CHECK_NEAR(flow.closure_state().thermodynamic_pressure_pa,
                    open_domain ? pressure : pressure * 310.0 / 300.0,
                    1.0e-12 * pressure);
  bool stale_rejected = false;
  try {
    static_cast<void>(closure_source.committed_step());
  } catch (const hundun::runtime::Error &error) {
    stale_rejected = std::string_view(error.what()) ==
                     "ideal-gas closure diagnostic source is stale";
  }
  HUNDUN_CHECK(stale_rejected);
  if (exhaust_generation) {
    const hundun::test::IdealGasStateSnapshot before_failure{
        state.snapshot(hundun::flow::FlowLayer::history),
        state.snapshot(hundun::flow::FlowLayer::committed),
        state.snapshot(hundun::flow::FlowLayer::trial), state.metadata(),
        flow.closure_state()};
    hundun::flow::test::IdealGasClosureTestAccess::set_uniform_enthalpy_rate(
        flow, -cp * 1000.0 / dt);
    const auto failed = flow.attempt(state, 0.0, bdf2, {}, {});
    HUNDUN_CHECK(failed.flow().flow().disposition ==
                 hundun::flow::StepAttemptDisposition::recoverable_failure);
    HUNDUN_CHECK(failed.flow().flow().reason ==
                 hundun::flow::StepFailureReason::density_closure_failure);
    HUNDUN_CHECK(failed.closure_report_available());
    HUNDUN_CHECK(failed.closure_report().stage() ==
                 hundun::flow::IdealGasClosureStage::predictor);
    const hundun::test::IdealGasStateSnapshot after_failure{
        state.snapshot(hundun::flow::FlowLayer::history),
        state.snapshot(hundun::flow::FlowLayer::committed),
        state.snapshot(hundun::flow::FlowLayer::trial), state.metadata(),
        flow.closure_state()};
    HUNDUN_CHECK(hundun::test::ideal_gas_state_bitwise_equal(before_failure,
                                                             after_failure));
    hundun::flow::test::IdealGasClosureTestAccess::exhaust_source_generation(
        flow);
    bool exhausted_rejected = false;
    try {
      static_cast<void>(flow.attempt(state, 0.0, bdf2, {}, {}));
    } catch (const hundun::runtime::Error &) {
      exhausted_rejected = true;
    }
    HUNDUN_CHECK(exhausted_rejected);
  }
}

} // namespace

#ifndef HUNDUN_TASK21_IDEAL_GAS_FIXTURE_ONLY
int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  return hundun::test::run([&] {
    auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
    HUNDUN_CHECK(mpi.size() == 1 || mpi.size() == 2 || mpi.size() == 4);
    const bool open_domain =
        argc > 1 && std::string_view(argv[1]) == "--open-plug";
    run(mpi, open_domain, argc == 1);
  });
}
#endif
