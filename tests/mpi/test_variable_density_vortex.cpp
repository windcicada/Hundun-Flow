// SPDX-License-Identifier: Apache-2.0

#include "flow/src/material_density_piso_test_access.hpp"
#include "hundun/boundary/basic_boundary.hpp"
#include "hundun/config/resolved_case.hpp"
#include "hundun/execution/execution.hpp"
#include "hundun/finite_volume/cell_centered_fvm.hpp"
#include "hundun/flow/material_density_piso.hpp"
#include "hundun/linear/conjugate_gradient.hpp"
#include "hundun/linear/preconditioners.hpp"
#include "hundun/mesh/mesh_geometry.hpp"
#include "hundun/mesh/mesh_topology.hpp"
#include "hundun/runtime/exchange_plan.hpp"
#include "hundun/runtime/field_access_plan.hpp"
#include "hundun/runtime/field_storage.hpp"
#include "hundun/runtime/halo_exchange.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/mpi_environment.hpp"
#include "hundun/runtime/structured_decomposition.hpp"
#include "tests/support/test_main.hpp"

#include <array>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <string_view>
#include <vector>

namespace {
constexpr double pi = 3.141592653589793238462643383279502884;

hundun::runtime::Int3 grid(int ranks) {
  return ranks == 1 ? hundun::runtime::Int3{1, 1, 1}
                    : ranks == 2 ? hundun::runtime::Int3{2, 1, 1}
                                 : hundun::runtime::Int3{4, 1, 1};
}

hundun::config::FlowCaseConfig periodic_case() {
  hundun::config::FlowCaseConfig config{};
  config.schema_version = 2;
  config.simulation_type =
      hundun::config::SimulationType::variable_density_flow;
  config.density_model = hundun::config::DensityModel::constant;
  config.physics.rho_ref_kg_per_m3 = 1.0;
  config.physics.inlet_consistency_rtol = 1.0e-12;
  constexpr std::array<hundun::config::PatchName, 6> names{
      hundun::config::PatchName::x_min, hundun::config::PatchName::x_max,
      hundun::config::PatchName::y_min, hundun::config::PatchName::y_max,
      hundun::config::PatchName::z_min, hundun::config::PatchName::z_max};
  for (std::size_t i = 0; i < names.size(); ++i) {
    config.boundaries[i].patch = names[i];
    config.boundaries[i].type = hundun::config::BoundaryType::periodic;
  }
  return config;
}

hundun::runtime::FieldDescriptor cell(const char *name, const char *unit,
                                      std::uint32_t components,
                                      bool conservative) {
  return {name, unit, "task20-vortex",
          hundun::runtime::FunctionSpace::cell_average,
          hundun::runtime::ScalarType::float64, components, 2, conservative,
          hundun::runtime::RestartPolicy::persistent,
          hundun::runtime::OutputPolicy::never};
}

hundun::runtime::FieldDescriptor face(const char *name,
                                      std::uint32_t components) {
  return {name, "m/s", "task20-vortex",
          hundun::runtime::FunctionSpace::face_value,
          hundun::runtime::ScalarType::float64, components, 0, false,
          hundun::runtime::RestartPolicy::persistent,
          hundun::runtime::OutputPolicy::never};
}

hundun::runtime::Real3 velocity(double x, double y) {
  return {std::sin(x) * std::cos(y), -std::cos(x) * std::sin(y), 0.0};
}

double density(double x, double y) {
  return 1.0 + 0.1 * std::sin(x) * std::sin(y);
}

double density_cell_average(double x, double y, double spacing) {
  const double factor = std::sin(0.5 * spacing) / (0.5 * spacing);
  return 1.0 + 0.1 * factor * factor * std::sin(x) * std::sin(y);
}

std::uint64_t bits(double value) {
  std::uint64_t result{};
  std::memcpy(&result, &value, sizeof(result));
  return result;
}
} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  return hundun::test::run([&] {
    auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
    HUNDUN_CHECK(mpi.size() == 1 || mpi.size() == 2 || mpi.size() == 4);

    const auto direct =
        hundun::flow::test::MaterialDensityPisoTestAccess::vortex_source(
            0.25 * pi, 0.5 * pi, 0.01);
    HUNDUN_CHECK_NEAR(
        direct.x,
        density(0.25 * pi, 0.5 * pi) * std::sin(0.25 * pi) *
                std::cos(0.25 * pi) +
            0.02 * std::sin(0.25 * pi) * std::cos(0.5 * pi),
        64.0 * std::numeric_limits<double>::epsilon());
    HUNDUN_CHECK(direct.z == 0.0);

    {
      constexpr double predictor_velocity_interpolation = 0.4;
      constexpr double predictor_momentum_interpolation = 0.72;
      constexpr double face_density = 1.2;
      constexpr double cell_momentum_n = 0.66;
      constexpr double face_momentum_n = 0.9;
      constexpr double cell_momentum_old = -0.15;
      constexpr double face_momentum_old = 0.05;
      constexpr double mobility = 0.08;
      constexpr double dt_s = 0.2;
      constexpr double alpha1 = -2.0;
      constexpr double alpha2 = 0.5;
      constexpr double pressure_correction = -0.03;
      const double value =
          hundun::flow::test::MaterialDensityPisoTestAccess::
              material_face_value(
                  predictor_velocity_interpolation,
                  predictor_momentum_interpolation, face_density,
                  cell_momentum_n, face_momentum_n, cell_momentum_old,
                  face_momentum_old, mobility, dt_s, alpha1, alpha2,
                  pressure_correction);
      const double delta_n = face_momentum_n - cell_momentum_n;
      const double delta_old = face_momentum_old - cell_momentum_old;
      const double time_term =
          (mobility / dt_s) *
          ((-alpha1) * delta_n + (-alpha2) * delta_old);
      const double expected = predictor_momentum_interpolation / face_density +
                              pressure_correction + time_term;
      HUNDUN_CHECK_NEAR(value, expected,
                        32.0 * std::numeric_limits<double>::epsilon());
      const double omitted = predictor_momentum_interpolation / face_density +
                             pressure_correction;
      const double plain_velocity_base = predictor_velocity_interpolation +
                                         pressure_correction + time_term;
      const double wrong_delta_n =
          face_density *
          (face_momentum_n / face_density -
           predictor_velocity_interpolation);
      const double wrong_delta_old =
          face_density *
          (face_momentum_old / face_density - (-0.1));
      const double wrong_discrepancy =
          predictor_momentum_interpolation / face_density +
          pressure_correction +
          (mobility / dt_s) *
              ((-alpha1) * wrong_delta_n + (-alpha2) * wrong_delta_old);
      HUNDUN_CHECK(std::abs(value - omitted) > 1.0e-3);
      HUNDUN_CHECK(std::abs(value - plain_velocity_base) > 1.0e-3);
      HUNDUN_CHECK(std::abs(value - wrong_discrepancy) > 1.0e-3);
    }

    struct Errors final {
      double rho{};
      double rho_u{};
      double rho_v{};
      double rho_w{};
      double pressure{};
      double pressure_mean{};
      double parity_pressure{};
      double parity_pressure_mean{};
    };
    const auto run_case = [&](int n) {
    const hundun::runtime::Int3 extent{n, n, 4};
    auto decomposition = hundun::runtime::StructuredDecomposition::create(
        mpi, extent, {true, true, true},
        hundun::runtime::DecompositionOptions{grid(mpi.size())});
    hundun::mesh::MeshTopology topology(decomposition);
    hundun::mesh::MeshGeometry geometry(
        topology, hundun::mesh::UniformBoxMapping(
                      {0.0, 0.0, 0.0}, {2.0 * pi, 2.0 * pi, 1.0}));
    auto boundaries = hundun::boundary::BoundaryRegistry::create(
        periodic_case(), topology);
    hundun::runtime::FieldRegistry registry;
    hundun::flow::FlowFieldIds fields;
    fields.density = registry.declare_field(cell("rho", "kg/m3", 1U, true));
    fields.velocity =
        registry.declare_field(cell("velocity", "m/s", 3U, false));
    fields.mechanical_pressure =
        registry.declare_field(cell("pi", "Pa", 1U, false));
    fields.face_velocity = registry.declare_field(face("face_velocity", 3U));
    fields.face_mass_flux =
        hundun::finite_volume::declare_face_mass_flux(registry);
    const auto rho_h =
        registry.declare_field(cell("rho_h", "J/m3", 1U, true));
    const auto initial_face_density =
        registry.declare_field(face("initial_face_density", 1U));
    fields.transported_cell_fields = {rho_h};
    registry.freeze();

    const auto box = decomposition.owned_box();
    const hundun::runtime::Int3 local{box.end.x - box.begin.x,
                                      box.end.y - box.begin.y,
                                      box.end.z - box.begin.z};
    auto halo = hundun::runtime::HaloExchange::create(
        decomposition,
        hundun::runtime::ExchangePlan::create(decomposition, local, 2));
    auto state = hundun::flow::FlowState::create(
        registry, {local, topology.local_face_count()}, fields,
        {0U, 0.0, 1.0e-4, 1.0e-4,
         hundun::flow::MomentumTimeOrder::bdf2});
    hundun::flow::FlowLayerValues exact;
    exact.density.resize(topology.owned_cell_count());
    exact.velocity.resize(topology.owned_cell_count() * 3U);
    exact.mechanical_pressure.assign(topology.owned_cell_count(), 0.0);
    exact.transported_cell_fields.resize(1U);
    exact.transported_cell_fields[0].resize(topology.owned_cell_count());
    for (hundun::mesh::LocalCellId cell_id = 0;
         cell_id < topology.owned_cell_count(); ++cell_id) {
      const auto centre = geometry.cell_center_m(cell_id);
      const auto u = velocity(centre.x, centre.y);
      const double rho =
          n == 8 ? density(centre.x, centre.y)
                 : density_cell_average(
                       centre.x, centre.y,
                       2.0 * pi / static_cast<double>(n));
      exact.density[cell_id] = rho;
      exact.velocity[cell_id * 3U] = u.x;
      exact.velocity[cell_id * 3U + 1U] = u.y;
      exact.velocity[cell_id * 3U + 2U] = 0.0;
      exact.transported_cell_fields[0][cell_id] = rho;
    }
    exact.face_velocity.resize(topology.local_face_count() * 3U);
    exact.face_mass_flux.resize(topology.local_face_count());
    for (hundun::mesh::LocalFaceId face_id = 0;
         face_id < topology.local_face_count(); ++face_id) {
      const auto centre = geometry.face_center_m(face_id);
      const auto u = velocity(centre.x, centre.y);
      exact.face_velocity[face_id * 3U] = u.x;
      exact.face_velocity[face_id * 3U + 1U] = u.y;
      exact.face_velocity[face_id * 3U + 2U] = 0.0;
      const auto area =
          geometry.face_area_vector_m2(face_id, hundun::mesh::FaceSide::owner);
      exact.face_mass_flux[face_id] =
          u.x * area.x + u.y * area.y + u.z * area.z;
    }
    constexpr hundun::runtime::PhaseId initial_phase = 2020U;
    constexpr hundun::runtime::ActorId initial_actor = 2020U;
    hundun::runtime::FieldAccessPlan initial_access(registry);
    initial_access.declare_access(initial_phase, initial_actor, fields.density,
                                  hundun::runtime::AccessMode::read_write);
    initial_access.declare_access(initial_phase, initial_actor,
                                  fields.face_mass_flux,
                                  hundun::runtime::AccessMode::read_write);
    initial_access.declare_access(initial_phase, initial_actor,
                                  initial_face_density,
                                  hundun::runtime::AccessMode::read_write);
    initial_access.freeze();
    hundun::runtime::FieldStorage initial_storage(
        registry, {local, topology.local_face_count()});
    {
      auto rho = initial_storage.acquire_write<double>(
          initial_access, initial_phase, initial_actor, fields.density);
      for (hundun::mesh::LocalCellId cell_id = 0;
           cell_id < topology.owned_cell_count(); ++cell_id) {
        const auto global = topology.global_cell(cell_id);
        rho(global.x - box.begin.x, global.y - box.begin.y,
            global.z - box.begin.z, 0) = exact.density[cell_id];
      }
      auto direction = initial_storage.acquire_face_write<double>(
          initial_access, initial_phase, initial_actor, fields.face_mass_flux);
      for (hundun::mesh::LocalFaceId face_id = 0;
           face_id < topology.local_face_count(); ++face_id)
        direction(face_id, 0) = exact.face_mass_flux[face_id];
    }
    halo.exchange(initial_storage, fields.density);
    {
      const auto direction = hundun::finite_volume::FaceMassFlux::acquire(
          registry, initial_storage, initial_access, initial_phase,
          initial_actor, fields.face_mass_flux, topology);
      const auto rho = initial_storage.acquire_read<double>(
          initial_access, initial_phase, initial_actor, fields.density);
      auto rho_face = initial_storage.acquire_face_write<double>(
          initial_access, initial_phase, initial_actor, initial_face_density);
      auto fvm = hundun::finite_volume::CellCenteredFvmOperators::create(
          topology, geometry);
      fvm.reconstruct_transport_faces(
          hundun::finite_volume::FiniteVolumeQuantity::density(), boundaries,
          direction, rho, rho_face);
      for (hundun::mesh::LocalFaceId face_id = 0;
           face_id < topology.local_face_count(); ++face_id) {
        exact.face_mass_flux[face_id] *= rho_face(face_id, 0);
      }
    }
    state.seed_accepted_layers(exact, exact);

    hundun::execution::CpuReferenceContext execution;
    hundun::linear::ConjugateGradientSolver momentum_solver(execution, mpi);
    hundun::linear::ConjugateGradientSolver pressure_solver(execution, mpi);
    hundun::linear::JacobiPreconditioner mx(execution), my(execution),
        mz(execution), pressure_preconditioner(execution);
    hundun::flow::MaterialDensityTransportSpec specification;
    specification.enthalpy_density = rho_h;
    specification.enthalpy_diffusivity_kg_per_m_s = 0.0;
    auto flow = hundun::flow::FixedStepMaterialDensityFlow::create(
        decomposition, topology, geometry, boundaries, mpi, execution, halo,
        momentum_solver, {&mx, &my, &mz}, pressure_solver,
        pressure_preconditioner, registry, fields, specification);
    hundun::flow::test::MaterialDensityPisoTestAccess::enable_vortex_source(
        flow, true);
    const auto stencil = hundun::flow::make_momentum_time_stencil(
        hundun::flow::MomentumTimeOrder::bdf2, 1.0e-4, 1.0e-4);
    const auto report = flow.attempt(state, 0.01, stencil, {}, {});
    if (mpi.rank() == 0 &&
        report.flow().disposition !=
            hundun::flow::StepAttemptDisposition::committed) {
      std::cerr << "Task20 vortex failure disposition="
                << static_cast<int>(report.flow().disposition)
                << " reason=" << static_cast<int>(report.flow().reason)
                << " correctors=" << report.flow().pressure_corrector_count
                << " pressure-solve="
                << static_cast<int>(report.flow().pressure[0].reason)
                << " pressure-final="
                << report.flow().pressure[0].final_residual
                << " pressure-initial="
                << report.flow().pressure[0].initial_residual
                << " pressure-iterations="
                << report.flow().pressure[0].iterations
                << " material="
                << static_cast<int>(report.material_failure_reason())
                << " continuity="
                << report.flow().final_continuity_normalized_l2
                << " pressure=" << report.final_pressure_normalized_residual()
                << " momentum=" << report.flow().final_momentum_normalized_l2[0]
                << ',' << report.flow().final_momentum_normalized_l2[1] << ','
                << report.flow().final_momentum_normalized_l2[2]
                << " mass-defect="
                << report.flow().final_mass_relative_conservation_defect
                << " material-mass="
                << (report.material_report_available()
                        ? report.material_report()
                              .mass_relative_conservation_defect()
                        : -1.0)
                << " momentum-defect="
                << report.flow().final_momentum_relative_conservation_defect[0]
                << ','
                << report.flow().final_momentum_relative_conservation_defect[1]
                << ','
                << report.flow().final_momentum_relative_conservation_defect[2]
                << '\n';
    }
    HUNDUN_CHECK(report.flow().disposition ==
                 hundun::flow::StepAttemptDisposition::committed);
    HUNDUN_CHECK(report.flow().pressure_corrector_count == 2U);
    HUNDUN_CHECK(report.final_continuity_residual_available());
    HUNDUN_CHECK(report.flow().final_continuity_normalized_l2 <= 1.0e-10);
    HUNDUN_CHECK(report.final_pressure_residual_available());
    HUNDUN_CHECK(report.final_pressure_normalized_residual() <= 1.0);
    for (double residual : report.flow().final_momentum_normalized_l2)
      HUNDUN_CHECK(residual <= 1.0e-9);
    HUNDUN_CHECK(report.material_report().disposition() ==
                 hundun::flow::MaterialTransportDisposition::finalized);
    HUNDUN_CHECK(report.flux_provenance() ==
                 hundun::flow::MaterialFluxProvenance::final_corrected);
    HUNDUN_CHECK(report.shared_face_mass_flux_field() ==
                 fields.face_mass_flux);
    HUNDUN_CHECK(report.material_report().density_normalized_l2() <= 1.0e-10);
    for (const auto residual :
         report.material_report().transport_normalized_l2())
      HUNDUN_CHECK(residual <= 1.0e-9);
    HUNDUN_CHECK(report.flow().final_mass_relative_conservation_defect <=
                 5.0e-12);
    HUNDUN_CHECK(
        report.material_report().mass_relative_conservation_defect() <=
        5.0e-12);
    for (const auto defect :
         report.flow().final_momentum_relative_conservation_defect)
      HUNDUN_CHECK(defect <= 5.0e-11);
    HUNDUN_CHECK(report.material_report().minimum_density_kg_per_m3() > 0.0);
    const auto committed = state.snapshot(hundun::flow::FlowLayer::committed);
    const auto &finalizer_input =
        hundun::flow::test::MaterialDensityPisoTestAccess::
            finalizer_flux_evidence(flow);
    HUNDUN_CHECK(finalizer_input.size() == committed.face_mass_flux.size());
    for (std::size_t face_index = 0; face_index < finalizer_input.size();
         ++face_index)
      HUNDUN_CHECK(bits(finalizer_input[face_index]) ==
                   bits(committed.face_mass_flux[face_index]));

    const std::size_t global_cells = static_cast<std::size_t>(n) *
                                     static_cast<std::size_t>(n) * 4U;
    std::vector<double> cell_cover(global_cells, 0.0);
    std::vector<double> global_fields(global_cells * 5U, 0.0);
    for (hundun::mesh::LocalCellId cell_id = 0;
         cell_id < topology.owned_cell_count(); ++cell_id) {
      const auto global =
          static_cast<std::size_t>(topology.global_cell_id(cell_id));
      HUNDUN_CHECK(global < global_cells);
      cell_cover[global] = 1.0;
      global_fields[global * 5U] = committed.density[cell_id];
      global_fields[global * 5U + 1U] =
          committed.density[cell_id] * committed.velocity[cell_id * 3U];
      global_fields[global * 5U + 2U] =
          committed.density[cell_id] * committed.velocity[cell_id * 3U + 1U];
      global_fields[global * 5U + 3U] =
          committed.density[cell_id] * committed.velocity[cell_id * 3U + 2U];
      global_fields[global * 5U + 4U] =
          committed.mechanical_pressure[cell_id];
    }
    mpi.allreduce_fp64_in_place(cell_cover.data(), cell_cover.size(),
                                hundun::runtime::Fp64ReductionOperation::sum);
    mpi.allreduce_fp64_in_place(global_fields.data(), global_fields.size(),
                                hundun::runtime::Fp64ReductionOperation::sum);
    HUNDUN_CHECK(std::all_of(cell_cover.begin(), cell_cover.end(),
                            [](double value) { return value == 1.0; }));
    double local_max_difference{};
    double local_infinity_norm{};
    for (hundun::mesh::LocalCellId cell_id = 0;
         cell_id < topology.owned_cell_count(); ++cell_id) {
      const auto global =
          static_cast<std::size_t>(topology.global_cell_id(cell_id));
      const std::array<double, 5> values{
          committed.density[cell_id],
          committed.density[cell_id] * committed.velocity[cell_id * 3U],
          committed.density[cell_id] * committed.velocity[cell_id * 3U + 1U],
          committed.density[cell_id] * committed.velocity[cell_id * 3U + 2U],
          committed.mechanical_pressure[cell_id]};
      for (std::size_t component = 0; component < values.size(); ++component) {
        local_max_difference = std::max(
            local_max_difference,
            std::abs(values[component] -
                     global_fields[global * values.size() + component]));
        local_infinity_norm =
            std::max(local_infinity_norm, std::abs(values[component]));
      }
    }
    double decomposition_values[2]{local_max_difference,
                                   local_infinity_norm};
    mpi.allreduce_fp64_in_place(
        decomposition_values, 2U,
        hundun::runtime::Fp64ReductionOperation::maximum);
    HUNDUN_CHECK(decomposition_values[0] <=
                 5.0e-12 * std::max(1.0, decomposition_values[1]));

    std::vector<double> face_cover(topology.global_face_count(), 0.0);
    std::vector<double> global_face_flux(topology.global_face_count(), 0.0);
    for (hundun::mesh::LocalFaceId face_id = 0;
         face_id < topology.local_face_count(); ++face_id) {
      if (topology.cell_ownership(topology.owner(face_id)) !=
          hundun::mesh::EntityOwnership::owned)
        continue;
      const auto global =
          static_cast<std::size_t>(topology.global_face_id(face_id));
      HUNDUN_CHECK(global < face_cover.size());
      face_cover[global] = 1.0;
      global_face_flux[global] = committed.face_mass_flux[face_id];
    }
    mpi.allreduce_fp64_in_place(face_cover.data(), face_cover.size(),
                                hundun::runtime::Fp64ReductionOperation::sum);
    mpi.allreduce_fp64_in_place(
        global_face_flux.data(), global_face_flux.size(),
        hundun::runtime::Fp64ReductionOperation::sum);
    HUNDUN_CHECK(std::all_of(face_cover.begin(), face_cover.end(),
                            [](double value) { return value == 1.0; }));
    std::array<double, 12> sums{};
    for (hundun::mesh::LocalCellId cell_id = 0;
         cell_id < topology.owned_cell_count(); ++cell_id) {
      const auto centre = geometry.cell_center_m(cell_id);
      const auto exact_u = velocity(centre.x, centre.y);
      const double exact_rho = density(centre.x, centre.y);
      const double volume = geometry.cell_volume_m3(cell_id);
      const double numeric_rho = committed.density[cell_id];
      const double numeric_u = committed.velocity[cell_id * 3U];
      const double numeric_v = committed.velocity[cell_id * 3U + 1U];
      const double numeric_w = committed.velocity[cell_id * 3U + 2U];
      const double rho_error = numeric_rho - exact_rho;
      const double rho_u_error = numeric_rho * numeric_u - exact_rho * exact_u.x;
      const double rho_v_error = numeric_rho * numeric_v - exact_rho * exact_u.y;
      const double rho_w_error = numeric_rho * numeric_w;
      const double pressure_value = committed.mechanical_pressure[cell_id];
      const auto global = topology.global_cell(cell_id);
      const double parity = ((global.x + global.y + global.z) & 1) == 0
                                ? 1.0
                                : -1.0;
      const double mutated_pressure = pressure_value + parity;
      sums[0] += volume * rho_error * rho_error;
      sums[1] += volume * exact_rho * exact_rho;
      sums[2] += volume * rho_u_error * rho_u_error;
      sums[3] += volume * (exact_rho * exact_u.x) * (exact_rho * exact_u.x);
      sums[4] += volume * rho_v_error * rho_v_error;
      sums[5] += volume * (exact_rho * exact_u.y) * (exact_rho * exact_u.y);
      sums[6] += volume * rho_w_error * rho_w_error;
      sums[7] += volume * pressure_value * pressure_value;
      sums[8] += volume;
      sums[9] += volume * pressure_value;
      sums[10] += volume * mutated_pressure * mutated_pressure;
      sums[11] += volume * mutated_pressure;
    }
    mpi.allreduce_fp64_in_place(sums.data(), sums.size(),
                                hundun::runtime::Fp64ReductionOperation::sum);
    return Errors{std::sqrt(sums[0] / sums[1]),
                  std::sqrt(sums[2] / sums[3]),
                  std::sqrt(sums[4] / sums[5]),
                  std::sqrt(sums[6] / sums[8]),
                  std::sqrt(sums[7] / sums[8]),
                  sums[9] / sums[8], std::sqrt(sums[10] / sums[8]),
                  sums[11] / sums[8]};
    };

    const bool full = argc > 1 && std::string_view(argv[1]) == "--full";
    if (!full) {
      static_cast<void>(run_case(8));
      return;
    }
    const std::array<int, 3> resolutions{16, 32, 64};
    std::array<Errors, 3> errors{};
    for (std::size_t index = 0; index < resolutions.size(); ++index)
      errors[index] = run_case(resolutions[index]);
    const auto require_order = [&](auto select) {
      for (std::size_t index = 0; index + 1U < errors.size(); ++index) {
        const double coarse = select(errors[index]);
        const double fine = select(errors[index + 1U]);
        HUNDUN_CHECK(std::isfinite(coarse) && coarse > 0.0);
        HUNDUN_CHECK(std::isfinite(fine) && fine > 0.0 && fine < coarse);
        HUNDUN_CHECK(std::log(coarse / fine) / std::log(2.0) >= 1.8);
      }
    };
    require_order([](const Errors &value) { return value.rho; });
    require_order([](const Errors &value) { return value.rho_u; });
    require_order([](const Errors &value) { return value.rho_v; });
    require_order([](const Errors &value) { return value.pressure; });
    for (std::size_t index = 0; index < errors.size(); ++index) {
      const double h = 2.0 * pi / static_cast<double>(resolutions[index]);
      HUNDUN_CHECK(errors[index].rho_w <= 1.0e-12);
      HUNDUN_CHECK(errors[index].pressure <= 2.0 * 1.1 * h * h);
      HUNDUN_CHECK(std::abs(errors[index].pressure_mean) <= 1.0e-12);
      HUNDUN_CHECK(errors[index].parity_pressure > 2.0 * 1.1 * h * h);
      HUNDUN_CHECK(std::abs(errors[index].parity_pressure_mean) <= 1.0e-12);
    }
  });
}
