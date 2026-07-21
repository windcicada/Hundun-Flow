// SPDX-License-Identifier: Apache-2.0

#include "hundun/boundary/basic_boundary.hpp"
#include "hundun/config/resolved_case.hpp"
#include "hundun/finite_volume/cell_centered_fvm.hpp"
#include "hundun/flow/flow_state.hpp"
#include "hundun/flow/material_density_transport.hpp"
#include "hundun/mesh/mesh_geometry.hpp"
#include "hundun/mesh/mesh_topology.hpp"
#include "hundun/runtime/exchange_plan.hpp"
#include "hundun/runtime/field_access_plan.hpp"
#include "hundun/runtime/field_registry.hpp"
#include "hundun/runtime/halo_exchange.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/mpi_environment.hpp"
#include "hundun/runtime/mpi_operation_error.hpp"
#include "hundun/runtime/structured_decomposition.hpp"
#include "tests/support/test_main.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

hundun::runtime::Int3 process_grid(int ranks) {
  if (ranks == 1)
    return {1, 1, 1};
  if (ranks == 2)
    return {2, 1, 1};
  if (ranks == 4)
    return {4, 1, 1};
  throw hundun::runtime::Error("unsupported density-wave rank count");
}

hundun::config::FlowCaseConfig periodic_case() {
  hundun::config::FlowCaseConfig config{};
  config.schema_version = 2;
  config.simulation_type =
      hundun::config::SimulationType::variable_density_flow;
  config.density_model = hundun::config::DensityModel::material;
  config.physics.rho_ref_kg_per_m3 = 1.0;
  config.physics.inlet_consistency_rtol = 1.0e-12;
  config.scalars.push_back({"alpha", 0.0});
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
                                      std::uint32_t components = 1U,
                                      bool conservative = true) {
  return {name,
          unit,
          "material-density-wave",
          hundun::runtime::FunctionSpace::cell_average,
          hundun::runtime::ScalarType::float64,
          components,
          2,
          conservative,
          hundun::runtime::RestartPolicy::persistent,
          hundun::runtime::OutputPolicy::never};
}

hundun::runtime::FieldDescriptor face(const char *name,
                                      std::uint32_t components) {
  return {name,
          "1",
          "material-density-wave",
          hundun::runtime::FunctionSpace::face_value,
          hundun::runtime::ScalarType::float64,
          components,
          0,
          false,
          hundun::runtime::RestartPolicy::persistent,
          hundun::runtime::OutputPolicy::never};
}

double exact_density(double x, double time) {
  constexpr double pi = 3.141592653589793238462643383279502884;
  return 1.0 + 0.2 * std::sin(2.0 * pi * (x - time));
}

double global_sum(const hundun::runtime::MpiContext &mpi, double local) {
  double result{};
  hundun::runtime::check_mpi_result(
      MPI_Allreduce(&local, &result, 1, MPI_DOUBLE, MPI_SUM, mpi.comm()),
      "MPI_Allreduce(density wave test sum)");
  return result;
}

struct WaveResult final {
  double l1{};
  double maximum_mass_error{};
  double minimum_density{};
  std::array<std::vector<std::pair<std::uint64_t, double>>, 4> fields;
};

std::vector<std::pair<std::uint64_t, double>> gather_owned_values(
    const hundun::runtime::MpiContext &mpi,
    const std::vector<std::pair<std::uint64_t, double>> &local) {
  const int local_count = static_cast<int>(local.size());
  std::vector<int> counts(static_cast<std::size_t>(mpi.size()));
  hundun::runtime::check_mpi_result(
      MPI_Allgather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT,
                    mpi.comm()),
      "MPI_Allgather(density wave decomposition counts)");
  std::vector<int> displacements(counts.size());
  int total{};
  for (std::size_t rank = 0; rank < counts.size(); ++rank) {
    displacements[rank] = total;
    total += counts[rank];
  }
  std::vector<std::uint64_t> local_ids(local.size());
  std::vector<double> local_values(local.size());
  for (std::size_t index = 0; index < local.size(); ++index) {
    local_ids[index] = local[index].first;
    local_values[index] = local[index].second;
  }
  std::vector<std::uint64_t> global_ids(static_cast<std::size_t>(total));
  std::vector<double> global_values(static_cast<std::size_t>(total));
  hundun::runtime::check_mpi_result(
      MPI_Allgatherv(local_ids.data(), local_count, MPI_UINT64_T,
                     global_ids.data(), counts.data(), displacements.data(),
                     MPI_UINT64_T, mpi.comm()),
      "MPI_Allgatherv(density wave decomposition IDs)");
  hundun::runtime::check_mpi_result(
      MPI_Allgatherv(local_values.data(), local_count, MPI_DOUBLE,
                     global_values.data(), counts.data(), displacements.data(),
                     MPI_DOUBLE, mpi.comm()),
      "MPI_Allgatherv(density wave decomposition values)");
  std::vector<std::pair<std::uint64_t, double>> result(
      static_cast<std::size_t>(total));
  for (std::size_t index = 0; index < result.size(); ++index)
    result[index] = {global_ids[index], global_values[index]};
  std::sort(result.begin(), result.end(),
            [](const auto &left, const auto &right) {
              return left.first < right.first;
            });
  HUNDUN_CHECK(std::adjacent_find(result.begin(), result.end(),
                                  [](const auto &left, const auto &right) {
                                    return left.first == right.first;
                                  }) == result.end());
  return result;
}

WaveResult run_wave(const hundun::runtime::MpiContext &mpi, int n) {
  using namespace hundun;
  const runtime::Int3 extent{n, 4, 4};
  auto decomposition = runtime::StructuredDecomposition::create(
      mpi, extent, {true, true, true},
      runtime::DecompositionOptions{process_grid(mpi.size())});
  mesh::MeshTopology topology(decomposition);
  mesh::MeshGeometry geometry(
      topology, mesh::UniformBoxMapping({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}));
  auto boundaries =
      boundary::BoundaryRegistry::create(periodic_case(), topology);
  runtime::FieldRegistry registry;
  flow::FlowFieldIds fields;
  fields.density = registry.declare_field(cell("rho", "kg/m3"));
  fields.velocity = registry.declare_field(cell("velocity", "1", 3U, false));
  fields.mechanical_pressure =
      registry.declare_field(cell("pi", "Pa", 1U, false));
  fields.face_velocity = registry.declare_field(face("u_face", 3U));
  fields.face_mass_flux = finite_volume::declare_face_mass_flux(registry);
  const auto rho_h = registry.declare_field(cell("rho_h", "J/m3"));
  const auto rho_phi = registry.declare_field(cell("rho_phi", "kg/m3"));
  fields.transported_cell_fields = {rho_h, rho_phi};
  registry.freeze();
  const auto box = decomposition.owned_box();
  const runtime::Int3 local{box.end.x - box.begin.x, box.end.y - box.begin.y,
                            box.end.z - box.begin.z};
  const double dt = 0.2 / static_cast<double>(n);
  auto state = flow::FlowState::create(
      registry, {local, topology.local_face_count()}, fields,
      {0U, 0.0, dt, dt, flow::MomentumTimeOrder::bdf2});
  flow::FlowLayerValues history;
  flow::FlowLayerValues committed;
  for (mesh::LocalCellId cell_id = 0; cell_id < topology.owned_cell_count();
       ++cell_id) {
    const double x = geometry.cell_center_m(cell_id).x;
    history.density.push_back(exact_density(x, -dt));
    committed.density.push_back(exact_density(x, 0.0));
    for (int component = 0; component < 3; ++component) {
      history.velocity.push_back(component == 0 ? 1.0 : 0.0);
      committed.velocity.push_back(component == 0 ? 1.0 : 0.0);
    }
    history.mechanical_pressure.push_back(0.0);
    committed.mechanical_pressure.push_back(0.0);
  }
  history.face_velocity.assign(topology.local_face_count() * 3U, 0.0);
  committed.face_velocity = history.face_velocity;
  history.face_mass_flux.assign(topology.local_face_count(), 0.0);
  committed.face_mass_flux = history.face_mass_flux;
  history.transported_cell_fields.resize(2U);
  committed.transported_cell_fields.resize(2U);
  for (const double rho : history.density) {
    history.transported_cell_fields[0].push_back(3.0 * rho);
    history.transported_cell_fields[1].push_back(0.4 * rho);
  }
  for (const double rho : committed.density) {
    committed.transported_cell_fields[0].push_back(3.0 * rho);
    committed.transported_cell_fields[1].push_back(0.4 * rho);
  }
  state.seed_accepted_layers(history, committed);

  runtime::FieldStorage flux_storage(
      registry, runtime::FieldLayoutSet{local, topology.local_face_count()});
  constexpr runtime::PhaseId phase = 1911U;
  constexpr runtime::ActorId actor = 1911U;
  runtime::FieldAccessPlan flux_access(registry);
  flux_access.declare_access(phase, actor, fields.face_mass_flux,
                             runtime::AccessMode::read_write);
  flux_access.freeze();
  auto halo = runtime::HaloExchange::create(
      decomposition, runtime::ExchangePlan::create(decomposition, local, 2));
  flow::MaterialDensityTransportSpec spec;
  spec.enthalpy_density = rho_h;
  spec.scalar_densities = {rho_phi};
  spec.scalar_diffusivities_kg_per_m_s = {0.0};
  auto transport = flow::MaterialDensityTransport::create(
      registry, decomposition, topology, geometry, boundaries, mpi, halo,
      fields, spec);
  const auto stencil =
      flow::make_momentum_time_stencil(flow::MomentumTimeOrder::bdf2, dt, dt);
  WaveResult result;
  result.minimum_density = 1.0;
  const int steps = n / 2;
  for (int step = 0; step < steps; ++step) {
    const double next_time = state.metadata().time_s + dt;
    auto flux_values = flux_storage.acquire_face_write<double>(
        flux_access, phase, actor, fields.face_mass_flux);
    for (mesh::LocalFaceId face_id = 0; face_id < topology.local_face_count();
         ++face_id) {
      const auto center = geometry.face_center_m(face_id);
      const auto area =
          geometry.face_area_vector_m2(face_id, mesh::FaceSide::owner);
      flux_values(face_id, 0) = exact_density(center.x, next_time) * area.x;
    }
    state.begin_attempt();
    auto flux = flow::MaterialFaceMassFlux::acquire(
        registry, flux_storage, flux_access, phase, actor,
        fields.face_mass_flux, topology,
        flow::MaterialFluxProvenance::final_corrected);
    const auto report = transport.finalize_trial(state, flux, stencil);
    if (report.disposition() != flow::MaterialTransportDisposition::finalized &&
        mpi.rank() == 0)
      std::cerr << "density-wave failure n=" << n << " step=" << step
                << " reason=" << static_cast<int>(report.reason())
                << " rho-res=" << report.density_normalized_l2()
                << " mass-defect=" << report.mass_relative_conservation_defect()
                << '\n';
    HUNDUN_CHECK(report.disposition() ==
                 flow::MaterialTransportDisposition::finalized);
    HUNDUN_CHECK(report.density_normalized_l2() <= 1.0e-10);
    HUNDUN_CHECK(report.mass_relative_conservation_defect() <= 5.0e-11);
    for (const auto value : report.transport_normalized_l2())
      HUNDUN_CHECK(value <= 1.0e-9);
    for (const auto value : report.transport_relative_conservation_defect())
      HUNDUN_CHECK(value <= 5.0e-11);
    const auto trial = state.snapshot(flow::FlowLayer::trial);
    double local_mass{};
    for (mesh::LocalCellId cell_id = 0; cell_id < topology.owned_cell_count();
         ++cell_id) {
      local_mass += geometry.cell_volume_m3(cell_id) * trial.density[cell_id];
      result.minimum_density =
          std::min(result.minimum_density, trial.density[cell_id]);
    }
    const double mass = global_sum(mpi, local_mass);
    result.maximum_mass_error =
        std::max(result.maximum_mass_error, std::abs(mass - 1.0));
    HUNDUN_CHECK(result.maximum_mass_error <= 5.0e-12);
    state.commit_attempt({static_cast<std::uint64_t>(step + 1), next_time, dt,
                          dt, flow::MomentumTimeOrder::bdf2});
  }
  const auto final = state.snapshot(flow::FlowLayer::committed);
  const double final_time = state.metadata().time_s;
  double local_error{};
  double local_volume{};
  for (mesh::LocalCellId cell_id = 0; cell_id < topology.owned_cell_count();
       ++cell_id) {
    const double volume = geometry.cell_volume_m3(cell_id);
    const double exact =
        exact_density(geometry.cell_center_m(cell_id).x, final_time);
    local_error += volume * std::abs(final.density[cell_id] - exact);
    local_volume += volume;
    HUNDUN_CHECK(final.density[cell_id] > 0.0);
    HUNDUN_CHECK(std::abs(final.transported_cell_fields[0][cell_id] -
                          3.0 * final.density[cell_id]) <= 1.0e-11);
    HUNDUN_CHECK(std::abs(final.transported_cell_fields[1][cell_id] -
                          0.4 * final.density[cell_id]) <= 1.0e-11);
    const auto global = topology.global_cell_id(cell_id);
    result.fields[0].push_back({global, final.density[cell_id]});
    result.fields[1].push_back(
        {global, final.transported_cell_fields[0][cell_id]});
    result.fields[2].push_back(
        {global, final.transported_cell_fields[1][cell_id]});
  }
  for (mesh::LocalFaceId face_id = 0; face_id < topology.local_face_count();
       ++face_id) {
    if (topology.face_ownership(face_id) == mesh::EntityOwnership::owned)
      result.fields[3].push_back(
          {topology.global_face_id(face_id), final.face_mass_flux[face_id]});
  }
  for (auto &field : result.fields)
    field = gather_owned_values(mpi, field);
  result.l1 = global_sum(mpi, local_error) / global_sum(mpi, local_volume);
  return result;
}

void require_decomposition_match(const WaveResult &serial,
                                 const WaveResult &distributed) {
  for (std::size_t field = 0; field < serial.fields.size(); ++field) {
    HUNDUN_CHECK(serial.fields[field].size() ==
                 distributed.fields[field].size());
    double infinity_norm{};
    double maximum_difference{};
    for (std::size_t index = 0; index < serial.fields[field].size(); ++index) {
      HUNDUN_CHECK(serial.fields[field][index].first ==
                   distributed.fields[field][index].first);
      infinity_norm =
          std::max({infinity_norm, std::abs(serial.fields[field][index].second),
                    std::abs(distributed.fields[field][index].second)});
      maximum_difference =
          std::max(maximum_difference,
                   std::abs(serial.fields[field][index].second -
                            distributed.fields[field][index].second));
    }
    HUNDUN_CHECK(maximum_difference <= 5.0e-12 * std::max(1.0, infinity_norm));
  }
}

void run(const hundun::runtime::MpiContext &mpi) {
  const auto serial_mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_SELF);
  const auto serial_coarse = run_wave(serial_mpi, 32);
  const auto coarse = run_wave(mpi, 32);
  require_decomposition_match(serial_coarse, coarse);
  const auto serial_medium = run_wave(serial_mpi, 64);
  const auto medium = run_wave(mpi, 64);
  require_decomposition_match(serial_medium, medium);
  const auto serial_fine = run_wave(serial_mpi, 128);
  const auto fine = run_wave(mpi, 128);
  require_decomposition_match(serial_fine, fine);
  const double first_order = std::log2(coarse.l1 / medium.l1);
  const double second_order = std::log2(medium.l1 / fine.l1);
  if (mpi.rank() == 0)
    std::cout << "density-wave l1=" << coarse.l1 << ',' << medium.l1 << ','
              << fine.l1 << " order=" << first_order << ',' << second_order
              << '\n';
  HUNDUN_CHECK(first_order >= 1.8);
  HUNDUN_CHECK(second_order >= 1.8);
  HUNDUN_CHECK(coarse.maximum_mass_error <= 5.0e-12);
  HUNDUN_CHECK(medium.maximum_mass_error <= 5.0e-12);
  HUNDUN_CHECK(fine.maximum_mass_error <= 5.0e-12);
  HUNDUN_CHECK(coarse.minimum_density > 0.0);
  HUNDUN_CHECK(medium.minimum_density > 0.0);
  HUNDUN_CHECK(fine.minimum_density > 0.0);
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  return hundun::test::run([&] {
    auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
    run(mpi);
  });
}
