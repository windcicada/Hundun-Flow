// SPDX-License-Identifier: Apache-2.0

#include "hundun/boundary/basic_boundary.hpp"
#include "hundun/finite_volume/cell_centered_fvm.hpp"
#include "hundun/flow/checkpoint_v2.hpp"
#include "hundun/mesh/mesh_geometry.hpp"
#include "hundun/mesh/mesh_topology.hpp"
#include "hundun/runtime/field_registry.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/mpi_environment.hpp"
#include "hundun/runtime/structured_decomposition.hpp"
#include "checkpoint_v2_test_access.hpp"
#include "tests/support/flow_state_equality.hpp"
#include "tests/support/test_main.hpp"

#include <array>
#include <filesystem>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

namespace {

template <class Function>
bool rejects(Function &&function) {
  try {
    function();
  } catch (const std::exception &) {
    return true;
  }
  return false;
}

hundun::runtime::Int3 grid(int ranks) { return {ranks, 1, 1}; }

hundun::config::FlowCaseConfig make_case(int ranks) {
  hundun::config::FlowCaseConfig config{};
  config.schema_version = 2;
  config.case_name = "checkpoint-v2";
  config.simulation_type =
      hundun::config::SimulationType::variable_density_flow;
  config.density_model = hundun::config::DensityModel::material;
  config.resources.expected_ranks = ranks;
  config.resources.process_grid = grid(ranks);
  config.mesh.cells = {8, 4, 3};
  config.mesh.origin_m = {0.0, 0.0, 0.0};
  config.mesh.length_m = {1.0, 1.0, 1.0};
  config.mesh.mapping = hundun::config::MeshMapping::uniform_box;
  config.time = {hundun::config::TimeMode::fixed, 4, 0.01, 0.001, 0.1,
                 0.5, 0.25, 1.25, 0.5, 8};
  config.physics.rho_ref_kg_per_m3 = 1.0;
  config.physics.dynamic_viscosity_pa_s = 0.01;
  config.physics.inlet_consistency_rtol = 1.0e-12;
  config.scalars = {{"alpha", 0.0}};
  constexpr std::array names{
      hundun::config::PatchName::x_min, hundun::config::PatchName::x_max,
      hundun::config::PatchName::y_min, hundun::config::PatchName::y_max,
      hundun::config::PatchName::z_min, hundun::config::PatchName::z_max};
  for (std::size_t index = 0; index < names.size(); ++index) {
    config.boundaries[index].patch = names[index];
    config.boundaries[index].type =
        hundun::config::BoundaryType::periodic;
  }
  return config;
}

hundun::runtime::FieldDescriptor cell(const char *name, std::uint32_t n,
                                      bool conservative) {
  return {name,
          "1",
          "checkpoint-v2",
          hundun::runtime::FunctionSpace::cell_average,
          hundun::runtime::ScalarType::float64,
          n,
          2,
          conservative,
          hundun::runtime::RestartPolicy::persistent,
          hundun::runtime::OutputPolicy::never};
}
hundun::runtime::FieldDescriptor face(const char *name, std::uint32_t n) {
  return {name,
          "1",
          "checkpoint-v2",
          hundun::runtime::FunctionSpace::face_value,
          hundun::runtime::ScalarType::float64,
          n,
          0,
          false,
          hundun::runtime::RestartPolicy::persistent,
          hundun::runtime::OutputPolicy::never};
}
hundun::runtime::FieldDescriptor physical_cell(const char *name,
                                               const char *unit,
                                               bool conservative) {
  auto result = cell(name, 1U, conservative);
  result.unit = unit;
  return result;
}

void run(const hundun::runtime::MpiContext &mpi) {
  HUNDUN_CHECK(
      hundun::test::checkpoint_v2_state_equality_oracle_is_mutation_sensitive());
  auto config = make_case(mpi.size());
  auto decomposition = hundun::runtime::StructuredDecomposition::create(
      mpi, config.mesh.cells, {true, true, true},
      hundun::runtime::DecompositionOptions{grid(mpi.size())});
  hundun::mesh::MeshTopology topology(decomposition);
  hundun::mesh::MeshGeometry geometry(
      topology, hundun::mesh::UniformBoxMapping(config.mesh.origin_m,
                                                config.mesh.length_m));
  auto boundaries =
      hundun::boundary::BoundaryRegistry::create(config, topology);
  hundun::runtime::FieldRegistry registry;
  hundun::flow::FlowFieldIds fields;
  fields.density = registry.declare_field(cell("rho", 1U, true));
  fields.velocity = registry.declare_field(cell("u", 3U, false));
  fields.mechanical_pressure =
      registry.declare_field(cell("pi", 1U, false));
  fields.face_velocity = registry.declare_field(face("uf", 3U));
  fields.face_mass_flux =
      hundun::finite_volume::declare_face_mass_flux(registry);
  fields.transported_cell_fields = {
      registry.declare_field(cell("rho_alpha", 1U, true))};
  registry.freeze();
  const auto metadata = hundun::flow::AcceptedStepMetadata{
      0U, 0.0, config.time.initial_dt_s, 0.0,
      hundun::flow::MomentumTimeOrder::backward_euler};
  auto make_state = [&] {
    return hundun::flow::FlowState::create(
        registry,
        {decomposition.local_extent(), topology.local_face_count()}, fields,
        metadata);
  };
  const std::size_t cells = topology.owned_cell_count();
  hundun::flow::FlowLayerValues history;
  history.density.resize(cells);
  history.velocity.resize(cells * 3U);
  history.mechanical_pressure.resize(cells);
  history.face_velocity.resize(topology.local_face_count() * 3U);
  history.face_mass_flux.resize(topology.local_face_count());
  history.transported_cell_fields = {std::vector<double>(cells)};
  for (std::size_t cell_id = 0; cell_id < cells; ++cell_id) {
    history.density[cell_id] =
        1.0 + 0.001 *
                  (static_cast<double>(mpi.rank()) +
                   static_cast<double>(cell_id));
    history.mechanical_pressure[cell_id] =
        -0.25 + static_cast<double>(cell_id);
    history.transported_cell_fields[0][cell_id] =
        history.density[cell_id] * 0.2;
    for (std::size_t component = 0; component < 3U; ++component)
      history.velocity[cell_id * 3U + component] =
          10.0 * static_cast<double>(mpi.rank()) +
          3.0 * static_cast<double>(cell_id) +
          static_cast<double>(component);
  }
  for (std::size_t face_id = 0; face_id < topology.local_face_count();
       ++face_id) {
    history.face_mass_flux[face_id] =
        (static_cast<double>(mpi.rank()) + static_cast<double>(face_id)) *
        0.01;
    for (std::size_t component = 0; component < 3U; ++component)
      history.face_velocity[face_id * 3U + component] =
          static_cast<double>(face_id + component) * 0.02;
  }
  auto committed = history;
  for (double &item : committed.density)
    item += 0.25;
  for (double &item : committed.transported_cell_fields[0])
    item += 0.125;

  auto source = make_state();
  source.seed_accepted_layers(history, committed);
  auto controller = hundun::flow::Bdf2RetryController::create(
      config.time, config.density_model, topology, geometry, mpi, source);
  const auto controller_state = controller.state();

  const auto directory =
      std::filesystem::temp_directory_path() /
      ("hundun-task23-checkpoint-" + std::to_string(mpi.size()));
  if (mpi.rank() == 0)
    std::filesystem::remove_all(directory);
  mpi.barrier();
  const auto written = hundun::flow::write_checkpoint_v2(
      mpi, decomposition, topology, geometry, boundaries, config, source,
      controller_state, std::nullopt, directory);
  HUNDUN_CHECK(written.disposition() ==
               hundun::flow::CheckpointV2Disposition::completed);
  HUNDUN_CHECK(written.file_count() ==
               static_cast<std::uint64_t>(mpi.size()) + 2U);
  HUNDUN_CHECK(std::filesystem::is_regular_file(directory / "COMPLETED"));

  auto destination = make_state();
  auto different = history;
  for (double &item : different.density)
    item += 9.0;
  destination.seed_accepted_layers(different, different);
  const auto restored = hundun::flow::read_checkpoint_v2(
      mpi, decomposition, topology, geometry, boundaries, config, destination,
      directory);
  HUNDUN_CHECK(restored.restored());
  HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
      history, destination.snapshot(hundun::flow::FlowLayer::history)));
  HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
      committed, destination.snapshot(hundun::flow::FlowLayer::committed)));
  HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
      committed, destination.snapshot(hundun::flow::FlowLayer::trial)));
  HUNDUN_CHECK(hundun::test::accepted_step_metadata_bitwise_equal(
      metadata, destination.metadata()));
  HUNDUN_CHECK(hundun::test::time_control_state_bitwise_equal(
      controller_state, restored.time_control_state()));

  using CheckpointAccess = hundun::flow::test::CheckpointV2TestAccess;
  HUNDUN_CHECK(
      hundun::flow::test::
          checkpoint_v2_deep_snapshot_oracle_is_mutation_sensitive(
              destination));
  for (std::size_t generation_layer = 0; generation_layer < 4U;
       ++generation_layer) {
    auto generation_state = make_state();
    generation_state.seed_accepted_layers(different, different);
    CheckpointAccess::force_generation(
        generation_state, generation_layer,
        std::numeric_limits<std::uint64_t>::max());
    const auto generation_before =
        CheckpointAccess::snapshot(generation_state);
    const auto generation_views =
        CheckpointAccess::density_views(generation_state);
    const auto rejected_entry = hundun::flow::read_checkpoint_v2(
        mpi, decomposition, topology, geometry, boundaries, config,
        generation_state, directory);
    HUNDUN_CHECK(!rejected_entry.restored());
    HUNDUN_CHECK(rejected_entry.report().reason() ==
                 hundun::flow::CheckpointV2FailureReason::state);
    HUNDUN_CHECK(rejected_entry.report().phase() ==
                 hundun::flow::CheckpointV2Phase::transaction_entry);
    HUNDUN_CHECK(
        rejected_entry.report().transaction_entry_status() ==
        hundun::flow::CheckpointV2CheckStatus::failed);
    HUNDUN_CHECK(hundun::flow::test::checkpoint_v2_deep_snapshot_equal(
        generation_before, CheckpointAccess::snapshot(generation_state)));
    for (const auto &view : generation_views)
      HUNDUN_CHECK(std::isfinite(view(0, 0, 0, 0)));
  }

  if (mpi.size() > 1) {
    auto partition_config = config;
    partition_config.resources.process_grid =
        hundun::runtime::Int3{1, mpi.size(), 1};
    auto partition_decomposition =
        hundun::runtime::StructuredDecomposition::create(
            mpi, partition_config.mesh.cells, {true, true, true},
            hundun::runtime::DecompositionOptions{
                *partition_config.resources.process_grid});
    hundun::mesh::MeshTopology partition_topology(partition_decomposition);
    hundun::mesh::MeshGeometry partition_geometry(
        partition_topology,
        hundun::mesh::UniformBoxMapping(partition_config.mesh.origin_m,
                                        partition_config.mesh.length_m));
    auto partition_boundaries = hundun::boundary::BoundaryRegistry::create(
        partition_config, partition_topology);
    auto partition_state = hundun::flow::FlowState::create(
        registry,
        {partition_decomposition.local_extent(),
         partition_topology.local_face_count()},
        fields, metadata);
    const auto partition_values =
        [&] {
          auto values = history;
          const auto partition_cells = partition_topology.owned_cell_count();
          values.density.assign(partition_cells, 2.0);
          values.velocity.assign(partition_cells * 3U, 0.0);
          values.mechanical_pressure.assign(partition_cells, 0.0);
          values.face_velocity.assign(
              partition_topology.local_face_count() * 3U, 0.0);
          values.face_mass_flux.assign(
              partition_topology.local_face_count(), 0.0);
          values.transported_cell_fields = {
              std::vector<double>(partition_cells, 0.4)};
          return values;
        }();
    partition_state.seed_accepted_layers(partition_values, partition_values);
    const auto partition_before = CheckpointAccess::snapshot(partition_state);
    const auto partition_views =
        CheckpointAccess::density_views(partition_state);
    const auto partition_read = hundun::flow::read_checkpoint_v2(
        mpi, partition_decomposition, partition_topology, partition_geometry,
        partition_boundaries, partition_config, partition_state, directory);
    HUNDUN_CHECK(!partition_read.restored());
    HUNDUN_CHECK(partition_read.report().reason() ==
                 hundun::flow::CheckpointV2FailureReason::layout);
    HUNDUN_CHECK(partition_read.report().partition_status() ==
                 hundun::flow::CheckpointV2CheckStatus::failed);
    HUNDUN_CHECK(
        hundun::flow::test::checkpoint_v2_failed_read_preserved_values(
            partition_before, CheckpointAccess::snapshot(partition_state)));
    for (const auto &view : partition_views)
      HUNDUN_CHECK(rejects(
          [&] { static_cast<void>(view(0, 0, 0, 0)); }));
  }

  const auto source_view =
      source.layer(hundun::flow::FlowLayer::committed)
          .view<double>(fields.density);
  const double source_first = source_view(0, 0, 0, 0);
  const auto duplicate_write = hundun::flow::write_checkpoint_v2(
      mpi, decomposition, topology, geometry, boundaries, config, source,
      controller_state, std::nullopt, directory);
  HUNDUN_CHECK(duplicate_write.disposition() ==
               hundun::flow::CheckpointV2Disposition::failed);
  HUNDUN_CHECK(source_view(0, 0, 0, 0) == source_first);

  CheckpointAccess::set_committed_density_ghost(destination, -0.0);
  CheckpointAccess::set_rollback_density_ghost(destination, 7.25);
  const auto before_deep = CheckpointAccess::snapshot(destination);
  const auto old_views = CheckpointAccess::density_views(destination);
  mpi.barrier();
  if (mpi.rank() == 0) {
    const auto manifest = directory / "manifest.v2.bin";
    std::fstream stream(manifest, std::ios::in | std::ios::out |
                                      std::ios::binary);
    HUNDUN_CHECK(static_cast<bool>(stream));
    stream.seekg(24);
    char byte{};
    stream.read(&byte, 1);
    HUNDUN_CHECK(static_cast<bool>(stream));
    byte ^= 1;
    stream.seekp(24);
    stream.write(&byte, 1);
    stream.flush();
    HUNDUN_CHECK(static_cast<bool>(stream));
  }
  mpi.barrier();
  const auto failed_read = hundun::flow::read_checkpoint_v2(
      mpi, decomposition, topology, geometry, boundaries, config, destination,
      directory);
  HUNDUN_CHECK(!failed_read.restored());
  HUNDUN_CHECK(failed_read.report().reason() ==
               hundun::flow::CheckpointV2FailureReason::file_integrity);
  HUNDUN_CHECK(
      hundun::flow::test::checkpoint_v2_failed_read_preserved_values(
          before_deep, CheckpointAccess::snapshot(destination)));
  for (const auto &view : old_views)
    HUNDUN_CHECK(
        rejects([&] { static_cast<void>(view(0, 0, 0, 0)); }));

  auto ideal_config = config;
  ideal_config.density_model = hundun::config::DensityModel::ideal_gas;
  ideal_config.scalars.clear();
  ideal_config.physics.cp_J_per_kg_K = 1000.0;
  ideal_config.physics.gas_constant_J_per_kg_K = 287.05;
  ideal_config.physics.thermodynamic_pressure_pa = 101325.0;
  auto ideal_boundaries =
      hundun::boundary::BoundaryRegistry::create(ideal_config, topology);
  hundun::runtime::FieldRegistry ideal_registry;
  hundun::flow::FlowFieldIds ideal_fields;
  ideal_fields.density = ideal_registry.declare_field(
      physical_cell("rho", "kg/m3", true));
  ideal_fields.velocity =
      ideal_registry.declare_field(cell("u", 3U, false));
  ideal_fields.mechanical_pressure =
      ideal_registry.declare_field(cell("pi", 1U, false));
  ideal_fields.face_velocity =
      ideal_registry.declare_field(face("uf", 3U));
  ideal_fields.face_mass_flux =
      hundun::finite_volume::declare_face_mass_flux(ideal_registry);
  const auto rho_h = ideal_registry.declare_field(
      physical_cell("rho_h", "J/m3", true));
  ideal_fields.transported_cell_fields = {rho_h};
  ideal_registry.freeze();
  auto ideal_state = hundun::flow::FlowState::create(
      ideal_registry,
      {decomposition.local_extent(), topology.local_face_count()},
      ideal_fields, metadata);
  constexpr double temperature = 300.0;
  const double ideal_density =
      101325.0 / (287.05 * temperature);
  hundun::flow::FlowLayerValues ideal_values;
  ideal_values.density.assign(cells, ideal_density);
  ideal_values.velocity.assign(cells * 3U, 0.0);
  ideal_values.mechanical_pressure.assign(cells, 0.0);
  ideal_values.face_velocity.assign(topology.local_face_count() * 3U, 0.0);
  ideal_values.face_mass_flux.assign(topology.local_face_count(), 0.0);
  ideal_values.transported_cell_fields = {std::vector<double>(
      cells, ideal_density * 1000.0 * temperature)};
  ideal_state.seed_accepted_layers(ideal_values, ideal_values);
  const hundun::flow::IdealGasClosureSpec ideal_spec{
      rho_h, 1000.0, 287.05, 101325.0};
  auto initial_closure = hundun::flow::IdealGasClosure::create(
      topology, geometry, ideal_boundaries, mpi, ideal_registry, ideal_fields,
      ideal_state, ideal_spec);
  auto persisted_closure = initial_closure.state();
  persisted_closure.revision = 4U;
  auto restored_closure = hundun::flow::IdealGasClosure::restore(
      topology, geometry, ideal_boundaries, mpi, ideal_registry, ideal_fields,
      ideal_state, ideal_spec, persisted_closure);
  HUNDUN_CHECK(hundun::test::ideal_gas_closure_state_bitwise_equal(
      restored_closure.state(), persisted_closure));
  auto invalid_closure = persisted_closure;
  invalid_closure.revision = std::numeric_limits<std::uint64_t>::max();
  HUNDUN_CHECK(rejects([&] {
    static_cast<void>(hundun::flow::IdealGasClosure::restore(
        topology, geometry, ideal_boundaries, mpi, ideal_registry,
        ideal_fields, ideal_state, ideal_spec, invalid_closure));
  }));

  auto ideal_controller = hundun::flow::Bdf2RetryController::create(
      ideal_config.time, ideal_config.density_model, topology, geometry, mpi,
      ideal_state);
  const auto ideal_controller_state = ideal_controller.state();
  const auto ideal_directory =
      std::filesystem::temp_directory_path() /
      ("hundun-task23-ideal-checkpoint-" + std::to_string(mpi.size()));
  if (mpi.rank() == 0)
    std::filesystem::remove_all(ideal_directory);
  mpi.barrier();
  const auto ideal_written = hundun::flow::write_checkpoint_v2(
      mpi, decomposition, topology, geometry, ideal_boundaries, ideal_config,
      ideal_state, ideal_controller_state, persisted_closure,
      ideal_directory);
  HUNDUN_CHECK(ideal_written.disposition() ==
               hundun::flow::CheckpointV2Disposition::completed);
  auto ideal_destination = hundun::flow::FlowState::create(
      ideal_registry,
      {decomposition.local_extent(), topology.local_face_count()},
      ideal_fields, metadata);
  ideal_destination.seed_accepted_layers(ideal_values, ideal_values);
  const auto ideal_read = hundun::flow::read_checkpoint_v2(
      mpi, decomposition, topology, geometry, ideal_boundaries, ideal_config,
      ideal_destination, ideal_directory);
  HUNDUN_CHECK(ideal_read.restored());
  HUNDUN_CHECK(ideal_read.ideal_gas_closure_state_available());
  HUNDUN_CHECK(hundun::test::ideal_gas_closure_state_bitwise_equal(
      ideal_read.ideal_gas_closure_state(), persisted_closure));
  HUNDUN_CHECK(hundun::test::time_control_state_bitwise_equal(
      ideal_read.time_control_state(), ideal_controller_state));
  auto resumed_closure = hundun::flow::IdealGasClosure::restore(
      topology, geometry, ideal_boundaries, mpi, ideal_registry, ideal_fields,
      ideal_destination, ideal_spec,
      ideal_read.ideal_gas_closure_state());
  HUNDUN_CHECK(hundun::test::ideal_gas_closure_state_bitwise_equal(
      resumed_closure.state(), persisted_closure));

  auto invalid_ideal_values = ideal_values;
  invalid_ideal_values.density.front() *= 1.01;
  auto invalid_ideal_state = hundun::flow::FlowState::create(
      ideal_registry,
      {decomposition.local_extent(), topology.local_face_count()},
      ideal_fields, metadata);
  invalid_ideal_state.seed_accepted_layers(invalid_ideal_values,
                                           invalid_ideal_values);
  const auto invalid_directory =
      std::filesystem::temp_directory_path() /
      ("hundun-task23-invalid-ideal-" + std::to_string(mpi.size()));
  if (mpi.rank() == 0)
    std::filesystem::remove_all(invalid_directory);
  mpi.barrier();
  const auto invalid_write = hundun::flow::write_checkpoint_v2(
      mpi, decomposition, topology, geometry, ideal_boundaries, ideal_config,
      invalid_ideal_state, ideal_controller_state, persisted_closure,
      invalid_directory);
  HUNDUN_CHECK(invalid_write.disposition() ==
               hundun::flow::CheckpointV2Disposition::failed);
  HUNDUN_CHECK(invalid_write.reason() ==
               hundun::flow::CheckpointV2FailureReason::state);
  HUNDUN_CHECK(!std::filesystem::exists(invalid_directory));

  mpi.barrier();
  if (mpi.rank() == 0) {
    std::filesystem::remove_all(directory);
    std::filesystem::remove_all(ideal_directory);
    std::filesystem::remove_all(invalid_directory);
  }
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  return hundun::test::run([&] {
    auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
    HUNDUN_CHECK(mpi.size() == 1 || mpi.size() == 2 || mpi.size() == 4);
    run(mpi);
  });
}
