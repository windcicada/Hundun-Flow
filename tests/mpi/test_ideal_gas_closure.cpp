// SPDX-License-Identifier: Apache-2.0

#include "flow/src/ideal_gas_closure_test_access.hpp"
#include "hundun/boundary/basic_boundary.hpp"
#include "hundun/config/resolved_case.hpp"
#include "hundun/finite_volume/cell_centered_fvm.hpp"
#include "hundun/flow/ideal_gas_closure.hpp"
#include "hundun/mesh/mesh_geometry.hpp"
#include "hundun/mesh/mesh_topology.hpp"
#include "hundun/runtime/field_registry.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/mpi_environment.hpp"
#include "hundun/runtime/mpi_operation_error.hpp"
#include "hundun/runtime/structured_decomposition.hpp"
#include "tests/support/flow_state_equality.hpp"
#include "tests/support/test_main.hpp"

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string_view>
#include <utility>
#include <vector>

namespace {

hundun::runtime::Int3 grid(int ranks) { return {ranks, 1, 1}; }

void require_fp64_delta(
    const hundun::runtime::Fp64ReductionCounters &before,
    const hundun::runtime::Fp64ReductionCounters &after,
    std::uint64_t calls, std::uint64_t scalars, std::uint64_t bytes) {
  HUNDUN_CHECK(after.collective_calls - before.collective_calls == calls);
  HUNDUN_CHECK(after.reduced_scalars - before.reduced_scalars == scalars);
  HUNDUN_CHECK(after.logical_payload_bytes - before.logical_payload_bytes ==
               bytes);
}

void require_candidate_oracle(
    const hundun::flow::IdealGasClosureReport &report) {
  HUNDUN_CHECK(hundun::flow::test::IdealGasClosureTestAccess::
                   candidate_pressure_mutation_rejected(report, true));
  HUNDUN_CHECK(hundun::flow::test::IdealGasClosureTestAccess::
                   candidate_pressure_mutation_rejected(report, false));
}

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

void run(const hundun::runtime::MpiContext &mpi) {
  // Deliberately uneven on 2/4 ranks: rank-local layout identity is allowed
  // to differ while the global box still has one exact owner per cell.
  constexpr hundun::runtime::Int3 extent{17, 4, 4};
  auto decomposition = hundun::runtime::StructuredDecomposition::create(
      mpi, extent, {true, true, true},
      hundun::runtime::DecompositionOptions{grid(mpi.size())});
  hundun::mesh::MeshTopology topology(decomposition);
  hundun::mesh::MeshGeometry geometry(
      topology,
      hundun::mesh::UniformBoxMapping({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}));
  auto alternate_decomposition =
      hundun::runtime::StructuredDecomposition::create(
          mpi, {18, 4, 4}, {true, true, true},
          hundun::runtime::DecompositionOptions{grid(mpi.size())});
  hundun::mesh::MeshTopology alternate_topology(alternate_decomposition);
  hundun::mesh::MeshGeometry alternate_geometry(
      alternate_topology,
      hundun::mesh::UniformBoxMapping({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}));
  auto boundaries =
      hundun::boundary::BoundaryRegistry::create(periodic_case(), topology);

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

  auto state = hundun::flow::FlowState::create(
      registry, {decomposition.local_extent(), topology.local_face_count()},
      fields,
      {0U, 0.0, 1.0e-3, 0.0, hundun::flow::MomentumTimeOrder::backward_euler});
  constexpr double cp = 1000.0;
  constexpr double gas_constant = 287.05;
  constexpr double temperature = 300.0;
  constexpr double pressure = 101325.0;
  const double density = pressure / (gas_constant * temperature);
  hundun::flow::FlowLayerValues initial;
  initial.density.assign(topology.owned_cell_count(), density);
  initial.velocity.assign(topology.owned_cell_count() * 3U, 0.0);
  initial.mechanical_pressure.assign(topology.owned_cell_count(), 0.0);
  initial.face_velocity.assign(topology.local_face_count() * 3U, 0.0);
  initial.face_mass_flux.assign(topology.local_face_count(), 0.0);
  initial.transported_cell_fields = {std::vector<double>(
      topology.owned_cell_count(), density * cp * temperature)};
  state.seed_accepted_layers(initial, initial);

  const auto expect_create_rejected = [&](hundun::flow::IdealGasClosureSpec spec) {
    bool rejected = false;
    try {
      static_cast<void>(hundun::flow::IdealGasClosure::create(
          topology, geometry, boundaries, mpi, registry, fields, state, spec));
    } catch (const hundun::runtime::Error &) {
      rejected = true;
    }
    HUNDUN_CHECK(rejected);
  };
  for (const int target : {0, mpi.size() - 1}) {
    auto invalid_cp = hundun::flow::IdealGasClosureSpec{
        rho_h, cp, gas_constant, pressure};
    if (mpi.rank() == target)
      invalid_cp.cp_J_per_kg_K = 0.0;
    expect_create_rejected(invalid_cp);

    auto invalid_r = hundun::flow::IdealGasClosureSpec{
        rho_h, cp, gas_constant, pressure};
    if (mpi.rank() == target)
      invalid_r.gas_constant_J_per_kg_K =
          std::numeric_limits<double>::quiet_NaN();
    expect_create_rejected(invalid_r);

    auto invalid_pressure = hundun::flow::IdealGasClosureSpec{
        rho_h, cp, gas_constant, pressure};
    if (mpi.rank() == target)
      invalid_pressure.configured_thermodynamic_pressure_pa =
          std::numeric_limits<double>::infinity();
    expect_create_rejected(invalid_pressure);

    auto invalid_enthalpy = hundun::flow::IdealGasClosureSpec{
        rho_h, cp, gas_constant, pressure};
    if (mpi.rank() == target)
      invalid_enthalpy.enthalpy_density = fields.density;
    expect_create_rejected(invalid_enthalpy);

    for (const bool geometry_mismatch : {false, true}) {
      const auto &candidate_topology =
          !geometry_mismatch && mpi.rank() == target ? alternate_topology
                                                     : topology;
      const auto &candidate_geometry =
          mpi.rank() == target ? alternate_geometry : geometry;
      bool collaborator_rejected = false;
      try {
        static_cast<void>(hundun::flow::IdealGasClosure::create(
            candidate_topology, candidate_geometry, boundaries, mpi, registry,
            fields, state, {rho_h, cp, gas_constant, pressure}));
      } catch (const hundun::runtime::Error &) {
        collaborator_rejected = true;
      }
      HUNDUN_CHECK(collaborator_rejected);
    }

    const std::size_t face_count =
        topology.local_face_count() + (mpi.rank() == target ? 1U : 0U);
    auto wrong_layout_state = hundun::flow::FlowState::create(
        registry, {decomposition.local_extent(), face_count}, fields,
        {0U, 0.0, 1.0e-3, 0.0,
         hundun::flow::MomentumTimeOrder::backward_euler});
    auto wrong_initial = initial;
    wrong_initial.face_velocity.resize(face_count * 3U, 0.0);
    wrong_initial.face_mass_flux.resize(face_count, 0.0);
    wrong_layout_state.seed_accepted_layers(wrong_initial, wrong_initial);
    bool layout_rejected = false;
    try {
      static_cast<void>(hundun::flow::IdealGasClosure::create(
          topology, geometry, boundaries, mpi, registry, fields,
          wrong_layout_state, {rho_h, cp, gas_constant, pressure}));
    } catch (const hundun::runtime::Error &) {
      layout_rejected = true;
    }
    HUNDUN_CHECK(layout_rejected);

    auto wrong_cell_extent = decomposition.local_extent();
    if (mpi.rank() == target)
      ++wrong_cell_extent.x;
    const auto wrong_cell_count =
        static_cast<std::size_t>(wrong_cell_extent.x) *
        static_cast<std::size_t>(wrong_cell_extent.y) *
        static_cast<std::size_t>(wrong_cell_extent.z);
    auto wrong_cell_state = hundun::flow::FlowState::create(
        registry, {wrong_cell_extent, topology.local_face_count()}, fields,
        {0U, 0.0, 1.0e-3, 0.0,
         hundun::flow::MomentumTimeOrder::backward_euler});
    auto wrong_cell_initial = initial;
    wrong_cell_initial.density.resize(wrong_cell_count, density);
    wrong_cell_initial.velocity.resize(wrong_cell_count * 3U, 0.0);
    wrong_cell_initial.mechanical_pressure.resize(wrong_cell_count, 0.0);
    wrong_cell_initial.transported_cell_fields[0].resize(
        wrong_cell_count, density * cp * temperature);
    wrong_cell_state.seed_accepted_layers(wrong_cell_initial,
                                          wrong_cell_initial);
    bool cell_layout_rejected = false;
    try {
      static_cast<void>(hundun::flow::IdealGasClosure::create(
          topology, geometry, boundaries, mpi, registry, fields,
          wrong_cell_state, {rho_h, cp, gas_constant, pressure}));
    } catch (const hundun::runtime::Error &) {
      cell_layout_rejected = true;
    }
    HUNDUN_CHECK(cell_layout_rejected);
    if (mpi.size() == 1)
      break;
  }

  using CreateFault = hundun::flow::test::IdealGasCreateFault;
  using TestAccess = hundun::flow::test::IdealGasClosureTestAccess;
  for (const int target : {0, mpi.size() - 1}) {
    for (const auto fault :
         {CreateFault::mode_disagreement, CreateFault::ownership_gap,
          CreateFault::ownership_overlap, CreateFault::ownership_swap,
          CreateFault::local_preparation}) {
      bool rejected = false;
      try {
        static_cast<void>(TestAccess::create(
            topology, geometry, boundaries, mpi, registry, fields, state,
            {rho_h, cp, gas_constant, pressure}, fault, target));
      } catch (const hundun::runtime::Error &) {
        rejected = true;
      }
      HUNDUN_CHECK(rejected);
      const auto successor = hundun::flow::IdealGasClosure::create(
          topology, geometry, boundaries, mpi, registry, fields, state,
          {rho_h, cp, gas_constant, pressure});
      HUNDUN_CHECK(successor.state().revision == 0U);
    }
    for (const auto fault :
         {CreateFault::sum_reduction, CreateFault::maximum_reduction}) {
      bool typed = false;
      try {
        static_cast<void>(TestAccess::create(
            topology, geometry, boundaries, mpi, registry, fields, state,
            {rho_h, cp, gas_constant, pressure}, fault, target));
      } catch (const hundun::runtime::MpiOperationError &error) {
        typed = std::string_view(error.what()).find(
                    fault == CreateFault::sum_reduction
                        ? "ideal-gas create sum reduction"
                        : "ideal-gas create maximum reduction") !=
                std::string_view::npos;
      }
      HUNDUN_CHECK(typed);
      const auto successor = hundun::flow::IdealGasClosure::create(
          topology, geometry, boundaries, mpi, registry, fields, state,
          {rho_h, cp, gas_constant, pressure});
      HUNDUN_CHECK(successor.state().revision == 0U);
    }
    if (mpi.size() == 1)
      break;
  }
  HUNDUN_CHECK(
      TestAccess::final_gate_rank_marker_is_collision_free(mpi.size()));

  expect_create_rejected(
      {fields.density, cp, gas_constant, pressure});
  expect_create_rejected(
      {std::numeric_limits<hundun::runtime::FieldId>::max(), cp,
       gas_constant, pressure});
  {
    for (const int target : {0, mpi.size() - 1}) {
      auto candidate_fields = fields;
      if (mpi.rank() == target)
        candidate_fields.transported_cell_fields = {fields.density};
      bool rejected = false;
      try {
        static_cast<void>(hundun::flow::IdealGasClosure::create(
            topology, geometry, boundaries, mpi, registry, candidate_fields,
            state, {rho_h, cp, gas_constant, pressure}));
      } catch (const hundun::runtime::Error &) {
        rejected = true;
      }
      HUNDUN_CHECK(rejected);
      if (mpi.size() == 1)
        break;
    }
  }
  {
    hundun::runtime::FieldRegistry wrong_registry;
    static_cast<void>(wrong_registry.declare_field(
        cell("rho", "kg/m3", 1U, true)));
    static_cast<void>(wrong_registry.declare_field(
        cell("velocity", "m/s", 3U, false)));
    static_cast<void>(wrong_registry.declare_field(
        cell("pi", "Pa", 1U, false)));
    static_cast<void>(wrong_registry.declare_field(
        face("face_velocity", "m/s", 3U)));
    static_cast<void>(
        hundun::finite_volume::declare_face_mass_flux(wrong_registry));
    static_cast<void>(wrong_registry.declare_field(
        cell("rho_h", "J/m3", 1U, true)));
    wrong_registry.freeze();
    for (const int target : {0, mpi.size() - 1}) {
      const auto &candidate_registry =
          mpi.rank() == target ? wrong_registry : registry;
      bool rejected = false;
      try {
        static_cast<void>(hundun::flow::IdealGasClosure::create(
            topology, geometry, boundaries, mpi, candidate_registry, fields,
            state, {rho_h, cp, gas_constant, pressure}));
      } catch (const hundun::runtime::Error &) {
        rejected = true;
      }
      HUNDUN_CHECK(rejected);
      if (mpi.size() == 1)
        break;
    }
  }
  {
    state.begin_attempt();
    expect_create_rejected({rho_h, cp, gas_constant, pressure});
    state.rollback_attempt();
  }
  const auto expect_state_rejected =
      [&](hundun::flow::FlowLayerValues history,
          hundun::flow::FlowLayerValues committed) {
        auto candidate = hundun::flow::FlowState::create(
            registry,
            {decomposition.local_extent(), topology.local_face_count()},
            fields,
            {0U, 0.0, 1.0e-3, 0.0,
             hundun::flow::MomentumTimeOrder::backward_euler});
        bool rejected = false;
        try {
          candidate.seed_accepted_layers(history, committed);
          static_cast<void>(hundun::flow::IdealGasClosure::create(
              topology, geometry, boundaries, mpi, registry, fields,
              candidate, {rho_h, cp, gas_constant, pressure}));
        } catch (const hundun::runtime::Error &) {
          rejected = true;
        }
        HUNDUN_CHECK(rejected);
      };
  for (const int target : {0, mpi.size() - 1}) {
    // FlowState itself rejects non-finite/non-positive seed values locally.
    // The closure matrix starts from finite states so all ranks reach the
    // same create preflight before the EOS/history contradiction is agreed.
    for (std::uint8_t mutation = 4U; mutation < 7U; ++mutation) {
      auto history = initial;
      auto committed = initial;
      if (mpi.rank() == target) {
        if (mutation == 4U)
          committed.transported_cell_fields[0][0] *= 1.01;
        else if (mutation == 5U)
          history.transported_cell_fields[0][0] *= 1.01;
        else
          history.density[0] *= 1.01;
      }
      expect_state_rejected(std::move(history), std::move(committed));
    }
    if (mpi.size() == 1)
      break;
  }
  const auto make_candidate_state = [&] {
    auto candidate = hundun::flow::FlowState::create(
        registry, {decomposition.local_extent(), topology.local_face_count()},
        fields,
        {0U, 0.0, 1.0e-3, 0.0,
         hundun::flow::MomentumTimeOrder::backward_euler});
    candidate.seed_accepted_layers(initial, initial);
    return candidate;
  };
  const auto commit_trial = [&](hundun::flow::FlowState &candidate,
                                std::uint64_t step) {
    candidate.commit_attempt(
        {step, static_cast<double>(step) * 1.0e-3, 1.0e-3,
         step == 1U ? 0.0 : 1.0e-3,
         step == 1U ? hundun::flow::MomentumTimeOrder::backward_euler
                    : hundun::flow::MomentumTimeOrder::bdf2});
  };
  const auto expect_layer_value_rejected =
      [&](int target, hundun::runtime::FieldId field, double value,
          bool history_layer) {
        auto candidate = make_candidate_state();
        candidate.begin_attempt();
        if (mpi.rank() == target) {
          auto changed = candidate.trial_layer().view<double>(field);
          changed(0, 0, 0, 0) = value;
        }
        commit_trial(candidate, 1U);
        if (history_layer) {
          candidate.begin_attempt();
          if (mpi.rank() == target) {
            auto restored = candidate.trial_layer().view<double>(field);
            restored(0, 0, 0, 0) =
                field == fields.density ? density : density * cp * temperature;
          }
          commit_trial(candidate, 2U);
        }
        bool rejected = false;
        try {
          static_cast<void>(hundun::flow::IdealGasClosure::create(
              topology, geometry, boundaries, mpi, registry, fields, candidate,
              {rho_h, cp, gas_constant, pressure}));
        } catch (const hundun::runtime::Error &) {
          rejected = true;
        }
        HUNDUN_CHECK(rejected);
        auto successor_state = make_candidate_state();
        const auto successor = hundun::flow::IdealGasClosure::create(
            topology, geometry, boundaries, mpi, registry, fields,
            successor_state, {rho_h, cp, gas_constant, pressure});
        HUNDUN_CHECK(successor.state().revision == 0U);
      };
  const auto expect_derived_value_rejected =
      [&](int target, double candidate_cp, double invalid_rho,
          double invalid_rho_h, bool history_layer) {
        auto candidate = make_candidate_state();
        const double base_rho_h = density * candidate_cp * temperature;
        const auto write_layer = [&](double target_rho,
                                     double target_rho_h) {
          candidate.begin_attempt();
          auto rho_values =
              candidate.trial_layer().view<double>(fields.density);
          auto rho_h_values = candidate.trial_layer().view<double>(rho_h);
          const auto extent = decomposition.local_extent();
          for (int k = 0; k < extent.z; ++k)
            for (int j = 0; j < extent.y; ++j)
              for (int i = 0; i < extent.x; ++i) {
                rho_values(i, j, k, 0) = density;
                rho_h_values(i, j, k, 0) = base_rho_h;
              }
          if (mpi.rank() == target) {
            rho_values(0, 0, 0, 0) = target_rho;
            rho_h_values(0, 0, 0, 0) = target_rho_h;
          }
        };
        write_layer(density, base_rho_h);
        commit_trial(candidate, 1U);
        write_layer(density, base_rho_h);
        commit_trial(candidate, 2U);
        if (history_layer) {
          write_layer(invalid_rho, invalid_rho_h);
          commit_trial(candidate, 3U);
          write_layer(density, base_rho_h);
          commit_trial(candidate, 4U);
        } else {
          write_layer(invalid_rho, invalid_rho_h);
          commit_trial(candidate, 3U);
        }

        const auto history_before =
            candidate.snapshot(hundun::flow::FlowLayer::history);
        const auto committed_before =
            candidate.snapshot(hundun::flow::FlowLayer::committed);
        const auto trial_before =
            candidate.snapshot(hundun::flow::FlowLayer::trial);
        const auto metadata_before = candidate.metadata();
        bool rejected = false;
        int failing_rank = -1;
        try {
          static_cast<void>(hundun::flow::IdealGasClosure::create(
              topology, geometry, boundaries, mpi, registry, fields, candidate,
              {rho_h, candidate_cp, gas_constant, pressure}));
        } catch (const hundun::runtime::Error &error) {
          rejected = true;
          failing_rank = TestAccess::preflight_failure_rank(error);
        }
        HUNDUN_CHECK(rejected && failing_rank == target);
        HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
            history_before,
            candidate.snapshot(hundun::flow::FlowLayer::history)));
        HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
            committed_before,
            candidate.snapshot(hundun::flow::FlowLayer::committed)));
        HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
            trial_before, candidate.snapshot(hundun::flow::FlowLayer::trial)));
        HUNDUN_CHECK(hundun::test::accepted_step_metadata_bitwise_equal(
            metadata_before, candidate.metadata()));

        auto successor_state = make_candidate_state();
        const auto successor = hundun::flow::IdealGasClosure::create(
            topology, geometry, boundaries, mpi, registry, fields,
            successor_state, {rho_h, cp, gas_constant, pressure});
        HUNDUN_CHECK(successor.state().revision == 0U);
      };
  for (const int target : {0, mpi.size() - 1}) {
    for (const bool history_layer : {false, true}) {
      for (const auto &[field, value] :
           {std::pair{fields.density, std::numeric_limits<double>::quiet_NaN()},
            std::pair{fields.density, 0.0},
            std::pair{rho_h, std::numeric_limits<double>::quiet_NaN()},
            std::pair{rho_h, 0.0},
            std::pair{rho_h, std::numeric_limits<double>::denorm_min()},
            std::pair{rho_h, std::numeric_limits<double>::max()}})
        expect_layer_value_rejected(target, field, value, history_layer);
      for (const auto &[candidate_cp, invalid_rho, invalid_rho_h] :
           {std::tuple{cp, std::numeric_limits<double>::denorm_min(),
                       std::numeric_limits<double>::max()},
            std::tuple{cp, std::numeric_limits<double>::max(),
                       std::numeric_limits<double>::denorm_min()},
            std::tuple{0.5, 1.0, std::numeric_limits<double>::max()},
            std::tuple{2.0, 1.0,
                       std::numeric_limits<double>::denorm_min()}})
        expect_derived_value_rejected(target, candidate_cp, invalid_rho,
                                      invalid_rho_h, history_layer);
    }

    {
      auto moved_from = make_candidate_state();
      auto moved_to = std::move(moved_from);
      const auto &candidate = mpi.rank() == target ? moved_from : moved_to;
      bool rejected = false;
      int failing_rank = -1;
      try {
        static_cast<void>(hundun::flow::IdealGasClosure::create(
            topology, geometry, boundaries, mpi, registry, fields, candidate,
            {rho_h, cp, gas_constant, pressure}));
      } catch (const hundun::runtime::Error &error) {
        rejected = true;
        failing_rank = TestAccess::preflight_failure_rank(error);
      }
      HUNDUN_CHECK(rejected && failing_rank == target);
      const auto successor = hundun::flow::IdealGasClosure::create(
          topology, geometry, boundaries, mpi, registry, fields, moved_to,
          {rho_h, cp, gas_constant, pressure});
      HUNDUN_CHECK(successor.state().revision == 0U);
    }
    {
      auto attempt_state = make_candidate_state();
      auto moved_to = std::move(attempt_state);
      auto attempt_closure = hundun::flow::IdealGasClosure::create(
          topology, geometry, boundaries, mpi, registry, fields, moved_to,
          {rho_h, cp, gas_constant, pressure});
      moved_to.begin_attempt();
      bool rejected = false;
      int failing_rank = -1;
      try {
        TestAccess::begin_attempt(
            attempt_closure, mpi.rank() == target ? attempt_state : moved_to,
            99U);
      } catch (const hundun::runtime::Error &error) {
        rejected = true;
        failing_rank = TestAccess::preflight_failure_rank(error);
      }
      HUNDUN_CHECK(rejected && failing_rank == target);
      moved_to.rollback_attempt();
      moved_to.begin_attempt();
      TestAccess::begin_attempt(attempt_closure, moved_to, 100U);
      TestAccess::rollback(attempt_closure);
      moved_to.rollback_attempt();
    }
    if (mpi.size() == 1)
      break;
  }

  const auto before = mpi.fp64_reduction_counters();
  auto closure = hundun::flow::IdealGasClosure::create(
      topology, geometry, boundaries, mpi, registry, fields, state,
      {rho_h, cp, gas_constant, pressure});
  const auto after = mpi.fp64_reduction_counters();
  HUNDUN_CHECK(after.collective_calls - before.collective_calls == 2U);
  HUNDUN_CHECK(after.reduced_scalars - before.reduced_scalars == 14U);
  HUNDUN_CHECK(after.logical_payload_bytes - before.logical_payload_bytes ==
               112U);
  const auto closure_state = closure.state();
  HUNDUN_CHECK(closure_state.mode ==
               hundun::flow::IdealGasPressureMode::closed_dynamic);
  HUNDUN_CHECK(closure_state.revision == 0U);
  HUNDUN_CHECK(closure_state.target_mass_kg.has_value());
  HUNDUN_CHECK_NEAR(closure_state.thermodynamic_pressure_pa, pressure,
                    1.0e-12 * pressure);
  HUNDUN_CHECK_NEAR(*closure_state.target_mass_kg, density, 5.0e-12 * density);
  HUNDUN_CHECK(hundun::flow::test::IdealGasClosureTestAccess::
                   post_store_rank_marker_is_collision_free(mpi.size()));

  bool inactive_rejected = false;
  try {
    hundun::flow::test::IdealGasClosureTestAccess::begin_attempt(closure, state,
                                                                 1U);
  } catch (const hundun::runtime::Error &) {
    inactive_rejected = true;
  }
  HUNDUN_CHECK(inactive_rejected);
  state.begin_attempt();
  bool zero_identity_rejected = false;
  try {
    hundun::flow::test::IdealGasClosureTestAccess::begin_attempt(closure, state,
                                                                 0U);
  } catch (const hundun::runtime::Error &) {
    zero_identity_rejected = true;
  }
  HUNDUN_CHECK(zero_identity_rejected);
  state.rollback_attempt();
  if (mpi.size() > 1)
    for (const int target : {0, mpi.size() - 1}) {
      state.begin_attempt();
      bool identity_rejected = false;
      try {
        hundun::flow::test::IdealGasClosureTestAccess::begin_attempt(
            closure, state, mpi.rank() == target ? 8U : 7U);
      } catch (const hundun::runtime::Error &) {
        identity_rejected = true;
      }
      HUNDUN_CHECK(identity_rejected);
      state.rollback_attempt();
    }
  state.begin_attempt();
  hundun::flow::test::IdealGasClosureTestAccess::begin_attempt(closure, state,
                                                               9U);
  bool overlap_rejected = false;
  try {
    hundun::flow::test::IdealGasClosureTestAccess::begin_attempt(closure, state,
                                                                 10U);
  } catch (const hundun::runtime::Error &) {
    overlap_rejected = true;
  }
  HUNDUN_CHECK(overlap_rejected);
  hundun::flow::test::IdealGasClosureTestAccess::rollback(closure);
  state.rollback_attempt();

  // Isolate the closure from the surrounding PISO arithmetic so the public
  // MpiContext counters prove each frozen P/C1/F prefix directly.
  state.begin_attempt();
  const auto success_before = mpi.fp64_reduction_counters();
  hundun::flow::test::IdealGasClosureTestAccess::begin_attempt(closure, state,
                                                               1U);
  const auto predictor =
      hundun::flow::test::IdealGasClosureTestAccess::evaluate(
          closure, state, hundun::flow::IdealGasClosureStage::predictor);
  require_fp64_delta(success_before, mpi.fp64_reduction_counters(), 3U, 19U,
                     152U);
  HUNDUN_CHECK(predictor.collective_count() == 5U);
  HUNDUN_CHECK(predictor.candidate_pressure_available());
  HUNDUN_CHECK_NEAR(predictor.candidate_pressure_pa(), pressure,
                    1.0e-12 * pressure);
  require_candidate_oracle(predictor);
  const auto provisional =
      hundun::flow::test::IdealGasClosureTestAccess::evaluate(
          closure, state, hundun::flow::IdealGasClosureStage::provisional);
  require_fp64_delta(success_before, mpi.fp64_reduction_counters(), 6U, 38U,
                     304U);
  HUNDUN_CHECK(provisional.collective_count() == 9U);
  HUNDUN_CHECK(provisional.candidate_pressure_available());
  HUNDUN_CHECK_NEAR(provisional.candidate_pressure_pa(), pressure,
                    1.0e-12 * pressure);
  require_candidate_oracle(provisional);
  const auto final = hundun::flow::test::IdealGasClosureTestAccess::evaluate(
      closure, state, hundun::flow::IdealGasClosureStage::final);
  require_fp64_delta(success_before, mpi.fp64_reduction_counters(), 10U, 59U,
                     472U);
  HUNDUN_CHECK(final.collective_count() == 14U);
  HUNDUN_CHECK(final.candidate_pressure_available());
  HUNDUN_CHECK_NEAR(final.candidate_pressure_pa(), pressure,
                    1.0e-12 * pressure);
  require_candidate_oracle(final);
  HUNDUN_CHECK(final.disposition() ==
               hundun::flow::IdealGasClosureDisposition::closed);
  hundun::flow::test::IdealGasClosureTestAccess::rollback(closure);
  state.rollback_attempt();

  const auto expect_derived_failure =
      [&](hundun::flow::IdealGasClosureStage failed_stage,
          std::uint64_t calls, std::uint64_t scalars, std::uint64_t bytes,
          std::uint64_t collectives, int target) {
        state.begin_attempt();
        hundun::flow::test::IdealGasClosureTestAccess::set_stage_failure(
            closure, failed_stage,
            hundun::flow::IdealGasClosureFailureReason::non_finite_temperature,
            target);
        const auto counters = mpi.fp64_reduction_counters();
        hundun::flow::test::IdealGasClosureTestAccess::begin_attempt(
            closure, state, 2U);
        hundun::flow::IdealGasClosureReport failed =
            hundun::flow::test::IdealGasClosureTestAccess::evaluate(
                closure, state,
                hundun::flow::IdealGasClosureStage::predictor);
        if (failed_stage != hundun::flow::IdealGasClosureStage::predictor)
          failed = hundun::flow::test::IdealGasClosureTestAccess::evaluate(
              closure, state,
              hundun::flow::IdealGasClosureStage::provisional);
        if (failed_stage == hundun::flow::IdealGasClosureStage::final)
          failed = hundun::flow::test::IdealGasClosureTestAccess::evaluate(
              closure, state, hundun::flow::IdealGasClosureStage::final);
        require_fp64_delta(counters, mpi.fp64_reduction_counters(), calls,
                           scalars, bytes);
        HUNDUN_CHECK(failed.collective_count() == collectives);
        HUNDUN_CHECK(failed.reason() ==
                     hundun::flow::IdealGasClosureFailureReason::
                         non_finite_temperature);
        HUNDUN_CHECK(failed.lowest_failing_rank() == target);
        if (failed_stage == hundun::flow::IdealGasClosureStage::predictor) {
          HUNDUN_CHECK(!failed.candidate_pressure_available());
          bool candidate_rejected = false;
          try {
            static_cast<void>(failed.candidate_pressure_pa());
          } catch (const hundun::runtime::Error &) {
            candidate_rejected = true;
          }
          HUNDUN_CHECK(candidate_rejected);
        } else {
          HUNDUN_CHECK(failed.candidate_pressure_available());
          HUNDUN_CHECK_NEAR(failed.candidate_pressure_pa(), pressure,
                            1.0e-12 * pressure);
        }
        require_candidate_oracle(failed);
        hundun::flow::test::IdealGasClosureTestAccess::rollback(closure);
        state.rollback_attempt();
      };
  for (const int target : {0, mpi.size() - 1}) {
    expect_derived_failure(hundun::flow::IdealGasClosureStage::predictor, 1U,
                           1U, 8U, 3U, target);
    expect_derived_failure(hundun::flow::IdealGasClosureStage::provisional,
                           4U, 20U, 160U, 7U, target);
    expect_derived_failure(hundun::flow::IdealGasClosureStage::final, 7U, 39U,
                           312U, 11U, target);
    if (mpi.size() == 1)
      break;
  }

  using MetricFault = hundun::flow::test::IdealGasMetricGateFault;
  const std::array metric_faults{MetricFault::eos, MetricFault::rho_remap,
                                 MetricFault::rho_h_remap, MetricFault::mass,
                                 MetricFault::enthalpy};
  const std::array metric_reasons{
      hundun::flow::IdealGasClosureFailureReason::eos_residual,
      hundun::flow::IdealGasClosureFailureReason::remap_residual,
      hundun::flow::IdealGasClosureFailureReason::remap_residual,
      hundun::flow::IdealGasClosureFailureReason::mass_conservation,
      hundun::flow::IdealGasClosureFailureReason::enthalpy_conservation};
  std::uint64_t metric_identity = 3U;
  for (const int target : {0, mpi.size() - 1}) {
    for (std::size_t fault_index = 0; fault_index < metric_faults.size();
         ++fault_index) {
      state.begin_attempt();
      const auto metric_before = mpi.fp64_reduction_counters();
      hundun::flow::test::IdealGasClosureTestAccess::begin_attempt(
          closure, state, metric_identity++);
      static_cast<void>(hundun::flow::test::IdealGasClosureTestAccess::evaluate(
          closure, state, hundun::flow::IdealGasClosureStage::predictor));
      static_cast<void>(hundun::flow::test::IdealGasClosureTestAccess::evaluate(
          closure, state, hundun::flow::IdealGasClosureStage::provisional));
      hundun::flow::test::IdealGasClosureTestAccess::set_metric_gate_failure(
          closure, metric_faults[fault_index], target);
      const auto metric_failed =
          hundun::flow::test::IdealGasClosureTestAccess::evaluate(
              closure, state, hundun::flow::IdealGasClosureStage::final);
      require_fp64_delta(metric_before, mpi.fp64_reduction_counters(), 9U, 57U,
                         456U);
      HUNDUN_CHECK(metric_failed.collective_count() == 13U);
      HUNDUN_CHECK(metric_failed.reason() == metric_reasons[fault_index]);
      HUNDUN_CHECK(metric_failed.lowest_failing_rank() == target);
      HUNDUN_CHECK(metric_failed.candidate_pressure_available());
      HUNDUN_CHECK_NEAR(metric_failed.candidate_pressure_pa(), pressure,
                        1.0e-12 * pressure);
      HUNDUN_CHECK(!metric_failed.final_metrics_available());
      require_candidate_oracle(metric_failed);
      hundun::flow::test::IdealGasClosureTestAccess::rollback(closure);
      state.rollback_attempt();
    }
    if (mpi.size() == 1)
      break;
  }

  for (const int target : {0, mpi.size() - 1}) {
    for (const bool enthalpy_density : {false, true}) {
      state.begin_attempt();
      const auto post_before = mpi.fp64_reduction_counters();
      hundun::flow::test::IdealGasClosureTestAccess::begin_attempt(
          closure, state, 4U);
      static_cast<void>(hundun::flow::test::IdealGasClosureTestAccess::evaluate(
          closure, state, hundun::flow::IdealGasClosureStage::predictor));
      static_cast<void>(hundun::flow::test::IdealGasClosureTestAccess::evaluate(
          closure, state, hundun::flow::IdealGasClosureStage::provisional));
      hundun::flow::test::IdealGasClosureTestAccess::set_post_store_corruption(
          closure, target, enthalpy_density);
      const auto post_failed =
          hundun::flow::test::IdealGasClosureTestAccess::evaluate(
              closure, state, hundun::flow::IdealGasClosureStage::final);
      require_fp64_delta(post_before, mpi.fp64_reduction_counters(), 10U, 59U,
                         472U);
      HUNDUN_CHECK(post_failed.collective_count() == 14U);
      HUNDUN_CHECK(post_failed.reason() ==
                   hundun::flow::IdealGasClosureFailureReason::eos_residual);
      HUNDUN_CHECK(post_failed.lowest_failing_rank() == target);
      HUNDUN_CHECK(post_failed.candidate_pressure_available());
      HUNDUN_CHECK_NEAR(post_failed.candidate_pressure_pa(), pressure,
                        1.0e-12 * pressure);
      require_candidate_oracle(post_failed);
      hundun::flow::test::IdealGasClosureTestAccess::rollback(closure);
      state.rollback_attempt();
    }
    if (mpi.size() == 1)
      break;
  }

  auto moved = std::move(closure);
  bool moved_from_rejected = false;
  try {
    static_cast<void>(closure.state());
  } catch (const hundun::runtime::Error &error) {
    moved_from_rejected = std::string_view(error.what()) ==
                          "ideal-gas closure has been moved from";
  }
  HUNDUN_CHECK(moved_from_rejected);
  HUNDUN_CHECK(moved.state().revision == 0U);

  auto open_decomposition = hundun::runtime::StructuredDecomposition::create(
      mpi, extent, {false, false, false},
      hundun::runtime::DecompositionOptions{grid(mpi.size())});
  hundun::mesh::MeshTopology open_topology(open_decomposition);
  hundun::mesh::MeshGeometry open_geometry(
      open_topology,
      hundun::mesh::UniformBoxMapping({0.0, 0.0, 0.0}, {1.0, 0.25, 0.25}));
  auto open_boundaries =
      hundun::boundary::BoundaryRegistry::create(open_case(), open_topology);
  auto open_state = hundun::flow::FlowState::create(
      registry,
      {open_decomposition.local_extent(), open_topology.local_face_count()},
      fields,
      {0U, 0.0, 6.25e-3, 0.0, hundun::flow::MomentumTimeOrder::backward_euler});
  auto open_initial = initial;
  open_initial.density.assign(open_topology.owned_cell_count(), density);
  open_initial.velocity.assign(open_topology.owned_cell_count() * 3U, 0.0);
  for (std::size_t cell_index = 0;
       cell_index < open_topology.owned_cell_count(); ++cell_index)
    open_initial.velocity[cell_index * 3U] = 1.0;
  open_initial.mechanical_pressure.assign(open_topology.owned_cell_count(),
                                          0.0);
  open_initial.face_velocity.assign(open_topology.local_face_count() * 3U, 0.0);
  open_initial.face_mass_flux.assign(open_topology.local_face_count(), 0.0);
  open_initial.transported_cell_fields = {std::vector<double>(
      open_topology.owned_cell_count(), density * cp * temperature)};
  open_state.seed_accepted_layers(open_initial, open_initial);
  const auto expect_open_spec_rejected =
      [&](hundun::flow::IdealGasClosureSpec spec) {
        bool rejected = false;
        try {
          static_cast<void>(hundun::flow::IdealGasClosure::create(
              open_topology, open_geometry, open_boundaries, mpi, registry,
              fields, open_state, spec));
        } catch (const hundun::runtime::Error &) {
          rejected = true;
        }
        HUNDUN_CHECK(rejected);
      };
  expect_open_spec_rejected({rho_h, cp * 1.01, gas_constant, pressure});
  expect_open_spec_rejected({rho_h, cp, gas_constant * 1.01, pressure});
  expect_open_spec_rejected({rho_h, cp, gas_constant, pressure * 1.01});
  for (const int target : {0, mpi.size() - 1}) {
    for (const auto fault :
         {CreateFault::inlet_cp, CreateFault::inlet_gas_constant,
          CreateFault::inlet_pressure}) {
      bool rejected = false;
      try {
        static_cast<void>(TestAccess::create(
            open_topology, open_geometry, open_boundaries, mpi, registry,
            fields, open_state, {rho_h, cp, gas_constant, pressure}, fault,
            target));
      } catch (const hundun::runtime::Error &) {
        rejected = true;
      }
      HUNDUN_CHECK(rejected);
      const auto successor = hundun::flow::IdealGasClosure::create(
          open_topology, open_geometry, open_boundaries, mpi, registry, fields,
          open_state, {rho_h, cp, gas_constant, pressure});
      HUNDUN_CHECK(successor.state().revision == 0U);
    }
    if (mpi.size() == 1)
      break;
  }
  const auto open_before = mpi.fp64_reduction_counters();
  auto open_closure = hundun::flow::IdealGasClosure::create(
      open_topology, open_geometry, open_boundaries, mpi, registry, fields,
      open_state, {rho_h, cp, gas_constant, pressure});
  const auto open_after = mpi.fp64_reduction_counters();
  HUNDUN_CHECK(open_after.collective_calls - open_before.collective_calls ==
               2U);
  HUNDUN_CHECK(open_after.reduced_scalars - open_before.reduced_scalars == 14U);
  const auto open_closure_state = open_closure.state();
  HUNDUN_CHECK(open_closure_state.mode ==
               hundun::flow::IdealGasPressureMode::open_fixed);
  HUNDUN_CHECK(!open_closure_state.target_mass_kg.has_value());
  HUNDUN_CHECK(open_closure_state.thermodynamic_pressure_pa == pressure);
  HUNDUN_CHECK(open_closure_state.revision == 0U);
  HUNDUN_CHECK(hundun::flow::test::IdealGasClosureTestAccess::
                   same_rank_reason_precedence_is_enum_order());
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
