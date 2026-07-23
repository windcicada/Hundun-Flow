// SPDX-License-Identifier: Apache-2.0

#include "hundun/boundary/basic_boundary.hpp"
#include "hundun/execution/execution.hpp"
#include "hundun/finite_volume/cell_centered_fvm.hpp"
#include "hundun/flow/adaptive_time_control.hpp"
#include "hundun/linear/conjugate_gradient.hpp"
#include "hundun/linear/preconditioners.hpp"
#include "hundun/mesh/mesh_geometry.hpp"
#include "hundun/mesh/mesh_topology.hpp"
#include "hundun/runtime/exchange_plan.hpp"
#include "hundun/runtime/field_registry.hpp"
#include "hundun/runtime/halo_exchange.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/mpi_environment.hpp"
#include "hundun/runtime/structured_decomposition.hpp"
#include "tests/support/flow_state_equality.hpp"
#include "tests/support/test_main.hpp"

#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <string>

namespace {

hundun::runtime::Int3 grid(int ranks) {
  if (ranks == 1)
    return {1, 1, 1};
  if (ranks == 2)
    return {2, 1, 1};
  if (ranks == 4)
    return {2, 2, 1};
  throw hundun::runtime::Error("unsupported Task22 rank count");
}

hundun::runtime::FieldDescriptor cell(const char *name,
                                      std::uint32_t components) {
  return {name, "1", "task22", hundun::runtime::FunctionSpace::cell_average,
          hundun::runtime::ScalarType::float64, components, 2, false,
          hundun::runtime::RestartPolicy::persistent,
          hundun::runtime::OutputPolicy::never};
}
hundun::runtime::FieldDescriptor face(const char *name,
                                      std::uint32_t components) {
  return {name, "1", "task22", hundun::runtime::FunctionSpace::face_value,
          hundun::runtime::ScalarType::float64, components, 0, false,
          hundun::runtime::RestartPolicy::persistent,
          hundun::runtime::OutputPolicy::never};
}

hundun::config::FlowCaseConfig periodic_case() {
  hundun::config::FlowCaseConfig config{};
  config.schema_version = 2;
  config.simulation_type =
      hundun::config::SimulationType::variable_density_flow;
  config.density_model = hundun::config::DensityModel::constant;
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

void run(const hundun::runtime::MpiContext &mpi) {
  constexpr hundun::runtime::Int3 extent{8, 6, 4};
  auto decomposition = hundun::runtime::StructuredDecomposition::create(
      mpi, extent, {true, true, true},
      hundun::runtime::DecompositionOptions{grid(mpi.size())});
  const hundun::mesh::MeshTopology topology(decomposition);
  const hundun::mesh::MeshGeometry geometry(
      topology,
      hundun::mesh::UniformBoxMapping({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}));
  auto boundaries =
      hundun::boundary::BoundaryRegistry::create(periodic_case(), topology);

  hundun::runtime::FieldRegistry registry;
  hundun::flow::FlowFieldIds fields;
  fields.density = registry.declare_field(cell("rho", 1U));
  fields.velocity = registry.declare_field(cell("u", 3U));
  fields.mechanical_pressure = registry.declare_field(cell("pi", 1U));
  fields.face_velocity = registry.declare_field(face("uf", 3U));
  fields.face_mass_flux =
      hundun::finite_volume::declare_face_mass_flux(registry);
  fields.transported_cell_fields.push_back(
      registry.declare_field(cell("alpha", 1U)));
  registry.freeze();

  const auto box = decomposition.owned_box();
  const hundun::runtime::Int3 local{box.end.x - box.begin.x,
                                    box.end.y - box.begin.y,
                                    box.end.z - box.begin.z};
  auto state = hundun::flow::FlowState::create(
      registry, {local, topology.local_face_count()}, fields,
      {0U, 0.0, 0.01, 0.0,
       hundun::flow::MomentumTimeOrder::backward_euler});
  const auto cells = topology.owned_cell_count();
  hundun::flow::FlowLayerValues initial;
  initial.density.assign(cells, 1.0);
  initial.velocity.assign(cells * 3U, 0.0);
  initial.mechanical_pressure.assign(cells, 0.0);
  initial.face_velocity.assign(topology.local_face_count() * 3U, 0.0);
  initial.face_mass_flux.assign(topology.local_face_count(), 0.0);
  initial.transported_cell_fields = {std::vector<double>(cells, 1.0)};
  state.seed_accepted_layers(initial, initial);
  auto other_state = hundun::flow::FlowState::create(
      registry, {local, topology.local_face_count()}, fields,
      {0U, 0.0, 0.01, 0.0,
       hundun::flow::MomentumTimeOrder::backward_euler});
  other_state.seed_accepted_layers(initial, initial);

  auto halo = hundun::runtime::HaloExchange::create(
      decomposition,
      hundun::runtime::ExchangePlan::create(decomposition, local, 2));
  hundun::execution::CpuReferenceContext execution;
  hundun::linear::ConjugateGradientSolver momentum_solver(execution, mpi);
  hundun::linear::ConjugateGradientSolver pressure_solver(execution, mpi);
  hundun::linear::JacobiPreconditioner mx(execution), my(execution),
      mz(execution), pressure_preconditioner(execution);
  auto facade = hundun::flow::FixedStepConstantDensityFlow::create(
      decomposition, topology, geometry, boundaries, mpi, execution, halo,
      momentum_solver, {&mx, &my, &mz}, pressure_solver,
      pressure_preconditioner,
      {{fields.transported_cell_fields.front(),
        hundun::finite_volume::FiniteVolumeQuantity::scalar(0U), 0.0}});
  const hundun::config::FlowTimeConfig time{
      hundun::config::TimeMode::adaptive, 2, 0.01, 0.00125, 0.02,
      0.5, 0.25, 1.25, 0.5, 8};
  const int mismatch_rank = mpi.size() == 1 ? 0 : 1;
  const auto expect_factory_failure =
      [&](const hundun::config::FlowTimeConfig &candidate,
          hundun::config::DensityModel model, const char *prefix) {
        bool rejected = false;
        try {
          static_cast<void>(hundun::flow::Bdf2RetryController::create(
              candidate, model, topology, geometry, mpi, state));
        } catch (const hundun::runtime::Error &error) {
          rejected = true;
          HUNDUN_CHECK(std::string(error.what()).find(prefix) !=
                       std::string::npos);
          HUNDUN_CHECK(
              std::string(error.what()).find(std::to_string(mismatch_rank)) !=
              std::string::npos);
        }
        HUNDUN_CHECK(rejected);
        HUNDUN_CHECK(state.metadata().step == 0U);
      };
  {
    auto candidate = time;
    if (mpi.rank() == mismatch_rank)
      candidate.mode = static_cast<hundun::config::TimeMode>(255);
    expect_factory_failure(candidate, hundun::config::DensityModel::constant,
                           "time-control.create.invalid-config");
  }
  {
    auto candidate = time;
    if (mpi.rank() == mismatch_rank)
      candidate.steps = -1;
    expect_factory_failure(candidate, hundun::config::DensityModel::constant,
                           "time-control.create.invalid-config");
  }
  {
    auto candidate = time;
    if (mpi.rank() == mismatch_rank)
      candidate.initial_dt_s = std::numeric_limits<double>::infinity();
    expect_factory_failure(candidate, hundun::config::DensityModel::constant,
                           "time-control.create.invalid-config");
  }
  {
    auto candidate = time;
    if (mpi.rank() == mismatch_rank)
      candidate.cfl_target = 0.6;
    expect_factory_failure(candidate, hundun::config::DensityModel::constant,
                           "time-control.create.invalid-config");
  }
  {
    auto candidate = time;
    if (mpi.size() > 1) {
      if (mpi.rank() == mismatch_rank)
        ++candidate.steps;
      expect_factory_failure(candidate, hundun::config::DensityModel::constant,
                             "time-control.create.agreement");
    }
  }
  {
    auto model = hundun::config::DensityModel::constant;
    if (mpi.rank() == mismatch_rank)
      model = static_cast<hundun::config::DensityModel>(255);
    expect_factory_failure(time, model, "time-control.create.identity");
  }
  auto controller = hundun::flow::Bdf2RetryController::create(
      time, hundun::config::DensityModel::constant, topology, geometry, mpi,
      state);
  {
    auto wrong_model = hundun::flow::Bdf2RetryController::create(
        time, hundun::config::DensityModel::material, topology, geometry, mpi,
        state);
    const auto metadata_before = state.metadata();
    auto rejected = wrong_model.advance(state, facade, 1.0, 0.0, {}, {});
    HUNDUN_CHECK(rejected.disposition() ==
                 hundun::flow::TimeAdvanceDisposition::non_retryable_failure);
    HUNDUN_CHECK(rejected.attempt_count() == 0U);
    HUNDUN_CHECK(rejected.lowest_failing_rank() == 0);
    HUNDUN_CHECK(state.metadata().step == metadata_before.step);
  }
  {
    constexpr std::array<int, 4> control_fields{0, 1, 2, 3};
    if (mpi.size() > 1) {
      for (const int field : control_fields) {
        hundun::linear::SolveControl momentum{};
        if (mpi.rank() == mismatch_rank) {
          if (field == 0)
            momentum.atol *= 2.0;
          else if (field == 1)
            momentum.rtol *= 2.0;
          else if (field == 2)
            ++momentum.max_iterations;
          else
            ++momentum.residual_recompute_interval;
        }
        const auto rejected =
            controller.advance(state, facade, 1.0, 0.0, momentum, {});
        HUNDUN_CHECK(
            rejected.disposition() ==
            hundun::flow::TimeAdvanceDisposition::non_retryable_failure);
        HUNDUN_CHECK(rejected.attempt_count() == 0U);
        HUNDUN_CHECK(rejected.lowest_failing_rank() == mismatch_rank);
        HUNDUN_CHECK(state.metadata().step == 0U);
      }
    }
  }
  {
    auto rejected =
        controller.advance(other_state, facade, 1.0, 0.0, {}, {});
    HUNDUN_CHECK(rejected.disposition() ==
                 hundun::flow::TimeAdvanceDisposition::non_retryable_failure);
    HUNDUN_CHECK(rejected.attempt_count() == 0U);
    HUNDUN_CHECK(rejected.lowest_failing_rank() == 0);
    HUNDUN_CHECK(other_state.metadata().step == 0U);
    HUNDUN_CHECK(state.metadata().step == 0U);
  }
  HUNDUN_CHECK(controller.state().accepted_step == 0U);
  const auto before = state.snapshot(hundun::flow::FlowLayer::committed);
  auto report = controller.advance(state, facade, 1.0, 0.0, {}, {});
  HUNDUN_CHECK(report.disposition() ==
               hundun::flow::TimeAdvanceDisposition::committed);
  HUNDUN_CHECK(report.attempt_count() == 1U);
  HUNDUN_CHECK(report.attempt(0).order ==
               hundun::flow::MomentumTimeOrder::backward_euler);
  HUNDUN_CHECK(report.accepted_dt_s() == 0.01);
  HUNDUN_CHECK(report.stability_metrics_available());
  HUNDUN_CHECK(report.convective_rate_per_s() == 0.0);
  HUNDUN_CHECK(report.diffusive_rate_per_s() == 0.0);
  HUNDUN_CHECK(controller.state().accepted_step == 1U);
  HUNDUN_CHECK(state.metadata().step == 1U);
  HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
      before, state.snapshot(hundun::flow::FlowLayer::committed)));

  auto restored = hundun::flow::Bdf2RetryController::restore(
      time, hundun::config::DensityModel::constant, topology, geometry, mpi,
      state, controller.state());
  HUNDUN_CHECK(restored.state().state_seal == controller.state().state_seal);
  auto second = restored.advance(state, facade, 1.0, 0.0, {}, {});
  HUNDUN_CHECK(second.disposition() ==
               hundun::flow::TimeAdvanceDisposition::committed);
  HUNDUN_CHECK(second.attempt(0).order ==
               hundun::flow::MomentumTimeOrder::bdf2);
  HUNDUN_CHECK(state.metadata().step == 2U);

  auto other_controller = hundun::flow::Bdf2RetryController::restore(
      time, hundun::config::DensityModel::constant, topology, geometry, mpi,
      state, restored.state());
  bool cross_controller_rejected = false;
  try {
    static_cast<void>(other_controller.diagnostic_source(state, second));
  } catch (const hundun::runtime::Error &) {
    cross_controller_rejected = true;
  }
  HUNDUN_CHECK(cross_controller_rejected);

  const auto out_of_band_stencil =
      hundun::flow::make_momentum_time_stencil(
          hundun::flow::MomentumTimeOrder::bdf2, 0.01,
          state.metadata().dt_s);
  const auto out_of_band =
      facade.attempt(state, 1.0, 0.0, out_of_band_stencil, {}, {});
  HUNDUN_CHECK(out_of_band.disposition ==
               hundun::flow::StepAttemptDisposition::committed);
  const auto controller_before_rejection = restored.state();
  const auto rejected_after_out_of_band =
      restored.advance(state, facade, 1.0, 0.0, {}, {});
  HUNDUN_CHECK(rejected_after_out_of_band.disposition() ==
               hundun::flow::TimeAdvanceDisposition::non_retryable_failure);
  HUNDUN_CHECK(rejected_after_out_of_band.attempt_count() == 0U);
  HUNDUN_CHECK(restored.state().state_seal ==
               controller_before_rejection.state_seal);
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
  return hundun::test::run([&] { run(mpi); });
}
