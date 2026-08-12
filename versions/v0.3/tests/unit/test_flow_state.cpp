// SPDX-License-Identifier: Apache-2.0

#include "hundun/fvm_cell_centered.hpp"
#include "hundun/flow_state.hpp"
#include "hundun/rt_field_access_plan.hpp"
#include "hundun/rt_field_registry.hpp"
#include "tests/support/allocation_attempt_guard.hpp"
#include "tests/support/test_main.hpp"

#include <type_traits>

namespace {

hundun::runtime::FieldDescriptor cell_field(const char *name,
                                            std::uint32_t components) {
  return {name,
          "1",
          "task18",
          hundun::runtime::FunctionSpace::cell_average,
          hundun::runtime::ScalarType::float64,
          components,
          2,
          false,
          hundun::runtime::RestartPolicy::persistent,
          hundun::runtime::OutputPolicy::never};
}

hundun::runtime::FieldDescriptor face_field(const char *name,
                                            std::uint32_t components) {
  return {name,
          "1",
          "task18",
          hundun::runtime::FunctionSpace::face_value,
          hundun::runtime::ScalarType::float64,
          components,
          0,
          false,
          hundun::runtime::RestartPolicy::persistent,
          hundun::runtime::OutputPolicy::never};
}

void run_flow_state_red_contract() {
  static_assert(!std::is_copy_constructible_v<hundun::flow::FlowState>);
  static_assert(std::is_move_constructible_v<hundun::flow::FlowState>);

  hundun::runtime::FieldRegistry registry;
  hundun::flow::FlowFieldIds ids;
  ids.density = registry.declare_field(cell_field("rho", 1U));
  ids.velocity = registry.declare_field(cell_field("velocity", 3U));
  ids.mechanical_pressure = registry.declare_field(cell_field("pi", 1U));
  ids.face_velocity = registry.declare_field(face_field("u_face", 3U));
  ids.face_mass_flux = hundun::finite_volume::declare_face_mass_flux(registry);
  registry.freeze();

  auto state = hundun::flow::FlowState::create(
      registry, {{2, 2, 1}, 12U}, ids,
      {4U, 0.4, 0.1, 0.1, hundun::flow::MomentumTimeOrder::backward_euler});
  HUNDUN_CHECK(state.metadata().step == 4U);
  HUNDUN_CHECK(
      state.layer(hundun::flow::FlowLayer::committed).layout_set().face_count ==
      12U);

  hundun::flow::FlowLayerValues history;
  history.density.assign(4U, 2.0);
  history.velocity.assign(12U, -1.0);
  history.mechanical_pressure.assign(4U, 3.0);
  history.face_velocity.assign(36U, -2.0);
  history.face_mass_flux.assign(12U, -4.0);
  auto committed = history;
  committed.velocity.assign(12U, 5.0);
  committed.mechanical_pressure.assign(4U, 6.0);
  committed.face_mass_flux.assign(12U, 7.0);
  state.seed_accepted_layers(history, committed);

  hundun::runtime::FieldAccessPlan stale_plan(registry);
  stale_plan.declare_access(1U, 1U, ids.velocity,
                            hundun::runtime::AccessMode::read_write);
  stale_plan.freeze();
  auto stale = state.trial_layer().acquire_write<double>(stale_plan, 1U, 1U,
                                                         ids.velocity);
  state.begin_attempt();
  bool rejected = false;
  try {
    static_cast<void>(stale(0, 0, 0, 0));
  } catch (const hundun::runtime::Error &) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);

  auto trial_values = state.snapshot(hundun::flow::FlowLayer::trial);
  HUNDUN_CHECK(trial_values.velocity == committed.velocity);
  trial_values.velocity.assign(12U, 9.0);
  trial_values.mechanical_pressure.assign(4U, 10.0);
  state.rollback_attempt();
  HUNDUN_CHECK(state.snapshot(hundun::flow::FlowLayer::committed).velocity ==
               committed.velocity);
  HUNDUN_CHECK(state.metadata().step == 4U);

  std::size_t warmed_copy_allocation_attempts = 0U;
  {
    hundun::test::allocation_probe::AllocationAttemptGuard guard;
    state.begin_attempt();
    state.rollback_attempt();
    warmed_copy_allocation_attempts = guard.attempts();
  }
  HUNDUN_CHECK(warmed_copy_allocation_attempts == 0U);

  state.begin_attempt();
  auto trial_velocity = state.trial_layer().acquire_write<double>(
      stale_plan, 1U, 1U, ids.velocity);
  for (int j = 0; j < 2; ++j) {
    for (int i = 0; i < 2; ++i) {
      for (int component = 0; component < 3; ++component) {
        trial_velocity(i, j, 0, component) = 9.0;
      }
    }
  }
  state.commit_attempt(
      {5U, 0.5, 0.1, 0.1, hundun::flow::MomentumTimeOrder::backward_euler});
  HUNDUN_CHECK(state.snapshot(hundun::flow::FlowLayer::history).velocity ==
               committed.velocity);
  HUNDUN_CHECK(state.snapshot(hundun::flow::FlowLayer::committed).velocity ==
               std::vector<double>(12U, 9.0));
  HUNDUN_CHECK(state.metadata().step == 5U);
  rejected = false;
  try {
    static_cast<void>(trial_velocity(0, 0, 0, 0));
  } catch (const hundun::runtime::Error &) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);
}

} // namespace

int main() {
  return hundun::test::run([] { run_flow_state_red_contract(); });
}
