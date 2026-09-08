// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "../../src/models_spray_transaction_detail.hpp"
#include "hundun/v04_field.hpp"

#include <mpi.h>
#include <array>
#include <iostream>

using namespace hundun::v04;
using namespace hundun::v04::spray;
using namespace hundun::v04::spray::detail;

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  int rank = 0, ranks = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &ranks);
  bool passed = true;
  auto expect = [&](bool condition, const char* message) {
    if (!condition) std::cerr << "rank " << rank << ": " << message << '\n';
    passed &= condition;
  };
  FieldRegistry registry;
  FieldId mass{};
  FieldSchema schema;
  StateLayers layers;
  ArenaLayout layout;
  Status status = registry.declare_field("mass", 1U, 0U, mass);
  if (status) status = registry.freeze(schema);
  const ArenaFieldRequest request{mass, {1, 1, 1}, {0U}, FieldLifetime::state_layer};
  if (status) status = ArenaLayout::compile(schema, {&request, 1U}, layout);
  if (status) status = StateLayers::allocate(layout, layers);
  expect(static_cast<bool>(status), "field fixture");
  AttemptTransaction gas;
  expect(static_cast<bool>(AttemptTransaction::create(1U, 0U, 1U, gas)), "field transaction");
  FieldView accepted;
  expect(static_cast<bool>(layers.view(StateRole::accepted_n, mass, accepted)), "accepted view");
  accepted.unchecked({0, 0, 0}, 0U) = 10.0;

  ParcelContainer parcels;
  expect(static_cast<bool>(parcels.reserve(8U)), "parcel capacity");
  SprayParcelState parcel;
  parcel.id = {static_cast<std::uint64_t>(rank + 1), 1U};
  parcel.droplet_mass_kg = 1.0;
  parcel.droplet_diameter_m = 0.01;
  parcel.multiplicity = 2.0;
  parcel.temperature_k = 300.0;
  parcel.liquid_material_fingerprint = 7U;
  expect(static_cast<bool>(parcels.begin_trial()) &&
         static_cast<bool>(parcels.stage_add(parcel)) &&
         static_cast<bool>(parcels.commit_trial()), "initial parcel");
  DeterministicInjector injector;
  InjectorSpec spec;
  spec.injector_id = rank + 1;
  spec.mass_flow_rate_kg_per_s = 0.25;
  spec.represented_mass_per_parcel_kg = 1.0;
  spec.droplet_mass_kg = 0.5;
  spec.droplet_diameter_m = 0.01;
  spec.temperature_k = 300.0;
  spec.liquid_material_fingerprint = 7U;
  expect(static_cast<bool>(injector.reserve(2U)) &&
         static_cast<bool>(injector.configure(spec)), "injector fixture");
  std::array<DeterministicInjector*, 1U> injectors{{&injector}};
  auto begin = [&]() {
    Status started = gas.begin(layers);
    if (started) started = gas.revise_trial(mass);
    FieldView trial;
    if (started) started = layers.view(StateRole::trial, mass, trial);
    if (started) trial.unchecked({0, 0, 0}, 0U) = 11.0;
    if (started) started = parcels.begin_trial();
    auto changed = parcel;
    changed.droplet_mass_kg = 0.5;
    changed.owner_global_cell = 5U;
    if (started) started = parcels.stage_update(changed);
    if (started) started = parcels.stage_tab_state(parcel.id, 0.2, 1.25);
    expect(static_cast<bool>(started), "begin isolated field and parcel trial");
    expect(injector.begin_trial(17U, 1.0).succeeded(), "injector residual trial");
  };
  auto unchanged = [&]() {
    FieldView current;
    SprayParcelState state;
    double deformation = 99.0, rate = 99.0;
    return static_cast<bool>(layers.view(StateRole::accepted_n, mass, current)) &&
        current.unchecked({0, 0, 0}, 0U) == 10.0 &&
        parcels.committed_at(0, state) && state.droplet_mass_kg == 1.0 &&
        state.owner_global_cell == 0U &&
        parcels.committed_tab_state(parcel.id, deformation, rate) &&
        deformation == 0.0 && rate == 0.0 &&
        injector.committed_state() == InjectorCommittedState{};
  };
  SprayAttemptFinishReport report;
  begin();
  const Status injected = rank == ranks - 1
      ? Status{StatusCode::numerical_failure, 777U} : Status{};
  expect(!finish_spray_attempt(MPI_COMM_WORLD, gas, parcels,
      {injectors.data(), injectors.size()}, injected, report), "one-rank gas failure rejects all");
  expect(!report.committed && report.lowest_failing_rank == ranks - 1 && unchanged(),
         "gas, parcel owner, TAB and injector residual all unchanged");
  expect(!gas.active() && !parcels.trial_active() && !injector.trial_active(),
         "rejected trial state is withdrawn");

  begin();
  if (rank == 0) expect(!parcels.stage_add(parcel), "duplicate parcel poisons local trial");
  expect(!finish_spray_attempt(MPI_COMM_WORLD, gas, parcels,
      {injectors.data(), injectors.size()}, {}, report) && unchanged(),
      "parcel preflight failure rejects an otherwise complete gas trial");

  begin();
  if (rank == ranks - 1) expect(static_cast<bool>(injector.rollback_trial()), "missing injector trial fixture");
  expect(!finish_spray_attempt(MPI_COMM_WORLD, gas, parcels,
      {injectors.data(), injectors.size()}, {}, report) && unchanged(),
      "missing injector candidate rejects complete parcel and gas trials");

  begin();
  std::array<DeterministicInjector*, 2U> repeated{{&injector, &injector}};
  expect(!finish_spray_attempt(MPI_COMM_WORLD, gas, parcels,
      {repeated.data(), repeated.size()}, {}, report) && unchanged(),
      "duplicate participant rejects before any publication");

  begin();
  expect(static_cast<bool>(finish_spray_attempt(MPI_COMM_WORLD, gas, parcels,
      {injectors.data(), injectors.size()}, {}, report)) && report.committed,
      "retry accepts every participant together");
  FieldView current;
  SprayParcelState result;
  double deformation = 0, rate = 0;
  expect(static_cast<bool>(layers.view(StateRole::accepted_n, mass, current)) &&
      current.unchecked({0, 0, 0}, 0U) == 11.0 && parcels.committed_at(0, result) &&
      result.droplet_mass_kg == 0.5 && result.owner_global_cell == 5U &&
      parcels.committed_tab_state(parcel.id, deformation, rate) &&
      deformation == 0.2 && rate == 1.25 &&
      injector.committed_state().residual_mass_kg == 0.25,
      "all committed authorities advance exactly once");
  int local = passed ? 1 : 0, global = 0;
  MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Finalize();
  return global ? 0 : 1;
}
