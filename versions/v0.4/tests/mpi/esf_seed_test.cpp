// SPDX-License-Identifier: Apache-2.0
#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "hundun/v04_app.hpp"
using namespace hundun::v04;
namespace {
bool all(bool value) {
  int result = value;
  MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  return result != 0;
}
bool run(const char* path, int rank) {
  ValidatedModel model;
  auto status = CaseCompiler::load_and_compile(MPI_COMM_WORLD, path, model);
  if (!all(bool(status))) {
    if (rank == 0)
      std::cerr << "seed_case status=" << unsigned(status.code) << '/'
                << status.detail << '\n';
    return false;
  }
  model.time.scheme = TimeScheme::cn_be;
  model.solver.coupling = CouplingKind::outer_corrected;
  model.pressure_reference = PressureReferenceKind::boundary_absolute;
  model.legacy_time_fingerprint = model.fingerprint;
  model.fingerprint ^= 0x45534653454544ULL;
  model.boundaries[2].flow_kind = BoundaryKind::symmetry;
  model.boundaries[2].scalars = {{"A", ScalarBoundaryKind::zero_gradient}};
  auto& outlet = model.boundaries[3];
  outlet.flow_kind = BoundaryKind::pressure_outlet;
  outlet.pressure = 101325;
  outlet.allow_backflow = true;
  outlet.backflow_temperature = 300;
  outlet.scalars = {{"A", ScalarBoundaryKind::zero_gradient}};
  outlet.scalars[0].backflow_kind = ScalarBoundaryKind::dirichlet;
  outlet.scalars[0].backflow_value = .25;
  auto create = [&](ProductDriver& d) {
    CompiledCasePlan plan;
    auto s = ProductCompiler::compile(MPI_COMM_WORLD, model, path, plan);
    if (s) s = ProductDriver::create(MPI_COMM_WORLD, std::move(plan), d);
    return s;
  };
  ProductDriver source;
  status = create(source);
  DriverInitialState initial;
  initial.pressure_reference = 101325;
  initial.temperature = 300;
  double y = .25;
  initial.transported_scalars = {&y, 1};
  if (status) status = source.initialize(initial);
  DriverStepReport step;
  if (status) status = source.advance({1, 1, 1, 1, 1}, step);
  RestartSnapshot state;
  if (status) status = source.committed_restart_snapshot(state);
  if (!all(bool(status))) {
    if (rank == 0)
      std::cerr << "seed_initial status=" << unsigned(status.code) << '/'
                << status.detail << '\n';
    return false;
  }
  RestartImage image;
  image.global_cells = state.global_cells;
  image.patch = state.patch;
  image.plan = state.plan;
  image.schema = state.schema;
  image.geometry = state.geometry;
  image.time = .004;
  image.dt = .001;
  image.step = 4;
  image.controller_state = 1;
  image.pressure_reference = 101325;
  image.backward_euler_recovery = true;
  image.final_mass_flux_revision = 17;
  std::vector<double> expected_cache;
  bool auxiliary_seen=false;
  for (std::size_t i = 0; i < state.fields.size; ++i) {
    const auto& f = state.fields.data[i];
    const auto v = f.values;
    RestartImageField copied;
    copied.role = f.role;
    copied.field = v.field;
    copied.components = v.components;
    for (int z = 0; z < v.interior.z; ++z)
      for (int j = 0; j < v.interior.y; ++j)
        for (int x = 0; x < v.interior.x; ++x)
          for (unsigned c = 0; c < v.components; ++c)
            copied.values.push_back(v.unchecked({x, j, z}, c));
    if (f.role == RestartFieldRole::stochastic_transport) {
      expected_cache = copied.values;
      for (double& value : copied.values) value *= 3;
    }
    if (f.role == RestartFieldRole::stochastic_auxiliary) {
      // Distinct signed auxiliary state survives both exact and explicit
      // method recovery; physical PDF values remain the source authority.
      if (copied.components != 3) return false;
      auxiliary_seen=true;
      for (std::size_t c=0;c<copied.values.size();c+=3) {
        copied.values[c]=-.0001;
        copied.values[c+1]=1.0001;
      }
    }
    image.fields.push_back(std::move(copied));
  }
  unsigned axis = 0;
  for (const auto v : {state.final_mass_flux.x, state.final_mass_flux.y,
                       state.final_mass_flux.z}) {
    for (int z = 0; z < v.extents.z; ++z)
      for (int j = 0; j < v.extents.y; ++j)
        for (int x = 0; x < v.extents.x; ++x)
          image.final_mass_flux[axis].push_back(v.unchecked({x, j, z}));
    ++axis;
  }
  if (!all(!expected_cache.empty() && auxiliary_seen)) return false;
  std::vector<double> rebuilt_cache(expected_cache.size());
  const RestartImageField *enthalpy = nullptr, *species = nullptr;
  for (const auto& f : image.fields) {
    if (f.role == RestartFieldRole::enthalpy) enthalpy = &f;
    if (f.role == RestartFieldRole::independent_species) species = &f;
  }
  if (!all(enthalpy && species)) return false;
  for (std::size_t i = 0; i < enthalpy->values.size(); ++i) {
    // Fixture: common cp=1000 and formation-enthalpy difference=100000.
    const double t =
        298.15 + (enthalpy->values[i] - 1e5 * species->values[i]) / 1000.;
    const double mu = 1e-5 * std::sqrt(t / 300.);
    rebuilt_cache[4 * i] = mu / .7;
    rebuilt_cache[4 * i + 2] = mu;
  }
  double max_error = 0;
  for (bool rebuild : {false, true}) {
    ProductDriver restored;
    status = create(restored);
    if (status)
      status = restored.initialize_restart(
          image, RestartStorageCompatibility::strict,
          rebuild ? RestartHistoryPolicy::rebuild_method_history
                  : RestartHistoryPolicy::require_compatible);
    RestartSnapshot loaded;
    if (status) status = restored.committed_restart_snapshot(loaded);
    if (!all(bool(status))) {
      if (rank == 0)
        std::cerr << "seed_restore status=" << unsigned(status.code) << '/'
                  << status.detail << '\n';
      return false;
    }
    double error = 0;
    bool unchanged = loaded.fields.size == image.fields.size();
    for (std::size_t i = 0; i < loaded.fields.size; ++i) {
      const auto& f = loaded.fields.data[i];
      const auto v = f.values;
      std::size_t index = 0;
      for (int z = 0; z < v.interior.z; ++z)
        for (int j = 0; j < v.interior.y; ++j)
          for (int x = 0; x < v.interior.x; ++x)
            for (unsigned c = 0; c < v.components; ++c) {
              if (f.role != RestartFieldRole::stochastic_transport) {
                unchanged =
                    unchanged && i < image.fields.size() &&
                    index < image.fields[i].values.size() &&
                    v.unchecked({x, j, z}, c) == image.fields[i].values[index];
                ++index;
                continue;
              }
              const double target =
                  rebuild ? rebuilt_cache[index] : expected_cache[index] * 3;
              ++index;
              error =
                  std::max(error, std::abs(v.unchecked({x, j, z}, c) - target) /
                                      std::max(1e-30, std::abs(target)));
            }
    }
    MPI_Allreduce(MPI_IN_PLACE, &error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    if (!all(unchanged)) return false;
    max_error = std::max(error, max_error);
    if (rank == 0)
      std::cout << "ESF_seed rebuild=" << rebuild
                << " cache_relative_error=" << error << '\n';
  }
  return all(max_error < 1e-12);
}
}  // namespace
int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  int rank{};
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  const bool okay = argc == 2 && run(argv[1], rank);
  MPI_Finalize();
  return okay ? 0 : 1;
}
