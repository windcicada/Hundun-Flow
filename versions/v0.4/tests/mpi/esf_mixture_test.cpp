// SPDX-License-Identifier: Apache-2.0
#include "../support/product_fixture.hpp"
#include "../support/dense_solve.hpp"
#include "hundun/v04_app.hpp"
#include "hundun/v04_portable.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <vector>

using namespace hundun::v04;
namespace {
// Three inert isomers isolate transport from chemical and thermal changes.
class Gas final : public portable::GasQueryProvider,
                  public portable::GasAdvanceProvider {
public:
  portable::GasIdentity identity;
  combustion::ChemistryIdentity closure;
  double transport_error{}, expected_transport{}, expected_low{}, expected_high{};
  unsigned transport_samples{};
  static constexpr double mw = 28.96546;
  static constexpr double cp = 3.5 * kUniversalGasConstant / mw;
  Gas() {
    identity.mechanism_sha256 = std::string(64, 'a');
    identity.phase = "inert-mixture";
    identity.species_names = {"A", "B", "C"};
    identity.element_names = {"X"};
    identity.element_counts = {1, 1, 1};
    identity.molecular_weights_kg_per_kmol = {mw, mw, mw};
    identity.enthalpy_reference = "nasa7-total-enthalpy";
    identity.composition_fingerprint = 3101;
    closure.element_count = 1;
    closure.species = {{mw, {1}, 0}, {mw, {1}, 0}, {mw, {1}, 0}};
    closure.fingerprint = combustion::chemistry_identity_fingerprint(closure);
    identity.closure_fingerprint = closure.fingerprint;
  }
  const portable::GasIdentity &gas_identity() const noexcept override {
    return identity;
  }
  portable::Status query_gas(const portable::GasQuery &q,
                             portable::GasQueryOutput &out) noexcept override {
    if (q.species_count != 3 || out.capacity != 3)
      return portable::Status::invalid_input;
    const double temperature =
        q.coordinates == portable::GasStateCoordinates::pressure_enthalpy
            ? q.enthalpy_j_per_kg / cp : q.temperature_k;
    out.sample = {q.revision, q.composition_fingerprint, q.pressure_pa,
                  temperature, q.pressure_pa * mw /
                      (kUniversalGasConstant * temperature),
                  cp * temperature, cp, 1.8e-5, .026};
    for (unsigned s = 0; s < 3; ++s) {
      out.diffusivities_m2_per_s[s] = 1e-5;
      out.species_enthalpies_j_per_kg[s] = cp * temperature;
      out.net_mass_rates_kg_per_m3_s[s] = 0;
    }
    return portable::Status::success;
  }
  portable::Status advance_gas(const portable::GasAdvanceQuery &q,
                               portable::GasAdvanceOutput &out) noexcept override {
    double d[3], h[3], w[3];
    portable::GasQueryOutput sample{{}, d, h, w, 3};
    const auto status = query_gas(q.state, sample);
    if (status != portable::Status::success) return status;
    if (q.start_time_s == 0 && q.state.mass_fractions[0] > .15 &&
        q.state.mass_fractions[0] < .25) {
      ++transport_samples;
      transport_error = std::max(transport_error,
          std::min(std::abs(q.state.mass_fractions[0]-expected_low),
                   std::abs(q.state.mass_fractions[0]-expected_high)));
    }
    for (unsigned s = 0; s < 3; ++s) {
      out.final_mass_fractions[s] = q.state.mass_fractions[s];
      out.integrated_species_density_delta_kg_per_m3[s] = 0;
    }
    out.final_sample = sample.sample;
    out.completed_duration_s = q.duration_s;
    out.internal_step_count = 1;
    out.integrated_heat_release_j_per_m3 = 0;
    return portable::Status::success;
  }
};

std::array<double, 3> composition(int x, bool zero) {
  auto y = x < 4 ? std::array<double, 3>{0, .6, .4}
       : x == 4 ? std::array<double, 3>{.2, .7, .1}
                : std::array<double, 3>{.5, .5, 0};
  if (!zero) for (auto &value : y) value = .001 + .997 * value;
  return y;
}

bool run(unsigned fields, bool zero) {
  const Int3 cells{12, 8, 8};
  constexpr double dt = 1e-3, pressure = 1e5, temperature = 300;
  const double density = pressure * Gas::mw / (kUniversalGasConstant * temperature);
  auto model = test::product_model(cells);
  model.turbulence = TurbulenceKind::none;
  model.time.control = TimeControlKind::fixed;
  model.time.initial_dt = model.time.minimum_dt = model.time.maximum_dt = dt;
  model.time.maximum_retries = 1;
  model.pressure_reference = PressureReferenceKind::boundary_absolute;
  auto &inlet = model.boundaries[0];
  inlet.flow_kind = BoundaryKind::velocity_inlet;
  inlet.velocity = {1, 0, 0}; inlet.direction = {1, 0, 0};
  inlet.temperature = temperature;
  const auto inlet_y = composition(0, zero);
  inlet.scalars = {{"A", ScalarBoundaryKind::dirichlet, inlet_y[0]},
                   {"B", ScalarBoundaryKind::dirichlet, inlet_y[1]}};
  auto &outlet = model.boundaries[1];
  outlet.flow_kind = BoundaryKind::pressure_outlet;
  outlet.pressure = pressure; outlet.backflow_temperature = temperature;
  outlet.scalars = {{"A", ScalarBoundaryKind::zero_gradient},
                    {"B", ScalarBoundaryKind::zero_gradient}};
  const auto species = model.thermophysics.species.front();
  model.thermophysics.species.assign(3, species);
  for (unsigned s = 0; s < 3; ++s)
    model.thermophysics.species[s].stable_name = std::string(1, char('A' + s));
  model.transported_scalars = {{"A", TransportedScalarRole::species, .7, .7},
                               {"B", TransportedScalarRole::species, .7, .7}};
  Gas gas;
  // The accepted mean has its own deterministic conservation/reconciliation
  // oracle; chemistry receives the joint implicit statistical transport.
  const double transported = .2 - dt * .2 / (2.0 / cells.x);
  gas.expected_transport = zero ? transported : .001 + .997 * transported;
  const double volume=1./(cells.x*cells.y*cells.z);
  const double mix=dt*1.8e-5/(density*std::pow(volume,2./3.));
  const double dx=2./cells.x,courant=dt/dx;
  const double inlet_diffusion=2*dt*(.026/Gas::cp)/(density*dx*dx);
  // Constant h contributes VLS ratio=0: all interior total conductances
  // equal rho*U*A/2. The resulting one-dimensional matrix is upwind BE.
  const auto matrix=[&](double mixing) {
    std::vector<std::vector<double>> a(cells.x,std::vector<double>(cells.x));
    for(int x=0;x<cells.x;++x) {
      a[x][x]=1+courant+mixing;
      if(x)a[x][x-1]=-courant;
    }
    a[0][0]+=inlet_diffusion;
    return a;
  };
  std::vector<double> old(cells.x),rhs(cells.x);
  for(int x=0;x<cells.x;++x)old[x]=composition(x,zero)[0];
  rhs=old;rhs[0]+=(courant+inlet_diffusion)*inlet_y[0];
  const auto mean=test::dense_solve(matrix(0),rhs);
  std::array<double,3> lower{1,1,1},upper{};
  for(int x=0;x<cells.x;++x)for(unsigned s=0;s<3;++s) {
    lower[s]=std::min(lower[s],composition(x,zero)[s]);
    upper[s]=std::max(upper[s],composition(x,zero)[s]);
  }
  std::array<double,2> endpoint{};
  for(unsigned sign=0;sign<2;++sign) {
    rhs=old;
    for(int x=0;x<cells.x;++x) {
      rhs[x]+=mix*mean[x];
      if(x==0 || x==cells.x-1)continue;
      const auto y=composition(x,zero);
      std::array<double,3> noise{};double alpha=1;
      for(unsigned s=0;s<3;++s) {
        noise[s]=(sign?1:-1)*std::sqrt(2*(1.8e-5/.7)*dt/density)*
            (composition(x+1,zero)[s]-composition(x-1,zero)[s])/(2*dx);
        if(noise[s]>0)alpha=std::min(alpha,(upper[s]-y[s])/noise[s]);
        if(noise[s]<0)alpha=std::min(alpha,(lower[s]-y[s])/noise[s]);
      }
      rhs[x]+=alpha*noise[0];
    }
    rhs[0]+=(courant+inlet_diffusion)*inlet_y[0];
    endpoint[sign]=test::dense_solve(matrix(mix),rhs)[4];
  }
  gas.expected_low=std::min(endpoint[0],endpoint[1]);
  gas.expected_high=std::max(endpoint[0],endpoint[1]);
  model.reaction.mode = ReactionMode::esf_tpdf;
  model.reaction.mechanism_sha256 = gas.identity.mechanism_sha256;
  model.reaction.phase = gas.identity.phase;
  model.reaction.esf = EsfSpec{};
  model.reaction.esf->fields = fields;
  CompiledCasePlan plan;
  auto status = ProductCompiler::compile(MPI_COMM_WORLD, model, {}, plan,
                                         {&gas, &gas, &gas.closure});
  ProductDriver driver;
  if (status) status = ProductDriver::create(MPI_COMM_WORLD, std::move(plan), driver);
  RestartExpected expected;
  if (status) status = driver.restart_expected(expected);
  RestartImage start;
  start.global_cells = expected.global_cells; start.patch = expected.target_patch;
  start.plan = expected.plan; start.schema = expected.schema; start.geometry = expected.geometry;
  start.dt = dt; start.step = 4; start.controller_state = 1;
  start.pressure_reference = pressure; start.backward_euler_recovery = true;
  const auto n = start.patch.cells;
  unsigned independent = 0;
  for (std::size_t f = 0; f < expected.fields.size && status; ++f) {
    const auto &spec = expected.fields.data[f];
    RestartImageField field;
    field.role = spec.role; field.field = spec.field; field.components = spec.components;
    field.values.resize(std::size_t(n.x) * n.y * n.z * spec.components);
    for (int z = 0; z < n.z; ++z) for (int y = 0; y < n.y; ++y) for (int x = 0; x < n.x; ++x) {
      const auto i = std::size_t((z * n.y + y) * n.x + x) * spec.components;
      const auto Y = composition(x + start.patch.begin.x, zero);
      if (spec.role == RestartFieldRole::velocity) field.values[i] = 1;
      if (spec.role == RestartFieldRole::pressure_absolute) field.values[i] = pressure;
      if (spec.role == RestartFieldRole::enthalpy) field.values[i] = Gas::cp * temperature;
      if (spec.role == RestartFieldRole::independent_species) field.values[i] = Y[independent];
      if (spec.role == RestartFieldRole::stochastic_field ||
          spec.role == RestartFieldRole::stochastic_auxiliary) {
        std::copy(Y.begin(), Y.end(), field.values.begin() + i);
        field.values[i + 3] = Gas::cp * temperature;
      }
      if (spec.role == RestartFieldRole::stochastic_transport) {
        field.values[i] = .026 / Gas::cp;
        field.values[i + 2] = 1.8e-5;
      }
    }
    if (spec.role == RestartFieldRole::independent_species) ++independent;
    start.fields.push_back(std::move(field));
  }
  for (unsigned a = 0; a < 3; ++a) {
    auto extent = n; (a == 0 ? extent.x : a == 1 ? extent.y : extent.z)++;
    start.final_mass_flux[a].assign(std::size_t(extent.x) * extent.y * extent.z,
                                  a == 0 ? density / 128 : 0);
  }
  if (status) status = driver.initialize_restart(start);
  if (!status) std::cerr << "ESF mixture initialize " << unsigned(status.code) << '/' << status.detail << '\n';
  DriverStepReport report;
  if (status) status = driver.advance({1, 1, 1, 1, 1}, report);
  if (!status) {
    const auto &failure = report.thermophysical_predictor.failure;
    std::cerr << "ESF mean predictor stage=" << report.failed_stage
        << " reason=" << unsigned(failure.reason) << " observed=" << failure.observed_value
        << " lower=" << failure.allowed_lower << " upper=" << failure.allowed_upper << '\n';
  }
  RestartSnapshot result;
  if (status) status = driver.committed_restart_snapshot(result);
  bool passed = status && report.accepted && result.step == 5;
  double error[4]{};
  std::vector<double> mean_a(std::size_t(n.x)*n.y*n.z,0);
  if (status) for (std::size_t f = 0; f < result.fields.size; ++f) {
    const auto &field = result.fields.data[f];
    if (field.role != RestartFieldRole::stochastic_field) continue;
    for (int z = 0; z < n.z; ++z) for (int y = 0; y < n.y; ++y) for (int x = 0; x < n.x; ++x) {
      double sum = 0;
      for (unsigned s = 0; s < 3; ++s) {
        const double value = field.values.unchecked({x, y, z}, s);
        passed &= value >= 0 && value <= 1;
        sum += value;
      }
      error[0] = std::max(error[0], std::abs(sum - 1));
      error[1] = std::max(error[1], std::abs(field.values.unchecked({x, y, z}, 3) / Gas::cp - temperature));
      mean_a[std::size_t(x)+std::size_t(n.x)*(y+std::size_t(n.y)*z)]+=
          field.values.unchecked({x,y,z},0)/fields;
    }
  }
  for(int z=0;z<n.z;++z) for(int y=0;y<n.y;++y) for(int x=0;x<n.x;++x)
    if(x+start.patch.begin.x==4)
      error[3]=std::max(error[3],std::abs(mean_a[std::size_t(x)+std::size_t(n.x)*
          (y+std::size_t(n.y)*z)]-gas.expected_transport));
  // Inspect the public chemistry request at the transport/chemistry boundary.
  // The later pressure/mean reconciliation has its own accepted-state checks.
  error[2] = gas.transport_error;
  unsigned samples = gas.transport_samples;
  MPI_Allreduce(MPI_IN_PLACE, &samples, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
  passed &= samples == fields * unsigned(cells.y * cells.z);
  MPI_Allreduce(MPI_IN_PLACE, error, 4, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  passed &= error[0] < 2e-12 && error[1] < 1e-8 && error[2] < 2e-12 && error[3] < 2e-12;
  int all = passed ? 1 : 0, rank{};
  MPI_Allreduce(MPI_IN_PLACE, &all, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) std::cout << "ESF mixture fields=" << fields << " status="
      << unsigned(status.code) << '/' << status.detail << " closure=" << error[0]
      << " temperature=" << error[1] << " face_update=" << error[2]
      << " final_update=" << error[3] << " passed=" << all << '\n';
  return all != 0;
}
} // namespace
int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  // Both the physical simplex boundary and its positive interior are covered.
  const bool zero = argc == 2 && std::string(argv[1]) == "--zero";
  bool passed = run(2, zero);
  passed = run(4, zero) && passed;
  MPI_Finalize();
  return passed ? 0 : 1;
}
