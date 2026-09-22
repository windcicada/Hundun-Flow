// SPDX-License-Identifier: Apache-2.0
#include "../../src/models_chemistry_adapter_detail.hpp"
#include "../../src/models_esf_detail.hpp"
#include "../../src/core_product_freeze_detail.hpp"
#include "../support/product_fixture.hpp"
#include "../support/dense_solve.hpp"
#include "hundun/v04_app.hpp"
#include "hundun/v04_portable.hpp"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <limits>
#include <mpi.h>
#include <unistd.h>
#include <vector>

using namespace hundun::v04;
namespace {
// Independent A -> B, equal molecular weights/cp/formation enthalpy. The
// source is -2*rho*Y_A; it cannot affect total gas mass or total enthalpy.
class IsomerGas final : public portable::GasQueryProvider,
                        public portable::GasAdvanceProvider {
public:
  IsomerGas() {
    identity_.mechanism_sha256 = std::string(64, 'a');
    identity_.phase = "synthetic-gas";
    identity_.species_names = {"A", "B"};
    identity_.element_names = {"X"};
    identity_.element_counts = {1, 1};
    identity_.molecular_weights_kg_per_kmol = {28.96546, 28.96546};
    identity_.enthalpy_reference = "nasa7-total-enthalpy";
    identity_.composition_fingerprint = 101;
    closure.element_count=1;
    closure.species={{28.96546,{1},0.},{28.96546,{1},0.}};
    closure.fingerprint=combustion::chemistry_identity_fingerprint(closure);
    identity_.closure_fingerprint = closure.fingerprint;
  }
  const portable::GasIdentity &gas_identity() const noexcept override {
    return identity_;
  }
  portable::Status query_gas(const portable::GasQuery &q,
                             portable::GasQueryOutput &out) noexcept override {
    ++queries;
    if (fail || (fail_at != 0 && queries == fail_at) || q.species_count != 2 || out.capacity != 2)
      return portable::Status::provider_failure;
    const double R = kUniversalGasConstant / 28.96546;
    const double cp = 3.5 * R;
    const double T =
        q.coordinates == portable::GasStateCoordinates::pressure_enthalpy
            ? q.enthalpy_j_per_kg / cp
            : q.temperature_k;
    const double rho = q.pressure_pa / (R * T);
    out.sample = {q.revision,
                  q.composition_fingerprint,
                  q.pressure_pa,
                  T,
                  rho,
                  cp * T,
                  cp,
                  1.8e-5,
                  0.026};
    for (int i = 0; i < 2; ++i) {
      out.diffusivities_m2_per_s[i] = 1e-5;
      out.species_enthalpies_j_per_kg[i] = cp * T;
    }
    out.net_mass_rates_kg_per_m3_s[0] = reverse
        ? 2 * rho * q.mass_fractions[1] : -2 * rho * q.mass_fractions[0];
    out.net_mass_rates_kg_per_m3_s[1] = -out.net_mass_rates_kg_per_m3_s[0];
    if (redirect_output) out.net_mass_rates_kg_per_m3_s = nullptr;
    return portable::Status::success;
  }
  portable::Status advance_gas(const portable::GasAdvanceQuery& q,
      portable::GasAdvanceOutput& out) noexcept override {
    if(++advances==fail_advance_at) return portable::Status::provider_failure;
    double d[2],h[2],w[2],y[2];
    portable::GasQueryOutput sample{{},d,h,w,2};
    auto status=query_gas(q.state,sample);
    if(status!=portable::Status::success) return status;
    const double rho=sample.sample.density_kg_per_m3;
    y[0]=reverse ? 1-(1-q.state.mass_fractions[0])*std::exp(-2*q.duration_s)
                 : q.state.mass_fractions[0]*std::exp(-2*q.duration_s);
    y[1]=1-y[0];
    auto final=q.state; final.mass_fractions=y;
    status=query_gas(final,sample);
    if(status!=portable::Status::success) return status;
    for(unsigned j=0;j<2;++j) {
      out.final_mass_fractions[j]=y[j];
      out.integrated_species_density_delta_kg_per_m3[j]=rho*(y[j]-q.state.mass_fractions[j]);
    }
    out.final_sample=sample.sample; out.completed_duration_s=q.duration_s;
    out.internal_step_count=1; out.integrated_heat_release_j_per_m3=0;
    return portable::Status::success;
  }
  combustion::ChemistryIdentity closure;
  std::uint64_t advances{},fail_advance_at{};
  bool reverse{};
  bool fail{};
  bool redirect_output{};
  std::uint64_t queries{}, fail_at{};

private:
  portable::GasIdentity identity_;
};
class EsfGas final : public portable::GasQueryProvider,
                     public portable::GasAdvanceProvider {
public:
  chemistry::detail::AnalyticIsomerBackend backend;
  bool fail_interval{}, zero_progress{}, hold_composition{};
  unsigned interval_calls{};
  std::uint64_t advance_calls{}, fail_query_after_intervals{};
  double minimum_advanced_enthalpy{std::numeric_limits<double>::infinity()};
  double maximum_advanced_enthalpy{-std::numeric_limits<double>::infinity()};
  const portable::GasIdentity &gas_identity() const noexcept override {
    return backend.gas_identity();
  }
  portable::Status query_gas(const portable::GasQuery &q,
                             portable::GasQueryOutput &out) noexcept override {
    if (fail_query_after_intervals != 0 &&
        advance_calls >= fail_query_after_intervals)
      return portable::Status::provider_failure;
    const auto status = backend.query_gas(q, out);
    if (status == portable::Status::success && zero_progress)
      for (std::size_t s = 0; s < q.species_count; ++s)
        out.net_mass_rates_kg_per_m3_s[s] = 0;
    return status;
  }
  portable::Status
  advance_gas(const portable::GasAdvanceQuery &q,
              portable::GasAdvanceOutput &out) noexcept override {
    ++advance_calls;
    minimum_advanced_enthalpy = std::min(minimum_advanced_enthalpy,
                                        q.state.enthalpy_j_per_kg);
    maximum_advanced_enthalpy = std::max(maximum_advanced_enthalpy,
                                        q.state.enthalpy_j_per_kg);
    if (fail_interval && ++interval_calls == 3)
      return portable::Status::provider_failure;
    // A conservative identity interval isolates the continuation witness from
    // changing eta. Rate queries remain available to construct actual eta/R.
    if (hold_composition) {
      double d[2], h[2], w[2];
      portable::GasQueryOutput sample{{}, d, h, w, 2};
      const auto status = backend.query_gas(q.state, sample);
      if (status != portable::Status::success)
        return status;
      for (unsigned i = 0; i < 2; ++i) {
        out.final_mass_fractions[i] = q.state.mass_fractions[i];
        out.integrated_species_density_delta_kg_per_m3[i] = 0;
      }
      out.final_sample = sample.sample;
      out.completed_duration_s = q.duration_s;
      out.internal_step_count = 1;
      out.integrated_heat_release_j_per_m3 = 0;
      return portable::Status::success;
    }
    return backend.advance_gas(q, out);
  }
  ProductCouplingBindings bindings() noexcept {
    return {this, this, &backend.closure_identity()};
  }
};
std::vector<double> physical_values(const RestartSnapshot &s,
                                    double reference_density = 0) {
  std::vector<double> values{s.time,
                             s.dt,
                             s.pressure_reference,
                             s.previous_pressure_reference,
                             s.closed_mass_target,
                             double(s.step)};
  for (auto fields : {s.fields, s.previous_fields, s.accepted_rate_fields,
                      s.previous_rate_fields})
    for (std::size_t f = 0; f < fields.size; ++f) {
      const auto v = fields.data[f].values;
      for (int z = 0; z < v.interior.z; ++z)
        for (int y = 0; y < v.interior.y; ++y)
          for (int x = 0; x < v.interior.x; ++x)
            for (std::uint8_t c = 0; c < v.components; ++c) {
              double scale = 1.0;
              if (reference_density > 0) {
                const double speed =
                    std::sqrt(s.pressure_reference / reference_density);
                switch (fields.data[f].role) {
                case RestartFieldRole::velocity:
                  scale = speed;
                  break;
                case RestartFieldRole::pressure_perturbation:
                case RestartFieldRole::pressure_absolute:
                  scale = s.pressure_reference;
                  break;
                case RestartFieldRole::enthalpy:
                  scale = speed * speed;
                  break;
                case RestartFieldRole::enthalpy_nonadvective_rate:
                  scale = s.pressure_reference / s.dt;
                  break;
                case RestartFieldRole::scalar_nonadvective_rate:
                  scale = reference_density / s.dt;
                  break;
                default:
                  break;
                }
              }
              values.push_back(v.unchecked({x, y, z}, c) / scale);
            }
    }
  for (auto flux : {s.final_mass_flux, s.previous_mass_flux})
    for (auto v : {flux.x, flux.y, flux.z})
      for (int z = 0; z < v.extents.z; ++z)
        for (int y = 0; y < v.extents.y; ++y)
          for (int x = 0; x < v.extents.x; ++x)
            values.push_back(
                v.unchecked({x, y, z}) /
                (reference_density > 0
                     ? std::sqrt(s.pressure_reference * reference_density)
                     : 1.0));
  for (std::size_t i = 0; i < s.cell_records.values.size;) {
    const auto offset = i % s.cell_records.record_bytes;
    // Rollback and exact restore compare every byte. For subsequent solves,
    // compare the six floating history quantities numerically; counters and
    // signs remain byte exact, including uint64 values above 2^53.
    if (reference_density > 0 && s.cell_records.record_bytes == 120 &&
        ((offset >= 24 && offset < 56) || (offset >= 80 && offset < 96))) {
      std::uint64_t bits = 0;
      for (unsigned b = 0; b < 8; ++b)
        bits |= std::uint64_t(s.cell_records.values.data[i + b]) << (8 * b);
      double value;
      std::memcpy(&value, &bits, 8);
      values.push_back(value);
      i += 8;
    } else
      values.push_back(s.cell_records.values.data[i++]);
  }
  return values;
}
bool collective(bool okay) {
  int a = okay, b = 0;
  MPI_Allreduce(&a, &b, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  return b;
}

bool cn_mean_reaction(ValidatedModel model, int rank, bool pasr = false,
                      double initial_y = .25, bool cn = true, bool reverse = false,
                      bool interval=false) {
  IsomerGas gas;
  gas.reverse = reverse;
  model.time.scheme = cn ? TimeScheme::cn_be : TimeScheme::backward_euler;
  model.solver.coupling = cn ? CouplingKind::outer_corrected : CouplingKind::piso;
  model.legacy_time_fingerprint = cn ? model.fingerprint + 1 : 0;
  model.reaction.mode = pasr ? ReactionMode::pasr_algebraic_v1 : ReactionMode::finite_rate_mean;
  model.reaction.mixing_c_z = .001;
  model.pressure_reference = PressureReferenceKind::boundary_absolute;
  if (cn) model.solver.cold_stopping = ColdStoppingSpec{1.0, 1e-7, 1e-10, 1e-10};
  else model.solver.cold_stopping.reset();
  model.boundaries[0].flow_kind = BoundaryKind::no_slip_wall;
  model.boundaries[0].thermal_kind = BoundaryKind::adiabatic_wall;
  model.boundaries[0].scalars = {{"A", ScalarBoundaryKind::zero_gradient}};
  auto& outlet = model.boundaries[1];
  outlet.flow_kind = BoundaryKind::pressure_outlet;
  outlet.pressure = 101325;
  outlet.allow_backflow = true;
  outlet.backflow_temperature = 300;
  outlet.scalars = {{"A", ScalarBoundaryKind::zero_gradient}};
  outlet.scalars[0].backflow_kind = ScalarBoundaryKind::dirichlet;
  outlet.scalars[0].backflow_value = initial_y;
  CompiledCasePlan plan;
  ProductCouplingBindings bindings{&gas};
  if(interval) { bindings.gas_advance=&gas; bindings.chemistry_identity=&gas.closure; }
  auto status = ProductCompiler::compile(MPI_COMM_WORLD, model, {}, plan, bindings);
  ProductDriver driver;
  if (status) status = ProductDriver::create(MPI_COMM_WORLD, std::move(plan), driver);
  DriverInitialState initial;
  initial.transported_scalars = {&initial_y, 1};
  if (status) status = driver.initialize(initial);
  DriverStepReport report;
  double maximum_error = 0;
  for (unsigned step = 1; step <= 3 && status; ++step) {
    if (cn && step == 3) {
      RestartSnapshot before, after;
      status = driver.committed_restart_snapshot(before);
      if (!collective(bool(status))) break;
      const auto saved = physical_values(before);
      std::uint64_t local_cells = 0;
      for (std::size_t i = 0; i < before.fields.size; ++i) {
        if (before.fields.data[i].role != RestartFieldRole::independent_species)
          continue;
        const auto n = before.fields.data[i].values.interior;
        local_cells = std::uint64_t(n.x) * n.y * n.z;
        break;
      }
      // Fail both accepted-state preparation and endpoint preparation after
      // one cell has overwritten its source. Each retry starts from this same
      // committed snapshot, including all method history and rate fields.
      for (const auto failure_query : {UINT64_C(2), local_cells + 2}) {
        gas.queries = 0;
        gas.fail_at = rank == 0 ? failure_query : 0;
        DriverStepReport failed;
        const auto failure = driver.advance({1, 1, 1, 1, 1}, failed);
        const bool injected = rank != 0 || gas.queries == gas.fail_at;
        gas.fail_at = 0;
        status = driver.committed_restart_snapshot(after);
        const bool rolled_back = collective(injected && !failure &&
            !failed.accepted && status && saved == physical_values(after));
        if (rank == 0)
          std::cerr << "CN_BE partial source rollback query=" << failure_query
                    << " status=" << unsigned(failure.code) << '/' << failure.detail
                    << " unchanged=" << rolled_back << '\n';
        if (!rolled_back) {
          status = {StatusCode::numerical_failure, 99002};
          break;
        }
      }
      if (status && interval) {
        gas.fail_advance_at=rank==0 ? gas.advances+3 : 0;
        DriverStepReport failed;
        const auto failure=driver.advance({1,1,1,1,1},failed);
        const bool injected=rank!=0 || gas.advances==gas.fail_advance_at;
        gas.fail_advance_at=0;
        status=driver.committed_restart_snapshot(after);
        if(!collective(injected && !failure && !failed.accepted && status &&
                       saved==physical_values(after)))
          status={StatusCode::numerical_failure,99002};
        if(rank==0) std::cerr << "CN_BE partial interval rollback status="
            << unsigned(failure.code) << '/' << failure.detail << " passed=" << bool(status) << '\n';
      }
      if (!status) break;
    }
    status = driver.advance({1, 1, 1, 1, 1}, report);
    RestartSnapshot state;
    if (status) status = driver.committed_restart_snapshot(state);
    if (!status) break;
    // The existing finite-rate model freezes its explicit derivative on the
    // accepted state. Outer corrections solve one target interval repeatedly.
    const auto extent = model.mesh.upper;
    const auto lower = model.mesh.lower;
    const auto cells = model.mesh.exact_cells;
    const double volume = (extent.x-lower.x)*(extent.y-lower.y)*(extent.z-lower.z) /
        (double(cells.x)*cells.y*cells.z);
    const double tau_mix = .001 * std::pow(volume, 2.0/3.0) / (2*1e-5);
    const double kappa = pasr ? .5/(.5+tau_mix) : 1.;
    const double decay = interval ? std::exp(-2*model.time.initial_dt*step)
        : std::pow(1 - 2*kappa*model.time.initial_dt, step);
    const double expected = reverse ? 1 - (1 - initial_y)*decay : initial_y*decay;
    bool found = false;
    for (std::size_t i = 0; i < state.fields.size; ++i) {
      const auto& f = state.fields.data[i];
      if (f.role != RestartFieldRole::independent_species) continue;
      found = true;
      const auto n = f.values.interior;
      for (int z = 0; z < n.z; ++z)
        for (int y = 0; y < n.y; ++y)
          for (int x = 0; x < n.x; ++x)
            maximum_error = std::max(maximum_error,
                std::abs(f.values.unchecked({x, y, z}, 0) - expected));
    }
    if (!collective(found && report.accepted && report.piso.cold.active == cn &&
        state.step == step && report.conservation.valid &&
        std::abs(report.conservation.total_energy_balance_defect) < 1e-6)) {
      status = {StatusCode::numerical_failure, 99001};
      break;
    }
    if (interval && cn) {
      const auto& balance=report.conservation;
      bool budget=balance.composition_valid && balance.composition_revision!=0 &&
          balance.composition_duration==model.time.initial_dt &&
          !balance.composition_after_parcel_exchange &&
          balance.species_balance.size()==2 && balance.element_balance.size()==1;
      if(budget) {
        const auto& a=balance.species_balance[0];
        const auto& b=balance.species_balance[1];
        const auto& element=balance.element_balance[0];
        const double mass=a.current_inventory+b.current_inventory;
        const double previous=reverse
            ? 1-(1-initial_y)*std::exp(-2*model.time.initial_dt*(step-1))
            : initial_y*std::exp(-2*model.time.initial_dt*(step-1));
        const double chemistry=mass*(expected-previous)/model.time.initial_dt;
        const double tolerance=2e-9*std::max(1.,std::abs(chemistry));
        budget=a.name=="A" && b.name=="B" && element.name=="X" && mass>0 &&
            std::abs(a.current_inventory/mass-expected)<2e-9 &&
            std::abs(a.accepted_inventory/mass-previous)<2e-9 &&
            (reverse ? a.chemistry_source>0 : a.chemistry_source<0) &&
            std::abs(a.chemistry_source-chemistry)<tolerance &&
            std::abs(a.chemistry_source+b.chemistry_source)<tolerance &&
            std::abs(a.defect)<tolerance && std::abs(b.defect)<tolerance &&
            std::abs(a.transport_outflow)<tolerance &&
            std::abs(element.chemistry_source)<tolerance &&
            std::abs(element.defect)<tolerance;
      }
      if(!collective(budget)) {
        if(rank==0)std::cerr << "mean reaction composition budget failed\n";
        status={StatusCode::numerical_failure,99005};break;
      }
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &maximum_error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  const bool okay = collective(status && maximum_error < 2e-9 && report.accepted_step == 3);
  if (rank == 0)
    std::cerr << "source_reaction cn=" << cn << " interval=" << interval << " pasr=" << pasr << " initial_Y=" << initial_y << " status=" << unsigned(status.code) << '/' << status.detail
              << " step=" << report.accepted_step << " Y_error=" << maximum_error
              << " energy_defect=" << report.conservation.total_energy_balance_defect
              << " passed=" << okay << '\n';
  return okay;
}
// Exercise a changing transport state: formation-enthalpy heat release drives
// expansion, and WALE supplies an active SGS contribution on the next step.
bool cn_pasr_sgs_retry(ValidatedModel model, const std::filesystem::path& assets,
                       int rank, bool interval=false) {
  class FaultGas final : public portable::GasQueryProvider,
                         public portable::GasAdvanceProvider {
  public:
    chemistry::detail::AnalyticIsomerBackend backend;
    std::uint64_t calls{}, fail_at{}, advance_calls{}, fail_advance_at{};
    portable::Status advance_gas(const portable::GasAdvanceQuery& query,
        portable::GasAdvanceOutput& output) noexcept override {
      if (++advance_calls==fail_advance_at) return portable::Status::provider_failure;
      return backend.advance_gas(query,output);
    }
    const portable::GasIdentity& gas_identity() const noexcept override {
      return backend.gas_identity();
    }
    portable::Status query_gas(const portable::GasQuery& query,
                              portable::GasQueryOutput& output) noexcept override {
      if (++calls == fail_at) return portable::Status::provider_failure;
      return backend.query_gas(query, output);
    }
  } gas, reference_gas;
  if (interval) model.reaction.mode=ReactionMode::finite_rate_mean;
  model.time.scheme = TimeScheme::cn_be;
  model.solver.coupling = CouplingKind::outer_corrected;
  model.legacy_time_fingerprint = model.fingerprint + 1;
  model.solver.cold_stopping = ColdStoppingSpec{1., 1e-7, 1e-10, 1e-10};
  model.pressure_reference = PressureReferenceKind::boundary_absolute;
  model.turbulence = TurbulenceKind::wale;
  model.boundaries[0].flow_kind = BoundaryKind::no_slip_wall;
  model.boundaries[0].thermal_kind = BoundaryKind::adiabatic_wall;
  model.boundaries[0].scalars = {{"A", ScalarBoundaryKind::zero_gradient}};
  auto& outlet = model.boundaries[1];
  outlet.flow_kind = BoundaryKind::pressure_outlet;
  outlet.pressure = 101325;
  outlet.allow_backflow = true;
  outlet.backflow_temperature = 300;
  outlet.scalars = {{"A", ScalarBoundaryKind::zero_gradient}};
  outlet.scalars[0].backflow_kind = ScalarBoundaryKind::dirichlet;
  outlet.scalars[0].backflow_value = .25;
  double fraction = .25;
  DriverInitialState initial;
  initial.transported_scalars = {&fraction, 1};
  const auto create = [&](FaultGas& provider, ProductDriver& driver) {
    CompiledCasePlan plan;
    ProductCouplingBindings bindings{&provider};
    if (interval) { bindings.gas_advance=&provider; bindings.chemistry_identity=&provider.backend.closure_identity(); }
    auto status = ProductCompiler::compile(MPI_COMM_WORLD, model, assets, plan, bindings);
    if (status) status = ProductDriver::create(MPI_COMM_WORLD, std::move(plan), driver);
    if (status) status = driver.initialize(initial);
    return status;
  };
  ProductDriver driver, reference;
  auto status = create(gas, driver);
  if (status) status = create(reference_gas, reference);
  DriverStepReport report, reference_report;
  for (unsigned step = 0; step < 2 && status; ++step) {
    status = driver.advance({1,1,1,1,1}, report);
    if (status) status = reference.advance({1,1,1,1,1}, reference_report);
  }
  RestartSnapshot before, after;
  if (status) status = driver.committed_restart_snapshot(before);
  if (!collective(bool(status))) return false;
  RestartSnapshot reference_before;
  status = reference.committed_restart_snapshot(reference_before);
  const auto saved = physical_values(before);
  const bool same_before = collective(status && saved == physical_values(reference_before));
  if (rank == 0) std::cerr << "PaSR WALE pre-failure equal=" << same_before << '\n';
  if (!same_before) return false;
  const auto local = before.patch.cells;
  gas.fail_at = !interval && rank == 0 ? gas.calls + std::uint64_t(local.x)*local.y*local.z + 2 : 0;
  gas.fail_advance_at = interval && rank == 0 ? gas.advance_calls+3 : 0;
  DriverStepReport failed;
  const auto failure = driver.advance({1,1,1,1,1}, failed);
  const bool injected = rank != 0 || (interval ? gas.advance_calls==gas.fail_advance_at : gas.calls==gas.fail_at);
  gas.fail_at = 0; gas.fail_advance_at=0;
  status = driver.committed_restart_snapshot(after);
  if (!collective(injected && !failure && !failed.accepted && status &&
                  saved == physical_values(after))) return false;
  status = driver.advance({1,1,1,1,1}, report);
  if (status) status = reference.advance({1,1,1,1,1}, reference_report);
  RestartSnapshot actual, expected;
  if (status) status = driver.committed_restart_snapshot(actual);
  if (status) status = reference.committed_restart_snapshot(expected);
  if (!collective(bool(status))) return false;
  const double density = 101325 * model.thermophysics.species[0].molecular_weight /
      (kUniversalGasConstant * 300);
  const auto a = physical_values(actual, density), b = physical_values(expected, density);
  double error = 0, maximum_nut = 0, chemistry_error=0;
  // This heated open-domain/WALE state includes spatial transport. Distance
  // from a uniform reactor is diagnostic; cn_mean_reaction tests the exact
  // homogeneous interval oracle independently, including reverse chemistry.
  if(interval) for(std::size_t f=0;f<actual.fields.size;++f) {
    if(actual.fields.data[f].role!=RestartFieldRole::independent_species) continue;
    const auto field=actual.fields.data[f].values;
    for(int z=0;z<field.interior.z;++z) for(int y=0;y<field.interior.y;++y)
      for(int x=0;x<field.interior.x;++x)
        chemistry_error=std::max(chemistry_error,std::abs(field.unchecked({x,y,z},0)-
            fraction*std::exp(-2*3*model.time.initial_dt)));
  }
  MPI_Allreduce(MPI_IN_PLACE,&chemistry_error,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  bool gradients_valid = true;
  if (!collective(a.size() == b.size())) return false;
  std::size_t worst = 0;
  for (std::size_t i = 0; i < a.size(); ++i) {
    const double difference = std::abs(a[i]-b[i])/std::max({1., std::abs(a[i]), std::abs(b[i])});
    if (difference > error) { error = difference; worst = i; }
  }
  double primitive_error = 0;
  std::size_t primitive_count = 6;
  for (std::size_t i = 0; i < actual.fields.size; ++i) {
    const auto field = actual.fields.data[i].values;
    primitive_count += std::size_t(field.interior.x)*field.interior.y*field.interior.z*field.components;
  }
  for (std::size_t i = 0; i < primitive_count; ++i)
    primitive_error = std::max(primitive_error, std::abs(a[i]-b[i]));
  if (rank == 0) std::cerr << "PaSR WALE primitive_error=" << primitive_error << '\n';
  if (error > 1e-12) std::cerr << "PaSR WALE worst rank=" << rank << " index=" << worst
      << " actual=" << a[worst] << " expected=" << b[worst]
      << " time=" << actual.time << '/' << expected.time << " dt=" << actual.dt << '/' << expected.dt << '\n';
  const auto n = model.mesh.exact_cells;
  const double spacing[3]{(model.mesh.upper.x-model.mesh.lower.x)/n.x,
      (model.mesh.upper.y-model.mesh.lower.y)/n.y,
      (model.mesh.upper.z-model.mesh.lower.z)/n.z};
  for (std::size_t f = 0; f < actual.fields.size; ++f) {
    if (actual.fields.data[f].role != RestartFieldRole::velocity) continue;
    const auto velocity = actual.fields.data[f].values;
    for (int z = 1; z+1 < local.z; ++z)
      for (int y = 1; y+1 < local.y; ++y)
        for (int x = 1; x+1 < local.x; ++x) {
          VelocityGradient gradient;
          for (unsigned direction = 0; direction < 3; ++direction) {
            Int3 minus{x,y,z}, plus{x,y,z};
            if (direction == 0) { --minus.x; ++plus.x; }
            if (direction == 1) { --minus.y; ++plus.y; }
            if (direction == 2) { --minus.z; ++plus.z; }
            for (unsigned component = 0; component < 3; ++component)
              gradient.value[3*component+direction] =
                  (velocity.unchecked(plus,component)-velocity.unchecked(minus,component))/(2*spacing[direction]);
          }
          double nut{};
          gradients_valid &= bool(wale_kinematic_viscosity(gradient,
              std::cbrt(spacing[0]*spacing[1]*spacing[2]), TurbulencePlanSpec{}.wale_coefficient, nut));
          maximum_nut = std::max(maximum_nut, nut);
        }
  }
  MPI_Allreduce(MPI_IN_PLACE, &error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &maximum_nut, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  const bool passed = collective(gradients_valid && error <= 1e-12 && maximum_nut > 0 &&
      report.accepted && reference_report.accepted && report.conservation.valid &&
      reference_report.conservation.valid && report.accepted_step == 3);
  if (rank == 0)
    std::cerr << "CN_BE WALE interval=" << interval << " uniform_reactor_distance=" << chemistry_error
              << " retry error=" << error << " max_nut=" << maximum_nut
              << " energy_defect=" << report.conservation.total_energy_balance_defect
              << " passed=" << passed << '\n';
  return passed;
}
class FoldWitness final : public ProductTcrFoldProvider {
public:
  int defect{};
  PlanFingerprint identity{UINT64_C(9007199254741011)};
  PlanFingerprint fingerprint() const noexcept override { return identity; }
  Status query(const ProductTcrFoldQuery &q,
               ProductTcrFoldEvidence &out) const noexcept override {
    out = {};
    // The provider failure occurs after the first local cell, while identity
    // defects continue to exercise the origin's continuation witness.
    if (q.global_cell != (defect == 4 ? 1U : 0U) || !q.initialized ||
        q.history_revision.accepted_step != 1)
      return {};
    if (defect == 4)
      return {StatusCode::numerical_failure, 12345};
    // Manufactured continuation at eta=1/4: R rises from 1 to 4/3 and
    // returns to 1 on the upper branch. It is a contract witness, not CFD
    // evidence for an automatically detected physical fold.
    out = {true,
           identity,
           q.global_cell,
           q.history_revision,
           q.input_revision,
           .25,
           4. / 3,
           1};
    if (defect == 1)
      ++out.global_cell;
    if (defect == 2)
      ++out.input_revision.input_revision;
    if (defect == 3)
      ++out.source_identity;
    if (defect == 5)
      ++out.history_revision.input_revision;
    if (defect == 6)
      out.rate_ratio = 1;
    return {};
  }
};
bool native_fold_continuation(const ValidatedModel &model,
                              const std::filesystem::path &case_root,
                              const std::filesystem::path &root,
                              DriverInitialState initial, int rank) {
  EsfGas gas, restored_gas;
  gas.hold_composition = restored_gas.hold_composition = true;
  FoldWitness witness, restored_witness;
  ProductDriver driver, restored;
  auto create = [&](EsfGas &backend, FoldWitness &fold, ProductDriver &out) {
    auto bindings = backend.bindings();
    bindings.tcr_fold = &fold;
    CompiledCasePlan plan;
    auto status = ProductCompiler::compile(MPI_COMM_WORLD, model, case_root,
                                           plan, bindings);
    if (status)
      status = ProductDriver::create(MPI_COMM_WORLD, std::move(plan), out);
    return status;
  };
  auto status = create(gas, witness, driver);
  if (status)
    status = driver.initialize(initial);
  DriverStepReport report;
  if (status)
    status = driver.advance({1, 1, 1, 1, 1}, report);
  if (!collective(bool(status) && report.accepted))
    return false;
  RestartSnapshot before;
  status = driver.committed_restart_snapshot(before);
  if (!collective(bool(status)))
    return false;
  const auto saved = physical_values(before);
  const auto variance_at_origin = [](const RestartSnapshot &snapshot) {
    double mean = 0, variance = 0;
    unsigned n = 0;
    for (std::size_t i = 0; i < snapshot.fields.size; ++i)
      if (snapshot.fields.data[i].role == RestartFieldRole::stochastic_field) {
        mean += snapshot.fields.data[i].values.unchecked({0, 0, 0}, 0);
        ++n;
      }
    mean /= n;
    for (std::size_t i = 0; i < snapshot.fields.size; ++i)
      if (snapshot.fields.data[i].role == RestartFieldRole::stochastic_field) {
        const double delta =
            snapshot.fields.data[i].values.unchecked({0, 0, 0}, 0) - mean;
        variance += delta * delta / n;
      }
    return variance;
  };
  const double old_variance = variance_at_origin(before);
  double viscosity = 0;
  for (std::size_t i = 0; i < before.fields.size; ++i)
    if (before.fields.data[i].role == RestartFieldRole::stochastic_transport)
      viscosity = before.fields.data[i].values.unchecked({0, 0, 0}, 2) +
                  before.fields.data[i].values.unchecked({0, 0, 0}, 3);
  const double rho = initial.pressure_reference * 28 /
                     (kUniversalGasConstant * initial.temperature);
  const double delta =
      std::cbrt((model.mesh.upper.x - model.mesh.lower.x) *
                (model.mesh.upper.y - model.mesh.lower.y) *
                (model.mesh.upper.z - model.mesh.lower.z) /
                (double(before.global_cells.x) * before.global_cells.y *
                 before.global_cells.z));
  // At kappa=1, uniform fields, zero advection/turbulence and identity
  // chemistry, the target variance follows this exact IEM decay.
  const double mixing=before.dt*viscosity/(rho*model.reaction.mixing_c_z*delta*delta);
  const double expected_variance=old_variance *
      (model.reaction.esf->tcr.mode==TcrMode::experimental
          ? std::exp(-2*mixing) : 1/std::pow(1+mixing,2));

  for (int defect = 1; defect <= 6; ++defect) {
    witness.defect = rank == 0 ? defect : 0;
    const auto chemical_calls = gas.advance_calls;
    status = driver.advance({1, 1, 1, 1, 1}, report);
    bool rejected = !status && !report.accepted && gas.advance_calls == chemical_calls;
    RestartSnapshot after;
    if (rejected)
      status = driver.committed_restart_snapshot(after);
    rejected &= bool(status) && physical_values(after) == saved;
    if (!collective(rejected)) {
      if (rank == 0)
        std::cerr << "fold witness defect did not roll back: " << defect
                  << '\n';
      return false;
    }
  }
  witness.defect = 0;
  status = driver.advance({1, 1, 1, 1, 1}, report);
  if (!collective(bool(status) && report.accepted))
    return false;
  RestartSnapshot folded;
  status = driver.committed_restart_snapshot(folded);
  if (!collective(bool(status)))
    return false;
  const auto u64 = [](const std::uint8_t *p, unsigned n = 8) {
    std::uint64_t v = 0;
    for (unsigned b = 0; b < n; ++b)
      v |= std::uint64_t(p[b]) << (8 * b);
    return v;
  };
  const auto real = [&](const std::uint8_t *p) {
    const auto bits = u64(p);
    double v;
    std::memcpy(&v, &bits, 8);
    return v;
  };
  bool valid = true;
  const auto cells = folded.patch.cells;
  for (int z = 0; z < cells.z; ++z)
    for (int y = 0; y < cells.y; ++y)
      for (int x = 0; x < cells.x; ++x) {
        const auto i = std::size_t(x) +
                       std::size_t(cells.x) * (y + std::size_t(cells.y) * z);
        const auto *row = folded.cell_records.values.data + 120 * i;
        const bool target = x + folded.patch.begin.x == 0 &&
                            y + folded.patch.begin.y == 0 &&
                            z + folded.patch.begin.z == 0;
        valid &= u64(row) == 2 && u64(row + 72) == (target ? 1 : 0);
        if (target)
          valid &=
              u64(row + 56, 4) == 2 && std::abs(real(row + 40) - .5) < 2e-12 &&
              std::abs(real(row + 48) - 1) < 2e-12 && real(row + 80) == .25 &&
              real(row + 88) == 4. / 3 && u64(row + 96) == 1;
      }
  if (!collective(valid))
    return false;
  const bool owns_origin = folded.patch.begin.x == 0 &&
                           folded.patch.begin.y == 0 &&
                           folded.patch.begin.z == 0;
  if (!collective(!owns_origin || std::abs(variance_at_origin(folded) -
                                           expected_variance) < 2e-12)) {
    if (rank == 0)
      std::cerr << "fold departure did not drive native IEM variance\n";
    return false;
  }
  const auto folded_values = physical_values(folded);
  status = RestartWriter::write(MPI_COMM_WORLD, root, folded);
  if (status)
    status = create(restored_gas, restored_witness, restored);
  RestartExpected expected;
  if (status)
    status = restored.restart_expected(expected);
  RestartImage image;
  if (status)
    status = RestartReader::load(MPI_COMM_WORLD, root, expected, image);
  if (status)
    status = restored.initialize_restart(image);
  RestartSnapshot after;
  if (status)
    status = restored.committed_restart_snapshot(after);
  if (!collective(bool(status) && physical_values(after) == folded_values))
    return false;
  status = restored.advance({1, 1, 1, 1, 1}, report);
  if (!collective(bool(status) && report.accepted))
    return false;
  status = restored.committed_restart_snapshot(after);
  if (!collective(bool(status)))
    return false;
  for (std::size_t i = 0; i < after.cell_records.values.size; i += 120)
    valid &= u64(after.cell_records.values.data + i) == 3 &&
             u64(after.cell_records.values.data + i + 72) ==
                 u64(folded.cell_records.values.data + i + 72);
  if (!collective(valid))
    return false;
  // A different evidence authority cannot restore this accepted branch history.
  ProductDriver wrong;
  restored_witness.identity += 2;
  status = create(restored_gas, restored_witness, wrong);
  if (status)
    status = wrong.initialize_restart(image);
  return collective(!status);
}
// Check the state actually supplied to chemistry, before the final mean-field
// reconciliation can hide a missing stochastic-field wall heat source.
bool thermal_esf(ValidatedModel model, const std::filesystem::path& assets,
                 int rank, double heat_flux) {
  model.fingerprint += 43;
  model.boundaries[2].flow_kind = BoundaryKind::no_slip_wall;
  model.boundaries[2].thermal_kind = BoundaryKind::heat_flux_wall;
  model.boundaries[2].heat_flux = heat_flux;
  model.boundaries[3].flow_kind = BoundaryKind::no_slip_wall;
  model.boundaries[3].thermal_kind = BoundaryKind::adiabatic_wall;
  for (unsigned face : {2U, 3U})
    for (const auto& scalar : model.transported_scalars)
      model.boundaries[face].scalars.push_back(
          {scalar.stable_name, ScalarBoundaryKind::zero_gradient, 0.0});
  EsfGas gas;
  CompiledCasePlan plan;
  ProductDriver driver;
  auto status = ProductCompiler::compile(MPI_COMM_WORLD, model, assets, plan,
                                         gas.bindings());
  if (status) status = ProductDriver::create(MPI_COMM_WORLD, std::move(plan), driver);
  double fraction = 0.25;
  DriverInitialState initial;
  initial.pressure_reference = 101325;
  initial.temperature = 300;
  initial.transported_scalars = {&fraction, 1};
  if (status) status = driver.initialize(initial);
  double h0{}, cp{}, gas_constant{};
  ThermodynamicsPlan thermo;
  if (status) status = ThermodynamicsPlan::compile(model.thermophysics,
      {model.transported_scalars.data(), model.transported_scalars.size()}, thermo);
  if (status) status = thermo.mixture_enthalpy(initial.temperature,
      {&fraction, 1}, h0, cp, gas_constant);
  const double rho = initial.pressure_reference *
      model.thermophysics.species[0].molecular_weight /
      (kUniversalGasConstant * initial.temperature);
  TransportPlan transport;
  MolecularTransportState material;
  if (status) status=TransportPlan::compile(model.thermophysics,thermo,transport);
  if (status) status=transport.evaluate(initial.temperature,{&fraction,1},material);
  double gamma=material.conductivity/cp;
  if (status && transport.has_effective_enthalpy_transport())
    status=transport.effective_enthalpy_transport(material.viscosity,material.viscosity,
                                                cp,material.conductivity,gamma);
  DriverStepReport report;
  if (status) status = driver.advance({1,1,1,1,1}, report);
  double observed = heat_flux > 0 ? gas.minimum_advanced_enthalpy : gas.maximum_advanced_enthalpy;
  MPI_Allreduce(MPI_IN_PLACE, &observed, 1, MPI_DOUBLE,
                heat_flux > 0 ? MPI_MIN : MPI_MAX, MPI_COMM_WORLD);
  const double dy = (model.mesh.upper.y - model.mesh.lower.y) / model.mesh.exact_cells.y;
  // The transported mean and each identical h field solve the same implicit
  // wall diffusion problem. Solve for delta h to retain FP64 precision.
  const int ny=model.mesh.exact_cells.y;
  const double d=model.time.initial_dt*gamma/(rho*dy*dy);
  std::vector<std::vector<double>> matrix(ny,std::vector<double>(ny));
  std::vector<double> rhs(ny);
  for(int j=0;j<ny;++j) {
    matrix[j][j]=1;
    if(j) {matrix[j][j]+=d;matrix[j][j-1]=-d;}
    if(j+1<ny) {matrix[j][j]+=d;matrix[j][j+1]=-d;}
  }
  rhs[0]=-model.time.initial_dt*heat_flux/(rho*dy);
  const double expected=h0+test::dense_solve(matrix,rhs)[0];
  const double heat_input = -heat_flux * (model.mesh.upper.x-model.mesh.lower.x) *
      (model.mesh.upper.z-model.mesh.lower.z);
  const bool passed = status && report.accepted && rho > 0 &&
      std::isfinite(observed) && std::abs(observed - expected) < 1e-8 &&
      report.conservation.valid &&
      std::abs(report.conservation.total_energy_balance_defect) <= 1e-6 * std::abs(heat_input);
  if (rank == 0)
    std::cerr << "esf_heat status=" << unsigned(status.code) << '/' << status.detail
              << " delta=" << observed-h0 << " expected=" << expected-h0
              << " energy_defect=" << report.conservation.total_energy_balance_defect
              << " passed=" << passed << '\n';
  return collective(passed);
}
// Exercise the production writer ceiling on a patch large enough that the
// face envelope and fixed metadata allowance expose an omitted field bank.
bool restart_capacity(ValidatedModel model, const std::filesystem::path& assets,
                      const std::filesystem::path& root, int rank) {
  model.mesh.exact_cells = {64, 64, 64};
  model.mesh.minimum_spacing = {1.0 / 64, 1.0 / 64, 1.0 / 64};
  model.mesh.limits.max_global_cells = 64U * 64U * 64U;
  model.mesh.limits.max_memory_bytes_per_rank = UINT64_C(1073741824);
  model.reaction.esf->fields = 2;
  model.reaction.esf->initial_species_offsets = {0.1, -0.1};
  EsfGas gas;
  gas.hold_composition = true;
  CompiledCasePlan plan;
  auto status = ProductCompiler::compile(MPI_COMM_WORLD, model, assets, plan,
                                         gas.bindings());
  std::size_t budget{};
  if (status) {
    const auto services = plan.io_services()->services();
    for (std::size_t i = 0; i < services.size; ++i)
      if (services.data[i].kind == RuntimeServiceKind::restart)
        budget = services.data[i].maximum_staging_bytes_per_rank;
  }
  ProductDriver driver;
  if (status) status = ProductDriver::create(MPI_COMM_WORLD, std::move(plan), driver);
  RestartExpected expected;
  if (status) status = driver.restart_expected(expected);
  double fraction = 0.25;
  DriverInitialState initial;
  initial.transported_scalars = {&fraction, 1};
  if (status) status = driver.initialize(initial);
  DriverStepReport step;
  if (status) status = driver.advance({1, 1, 1, 1, 1}, step);
  RestartSnapshot snapshot;
  if (status) status = driver.committed_restart_snapshot(snapshot);
  RestartWriteReport write;
  if (status && budget)
    status = RestartWriter::write(MPI_COMM_WORLD, root, snapshot, {1, &write, budget});
  bool passed = status && budget && write.rank_payload_bytes > 0 &&
                write.peak_bulk_staging_bytes <= budget;
  RestartImage image;
  if (passed && status)
    status = RestartReader::load(MPI_COMM_WORLD, root, expected, image);
  if (passed && status) {
    // Compare every serialized field, including independent field0, at both
    // state levels and both nonadvective rate levels.
    const auto equal = [](Span<const RestartFieldView> views,
                          const std::vector<RestartImageField>& fields) {
      if (views.size != fields.size()) return false;
      for (std::size_t f = 0; f < views.size; ++f) {
        const auto v = views.data[f].values;
        const auto& a = fields[f];
        if (views.data[f].role != a.role || v.field != a.field ||
            v.components != a.components) return false;
        std::size_t i{};
        for (int z = 0; z < v.interior.z; ++z)
          for (int y = 0; y < v.interior.y; ++y)
            for (int x = 0; x < v.interior.x; ++x)
              for (unsigned c = 0; c < v.components; ++c)
                if (i >= a.values.size() ||
                    a.values[i++] != v.unchecked({x, y, z}, c)) return false;
        if (i != a.values.size()) return false;
      }
      return true;
    };
    passed = equal(snapshot.fields, image.fields) &&
             equal(snapshot.previous_fields, image.previous_fields) &&
             equal(snapshot.accepted_rate_fields, image.accepted_rate_fields) &&
             equal(snapshot.previous_rate_fields, image.previous_rate_fields);
  }
  passed = collective(passed && bool(status));
  if (rank == 0) {
    std::cerr << "restart_capacity status=" << unsigned(status.code) << '/'
              << status.detail << " budget=" << budget
              << " payload=" << write.rank_payload_bytes
              << " peak=" << write.peak_bulk_staging_bytes
              << " passed=" << passed << '\n';
    std::error_code error;
    std::filesystem::remove_all(root, error);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  return passed;
}
} // namespace
int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int id = int(getpid());
  MPI_Bcast(&id, 1, MPI_INT, 0, MPI_COMM_WORLD);
  const auto root = std::filesystem::temp_directory_path() /
                    ("hundun-product-reaction-" + std::to_string(id));
  bool okay = true;
  {
    IsomerGas gas;
    auto model = test::product_model({8, 8, 8});
    model.turbulence = TurbulenceKind::none;
    model.time.control = TimeControlKind::fixed;
    model.time.initial_dt = model.time.minimum_dt = model.time.maximum_dt =
        0.001;
    model.time.maximum_retries = 1;
    model.thermophysics.species[0].stable_name = "A";
    auto b = model.thermophysics.species[0];
    b.stable_name = "B";
    model.thermophysics.species.push_back(b);
    model.transported_scalars.push_back(
        {"A", TransportedScalarRole::species, 1.0, 1.0});
    model.reaction.mode = ReactionMode::finite_rate_mean;
    model.reaction.mechanism_sha256 = gas.gas_identity().mechanism_sha256;
    model.reaction.phase = gas.gas_identity().phase;
    if (argc == 1 && (!cn_mean_reaction(model, rank, false,.25,true,false,true) ||
        !cn_mean_reaction(model, rank, false,0.,true,true,true) ||
        !cn_mean_reaction(model, rank, false, 1.) ||
        !cn_mean_reaction(model, rank, false, 0., true, true) ||
        !cn_mean_reaction(model, rank) ||
        !cn_mean_reaction(model, rank, true, 0., false) ||
        !cn_mean_reaction(model, rank, true) ||
        !cn_mean_reaction(model, rank, true, 0.))) {
      MPI_Finalize();
      return 1;
    }
    Status status;
    const std::filesystem::path case_root = argc > 1 ? argv[1] : "";
    if (!case_root.empty())
      status = CaseCompiler::load_and_compile(MPI_COMM_WORLD, case_root, model);
    if (status && !case_root.empty() && model.reaction.mode == ReactionMode::pasr_algebraic_v1 &&
        (!cn_pasr_sgs_retry(model, case_root, rank) ||
         !cn_pasr_sgs_retry(model, case_root, rank, true))) {
      MPI_Finalize();
      return 1;
    }
    if (argc > 2 && std::string_view(argv[2]) == "--restart-capacity") {
      const bool passed = status && restart_capacity(model, case_root, root, rank);
      MPI_Finalize();
      return passed ? 0 : 1;
    }
    EsfGas esf_gas;
    const bool esf = model.reaction.mode == ReactionMode::esf_tpdf;
    const bool external_esf =
        esf && model.reaction.representation !=
                   ReactionSpec::Representation::direct_cantera;
    const bool tcr = esf && model.reaction.esf->tcr.mode != TcrMode::off;
    if (external_esf && !tcr &&
        (!thermal_esf(model, case_root, rank, 10.0) ||
         !thermal_esf(model, case_root, rank, -10.0))) {
      MPI_Finalize();
      return 1;
    }
    CompiledCasePlan plan;
    if (status)
      status = ProductCompiler::compile(
          MPI_COMM_WORLD, model, case_root, plan,
          external_esf ? esf_gas.bindings()
                       : (case_root.empty() ? ProductCouplingBindings{&gas}
                                            : ProductCouplingBindings{}));
    okay = collective(static_cast<bool>(status));
    if (!okay && rank == 0)
      std::cerr << "reaction compile failed " << unsigned(status.code) << ":"
                << status.detail << '\n';
    ProductDriver driver;
    if (okay)
      status = ProductDriver::create(MPI_COMM_WORLD, std::move(plan), driver);
    if (okay)
      okay = collective(static_cast<bool>(status));
    double initial_fraction = 0.25;
    DriverInitialState initial;
    initial.transported_scalars = {&initial_fraction, 1};
    if (okay) {
      status = driver.initialize(initial);
      okay = collective(static_cast<bool>(status));
      if (!okay && rank == 0)
        std::cerr << "initialize " << unsigned(status.code) << ":"
                  << status.detail << "\n";
    }
    DriverStepReport report;
    if (okay) {
      status = driver.advance({1, 1, 1, 1, 1}, report);
      okay = collective(static_cast<bool>(status) && report.accepted);
      if (!okay && rank == 0)
        std::cerr << "advance " << unsigned(status.code) << ":" << status.detail
                  << "\n";
    }
    RestartSnapshot snap;
    if (okay)
      okay = collective(
          static_cast<bool>(driver.committed_restart_snapshot(snap)));
    double effective_rate = 2.0;
    if (model.reaction.mode == ReactionMode::pasr_algebraic_v1) {
      const double volume = (model.mesh.upper.x - model.mesh.lower.x) *
                            (model.mesh.upper.y - model.mesh.lower.y) *
                            (model.mesh.upper.z - model.mesh.lower.z);
      const double delta = std::cbrt(
          volume / (model.mesh.exact_cells.x * model.mesh.exact_cells.y *
                    model.mesh.exact_cells.z));
      const double tau_mix =
          model.reaction.mixing_c_z * delta * delta / (2 * 2e-5);
      effective_rate *= 0.5 / (0.5 + tau_mix);
    }
    double accepted_y = 0;
    if (okay) {
      bool found = false;
      for (std::size_t i = 0; i < snap.fields.size; ++i) {
        const auto &f = snap.fields.data[i];
        if (f.role != RestartFieldRole::independent_species)
          continue;
        accepted_y = f.values.unchecked({0, 0, 0}, 0);
        found = true;

        // Independent analytic exp(-2dt), with first-start BE/explicit-source
        // truncation bounded by 0.5*Y0*(2dt)^2 = 5e-7.
        okay &=
            std::abs(accepted_y - 0.25 * std::exp(-effective_rate * 0.001)) <
            0.5 * 0.25 * std::pow(effective_rate * 0.001, 2) + 1e-10;
        okay &= accepted_y < 0.25 && accepted_y > 0.249;
      }
      okay = collective(okay && found && snap.step == 1);
    }
    if (okay && !case_root.empty()) {
      double h = 0, pi = 0;
      for (std::size_t i = 0; i < snap.fields.size; ++i) {
        const auto &f = snap.fields.data[i];
        if (f.role == RestartFieldRole::enthalpy)
          h = f.values.unchecked({0, 0, 0}, 0);
        if (f.role == RestartFieldRole::pressure_perturbation)
          pi = f.values.unchecked({0, 0, 0}, 0);
      }
      const double R = kUniversalGasConstant / model.thermophysics.species[0].molecular_weight, cp = 1000.0,
                   heat = 100000.0;
      const double rho0 =
          initial.pressure_reference / (R * initial.temperature);
      const double h0 =
          cp * (initial.temperature - 298.15) + heat * initial_fraction;
      const double pressure = snap.pressure_reference + pi;
      const double temperature = 298.15 + (h - heat * accepted_y) / cp;
      const double expected_temperature =
          initial.temperature +
          heat * (initial_fraction - accepted_y) / (cp - R);
      okay =
          collective(temperature > initial.temperature &&
                     std::abs(temperature - expected_temperature) < 1e-7 &&
                     std::abs(rho0 * (h - h0) -
                              (pressure - initial.pressure_reference)) < 1e-5 &&
                     std::abs(pressure / (R * temperature) - rho0) < 1e-10);
    }
    if (okay && esf) {
      const double R = kUniversalGasConstant /
                       model.thermophysics.species[0].molecular_weight;
      const double rho = initial.pressure_reference / (R * initial.temperature);
      const double volume =
          (model.mesh.upper.x - model.mesh.lower.x) *
          (model.mesh.upper.y - model.mesh.lower.y) *
          (model.mesh.upper.z - model.mesh.lower.z) /
          (model.mesh.exact_cells.x * model.mesh.exact_cells.y *
           model.mesh.exact_cells.z);
      const double viscosity = 1e-5; // IEM Cd=2 uses mu(300 K), independently of Pr/Sc.
      const double tau = model.reaction.mixing_c_z *
                         std::pow(std::cbrt(volume), 2) / (2 * viscosity / rho);
      const double control =
          tcr && model.reaction.esf->tcr.mode == TcrMode::experimental
              ? 1.0 / 3.0
              : 1.0;
      const double dt = snap.dt,
                   relaxation = tcr && model.reaction.esf->tcr.mode == TcrMode::experimental
                       ? std::exp(-std::cbrt(control) * dt / (2 * tau))
                       : 1/(1+dt/(2*tau)),
                   decay = std::exp(-2*dt);
      double fraction_delta = 0, expected_variance = 0;
      for (double offset : model.reaction.esf->initial_species_offsets) {
        const double y0 = 0.25 + offset * relaxation, ye = y0 * decay;
        // Equal-weight Favre samples define the mean endpoint increment;
        // the common carrier density belongs to the flow source assembly.
        fraction_delta += (ye - y0) / 4;
        expected_variance += std::pow(offset * relaxation * decay, 2) / 4;
      }
      if (rank == 0 && !external_esf &&
          std::abs(accepted_y - (0.25 + fraction_delta)) >= 2e-12)
        std::cerr << std::scientific << "Cantera ESF source error="
                  << accepted_y - (0.25 + fraction_delta) << "\n";
      okay &= std::abs(accepted_y - (0.25 + fraction_delta)) < 2e-12;
      std::vector<ConstFieldView> ensemble;
      for (std::size_t i = 0; i < snap.fields.size; ++i)
        if (snap.fields.data[i].role == RestartFieldRole::stochastic_field)
          ensemble.push_back(snap.fields.data[i].values);
      okay &= ensemble.size() == 4;
      if (ensemble.size() == 4) {
        const auto cells = ensemble[0].interior;
        for (int z = 0; z < cells.z; ++z)
          for (int y = 0; y < cells.y; ++y)
            for (int x = 0; x < cells.x; ++x) {
              double mean = 0, variance = 0;
              for (auto field : ensemble) {
                okay &= field.components == 3;
                const double a = field.unchecked({x, y, z}, 0),
                             b = field.unchecked({x, y, z}, 1);
                mean += a / 4;
                variance += std::pow(a - accepted_y, 2) / 4;
                okay &= a >= 0 && b >= 0 && std::abs(a + b - 1) < 2e-12;
              }
              okay &= std::abs(mean - accepted_y) < 2e-12;
              if (rank == 0 && !external_esf && x == 0 && y == 0 && z == 0 &&
                  std::abs(variance - expected_variance) >= 2e-12)
                std::cerr << "Cantera ESF variance error="
                          << variance - expected_variance << "\n";
              okay &= std::abs(variance - expected_variance) < 2e-12;
            }
      }
      okay = collective(okay);
    }
    if (!okay && rank == 0)
      std::cerr << "first-step physical oracle failed, Y=" << accepted_y
                << "\n";
    if (okay && tcr) {
      // Public product record ABI: little-endian step/input revision,
      // initialized flag, eta/R, signed root and kappa. A -> B gives eta=1/4,
      // R=1, negative root=-1/2 and kappa=1/3, independently of field variance.
      const auto records = snap.cell_records;
      const auto u64 = [&](std::size_t offset) {
        std::uint64_t value = 0;
        for (unsigned i = 0; i < 8; ++i)
          value |= std::uint64_t(records.values.data[offset + i]) << (8 * i);
        return value;
      };
      const auto real = [&](std::size_t offset) {
        const auto bits = u64(offset);
        double value;
        std::memcpy(&value, &bits, 8);
        return value;
      };
      okay &= records.identity != 0 && records.record_bytes == 120U &&
              records.values.size > 0;
      if (okay)
        for (std::size_t cell = 0; cell < records.values.size;
             cell += records.record_bytes)
          okay &= u64(cell) == 1 && u64(cell + 8) == 2 &&
                  records.values.data[cell + 20] == 1 &&
                  std::abs(real(cell + 24) - .25) < 2e-12 &&
                  std::abs(real(cell + 32) - 1) < 2e-12 &&
                  std::abs(real(cell + 40) + .5) < 2e-12 &&
                  std::abs(real(cell + 48) - 1.0 / 3.0) < 2e-12;
      if (!okay && rank == 0)
        std::cerr << "TCR record oracle: identity=" << records.identity
                  << " width=" << records.record_bytes
                  << " bytes=" << records.values.size
                  << " step=" << (records.values.size ? u64(0) : 0)
                  << " input=" << (records.values.size ? u64(8) : 0)
                  << " eta=" << (records.values.size ? real(24) : 0)
                  << " R=" << (records.values.size ? real(32) : 0)
                  << " root=" << (records.values.size ? real(40) : 0)
                  << " kappa=" << (records.values.size ? real(48) : 0) << "\n";
      okay = collective(okay);
    }
    // Compare a real disk restart with uninterrupted continuation.
    IsomerGas restored_gas;
    EsfGas restored_esf;
    ProductDriver restored;
    if (okay) {
      status = RestartWriter::write(MPI_COMM_WORLD, root, snap);
      CompiledCasePlan restored_plan;
      if (status)
        status = ProductCompiler::compile(
            MPI_COMM_WORLD, model, case_root, restored_plan,
            external_esf
                ? restored_esf.bindings()
                : (case_root.empty() ? ProductCouplingBindings{&restored_gas}
                                     : ProductCouplingBindings{}));
      if (status)
        status = ProductDriver::create(MPI_COMM_WORLD, std::move(restored_plan),
                                       restored);
      RestartExpected expected;
      if (status)
        status = restored.restart_expected(expected);
      RestartImage image;
      if (status)
        status = RestartReader::load(MPI_COMM_WORLD, root, expected, image);
      if (status && esf) {
        // Native restore must reject a rank-local invalid field before any
        // primary state is installed, then accept the intact image.
        const std::size_t field = 3 + model.transported_scalars.size();
        const double value = image.fields[field].values[0];
        if (rank == 0)
          image.fields[field].values[0] =
              std::numeric_limits<double>::quiet_NaN();
        const auto rejected = restored.initialize_restart(image);
        okay = collective(!rejected);
        image.fields[field].values[0] = value;
        // Both accepted material histories participate in restore validation.
        for (auto* layer : {&image.fields,&image.previous_fields}) {
          auto& cache=*std::find_if(layer->begin(),layer->end(),[](const auto& field) {
            return field.role==RestartFieldRole::stochastic_transport;
          });
          for (unsigned component : {2U,3U}) {
            const double saved=cache.values[component];
            if (rank==0) cache.values[component]=component==2
                ? std::numeric_limits<double>::quiet_NaN() : -1e-5;
            const auto rejected_cache=restored.initialize_restart(image);
            const bool valid_rejection=collective(!rejected_cache);
            okay=okay && valid_rejection;
            cache.values[component]=saved;
          }
        }
      }
      if (status && tcr) {
        for (const auto offset : {0U, 64U}) {
          const auto old = image.cell_records[offset];
          if (rank == 0)
            image.cell_records[offset] ^= 1U; // wrong step or mapping identity
          const auto rejected = restored.initialize_restart(image);
          okay = collective(okay && !rejected);
          image.cell_records[offset] = old;
        }
      }
      if (status)
        status = restored.initialize_restart(image);
      if (status) {
        RestartSnapshot check;
        status = restored.committed_restart_snapshot(check);
        okay = collective(bool(status) &&
                          physical_values(check) == physical_values(snap));
      }
      if (!status && rank == 0)
        std::cerr << "restart status " << unsigned(status.code) << ":"
                  << status.detail << "\n";
      okay = collective(okay && bool(status));
      if (!okay && rank == 0)
        std::cerr << "restart exact payload/rejection failed\n";
    }
    if (okay) {
      // Match the linear initial guess of native restore. Restart transports
      // the physical state; pressure warm-start efficacy has its own same-state
      // solver regression in core_product_freeze_mpi_test.
      detail::suppress_pressure_correction_warm_start_once_for_test();
      status = driver.advance({1, 1, 1, 1, 1}, report);
      if (status)
        status = driver.committed_restart_snapshot(snap);
      const double density_scale =
          case_root.empty()
              ? 0.0
              : initial.pressure_reference *
                    model.thermophysics.species[0].molecular_weight /
                    (kUniversalGasConstant * initial.temperature);
      const auto continuous = physical_values(snap, density_scale);
      if (status)
        status = restored.advance({1, 1, 1, 1, 1}, report);
      RestartSnapshot resumed;
      if (status)
        status = restored.committed_restart_snapshot(resumed);
      const auto recovered = physical_values(resumed, density_scale);
      // Restart preserves stored accepted values exactly. Compare subsequent
      // pressure/enthalpy solves with matching linear initial guesses and the
      // original terminal tolerance.
      // SI fields are normalized by p_ref, sqrt(p_ref/rho), p_ref/rho,
      // p_ref/dt, rho/dt, and rho*sqrt(p_ref/rho), respectively. Comparing
      // tiny pi or zero-flow velocity against a unit Pa or m/s scale would
      // confuse solver-relative tolerance with an absolute tolerance.
      bool same = continuous.size() == recovered.size();
      for (std::size_t i = 0; i < continuous.size() && same; ++i) {
        const double tolerance =
            case_root.empty() ? 0.0
                              : model.solver.terminal.eos *
                                    std::max({1.0, std::abs(continuous[i]),
                                              std::abs(recovered[i])});
        same = std::abs(continuous[i] - recovered[i]) <= tolerance;
        if (!same)
          std::cerr << "restart difference rank=" << rank << " entry=" << i
                    << " delta=" << continuous[i] - recovered[i]
                    << " tolerance=" << tolerance << '\n';
      }
      if ((!status || !same) && rank == 0)
        std::cerr << "continuation " << unsigned(status.code) << ":"
                  << status.detail << " equal=" << same << "\n";
      okay = collective(bool(status) && same);
    }
    // One failing rank must reject every rank and preserve every accepted
    // primitive, previous layer, rate history and mass flux, not just Y_A.
    gas.fail = rank == 0;
    if (okay && case_root.empty()) {
      const auto accepted = physical_values(snap);
      status = driver.advance({1, 1, 1, 1, 1}, report);
      okay = collective(!status && !report.accepted);
      if (okay)
        okay = collective(bool(driver.committed_restart_snapshot(snap)) &&
                          snap.step == 2 && accepted == physical_values(snap));
    }
    if (okay && case_root.empty()) {
      RestartSnapshot before;
      status = restored.committed_restart_snapshot(before);
      const auto accepted = physical_values(before);
      restored_gas.redirect_output = rank == 0;
      if (status) status = restored.advance({1, 1, 1, 1, 1}, report);
      okay = collective(!status && !report.accepted);
      if (okay) okay = collective(bool(restored.committed_restart_snapshot(before)) &&
                                 accepted == physical_values(before));
    }
    if (okay && external_esf) {
      const auto accepted = physical_values(snap);
      esf_gas.fail_interval = rank == 0;
      status = driver.advance({1, 1, 1, 1, 1}, report);
      okay = collective(!status && !report.accepted);
      if (okay)
        okay = collective(bool(driver.committed_restart_snapshot(snap)) &&
                          accepted == physical_values(snap));
      if (!okay && rank == 0)
        std::cerr << "interval failure rollback failed\n";
      esf_gas.fail_interval = false;
      if (okay && tcr) {
        // Fail only after every cell has completed chemistry and staged its
        // TCR history, at the final ensemble reconciliation query.
        const auto cells = snap.patch.cells;
        esf_gas.fail_query_after_intervals =
            rank == 0 ? esf_gas.advance_calls +
                            model.reaction.esf->fields * std::uint64_t(cells.x) * cells.y * cells.z
                      : 0;
        status = driver.advance({1, 1, 1, 1, 1}, report);
        okay = collective(!status && !report.accepted &&
                          esf_gas.advance_calls >=
                              esf_gas.fail_query_after_intervals);
        if (okay)
          okay = collective(bool(driver.committed_restart_snapshot(snap)) &&
                            accepted == physical_values(snap));
        esf_gas.fail_query_after_intervals = 0;
        if (!okay && rank == 0)
          std::cerr << "late TCR transaction rollback failed\n";
      }
      if (okay) {
        status = driver.advance({1, 1, 1, 1, 1}, report);
        if (status)
          status = restored.advance({1, 1, 1, 1, 1}, report);
        RestartSnapshot retried, continuous;
        if (status)
          status = driver.committed_restart_snapshot(retried);
        if (status)
          status = restored.committed_restart_snapshot(continuous);
        const auto a = physical_values(
            retried, initial.pressure_reference *
                         model.thermophysics.species[0].molecular_weight /
                         (kUniversalGasConstant * initial.temperature));
        const auto b = physical_values(
            continuous, initial.pressure_reference *
                            model.thermophysics.species[0].molecular_weight /
                            (kUniversalGasConstant * initial.temperature));
        bool same = a.size() == b.size();
        for (std::size_t i = 0; i < a.size() && same; ++i)
          same = std::abs(a[i] - b[i]) <=
                 model.solver.terminal.eos *
                     std::max({1., std::abs(a[i]), std::abs(b[i])});
        if ((!status || !same) && rank == 0) {
          std::cerr << "retry comparison failed " << unsigned(status.code)
                    << ":" << status.detail << "\n";
          for (std::size_t i = 0; i < std::min(a.size(), b.size()); ++i)
            if (std::abs(a[i] - b[i]) >
                model.solver.terminal.eos *
                    std::max({1., std::abs(a[i]), std::abs(b[i])})) {
              std::cerr << "first mismatch " << i << " of " << a.size()
                        << " values " << a[i] << " " << b[i] << "\n";
              break;
            }
        }
        okay = collective(bool(status) && same);
      }
    }
    if (okay && tcr) {
      RestartSnapshot before, after;
      status = driver.committed_restart_snapshot(before);
      const auto saved = physical_values(before);
      const auto records = before.cell_records.values;
      const std::vector<std::uint8_t> saved_history(
          records.data, records.data + records.size);
      esf_gas.zero_progress = true;
      if (status)
        status = driver.advance({1, 1, 1, 1, 1}, report);
      const bool shadow = model.reaction.esf->tcr.mode == TcrMode::shadow;
      okay = collective(shadow ? bool(status) && report.accepted
                               : !status && !report.accepted);
      if (okay)
        okay = collective(bool(driver.committed_restart_snapshot(after)));
      if (okay && !shadow)
        okay = collective(physical_values(after) == saved);
      if (okay && shadow) {
        okay &= after.step == before.step + 1;
        for (std::size_t i = 0; i < saved_history.size(); ++i)
          if (i % 120 >= 20)
            okay &= saved_history[i] == after.cell_records.values.data[i];
        okay = collective(okay);
      }
      esf_gas.zero_progress = false;
      if (!okay && rank == 0)
        std::cerr << "TCR weak denominator mode contract failed\n";
    }
    if (okay && tcr) {
      EsfGas counter_gas;
      CompiledCasePlan counter_plan;
      ProductDriver counter_driver;
      status = ProductCompiler::compile(MPI_COMM_WORLD, model, case_root,
                                        counter_plan, counter_gas.bindings());
      if (status)
        status = ProductDriver::create(MPI_COMM_WORLD, std::move(counter_plan),
                                       counter_driver);
      RestartExpected expected;
      if (status)
        status = counter_driver.restart_expected(expected);
      RestartImage image;
      if (status)
        status = RestartReader::load(MPI_COMM_WORLD, root, expected, image);
      constexpr std::uint64_t large_revision = UINT64_C(9007199254740993);
      if (status) {
        for (std::size_t cell = 0; cell < image.cell_records.size();
             cell += 120)
          for (unsigned b = 0; b < 8; ++b)
            image.cell_records[cell + 8 + b] =
                std::uint8_t(large_revision >> (8 * b));
        status = counter_driver.initialize_restart(image);
      }
      if (status)
        status = counter_driver.advance({1, 1, 1, 1, 1}, report);
      RestartSnapshot result;
      if (status)
        status = counter_driver.committed_restart_snapshot(result);
      okay = collective(bool(status));
      if (okay) {
        for (std::size_t cell = 0; cell < result.cell_records.values.size;
             cell += 120) {
          std::uint64_t value = 0;
          for (unsigned b = 0; b < 8; ++b)
            value |=
                std::uint64_t(result.cell_records.values.data[cell + 8 + b])
                << (8 * b);
          okay &= value == large_revision + 1;
        }
        okay = collective(okay);
      }
      if (!okay && rank == 0)
        std::cerr << "TCR integer history continuation failed\n";
    }
    if (okay && tcr && external_esf) {
      okay = native_fold_continuation(model, case_root, root / "fold", initial,
                                      rank);
      if (!okay && rank == 0)
        std::cerr << "native TCR fold continuation failed\n";
    }
    if (okay && esf && !tcr) {
      // Restart carries a zero-mean sinusoidal composition disturbance across
      // patch boundaries. The discrete Fourier eigenvalue is independent of
      // the production Halo and diffusion kernels.
      EsfGas spatial_gas;
      CompiledCasePlan spatial_plan;
      ProductDriver spatial;
      status = ProductCompiler::compile(
          MPI_COMM_WORLD, model, case_root, spatial_plan,
          external_esf ? spatial_gas.bindings() : ProductCouplingBindings{});
      if (status)
        status = ProductDriver::create(MPI_COMM_WORLD, std::move(spatial_plan),
                                       spatial);
      RestartExpected expected;
      if (status)
        status = spatial.restart_expected(expected);
      RestartImage image;
      if (status)
        status = RestartReader::load(MPI_COMM_WORLD, root, expected, image);
      std::vector<double> variance;
      if (status) {
        const auto cells = image.patch.cells;
        const auto count = std::size_t(cells.x) * cells.y * cells.z;
        variance.resize(count);
        // Manufacture a constant enthalpy coordinate explicitly. The
        // analytic Fourier oracle below perturbs composition alone; tiny
        // spatial h differences from the earlier pressure solve have their
        // own bounds and would also participate in the shared noise limiter.
        double uniform_h{};
        for(const auto& field:image.fields)
          if(field.role==RestartFieldRole::enthalpy) uniform_h=field.values[0];
        MPI_Bcast(&uniform_h,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
        for(auto& field:image.fields) {
          if(field.role==RestartFieldRole::enthalpy)
            std::fill(field.values.begin(),field.values.end(),uniform_h);
          if(field.role==RestartFieldRole::stochastic_field)
            for(std::size_t c=0;c<count;++c) field.values[3*c+2]=uniform_h;
        }
        const double pi = std::acos(-1.0), amplitude = 0.003;
        const double dx =
            (model.mesh.upper.x - model.mesh.lower.x) / image.global_cells.x;
        const double dy =
            (model.mesh.upper.y - model.mesh.lower.y) / image.global_cells.y;
        const double dz =
            (model.mesh.upper.z - model.mesh.lower.z) / image.global_cells.z;
        std::array<std::array<double,3>,4> increments{};
        const auto wiener=esf::detail::balanced_wiener(4,image.dt,
            {model.reaction.esf->seed,image.step,1,0,0,1},increments.data(),increments.size());
        if (wiener!=portable::Status::success) return 1;
        const double maximum_perturbation=amplitude*std::cos(pi/image.global_cells.x);
        const double rho = initial.pressure_reference *
                           model.thermophysics.species[0].molecular_weight /
                           (kUniversalGasConstant * initial.temperature);
        const std::size_t start = 3 + model.transported_scalars.size();
        // Fix the manufactured mean, material and per-field base globally.
        // This isolates bounded stochastic forcing followed by one joint
        // diffusion/IEM solve; each field retains the same accepted history.
        std::array<double,4> base{};
        for(std::size_t f=0;f<4;++f)base[f]=image.fields[start+f].values[0];
        MPI_Bcast(base.data(),4,MPI_DOUBLE,0,MPI_COMM_WORLD);
        double mean=0;for(double q:base)mean+=q/4;
        std::array<double,4> cache{};
        const auto material=std::find_if(image.fields.begin(),image.fields.end(),[](const auto& field) {
          return field.role==RestartFieldRole::stochastic_transport;
        });
        std::copy_n(material->values.begin(),4,cache.begin());
        MPI_Bcast(cache.data(),4,MPI_DOUBLE,0,MPI_COMM_WORLD);
        const double mix=image.dt*(cache[2]+cache[3])/
            (model.reaction.mixing_c_z*rho*std::pow(dx*dy*dz,2./3.));
        const double diffusion=image.dt*cache[0]/(rho*dx*dx);
        const int nx=image.global_cells.x;
        std::vector<std::vector<double>> matrix(nx,std::vector<double>(nx));
        for(int x=0;x<nx;++x) {
          matrix[x][x]=1+2*diffusion+mix;
          matrix[x][(x+nx-1)%nx]-=diffusion;
          matrix[x][(x+1)%nx]-=diffusion;
        }
        std::array<std::vector<double>,4> advanced;
        for(std::size_t f=0;f<4;++f) {
          std::vector<double> rhs(nx);
          const double sign=f%2 ? -1 : 1;
          for(int x=0;x<nx;++x) {
            const double perturbation=amplitude*std::sin(2*pi*(x+.5)/nx);
            const double gradient=sign*amplitude*std::sin(2*pi/nx)/dx*
                std::cos(2*pi*(x+.5)/nx);
            const double raw_noise=std::sqrt(2*(cache[2]/.7+cache[3]/.5)/rho)*
                gradient*increments[f][0];
            const double noise=std::max(-maximum_perturbation-sign*perturbation,
                std::min(maximum_perturbation-sign*perturbation,raw_noise));
            rhs[x]=base[f]+sign*perturbation+noise+mix*mean;
          }
          advanced[f]=test::dense_solve(matrix,rhs);
        }
        for (int z=0;z<cells.z;++z) for(int y=0;y<cells.y;++y) for(int x=0;x<cells.x;++x) {
          const auto cell=std::size_t(x)+std::size_t(cells.x)*(y+std::size_t(cells.y)*z);
          const int gx=image.patch.begin.x+x;
          const double perturbation=amplitude*std::sin(2*pi*(gx+.5)/nx);
          double shifted_mean=0,moment=0;
          for(std::size_t f=0;f<4;++f) {
            auto& values=image.fields[start+f].values;
            values[3*cell]=base[f]+(f%2 ? -1 : 1)*perturbation;
            values[3*cell+1]=1-values[3*cell];
            shifted_mean+=advanced[f][gx]/4;
          }
          for(std::size_t f=0;f<4;++f)moment+=std::pow(advanced[f][gx]-shifted_mean,2)/4;
          variance[cell]=moment*std::exp(-4*image.dt);
          for(auto& field:image.fields) {
            if(field.role==RestartFieldRole::independent_species)field.values[cell]=mean;
            if(field.role==RestartFieldRole::stochastic_transport)
              std::copy(cache.begin(),cache.end(),field.values.begin()+field.components*cell);
          }
        }
        status = spatial.initialize_restart(image);
      }
      if (status)
        status = spatial.advance({1, 1, 1, 1, 1}, report);
      RestartSnapshot result;
      if (status)
        status = spatial.committed_restart_snapshot(result);
      bool matches = bool(status);
      unsigned noise_mismatches{};
      if (status) {
        const auto cells = image.patch.cells;
        const std::size_t start = 3 + model.transported_scalars.size();
        for (int z = 0; z < cells.z; ++z)
          for (int y = 0; y < cells.y; ++y)
            for (int x = 0; x < cells.x; ++x) {
              const auto cell =
                  std::size_t(x) +
                  std::size_t(cells.x) * (y + std::size_t(cells.y) * z);
              double mean = 0, moment = 0;
              for (std::size_t f = 0; f < 4; ++f)
                mean += result.fields.data[start + f].values.unchecked(
                            {x, y, z}, 0) /
                        4;
              for (std::size_t f = 0; f < 4; ++f)
                moment +=
                    std::pow(result.fields.data[start + f].values.unchecked(
                                 {x, y, z}, 0) -
                                 mean,
                             2) /
                    4;
              if(!(std::abs(moment-variance[cell])<3e-12) && noise_mismatches++<8 && rank==0)
                std::cerr<<"ESF noise variance x="<<x+image.patch.begin.x
                    <<" y="<<y<<" z="<<z<<" actual="<<std::setprecision(17)<<moment<<" expected="<<variance[cell]
                    <<" step="<<image.step<<" dt="<<image.dt<<'\n';
              matches &= std::abs(moment - variance[cell]) < 3e-12;
            }
      }
      okay = collective(matches);
      if (!okay && rank == 0)
        std::cerr << "ESF spatial continuation failed " << unsigned(status.code)
                  << ":" << status.detail << "\n";
    }
    if (rank == 0)
      std::cout << (okay ? "PASS" : "FAIL")
                << ": product reaction changes accepted species; provider "
                   "failure rolls back full step\n";
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0) {
    std::error_code error;
    std::filesystem::remove_all(root, error);
  }
  MPI_Finalize();
  return okay ? 0 : 1;
}
