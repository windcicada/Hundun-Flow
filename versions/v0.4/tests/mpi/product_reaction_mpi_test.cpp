// SPDX-License-Identifier: Apache-2.0
#include "../../src/models_chemistry_adapter_detail.hpp"
#include "../support/product_fixture.hpp"
#include "hundun/v04_app.hpp"
#include "hundun/v04_portable.hpp"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <mpi.h>
#include <unistd.h>
#include <vector>

using namespace hundun::v04;
namespace {
// Independent A -> B, equal molecular weights/cp/formation enthalpy. The
// source is -2*rho*Y_A; it cannot affect total gas mass or total enthalpy.
class IsomerGas final : public portable::GasQueryProvider {
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
    identity_.closure_fingerprint = 102;
  }
  const portable::GasIdentity &gas_identity() const noexcept override {
    return identity_;
  }
  portable::Status query_gas(const portable::GasQuery &q,
                             portable::GasQueryOutput &out) noexcept override {
    if (fail || q.species_count != 2 || out.capacity != 2)
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
    out.net_mass_rates_kg_per_m3_s[0] = -2 * rho * q.mass_fractions[0];
    out.net_mass_rates_kg_per_m3_s[1] = -out.net_mass_rates_kg_per_m3_s[0];
    if (redirect_output) out.net_mass_rates_kg_per_m3_s = nullptr;
    return portable::Status::success;
  }
  bool fail{};
  bool redirect_output{};

private:
  portable::GasIdentity identity_;
};
class EsfGas final : public portable::GasQueryProvider,
                     public portable::GasAdvanceProvider {
public:
  chemistry::detail::AnalyticIsomerBackend backend;
  bool fail_half{}, zero_progress{};
  unsigned half_calls{};
  std::uint64_t advance_calls{}, fail_query_after_halves{};
  const portable::GasIdentity &gas_identity() const noexcept override {
    return backend.gas_identity();
  }
  portable::Status query_gas(const portable::GasQuery &q,
                             portable::GasQueryOutput &out) noexcept override {
    if (fail_query_after_halves != 0 &&
        advance_calls >= fail_query_after_halves)
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
    if (fail_half && ++half_calls == 3)
      return portable::Status::provider_failure;
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
    Status status;
    const std::filesystem::path case_root = argc > 1 ? argv[1] : "";
    if (!case_root.empty())
      status = CaseCompiler::load_and_compile(MPI_COMM_WORLD, case_root, model);
    EsfGas esf_gas;
    const bool esf = model.reaction.mode == ReactionMode::esf_tpdf;
    const bool tcr = esf && model.reaction.esf->tcr.mode != TcrMode::off;
    CompiledCasePlan plan;
    if (status)
      status = ProductCompiler::compile(
          MPI_COMM_WORLD, model, case_root, plan,
          esf ? esf_gas.bindings()
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
      const double gamma = 1e-5 / 0.7; // Synthetic fixture: mu(300 K)/Pr.
      const double tau = model.reaction.mixing_c_z *
                         std::pow(std::cbrt(volume), 2) / (2 * gamma / rho);
      const double control =
          tcr && model.reaction.esf->tcr.mode == TcrMode::experimental
              ? 1.0 / 3.0
              : 1.0;
      const double dt = snap.dt,
                   relaxation = std::exp(-std::cbrt(control) * dt / (2 * tau)),
                   half = std::exp(-dt);
      double density_delta = 0, expected_variance = 0;
      for (double offset : model.reaction.esf->initial_species_offsets) {
        const double y0 = 0.25 + offset * relaxation, ym = y0 * half,
                     ye = ym * half;
        const double t0 = initial.temperature - 100 * offset * relaxation;
        const double tm = t0 + 100 * (y0 - ym);
        density_delta += (initial.pressure_reference / (R * t0) * (ym - y0) +
                          initial.pressure_reference / (R * tm) * (ye - ym)) /
                         4;
        expected_variance += std::pow(offset * relaxation * half * half, 2) / 4;
      }
      okay &= std::abs(accepted_y - (0.25 + density_delta / rho)) < 2e-12;
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
            esf ? restored_esf.bindings()
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
      // Restart preserves the stored accepted values exactly. Pressure
      // correction warm starts are deliberately not persisted; a reacting
      // pressure/enthalpy solve is equivalent within its terminal tolerance.
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
    if (okay && esf) {
      const auto accepted = physical_values(snap);
      esf_gas.fail_half = rank == 0;
      status = driver.advance({1, 1, 1, 1, 1}, report);
      okay = collective(!status && !report.accepted);
      if (okay)
        okay = collective(bool(driver.committed_restart_snapshot(snap)) &&
                          accepted == physical_values(snap));
      if (!okay && rank == 0)
        std::cerr << "half failure rollback failed\n";
      esf_gas.fail_half = false;
      if (okay && tcr) {
        // Fail only after every cell has completed chemistry and staged its
        // TCR history, at the final ensemble reconciliation query.
        const auto cells = snap.patch.cells;
        esf_gas.fail_query_after_halves =
            rank == 0 ? esf_gas.advance_calls +
                            8 * std::uint64_t(cells.x) * cells.y * cells.z
                      : 0;
        status = driver.advance({1, 1, 1, 1, 1}, report);
        okay = collective(!status && !report.accepted &&
                          esf_gas.advance_calls >=
                              esf_gas.fail_query_after_halves);
        if (okay)
          okay = collective(bool(driver.committed_restart_snapshot(snap)) &&
                            accepted == physical_values(snap));
        esf_gas.fail_query_after_halves = 0;
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
    if (okay && esf && !tcr) {
      // Restart carries a zero-mean sinusoidal composition disturbance across
      // patch boundaries. The discrete Fourier eigenvalue is independent of
      // the production Halo and diffusion kernels.
      EsfGas spatial_gas;
      CompiledCasePlan spatial_plan;
      ProductDriver spatial;
      status = ProductCompiler::compile(MPI_COMM_WORLD, model, case_root,
                                        spatial_plan, spatial_gas.bindings());
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
        const double pi = std::acos(-1.0), amplitude = 0.003;
        const double dx =
            (model.mesh.upper.x - model.mesh.lower.x) / image.global_cells.x;
        const double dy =
            (model.mesh.upper.y - model.mesh.lower.y) / image.global_cells.y;
        const double dz =
            (model.mesh.upper.z - model.mesh.lower.z) / image.global_cells.z;
        const double eigenvalue =
            4 * std::pow(std::sin(pi / image.global_cells.x), 2) / (dx * dx);
        const double rho = initial.pressure_reference *
                           model.thermophysics.species[0].molecular_weight /
                           (kUniversalGasConstant * initial.temperature);
        const std::size_t start = 3 + model.transported_scalars.size();
        for (int z = 0; z < cells.z; ++z)
          for (int y = 0; y < cells.y; ++y)
            for (int x = 0; x < cells.x; ++x) {
              const auto cell =
                  std::size_t(x) +
                  std::size_t(cells.x) * (y + std::size_t(cells.y) * z);
              const double gamma = image.fields.back().values[2 * cell];
              const double tau = model.reaction.mixing_c_z *
                                 std::pow(std::cbrt(dx * dy * dz), 2) /
                                 (2 * gamma / rho);
              const double perturbation =
                  amplitude *
                  std::sin(2 * pi * (image.patch.begin.x + x + 0.5) /
                           image.global_cells.x);
              double mean = 0, moment = 0;
              for (std::size_t f = 0; f < 4; ++f)
                mean += image.fields[start + f].values[3 * cell] / 4;
              for (std::size_t f = 0; f < 4; ++f) {
                const double sign = f % 2 ? -1 : 1;
                auto &values = image.fields[start + f].values;
                const double transported_deviation =
                    values[3 * cell] - mean +
                    sign * perturbation *
                        (1 - image.dt * gamma / rho * eigenvalue);
                moment += transported_deviation * transported_deviation / 4;
                values[3 * cell] += sign * perturbation;
                values[3 * cell + 1] -= sign * perturbation;
              }
              variance[cell] =
                  moment * std::exp(-image.dt / tau - 4 * image.dt);
            }
        status = spatial.initialize_restart(image);
      }
      if (status)
        status = spatial.advance({1, 1, 1, 1, 1}, report);
      RestartSnapshot result;
      if (status)
        status = spatial.committed_restart_snapshot(result);
      bool matches = bool(status);
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
