// SPDX-License-Identifier: Apache-2.0
#include "../support/product_fixture.hpp"
#include "hundun/v04_app.hpp"
#include "hundun/v04_portable.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
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
    CompiledCasePlan plan;
    if (status)
      status = ProductCompiler::compile(MPI_COMM_WORLD, model, case_root, plan,
                                        case_root.empty()
                                            ? ProductCouplingBindings{&gas}
                                            : ProductCouplingBindings{});
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
    // Compare a real disk restart with uninterrupted BDF2 continuation.
    IsomerGas restored_gas;
    ProductDriver restored;
    if (okay) {
      status = RestartWriter::write(MPI_COMM_WORLD, root, snap);
      CompiledCasePlan restored_plan;
      if (status)
        status = ProductCompiler::compile(
            MPI_COMM_WORLD, model, case_root, restored_plan,
            case_root.empty() ? ProductCouplingBindings{&restored_gas}
                              : ProductCouplingBindings{});
      if (status)
        status = ProductDriver::create(MPI_COMM_WORLD, std::move(restored_plan),
                                       restored);
      RestartExpected expected;
      if (status)
        status = restored.restart_expected(expected);
      RestartImage image;
      if (status)
        status = RestartReader::load(MPI_COMM_WORLD, root, expected, image);
      if (status)
        status = restored.initialize_restart(image);
      if (status) {
        RestartSnapshot check;
        status = restored.committed_restart_snapshot(check);
        okay = collective(bool(status) &&
                          physical_values(check) == physical_values(snap));
      }
      okay = collective(okay && bool(status));
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
