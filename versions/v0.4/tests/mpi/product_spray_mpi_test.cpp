// SPDX-License-Identifier: Apache-2.0
#include "../../src/models_chemistry_adapter_detail.hpp"
#include "hundun/v04_app.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <iostream>
#include <mpi.h>
#include <unistd.h>
#include <vector>
using namespace hundun::v04;
namespace {
class Gas final : public portable::GasQueryProvider,
                  public portable::GasAdvanceProvider {
public:
  chemistry::detail::AnalyticIsomerBackend backend;
  std::uint64_t calls{}, fail_at{};
  bool fail_half{};
  unsigned half_calls{};
  portable::Status
  advance_gas(const portable::GasAdvanceQuery &q,
              portable::GasAdvanceOutput &out) noexcept override {
    if (fail_half && ++half_calls == 2)
      return portable::Status::provider_failure;
    return backend.advance_gas(q, out);
  }
  ProductCouplingBindings bindings() noexcept {
    return {this, this, &backend.closure_identity()};
  }
  const portable::GasIdentity &gas_identity() const noexcept override {
    return backend.gas_identity();
  }
  portable::Status query_gas(const portable::GasQuery &q,
                             portable::GasQueryOutput &out) noexcept override {
    ++calls;
    if (fail_at && calls >= fail_at)
      return portable::Status::provider_failure;
    return backend.query_gas(q, out);
  }
};
bool agree(bool local) {
  int a = local, b = 0;
  MPI_Allreduce(&a, &b, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  return b;
}
std::uint64_t u64(const std::uint8_t *p, unsigned n = 8) {
  std::uint64_t v = 0;
  for (unsigned i = 0; i < n; ++i)
    v |= std::uint64_t(p[i]) << (8 * i);
  return v;
}
double real(const std::uint8_t *p) {
  auto b = u64(p);
  double v;
  std::memcpy(&v, &b, 8);
  return v;
}
std::vector<double> snapshot(const RestartSnapshot &s,
                             bool normalized = false) {
  std::vector<double> v{s.time,
                        s.dt,
                        s.pressure_reference,
                        s.previous_pressure_reference,
                        s.closed_mass_target,
                        double(s.step)};
  for (auto fields : {s.fields, s.previous_fields, s.accepted_rate_fields,
                      s.previous_rate_fields})
    for (std::size_t f = 0; f < fields.size; ++f) {
      auto q = fields.data[f].values;
      double scale = 1.;
      if (normalized) {
        const double pressure = 101325.,
                     density = pressure * 28 / (kUniversalGasConstant * 400);
        switch (fields.data[f].role) {
        case RestartFieldRole::velocity:
          scale = std::sqrt(pressure / density);
          break;
        case RestartFieldRole::pressure_perturbation:
        case RestartFieldRole::pressure_absolute:
          scale = pressure;
          break;
        case RestartFieldRole::enthalpy:
          scale = pressure / density;
          break;
        case RestartFieldRole::enthalpy_nonadvective_rate:
          scale = pressure / s.dt;
          break;
        case RestartFieldRole::scalar_nonadvective_rate:
          scale = density / s.dt;
          break;
        default:
          break;
        }
      }
      for (int z = 0; z < q.interior.z; ++z)
        for (int y = 0; y < q.interior.y; ++y)
          for (int x = 0; x < q.interior.x; ++x)
            for (unsigned c = 0; c < q.components; ++c)
              v.push_back(q.unchecked({x, y, z}, c) / scale);
    }
  for (auto flux : {s.final_mass_flux, s.previous_mass_flux})
    for (auto q : {flux.x, flux.y, flux.z})
      for (int z = 0; z < q.extents.z; ++z)
        for (int y = 0; y < q.extents.y; ++y)
          for (int x = 0; x < q.extents.x; ++x)
            v.push_back(q.unchecked({x, y, z}));
  std::size_t offset = 0;
  for (std::size_t cell = 0; cell < s.cell_records.variable_cell_bytes.size;
       ++cell) {
    const auto length = s.cell_records.variable_cell_bytes.data[cell];
    const auto *row = s.cell_records.values.data + offset;
    offset += length;
    const auto tcr = length >= 24 ? u64(row + 16, 4) : 0;
    for (std::size_t i = 0; i < length;) {
      const auto history_offset = i >= 24 ? i - 24 : SIZE_MAX;
      if (normalized && tcr == 120 &&
          (history_offset == 24 || history_offset == 32 ||
           history_offset == 40 || history_offset == 48 ||
           history_offset == 80 || history_offset == 88)) {
        v.push_back(real(row + i));
        i += 8;
      } else
        v.push_back(row[i++]);
    }
  }
  for (std::size_t i = 0; i < s.cell_records.variable_cell_bytes.size; ++i)
    v.push_back(s.cell_records.variable_cell_bytes.data[i]);
  return v;
}
bool inventory(const RestartSnapshot &s, int rank, bool wall, bool ibm) {
  ConstFieldView h, p, u, Y, passive;
  for (std::size_t i = 0; i < s.fields.size; ++i) {
    auto f = s.fields.data[i];
    switch (f.role) {
    case RestartFieldRole::enthalpy:
      h = f.values;
      break;
    case RestartFieldRole::pressure_perturbation:
      p = f.values;
      break;
    case RestartFieldRole::velocity:
      u = f.values;
      break;
    case RestartFieldRole::independent_species:
      Y = f.values;
      break;
    case RestartFieldRole::transported_scalar:
      passive = f.values;
      break;
    default:
      break;
    }
  }
  // Independent constant-cp EOS from the fixture: equal MW=28, h_A-h_B=1e5.
  std::array<double, 8> local{}, global{};
  const double V =
      1. / (s.global_cells.x * s.global_cells.y * s.global_cells.z);
  for (int z = 0; z < h.interior.z; ++z)
    for (int y = 0; y < h.interior.y; ++y)
      for (int x = 0; x < h.interior.x; ++x) {
        Int3 c{x, y, z};
        const double px = (x + s.patch.begin.x + .5) / s.global_cells.x,
                     py = (y + s.patch.begin.y + .5) / s.global_cells.y,
                     pz = (z + s.patch.begin.z + .5) / s.global_cells.z;
        if (ibm && px > .375 && px < .625 && py > .375 && py < .625 &&
            pz > .375 && pz < .625)
          continue;
        const double H = h.unchecked(c, 0),
                     P = s.pressure_reference + p.unchecked(c, 0);
        const double T = 298.15 + (H - 1e5 * Y.unchecked(c, 0)) / 1000,
                     rho = P * 28 / (kUniversalGasConstant * T);
        double K = 0;
        for (unsigned d = 0; d < 3; ++d) {
          const auto speed = u.unchecked(c, d);
          K += .5 * speed * speed;
          local[d + 2] += rho * V * speed;
        }
        local[0] += rho * V;
        local[1] += (rho * (H + K) - P) * V;
        local[5] += rho * V;
        if (passive.base)
          local[7] += rho * V * passive.unchecked(c, 0);
      }
  bool ensemble_ok = true;
  std::size_t fields = 0;
  for (std::size_t f = 0; f < s.fields.size; ++f)
    if (s.fields.data[f].role == RestartFieldRole::stochastic_field)
      ++fields;
  if (fields)
    for (int z = 0; z < h.interior.z; ++z)
      for (int y = 0; y < h.interior.y; ++y)
        for (int x = 0; x < h.interior.x; ++x) {
          const Int3 c{x, y, z};
          double mean_y = 0, mean_h = 0;
          for (std::size_t f = 0; f < s.fields.size; ++f)
            if (s.fields.data[f].role == RestartFieldRole::stochastic_field) {
              const auto field = s.fields.data[f].values;
              mean_y += field.unchecked(c, 0) / fields;
              mean_h += field.unchecked(c, 2) / fields;
            }
          ensemble_ok &= std::abs(mean_y - Y.unchecked(c, 0)) < 2e-12 &&
                         std::abs(mean_h - h.unchecked(c, 0)) <
                             2e-12 * std::max(1., std::abs(h.unchecked(c, 0)));
        }
  if (!agree(ensemble_ok))
    return false;
  auto records = s.cell_records;
  std::size_t offset = 0;
  bool valid = true;
  for (std::size_t c = 0; c < records.variable_cell_bytes.size; ++c) {
    const auto bytes = records.variable_cell_bytes.data[c];
    if (!bytes)
      continue;
    const auto *base = records.values.data + offset;
    offset += bytes;
    const auto count = u64(base + 8, 4), injectors = u64(base + 12, 4),
               tcr = u64(base + 16, 4);
    const auto *row = base + 24 + tcr;
    for (std::uint64_t i = 0; i < count; ++i, row += 144) {
      const double mass = real(row + 64) * real(row + 80), T = real(row + 88);
      double K = 0;
      for (unsigned d = 0; d < 3; ++d) {
        const auto speed = real(row + 40 + 8 * d);
        K += .5 * speed * speed;
        local[d + 2] += mass * speed;
      }
      local[0] += mass;
      local[1] += mass * (-1e5 + 1000 * (T - 298.15) + K);
      local[6] += 1;
      if (ibm)
        valid &= real(row + 16) < .375 && real(row + 16) > .3748 &&
                 real(row + 40) < 0;
      if (wall)
        valid &=
            real(row + 16) < 1 && real(row + 16) > .9998 && real(row + 40) < 0;
    }
    for (std::uint64_t i = 0; i < injectors; ++i, row += 24)
      if (u64(row) != UINT64_C(9007199254740997) || u64(row + 16) != s.step ||
          real(row + 8) != 0)
        valid = false;
  }
  if (!agree(valid))
    return false;
  MPI_Allreduce(local.data(), global.data(), 8, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  const double initial_mass =
      101325 * 28 / (kUniversalGasConstant * 400) * (ibm ? 63. / 64 : 1);
  const double injected = 1e-4 * s.step,
               expected_energy = initial_mass * 102850 -
                                 101325 * (ibm ? 63. / 64 : 1) +
                                 injected * (-1e5 + 2);
  const double dm = global[0] - initial_mass - injected,
               dE = global[1] - expected_energy;
  if (rank == 0)
    std::cerr.precision(17),
        std::cerr << "step=" << s.step << " mass_error=" << dm
                  << " energy_error=" << dE
                  << " momentum_error=" << global[2] - 2 * injected
                  << " gas_target_error=" << global[5] - s.closed_mass_target
                  << " parcels=" << global[6] << '\n';
  return std::abs(dm) < 1e-10 && std::abs(dE) < 1e-5 &&
         (wall || ibm || std::abs(global[2] - 2 * injected) < 1e-9) &&
         std::abs(global[3]) < 1e-9 && std::abs(global[4]) < 1e-9 &&
         std::abs(global[5] - s.closed_mass_target) < 1e-10 &&
         global[6] == s.step &&
         (!passive.base || std::abs(global[7] + 2 * initial_mass) < 1e-10);
}
} // namespace
int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  bool ok = argc == 2;
  {
    ValidatedModel model;
    Status status;
    if (ok)
      status = CaseCompiler::load_and_compile(MPI_COMM_WORLD, argv[1], model);
    Gas reference_gas, retry_gas;
    ProductDriver reference, retry;
    auto create = [&](Gas &gas, ProductDriver &driver) {
      CompiledCasePlan plan;
      auto r = ProductCompiler::compile(MPI_COMM_WORLD, model, argv[1], plan,
                                        gas.bindings());
      if (r)
        r = ProductDriver::create(MPI_COMM_WORLD, std::move(plan), driver);
      std::array<double, 2> initial_scalars{.01, -2};
      DriverInitialState init;
      init.temperature = 400;
      init.transported_scalars = {initial_scalars.data(),
                                  model.transported_scalars.size()};
      if (r)
        r = driver.initialize(init);
      return r;
    };
    if (status)
      status = create(reference_gas, reference);
    if (status)
      status = create(retry_gas, retry);
    ok = agree(ok && bool(status));
    DriverStepReport report;
    if (ok) {
      status = reference.advance({1, 1, 1, 1, 1}, report);
      ok = agree(bool(status) && report.accepted);
    }
    if (ok) {
      status = retry.advance({1, 1, 1, 1, 1}, report);
      ok = agree(bool(status) && report.accepted);
    }
    if (ok && model.reaction.esf) {
      RestartSnapshot first;
      status = reference.committed_restart_snapshot(first);
      ok = agree(bool(status));
      if (ok) {
        std::size_t offset = 0;
        int enriched = 0;
        bool valid = true;
        for (std::size_t cell = 0;
             cell < first.cell_records.variable_cell_bytes.size; ++cell) {
          const auto bytes = first.cell_records.variable_cell_bytes.data[cell];
          const auto *row = first.cell_records.values.data + offset;
          offset += bytes;
          if (bytes < 144 || u64(row + 16, 4) != 120) {
            valid = false;
            continue;
          }
          const auto eta = real(row + 24 + 24);
          valid &= eta >= .01 - 1e-13;
          if (eta > .01 + 1e-13)
            ++enriched;
        }
        int total = 0;
        MPI_Allreduce(&enriched, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        // Initial gas/field mean is exactly .01 everywhere. The actual
        // trilinear deposition enriches eight interior or four wall-adjacent
        // cells before TCR observes eta.
        const int expected =
            model.boundaries[0].flow_kind == BoundaryKind::slip ? 4 : 8;
        ok = agree(valid && total == expected);
        if (!ok && rank == 0)
          std::cerr << "TCR post-source enriched cells=" << total << '\n';
      }
    }
    if (ok) {
      reference_gas.calls = 0;
      status = reference.advance({1, 1, 1, 1, 1}, report);
      ok = agree(bool(status) && report.accepted);
      if (ok) {
        RestartSnapshot s;
        status = reference.committed_restart_snapshot(s);
        ok = agree(bool(status));
        if (ok)
          ok = agree(inventory(
              s, rank, model.boundaries[0].flow_kind == BoundaryKind::slip,
              model.immersed_boundary.has_value()));
      }
    }
    if (ok) {
      RestartSnapshot before;
      status = retry.committed_restart_snapshot(before);
      auto saved = snapshot(before);
      retry_gas.calls = 0;
      retry_gas.fail_at = rank == 0 ? reference_gas.calls - 1 : 0;
      status = retry.advance({1, 1, 1, 1, 1}, report);
      ok = agree(!status && !report.accepted);
      if (rank == 0)
        std::cerr << "late failure " << unsigned(status.code) << ':'
                  << status.detail << " query=" << retry_gas.calls << '/'
                  << reference_gas.calls << '\n';
      if (ok) {
        status = retry.committed_restart_snapshot(before);
        auto current = snapshot(before);
        bool same = current == saved;
        if (!same) {
          for (std::size_t i = 0; i < std::min(saved.size(), current.size());
               ++i)
            if (saved[i] != current[i]) {
              std::cerr << "rollback rank=" << rank << " index=" << i
                        << " before=" << saved[i] << " after=" << current[i]
                        << '\n';
              break;
            }
        }
        ok = agree(bool(status) && same);
      }
      retry_gas.fail_at = 0;
    }
    if (ok) {
      status = retry.advance({1, 1, 1, 1, 1}, report);
      ok = agree(bool(status) && report.accepted);
    }
    if (ok) {
      RestartSnapshot a, b;
      status = reference.committed_restart_snapshot(a);
      if (status)
        status = retry.committed_restart_snapshot(b);
      auto av = snapshot(a), bv = snapshot(b);
      bool same = av.size() == bv.size();
      // A rejected pressure solve invalidates its warm start. Accepted state
      // rollback is byte exact above; a new solve is compared at solver
      // accuracy.
      for (std::size_t i = 0; i < std::min(av.size(), bv.size()); ++i)
        same &= std::abs(av[i] - bv[i]) <=
                2e-12 * std::max({1., std::abs(av[i]), std::abs(bv[i])});
      same &=
          a.cell_records.values.size == b.cell_records.values.size &&
          std::memcmp(a.cell_records.values.data, b.cell_records.values.data,
                      a.cell_records.values.size) == 0;
      if (!same) {
        for (std::size_t i = 0; i < std::min(av.size(), bv.size()); ++i)
          if (av[i] != bv[i]) {
            std::cerr.precision(17);
            std::cerr << "retry rank=" << rank << " index=" << i
                      << " control=" << av[i] << " retry=" << bv[i] << '\n';
            break;
          }
      }
      ok = agree(bool(status) && same);
    }
    if (ok) {
      status = reference.advance({1, 1, 1, 1, 1}, report);
      ok = agree(bool(status) && report.accepted);
    }
    if (ok) {
      RestartSnapshot s;
      status = reference.committed_restart_snapshot(s);
      ok = agree(bool(status));
      if (ok)
        ok = agree(inventory(
            s, rank, model.boundaries[0].flow_kind == BoundaryKind::slip,
            model.immersed_boundary.has_value()));
    }
    if (ok) {
      int pid = int(getpid());
      MPI_Bcast(&pid, 1, MPI_INT, 0, MPI_COMM_WORLD);
      const auto root = std::filesystem::temp_directory_path() /
                        ("hundun-native-spray-" + std::to_string(pid));
      RestartSnapshot accepted;
      status = reference.committed_restart_snapshot(accepted);
      if (status)
        status = RestartWriter::write(MPI_COMM_WORLD, root, accepted);
      Gas restored_gas;
      ProductDriver restored;
      CompiledCasePlan restored_plan;
      if (status)
        status =
            ProductCompiler::compile(MPI_COMM_WORLD, model, argv[1],
                                     restored_plan, restored_gas.bindings());
      if (status)
        status = ProductDriver::create(MPI_COMM_WORLD, std::move(restored_plan),
                                       restored);
      RestartExpected expected;
      RestartImage image;
      if (status)
        status = restored.restart_expected(expected);
      if (status)
        status = RestartReader::load(MPI_COMM_WORLD, root, expected, image);
      ok = agree(bool(status));
      if (ok) {
        auto intact = image.cell_records;
        std::size_t offset = 0;
        bool changed = false;
        for (auto bytes : image.cell_record_lengths) {
          if (bytes && u64(image.cell_records.data() + offset + 8, 4) > 0) {
            // Change only diameter, retaining positive/finite parcel values.
            // The native material/EOS restore preflight must reject this.
            const auto index = offset + 24 +
                               u64(image.cell_records.data() + offset + 16, 4) +
                               72;
            const auto d = 2 * real(image.cell_records.data() + index);
            std::uint64_t bits;
            std::memcpy(&bits, &d, 8);
            for (unsigned b = 0; b < 8; ++b)
              image.cell_records[index + b] = (bits >> (8 * b)) & 255;
            changed = true;
            break;
          }
          offset += bytes;
        }
        int a = changed, b = 0;
        MPI_Allreduce(&a, &b, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        auto rejected = restored.initialize_restart(image);
        ok = agree(b > 0 && !rejected);
        image.cell_records = std::move(intact);
        if (ok)
          status = restored.initialize_restart(image);
        RestartSnapshot actual;
        if (ok && status)
          status = restored.committed_restart_snapshot(actual);
        if (ok)
          ok = agree(bool(status) && snapshot(actual) == snapshot(accepted));
      }
      if (ok && model.reaction.esf) {
        const auto saved = snapshot(accepted);
        reference_gas.half_calls = 0;
        reference_gas.fail_half = rank == 0;
        auto failure = reference.advance({1, 1, 1, 1, 1}, report);
        ok = agree(!failure && !report.accepted);
        if (rank == 0)
          std::cerr << "second half status=" << unsigned(failure.code) << ':'
                    << failure.detail << " accepted=" << report.accepted
                    << " halves=" << reference_gas.half_calls << '\n';
        reference_gas.fail_half = false;
        RestartSnapshot after;
        if (ok)
          status = reference.committed_restart_snapshot(after);
        if (ok)
          ok = agree(bool(status) && snapshot(after) == saved);
        if (ok)
          status = reference.advance({1, 1, 1, 1, 1}, report);
        if (ok)
          ok = agree(bool(status) && report.accepted);
        if (ok)
          status = restored.advance({1, 1, 1, 1, 1}, report);
        if (ok)
          ok = agree(bool(status) && report.accepted);
        RestartSnapshot control;
        if (ok)
          status = reference.committed_restart_snapshot(after);
        if (ok && status)
          status = restored.committed_restart_snapshot(control);
        if (ok) {
          // Restart drops numerical warm starts. Compare the same dimensionless
          // SI scales as the native reacting restart oracle. TCR numerical
          // history uses doubles; parcel values and integer records stay exact.
          const auto a = snapshot(after, true), b = snapshot(control, true);
          bool same = a.size() == b.size();
          for (std::size_t i = 0; i < std::min(a.size(), b.size()); ++i)
            same &= std::abs(a[i] - b[i]) <=
                    2e-12 * std::max({1., std::abs(a[i]), std::abs(b[i])});
          if (!same)
            for (std::size_t i = 0; i < std::min(a.size(), b.size()); ++i)
              if (std::abs(a[i] - b[i]) >
                  2e-12 * std::max({1., std::abs(a[i]), std::abs(b[i])})) {
                std::cerr.precision(17);
                std::cerr << "second half retry rank=" << rank << " index=" << i
                          << " after=" << a[i] << " control=" << b[i] << '\n';
                break;
              }
          ok = agree(bool(status) && same);
        }
        if (!ok && rank == 0)
          std::cerr << "P8 second-half rejection/retry diverged\n";
      }
      MPI_Barrier(MPI_COMM_WORLD);
      if (rank == 0)
        std::filesystem::remove_all(root);
    }
    if (!ok && rank == 0)
      for (const auto &failure :
           report.pressure_energy_globalization.candidate_evaluation_status)
        if (!failure)
          std::cerr << "candidate rejection " << unsigned(failure.code) << ':'
                    << failure.detail << '\n';
    if (!ok && rank == 0)
      std::cerr << "native spray product failure " << unsigned(status.code)
                << ':' << status.detail << '\n';
  }
  MPI_Finalize();
  return ok ? 0 : 1;
}
