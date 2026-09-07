// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_app.hpp"
#include "../support/product_fixture.hpp"
#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string_view>
#include <vector>

using namespace hundun::v04;

bool agree(bool valid) {
  const int local = valid ? 1 : 0;
  int global = 0;
  return MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD) == MPI_SUCCESS && global != 0;
}

bool run(int rank) {
  auto model = test::product_model({16, 16, 16});
  model.mesh.lower = {-2.0, -2.0, -2.0};
  model.mesh.upper = {2.0, 2.0, 2.0};
  model.mesh.minimum_spacing = {0.25, 0.25, 0.25};
  model.time.control = TimeControlKind::fixed;
  model.time.initial_dt = model.time.minimum_dt = model.time.maximum_dt = 1.0e-5;
  model.immersed_boundary = ImmersedBoundarySpec{
      "cylinder_ascii.stl", ImmersedFluidSide::outside};
  CompiledCasePlan plan;
  Status status = ProductCompiler::compile(MPI_COMM_WORLD, model,
      std::filesystem::path{HUNDUN_V04_SOURCE_ROOT} / "tests" / "data", plan);
  ProductDriver driver;
  if (status) status = ProductDriver::create(MPI_COMM_WORLD, std::move(plan), driver);
  DriverInitialState initial;
  initial.temperature = 300.0;
  initial.velocity = {1.0, 0.0, 0.0};
  if (status) status = driver.initialize(initial);
  DriverCellTraceWindow window;
  window.cells = {Int3{11, 11, 8}, Int3{11, 8, 8}};
  window.count = 2U;
  window.first_step = 1U;
  window.last_step = 4U;
  if (status) status = driver.set_cell_trace_window(window);
  CommittedOutputSnapshot snapshot;
  if (status) status = driver.committed_output_snapshot(snapshot);
  if (!status) {
    std::cerr << "rank=" << rank << " placeholder setup="
              << unsigned(status.code) << '/' << status.detail << '\n';
    return false;
  }
  const Int3 cells = snapshot.patch.cells;
  const std::size_t count = std::size_t(cells.x) * cells.y * cells.z;
  if (!agree(snapshot.cell_activity.size == count)) return false;
  const auto enthalpy = [](const CommittedOutputSnapshot& output) {
    for (std::size_t i = 0U; i < output.fields.size; ++i)
      if (output.fields.data[i].stable_name == "h") return output.fields.data[i].values;
    return ConstFieldView{};
  };
  std::vector<double> initial_h(count);
  const auto h = enthalpy(snapshot);
  std::size_t flat = 0U;
  for (int z = 0; z < cells.z; ++z)
    for (int y = 0; y < cells.y; ++y)
      for (int x = 0; x < cells.x; ++x, ++flat)
        initial_h[flat] = h.unchecked({x, y, z}, 0U);
  double maximum_drift = 0.0;
  for (unsigned step_index = 0U; step_index < 4U; ++step_index) {
    DriverStepReport step;
    status = driver.advance({1, 1, 1, 1, 1}, step);
    for (std::size_t k = 0; k < step.cell_trace.count; ++k) {
      const auto& s = step.cell_trace.samples[k];
      if (std::getenv("HUNDUN_TEST_PRINT_CELL_TRACE") == nullptr) continue;
      std::cerr << std::setprecision(17) << "placeholder_trace step=" << s.step
                << " rank=" << rank << " cell=" << s.global_index.x << ','
                << s.global_index.y << ',' << s.global_index.z << " stage=" << s.stage
                << " active=" << unsigned(s.active) << " rho=" << s.rho
                << " h=" << s.h << " T=" << s.temperature << " p=" << s.pressure
                << " rate=" << s.rate << " rho_n=" << s.rho_accepted
                << " rho_nm1=" << s.rho_previous << " h_n=" << s.h_accepted
                << " h_nm1=" << s.h_previous << '\n';
    }
    if (!status) {
      if (rank == 0) {
        std::cerr << "placeholder advance=" << unsigned(status.code)
                  << '/' << status.detail << " stage=" << step.failed_stage << '\n';
        (void)write_numerical_failure(std::cerr, step.numerical_failure);
      }
      return false;
    }
    if (!agree(static_cast<bool>(driver.committed_output_snapshot(snapshot)))) return false;
    bool trace_valid = step.cell_trace.dropped == 0U;
    for (std::size_t k = 0U; k < step.cell_trace.count; ++k) {
      const auto& s = step.cell_trace.samples[k];
      if (s.attempt != step.attempts) trace_valid = false;
      if (s.stage >= 12U && (s.h != s.h_accepted ||
          std::abs(s.rho - s.rho_accepted) > 4e-15 || s.rate != 0.0)) trace_valid = false;
    }
    if (!agree(trace_valid)) return false;
    const auto current_h = enthalpy(snapshot);
    flat = 0U;
    for (int z = 0; z < cells.z; ++z)
      for (int y = 0; y < cells.y; ++y)
        for (int x = 0; x < cells.x; ++x, ++flat)
          if (snapshot.cell_activity.data[flat] == 0U)
            maximum_drift = std::max(maximum_drift,
                std::abs(current_h.unchecked({x, y, z}, 0U) - initial_h[flat]));
  }
  double global_drift = 0.0;
  MPI_Allreduce(&maximum_drift, &global_drift, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  if (rank == 0)
    std::cerr << std::setprecision(17) << "stationary-placeholder-h-drift="
              << global_drift << '\n';
  return global_drift == 0.0;
}

bool native_placeholder_perturbation(int rank, double solid_temperature,
                                     std::vector<double>& fluid_state) {
  auto model = test::product_model({16, 16, 16});
  model.mesh.lower = {-2, -2, -2};
  model.mesh.upper = {2, 2, 2};
  model.mesh.minimum_spacing = {0.25, 0.25, 0.25};
  model.time.control = TimeControlKind::fixed;
  model.time.initial_dt = model.time.minimum_dt = model.time.maximum_dt = 1e-5;
  model.immersed_boundary = ImmersedBoundarySpec{"cylinder_ascii.stl", ImmersedFluidSide::outside};
  auto& spec = model.thermophysics;
  spec.minimum_temperature = 273.15;
  spec.maximum_temperature = 6000.0;
  spec.maximum_temperature_iterations = 200U;
  auto& air = spec.species[0];
  air.molecular_weight = 28.850334;
  air.nasa7_low = {3.5838100068, -7.2700635412e-4, 1.67056387003e-6,
                   -1.091801341e-10, -4.317787988e-13, -1050.5394088, 3.1124135035};
  air.nasa7_high = {3.1013370688, 1.24138813631e-3, -4.1882038804e-7,
                    6.641656204e-11, -3.9127843272e-15, -985.27467132, 5.3560174057};
  air.transport_law = TransportLaw::coast_native_air;
  air.viscosity_reference = air.conductivity = 0.0;
  ThermodynamicsPlan thermo;
  Status status = ThermodynamicsPlan::compile(spec, {}, thermo);
  if (!status && rank == 0) std::cerr << "native thermo compile=" << unsigned(status.code) << '/' << status.detail << '\n';
  CompiledCasePlan plan;
  if (status) status = ProductCompiler::compile(MPI_COMM_WORLD, model,
      std::filesystem::path{HUNDUN_V04_SOURCE_ROOT} / "tests" / "data", plan);
  ProductDriver driver;
  if (status) status = ProductDriver::create(MPI_COMM_WORLD, std::move(plan), driver);
  RestartExpected expected;
  if (status) status = driver.restart_expected(expected);
  if (status) status = driver.initialize({});
  CommittedOutputSnapshot output;
  if (status) status = driver.committed_output_snapshot(output);
  if (!agree(static_cast<bool>(status))) {
    if (rank == 0) std::cerr << "native setup=" << unsigned(status.code) << '/' << status.detail << '\n';
    return false;
  }
  RestartImage image;
  image.global_cells = expected.global_cells;
  image.patch = expected.target_patch;
  image.plan = expected.plan; image.schema = expected.schema; image.geometry = expected.geometry;
  image.time = image.dt = 1e-5; image.step = image.controller_state = 1U;
  image.pressure_reference = 101325.0;
  image.backward_euler_recovery = true; // A genuine current-state-only V1 fixture.
  const Int3 n = image.patch.cells;
  const auto activity = output.cell_activity;
  const std::size_t count = std::size_t(n.x) * n.y * n.z;
  for (std::size_t k = 0U; k < expected.fields.size; ++k) {
    const auto descriptor = expected.fields.data[k];
    RestartImageField field;
    field.role = descriptor.role; field.field = descriptor.field; field.components = descriptor.components;
    field.values.assign(count * field.components, 0.0);
    for (int z = 0; z < n.z; ++z)
      for (int y = 0; y < n.y; ++y)
        for (int x = 0; x < n.x; ++x) {
          const std::size_t flat = std::size_t(x) + std::size_t(n.x) * (y + std::size_t(n.y) * z);
          if (field.role == RestartFieldRole::velocity && activity.data[flat] != 0U)
            field.values[flat * 3U] = 0.1;
          if (field.role == RestartFieldRole::enthalpy) {
            const double t = activity.data[flat] == 0U ? solid_temperature :
                300.0 + 0.2 * std::cos(2.0 * 3.141592653589793 *
                    (image.patch.begin.x + x + 0.5) / 16.0);
            double h = 0.0, cp = 0.0, gas = 0.0;
            status = thermo.mixture_enthalpy(t, {}, h, cp, gas);
            if (!status) return false;
            field.values[flat] = h;
          }
        }
    image.fields.push_back(std::move(field));
  }
  image.final_mass_flux[0].assign(std::size_t(n.x + 1) * n.y * n.z, 0.0);
  image.final_mass_flux[1].assign(std::size_t(n.x) * (n.y + 1) * n.z, 0.0);
  image.final_mass_flux[2].assign(std::size_t(n.x) * n.y * (n.z + 1), 0.0);
  CompiledCasePlan restart_plan;
  status = ProductCompiler::compile(MPI_COMM_WORLD, model,
      std::filesystem::path{HUNDUN_V04_SOURCE_ROOT} / "tests" / "data", restart_plan);
  ProductDriver restarted_driver;
  if (status) status = ProductDriver::create(MPI_COMM_WORLD, std::move(restart_plan), restarted_driver);
  if (status) status = restarted_driver.initialize_restart(image);
  DriverStepReport step;
  if (status) status = restarted_driver.advance({1, 1, 1, 1, 1}, step);
  if (!status) {
    if (rank == 0) std::cerr << "native-placeholder variant=" << solid_temperature
        << " status=" << unsigned(status.code) << '/' << status.detail
        << " stage=" << step.failed_stage << '\n';
    return false;
  }
  status = restarted_driver.committed_output_snapshot(output);
  if (!agree(static_cast<bool>(status))) return false;
  fluid_state.clear();
  for (std::size_t k = 0U; k < output.fields.size; ++k) {
    const auto field = output.fields.data[k].values;
    for (int z = 0; z < n.z; ++z)
      for (int y = 0; y < n.y; ++y)
        for (int x = 0; x < n.x; ++x) {
          const std::size_t flat = std::size_t(x) + std::size_t(n.x) * (y + std::size_t(n.y) * z);
          if (output.cell_activity.data[flat] == 0U) continue;
          for (std::uint8_t c = 0U; c < field.components; ++c)
            fluid_state.push_back(field.unchecked({x, y, z}, c));
        }
  }
  return true;
}

// Manufacture a checksummed complete-history fixture through public snapshots
// and Writer. One synthetic accepted source rate exhausts the deliberately
// iteration-limited CLI inversion. This is fault injection, not a physical run or
// a conversion of an existing complete image into a history-less image.
bool seed_cli_failure(const char* case_root, const char* restart_root, int target) {
  ValidatedModel model;
  Status status = CaseCompiler::load_and_compile(MPI_COMM_WORLD, case_root, model);
  CompiledCasePlan plan;
  if (status) status = ProductCompiler::compile(MPI_COMM_WORLD, model, case_root, plan);
  ProductDriver driver;
  if (status) status = ProductDriver::create(MPI_COMM_WORLD, std::move(plan), driver);
  DriverInitialState initial;
  initial.temperature = 300.0;
  if (status) status = driver.initialize(initial);
  DriverStepReport uniform;
  if (status) status = driver.advance({1, 1, 1, 1, 1}, uniform);
  CommittedOutputSnapshot output;
  RestartSnapshot base;
  if (status) status = driver.committed_output_snapshot(output);
  if (status) status = driver.committed_restart_snapshot(base);
  if (!agree(bool(status))) {
    std::cerr << "seed capture status=" << unsigned(status.code) << '/' << status.detail << '\n';
    return false;
  }
  std::vector<std::vector<double>> storage(base.accepted_rate_fields.size);
  std::vector<RestartFieldView> fields(base.accepted_rate_fields.size);
  const Int3 n = base.patch.cells;
  const std::size_t count = std::size_t(n.x) * n.y * n.z;
  for (std::size_t f = 0U; f < base.accepted_rate_fields.size; ++f) {
    const auto source = base.accepted_rate_fields.data[f];
    storage[f].resize(count * source.values.components);
    auto view = source.values;
    view.base = storage[f].data();
    view.storage_identity = 1000U + f;
    view.ghosts = {};
    view.stride_y = n.x; view.stride_z = std::size_t(n.x) * n.y;
    view.component_stride = count;
    for (int z = 0; z < n.z; ++z)
      for (int y = 0; y < n.y; ++y)
        for (int x = 0; x < n.x; ++x) {
          const std::size_t flat = std::size_t(x) + std::size_t(n.x) * (y + std::size_t(n.y) * z);
          for (std::uint8_t c = 0U; c < view.components; ++c)
            storage[f][flat + c * count] = source.values.unchecked({x,y,z}, c);
          if (source.role == RestartFieldRole::enthalpy_nonadvective_rate && output.cell_activity.data[flat] != 0U &&
              base.patch.begin.x + x == (target == 0 ? 3 : 12) &&
              base.patch.begin.y + y == 2 && base.patch.begin.z + z == (target == 0 ? 3 : 12))
            storage[f][flat] = 1.0e8;
        }
    fields[f] = {source.role, view};
  }
  RestartSnapshot seed = base;
  seed.accepted_rate_fields = {fields.data(), fields.size()};
  status = RestartWriter::write(MPI_COMM_WORLD, restart_root, seed);
  if (!status) std::cerr << "seed writer status=" << unsigned(status.code) << '/' << status.detail << '\n';
  int rank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (status && rank == 0)
    std::cout << "seed plan=" << seed.plan << " schema=" << seed.schema << '\n';
  return agree(bool(status));
}

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 2;
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (argc == 5 && std::string_view(argv[1]) == "--seed-cli-failure") {
    const bool okay = seed_cli_failure(argv[2], argv[3], std::atoi(argv[4]));
    MPI_Finalize();
    return okay ? 0 : 1;
  }
  bool passed = run(rank);
  if (agree(passed)) {
    std::vector<double> reference, perturbed;
    passed = native_placeholder_perturbation(rank, 280.0, reference);
    if (agree(passed)) passed = native_placeholder_perturbation(rank, 370.0, perturbed);
    if (agree(passed)) {
      double error = 0.0;
      if (reference.size() != perturbed.size()) passed = false;
      else for (std::size_t k = 0U; k < reference.size(); ++k)
        error = std::max(error, std::abs(reference[k] - perturbed[k]) / std::max(1.0, std::abs(reference[k])));
      double global_error = 0.0;
      MPI_Allreduce(&error, &global_error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      if (rank == 0) std::cerr << "native-fluid-placeholder-relative-response=" << global_error << '\n';
      passed = passed && global_error < 1e-10;
    }
  }
  const int local = passed ? 1 : 0;
  int global = 0;
  MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Finalize();
  return global ? 0 : 1;
}
