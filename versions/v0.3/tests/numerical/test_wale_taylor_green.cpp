// SPDX-License-Identifier: Apache-2.0

#include "tests/support/flow_immersed_test_access.hpp"
#include "tests/support/stage3_scientific_row.hpp"
#include "tests/support/test_main.hpp"

#include "hundun/bc_basic_boundary.hpp"
#include "hundun/exec_execution.hpp"
#include "hundun/flow_immersed.hpp"
#include "hundun/fvm_cell_centered.hpp"
#include "hundun/les_wale.hpp"
#include "hundun/lin_conjugate_gradient.hpp"
#include "hundun/lin_preconditioners.hpp"
#include "hundun/mesh_geometry.hpp"
#include "hundun/mesh_topology.hpp"
#include "hundun/rt_exchange_plan.hpp"
#include "hundun/rt_halo_exchange.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "hundun/rt_structured_decomposition.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <string>
#include <variant>
#include <vector>

namespace {

using namespace hundun;

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kRho = 1.0;
constexpr double kNu = 0.01;
constexpr double kFinalTime = 1.0e-4;
// A uniform spanwise Galilean velocity leaves the 2-D Taylor--Green solution
// unchanged while giving the third momentum equation a non-degenerate scale.
constexpr double kSpanwiseVelocity = 0.125;

runtime::Int3 process_grid(int ranks) {
  if (ranks == 1)
    return {1, 1, 1};
  if (ranks == 2)
    return {2, 1, 1};
  if (ranks == 4)
    return {2, 2, 1};
  throw runtime::Error("unsupported WALE Taylor-Green rank count");
}

config::FlowCaseConfig periodic_case() {
  config::FlowCaseConfig result{};
  result.schema_version = 2;
  result.simulation_type = config::SimulationType::variable_density_flow;
  result.density_model = config::DensityModel::constant;
  result.physics.rho_ref_kg_per_m3 = kRho;
  result.physics.dynamic_viscosity_pa_s = kRho * kNu;
  result.physics.inlet_consistency_rtol = 1.0e-12;
  constexpr std::array<config::PatchName, 6> names{
      config::PatchName::x_min, config::PatchName::x_max,
      config::PatchName::y_min, config::PatchName::y_max,
      config::PatchName::z_min, config::PatchName::z_max};
  for (std::size_t index = 0U; index < names.size(); ++index) {
    result.boundaries[index].patch = names[index];
    result.boundaries[index].type = config::BoundaryType::periodic;
  }
  return result;
}

runtime::FieldDescriptor cell_field(const char *name,
                                    std::uint32_t components) {
  return {name,
          "1",
          "stage3_s1_wale_tgv",
          runtime::FunctionSpace::cell_average,
          runtime::ScalarType::float64,
          components,
          2,
          false,
          runtime::RestartPolicy::persistent,
          runtime::OutputPolicy::never};
}

runtime::FieldDescriptor face_field(const char *name,
                                    std::uint32_t components) {
  return {name,
          "1",
          "stage3_s1_wale_tgv",
          runtime::FunctionSpace::face_value,
          runtime::ScalarType::float64,
          components,
          0,
          false,
          runtime::RestartPolicy::persistent,
          runtime::OutputPolicy::never};
}

double average_sin(int cell, int cells, int wave = 1) noexcept {
  const double lower = 2.0 * kPi * static_cast<double>(cell) / cells;
  const double upper = 2.0 * kPi * static_cast<double>(cell + 1) / cells;
  return (std::cos(wave * lower) - std::cos(wave * upper)) /
         (wave * (upper - lower));
}

double average_cos(int cell, int cells, int wave = 1) noexcept {
  const double lower = 2.0 * kPi * static_cast<double>(cell) / cells;
  const double upper = 2.0 * kPi * static_cast<double>(cell + 1) / cells;
  return (std::sin(wave * upper) - std::sin(wave * lower)) /
         (wave * (upper - lower));
}

runtime::Real3 exact_velocity_average(runtime::Int3 cell, int cells,
                                      double time_s) noexcept {
  const double decay = std::exp(-2.0 * kNu * time_s);
  return {average_sin(cell.x, cells) * average_cos(cell.y, cells) * decay,
          -average_cos(cell.x, cells) * average_sin(cell.y, cells) * decay,
          kSpanwiseVelocity};
}

double exact_pressure_average(runtime::Int3 cell, int cells,
                              double time_s) noexcept {
  return 0.25 * kRho *
         (average_cos(cell.x, cells, 2) + average_cos(cell.y, cells, 2)) *
         std::exp(-4.0 * kNu * time_s);
}

runtime::Real3 exact_velocity_point(runtime::Real3 point,
                                    double time_s) noexcept {
  const double decay = std::exp(-2.0 * kNu * time_s);
  return {std::sin(point.x) * std::cos(point.y) * decay,
          -std::cos(point.x) * std::sin(point.y) * decay, kSpanwiseVelocity};
}

flow::FlowLayerValues exact_layer(const mesh::MeshTopology &topology,
                                  const mesh::MeshGeometry &geometry, int cells,
                                  double time_s) {
  flow::FlowLayerValues result;
  result.density.assign(topology.owned_cell_count(), kRho);
  result.velocity.resize(topology.owned_cell_count() * 3U);
  result.mechanical_pressure.resize(topology.owned_cell_count());
  result.face_velocity.resize(topology.local_face_count() * 3U);
  result.face_mass_flux.resize(topology.local_face_count());
  for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count();
       ++cell) {
    const auto exact =
        exact_velocity_average(topology.global_cell(cell), cells, time_s);
    result.velocity[cell * 3U] = exact.x;
    result.velocity[cell * 3U + 1U] = exact.y;
    result.velocity[cell * 3U + 2U] = exact.z;
    result.mechanical_pressure[cell] =
        exact_pressure_average(topology.global_cell(cell), cells, time_s);
  }
  for (mesh::LocalFaceId face = 0U; face < topology.local_face_count();
       ++face) {
    const auto exact =
        exact_velocity_point(geometry.face_center_m(face), time_s);
    result.face_velocity[face * 3U] = exact.x;
    result.face_velocity[face * 3U + 1U] = exact.y;
    result.face_velocity[face * 3U + 2U] = exact.z;
    const auto area = geometry.face_area_vector_m2(face, mesh::FaceSide::owner);
    result.face_mass_flux[face] =
        kRho * (exact.x * area.x + exact.y * area.y + exact.z * area.z);
  }
  return result;
}

struct LevelResult final {
  test::stage3::ScientificRow row;
  test::stage3::CellAverageSnapshot fields;
  double velocity_l2{};
  double pressure_l2{};
  double nu_t_l2{};
  std::uint64_t density_revision{};
  std::uint64_t wale_evaluation_count{};
};

std::size_t global_offset(runtime::Int3 cell, int cells) noexcept {
  return (static_cast<std::size_t>(cell.z) * static_cast<std::size_t>(cells) +
          static_cast<std::size_t>(cell.y)) *
             static_cast<std::size_t>(cells) +
         static_cast<std::size_t>(cell.x);
}

LevelResult run_level(const runtime::MpiContext &mpi, int cells,
                      std::string row_id) {
  HUNDUN_CHECK(cells == 12 || cells == 24 || cells == 48);
  const int ratio = cells / 12;
  const std::uint64_t steps = static_cast<std::uint64_t>(ratio);
  const double dt = kFinalTime / static_cast<double>(steps);
  const runtime::Int3 extent{cells, cells, cells};
  const auto grid = process_grid(mpi.size());
  auto decomposition = runtime::StructuredDecomposition::create(
      mpi, extent, {true, true, true}, runtime::DecompositionOptions{grid});
  const mesh::MeshTopology topology(decomposition);
  const mesh::MeshGeometry geometry(
      topology, mesh::UniformBoxMapping({0.0, 0.0, 0.0},
                                        {2.0 * kPi, 2.0 * kPi, 2.0 * kPi}));
  const auto boundaries =
      boundary::BoundaryRegistry::create(periodic_case(), topology);

  runtime::FieldRegistry registry;
  flow::FlowFieldIds fields;
  fields.density = registry.declare_field(cell_field("rho", 1U));
  fields.velocity = registry.declare_field(cell_field("velocity", 3U));
  fields.mechanical_pressure = registry.declare_field(cell_field("pi", 1U));
  fields.face_velocity =
      registry.declare_field(face_field("face_velocity", 3U));
  fields.face_mass_flux = finite_volume::declare_face_mass_flux(registry);
  registry.freeze();
  auto state = flow::FlowState::create(
      registry, {decomposition.local_extent(), topology.local_face_count()},
      fields, {0U, 0.0, dt, 0.0, flow::MomentumTimeOrder::backward_euler});
  state.seed_accepted_layers(exact_layer(topology, geometry, cells, -dt),
                             exact_layer(topology, geometry, cells, 0.0));

  execution::CpuReferenceContext execution;
  std::vector<mesh::GlobalCellId> active;
  active.reserve(topology.local_cell_count());
  for (mesh::LocalCellId cell = 0U; cell < topology.local_cell_count(); ++cell)
    active.push_back(topology.global_cell_id(cell));
  auto wale =
      les::WaleModel::create({0.5, 0.9, 0.7}, topology, geometry,
                             topology.owned_cell_count(), active, execution);
  auto halo = runtime::HaloExchange::create(
      decomposition, runtime::ExchangePlan::create(
                         decomposition, decomposition.local_extent(), 2));
  linear::ConjugateGradientSolver momentum_solver(execution, mpi);
  linear::ConjugateGradientSolver pressure_solver(execution, mpi);
  linear::JacobiPreconditioner mx(execution), my(execution), mz(execution),
      pressure_preconditioner(execution);
  auto facade = flow::FixedStepImmersedFlow::create(
      decomposition, topology, geometry, boundaries, nullptr, nullptr, nullptr,
      nullptr, &wale, mpi, execution, halo, momentum_solver, {&mx, &my, &mz},
      pressure_solver, pressure_preconditioner);
  const flow::ImmersedFlowPhysics physics{config::DensityModel::constant,
                                          kRho,
                                          kRho * kNu,
                                          std::nullopt,
                                          std::nullopt,
                                          std::nullopt};
  linear::SolveControl control;
  control.max_iterations = 5000U;
  control.residual_recompute_interval = 20U;
  double maximum_continuity = 0.0;
  double maximum_conservation = 0.0;
  les::WaleSummary final_wale;
  std::uint64_t final_wale_evaluation_count{};
  for (std::uint64_t step = 0U; step < steps; ++step) {
    const auto order = step == 0U ? flow::MomentumTimeOrder::backward_euler
                                  : flow::MomentumTimeOrder::bdf2;
    const auto stencil =
        flow::make_momentum_time_stencil(order, dt, step == 0U ? 0.0 : dt);
    const auto report =
        facade.attempt(state, physics, stencil, control, control);
    const auto *base = std::get_if<flow::StepAttemptReport>(&report.base);
    HUNDUN_CHECK(base != nullptr);
    if (base->disposition != flow::StepAttemptDisposition::committed &&
        mpi.rank() == 0) {
      std::cerr << "WALE_TGV_FAILURE cells=" << cells << " step=" << step
                << " reason=" << static_cast<int>(base->reason)
                << " rank=" << base->lowest_failing_rank
                << " momentum=" << base->final_momentum_normalized_l2[0] << ','
                << base->final_momentum_normalized_l2[1] << ','
                << base->final_momentum_normalized_l2[2]
                << " continuity=" << base->final_continuity_normalized_l2
                << " pressure=" << base->final_pressure_residual_l2 << '\n';
    }
    HUNDUN_CHECK(base->disposition == flow::StepAttemptDisposition::committed);
    HUNDUN_CHECK(base->pressure_corrector_count == 2U);
    HUNDUN_CHECK(report.wale.has_value());
    HUNDUN_CHECK(!report.force.has_value());
    final_wale = *report.wale;
    final_wale_evaluation_count =
        flow::test::ImmersedFlowTestAccess::wale_evaluation_count(facade);
    HUNDUN_CHECK(final_wale_evaluation_count == 1U);
    HUNDUN_CHECK(flow::test::ImmersedFlowTestAccess::wale_coefficient_identity(
                     facade) == final_wale.identity);
    maximum_continuity =
        std::max(maximum_continuity, base->final_continuity_normalized_l2);
    maximum_conservation = std::max(
        maximum_conservation, base->final_mass_relative_conservation_defect);
    for (double value : base->final_momentum_relative_conservation_defect)
      maximum_conservation = std::max(maximum_conservation, value);
  }
  HUNDUN_CHECK(state.metadata().time_s == kFinalTime);
  HUNDUN_CHECK(final_wale.identity.value != 0U);
  HUNDUN_CHECK(final_wale.l2_nu_t_m2_per_s > 0.0);

  const auto committed = state.snapshot(flow::FlowLayer::committed);
  double sums[4]{};
  for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count();
       ++cell) {
    const double volume = geometry.cell_volume_m3(cell);
    const auto exact =
        exact_velocity_average(topology.global_cell(cell), cells, kFinalTime);
    const double pressure_exact =
        exact_pressure_average(topology.global_cell(cell), cells, kFinalTime);
    const std::array<double, 3> exact_components{exact.x, exact.y, exact.z};
    for (std::size_t component = 0U; component < 3U; ++component) {
      const double difference = committed.velocity[cell * 3U + component] -
                                exact_components[component];
      sums[0] += volume * difference * difference;
    }
    sums[1] += volume * committed.mechanical_pressure[cell];
    sums[2] += volume * pressure_exact;
    sums[3] += volume;
  }
  mpi.allreduce_fp64_in_place(sums, 4U, runtime::Fp64ReductionOperation::sum);
  const double numerical_pressure_mean = sums[1] / sums[3];
  const double exact_pressure_mean = sums[2] / sums[3];
  double pressure_error_square = 0.0;
  for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count();
       ++cell) {
    const double volume = geometry.cell_volume_m3(cell);
    const double exact =
        exact_pressure_average(topology.global_cell(cell), cells, kFinalTime) -
        exact_pressure_mean;
    const double numerical =
        committed.mechanical_pressure[cell] - numerical_pressure_mean;
    const double difference = numerical - exact;
    pressure_error_square += volume * difference * difference;
  }
  mpi.allreduce_fp64_in_place(&pressure_error_square, 1U,
                              runtime::Fp64ReductionOperation::sum);
  const double velocity_l2 = std::sqrt(sums[0] / sums[3]);
  const double pressure_l2 = std::sqrt(pressure_error_square / sums[3]);
  HUNDUN_CHECK(velocity_l2 > 0.0 && std::isfinite(velocity_l2));
  HUNDUN_CHECK(pressure_l2 > 0.0 && std::isfinite(pressure_l2));

  double nu_t_square =
      final_wale.l2_nu_t_m2_per_s * final_wale.l2_nu_t_m2_per_s;
  mpi.allreduce_fp64_in_place(&nu_t_square, 1U,
                              runtime::Fp64ReductionOperation::sum);
  const double nu_t_l2 = std::sqrt(nu_t_square * geometry.cell_volume_m3(0U));
  HUNDUN_CHECK(nu_t_l2 > 0.0 && std::isfinite(nu_t_l2));

  const auto global_count = static_cast<std::size_t>(cells) *
                            static_cast<std::size_t>(cells) *
                            static_cast<std::size_t>(cells);
  test::stage3::CellAverageSnapshot snapshot;
  snapshot.cells = extent;
  snapshot.velocity.assign(global_count * 3U, 0.0);
  snapshot.pressure.assign(global_count, 0.0);
  std::vector<double> coverage(global_count, 0.0);
  for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count();
       ++cell) {
    const auto offset = global_offset(topology.global_cell(cell), cells);
    HUNDUN_CHECK(offset < global_count);
    HUNDUN_CHECK(coverage[offset] == 0.0);
    coverage[offset] = 1.0;
    for (std::size_t component = 0U; component < 3U; ++component)
      snapshot.velocity[offset * 3U + component] =
          committed.velocity[cell * 3U + component];
    snapshot.pressure[offset] =
        committed.mechanical_pressure[cell] - numerical_pressure_mean;
  }
  mpi.allreduce_fp64_in_place(snapshot.velocity.data(),
                              snapshot.velocity.size(),
                              runtime::Fp64ReductionOperation::sum);
  mpi.allreduce_fp64_in_place(snapshot.pressure.data(),
                              snapshot.pressure.size(),
                              runtime::Fp64ReductionOperation::sum);
  mpi.allreduce_fp64_in_place(coverage.data(), coverage.size(),
                              runtime::Fp64ReductionOperation::sum);
  HUNDUN_CHECK(std::all_of(coverage.begin(), coverage.end(),
                           [](double value) { return value == 1.0; }));

  test::stage3::ScientificRow row{std::move(row_id),
                                  extent,
                                  mpi.size(),
                                  grid,
                                  steps,
                                  dt,
                                  kFinalTime,
                                  {true, velocity_l2},
                                  {true, pressure_l2},
                                  {true, nu_t_l2},
                                  {true, maximum_continuity},
                                  {true, maximum_conservation},
                                  {false, 0.0},
                                  {false, 0.0},
                                  {true, final_wale.identity.value},
                                  "pass"};
  HUNDUN_CHECK(test::stage3::validate_scientific_row(row));
  return {std::move(row),
          std::move(snapshot),
          velocity_l2,
          pressure_l2,
          nu_t_l2,
          state.metadata().step,
          final_wale_evaluation_count};
}

void print_row(const runtime::MpiContext &mpi,
               const test::stage3::ScientificRow &row) {
  if (mpi.rank() == 0)
    std::cout << test::stage3::serialize_scientific_row(row) << '\n';
}

struct FieldDifference final {
  double velocity_l2{};
  double pressure_l2{};
};

FieldDifference
field_difference(const test::stage3::CellAverageSnapshot &left,
                 const test::stage3::CellAverageSnapshot &right) {
  HUNDUN_CHECK(left.cells.x == right.cells.x && left.cells.y == right.cells.y &&
               left.cells.z == right.cells.z);
  HUNDUN_CHECK(left.velocity.size() == right.velocity.size());
  HUNDUN_CHECK(left.pressure.size() == right.pressure.size());
  HUNDUN_CHECK(!left.velocity.empty() && !left.pressure.empty());
  double velocity_square = 0.0;
  for (std::size_t index = 0U; index < left.velocity.size(); ++index) {
    const double difference = left.velocity[index] - right.velocity[index];
    velocity_square += difference * difference;
  }
  double pressure_square = 0.0;
  for (std::size_t index = 0U; index < left.pressure.size(); ++index) {
    const double difference = left.pressure[index] - right.pressure[index];
    pressure_square += difference * difference;
  }
  return {
      std::sqrt(velocity_square / static_cast<double>(left.pressure.size())),
      std::sqrt(pressure_square / static_cast<double>(left.pressure.size()))};
}

void require_cell_average_self_convergence(
    const std::array<LevelResult, 3> &levels) {
  const auto medium_on_coarse = test::stage3::restrict_cell_averages(
      levels[1].fields, levels[0].fields.cells,
      test::stage3::CellRestriction::cell_average);
  const auto fine_on_medium = test::stage3::restrict_cell_averages(
      levels[2].fields, levels[1].fields.cells,
      test::stage3::CellRestriction::cell_average);
  const auto coarse_medium =
      field_difference(levels[0].fields, medium_on_coarse);
  const auto medium_fine = field_difference(levels[1].fields, fine_on_medium);
  for (const std::array<double, 2> errors :
       {std::array<double, 2>{coarse_medium.velocity_l2,
                              medium_fine.velocity_l2},
        std::array<double, 2>{coarse_medium.pressure_l2,
                              medium_fine.pressure_l2}}) {
    HUNDUN_CHECK(std::isfinite(errors[0]) && std::isfinite(errors[1]));
    HUNDUN_CHECK(errors[0] > errors[1] && errors[1] > 0.0);
    const double order = std::log(errors[0] / errors[1]) / std::log(2.0);
    HUNDUN_CHECK(std::isfinite(order) && order >= 1.8);
  }
}

void require_decomposition_match(const runtime::MpiContext &mpi,
                                 const LevelResult &distributed,
                                 const LevelResult &single_rank) {
  HUNDUN_CHECK(distributed.row.cells.x == single_rank.row.cells.x &&
               distributed.row.cells.y == single_rank.row.cells.y &&
               distributed.row.cells.z == single_rank.row.cells.z);
  HUNDUN_CHECK(distributed.row.steps == single_rank.row.steps);
  HUNDUN_CHECK(distributed.row.dt_s == single_rank.row.dt_s);
  HUNDUN_CHECK(distributed.row.final_time_s == single_rank.row.final_time_s);
  constexpr std::size_t field_count = 5U;
  std::array<double, field_count * 2U> values{};
  const auto compare = [&](std::size_t field, double actual, double reference) {
    values[field] = std::max(values[field], std::abs(actual - reference));
    values[field + field_count] =
        std::max(values[field + field_count], std::abs(reference));
  };
  HUNDUN_CHECK(distributed.fields.velocity.size() ==
               single_rank.fields.velocity.size());
  HUNDUN_CHECK(distributed.fields.pressure.size() ==
               single_rank.fields.pressure.size());
  for (std::size_t index = 0U; index < distributed.fields.velocity.size();
       ++index)
    compare(0U, distributed.fields.velocity[index],
            single_rank.fields.velocity[index]);
  for (std::size_t index = 0U; index < distributed.fields.pressure.size();
       ++index)
    compare(1U, distributed.fields.pressure[index],
            single_rank.fields.pressure[index]);
  compare(2U, distributed.velocity_l2, single_rank.velocity_l2);
  compare(3U, distributed.pressure_l2, single_rank.pressure_l2);
  compare(4U, distributed.nu_t_l2, single_rank.nu_t_l2);
  mpi.allreduce_fp64_in_place(values.data(), values.size(),
                              runtime::Fp64ReductionOperation::maximum);
  for (std::size_t field = 0U; field < field_count; ++field)
    HUNDUN_CHECK(values[field] <=
                 5.0e-12 * std::max(1.0, values[field + field_count]));
}

void run_smoke(const runtime::MpiContext &mpi) {
  HUNDUN_CHECK(mpi.size() == 1);
  const auto result = run_level(mpi, 12, "wale-tgv-n12-smoke-r1");
  print_row(mpi, result.row);
}

void run_formal_n24(const runtime::MpiContext &mpi) {
  const auto result =
      run_level(mpi, 24, "wale-tgv-n24-r" + std::to_string(mpi.size()));
  if (mpi.size() > 1) {
    auto self = runtime::MpiContext::duplicate(MPI_COMM_SELF);
    const auto single_rank =
        run_level(self, 24, "wale-tgv-n24-decomposition-reference");
    require_decomposition_match(mpi, result, single_rank);
  }
  print_row(mpi, result.row);
}

void run_formal_convergence(const runtime::MpiContext &mpi) {
  HUNDUN_CHECK(mpi.size() == 1);
  std::array<LevelResult, 3> levels{run_level(mpi, 12, "wale-tgv-n12-r1"),
                                    run_level(mpi, 24, "wale-tgv-n24-r1"),
                                    run_level(mpi, 48, "wale-tgv-n48-r1")};
  require_cell_average_self_convergence(levels);
  test::stage3::TgvConvergenceInput input;
  for (std::size_t level = 0U; level < levels.size(); ++level) {
    input.velocity_l2[level] = levels[level].velocity_l2;
    input.pressure_l2[level] = levels[level].pressure_l2;
    input.nu_t_l2[level] = levels[level].nu_t_l2;
    input.final_time_s[level] = levels[level].row.final_time_s;
    input.density_revision[level] = levels[level].density_revision;
    input.wale_density_revision[level] = levels[level].density_revision;
    input.wale_evaluation_count[level] = levels[level].wale_evaluation_count;
  }
  const auto convergence = test::stage3::assess_tgv_convergence(input);
  HUNDUN_CHECK(convergence.accepted);
  auto aggregate = levels.back().row;
  aggregate.row_id = "wale-tgv-convergence-r1";
  print_row(mpi, aggregate);
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  enum class Mode { invalid, smoke, formal_convergence, formal_n24 };
  Mode mode{Mode::invalid};
  if (argc == 3 && std::string(argv[1]) == "smoke" &&
      std::string(argv[2]) == "12")
    mode = Mode::smoke;
  else if (argc == 3 && std::string(argv[1]) == "formal" &&
           std::string(argv[2]) == "convergence")
    mode = Mode::formal_convergence;
  else if (argc == 3 && std::string(argv[1]) == "formal" &&
           std::string(argv[2]) == "24")
    mode = Mode::formal_n24;
  if (mode == Mode::invalid)
    return 2;
  auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
  return hundun::test::run([&] {
    if (mode == Mode::smoke)
      run_smoke(mpi);
    else if (mode == Mode::formal_convergence)
      run_formal_convergence(mpi);
    else
      run_formal_n24(mpi);
  });
}
