// SPDX-License-Identifier: Apache-2.0

#include "hundun/flow_passive_scalar.hpp"

#include "hundun/mesh_uniform_structured.hpp"
#include "hundun/rt_exchange_plan.hpp"
#include "hundun/rt_field_registry.hpp"
#include "hundun/rt_field_storage.hpp"
#include "hundun/rt_halo_exchange.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "hundun/rt_structured_decomposition.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

namespace {

using hundun::mesh::UniformStructuredMesh;
using hundun::runtime::ExchangePlan;
using hundun::runtime::FieldDescriptor;
using hundun::runtime::FieldId;
using hundun::runtime::FieldRegistry;
using hundun::runtime::FieldStorage;
using hundun::runtime::FunctionSpace;
using hundun::runtime::HaloExchange;
using hundun::runtime::Int3;
using hundun::runtime::MpiContext;
using hundun::runtime::MpiEnvironment;
using hundun::runtime::OutputPolicy;
using hundun::runtime::Real3;
using hundun::runtime::RestartPolicy;
using hundun::runtime::ScalarType;
using hundun::runtime::StructuredDecomposition;
using hundun::solver::global_l1_error;
using hundun::solver::global_mass;
using hundun::solver::mc_limiter;
using hundun::solver::PassiveScalarSolver;

constexpr std::array<bool, 3> kPeriodic{true, true, true};
constexpr double kPi = 3.141592653589793238462643383279502884;

static_assert(noexcept(mc_limiter(1.0, 1.0)));
static_assert(std::is_constructible_v<PassiveScalarSolver, const MpiContext &,
                                      const StructuredDecomposition &,
                                      const UniformStructuredMesh &,
                                      HaloExchange &, Real3, double>);

using AdvanceSignature = void (PassiveScalarSolver::*)(FieldStorage &, FieldId,
                                                       FieldId, double);
using MassSignature = double (*)(const MpiContext &,
                                 const UniformStructuredMesh &,
                                 const FieldStorage &, FieldId);
using L1Signature = double (*)(const MpiContext &,
                               const UniformStructuredMesh &,
                               const FieldStorage &, FieldId, FieldId);

[[maybe_unused]] constexpr AdvanceSignature kAdvanceSignature =
    &PassiveScalarSolver::advance_ssprk2;
[[maybe_unused]] constexpr MassSignature kMassSignature = &global_mass;
[[maybe_unused]] constexpr L1Signature kL1Signature = &global_l1_error;

FieldId declare_scalar(FieldRegistry &registry, const std::string &name,
                       int ghost_width = 2) {
  return registry.declare_field(FieldDescriptor{
      name, "1", "passive_scalar_test", FunctionSpace::cell_average,
      ScalarType::float64, 1U, ghost_width, true, RestartPolicy::persistent,
      OutputPolicy::selected});
}

std::size_t flat_index(Int3 extent, int i, int j, int k) {
  return (static_cast<std::size_t>(k) * static_cast<std::size_t>(extent.y) +
          static_cast<std::size_t>(j)) *
             static_cast<std::size_t>(extent.x) +
         static_cast<std::size_t>(i);
}

void check_limiter() {
  HUNDUN_CHECK_NEAR(mc_limiter(1.0, 1.0), 1.0, 0.0);
  HUNDUN_CHECK_NEAR(mc_limiter(-1.0, -1.0), -1.0, 0.0);
  HUNDUN_CHECK_NEAR(mc_limiter(1.0, -1.0), 0.0, 0.0);
  HUNDUN_CHECK_NEAR(mc_limiter(1.0, 3.0), 2.0, 0.0);
  HUNDUN_CHECK_NEAR(mc_limiter(3.0, 1.0), 2.0, 0.0);

  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double infinity = std::numeric_limits<double>::infinity();
  HUNDUN_CHECK(mc_limiter(nan, 1.0) == 0.0);
  HUNDUN_CHECK(mc_limiter(1.0, nan) == 0.0);
  HUNDUN_CHECK(mc_limiter(infinity, 1.0) == 0.0);
  HUNDUN_CHECK(mc_limiter(-infinity, 1.0) == 0.0);
  HUNDUN_CHECK(mc_limiter(1.0, infinity) == 0.0);
  HUNDUN_CHECK(mc_limiter(1.0, -infinity) == 0.0);
}

struct ConvergenceResult {
  double error{};
  double relative_mass_error{};
};

ConvergenceResult run_sine_resolution(const MpiContext &context, int nx) {
  const Int3 cells{nx, 4, 4};
  auto decomposition =
      StructuredDecomposition::create(context, cells, kPeriodic);
  UniformStructuredMesh mesh(cells, Real3{0.0, 0.0, 0.0},
                             Real3{1.0, 0.125, 0.125}, decomposition);

  FieldRegistry registry;
  const FieldId scalar = declare_scalar(registry, "scalar");
  const FieldId stage = declare_scalar(registry, "stage");
  const FieldId reference = declare_scalar(registry, "reference");
  registry.freeze();
  FieldStorage storage(registry, decomposition.local_extent());

  auto scalar_view = storage.view<double>(scalar);
  auto reference_view = storage.view<double>(reference);
  const Int3 local = decomposition.local_extent();
  for (int k = 0; k < local.z; ++k) {
    for (int j = 0; j < local.y; ++j) {
      for (int i = 0; i < local.x; ++i) {
        const double x = mesh.cell_center(Int3{i, j, k}).x;
        const double value = 1.0 + 0.2 * std::sin(2.0 * kPi * x);
        scalar_view(i, j, k, 0) = value;
        reference_view(i, j, k, 0) = value;
      }
    }
  }

  auto halo = HaloExchange::create(
      decomposition,
      ExchangePlan::create(decomposition, decomposition.local_extent(), 2));
  PassiveScalarSolver solver(context, decomposition, mesh, halo,
                             Real3{1.0, 0.0, 0.0}, 0.0);
  const double initial_mass = global_mass(context, mesh, storage, scalar);
  const int steps = static_cast<int>(std::ceil(static_cast<double>(nx) / 0.35));
  const double dt = 1.0 / static_cast<double>(steps);
  for (int step = 0; step < steps; ++step) {
    solver.advance_ssprk2(storage, scalar, stage, dt);
  }

  const double final_mass = global_mass(context, mesh, storage, scalar);
  return ConvergenceResult{
      global_l1_error(context, mesh, storage, scalar, reference),
      std::abs(final_mass - initial_mass) / std::abs(initial_mass)};
}

double independent_limiter(double left, double right) {
  if (!std::isfinite(left) || !std::isfinite(right) || left == 0.0 ||
      right == 0.0 || std::signbit(left) != std::signbit(right)) {
    return 0.0;
  }
  const double sign = left > 0.0 ? 1.0 : -1.0;
  return sign * std::min({2.0 * std::abs(left), 0.5 * std::abs(left + right),
                          2.0 * std::abs(right)});
}

int periodic_index(int coordinate, int extent) {
  const int reduced = coordinate % extent;
  return reduced < 0 ? reduced + extent : reduced;
}

double independent_sample(const std::vector<double> &values, Int3 extent, int i,
                          int j, int k) {
  return values[flat_index(extent, periodic_index(i, extent.x),
                           periodic_index(j, extent.y),
                           periodic_index(k, extent.z))];
}

double independent_face_flux(const std::vector<double> &values, Int3 extent,
                             Int3 lower, Int3 direction, double velocity) {
  if (velocity == 0.0) {
    return 0.0;
  }
  const auto sample = [&](int offset) {
    return independent_sample(values, extent, lower.x + offset * direction.x,
                              lower.y + offset * direction.y,
                              lower.z + offset * direction.z);
  };
  const double minus_one = sample(-1);
  const double lower_value = sample(0);
  const double upper_value = sample(1);
  const double plus_two = sample(2);
  const double lower_slope =
      independent_limiter(lower_value - minus_one, upper_value - lower_value);
  const double upper_slope =
      independent_limiter(upper_value - lower_value, plus_two - upper_value);
  const double left_state = lower_value + 0.5 * lower_slope;
  const double right_state = upper_value - 0.5 * upper_slope;
  return velocity * (velocity > 0.0 ? left_state : right_state);
}

double independent_spatial(const std::vector<double> &values, Int3 extent,
                           Real3 spacing, Real3 velocity, int i, int j, int k) {
  const Int3 cell{i, j, k};
  const double x_plus =
      independent_face_flux(values, extent, cell, Int3{1, 0, 0}, velocity.x);
  const double x_minus = independent_face_flux(
      values, extent, Int3{i - 1, j, k}, Int3{1, 0, 0}, velocity.x);
  const double y_plus =
      independent_face_flux(values, extent, cell, Int3{0, 1, 0}, velocity.y);
  const double y_minus = independent_face_flux(
      values, extent, Int3{i, j - 1, k}, Int3{0, 1, 0}, velocity.y);
  const double z_plus =
      independent_face_flux(values, extent, cell, Int3{0, 0, 1}, velocity.z);
  const double z_minus = independent_face_flux(
      values, extent, Int3{i, j, k - 1}, Int3{0, 0, 1}, velocity.z);
  return -(x_plus - x_minus) / spacing.x - (y_plus - y_minus) / spacing.y -
         (z_plus - z_minus) / spacing.z;
}

double check_independent_oracle(const MpiContext &context) {
  const Int3 cells{7, 6, 5};
  const Real3 lengths{1.4, 0.9, 0.7};
  const Real3 velocity{0.37, -0.29, 0.23};
  constexpr double dt = 0.02;
  auto decomposition =
      StructuredDecomposition::create(context, cells, kPeriodic);
  UniformStructuredMesh mesh(cells, Real3{-0.2, 0.1, -0.3}, lengths,
                             decomposition);

  FieldRegistry registry;
  const FieldId scalar = declare_scalar(registry, "oracle_scalar");
  const FieldId stage = declare_scalar(registry, "oracle_stage");
  registry.freeze();
  FieldStorage storage(registry, cells);
  auto scalar_view = storage.view<double>(scalar);

  const std::size_t count = static_cast<std::size_t>(cells.x) *
                            static_cast<std::size_t>(cells.y) *
                            static_cast<std::size_t>(cells.z);
  std::vector<double> initial(count);
  for (int k = 0; k < cells.z; ++k) {
    for (int j = 0; j < cells.y; ++j) {
      for (int i = 0; i < cells.x; ++i) {
        const double x =
            2.0 * kPi * static_cast<double>(i) / static_cast<double>(cells.x);
        const double y =
            2.0 * kPi * static_cast<double>(j) / static_cast<double>(cells.y);
        const double z =
            2.0 * kPi * static_cast<double>(k) / static_cast<double>(cells.z);
        const double value = 1.4 + 0.13 * std::sin(x) + 0.19 * std::cos(y) +
                             0.11 * std::sin(z) +
                             0.07 * std::cos(x - 0.5 * y + z);
        initial[flat_index(cells, i, j, k)] = value;
        scalar_view(i, j, k, 0) = value;
      }
    }
  }

  const Real3 spacing{lengths.x / static_cast<double>(cells.x),
                      lengths.y / static_cast<double>(cells.y),
                      lengths.z / static_cast<double>(cells.z)};
  std::vector<double> oracle_stage(count);
  for (int k = 0; k < cells.z; ++k) {
    for (int j = 0; j < cells.y; ++j) {
      for (int i = 0; i < cells.x; ++i) {
        const std::size_t index = flat_index(cells, i, j, k);
        oracle_stage[index] =
            initial[index] + dt * independent_spatial(initial, cells, spacing,
                                                      velocity, i, j, k);
      }
    }
  }
  std::vector<double> expected(count);
  for (int k = 0; k < cells.z; ++k) {
    for (int j = 0; j < cells.y; ++j) {
      for (int i = 0; i < cells.x; ++i) {
        const std::size_t index = flat_index(cells, i, j, k);
        expected[index] =
            0.5 * initial[index] +
            0.5 * (oracle_stage[index] +
                   dt * independent_spatial(oracle_stage, cells, spacing,
                                            velocity, i, j, k));
      }
    }
  }

  auto halo = HaloExchange::create(
      decomposition,
      ExchangePlan::create(decomposition, decomposition.local_extent(), 2));
  PassiveScalarSolver solver(context, decomposition, mesh, halo, velocity, 0.0);
  solver.advance_ssprk2(storage, scalar, stage, dt);

  double maximum_difference = 0.0;
  for (int k = 0; k < cells.z; ++k) {
    for (int j = 0; j < cells.y; ++j) {
      for (int i = 0; i < cells.x; ++i) {
        maximum_difference = std::max(
            maximum_difference, std::abs(scalar_view(i, j, k, 0) -
                                         expected[flat_index(cells, i, j, k)]));
      }
    }
  }
  HUNDUN_CHECK(maximum_difference < 2.0e-14);
  return maximum_difference;
}

double check_constant_state(const MpiContext &context) {
  const Int3 cells{9, 7, 5};
  auto decomposition =
      StructuredDecomposition::create(context, cells, kPeriodic);
  UniformStructuredMesh mesh(cells, Real3{0.0, 0.0, 0.0}, Real3{1.0, 0.8, 0.6},
                             decomposition);
  FieldRegistry registry;
  const FieldId scalar = declare_scalar(registry, "constant_scalar");
  const FieldId stage = declare_scalar(registry, "constant_stage");
  registry.freeze();
  FieldStorage storage(registry, cells);
  auto scalar_view = storage.view<double>(scalar);
  constexpr double constant = 3.25;
  for (int k = 0; k < cells.z; ++k) {
    for (int j = 0; j < cells.y; ++j) {
      for (int i = 0; i < cells.x; ++i) {
        scalar_view(i, j, k, 0) = constant;
      }
    }
  }
  auto halo = HaloExchange::create(
      decomposition,
      ExchangePlan::create(decomposition, decomposition.local_extent(), 2));
  PassiveScalarSolver solver(context, decomposition, mesh, halo,
                             Real3{0.31, -0.27, 0.19}, 0.0);
  for (int step = 0; step < 12; ++step) {
    solver.advance_ssprk2(storage, scalar, stage, 0.01);
  }
  double maximum_difference = 0.0;
  for (int k = 0; k < cells.z; ++k) {
    for (int j = 0; j < cells.y; ++j) {
      for (int i = 0; i < cells.x; ++i) {
        maximum_difference = std::max(
            maximum_difference, std::abs(scalar_view(i, j, k, 0) - constant));
      }
    }
  }
  HUNDUN_CHECK(maximum_difference == 0.0);
  return maximum_difference;
}

void check_serial_numerics(const MpiContext &context) {
  const ConvergenceResult coarse = run_sine_resolution(context, 32);
  const ConvergenceResult medium = run_sine_resolution(context, 64);
  const ConvergenceResult fine = run_sine_resolution(context, 128);
  const double order_coarse_medium =
      std::log(coarse.error / medium.error) / std::log(2.0);
  const double order_medium_fine =
      std::log(medium.error / fine.error) / std::log(2.0);
  const double maximum_mass_error =
      std::max({coarse.relative_mass_error, medium.relative_mass_error,
                fine.relative_mass_error});
  HUNDUN_CHECK(order_coarse_medium > 1.70);
  HUNDUN_CHECK(order_medium_fine > 1.70);
  HUNDUN_CHECK(maximum_mass_error < 1.0e-12);

  const double constant_difference = check_constant_state(context);
  std::cout << "PASSIVE_SCALAR_SERIAL"
            << " error32=" << coarse.error << " error64=" << medium.error
            << " error128=" << fine.error
            << " order32_64=" << order_coarse_medium
            << " order64_128=" << order_medium_fine
            << " relative_mass_error=" << maximum_mass_error
            << " constant_max_difference=" << constant_difference << '\n';
}

void check_oracle_numerics(const MpiContext &context) {
  const double oracle_difference = check_independent_oracle(context);
  std::cout << "PASSIVE_SCALAR_ORACLE"
            << " max_difference=" << oracle_difference << '\n';
}

} // namespace

int main(int argc, char **argv) {
  return hundun::test::run([&] {
    HUNDUN_CHECK(argc == 2);
    const std::string mode(argv[1]);
    if (mode == "limiter") {
      check_limiter();
      return;
    }
    HUNDUN_CHECK(mode == "serial" || mode == "oracle");
    MpiEnvironment environment(argc, argv);
    auto context = MpiContext::duplicate(MPI_COMM_SELF);
    check_limiter();
    if (mode == "oracle") {
      check_oracle_numerics(context);
      return;
    }
    check_serial_numerics(context);
  });
}
