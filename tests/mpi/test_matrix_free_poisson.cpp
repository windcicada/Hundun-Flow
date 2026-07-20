// SPDX-License-Identifier: Apache-2.0

#include "hundun/execution/execution.hpp"
#include "hundun/finite_volume/matrix_free_poisson.hpp"
#include "hundun/linear/ghosted_vector_halo.hpp"
#include "hundun/linear/bicgstab.hpp"
#include "hundun/linear/conjugate_gradient.hpp"
#include "hundun/linear/preconditioners.hpp"
#include "hundun/mesh/mesh_geometry.hpp"
#include "hundun/mesh/mesh_topology.hpp"
#include "hundun/runtime/error.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/mpi_environment.hpp"
#include "hundun/runtime/structured_decomposition.hpp"
#include "finite_volume/src/matrix_free_poisson_detail.hpp"
#include "tests/support/allocation_attempt_guard.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <iostream>
#include <string>
#include <vector>

namespace {

using hundun::execution::Buffer;
using hundun::execution::CpuReferenceContext;
using hundun::finite_volume::MatrixFreePoissonOperator;
using hundun::finite_volume::PoissonBoundarySpec;
using hundun::finite_volume::PoissonConstraint;
using hundun::finite_volume::PoissonSolverFamily;
using hundun::finite_volume::PressureConstraintMode;
using hundun::linear::BiCGStabSolver;
using hundun::linear::ConjugateGradientSolver;
using hundun::linear::JacobiPreconditioner;
using hundun::linear::SolveControl;
using hundun::linear::SolveReport;
using hundun::linear::SolveTerminationReason;
using hundun::mesh::EntityOwnership;
using hundun::mesh::AnalyticWarpedBoxMapping;
using hundun::mesh::FaceSide;
using hundun::mesh::LocalCellId;
using hundun::mesh::LocalFaceId;
using hundun::mesh::MeshGeometry;
using hundun::mesh::MeshTopology;
using hundun::mesh::UniformBoxMapping;
using hundun::runtime::DecompositionOptions;
using hundun::runtime::Int3;
using hundun::runtime::MpiContext;
using hundun::runtime::MpiEnvironment;
using hundun::runtime::Real3;
using hundun::runtime::StructuredDecomposition;

constexpr Int3 kExtent{5, 4, 3};
constexpr Real3 kOrigin{-0.5, 0.25, 1.5};
constexpr Real3 kLength{2.0, 3.0, 4.5};
namespace poisson_detail = hundun::finite_volume::detail;
namespace allocation_probe = hundun::test::allocation_probe;

double dot(Real3 left, Real3 right) {
  return left.x * right.x + left.y * right.y + left.z * right.z;
}

double norm_squared(Real3 value) { return dot(value, value); }

std::uint64_t bits(double value) {
  std::uint64_t result = 0U;
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

Int3 process_grid_for(int ranks) {
  switch (ranks) {
    case 1:
      return {1, 1, 1};
    case 2:
      return {2, 1, 1};
    case 4:
      return {2, 2, 1};
    default:
      throw hundun::runtime::Error("unsupported Poisson test rank count");
  }
}

std::size_t bytes_for(std::size_t count) {
  HUNDUN_CHECK(count <=
               std::numeric_limits<std::size_t>::max() / sizeof(double));
  return count * sizeof(double);
}

template <class Function>
void expect_error(Function&& function) {
  bool rejected = false;
  try {
    function();
  } catch (const hundun::runtime::Error& error) {
    rejected = !std::string(error.what()).empty();
  }
  HUNDUN_CHECK(rejected);
}

std::vector<double> face_coefficients(const MeshTopology& topology,
                                      double scale) {
  std::vector<double> result(topology.local_face_count());
  for (LocalFaceId face = 0; face < topology.local_face_count(); ++face) {
    const auto global_id = topology.global_face_id(face);
    const auto periodic_pair = topology.periodic_pair(face);
    const auto coefficient_id =
        periodic_pair.has_value() ? std::min(global_id, *periodic_pair)
                                  : global_id;
    result[face] = 1.0 +
                   scale * static_cast<double>(coefficient_id + 1U);
  }
  return result;
}

void copy_to_buffer(const std::vector<double>& values, Buffer& buffer) {
  auto view = buffer.view(0U, values.size());
  for (std::size_t index = 0; index < values.size(); ++index) {
    view[index] = values[index];
  }
}

std::vector<double> owned_values(const MeshTopology& topology, double shift) {
  std::vector<double> result(topology.owned_cell_count());
  for (LocalCellId cell = 0; cell < topology.owned_cell_count(); ++cell) {
    const double id = static_cast<double>(topology.global_cell_id(cell));
    result[cell] = shift + 0.125 * id + 0.015625 * id * id;
  }
  return result;
}

struct OracleResult final {
  std::vector<double> applied;
  std::vector<double> diagonal;
};

OracleResult independent_oracle(const MeshTopology& topology,
                                const MeshGeometry& geometry,
                                const std::vector<double>& gamma,
                                const std::vector<double>& local_values,
                                PoissonBoundarySpec boundary) {
  OracleResult result{std::vector<double>(topology.owned_cell_count(), 0.0),
                      std::vector<double>(topology.owned_cell_count(), 0.0)};
  for (LocalFaceId face = 0; face < topology.local_face_count(); ++face) {
    const LocalCellId owner = topology.owner(face);
    const auto neighbour = topology.neighbour(face);
    const bool owner_owned =
        topology.cell_ownership(owner) == EntityOwnership::owned;
    const bool neighbour_owned =
        neighbour.has_value() && topology.cell_ownership(*neighbour) ==
                                     EntityOwnership::owned;
    const bool periodic = topology.periodic_pair(face).has_value();
    const auto patch = topology.patch_id(face);
    const bool reference =
        patch.has_value() &&
        boundary.mode == PressureConstraintMode::pressure_reference_patch &&
        boundary.pressure_reference_patch_id == patch;
    if (!neighbour.has_value() && !reference) {
      continue;
    }
    const Real3 area =
        geometry.face_area_vector_m2(face, FaceSide::owner);
    const Real3 displacement = geometry.face_displacement_m(face);
    const double coefficient =
        gamma[face] * norm_squared(area) / dot(area, displacement);
    if (reference) {
      if (owner_owned) {
        const double volume = geometry.cell_volume_m3(owner);
        result.applied[owner] += coefficient * local_values[owner] / volume;
        result.diagonal[owner] += coefficient / volume;
      }
      continue;
    }
    HUNDUN_CHECK(neighbour.has_value());
    if (periodic) {
      if (owner_owned) {
        const double volume = geometry.cell_volume_m3(owner);
        result.applied[owner] +=
            coefficient * (local_values[owner] - local_values[*neighbour]) /
            volume;
        if (owner != *neighbour) {
          result.diagonal[owner] += coefficient / volume;
        }
      }
      continue;
    }
    if (owner_owned) {
      const double volume = geometry.cell_volume_m3(owner);
      result.applied[owner] +=
          coefficient * (local_values[owner] - local_values[*neighbour]) /
          volume;
      result.diagonal[owner] += coefficient / volume;
    }
    if (neighbour_owned) {
      const double volume = geometry.cell_volume_m3(*neighbour);
      result.applied[*neighbour] +=
          coefficient *
          (local_values[*neighbour] - local_values[owner]) / volume;
      result.diagonal[*neighbour] += coefficient / volume;
    }
  }
  return result;
}

void check_near(double actual, double expected, double factor = 256.0) {
  const double scale = std::max({1.0, std::abs(actual), std::abs(expected)});
  HUNDUN_CHECK(std::abs(actual - expected) <=
               factor * std::numeric_limits<double>::epsilon() * scale);
}

void run_operator_case(const MpiContext& mpi, PoissonBoundarySpec boundary,
                       bool periodic_x) {
  CpuReferenceContext execution;
  auto decomposition = StructuredDecomposition::create(
      mpi, kExtent, {periodic_x, false, false},
      DecompositionOptions{process_grid_for(mpi.size())});
  const MeshTopology topology(decomposition);
  const MeshGeometry geometry(topology, UniformBoxMapping(kOrigin, kLength));
  std::vector<double> gamma = face_coefficients(topology, 0.0025);
  Buffer gamma_buffer(execution, bytes_for(gamma.size()));
  copy_to_buffer(gamma, gamma_buffer);

  auto linear_operator = MatrixFreePoissonOperator::create(
      decomposition, topology, geometry, execution,
      gamma_buffer.view(0U, gamma.size()), boundary);
  HUNDUN_CHECK(linear_operator.revision() != 0U);
  HUNDUN_CHECK(linear_operator.domain_layout() ==
               hundun::linear::VectorLayout::from_topology(topology));
  HUNDUN_CHECK(linear_operator.range_layout() ==
               linear_operator.domain_layout());
  HUNDUN_CHECK(&linear_operator.context() == &execution);
  HUNDUN_CHECK(linear_operator.has_diagonal());
  HUNDUN_CHECK(linear_operator.constraint_mode() == boundary.mode);
  HUNDUN_CHECK(linear_operator.pressure_reference_patch_id() ==
               boundary.pressure_reference_patch_id);
  HUNDUN_CHECK(linear_operator.solver_family() ==
               PoissonSolverFamily::conjugate_gradient);

  const std::vector<double> owned = owned_values(topology, -0.75);
  Buffer x_buffer(execution, bytes_for(owned.size()));
  Buffer y_buffer(execution, bytes_for(owned.size()));
  Buffer diagonal_buffer(execution, bytes_for(owned.size()));
  copy_to_buffer(owned, x_buffer);
  auto y = y_buffer.view(0U, owned.size());
  auto diagonal = diagonal_buffer.view(0U, owned.size());
  for (std::size_t index = 0; index < owned.size(); ++index) {
    y[index] = -777.0;
    diagonal[index] = -888.0;
  }
  auto apply_event =
      linear_operator.apply(x_buffer.view(0U, owned.size()), y);
  HUNDUN_CHECK(apply_event.ready());
  apply_event.wait();
  auto diagonal_event = linear_operator.diagonal(diagonal);
  HUNDUN_CHECK(diagonal_event.ready());
  diagonal_event.wait();

  hundun::linear::GhostedVector oracle_values(
      execution, hundun::linear::VectorLayout::from_topology(topology));
  auto oracle_local = oracle_values.local_view();
  for (std::size_t index = 0; index < owned.size(); ++index) {
    oracle_local[index] = owned[index];
  }
  auto oracle_halo = hundun::linear::GhostedVectorHalo::create(
      decomposition, topology, execution);
  oracle_halo.exchange(oracle_values);
  std::vector<double> local_values(topology.local_cell_count());
  for (std::size_t index = 0; index < local_values.size(); ++index) {
    local_values[index] = oracle_local[index];
  }
  const OracleResult oracle = independent_oracle(
      topology, geometry, gamma, local_values, boundary);
  for (std::size_t index = 0; index < owned.size(); ++index) {
    check_near(y[index], oracle.applied[index]);
    check_near(diagonal[index], oracle.diagonal[index]);
    HUNDUN_CHECK(diagonal[index] >= 0.0);
  }

  if (boundary.mode == PressureConstraintMode::constant_nullspace &&
      periodic_x) {
    poisson_detail::set_poisson_test_options(
        {{}, false, true});
    poisson_detail::reset_poisson_apply_trace();
    {
      allocation_probe::AllocationAttemptGuard allocation_guard;
      linear_operator.apply(x_buffer.view(0U, owned.size()), y).wait();
      HUNDUN_CHECK(allocation_guard.attempts() == 0U);
    }
    const auto trace = poisson_detail::poisson_apply_trace();
    HUNDUN_CHECK(trace.count == 4U);
    HUNDUN_CHECK(trace.phases[0] ==
                 poisson_detail::PoissonApplyPhase::halo_begin);
    HUNDUN_CHECK(trace.phases[1] ==
                 poisson_detail::PoissonApplyPhase::local_rows);
    HUNDUN_CHECK(trace.phases[2] ==
                 poisson_detail::PoissonApplyPhase::halo_wait);
    HUNDUN_CHECK(trace.phases[3] ==
                 poisson_detail::PoissonApplyPhase::partition_rows);

    for (std::size_t index = 0; index < owned.size(); ++index) {
      y[index] = -919.25;
    }
    poisson_detail::set_poisson_test_options(
        {poisson_detail::PoissonApplyPhase::local_rows, false, true});
    poisson_detail::reset_poisson_apply_trace();
    expect_error([&] {
      linear_operator.apply(x_buffer.view(0U, owned.size()), y).wait();
    });
    const auto failed_trace = poisson_detail::poisson_apply_trace();
    HUNDUN_CHECK(failed_trace.count == 3U);
    HUNDUN_CHECK(failed_trace.phases[0] ==
                 poisson_detail::PoissonApplyPhase::halo_begin);
    HUNDUN_CHECK(failed_trace.phases[1] ==
                 poisson_detail::PoissonApplyPhase::local_rows);
    HUNDUN_CHECK(failed_trace.phases[2] ==
                 poisson_detail::PoissonApplyPhase::halo_wait);
    for (std::size_t index = 0; index < owned.size(); ++index) {
      HUNDUN_CHECK(y[index] == -919.25);
    }
    poisson_detail::set_poisson_test_options({});
    {
      allocation_probe::AllocationAttemptGuard allocation_guard;
      linear_operator.diagonal(diagonal).wait();
      HUNDUN_CHECK(allocation_guard.attempts() == 0U);
    }

    for (std::size_t index = 0; index < owned.size(); ++index) {
      y[index] = -451.75;
    }
    poisson_detail::set_poisson_test_options(
        {{}, false, false, true});
    expect_error([&] {
      linear_operator.apply(x_buffer.view(0U, owned.size()), y).wait();
    });
    poisson_detail::set_poisson_test_options({});
    for (std::size_t index = 0; index < owned.size(); ++index) {
      HUNDUN_CHECK(y[index] == -451.75);
    }
  }

  if (boundary.mode == PressureConstraintMode::constant_nullspace) {
    for (std::size_t index = 0; index < owned.size(); ++index) {
      x_buffer.view(0U, owned.size())[index] = 3.25;
    }
    linear_operator.apply(x_buffer.view(0U, owned.size()), y).wait();
    for (std::size_t index = 0; index < owned.size(); ++index) {
      check_near(y[index], 0.0, 1024.0);
    }
  }

  const std::uint64_t initial_revision = linear_operator.revision();
  const std::vector<double> replacement = face_coefficients(topology, 0.005);
  copy_to_buffer(replacement, gamma_buffer);
  HUNDUN_CHECK(linear_operator.replace_face_coefficients(
                   gamma_buffer.view(0U, replacement.size())) ==
               initial_revision + 1U);
  HUNDUN_CHECK(linear_operator.revision() == initial_revision + 1U);
  HUNDUN_CHECK(linear_operator.replace_face_coefficients(
                   gamma_buffer.view(0U, replacement.size())) ==
               initial_revision + 2U);

  const std::uint64_t before_failure = linear_operator.revision();
  gamma_buffer.view(0U, replacement.size())[0] =
      std::numeric_limits<double>::quiet_NaN();
  expect_error([&] {
    static_cast<void>(linear_operator.replace_face_coefficients(
        gamma_buffer.view(0U, replacement.size())));
  });
  HUNDUN_CHECK(linear_operator.revision() == before_failure);
  gamma_buffer.view(0U, replacement.size())[0] = replacement[0];
  poisson_detail::set_poisson_test_options({{}, true, false});
  expect_error([&] {
    static_cast<void>(linear_operator.replace_face_coefficients(
        gamma_buffer.view(0U, replacement.size())));
  });
  poisson_detail::set_poisson_test_options({});
  HUNDUN_CHECK(linear_operator.revision() == before_failure);
}

void test_construction_rejections(const MpiContext& mpi) {
  CpuReferenceContext execution;
  auto decomposition = StructuredDecomposition::create(
      mpi, kExtent, {false, false, false},
      DecompositionOptions{process_grid_for(mpi.size())});
  const MeshTopology topology(decomposition);
  const MeshGeometry geometry(topology, UniformBoxMapping(kOrigin, kLength));
  const std::vector<double> gamma = face_coefficients(topology, 0.001);
  Buffer gamma_buffer(execution, bytes_for(gamma.size()));
  copy_to_buffer(gamma, gamma_buffer);
  expect_error([&] {
    static_cast<void>(MatrixFreePoissonOperator::create(
        decomposition, topology, geometry, execution,
        gamma_buffer.view(0U, gamma.size()),
        {PressureConstraintMode::constant_nullspace, 0U}));
  });
  expect_error([&] {
    static_cast<void>(MatrixFreePoissonOperator::create(
        decomposition, topology, geometry, execution,
        gamma_buffer.view(0U, gamma.size()),
        {PressureConstraintMode::pressure_reference_patch, std::nullopt}));
  });
  expect_error([&] {
    static_cast<void>(MatrixFreePoissonOperator::create(
        decomposition, topology, geometry, execution,
        gamma_buffer.view(0U, gamma.size() - 1U), {}));
  });
}

void test_extent_one_periodic_operator(const MpiContext& mpi) {
  if (mpi.size() != 1) {
    return;
  }
  CpuReferenceContext execution;
  auto decomposition = StructuredDecomposition::create(
      mpi, {1, 1, 1}, {true, true, true},
      DecompositionOptions{Int3{1, 1, 1}});
  const MeshTopology topology(decomposition);
  const MeshGeometry geometry(
      topology, UniformBoxMapping({0.0, 0.0, 0.0}, {2.0, 3.0, 4.0}));
  std::vector<double> gamma(topology.local_face_count(), 1.0);
  Buffer gamma_buffer(execution, bytes_for(gamma.size()));
  copy_to_buffer(gamma, gamma_buffer);
  auto linear_operator = MatrixFreePoissonOperator::create(
      decomposition, topology, geometry, execution,
      gamma_buffer.view(0U, gamma.size()), {});
  Buffer input(execution, sizeof(double));
  Buffer output(execution, sizeof(double));
  Buffer diagonal(execution, sizeof(double));
  input.view(0U, 1U)[0] = 7.0;
  output.view(0U, 1U)[0] = -1.0;
  diagonal.view(0U, 1U)[0] = -1.0;
  linear_operator.apply(input.view(0U, 1U), output.view(0U, 1U)).wait();
  linear_operator.diagonal(diagonal.view(0U, 1U)).wait();
  HUNDUN_CHECK(output.view(0U, 1U)[0] == 0.0);
  HUNDUN_CHECK(diagonal.view(0U, 1U)[0] == 0.0);
}

double global_volume_weighted_mean(const MpiContext& mpi,
                                   const MeshTopology& topology,
                                   const MeshGeometry& geometry,
                                   hundun::execution::VectorView<double> values) {
  double sums[2]{};
  for (LocalCellId cell = 0; cell < topology.owned_cell_count(); ++cell) {
    const double volume = geometry.cell_volume_m3(cell);
    sums[0] += volume * values[cell];
    sums[1] += volume;
  }
  double global[2]{};
  HUNDUN_CHECK(MPI_Allreduce(sums, global, 2, MPI_DOUBLE, MPI_SUM,
                             mpi.comm()) == MPI_SUCCESS);
  return global[0] / global[1];
}

void test_constraint_case(const MpiContext& mpi, bool warped) {
  CpuReferenceContext execution;
  auto decomposition = StructuredDecomposition::create(
      mpi, kExtent, {true, true, true},
      DecompositionOptions{process_grid_for(mpi.size())});
  const MeshTopology topology(decomposition);
  const MeshGeometry geometry =
      warped
          ? MeshGeometry(topology,
                         AnalyticWarpedBoxMapping(
                             kOrigin, kLength, {0.02, -0.015, 0.01}))
          : MeshGeometry(topology, UniformBoxMapping(kOrigin, kLength));
  auto constraint = PoissonConstraint::create(
      topology, geometry, execution, mpi,
      PressureConstraintMode::constant_nullspace);
  HUNDUN_CHECK(constraint.mode() ==
               PressureConstraintMode::constant_nullspace);
  Buffer values_buffer(execution, bytes_for(topology.owned_cell_count()));
  auto values = values_buffer.view(0U, topology.owned_cell_count());
  for (LocalCellId cell = 0; cell < topology.owned_cell_count(); ++cell) {
    values[cell] = 4.0 +
                   0.03125 * static_cast<double>(topology.global_cell_id(cell));
  }
  const auto before = mpi.fp64_reduction_counters();
  constraint.project_rhs(values);
  const auto after_project = mpi.fp64_reduction_counters();
  HUNDUN_CHECK(after_project.collective_calls - before.collective_calls == 2U);
  HUNDUN_CHECK(after_project.reduced_scalars - before.reduced_scalars == 3U);
  HUNDUN_CHECK(after_project.logical_payload_bytes -
                   before.logical_payload_bytes ==
               3U * sizeof(double));
  check_near(global_volume_weighted_mean(mpi, topology, geometry, values),
             0.0, 4096.0);

  for (LocalCellId cell = 0; cell < topology.owned_cell_count(); ++cell) {
    values[cell] += 7.25;
  }
  constraint.normalize_solution(values);
  check_near(global_volume_weighted_mean(mpi, topology, geometry, values),
             0.0, 4096.0);

  std::vector<std::uint64_t> original_bits(topology.owned_cell_count());
  for (LocalCellId cell = 0; cell < topology.owned_cell_count(); ++cell) {
    values[cell] = -2.0 + 0.5 * static_cast<double>(cell);
    original_bits[cell] = bits(values[cell]);
  }
  const int failing_rank = mpi.size() > 1 ? 1 : 0;
  if (mpi.rank() == failing_rank) {
    values[0] = std::numeric_limits<double>::quiet_NaN();
    original_bits[0] = bits(values[0]);
  }
  expect_error([&] { constraint.project_rhs(values); });
  for (LocalCellId cell = 0; cell < topology.owned_cell_count(); ++cell) {
    HUNDUN_CHECK(bits(values[cell]) == original_bits[cell]);
  }

  auto reference = PoissonConstraint::create(
      topology, geometry, execution, mpi,
      PressureConstraintMode::pressure_reference_patch);
  for (LocalCellId cell = 0; cell < topology.owned_cell_count(); ++cell) {
    values[cell] = 100.0 + static_cast<double>(cell);
    original_bits[cell] = bits(values[cell]);
  }
  const auto before_noop = mpi.fp64_reduction_counters();
  reference.project_rhs(values);
  reference.normalize_solution(values);
  const auto after_noop = mpi.fp64_reduction_counters();
  HUNDUN_CHECK(after_noop.collective_calls == before_noop.collective_calls);
  HUNDUN_CHECK(after_noop.reduced_scalars == before_noop.reduced_scalars);
  for (LocalCellId cell = 0; cell < topology.owned_cell_count(); ++cell) {
    HUNDUN_CHECK(bits(values[cell]) == original_bits[cell]);
  }
  Buffer stale_buffer(execution, sizeof(double));
  auto stale = stale_buffer.view(0U, 1U);
  stale_buffer.reallocate(2U * sizeof(double));
  reference.project_rhs(stale);
  reference.normalize_solution(stale);
}

void test_collective_create_rejection(const MpiContext& mpi) {
  CpuReferenceContext execution;
  auto decomposition = StructuredDecomposition::create(
      mpi, kExtent, {true, false, false},
      DecompositionOptions{process_grid_for(mpi.size())});
  const MeshTopology topology(decomposition);
  const MeshGeometry geometry(topology, UniformBoxMapping(kOrigin, kLength));
  const std::vector<double> gamma = face_coefficients(topology, 0.001);
  Buffer gamma_buffer(execution, bytes_for(gamma.size()));
  copy_to_buffer(gamma, gamma_buffer);
  const int failing_rank = mpi.size() > 1 ? 1 : 0;
  if (mpi.rank() == failing_rank) {
    gamma_buffer.view(0U, gamma.size())[0] =
        std::numeric_limits<double>::quiet_NaN();
  }
  expect_error([&] {
    static_cast<void>(MatrixFreePoissonOperator::create(
        decomposition, topology, geometry, execution,
        gamma_buffer.view(0U, gamma.size()), {}));
  });
  mpi.barrier();
}

void test_warped_explicit_correction(const MpiContext& mpi) {
  CpuReferenceContext execution;
  auto decomposition = StructuredDecomposition::create(
      mpi, kExtent, {true, false, false},
      DecompositionOptions{process_grid_for(mpi.size())});
  const MeshTopology topology(decomposition);
  const MeshGeometry geometry(
      topology,
      AnalyticWarpedBoxMapping(kOrigin, kLength, {0.02, -0.015, 0.01}));
  const std::vector<double> gamma = face_coefficients(topology, 0.0015);
  Buffer gamma_buffer(execution, bytes_for(gamma.size()));
  copy_to_buffer(gamma, gamma_buffer);
  auto linear_operator = MatrixFreePoissonOperator::create(
      decomposition, topology, geometry, execution,
      gamma_buffer.view(0U, gamma.size()), {});
  HUNDUN_CHECK(linear_operator.solver_family() ==
               PoissonSolverFamily::bicgstab);

  Buffer gx_buffer(execution, bytes_for(gamma.size()));
  Buffer gy_buffer(execution, bytes_for(gamma.size()));
  Buffer gz_buffer(execution, bytes_for(gamma.size()));
  auto gx = gx_buffer.view(0U, gamma.size());
  auto gy = gy_buffer.view(0U, gamma.size());
  auto gz = gz_buffer.view(0U, gamma.size());
  for (LocalFaceId face = 0; face < topology.local_face_count(); ++face) {
    const double id = static_cast<double>(topology.global_face_id(face));
    gx[face] = 0.125 + 0.0005 * id;
    gy[face] = -0.25 + 0.00025 * id;
    gz[face] = 0.375 - 0.000125 * id;
  }
  Buffer rhs_buffer(execution, bytes_for(topology.owned_cell_count()));
  auto rhs = rhs_buffer.view(0U, topology.owned_cell_count());
  std::vector<double> expected(topology.owned_cell_count());
  for (LocalCellId cell = 0; cell < topology.owned_cell_count(); ++cell) {
    rhs[cell] = -1.0 + 0.0625 * static_cast<double>(cell);
    expected[cell] = rhs[cell];
  }
  for (LocalFaceId face = 0; face < topology.local_face_count(); ++face) {
    const LocalCellId owner = topology.owner(face);
    const auto neighbour = topology.neighbour(face);
    const bool owner_owned =
        topology.cell_ownership(owner) == EntityOwnership::owned;
    const bool neighbour_owned =
        neighbour.has_value() && topology.cell_ownership(*neighbour) ==
                                     EntityOwnership::owned;
    const bool periodic = topology.periodic_pair(face).has_value();
    if (!neighbour.has_value() || (periodic && !owner_owned) ||
        (!owner_owned && !neighbour_owned)) {
      continue;
    }
    const Real3 area =
        geometry.face_area_vector_m2(face, FaceSide::owner);
    const Real3 displacement = geometry.face_displacement_m(face);
    const double factor = norm_squared(area) / dot(area, displacement);
    const Real3 nonorthogonal{area.x - factor * displacement.x,
                              area.y - factor * displacement.y,
                              area.z - factor * displacement.z};
    const Real3 gradient{gx[face], gy[face], gz[face]};
    const double correction = gamma[face] * dot(gradient, nonorthogonal);
    if (owner_owned) {
      expected[owner] += correction / geometry.cell_volume_m3(owner);
    }
    if (!periodic && neighbour_owned) {
      expected[*neighbour] -=
          correction / geometry.cell_volume_m3(*neighbour);
    }
  }
  const std::uint64_t revision = linear_operator.revision();
  linear_operator
      .accumulate_explicit_nonorthogonal_rhs(gx, gy, gz, rhs)
      .wait();
  HUNDUN_CHECK(linear_operator.revision() == revision);
  for (LocalCellId cell = 0; cell < topology.owned_cell_count(); ++cell) {
    check_near(rhs[cell], expected[cell], 1024.0);
  }
  {
    allocation_probe::AllocationAttemptGuard allocation_guard;
    linear_operator
        .accumulate_explicit_nonorthogonal_rhs(gx, gy, gz, rhs)
        .wait();
    HUNDUN_CHECK(allocation_guard.attempts() == 0U);
  }

  std::vector<std::uint64_t> before_bits(topology.owned_cell_count());
  for (LocalCellId cell = 0; cell < topology.owned_cell_count(); ++cell) {
    rhs[cell] = 19.0 + static_cast<double>(cell);
    before_bits[cell] = bits(rhs[cell]);
  }
  gx[0] = std::numeric_limits<double>::quiet_NaN();
  expect_error([&] {
    linear_operator
        .accumulate_explicit_nonorthogonal_rhs(gx, gy, gz, rhs)
        .wait();
  });
  for (LocalCellId cell = 0; cell < topology.owned_cell_count(); ++cell) {
    HUNDUN_CHECK(bits(rhs[cell]) == before_bits[cell]);
  }
}

double independent_residual_norm(
    const MpiContext& mpi, MatrixFreePoissonOperator& linear_operator,
    hundun::execution::VectorView<const double> rhs,
    hundun::execution::VectorView<const double> solution,
    hundun::execution::VectorView<double> work) {
  linear_operator.apply(solution, work).wait();
  double local_square = 0.0;
  for (std::size_t cell = 0; cell < work.size(); ++cell) {
    const double residual = rhs[cell] - work[cell];
    local_square += residual * residual;
  }
  double global_square = 0.0;
  HUNDUN_CHECK(MPI_Allreduce(&local_square, &global_square, 1, MPI_DOUBLE,
                             MPI_SUM, mpi.comm()) == MPI_SUCCESS);
  return std::sqrt(global_square);
}

void check_solver_report(const SolveReport& report, double independent,
                         bool warped, int rank) {
  if (report.reason != SolveTerminationReason::converged) {
    std::cerr << "POISSON_SOLVER_REPORT rank=" << rank
              << " branch=" << (warped ? "warped" : "uniform")
              << " reason=" << static_cast<int>(report.reason)
              << " iterations=" << report.iterations
              << " initial=" << report.initial_residual
              << " recursive=" << report.recursive_residual
              << " final=" << report.final_residual
              << " independent=" << independent
              << " matvec=" << report.matvec_count
              << " preconditioner=" << report.preconditioner_apply_count
              << " reductions=" << report.global_reduction_count
              << " failing_rank=" << report.lowest_failing_rank << '\n';
  }
  HUNDUN_CHECK(report.reason == SolveTerminationReason::converged);
  HUNDUN_CHECK(report.iterations > 0U);
  HUNDUN_CHECK(report.matvec_count >= report.iterations + 2U);
  HUNDUN_CHECK(report.preconditioner_apply_count >= report.iterations);
  HUNDUN_CHECK(report.global_reduction_count > 0U);
  HUNDUN_CHECK(report.lowest_failing_rank == -1);
  const double tolerance =
      64.0 * std::numeric_limits<double>::epsilon() *
      std::max(1.0, independent);
  HUNDUN_CHECK(std::abs(report.final_residual - independent) <= tolerance);
}

void test_solver_integration(const MpiContext& mpi, bool warped) {
  CpuReferenceContext execution;
  constexpr Int3 extent{6, 5, 4};
  auto decomposition = StructuredDecomposition::create(
      mpi, extent, {true, true, true},
      DecompositionOptions{process_grid_for(mpi.size())});
  const MeshTopology topology(decomposition);
  const MeshGeometry geometry =
      warped
          ? MeshGeometry(topology,
                         AnalyticWarpedBoxMapping(
                             {0.0, 0.0, 0.0}, {1.0, 1.25, 0.75},
                             {0.02, -0.015, 0.01}))
          : MeshGeometry(topology,
                         UniformBoxMapping({0.0, 0.0, 0.0},
                                           {1.0, 1.25, 0.75}));
  const std::vector<double> gamma = face_coefficients(topology, 0.0001);
  Buffer gamma_buffer(execution, bytes_for(gamma.size()));
  copy_to_buffer(gamma, gamma_buffer);
  auto linear_operator = MatrixFreePoissonOperator::create(
      decomposition, topology, geometry, execution,
      gamma_buffer.view(0U, gamma.size()), {});
  HUNDUN_CHECK(linear_operator.solver_family() ==
               (warped ? PoissonSolverFamily::bicgstab
                       : PoissonSolverFamily::conjugate_gradient));
  auto constraint = PoissonConstraint::create(
      topology, geometry, execution, mpi,
      PressureConstraintMode::constant_nullspace);

  const std::size_t count = topology.owned_cell_count();
  Buffer exact_buffer(execution, bytes_for(count));
  Buffer rhs_buffer(execution, bytes_for(count));
  Buffer solution_buffer(execution, bytes_for(count));
  Buffer residual_buffer(execution, bytes_for(count));
  auto exact = exact_buffer.view(0U, count);
  auto rhs = rhs_buffer.view(0U, count);
  auto solution = solution_buffer.view(0U, count);
  auto residual = residual_buffer.view(0U, count);
  for (LocalCellId cell = 0; cell < count; ++cell) {
    const auto global = topology.global_cell(cell);
    exact[cell] =
        std::sin(2.0 * 3.14159265358979323846 *
                 (static_cast<double>(global.x) + 0.5) /
                 static_cast<double>(extent.x)) +
        0.25 * std::cos(2.0 * 3.14159265358979323846 *
                        (static_cast<double>(global.y) + 0.5) /
                        static_cast<double>(extent.y));
    solution[cell] = 0.0;
  }
  constraint.normalize_solution(exact);
  linear_operator.apply(exact, rhs).wait();
  constraint.project_rhs(rhs);
  constraint.normalize_solution(solution);

  JacobiPreconditioner jacobi(execution);
  jacobi.update(linear_operator, linear_operator.revision());
  const SolveControl control{1.0e-12, 1.0e-10, 500U, 20U};
  const auto before = mpi.fp64_reduction_counters();
  SolveReport report;
  if (warped) {
    BiCGStabSolver solver(execution, mpi);
    report = solver.solve(linear_operator, jacobi, rhs, solution, control);
  } else {
    ConjugateGradientSolver solver(execution, mpi);
    report = solver.solve(linear_operator, jacobi, rhs, solution, control);
  }
  const auto after = mpi.fp64_reduction_counters();
  HUNDUN_CHECK(after.collective_calls - before.collective_calls ==
               report.global_reduction_count);
  constraint.normalize_solution(solution);
  const double independent = independent_residual_norm(
      mpi, linear_operator, rhs, solution, residual);
  check_solver_report(report, independent, warped, mpi.rank());

  double local_max_error = 0.0;
  for (LocalCellId cell = 0; cell < count; ++cell) {
    local_max_error =
        std::max(local_max_error, std::abs(solution[cell] - exact[cell]));
  }
  double global_max_error = 0.0;
  HUNDUN_CHECK(MPI_Allreduce(&local_max_error, &global_max_error, 1,
                             MPI_DOUBLE, MPI_MAX, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(global_max_error <= 2.0e-8);

  const std::uint64_t old_revision = linear_operator.revision();
  auto replacement_gamma = gamma_buffer.view(0U, gamma.size());
  for (std::size_t face = 0; face < gamma.size(); ++face) {
    replacement_gamma[face] *= 1.125;
  }
  HUNDUN_CHECK(linear_operator.replace_face_coefficients(
                   gamma_buffer.view(0U, gamma.size())) ==
               old_revision + 1U);
  linear_operator.apply(exact, rhs).wait();
  constraint.project_rhs(rhs);
  expect_error([&] { jacobi.apply(rhs, residual).wait(); });
  for (LocalCellId cell = 0; cell < count; ++cell) {
    solution[cell] = 0.0;
  }
  SolveReport refreshed_report;
  if (warped) {
    BiCGStabSolver solver(execution, mpi);
    refreshed_report =
        solver.solve(linear_operator, jacobi, rhs, solution, control);
  } else {
    ConjugateGradientSolver solver(execution, mpi);
    refreshed_report =
        solver.solve(linear_operator, jacobi, rhs, solution, control);
  }
  constraint.normalize_solution(solution);
  const double refreshed_independent = independent_residual_norm(
      mpi, linear_operator, rhs, solution, residual);
  check_solver_report(refreshed_report, refreshed_independent, warped,
                      mpi.rank());
  local_max_error = 0.0;
  for (LocalCellId cell = 0; cell < count; ++cell) {
    local_max_error =
        std::max(local_max_error, std::abs(solution[cell] - exact[cell]));
  }
  global_max_error = 0.0;
  HUNDUN_CHECK(MPI_Allreduce(&local_max_error, &global_max_error, 1,
                             MPI_DOUBLE, MPI_MAX, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(global_max_error <= 2.0e-8);
}

void run_tests(const MpiContext& mpi) {
  run_operator_case(mpi, {}, true);
  run_operator_case(
      mpi,
      {PressureConstraintMode::pressure_reference_patch, std::uint32_t{0}},
      false);
  test_construction_rejections(mpi);
  test_extent_one_periodic_operator(mpi);
  test_collective_create_rejection(mpi);
  test_constraint_case(mpi, false);
  test_constraint_case(mpi, true);
  test_warped_explicit_correction(mpi);
  test_solver_integration(mpi, false);
  test_solver_integration(mpi, true);
  mpi.barrier();
}

}  // namespace

int main(int argc, char** argv) {
  MpiEnvironment environment(argc, argv);
  MpiContext mpi = MpiContext::duplicate(MPI_COMM_WORLD);
  const int result = hundun::test::run([&] { run_tests(mpi); });
  if (result != EXIT_SUCCESS) {
    MPI_Abort(mpi.comm(), result);
  }
  return result;
}
