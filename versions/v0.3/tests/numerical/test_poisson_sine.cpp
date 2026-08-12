// SPDX-License-Identifier: Apache-2.0

#include "hundun/exec_execution.hpp"
#include "hundun/fvm_matrix_free_poisson.hpp"
#include "hundun/lin_conjugate_gradient.hpp"
#include "hundun/lin_preconditioners.hpp"
#include "hundun/mesh_geometry.hpp"
#include "hundun/mesh_topology.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "hundun/rt_structured_decomposition.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>

namespace {

using hundun::execution::Buffer;
using hundun::execution::CpuReferenceContext;
using hundun::finite_volume::MatrixFreePoissonOperator;
using hundun::finite_volume::PoissonConstraint;
using hundun::finite_volume::PoissonSolverFamily;
using hundun::finite_volume::PressureConstraintMode;
using hundun::linear::ConjugateGradientSolver;
using hundun::linear::JacobiPreconditioner;
using hundun::linear::SolveControl;
using hundun::linear::SolveTerminationReason;
using hundun::mesh::LocalCellId;
using hundun::mesh::LocalFaceId;
using hundun::mesh::MeshGeometry;
using hundun::mesh::MeshTopology;
using hundun::mesh::UniformBoxMapping;
using hundun::runtime::DecompositionOptions;
using hundun::runtime::Int3;
using hundun::runtime::MpiContext;
using hundun::runtime::MpiEnvironment;
using hundun::runtime::StructuredDecomposition;

constexpr double kPi = 3.141592653589793238462643383279502884;

std::size_t bytes_for(std::size_t count) {
  HUNDUN_CHECK(count <=
               std::numeric_limits<std::size_t>::max() / sizeof(double));
  return count * sizeof(double);
}

struct GridResult final {
  double l2_error{};
  double absolute_mean{};
  double residual{};
  std::uint64_t iterations{};
};

GridResult run_grid(const MpiContext& mpi, int cells) {
  CpuReferenceContext execution;
  const Int3 extent{cells, cells, cells};
  auto decomposition = StructuredDecomposition::create(
      mpi, extent, {true, true, true},
      DecompositionOptions{Int3{mpi.size(), 1, 1}});
  const MeshTopology topology(decomposition);
  const MeshGeometry geometry(
      topology,
      UniformBoxMapping({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}));
  Buffer gamma(execution, bytes_for(topology.local_face_count()));
  auto gamma_view = gamma.view(0U, topology.local_face_count());
  for (LocalFaceId face = 0; face < topology.local_face_count(); ++face) {
    gamma_view[face] = 1.0;
  }
  auto linear_operator = MatrixFreePoissonOperator::create(
      decomposition, topology, geometry, execution, gamma_view, {});
  HUNDUN_CHECK(linear_operator.solver_family() ==
               PoissonSolverFamily::conjugate_gradient);
  auto constraint = PoissonConstraint::create(
      topology, geometry, execution, mpi,
      PressureConstraintMode::constant_nullspace);

  const std::size_t owned_count = topology.owned_cell_count();
  Buffer rhs_buffer(execution, bytes_for(owned_count));
  Buffer solution_buffer(execution, bytes_for(owned_count));
  Buffer residual_buffer(execution, bytes_for(owned_count));
  auto rhs = rhs_buffer.view(0U, owned_count);
  auto solution = solution_buffer.view(0U, owned_count);
  auto residual = residual_buffer.view(0U, owned_count);
  for (LocalCellId cell = 0; cell < owned_count; ++cell) {
    const auto global = topology.global_cell(cell);
    const double x =
        (static_cast<double>(global.x) + 0.5) / static_cast<double>(cells);
    const double y =
        (static_cast<double>(global.y) + 0.5) / static_cast<double>(cells);
    const double z =
        (static_cast<double>(global.z) + 0.5) / static_cast<double>(cells);
    const double exact = std::sin(2.0 * kPi * x) *
                         std::sin(2.0 * kPi * y) *
                         std::sin(2.0 * kPi * z);
    rhs[cell] = 12.0 * kPi * kPi * exact;
    solution[cell] = 0.0;
  }
  constraint.project_rhs(rhs);
  constraint.normalize_solution(solution);

  JacobiPreconditioner jacobi(execution);
  jacobi.update(linear_operator, linear_operator.revision());
  ConjugateGradientSolver solver(execution, mpi);
  const SolveControl control{1.0e-12, 1.0e-10, 500U, 20U};
  const auto report =
      solver.solve(linear_operator, jacobi, rhs, solution, control);
  HUNDUN_CHECK(report.reason == SolveTerminationReason::converged);
  HUNDUN_CHECK(report.iterations > 0U);
  HUNDUN_CHECK(report.matvec_count >= report.iterations + 2U);
  HUNDUN_CHECK(report.preconditioner_apply_count >= report.iterations);
  HUNDUN_CHECK(report.global_reduction_count > 0U);
  constraint.normalize_solution(solution);

  linear_operator.apply(solution, residual).wait();
  double local_residual_square = 0.0;
  double local_rhs_square = 0.0;
  double local_error_integral = 0.0;
  double local_weighted_solution = 0.0;
  double local_volume = 0.0;
  for (LocalCellId cell = 0; cell < owned_count; ++cell) {
    const auto global = topology.global_cell(cell);
    const double x =
        (static_cast<double>(global.x) + 0.5) / static_cast<double>(cells);
    const double y =
        (static_cast<double>(global.y) + 0.5) / static_cast<double>(cells);
    const double z =
        (static_cast<double>(global.z) + 0.5) / static_cast<double>(cells);
    const double exact = std::sin(2.0 * kPi * x) *
                         std::sin(2.0 * kPi * y) *
                         std::sin(2.0 * kPi * z);
    const double difference = residual[cell] - rhs[cell];
    const double error = solution[cell] - exact;
    const double volume = geometry.cell_volume_m3(cell);
    local_residual_square += difference * difference;
    local_rhs_square += rhs[cell] * rhs[cell];
    local_error_integral += volume * error * error;
    local_weighted_solution += volume * solution[cell];
    local_volume += volume;
  }
  double values[5]{local_residual_square, local_rhs_square,
                   local_error_integral, local_weighted_solution,
                   local_volume};
  HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, values, 5, MPI_DOUBLE, MPI_SUM,
                             mpi.comm()) == MPI_SUCCESS);
  const double independent_residual = std::sqrt(values[0]);
  const double rhs_norm = std::sqrt(values[1]);
  const double tolerance = std::max(control.atol, control.rtol * rhs_norm);
  HUNDUN_CHECK(independent_residual <= tolerance);
  HUNDUN_CHECK(std::abs(report.final_residual - independent_residual) <=
               64.0 * std::numeric_limits<double>::epsilon() *
                   std::max(1.0, independent_residual));
  const double mean = values[3] / values[4];
  HUNDUN_CHECK(std::abs(mean) <= 1.0e-12);
  const double l2_error = std::sqrt(values[2]);
  HUNDUN_CHECK(std::isfinite(l2_error));
  HUNDUN_CHECK(l2_error > 0.0);
  if (mpi.rank() == 0) {
    std::cout << "POISSON_SINE cells=" << cells << " l2=" << l2_error
              << " mean=" << std::abs(mean)
              << " residual=" << independent_residual
              << " iterations=" << report.iterations << '\n';
  }
  return {l2_error, std::abs(mean), independent_residual,
          report.iterations};
}

void run_tests(const MpiContext& mpi) {
  const GridResult coarse = run_grid(mpi, 16);
  const GridResult medium = run_grid(mpi, 32);
  const GridResult fine = run_grid(mpi, 64);
  const double first_order = std::log(coarse.l2_error / medium.l2_error) /
                             std::log(2.0);
  const double second_order = std::log(medium.l2_error / fine.l2_error) /
                              std::log(2.0);
  HUNDUN_CHECK(first_order >= 1.8);
  HUNDUN_CHECK(second_order >= 1.8);
  HUNDUN_CHECK(coarse.absolute_mean <= 1.0e-12);
  HUNDUN_CHECK(medium.absolute_mean <= 1.0e-12);
  HUNDUN_CHECK(fine.absolute_mean <= 1.0e-12);
  if (mpi.rank() == 0) {
    std::cout << "POISSON_SINE_ORDER first=" << first_order
              << " second=" << second_order << '\n';
  }
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
