// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/execution/execution.hpp"
#include "hundun/linear/linear_system.hpp"
#include "hundun/mesh/mesh_geometry.hpp"
#include "hundun/mesh/mesh_topology.hpp"
#include "hundun/runtime/error.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/structured_decomposition.hpp"

#include <cstdint>
#include <memory>
#include <optional>
#include <string>

namespace hundun::finite_volume {

enum class PressureConstraintMode : std::uint8_t {
  constant_nullspace,
  pressure_reference_patch
};

enum class PoissonSolverFamily : std::uint8_t {
  conjugate_gradient,
  bicgstab
};

struct PoissonBoundarySpec final {
  PressureConstraintMode mode{PressureConstraintMode::constant_nullspace};
  std::optional<std::uint32_t> pressure_reference_patch_id{};
};

struct PoissonCoefficientReplacementResult final {
  bool accepted{};
  bool changed{};
  int lowest_failing_rank{-1};
  std::uint64_t revision{};
};

// Construction performs a collective preflight.  This public error preserves
// the lowest rank selected by that preflight without requiring callers to
// interpret diagnostic text.
class PoissonConstructionError final : public runtime::Error {
 public:
  PoissonConstructionError(int failing_rank, std::string message);
  int failing_rank() const noexcept;

 private:
  int failing_rank_;
};

class MatrixFreePoissonOperator final : public linear::LinearOperator {
 public:
  static MatrixFreePoissonOperator create(
      const runtime::StructuredDecomposition& decomposition,
      const mesh::MeshTopology& topology,
      const mesh::MeshGeometry& geometry,
      execution::ExecutionContext& execution_context,
      execution::VectorView<const double> gamma_by_local_face,
      PoissonBoundarySpec boundary);

  ~MatrixFreePoissonOperator() noexcept override;
  MatrixFreePoissonOperator(MatrixFreePoissonOperator&&) noexcept;
  MatrixFreePoissonOperator& operator=(MatrixFreePoissonOperator&&) = delete;
  MatrixFreePoissonOperator(const MatrixFreePoissonOperator&) = delete;
  MatrixFreePoissonOperator& operator=(const MatrixFreePoissonOperator&) =
      delete;

  linear::VectorLayout domain_layout() const override;
  linear::VectorLayout range_layout() const override;
  const execution::ExecutionContext& context() const override;
  std::uint64_t revision() const override;
  execution::ExecutionEvent apply(
      execution::VectorView<const double> x,
      execution::VectorView<double> y) const override;
  bool has_diagonal() const override;
  execution::ExecutionEvent diagonal(
      execution::VectorView<double> d) const override;

  PressureConstraintMode constraint_mode() const noexcept;
  std::optional<std::uint32_t> pressure_reference_patch_id() const noexcept;
  PoissonSolverFamily solver_family() const noexcept;

  // Local transactional replacement. Failure preserves the active
  // coefficients, diagonal and revision. Distributed callers that require
  // one common commit must use collectively_replace_face_coefficients().
  std::uint64_t replace_face_coefficients(
      execution::VectorView<const double> gamma_by_local_face);

  // Collective prepare/commit over mpi_context. No rank changes active state
  // until every rank has prepared a valid candidate and any required revision
  // increment has been validated.
  // A rejected result preserves the active coefficients, diagonal and
  // revision everywhere and carries the lowest failing rank. MPI operation
  // failures remain typed exceptions from the runtime layer.
  PoissonCoefficientReplacementResult collectively_replace_face_coefficients(
      execution::VectorView<const double> gamma_by_local_face,
      const runtime::MpiContext& mpi_context);

  execution::ExecutionEvent accumulate_explicit_nonorthogonal_rhs(
      execution::VectorView<const double> gradient_x_by_local_face,
      execution::VectorView<const double> gradient_y_by_local_face,
      execution::VectorView<const double> gradient_z_by_local_face,
      execution::VectorView<double> rhs) const;

 private:
  struct Impl;
  explicit MatrixFreePoissonOperator(std::unique_ptr<Impl>) noexcept;
  std::unique_ptr<Impl> impl_;
};

class PoissonConstraint final {
 public:
  static PoissonConstraint create(
      const mesh::MeshTopology& topology,
      const mesh::MeshGeometry& geometry,
      execution::ExecutionContext& execution_context,
      const runtime::MpiContext& mpi_context,
      PressureConstraintMode mode);

  ~PoissonConstraint() noexcept;
  PoissonConstraint(PoissonConstraint&&) noexcept;
  PoissonConstraint& operator=(PoissonConstraint&&) = delete;
  PoissonConstraint(const PoissonConstraint&) = delete;
  PoissonConstraint& operator=(const PoissonConstraint&) = delete;

  PressureConstraintMode mode() const noexcept;
  void project_rhs(execution::VectorView<double> rhs) const;
  void normalize_solution(execution::VectorView<double> x) const;

 private:
  struct Impl;
  explicit PoissonConstraint(std::unique_ptr<Impl>) noexcept;
  std::unique_ptr<Impl> impl_;
};

}  // namespace hundun::finite_volume
