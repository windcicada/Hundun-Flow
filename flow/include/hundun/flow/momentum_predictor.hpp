// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/boundary/basic_boundary.hpp"
#include "hundun/execution/execution.hpp"
#include "hundun/linear/linear_system.hpp"
#include "hundun/mesh/mesh_geometry.hpp"
#include "hundun/mesh/mesh_topology.hpp"
#include "hundun/runtime/field_access_plan.hpp"
#include "hundun/runtime/field_registry.hpp"
#include "hundun/runtime/field_storage.hpp"
#include "hundun/runtime/field_view.hpp"

#include <array>
#include <cstdint>
#include <memory>

namespace hundun::runtime {
class MpiContext;
}

namespace hundun::flow {

enum class MomentumTimeOrder : std::uint8_t { backward_euler = 1, bdf2 = 2 };

struct MomentumTimeStencil final {
  MomentumTimeOrder order;
  double dt_s;
  double previous_dt_s;
  double alpha0;
  double alpha1;
  double alpha2;
};

MomentumTimeStencil make_momentum_time_stencil(MomentumTimeOrder order,
                                               double dt_s,
                                               double previous_dt_s);

struct MomentumComponentEquation final {
  const linear::LinearOperator *linear_operator;
  linear::Preconditioner *preconditioner;
  execution::VectorView<const double> rhs;
  execution::VectorView<double> predictor;
  execution::VectorView<double> actual_diagonal;
};

struct MomentumPredictorReport final {
  std::array<linear::SolveReport, 3> components;
  bool all_converged() const noexcept;
};

class MomentumPredictor final {
public:
  explicit MomentumPredictor(const linear::LinearSolver &solver) noexcept;

  MomentumPredictorReport
  solve(const runtime::MpiContext &mpi,
        const std::array<MomentumComponentEquation, 3> &equations,
        const linear::SolveControl &control) const;

private:
  const linear::LinearSolver *solver_;
};

struct MomentumFaceHistory final {
  const runtime::FieldView<const double> &velocity_n;
  const runtime::FaceFieldView<const double> &face_velocity_n;
  const runtime::FieldView<const double> *velocity_n_minus_1;
  const runtime::FaceFieldView<const double> *face_velocity_n_minus_1;
};

struct MaterialMomentumFaceHistory final {
  const runtime::FieldView<const double> &density_n;
  const runtime::FieldView<const double> &velocity_n;
  const runtime::FaceFieldView<const double> &face_density_n;
  const runtime::FaceFieldView<const double> &face_velocity_n;
  const runtime::FieldView<const double> *density_n_minus_1;
  const runtime::FieldView<const double> *velocity_n_minus_1;
  const runtime::FaceFieldView<const double> *face_density_n_minus_1;
  const runtime::FaceFieldView<const double> *face_velocity_n_minus_1;
};

class TimeConsistentFaceVelocity final {
public:
  static TimeConsistentFaceVelocity create(const mesh::MeshTopology &topology,
                                           const mesh::MeshGeometry &geometry);

  ~TimeConsistentFaceVelocity() noexcept;
  TimeConsistentFaceVelocity(TimeConsistentFaceVelocity &&) noexcept;
  TimeConsistentFaceVelocity &operator=(TimeConsistentFaceVelocity &&) = delete;
  TimeConsistentFaceVelocity(const TimeConsistentFaceVelocity &) = delete;
  TimeConsistentFaceVelocity &
  operator=(const TimeConsistentFaceVelocity &) = delete;

  void assemble_constant_density(
      const boundary::BoundaryRegistry &boundaries, double rho_ref,
      const MomentumTimeStencil &stencil,
      const runtime::FieldView<const double> &predictor_velocity,
      const runtime::FieldView<const double> &mechanical_pressure,
      const runtime::FieldView<const double> &pressure_gradient,
      const runtime::FieldView<const double> &actual_diagonal,
      const MomentumFaceHistory &history,
      const runtime::FaceFieldView<double> &trial_face_velocity,
      const runtime::FieldRegistry &registry, runtime::FieldStorage &storage,
      const runtime::FieldAccessPlan &access_plan, runtime::PhaseId phase,
      runtime::ActorId actor, runtime::FieldId face_mass_flux_field) const;

  void assemble_material_density(
      const boundary::BoundaryRegistry &boundaries,
      const MomentumTimeStencil &stencil,
      const runtime::FieldView<const double> &predictor_velocity,
      const runtime::FieldView<const double> &mechanical_pressure,
      const runtime::FieldView<const double> &pressure_gradient,
      const runtime::FieldView<const double> &actual_diagonal,
      const MaterialMomentumFaceHistory &history,
      const runtime::FaceFieldView<double> &trial_face_velocity) const;

private:
  struct Impl;
  explicit TimeConsistentFaceVelocity(std::unique_ptr<Impl>) noexcept;
  std::unique_ptr<Impl> impl_;
};

} // namespace hundun::flow
