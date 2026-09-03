// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/fvm_matrix_free_poisson.hpp"
#include "hundun/mesh_geometry.hpp"

#include <cstdint>
#include <optional>
#include <type_traits>

namespace {

using hundun::execution::ExecutionContext;
using hundun::execution::ExecutionEvent;
using hundun::execution::VectorView;
using hundun::finite_volume::MatrixFreePoissonOperator;
using hundun::finite_volume::PoissonBoundarySpec;
using hundun::finite_volume::PoissonCoefficientReplacementResult;
using hundun::finite_volume::PoissonConstructionError;
using hundun::finite_volume::PoissonConstraint;
using hundun::finite_volume::PoissonSolverFamily;
using hundun::finite_volume::PressureConstraintMode;
using hundun::linear::LinearOperator;
using hundun::linear::VectorLayout;
using hundun::mesh::LocalFaceId;
using hundun::mesh::MeshGeometry;
using hundun::mesh::MeshTopology;
using hundun::runtime::MpiContext;
using hundun::runtime::Real3;
using hundun::runtime::StructuredDecomposition;

using OperatorCreate = MatrixFreePoissonOperator (*)(
    const StructuredDecomposition&, const MeshTopology&, const MeshGeometry&,
    ExecutionContext&, VectorView<const double>, PoissonBoundarySpec);
using OperatorApply = ExecutionEvent (MatrixFreePoissonOperator::*)(
    VectorView<const double>, VectorView<double>) const;
using OperatorDiagonal = ExecutionEvent (MatrixFreePoissonOperator::*)(
    VectorView<double>) const;
using ReplaceFaceCoefficients = std::uint64_t (
    MatrixFreePoissonOperator::*)(VectorView<const double>);
using CollectivelyReplaceFaceCoefficients =
    PoissonCoefficientReplacementResult (MatrixFreePoissonOperator::*)(
        VectorView<const double>, const MpiContext&);
using AccumulateExplicitCorrection = ExecutionEvent (
    MatrixFreePoissonOperator::*)(VectorView<const double>,
                                  VectorView<const double>,
                                  VectorView<const double>,
                                  VectorView<double>) const;
using ConstraintCreate = PoissonConstraint (*)(
    const MeshTopology&, const MeshGeometry&, ExecutionContext&,
    const MpiContext&, PressureConstraintMode);
using ConstraintOperation = void (PoissonConstraint::*)(VectorView<double>)
    const;
using FaceDisplacement = Real3 (MeshGeometry::*)(LocalFaceId) const;

static_assert(std::is_enum_v<PressureConstraintMode>);
static_assert(std::is_same_v<std::underlying_type_t<PressureConstraintMode>,
                             std::uint8_t>);
static_assert(std::is_enum_v<PoissonSolverFamily>);
static_assert(std::is_same_v<std::underlying_type_t<PoissonSolverFamily>,
                             std::uint8_t>);
static_assert(std::is_final_v<PoissonBoundarySpec>);
static_assert(std::is_final_v<PoissonCoefficientReplacementResult>);
static_assert(std::is_final_v<PoissonConstructionError>);
static_assert(
    std::is_base_of_v<hundun::runtime::Error, PoissonConstructionError>);
static_assert(std::is_final_v<MatrixFreePoissonOperator>);
static_assert(std::is_same_v<decltype(&PoissonConstructionError::failing_rank),
                             int (PoissonConstructionError::*)()
                                 const noexcept>);
static_assert(std::is_final_v<PoissonConstraint>);
static_assert(std::is_base_of_v<LinearOperator, MatrixFreePoissonOperator>);

static_assert(!std::is_copy_constructible_v<MatrixFreePoissonOperator>);
static_assert(!std::is_copy_assignable_v<MatrixFreePoissonOperator>);
static_assert(std::is_nothrow_move_constructible_v<
              MatrixFreePoissonOperator>);
static_assert(!std::is_move_assignable_v<MatrixFreePoissonOperator>);
static_assert(std::is_nothrow_destructible_v<MatrixFreePoissonOperator>);
static_assert(!std::is_copy_constructible_v<PoissonConstraint>);
static_assert(!std::is_copy_assignable_v<PoissonConstraint>);
static_assert(std::is_nothrow_move_constructible_v<PoissonConstraint>);
static_assert(!std::is_move_assignable_v<PoissonConstraint>);
static_assert(std::is_nothrow_destructible_v<PoissonConstraint>);

static_assert(std::is_same_v<decltype(&MatrixFreePoissonOperator::create),
                             OperatorCreate>);
static_assert(std::is_same_v<decltype(&MatrixFreePoissonOperator::apply),
                             OperatorApply>);
static_assert(std::is_same_v<decltype(&MatrixFreePoissonOperator::diagonal),
                             OperatorDiagonal>);
static_assert(std::is_same_v<decltype(
                                 &MatrixFreePoissonOperator::
                                     replace_face_coefficients),
                             ReplaceFaceCoefficients>);
static_assert(std::is_same_v<
              decltype(&MatrixFreePoissonOperator::
                           collectively_replace_face_coefficients),
              CollectivelyReplaceFaceCoefficients>);
static_assert(std::is_same_v<
              decltype(PoissonCoefficientReplacementResult::accepted),
              bool>);
static_assert(std::is_same_v<
              decltype(PoissonCoefficientReplacementResult::changed),
              bool>);
static_assert(std::is_same_v<decltype(
                                 PoissonCoefficientReplacementResult::
                                     lowest_failing_rank),
                             int>);
static_assert(std::is_same_v<
              decltype(PoissonCoefficientReplacementResult::revision),
              std::uint64_t>);
static_assert(std::is_same_v<decltype(
                                 &MatrixFreePoissonOperator::
                                     accumulate_explicit_nonorthogonal_rhs),
                             AccumulateExplicitCorrection>);
static_assert(std::is_same_v<decltype(
                                 &MatrixFreePoissonOperator::domain_layout),
                             VectorLayout (MatrixFreePoissonOperator::*)()
                                 const>);
static_assert(std::is_same_v<decltype(
                                 &MatrixFreePoissonOperator::range_layout),
                             VectorLayout (MatrixFreePoissonOperator::*)()
                                 const>);
static_assert(std::is_same_v<decltype(&MatrixFreePoissonOperator::context),
                             const ExecutionContext& (
                                 MatrixFreePoissonOperator::*)() const>);
static_assert(std::is_same_v<decltype(&MatrixFreePoissonOperator::revision),
                             std::uint64_t (
                                 MatrixFreePoissonOperator::*)() const>);
static_assert(std::is_same_v<decltype(
                                 &MatrixFreePoissonOperator::has_diagonal),
                             bool (MatrixFreePoissonOperator::*)() const>);
static_assert(std::is_same_v<decltype(
                                 &MatrixFreePoissonOperator::constraint_mode),
                             PressureConstraintMode (
                                 MatrixFreePoissonOperator::*)()
                                 const noexcept>);
static_assert(std::is_same_v<decltype(
                                 &MatrixFreePoissonOperator::
                                     pressure_reference_patch_id),
                             std::optional<std::uint32_t> (
                                 MatrixFreePoissonOperator::*)()
                                 const noexcept>);
static_assert(std::is_same_v<decltype(
                                 &MatrixFreePoissonOperator::solver_family),
                             PoissonSolverFamily (
                                 MatrixFreePoissonOperator::*)()
                                 const noexcept>);

static_assert(std::is_same_v<decltype(&PoissonConstraint::create),
                             ConstraintCreate>);
static_assert(std::is_same_v<decltype(&PoissonConstraint::mode),
                             PressureConstraintMode (PoissonConstraint::*)()
                                 const noexcept>);
static_assert(std::is_same_v<decltype(&PoissonConstraint::project_rhs),
                             ConstraintOperation>);
static_assert(std::is_same_v<decltype(
                                 &PoissonConstraint::normalize_solution),
                             ConstraintOperation>);
static_assert(std::is_same_v<decltype(&MeshGeometry::face_displacement_m),
                             FaceDisplacement>);

}  // namespace

int main() {
  const PoissonBoundarySpec boundary;
  return boundary.mode == PressureConstraintMode::constant_nullspace &&
                 !boundary.pressure_reference_patch_id.has_value()
             ? 0
             : 1;
}
