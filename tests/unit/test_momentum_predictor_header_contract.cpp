// SPDX-License-Identifier: Apache-2.0

#include "hundun/flow/momentum_predictor.hpp"

#include <array>
#include <cstdint>
#include <memory>
#include <type_traits>

namespace {

using hundun::boundary::BoundaryRegistry;
using hundun::flow::MomentumComponentEquation;
using hundun::flow::MomentumFaceHistory;
using hundun::flow::MomentumPredictor;
using hundun::flow::MomentumPredictorReport;
using hundun::flow::MomentumTimeOrder;
using hundun::flow::MomentumTimeStencil;
using hundun::flow::TimeConsistentFaceVelocity;
using hundun::linear::LinearSolver;
using hundun::linear::SolveControl;
using hundun::mesh::MeshGeometry;
using hundun::mesh::MeshTopology;
using hundun::runtime::ActorId;
using hundun::runtime::FaceFieldView;
using hundun::runtime::FieldAccessPlan;
using hundun::runtime::FieldId;
using hundun::runtime::FieldRegistry;
using hundun::runtime::FieldStorage;
using hundun::runtime::FieldView;
using hundun::runtime::PhaseId;

using MakeStencil = MomentumTimeStencil (*)(MomentumTimeOrder, double, double);
using PredictorSolve = MomentumPredictorReport (MomentumPredictor::*)(
    const std::array<MomentumComponentEquation, 3> &,
    const SolveControl &) const;
using FaceCreate = TimeConsistentFaceVelocity (*)(const MeshTopology &,
                                                  const MeshGeometry &);
using FaceAssemble = void (TimeConsistentFaceVelocity::*)(
    const BoundaryRegistry &, double, const MomentumTimeStencil &,
    const FieldView<const double> &, const FieldView<const double> &,
    const FieldView<const double> &, const FieldView<const double> &,
    const MomentumFaceHistory &, const FaceFieldView<double> &,
    const FieldRegistry &, FieldStorage &, const FieldAccessPlan &, PhaseId,
    ActorId, FieldId) const;

static_assert(std::is_enum_v<MomentumTimeOrder>);
static_assert(
    std::is_same_v<std::underlying_type_t<MomentumTimeOrder>, std::uint8_t>);
static_assert(static_cast<int>(MomentumTimeOrder::backward_euler) == 1);
static_assert(static_cast<int>(MomentumTimeOrder::bdf2) == 2);
static_assert(std::is_final_v<MomentumTimeStencil>);
static_assert(std::is_final_v<MomentumComponentEquation>);
static_assert(std::is_final_v<MomentumPredictorReport>);
static_assert(std::is_final_v<MomentumFaceHistory>);
static_assert(std::is_final_v<MomentumPredictor>);
static_assert(std::is_final_v<TimeConsistentFaceVelocity>);
static_assert(
    std::is_same_v<decltype(&hundun::flow::make_momentum_time_stencil),
                   MakeStencil>);
static_assert(
    std::is_same_v<decltype(&MomentumPredictor::solve), PredictorSolve>);
static_assert(
    std::is_same_v<decltype(&TimeConsistentFaceVelocity::create), FaceCreate>);
static_assert(std::is_same_v<
              decltype(&TimeConsistentFaceVelocity::assemble_constant_density),
              FaceAssemble>);
static_assert(std::is_constructible_v<MomentumPredictor, const LinearSolver &>);
static_assert(!std::is_convertible_v<const LinearSolver &, MomentumPredictor>);
static_assert(!std::is_copy_constructible_v<TimeConsistentFaceVelocity>);
static_assert(!std::is_copy_assignable_v<TimeConsistentFaceVelocity>);
static_assert(std::is_nothrow_move_constructible_v<TimeConsistentFaceVelocity>);
static_assert(!std::is_move_assignable_v<TimeConsistentFaceVelocity>);
static_assert(std::is_nothrow_destructible_v<TimeConsistentFaceVelocity>);

} // namespace

int main() { return 0; }
