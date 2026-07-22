// SPDX-License-Identifier: Apache-2.0

#include "hundun/diagnostics/material_density_piso_diagnostics.hpp"
#include "hundun/flow/material_density_piso.hpp"

#include <array>
#include <cstdint>
#include <string_view>
#include <type_traits>
#include <utility>

namespace {

using History = hundun::flow::MaterialMomentumFaceHistory;
using FaceAssembler = hundun::flow::TimeConsistentFaceVelocity;
using Coupler = hundun::flow::PisoCoupler;
using Flow = hundun::flow::FixedStepMaterialDensityFlow;
using Report = hundun::flow::MaterialDensityStepAttemptReport;
using Source = hundun::flow::MaterialDensityFlowDiagnosticSource;

using AssembleMaterial = void (FaceAssembler::*)(
    const hundun::boundary::BoundaryRegistry &,
    const hundun::flow::MomentumTimeStencil &,
    const hundun::runtime::FieldView<const double> &,
    const hundun::runtime::FieldView<const double> &,
    const hundun::runtime::FieldView<const double> &,
    const hundun::runtime::FieldView<const double> &, const History &,
    const hundun::runtime::FaceFieldView<double> &) const;

using CorrectMaterial = hundun::flow::PressureCorrectionReport (Coupler::*)(
    hundun::flow::FlowState &,
    const hundun::flow::MomentumTimeStencil &,
    const hundun::runtime::FieldView<const double> &,
    const hundun::linear::SolveControl &) const;

static_assert(std::is_final_v<History>);
static_assert(std::is_aggregate_v<History>);
static_assert(std::is_same_v<decltype(History::density_n),
                             const hundun::runtime::FieldView<const double> &>);
static_assert(std::is_same_v<decltype(History::velocity_n),
                             const hundun::runtime::FieldView<const double> &>);
static_assert(std::is_same_v<
              decltype(History::face_density_n),
              const hundun::runtime::FaceFieldView<const double> &>);
static_assert(std::is_same_v<
              decltype(History::face_velocity_n),
              const hundun::runtime::FaceFieldView<const double> &>);
static_assert(std::is_same_v<
              decltype(History::density_n_minus_1),
              const hundun::runtime::FieldView<const double> *>);
static_assert(std::is_same_v<
              decltype(History::velocity_n_minus_1),
              const hundun::runtime::FieldView<const double> *>);
static_assert(std::is_same_v<
              decltype(History::face_density_n_minus_1),
              const hundun::runtime::FaceFieldView<const double> *>);
static_assert(std::is_same_v<
              decltype(History::face_velocity_n_minus_1),
              const hundun::runtime::FaceFieldView<const double> *>);
static_assert(std::is_same_v<decltype(&FaceAssembler::assemble_material_density),
                             AssembleMaterial>);
static_assert(std::is_same_v<decltype(&Coupler::correct_material_density),
                             CorrectMaterial>);

static_assert(std::is_nothrow_destructible_v<Flow>);
static_assert(std::is_nothrow_move_constructible_v<Flow>);
static_assert(!std::is_copy_constructible_v<Flow>);
static_assert(!std::is_copy_assignable_v<Flow>);
static_assert(!std::is_move_assignable_v<Flow>);
static_assert(std::is_copy_constructible_v<Report>);
static_assert(std::is_nothrow_move_constructible_v<Report>);
static_assert(std::is_nothrow_destructible_v<Source>);
static_assert(std::is_nothrow_move_constructible_v<Source>);
static_assert(!std::is_copy_constructible_v<Source>);
static_assert(!std::is_copy_assignable_v<Source>);
static_assert(!std::is_move_assignable_v<Source>);

static_assert(std::is_same_v<
              decltype(std::declval<const Report &>().material_field_count()),
              std::uint64_t>);
static_assert(std::is_same_v<
              decltype(std::declval<const Report &>()
                           .final_momentum_residual_availability()),
              const std::array<std::uint8_t, 3> &>);
static_assert(std::is_same_v<
              decltype(std::declval<const Source &>().fingerprint_field_id(0)),
              std::string_view>);
static_assert(std::is_same_v<
              decltype(std::declval<const Source &>().field_entity(0)),
              hundun::flow::MaterialDensityDiagnosticEntity>);

} // namespace

int main() { return 0; }
