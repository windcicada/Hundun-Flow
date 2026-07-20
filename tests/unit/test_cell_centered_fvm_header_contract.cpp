// SPDX-License-Identifier: Apache-2.0

#include "hundun/finite_volume/cell_centered_fvm.hpp"

#include <cstddef>
#include <cstdint>
#include <type_traits>

namespace {

using hundun::boundary::BoundaryRegistry;
using hundun::finite_volume::CellCenteredFvmOperators;
using hundun::finite_volume::FaceMassFlux;
using hundun::finite_volume::FiniteVolumeQuantity;
using hundun::finite_volume::FiniteVolumeQuantityKind;
using hundun::finite_volume::GradientScheme;
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

static_assert(std::is_enum_v<GradientScheme>);
static_assert(
    std::is_same_v<std::underlying_type_t<GradientScheme>, std::uint8_t>);
static_assert(std::is_enum_v<FiniteVolumeQuantityKind>);
static_assert(std::is_same_v<std::underlying_type_t<FiniteVolumeQuantityKind>,
                             std::uint8_t>);
static_assert(std::is_final_v<FiniteVolumeQuantity>);
static_assert(std::is_final_v<FaceMassFlux>);
static_assert(std::is_final_v<CellCenteredFvmOperators>);
static_assert(!std::is_copy_constructible_v<FaceMassFlux>);
static_assert(!std::is_copy_assignable_v<FaceMassFlux>);
static_assert(std::is_nothrow_move_constructible_v<FaceMassFlux>);
static_assert(!std::is_move_assignable_v<FaceMassFlux>);
static_assert(std::is_nothrow_destructible_v<FaceMassFlux>);
static_assert(!std::is_copy_constructible_v<CellCenteredFvmOperators>);
static_assert(!std::is_copy_assignable_v<CellCenteredFvmOperators>);
static_assert(std::is_nothrow_move_constructible_v<CellCenteredFvmOperators>);
static_assert(!std::is_move_assignable_v<CellCenteredFvmOperators>);

using FluxAcquire = FaceMassFlux (*)(const FieldRegistry &,
                                     const FieldStorage &,
                                     const FieldAccessPlan &, PhaseId, ActorId,
                                     FieldId, const MeshTopology &);
using OperatorCreate = CellCenteredFvmOperators (*)(const MeshTopology &,
                                                    const MeshGeometry &);
using ComputeGradient = void (CellCenteredFvmOperators::*)(
    GradientScheme, FiniteVolumeQuantity, const BoundaryRegistry &,
    const FieldView<const double> &, const FieldView<double> &) const;
using ReconstructTransport = void (CellCenteredFvmOperators::*)(
    FiniteVolumeQuantity, const BoundaryRegistry &, const FaceMassFlux &,
    const FieldView<const double> &, const FaceFieldView<double> &) const;
using ReconstructMomentum = void (CellCenteredFvmOperators::*)(
    const BoundaryRegistry &, const FaceMassFlux &,
    const FieldView<const double> &, const FaceFieldView<double> &) const;
using AssembleProvisionalMassFlux = void (CellCenteredFvmOperators::*)(
    const BoundaryRegistry &, const FieldView<const double> &,
    const FieldView<const double> &, const FieldRegistry &, FieldStorage &,
    const FieldAccessPlan &, PhaseId, ActorId, FieldId) const;
using AccumulateMass = void (CellCenteredFvmOperators::*)(
    const FaceMassFlux &, const FieldView<double> &) const;
using AccumulateConvective = void (CellCenteredFvmOperators::*)(
    const FaceMassFlux &, const FaceFieldView<const double> &,
    const FieldView<double> &) const;
using AccumulateDiffusive = void (CellCenteredFvmOperators::*)(
    FiniteVolumeQuantity, const BoundaryRegistry &,
    const FieldView<const double> &, const FieldView<const double> &,
    const FaceFieldView<const double> &, const FieldView<double> &) const;
using AccumulateViscous = void (CellCenteredFvmOperators::*)(
    const BoundaryRegistry &, const FieldView<const double> &,
    const FieldView<const double> &, double, const FieldView<double> &) const;
using FluxDescriptor = hundun::runtime::FieldDescriptor (*)();
using DeclareFlux = FieldId (*)(FieldRegistry &);
using RequireFlux = void (*)(const FieldRegistry &, FieldId);
using Limiter = double (*)(double, double) noexcept;

static_assert(std::is_same_v<decltype(&FaceMassFlux::acquire), FluxAcquire>);
static_assert(std::is_same_v<decltype(&CellCenteredFvmOperators::create),
                             OperatorCreate>);
static_assert(
    std::is_same_v<decltype(&CellCenteredFvmOperators::compute_gradient),
                   ComputeGradient>);
static_assert(std::is_same_v<
              decltype(&CellCenteredFvmOperators::reconstruct_transport_faces),
              ReconstructTransport>);
static_assert(std::is_same_v<
              decltype(&CellCenteredFvmOperators::reconstruct_momentum_faces),
              ReconstructMomentum>);
static_assert(
    std::is_same_v<
        decltype(&CellCenteredFvmOperators::assemble_provisional_mass_flux),
        AssembleProvisionalMassFlux>);
static_assert(std::is_same_v<
              decltype(&CellCenteredFvmOperators::accumulate_mass_residual),
              AccumulateMass>);
static_assert(
    std::is_same_v<
        decltype(&CellCenteredFvmOperators::accumulate_convective_residual),
        AccumulateConvective>);
static_assert(std::is_same_v<decltype(&CellCenteredFvmOperators::
                                          accumulate_scalar_diffusive_residual),
                             AccumulateDiffusive>);
static_assert(std::is_same_v<
              decltype(&CellCenteredFvmOperators::accumulate_viscous_residual),
              AccumulateViscous>);
static_assert(
    std::is_same_v<decltype(&hundun::finite_volume::face_mass_flux_descriptor),
                   FluxDescriptor>);
static_assert(
    std::is_same_v<decltype(&hundun::finite_volume::declare_face_mass_flux),
                   DeclareFlux>);
static_assert(std::is_same_v<
              decltype(&hundun::finite_volume::require_face_mass_flux_field),
              RequireFlux>);
static_assert(std::is_same_v<
              decltype(&hundun::finite_volume::monotonized_central), Limiter>);

} // namespace

int main() {
  const auto density = FiniteVolumeQuantity::density();
  const auto scalar = FiniteVolumeQuantity::scalar(3U);
  return density.kind == FiniteVolumeQuantityKind::density &&
                 density.scalar_index == 0U &&
                 scalar.kind == FiniteVolumeQuantityKind::scalar &&
                 scalar.scalar_index == 3U
             ? 0
             : 1;
}
