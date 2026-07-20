// SPDX-License-Identifier: Apache-2.0

#include "hundun/boundary/basic_boundary.hpp"
#include "hundun/finite_volume/poisson_boundary_adapter.hpp"

#include "tests/support/test_main.hpp"

#include <cstdint>
#include <type_traits>
#include <utility>

namespace {

using hundun::boundary::BoundaryDescriptor;
using hundun::boundary::BoundaryKind;
using hundun::boundary::BoundaryRegistry;
using hundun::boundary::FinalFluxAdmissibility;
using hundun::boundary::FinalFluxDecision;
using hundun::boundary::MassFluxRule;
using hundun::boundary::PressureRule;
using hundun::boundary::ScalarBoundaryValues;
using hundun::boundary::TransportRule;
using hundun::boundary::VelocityBoundaryValues;
using hundun::boundary::VelocityRule;

static_assert(
    std::is_same_v<std::underlying_type_t<BoundaryKind>, std::uint8_t>);
static_assert(
    std::is_same_v<std::underlying_type_t<VelocityRule>, std::uint8_t>);
static_assert(
    std::is_same_v<std::underlying_type_t<PressureRule>, std::uint8_t>);
static_assert(
    std::is_same_v<std::underlying_type_t<TransportRule>, std::uint8_t>);
static_assert(
    std::is_same_v<std::underlying_type_t<MassFluxRule>, std::uint8_t>);
static_assert(
    std::is_same_v<std::underlying_type_t<FinalFluxDecision>, std::uint8_t>);

static_assert(std::is_final_v<BoundaryDescriptor>);
static_assert(!std::is_copy_constructible_v<BoundaryDescriptor>);
static_assert(!std::is_copy_assignable_v<BoundaryDescriptor>);
static_assert(std::is_nothrow_move_constructible_v<BoundaryDescriptor>);
static_assert(std::is_nothrow_move_assignable_v<BoundaryDescriptor>);

static_assert(std::is_final_v<BoundaryRegistry>);
static_assert(!std::is_copy_constructible_v<BoundaryRegistry>);
static_assert(!std::is_copy_assignable_v<BoundaryRegistry>);
static_assert(std::is_nothrow_move_constructible_v<BoundaryRegistry>);
static_assert(!std::is_move_assignable_v<BoundaryRegistry>);
static_assert(std::is_nothrow_destructible_v<BoundaryRegistry>);

using CreateResult = decltype(BoundaryRegistry::create(
    std::declval<const hundun::config::FlowCaseConfig &>(),
    std::declval<const hundun::mesh::MeshTopology &>()));
static_assert(std::is_same_v<CreateResult, BoundaryRegistry>);

using VelocityResult =
    decltype(std::declval<const BoundaryRegistry &>().evaluate_velocity(
        std::uint32_t{}, hundun::runtime::Real3{}, hundun::runtime::Real3{}));
using PressureResult =
    decltype(std::declval<const BoundaryRegistry &>().evaluate_pressure(
        std::uint32_t{}, double{}));
using FluxResult =
    decltype(std::declval<const BoundaryRegistry &>()
                 .assess_final_pressure_outlet_flux(
                     std::declval<const hundun::mesh::MeshTopology &>(),
                     std::declval<const hundun::runtime::MpiContext &>(),
                     std::declval<const hundun::runtime::FaceFieldView<
                         const double> &>(),
                     std::uint64_t{}, double{}));
static_assert(std::is_same_v<VelocityResult, VelocityBoundaryValues>);
static_assert(std::is_same_v<PressureResult, ScalarBoundaryValues>);
static_assert(std::is_same_v<FluxResult, FinalFluxAdmissibility>);

using AdapterResult =
    decltype(hundun::finite_volume::make_poisson_boundary_spec(
        std::declval<const BoundaryRegistry &>()));
static_assert(
    std::is_same_v<AdapterResult, hundun::finite_volume::PoissonBoundarySpec>);

} // namespace

int main() {
  return hundun::test::run([] {
    HUNDUN_CHECK(static_cast<std::uint8_t>(BoundaryKind::periodic) == 0U);
    HUNDUN_CHECK(static_cast<std::uint8_t>(BoundaryKind::pressure_outlet) ==
                 4U);
    HUNDUN_CHECK(static_cast<std::uint8_t>(FinalFluxDecision::admissible) ==
                 0U);
  });
}
