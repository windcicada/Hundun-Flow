// SPDX-License-Identifier: Apache-2.0

#include "hundun/diag_immersed_static.hpp"
#include "hundun/flow_checkpoint_v3.hpp"
#include "tests/support/test_main.hpp"

#include <array>
#include <cstdint>
#include <vector>

namespace hundun::application {

bool stage3_dispatch_inventory_complete(
    const flow::CheckpointV3WriteModules &) noexcept;

} // namespace hundun::application

namespace {

using namespace hundun;

template <class Type> const Type *object(std::uintptr_t identity) noexcept {
  return reinterpret_cast<const Type *>(identity);
}

flow::CheckpointV3WriteModules modules(std::uint8_t profile) {
  flow::CheckpointV3WriteModules result;
  result.presence = static_cast<flow::CheckpointV3Presence>(profile);
  const bool ibm = profile == 1U || profile == 3U || profile == 4U ||
                   profile == 6U || profile == 7U || profile == 9U;
  const bool wale = profile == 2U || profile == 3U || profile == 5U ||
                    profile == 6U || profile == 8U || profile == 9U;
  const bool ideal = profile >= 7U && profile <= 9U;
  if (ibm) {
    result.surface = object<immersed::ImmersedSurface>(1U);
    result.query = object<immersed::SurfaceQuery>(2U);
    result.domain = object<immersed::ImmersedDomain>(3U);
    result.ghost_plan = object<immersed::GhostStencilPlan>(4U);
    result.wall_plan = object<immersed::WallQuadraturePlan>(5U);
    result.transform = object<immersed::LocalFlowPatternTransform>(6U);
    result.flow = object<flow::FixedStepImmersedFlow>(7U);
  }
  if (wale)
    result.wale = object<les::WaleModel>(8U);
  if (ideal)
    result.ideal_gas = object<flow::IdealGasClosure>(9U);
  return result;
}

std::vector<diagnostics::DiagnosticModuleKind>
expected_providers(std::uint8_t profile) {
  using K = diagnostics::DiagnosticModuleKind;
  const bool ibm = profile == 1U || profile == 3U || profile == 4U ||
                   profile == 6U || profile == 7U || profile == 9U;
  const bool wale = profile == 2U || profile == 3U || profile == 5U ||
                    profile == 6U || profile == 8U || profile == 9U;
  std::vector<K> result;
  if (ibm)
    result = {K::immersed_surface, K::ghost_stencil,
              K::local_flow_pattern, K::wall_force};
  if (wale)
    result.push_back(K::les);
  return result;
}

void check_legal_matrix() {
  for (std::uint8_t profile = 1U; profile <= 9U; ++profile) {
    const auto row = modules(profile);
    HUNDUN_CHECK(application::stage3_dispatch_inventory_complete(row));
    HUNDUN_CHECK(diagnostics::stage3_added_provider_inventory(row.presence) ==
                 expected_providers(profile));
  }
}

void check_mutations() {
  auto missing_ibm = modules(6U);
  missing_ibm.surface = nullptr;
  HUNDUN_CHECK(
      !application::stage3_dispatch_inventory_complete(missing_ibm));

  auto unknown = modules(10U);
  HUNDUN_CHECK(!application::stage3_dispatch_inventory_complete(unknown));

  for (std::uint8_t profile = 1U; profile <= 9U; ++profile) {
    const auto legal = modules(profile);
    auto incomplete = legal;
    if (legal.flow != nullptr)
      incomplete.flow = nullptr;
    else
      incomplete.wale = nullptr;
    HUNDUN_CHECK(
        !application::stage3_dispatch_inventory_complete(incomplete));
  }
}

} // namespace

int main() {
  return hundun::test::run([] {
    check_legal_matrix();
    check_mutations();
  });
}
