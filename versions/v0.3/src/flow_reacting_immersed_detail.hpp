// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "flow_reacting_boundary_detail.hpp"
#include "flow_reacting_transaction_detail.hpp"
#include "hundun/flow_material_density_transport.hpp"
#include "hundun/les_wale.hpp"

#include <optional>
#include <vector>

namespace hundun::flow::detail {

struct ReactingImmersedAuthority final {
  std::uint64_t geometry_fingerprint{};
  std::uint64_t wall_ownership_fingerprint{};
  double thermodynamic_pressure_pa{};
  runtime::FieldId shared_face_mass_flux_field{};
  MaterialFluxProvenance flux_provenance{MaterialFluxProvenance::predictor};
  std::optional<les::WaleCoefficientIdentity> wale_identity;
  std::uint32_t wale_evaluation_count{};
};

struct ReactingImmersedRequest final {
  ReactingImmersedAuthority authority;
  runtime::FieldId expected_face_mass_flux_field{};
  double expected_thermodynamic_pressure_pa{};
  ReactingWallBoundaryState wall;
  std::vector<double> molecular_diffusivity_m2_per_s;
  std::optional<double> wale_kinematic_diffusivity_m2_per_s;
  double turbulent_schmidt{1.0};
};

struct ReactingImmersedCoupling final {
  std::uint64_t geometry_fingerprint{};
  std::uint64_t wall_ownership_fingerprint{};
  std::vector<double> species_wall_flux_kg_per_m2_s;
  double wall_heat_flux_w_per_m2{};
  std::vector<double> effective_diffusivity_m2_per_s;
  ReactingSourceIdentity wall_source;
  std::optional<les::WaleCoefficientIdentity> wale_identity;
};

ReactingImmersedCoupling
compose_reacting_immersed_coupling(const ReactingImmersedRequest &);

} // namespace hundun::flow::detail
