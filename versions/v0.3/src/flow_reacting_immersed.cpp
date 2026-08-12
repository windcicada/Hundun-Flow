// SPDX-License-Identifier: Apache-2.0

#include "flow_reacting_immersed_detail.hpp"

#include <cmath>
#include <stdexcept>

namespace hundun::flow::detail {

ReactingImmersedCoupling
compose_reacting_immersed_coupling(const ReactingImmersedRequest &request) {
  const auto &authority = request.authority;
  if (authority.geometry_fingerprint == 0U ||
      authority.wall_ownership_fingerprint == 0U ||
      !std::isfinite(authority.thermodynamic_pressure_pa) ||
      authority.thermodynamic_pressure_pa <= 0.0 ||
      authority.thermodynamic_pressure_pa !=
          request.expected_thermodynamic_pressure_pa ||
      authority.shared_face_mass_flux_field !=
          request.expected_face_mass_flux_field ||
      authority.flux_provenance != MaterialFluxProvenance::final_corrected ||
      request.molecular_diffusivity_m2_per_s.empty() ||
      request.wall.species_flux_kg_per_m2_s.size() !=
          request.molecular_diffusivity_m2_per_s.size() ||
      !std::isfinite(request.turbulent_schmidt) ||
      request.turbulent_schmidt <= 0.0)
    throw std::invalid_argument("reacting immersed authority is invalid");
  for (double flux : request.wall.species_flux_kg_per_m2_s)
    if (flux != 0.0)
      throw std::invalid_argument("reacting IBM wall must be non-catalytic");
  for (double value : request.molecular_diffusivity_m2_per_s)
    if (!std::isfinite(value) || value <= 0.0)
      throw std::invalid_argument("reacting molecular diffusivity is invalid");

  if (request.wale_kinematic_diffusivity_m2_per_s) {
    if (!authority.wale_identity || authority.wale_evaluation_count != 1U ||
        !std::isfinite(*request.wale_kinematic_diffusivity_m2_per_s) ||
        *request.wale_kinematic_diffusivity_m2_per_s < 0.0)
      throw std::invalid_argument("reacting WALE authority is invalid");
  } else if (authority.wale_identity || authority.wale_evaluation_count != 0U) {
    throw std::invalid_argument("reacting molecular path has stray WALE state");
  }

  ReactingImmersedCoupling result;
  result.geometry_fingerprint = authority.geometry_fingerprint;
  result.wall_ownership_fingerprint = authority.wall_ownership_fingerprint;
  result.species_wall_flux_kg_per_m2_s =
      request.wall.species_flux_kg_per_m2_s;
  result.wall_heat_flux_w_per_m2 = request.wall.heat_flux_w_per_m2;
  result.effective_diffusivity_m2_per_s =
      request.molecular_diffusivity_m2_per_s;
  if (request.wale_kinematic_diffusivity_m2_per_s) {
    const double turbulent =
        *request.wale_kinematic_diffusivity_m2_per_s /
        request.turbulent_schmidt;
    for (double &value : result.effective_diffusivity_m2_per_s)
      value += turbulent;
  }
  result.wall_source = {"immersed-reacting-wall",
                        ReactingSourceKind::boundary};
  result.wale_identity = authority.wale_identity;
  return result;
}

} // namespace hundun::flow::detail
