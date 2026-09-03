// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "flow_reacting_transport_detail.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace hundun::flow::detail {

struct ReactingTransportFluxAccess final {
  static double value(const MaterialFaceMassFlux &flux,
                      mesh::LocalFaceId face) {
    return flux.value_for_reacting_transport(face);
  }
};

namespace {
bool positive(double value) noexcept {
  return std::isfinite(value) && value > 0.0;
}

void require_fractions(const std::vector<double> &values,
                       std::size_t species) {
  if (values.size() != species)
    throw std::invalid_argument("reacting transport species shape mismatch");
  double sum{};
  for (double value : values) {
    if (!std::isfinite(value) || value < 0.0)
      throw std::invalid_argument("reacting transport mass fraction is invalid");
    sum += value;
  }
  if (std::abs(sum - 1.0) > 64.0 * std::numeric_limits<double>::epsilon())
    throw std::invalid_argument("reacting transport fractions do not sum to one");
}
} // namespace

ReactingTransportAssembly
assemble_reacting_transport(const MaterialFaceMassFlux &mass_flux,
                            const ReactingTransportRequest &request) {
  chemistry::validate_composition_identity(request.composition);
  const std::size_t species = request.composition.species.size();
  if (mass_flux.provenance() != MaterialFluxProvenance::final_corrected ||
      mass_flux.field_id() != request.shared_face_mass_flux_field ||
      request.cells.empty())
    throw std::invalid_argument(
        "reacting transport requires accepted final face mass flux");

  for (const auto &cell : request.cells) {
    if (!positive(cell.density_kg_per_m3) ||
        cell.rho_y_kg_per_m3.size() != species ||
        !std::isfinite(cell.rho_h_tc_j_per_m3) ||
        !positive(cell.thermodynamics.cp_j_per_kg_k) ||
        !positive(cell.molecular_transport.conductivity_w_per_m_k) ||
        cell.molecular_transport.mixture_diffusivity_m2_per_s.size() != species ||
        (cell.wale_kinematic_diffusivity_m2_per_s &&
         (!std::isfinite(*cell.wale_kinematic_diffusivity_m2_per_s) ||
          *cell.wale_kinematic_diffusivity_m2_per_s < 0.0)))
      throw std::invalid_argument("reacting transport cell state is invalid");
    double density_sum{};
    for (std::size_t k = 0; k < species; ++k) {
      const double diffusivity =
          cell.molecular_transport.mixture_diffusivity_m2_per_s[k];
      if (!std::isfinite(cell.rho_y_kg_per_m3[k]) ||
          cell.rho_y_kg_per_m3[k] < 0.0 || !positive(diffusivity))
        throw std::invalid_argument("reacting transport cell state is invalid");
      density_sum += cell.rho_y_kg_per_m3[k];
    }
    if (std::abs(density_sum - cell.density_kg_per_m3) >
        64.0 * std::numeric_limits<double>::epsilon() *
            std::max(1.0, cell.density_kg_per_m3))
      throw std::invalid_argument("reacting transport density is stale");
  }

  ReactingTransportAssembly result;
  result.shared_face_mass_flux_field = mass_flux.field_id();
  result.flux_provenance = mass_flux.provenance();
  result.species_residual_kg_per_s.assign(
      species, std::vector<double>(request.cells.size(), 0.0));
  result.enthalpy_residual_w.assign(request.cells.size(), 0.0);
  result.face_fluxes.reserve(request.interior_faces.size());

  for (const auto &face : request.interior_faces) {
    if (face.owner_cell >= request.cells.size() ||
        face.neighbour_cell >= request.cells.size() ||
        face.owner_cell == face.neighbour_cell || !positive(face.area_m2) ||
        !positive(face.owner_to_neighbour_distance_m) ||
        !std::isfinite(face.muscl_owner_h_tc_j_per_kg) ||
        !std::isfinite(face.muscl_neighbour_h_tc_j_per_kg))
      throw std::invalid_argument("reacting transport face geometry is invalid");
    require_fractions(face.muscl_owner_mass_fractions, species);
    require_fractions(face.muscl_neighbour_mass_fractions, species);
    const auto &owner = request.cells[face.owner_cell];
    const auto &neighbour = request.cells[face.neighbour_cell];
    const double mdot =
        ReactingTransportFluxAccess::value(mass_flux, face.local_face);
    if (!std::isfinite(mdot))
      throw std::invalid_argument("reacting transport mass flux is non-finite");
    const double inverse_distance = 1.0 / face.owner_to_neighbour_distance_m;
    const double rho_face =
        0.5 * (owner.density_kg_per_m3 + neighbour.density_kg_per_m3);
    const double wale_face =
        0.5 * (owner.wale_kinematic_diffusivity_m2_per_s.value_or(0.0) +
               neighbour.wale_kinematic_diffusivity_m2_per_s.value_or(0.0));

    ReactingTransportFaceFlux assembled;
    assembled.local_face = face.local_face;
    assembled.species_advective_kg_per_s.resize(species);
    assembled.species_diffusive_kg_per_s.resize(species);
    std::vector<double> central_y(species);
    double raw_sum{};
    for (std::size_t k = 0; k < species; ++k) {
      const double owner_y =
          owner.rho_y_kg_per_m3[k] / owner.density_kg_per_m3;
      const double neighbour_y =
          neighbour.rho_y_kg_per_m3[k] / neighbour.density_kg_per_m3;
      central_y[k] = 0.5 * (owner_y + neighbour_y);
      const double diffusivity =
          0.5 * (owner.molecular_transport.mixture_diffusivity_m2_per_s[k] +
                 neighbour.molecular_transport.mixture_diffusivity_m2_per_s[k]) +
          wale_face;
      assembled.species_diffusive_kg_per_s[k] =
          -face.area_m2 * rho_face * diffusivity *
          (neighbour_y - owner_y) * inverse_distance;
      raw_sum += assembled.species_diffusive_kg_per_s[k];
      const auto &muscl = mdot >= 0.0 ? face.muscl_owner_mass_fractions
                                      : face.muscl_neighbour_mass_fractions;
      assembled.species_advective_kg_per_s[k] = mdot * muscl[k];
    }
    double corrected_prefix{};
    for (std::size_t k = 0; k + 1U < species; ++k) {
      assembled.species_diffusive_kg_per_s[k] -= central_y[k] * raw_sum;
      corrected_prefix += assembled.species_diffusive_kg_per_s[k];
    }
    // All raw species fluxes are evaluated above. This removes only the final
    // floating-point roundoff after applying correction velocity once.
    assembled.species_diffusive_kg_per_s.back() = -corrected_prefix;

    const double owner_h =
        owner.rho_h_tc_j_per_m3 / owner.density_kg_per_m3;
    const double neighbour_h =
        neighbour.rho_h_tc_j_per_m3 / neighbour.density_kg_per_m3;
    const double conductivity_over_cp =
        0.5 * (owner.molecular_transport.conductivity_w_per_m_k /
                   owner.thermodynamics.cp_j_per_kg_k +
               neighbour.molecular_transport.conductivity_w_per_m_k /
                   neighbour.thermodynamics.cp_j_per_kg_k);
    assembled.enthalpy_advective_w =
        mdot * (mdot >= 0.0 ? face.muscl_owner_h_tc_j_per_kg
                            : face.muscl_neighbour_h_tc_j_per_kg);
    assembled.enthalpy_diffusive_w =
        -face.area_m2 * (conductivity_over_cp + rho_face * wale_face) *
        (neighbour_h - owner_h) * inverse_distance;

    for (std::size_t k = 0; k < species; ++k) {
      const double total = assembled.species_advective_kg_per_s[k] +
                           assembled.species_diffusive_kg_per_s[k];
      result.species_residual_kg_per_s[k][face.owner_cell] += total;
      result.species_residual_kg_per_s[k][face.neighbour_cell] -= total;
    }
    const double total_h =
        assembled.enthalpy_advective_w + assembled.enthalpy_diffusive_w;
    result.enthalpy_residual_w[face.owner_cell] += total_h;
    result.enthalpy_residual_w[face.neighbour_cell] -= total_h;
    result.face_fluxes.push_back(std::move(assembled));
  }
  return result;
}

} // namespace hundun::flow::detail
