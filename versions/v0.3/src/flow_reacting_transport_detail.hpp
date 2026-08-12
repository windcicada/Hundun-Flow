// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/chem_composition.hpp"
#include "hundun/flow_material_density_transport.hpp"
#include "hundun/mesh_topology.hpp"

#include <optional>
#include <vector>

namespace hundun::flow::detail {

struct ReactingTransportCell final {
  double density_kg_per_m3{};
  std::vector<double> rho_y_kg_per_m3;
  double rho_h_tc_j_per_m3{};
  chemistry::ThermodynamicProperties thermodynamics;
  chemistry::TransportProperties molecular_transport;
  std::optional<double> wale_kinematic_diffusivity_m2_per_s;
};

struct ReactingTransportInteriorFace final {
  mesh::LocalFaceId local_face{};
  mesh::LocalCellId owner_cell{};
  mesh::LocalCellId neighbour_cell{};
  double area_m2{};
  double owner_to_neighbour_distance_m{};
  std::vector<double> muscl_owner_mass_fractions;
  std::vector<double> muscl_neighbour_mass_fractions;
  double muscl_owner_h_tc_j_per_kg{};
  double muscl_neighbour_h_tc_j_per_kg{};
};

struct ReactingTransportRequest final {
  runtime::FieldId shared_face_mass_flux_field{};
  chemistry::CompositionIdentity composition;
  std::vector<ReactingTransportCell> cells;
  std::vector<ReactingTransportInteriorFace> interior_faces;
};

struct ReactingTransportFaceFlux final {
  mesh::LocalFaceId local_face{};
  std::vector<double> species_advective_kg_per_s;
  std::vector<double> species_diffusive_kg_per_s;
  double enthalpy_advective_w{};
  double enthalpy_diffusive_w{};
};

struct ReactingTransportAssembly final {
  runtime::FieldId shared_face_mass_flux_field{};
  MaterialFluxProvenance flux_provenance{MaterialFluxProvenance::predictor};
  std::vector<ReactingTransportFaceFlux> face_fluxes;
  std::vector<std::vector<double>> species_residual_kg_per_s;
  std::vector<double> enthalpy_residual_w;
};

ReactingTransportAssembly
assemble_reacting_transport(const MaterialFaceMassFlux &,
                            const ReactingTransportRequest &);

} // namespace hundun::flow::detail
