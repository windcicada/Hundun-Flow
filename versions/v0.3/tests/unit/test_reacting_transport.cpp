// SPDX-License-Identifier: Apache-2.0

#include "flow_reacting_transport_detail.hpp"

#include "hundun/fvm_cell_centered.hpp"
#include "hundun/flow_material_density_transport.hpp"
#include "hundun/mesh_topology.hpp"
#include "hundun/rt_field_access_plan.hpp"
#include "hundun/rt_field_registry.hpp"
#include "hundun/rt_field_storage.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "hundun/rt_structured_decomposition.hpp"
#include "tests/support/test_main.hpp"

#include <cmath>
#include <stdexcept>

namespace {

using namespace hundun;

chemistry::CompositionIdentity composition() {
  chemistry::CompositionIdentity value;
  value.element_names = {"A", "B"};
  value.species = {{"A2", 2.0, {2, 0}},
                   {"B2", 4.0, {0, 2}},
                   {"AB", 3.0, {1, 1}}};
  value.fingerprint = chemistry::composition_identity_fingerprint(value);
  return value;
}

void run(const runtime::MpiContext &mpi) {
  auto decomposition = runtime::StructuredDecomposition::create(
      mpi, {2, 1, 1}, {false, false, false},
      runtime::DecompositionOptions{runtime::Int3{1, 1, 1}});
  mesh::MeshTopology topology(decomposition);
  runtime::FieldRegistry registry;
  const auto flux_field = finite_volume::declare_face_mass_flux(registry);
  registry.freeze();
  runtime::FieldStorage storage(
      registry, {decomposition.local_extent(), topology.local_face_count()});
  constexpr runtime::PhaseId phase = 4202U;
  constexpr runtime::ActorId actor = 4202U;
  runtime::FieldAccessPlan access(registry);
  access.declare_access(phase, actor, flux_field,
                        runtime::AccessMode::read_write);
  access.freeze();
  auto values = storage.acquire_face_write<double>(access, phase, actor,
                                                    flux_field);
  for (std::size_t face = 0; face < topology.local_face_count(); ++face)
    values(face, 0) = 0.0;

  mesh::LocalFaceId interior{};
  bool found = false;
  for (std::size_t face = 0; face < topology.local_face_count(); ++face) {
    if (topology.neighbour(face)) {
      interior = face;
      found = true;
      break;
    }
  }
  HUNDUN_CHECK(found);
  values(interior, 0) = 2.0;
  auto final_flux = flow::MaterialFaceMassFlux::acquire(
      registry, storage, access, phase, actor, flux_field, topology,
      flow::MaterialFluxProvenance::final_corrected);

  flow::detail::ReactingTransportRequest request;
  request.shared_face_mass_flux_field = flux_field;
  request.composition = composition();
  request.cells.resize(2U);
  request.cells[0].density_kg_per_m3 = 1.0;
  request.cells[0].rho_y_kg_per_m3 = {0.2, 0.3, 0.5};
  request.cells[0].rho_h_tc_j_per_m3 = 10.0;
  request.cells[0].thermodynamics = {300.0, 1.0, 2.0, 3.0};
  request.cells[0].molecular_transport = {1.0e-5, 2.0, {1.0e-5, 2.0e-5, 3.0e-5}};
  request.cells[1].density_kg_per_m3 = 2.0;
  request.cells[1].rho_y_kg_per_m3 = {0.8, 0.2, 1.0};
  request.cells[1].rho_h_tc_j_per_m3 = 28.0;
  request.cells[1].thermodynamics = {400.0, 2.0, 2.0, 3.0};
  request.cells[1].molecular_transport = {1.0e-5, 2.0, {1.0e-5, 2.0e-5, 3.0e-5}};

  flow::detail::ReactingTransportInteriorFace face;
  face.local_face = interior;
  face.owner_cell = topology.owner(interior);
  face.neighbour_cell = *topology.neighbour(interior);
  face.area_m2 = 1.0;
  face.owner_to_neighbour_distance_m = 0.5;
  face.muscl_owner_mass_fractions = {0.25, 0.25, 0.5};
  face.muscl_neighbour_mass_fractions = {0.35, 0.15, 0.5};
  face.muscl_owner_h_tc_j_per_kg = 11.0;
  face.muscl_neighbour_h_tc_j_per_kg = 13.0;
  request.interior_faces.push_back(face);

  const auto report =
      flow::detail::assemble_reacting_transport(final_flux, request);
  HUNDUN_CHECK(report.shared_face_mass_flux_field == flux_field);
  HUNDUN_CHECK(report.flux_provenance ==
               flow::MaterialFluxProvenance::final_corrected);
  HUNDUN_CHECK(report.face_fluxes.size() == 1U);
  const auto &assembled = report.face_fluxes.front();
  HUNDUN_CHECK(assembled.species_advective_kg_per_s.size() == 3U);
  HUNDUN_CHECK_NEAR(assembled.species_advective_kg_per_s[0], 0.5, 1.0e-15);
  HUNDUN_CHECK_NEAR(assembled.species_advective_kg_per_s[1], 0.5, 1.0e-15);
  HUNDUN_CHECK_NEAR(assembled.species_advective_kg_per_s[2], 1.0, 1.0e-15);
  double diffusion_sum = 0.0;
  for (double value : assembled.species_diffusive_kg_per_s)
    diffusion_sum += value;
  HUNDUN_CHECK(diffusion_sum == 0.0);
  HUNDUN_CHECK_NEAR(assembled.enthalpy_advective_w, 22.0, 1.0e-14);
  HUNDUN_CHECK_NEAR(assembled.enthalpy_diffusive_w, -8.0, 1.0e-14);

  for (std::size_t species = 0; species < 3U; ++species) {
    HUNDUN_CHECK_NEAR(report.species_residual_kg_per_s[species][0] +
                          report.species_residual_kg_per_s[species][1],
                      0.0, 0.0);
  }
  HUNDUN_CHECK_NEAR(report.enthalpy_residual_w[0] +
                        report.enthalpy_residual_w[1],
                    0.0, 0.0);
  for (std::size_t element = 0; element < 2U; ++element) {
    double element_residual = 0.0;
    for (std::size_t species = 0; species < 3U; ++species) {
      element_residual +=
          (report.species_residual_kg_per_s[species][0] +
           report.species_residual_kg_per_s[species][1]) *
          static_cast<double>(request.composition.species[species]
                                  .element_counts[element]) /
          request.composition.species[species]
              .molecular_weight_kg_per_kmol;
    }
    HUNDUN_CHECK(element_residual == 0.0);
  }

  auto predictor = flow::MaterialFaceMassFlux::acquire(
      registry, storage, access, phase, actor, flux_field, topology,
      flow::MaterialFluxProvenance::predictor);
  bool rejected = false;
  try {
    static_cast<void>(
        flow::detail::assemble_reacting_transport(predictor, request));
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);

  auto stale = request;
  stale.cells[0].density_kg_per_m3 = 1.1;
  rejected = false;
  try {
    static_cast<void>(
        flow::detail::assemble_reacting_transport(final_flux, stale));
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);

  auto clipped = request;
  clipped.interior_faces[0].muscl_owner_mass_fractions = {-0.1, 0.6, 0.5};
  rejected = false;
  try {
    static_cast<void>(
        flow::detail::assemble_reacting_transport(final_flux, clipped));
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);

  auto wale = request;
  wale.cells[0].wale_kinematic_diffusivity_m2_per_s = 1.0e-5;
  wale.cells[1].wale_kinematic_diffusivity_m2_per_s = 1.0e-5;
  const auto wale_report =
      flow::detail::assemble_reacting_transport(final_flux, wale);
  double wale_diffusion_sum = 0.0;
  for (double value :
       wale_report.face_fluxes.front().species_diffusive_kg_per_s)
    wale_diffusion_sum += value;
  HUNDUN_CHECK(wale_diffusion_sum == 0.0);
  HUNDUN_CHECK(wale_report.face_fluxes.front()
                   .species_diffusive_kg_per_s[0] !=
               assembled.species_diffusive_kg_per_s[0]);
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  return hundun::test::run([&] {
    auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
    run(mpi);
  });
}
