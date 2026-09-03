// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "flow_reacting_transport_detail.hpp"

#include "hundun/fvm_cell_centered.hpp"
#include "hundun/mesh_topology.hpp"
#include "hundun/rt_field_access_plan.hpp"
#include "hundun/rt_field_registry.hpp"
#include "hundun/rt_field_storage.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "hundun/rt_structured_decomposition.hpp"
#include "tests/support/test_main.hpp"

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
  const runtime::Int3 process_grid{mpi.size(), 1, 1};
  auto decomposition = runtime::StructuredDecomposition::create(
      mpi, {2, 1, 1}, {false, false, false},
      runtime::DecompositionOptions{process_grid});
  mesh::MeshTopology topology(decomposition);
  runtime::FieldRegistry registry;
  const auto flux_field = finite_volume::declare_face_mass_flux(registry);
  registry.freeze();
  runtime::FieldStorage storage(
      registry, {decomposition.local_extent(), topology.local_face_count()});
  constexpr runtime::PhaseId phase = 4203U;
  constexpr runtime::ActorId actor = 4203U;
  runtime::FieldAccessPlan access(registry);
  access.declare_access(phase, actor, flux_field,
                        runtime::AccessMode::read_write);
  access.freeze();
  auto values = storage.acquire_face_write<double>(access, phase, actor,
                                                    flux_field);
  for (std::size_t face = 0; face < topology.local_face_count(); ++face)
    values(face, 0) = 0.0;

  const auto global_face = topology.global_face_id(
      mesh::LogicalFace{mesh::FaceAxis::x, {1, 0, 0}});
  const auto local_face = topology.find_local_face(global_face);
  HUNDUN_CHECK(local_face.has_value());
  HUNDUN_CHECK(topology.neighbour(*local_face).has_value());
  values(*local_face, 0) = 2.0;
  auto final_flux = flow::MaterialFaceMassFlux::acquire(
      registry, storage, access, phase, actor, flux_field, topology,
      flow::MaterialFluxProvenance::final_corrected);

  flow::detail::ReactingTransportRequest request;
  request.shared_face_mass_flux_field = flux_field;
  request.composition = composition();
  request.cells.resize(topology.local_cell_count());
  for (std::size_t cell = 0; cell < topology.local_cell_count(); ++cell) {
    auto &state = request.cells[cell];
    if (topology.global_cell(cell).x == 0) {
      state.density_kg_per_m3 = 1.0;
      state.rho_y_kg_per_m3 = {0.2, 0.3, 0.5};
      state.rho_h_tc_j_per_m3 = 10.0;
      state.thermodynamics = {300.0, 1.0, 2.0, 3.0};
    } else {
      state.density_kg_per_m3 = 2.0;
      state.rho_y_kg_per_m3 = {0.8, 0.2, 1.0};
      state.rho_h_tc_j_per_m3 = 28.0;
      state.thermodynamics = {400.0, 2.0, 2.0, 3.0};
    }
    state.molecular_transport =
        {1.0e-5, 2.0, {1.0e-5, 2.0e-5, 3.0e-5}};
  }
  request.interior_faces.push_back(
      {*local_face,
       topology.owner(*local_face),
       *topology.neighbour(*local_face),
       1.0,
       0.5,
       {0.25, 0.25, 0.5},
       {0.35, 0.15, 0.5},
       11.0,
       13.0});

  const auto report =
      flow::detail::assemble_reacting_transport(final_flux, request);
  double reductions[4]{};
  for (std::size_t cell = 0; cell < topology.local_cell_count(); ++cell) {
    if (topology.cell_ownership(cell) != mesh::EntityOwnership::owned)
      continue;
    for (std::size_t species = 0; species < 3U; ++species)
      reductions[species] += report.species_residual_kg_per_s[species][cell];
    reductions[3] += report.enthalpy_residual_w[cell];
  }
  mpi.allreduce_fp64_in_place(reductions, 4U,
                              runtime::Fp64ReductionOperation::sum);
  for (double value : reductions)
    HUNDUN_CHECK(value == 0.0);
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  return hundun::test::run([&] {
    auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
    HUNDUN_CHECK(mpi.size() == 1 || mpi.size() == 2);
    run(mpi);
  });
}
