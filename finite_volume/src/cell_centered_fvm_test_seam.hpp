// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/finite_volume/cell_centered_fvm.hpp"
#include "hundun/mesh/mesh_topology.hpp"

namespace hundun::finite_volume::test {

enum class TopologySignatureMutationForTest {
  cell_global_id,
  cell_ownership,
  face_ownership,
  face_owner_local_id,
  face_owner_global_id,
  face_owner_ownership,
  face_neighbour_presence,
  face_neighbour_local_id,
  face_neighbour_global_id,
  face_neighbour_ownership,
  logical_face,
  face_patch_membership,
  periodic_pair,
  patch_stable_id,
  patch_name,
  patch_pairing_kind,
  patch_paired_id,
  patch_exact_membership
};

enum class FaceMassFluxConstructionFailureForTest { bad_alloc, length_error };

void mutate_next_topology_signature(TopologySignatureMutationForTest mutation);

void override_next_face_metrics(mesh::GlobalFaceId global_face, double skewness,
                                double non_orthogonality_degrees);

void force_next_least_squares_singular();

void fail_next_face_mass_flux_construction(
    FaceMassFluxConstructionFailureForTest failure);

#ifdef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS
void reverse_scalar_diffusion_nonorthogonal_contribution_for_test(
    bool enabled);

class PreparedFaceMassFluxForTest final {
public:
  static PreparedFaceMassFluxForTest create(
      const mesh::MeshTopology &topology);

  ~PreparedFaceMassFluxForTest() noexcept;
  PreparedFaceMassFluxForTest(PreparedFaceMassFluxForTest &&) noexcept;
  PreparedFaceMassFluxForTest &
  operator=(PreparedFaceMassFluxForTest &&) = delete;
  PreparedFaceMassFluxForTest(const PreparedFaceMassFluxForTest &) = delete;
  PreparedFaceMassFluxForTest &
  operator=(const PreparedFaceMassFluxForTest &) = delete;

  FaceMassFlux bind(const runtime::FieldRegistry &registry,
                    const runtime::FieldStorage &storage,
                    const runtime::FieldAccessPlan &access_plan,
                    runtime::PhaseId phase, runtime::ActorId actor,
                    runtime::FieldId field,
                    const mesh::MeshTopology &topology);

private:
  explicit PreparedFaceMassFluxForTest(
      FaceMassFlux::PreparedStatePtr state) noexcept;
  FaceMassFlux::PreparedStatePtr state_;
};
#endif

} // namespace hundun::finite_volume::test
