// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/fvm_cell_centered.hpp"
#include "hundun/mesh_topology.hpp"

#include <utility>

namespace hundun::finite_volume::detail {

void mutate_next_topology_signature_raw(int mutation);
void override_next_face_metrics_raw(mesh::GlobalFaceId global_face,
                                    double skewness,
                                    double non_orthogonality_degrees);
void force_next_least_squares_singular_raw();
void fail_next_face_mass_flux_construction_raw(int failure);
void reverse_scalar_diffusion_nonorthogonal_contribution_raw(bool enabled);

} // namespace hundun::finite_volume::detail

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

inline void
mutate_next_topology_signature(TopologySignatureMutationForTest mutation) {
  detail::mutate_next_topology_signature_raw(static_cast<int>(mutation));
}

inline void override_next_face_metrics(
    mesh::GlobalFaceId global_face, double skewness,
    double non_orthogonality_degrees) {
  detail::override_next_face_metrics_raw(global_face, skewness,
                                         non_orthogonality_degrees);
}

inline void force_next_least_squares_singular() {
  detail::force_next_least_squares_singular_raw();
}

inline void fail_next_face_mass_flux_construction(
    FaceMassFluxConstructionFailureForTest failure) {
  detail::fail_next_face_mass_flux_construction_raw(
      static_cast<int>(failure));
}

inline void reverse_scalar_diffusion_nonorthogonal_contribution_for_test(
    bool enabled) {
  detail::reverse_scalar_diffusion_nonorthogonal_contribution_raw(enabled);
}

class PreparedFaceMassFluxForTest final {
public:
  static PreparedFaceMassFluxForTest create(
      const mesh::MeshTopology &topology) {
    return PreparedFaceMassFluxForTest(FaceMassFlux::prepare(topology));
  }

  ~PreparedFaceMassFluxForTest() noexcept = default;
  PreparedFaceMassFluxForTest(PreparedFaceMassFluxForTest &&) noexcept =
      default;
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
                    const mesh::MeshTopology &topology) {
    if (!state_)
      throw runtime::Error(
          "prepared face mass flux test handle has been moved from");
    return FaceMassFlux::bind_prepared(*state_, registry, storage, access_plan,
                                       phase, actor, field, topology);
  }

private:
  explicit PreparedFaceMassFluxForTest(
      FaceMassFlux::PreparedStatePtr state) noexcept
      : state_(std::move(state)) {}

  FaceMassFlux::PreparedStatePtr state_;
};

} // namespace hundun::finite_volume::test
