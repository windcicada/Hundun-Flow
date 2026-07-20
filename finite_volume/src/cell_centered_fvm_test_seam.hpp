// SPDX-License-Identifier: Apache-2.0
#pragma once

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

} // namespace hundun::finite_volume::test
