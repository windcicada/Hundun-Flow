// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/mesh/mesh_topology.hpp"

namespace hundun::finite_volume::test {

void override_next_face_metrics(mesh::GlobalFaceId global_face, double skewness,
                                double non_orthogonality_degrees);

void force_next_least_squares_singular();

} // namespace hundun::finite_volume::test
