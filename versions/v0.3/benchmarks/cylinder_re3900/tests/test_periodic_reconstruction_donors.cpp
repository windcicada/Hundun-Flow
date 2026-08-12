// SPDX-License-Identifier: Apache-2.0

#include "ib_quadratic_reconstruction_detail.hpp"
#include "ib_deterministic_qr_detail.hpp"

#include "hundun/ib_quadratic_reconstruction.hpp"
#include "hundun/mesh_geometry.hpp"
#include "hundun/mesh_topology.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_structured_decomposition.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <algorithm>
#include <vector>

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  {
    auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
    HUNDUN_CHECK(mpi.size() == 1);
    auto decomposition = hundun::runtime::StructuredDecomposition::create(
        mpi, {8, 8, 8}, {false, false, true},
        hundun::runtime::DecompositionOptions{
            hundun::runtime::Int3{1, 1, 1}});
    hundun::mesh::MeshTopology topology(decomposition);
    hundun::mesh::MeshGeometry geometry(
        topology,
        hundun::mesh::UniformBoxMapping({0.0, 0.0, 0.0}, {8.0, 8.0, 8.0}));

    std::vector<hundun::runtime::Int3> donors;
    for (int k : {7, 0, 1})
      for (int j = 2; j <= 4; ++j)
        for (int i = 2; i <= 4; ++i)
          donors.push_back({i, j, k});

    const auto reconstruction =
        hundun::immersed::QuadraticReconstruction::create(
            {1.5, 3.5, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0},
            {0.0, 0.0, 1.0}, 1.0, {2, 3, 0}, donors, topology, geometry);
    const auto &logical_donors =
        hundun::immersed::detail::QuadraticReconstructionWeights::
            donor_stencil_images(reconstruction);

    HUNDUN_CHECK(std::any_of(logical_donors.begin(), logical_donors.end(),
                             [](const auto donor) { return donor.z == -1; }));
    HUNDUN_CHECK(std::none_of(logical_donors.begin(), logical_donors.end(),
                              [](const auto donor) { return donor.z == 7; }));
  }
  MPI_Finalize();
}
