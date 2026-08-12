// SPDX-License-Identifier: Apache-2.0

#include "flow_reacting_immersed_detail.hpp"

#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "tests/support/test_main.hpp"

#include <cstring>
#include <stdexcept>

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  return hundun::test::run([&] {
    auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
    hundun::flow::detail::ReactingImmersedRequest request;
    request.authority.geometry_fingerprint = 101U;
    request.authority.wall_ownership_fingerprint = 202U;
    request.authority.thermodynamic_pressure_pa = 101325.0;
    request.authority.shared_face_mass_flux_field = {7U};
    request.authority.flux_provenance =
        hundun::flow::MaterialFluxProvenance::final_corrected;
    request.expected_face_mass_flux_field = {7U};
    request.expected_thermodynamic_pressure_pa = 101325.0;
    request.wall.species_flux_kg_per_m2_s = {0.0, 0.0, 0.0};
    request.wall.heat_flux_w_per_m2 = 0.0;
    request.molecular_diffusivity_m2_per_s = {1.0e-5, 2.0e-5, 3.0e-5};

    const auto molecular =
        hundun::flow::detail::compose_reacting_immersed_coupling(request);
    HUNDUN_CHECK(std::memcmp(molecular.effective_diffusivity_m2_per_s.data(),
                             request.molecular_diffusivity_m2_per_s.data(),
                             3U * sizeof(double)) == 0);
    HUNDUN_CHECK(molecular.wall_source.kind ==
                 hundun::flow::detail::ReactingSourceKind::boundary);

    request.authority.wale_identity = hundun::les::WaleCoefficientIdentity{9U};
    request.authority.wale_evaluation_count = 1U;
    request.wale_kinematic_diffusivity_m2_per_s = 4.0e-5;
    request.turbulent_schmidt = 2.0;
    const auto wale =
        hundun::flow::detail::compose_reacting_immersed_coupling(request);
    HUNDUN_CHECK_NEAR(wale.effective_diffusivity_m2_per_s[0], 3.0e-5,
                      1.0e-18);
    HUNDUN_CHECK(wale.wale_identity == request.authority.wale_identity);

    bool rejected = false;
    auto mutation = request;
    mutation.authority.wale_evaluation_count = 3U;
    try {
      static_cast<void>(
          hundun::flow::detail::compose_reacting_immersed_coupling(mutation));
    } catch (const std::invalid_argument &) {
      rejected = true;
    }
    HUNDUN_CHECK(rejected);
    mutation = request;
    mutation.authority.flux_provenance =
        hundun::flow::MaterialFluxProvenance::predictor;
    rejected = false;
    try {
      static_cast<void>(
          hundun::flow::detail::compose_reacting_immersed_coupling(mutation));
    } catch (const std::invalid_argument &) {
      rejected = true;
    }
    HUNDUN_CHECK(rejected);

    double fingerprint = static_cast<double>(
        wale.geometry_fingerprint + wale.wall_ownership_fingerprint);
    mpi.allreduce_fp64_in_place(
        &fingerprint, 1U, hundun::runtime::Fp64ReductionOperation::sum);
    HUNDUN_CHECK(fingerprint == 303.0 * static_cast<double>(mpi.size()));
  });
}
