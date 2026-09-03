// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "app_evidence_detail.hpp"

#include <algorithm>
#include <array>
#include <cstdint>

namespace hundun::v04::detail {
namespace {

constexpr std::uint32_t kRuntimeEvidenceFailure = 10506U;

RuntimeConvectiveCflWinner runtime_cfl_winner(
    const ConvectiveCflFailureWitness& witness) noexcept {
  return {witness.valid,
          witness.global_cell,
          static_cast<std::int32_t>(witness.rank),
          witness.out,
          witness.absolute,
          witness.density_volume,
          witness.outgoing_mass_flow,
          witness.absolute_mass_flow};
}

void fingerprint_word(std::uint64_t word, std::uint64_t& hash) noexcept {
  constexpr std::uint64_t kPrime = UINT64_C(1099511628211);
  for (unsigned byte = 0U; byte < 8U; ++byte) {
    hash ^= static_cast<std::uint8_t>(word >> (8U * byte));
    hash *= kPrime;
  }
}

template <class Field>
void fingerprint_field_view(const Field& field,
                            std::uint64_t& hash) noexcept {
  fingerprint_word(static_cast<std::uint32_t>(field.interior.x), hash);
  fingerprint_word(static_cast<std::uint32_t>(field.interior.y), hash);
  fingerprint_word(static_cast<std::uint32_t>(field.interior.z), hash);
  fingerprint_word(static_cast<std::uint32_t>(field.ghosts.x), hash);
  fingerprint_word(static_cast<std::uint32_t>(field.ghosts.y), hash);
  fingerprint_word(static_cast<std::uint32_t>(field.ghosts.z), hash);
  fingerprint_word(field.components, hash);
  fingerprint_word(field.stride_y, hash);
  fingerprint_word(field.stride_z, hash);
  fingerprint_word(field.component_stride, hash);
  fingerprint_word(field.replica, hash);
  fingerprint_word(field.field, hash);
  fingerprint_word(field.revision, hash);
  fingerprint_word(field.storage_identity, hash);
  fingerprint_word(field.revision_domain, hash);
}

template <class Face>
void fingerprint_face_view(const Face& face, std::uint64_t& hash) noexcept {
  fingerprint_word(static_cast<std::uint32_t>(face.extents.x), hash);
  fingerprint_word(static_cast<std::uint32_t>(face.extents.y), hash);
  fingerprint_word(static_cast<std::uint32_t>(face.extents.z), hash);
  fingerprint_word(face.stride_y, hash);
  fingerprint_word(face.stride_z, hash);
  fingerprint_word(static_cast<std::uint8_t>(face.axis), hash);
  fingerprint_word(face.storage_identity, hash);
  fingerprint_word(face.revision_domain, hash);
}

void fingerprint_face_flux_view(ConstFaceFluxView flux,
                                std::uint64_t& hash) noexcept {
  fingerprint_word(flux.revision, hash);
  fingerprint_word(flux.certificate.revision(), hash);
  fingerprint_word(flux.certificate.authority(), hash);
  fingerprint_word(flux.certificate.storage(), hash);
  fingerprint_word(flux.certificate.revision_domain(), hash);
  fingerprint_face_view(flux.x, hash);
  fingerprint_face_view(flux.y, hash);
  fingerprint_face_view(flux.z, hash);
}

}  // namespace

Status runtime_committed_cfl(
    MPI_Comm communicator,
    const CommittedConvectiveCflCertificate& certificate,
    RuntimeCommittedConvectiveCflAudit& runtime) noexcept {
  if (communicator == MPI_COMM_NULL || !certificate.valid())
    return {StatusCode::invalid_plan, kRuntimeEvidenceFailure};
  int rank = 0;
  if (MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kRuntimeEvidenceFailure};
  constexpr std::uint64_t kOffset = UINT64_C(1469598103934665603);
  std::array<std::uint64_t, 2U> collective{kOffset, kOffset};
  fingerprint_word(static_cast<std::uint32_t>(rank), collective[0U]);
  fingerprint_word(static_cast<std::uint32_t>(rank), collective[1U]);
  fingerprint_field_view(certificate.density_view_identity.view,
                         collective[0U]);
  const ConstFaceFluxView flux = certificate.face_flux_view_identity.view;
  fingerprint_face_flux_view(flux, collective[1U]);
  if (MPI_Allreduce(MPI_IN_PLACE, collective.data(),
                    static_cast<int>(collective.size()), MPI_UINT64_T,
                    MPI_BXOR, communicator) != MPI_SUCCESS ||
      collective[0U] == 0U || collective[1U] == 0U)
    return {StatusCode::mpi_failure, kRuntimeEvidenceFailure};
  runtime = {};
  runtime.valid = certificate.valid();
  runtime.density_revision = certificate.density;
  runtime.final_flux_revision = certificate.final_flux;
  runtime.density_field = certificate.density_view_identity.view.field;
  runtime.density_view_collective = collective[0U];
  runtime.final_flux_view_collective = collective[1U];
  runtime.activity_collective = certificate.activity_collective;
  runtime.dt = certificate.dt;
  runtime.out_max = certificate.out_max;
  runtime.abs_max = certificate.absolute_max;
  runtime.limit = certificate.limit;
  runtime.out_winner = runtime_cfl_winner(certificate.out_winner);
  runtime.abs_winner = runtime_cfl_winner(certificate.absolute_winner);
  return {};
}

Status runtime_advective_cfl(
    MPI_Comm communicator,
    const MomentumAdvectiveCflCertificate& certificate,
    RuntimeAdvectiveCflAudit& runtime) noexcept {
  if (communicator == MPI_COMM_NULL || !certificate.valid())
    return {StatusCode::invalid_plan, kRuntimeEvidenceFailure};
  int rank = 0;
  if (MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kRuntimeEvidenceFailure};
  constexpr std::uint64_t kOffset = UINT64_C(1469598103934665603);
  std::array<std::uint64_t, 3U> collective{kOffset, kOffset, kOffset};
  for (std::uint64_t& fingerprint : collective)
    fingerprint_word(static_cast<std::uint32_t>(rank), fingerprint);
  fingerprint_word(certificate.time, collective[0U]);
  fingerprint_field_view(certificate.density_view_identity.view,
                         collective[1U]);
  const ConstFaceFluxView flux = certificate.face_flux_view_identity.view;
  fingerprint_face_flux_view(flux, collective[2U]);
  if (MPI_Allreduce(MPI_IN_PLACE, collective.data(),
                    static_cast<int>(collective.size()), MPI_UINT64_T,
                    MPI_BXOR, communicator) != MPI_SUCCESS ||
      std::any_of(collective.begin(), collective.end(),
                  [](std::uint64_t value) { return value == 0U; }))
    return {StatusCode::mpi_failure, kRuntimeEvidenceFailure};
  runtime = {true,
             certificate.plan,
             collective[0U],
             collective[1U],
             collective[2U],
             certificate.activity_collective,
             certificate.dt,
             certificate.out_max,
             certificate.absolute_max,
             certificate.limit};
  return {};
}

}  // namespace hundun::v04::detail
