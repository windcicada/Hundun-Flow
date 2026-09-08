// SPDX-License-Identifier: Apache-2.0
#pragma once
#include "hundun/v04_mesh.hpp"
#include "models_exchange_batch_detail.hpp"
#include <mpi.h>

namespace hundun::v04::portable {
// Production, bounded owner routing. First audit each complete physical segment
// at a deterministic auditor, then deliver rows only to the actual cell owner.
// No global cell census or global segment replication. The communicator is
// borrowed. All C++ buffers are allocated by configure, never by prepare.
class OwnerExchangeRoutingPlan {
public:
  OwnerExchangeRoutingPlan() = default;
  OwnerExchangeRoutingPlan(const OwnerExchangeRoutingPlan &) = delete;
  OwnerExchangeRoutingPlan &
  operator=(const OwnerExchangeRoutingPlan &) = delete;
  hundun::v04::Status configure(MPI_Comm, Int3 global, MeshPatch,
                                std::size_t maximum_local_segments,
                                std::size_t maximum_auditor_segments,
                                std::size_t species_count,
                                std::uint64_t maximum_bytes) noexcept;
  Status prepare(Revision, const ExchangeSegment *, std::size_t) noexcept;
  RoutedExchangeSegments candidates() const noexcept;
  void discard() noexcept {
    ++generation_;
    available_ = false;
  }
  std::uint64_t owned_bytes() const noexcept { return owned_bytes_; }

private:
  OwnerExchangeRoutingPlan &operator=(OwnerExchangeRoutingPlan &&) = default;
  int owner(std::uint64_t) const noexcept;
  int auditor(const ExchangeSegment &) const noexcept;
  Status agree(Status) const noexcept;
  Status exchange(const ExchangeSegment *, std::size_t, bool,
                  std::size_t &) noexcept;
  MPI_Comm comm_{MPI_COMM_NULL};
  Int3 global_{}, grid_{};
  int rank_{}, ranks_{};
  std::size_t local_capacity_{}, audit_capacity_{}, species_{};
  std::uint64_t owned_bytes_{}, generation_{};
  Revision revision_{};
  std::vector<int> send_counts_, receive_counts_, send_offsets_,
      receive_offsets_, cursor_;
  std::vector<std::uint64_t> send_, receive_;
  std::vector<ExchangeSegment> values_;
  std::size_t count_{};
  bool available_{};
};
} // namespace hundun::v04::portable
