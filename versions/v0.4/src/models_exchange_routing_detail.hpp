// SPDX-License-Identifier: Apache-2.0
#pragma once
#include "models_exchange_batch_detail.hpp"
#include <mpi.h>
namespace hundun::v04::portable {
struct ExchangeRoutingReport {
  Status status{Status::invalid_input};
  bool available{};
  const ExchangeSegment *segments{};
  std::size_t count{}, global_count{};
};
// Bounded reference routing of ordinary source candidates. This is not field
// Halo or a production coupling stencil. Each global cell has exactly one
// owner in the supplied owned-cell lists; only that owner receives its sources.
// prepare uses fixed HUNDUN buffers. MPI's internal resource policy is
// separate.
class ExchangeRoutingWorkspace {
public:
  ExchangeRoutingWorkspace(std::size_t local_cell_capacity,
                           std::size_t global_segment_capacity,
                           std::size_t rank_capacity = 4);
  ExchangeRoutingReport prepare(MPI_Comm, const ExchangeCell *, std::size_t,
                                const ExchangeSegment *, std::size_t) noexcept;
  // Includes removed parents and every created child, not only survivors.
  Status audit_ids(MPI_Comm, const spray::ParcelId *, std::size_t) noexcept;

private:
  std::size_t cell_capacity_{}, segment_capacity_{}, rank_capacity_{};
  std::vector<int> counts_, offsets_;
  std::vector<std::uint64_t> local_ids_, global_ids_, send_, receive_;
  std::vector<ExchangeSegment> candidates_;
};
} // namespace hundun::v04::portable
