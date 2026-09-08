// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#pragma once

#include "hundun/v04_mesh.hpp"
#include "hundun/v04_spray.hpp"

#include <mpi.h>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace hundun::v04::spray::detail {

struct ParcelMigrationValue {
  SprayParcelState parcel{};
  double tab_deformation{};
  double tab_deformation_rate_per_s{};
};

enum class MigrationDetail : std::uint32_t {
  invalid_plan = 26100U,
  invalid_value,
  invalid_owner,
  capacity_exceeded,
  duplicate_id,
  communication_failure,
  allocation_failure
};

struct ParcelMigrationReport {
  Status status{StatusCode::invalid_plan,
                static_cast<std::uint32_t>(MigrationDetail::invalid_plan)};
  bool available{};
  int lowest_failing_rank{-1};
  std::uint64_t input_parcels{};
  std::uint64_t output_parcels{};
  std::uint64_t outgoing_parcels{};
  std::uint64_t incoming_parcels{};
  std::uint64_t global_parcels{};
};

// Variable-size ownership exchange, separate from fixed-field Halo. configure
// allocates every buffer; prepare performs no C++ allocation. The communicator
// is borrowed and must remain alive until the plan is no longer used.
// Candidates are borrowed until the next prepare/discard/configure. This class
// never changes a committed parcel container or publishes a gas source.
// owner_global_cell must already be resolved from the physical endpoint by the
// mesh/trajectory module. This exchange validates its Cartesian range/owner,
// not position-to-cell consistency; physical geometry is not an input here.
class ParcelMigrationPlan {
 public:
  ParcelMigrationPlan() noexcept = default;
  ParcelMigrationPlan(const ParcelMigrationPlan&) = delete;
  ParcelMigrationPlan& operator=(const ParcelMigrationPlan&) = delete;
  ParcelMigrationPlan(ParcelMigrationPlan&&) = delete;
  ParcelMigrationPlan& operator=(ParcelMigrationPlan&&) = delete;

  Status configure(MPI_Comm communicator, Int3 global_cells, MeshPatch patch,
                   std::size_t maximum_local_parcels,
                   std::size_t maximum_ids_per_auditor) noexcept;
  Status prepare(Span<const ParcelMigrationValue> trial,
                 ParcelMigrationReport& report) noexcept;
  Span<const ParcelMigrationValue> candidates() const noexcept {
    return available_ ? Span<const ParcelMigrationValue>{candidates_.data(),
                                                         candidates_.size()}
                      : Span<const ParcelMigrationValue>{};
  }
  void discard() noexcept { available_ = false; candidates_.clear(); }

 private:
  int owner(std::uint64_t global_cell) const noexcept;
  int auditor(ParcelId id) const noexcept;
  Status consensus(Status local, int& failing_rank) const noexcept;
  Status counts(std::size_t receive_capacity, int lanes,
                std::size_t& receive_count, int& failing_rank) noexcept;

  MPI_Comm communicator_{MPI_COMM_NULL};
  Int3 global_cells_{};
  Int3 process_grid_{};
  int rank_{};
  int ranks_{};
  std::size_t local_capacity_{};
  std::size_t audit_capacity_{};
  std::vector<int> send_counts_, receive_counts_, send_offsets_, receive_offsets_, cursor_;
  std::vector<int> destinations_;
  std::vector<std::uint64_t> send_wire_, receive_wire_, audit_wire_;
  std::vector<ParcelId> audit_ids_;
  std::vector<ParcelMigrationValue> candidates_;
  bool available_{};
};

}  // namespace hundun::v04::spray::detail
