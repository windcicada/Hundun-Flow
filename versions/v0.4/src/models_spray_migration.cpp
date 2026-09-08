// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat:
// windcicada | Year.M: 2026.09

#include "mesh_focus_detail.hpp"
#include "models_spray_migration_detail.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <new>
#include <utility>

namespace hundun::v04::spray::detail {
namespace {
constexpr int kLanes = 18;
static_assert(sizeof(double) == sizeof(std::uint64_t) &&
                  std::numeric_limits<double>::is_iec559,
              "parcel wire requires IEEE binary64");

Status invalid(MigrationDetail detail) noexcept {
  return {StatusCode::invalid_plan, static_cast<std::uint32_t>(detail)};
}
Status mpi_failure() noexcept {
  return {StatusCode::mpi_failure,
          static_cast<std::uint32_t>(MigrationDetail::communication_failure)};
}
bool same(Int3 a, Int3 b) noexcept {
  return a.x == b.x && a.y == b.y && a.z == b.z;
}
Status collective(MPI_Comm comm, int rank, Status local,
                  int &failing_rank) noexcept {
  int failing = local ? std::numeric_limits<int>::max() : rank;
  int first = 0;
  if (MPI_Allreduce(&failing, &first, 1, MPI_INT, MPI_MIN, comm) != MPI_SUCCESS)
    return mpi_failure();
  failing_rank = first == std::numeric_limits<int>::max() ? -1 : first;
  if (failing_rank < 0)
    return {};
  std::uint32_t wire[]{static_cast<std::uint32_t>(local.code), local.detail};
  if (MPI_Bcast(wire, 2, MPI_UINT32_T, first, comm) != MPI_SUCCESS)
    return mpi_failure();
  return {static_cast<StatusCode>(wire[0]), wire[1]};
}
std::uint64_t bits(double value) noexcept {
  std::uint64_t out;
  std::memcpy(&out, &value, sizeof(out));
  return out;
}
double real(std::uint64_t value) noexcept {
  double out;
  std::memcpy(&out, &value, sizeof(out));
  return out;
}
void encode(const ParcelMigrationValue &value, std::uint64_t *wire) noexcept {
  const auto &p = value.parcel;
  wire[0] = p.id.high;
  wire[1] = p.id.low;
  for (int c = 0; c < 3; ++c) {
    wire[2 + c] = bits(p.position_m[c]);
    wire[5 + c] = bits(p.velocity_m_per_s[c]);
  }
  wire[8] = bits(p.droplet_mass_kg);
  wire[9] = bits(p.droplet_diameter_m);
  wire[10] = bits(p.multiplicity);
  wire[11] = bits(p.temperature_k);
  wire[12] = p.liquid_material_fingerprint;
  wire[13] = p.owner_global_cell;
  wire[14] = bits(p.age_s);
  wire[15] = bits(value.tab_deformation);
  wire[16] = bits(value.tab_deformation_rate_per_s);
  wire[17] = value.breakup_ordinal;
}
ParcelMigrationValue decode(const std::uint64_t *wire) noexcept {
  ParcelMigrationValue value;
  auto &p = value.parcel;
  p.id = {wire[0], wire[1]};
  for (int c = 0; c < 3; ++c) {
    p.position_m[c] = real(wire[2 + c]);
    p.velocity_m_per_s[c] = real(wire[5 + c]);
  }
  p.droplet_mass_kg = real(wire[8]);
  p.droplet_diameter_m = real(wire[9]);
  p.multiplicity = real(wire[10]);
  p.temperature_k = real(wire[11]);
  p.liquid_material_fingerprint = wire[12];
  p.owner_global_cell = wire[13];
  p.age_s = real(wire[14]);
  value.tab_deformation = real(wire[15]);
  value.tab_deformation_rate_per_s = real(wire[16]);
  value.breakup_ordinal = wire[17];
  return value;
}
bool valid(const ParcelMigrationValue &value) noexcept {
  return validate_parcel_state(value.parcel) == ParcelStateStatus::success &&
         std::isfinite(value.tab_deformation) &&
         std::isfinite(value.tab_deformation_rate_per_s);
}
int owner_coordinate(int index, int cells, int partitions) noexcept {
  const int ordinary = cells / partitions;
  const int extra = cells % partitions;
  const std::int64_t extended = static_cast<std::int64_t>(ordinary + 1) * extra;
  return index < extended
             ? index / (ordinary + 1)
             : extra + static_cast<int>((index - extended) / ordinary);
}
} // namespace

Status ParcelMigrationPlan::configure(MPI_Comm communicator, Int3 cells,
                                      MeshPatch patch,
                                      std::size_t maximum_local,
                                      std::size_t maximum_audit) noexcept {
  if (communicator == MPI_COMM_NULL)
    return invalid(MigrationDetail::invalid_plan);
  int rank = 0, ranks = 0, failing_rank = -1;
  if (MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(communicator, &ranks) != MPI_SUCCESS)
    return mpi_failure();
  MeshPatch authoritative;
  Status local =
      hundun::v04::detail::make_mesh_patch(rank, ranks, cells, authoritative);
  if (!local || !same(patch.begin, authoritative.begin) ||
      !same(patch.cells, authoritative.cells) ||
      !same(patch.process_grid, authoritative.process_grid) ||
      !same(patch.process_coord, authoritative.process_coord) ||
      maximum_local > static_cast<std::size_t>(INT32_MAX / kLanes) ||
      maximum_audit > static_cast<std::size_t>(INT32_MAX / 2))
    local = invalid(MigrationDetail::invalid_plan);
  Status status = collective(communicator, rank, local, failing_rank);
  if (!status)
    return status;
  int dimensions[]{cells.x, cells.y, cells.z}, minima[3]{}, maxima[3]{};
  if (MPI_Allreduce(dimensions, minima, 3, MPI_INT, MPI_MIN, communicator) !=
          MPI_SUCCESS ||
      MPI_Allreduce(dimensions, maxima, 3, MPI_INT, MPI_MAX, communicator) !=
          MPI_SUCCESS)
    return mpi_failure();
  if (!std::equal(minima, minima + 3, maxima))
    return invalid(MigrationDetail::invalid_plan);

  // Build a cold candidate so failed reconfiguration preserves the old plan.
  std::vector<int> send_counts, receive_counts, send_offsets, receive_offsets,
      cursor, destinations;
  std::vector<std::uint64_t> send_wire, receive_wire, audit_wire;
  std::vector<ParcelId> audit_ids;
  std::vector<ParcelMigrationValue> candidates;
  try {
    send_counts.resize(ranks);
    receive_counts.resize(ranks);
    send_offsets.resize(ranks);
    receive_offsets.resize(ranks);
    cursor.resize(ranks);
    destinations.resize(maximum_local);
    // One sentinel word also makes an empty MPI buffer a valid address.
    send_wire.resize(std::max<std::size_t>(1, maximum_local * kLanes));
    receive_wire.resize(std::max<std::size_t>(1, maximum_local * kLanes));
    audit_wire.resize(std::max<std::size_t>(1, maximum_audit * 2));
    audit_ids.reserve(maximum_audit);
    candidates.reserve(maximum_local);
  } catch (...) {
    local = {StatusCode::allocation_failure,
             static_cast<std::uint32_t>(MigrationDetail::allocation_failure)};
  }
  status = collective(communicator, rank, local, failing_rank);
  if (!status)
    return status;
  send_counts_.swap(send_counts);
  receive_counts_.swap(receive_counts);
  send_offsets_.swap(send_offsets);
  receive_offsets_.swap(receive_offsets);
  cursor_.swap(cursor);
  destinations_.swap(destinations);
  send_wire_.swap(send_wire);
  receive_wire_.swap(receive_wire);
  audit_wire_.swap(audit_wire);
  audit_ids_.swap(audit_ids);
  candidates_.swap(candidates);
  communicator_ = communicator;
  rank_ = rank;
  ranks_ = ranks;
  global_cells_ = cells;
  process_grid_ = patch.process_grid;
  local_capacity_ = maximum_local;
  audit_capacity_ = maximum_audit;
  available_ = false;
  return {};
}

int ParcelMigrationPlan::owner(std::uint64_t cell) const noexcept {
  const std::uint64_t x = cell % static_cast<std::uint64_t>(global_cells_.x);
  cell /= static_cast<std::uint64_t>(global_cells_.x);
  const std::uint64_t y = cell % static_cast<std::uint64_t>(global_cells_.y);
  const std::uint64_t z = cell / static_cast<std::uint64_t>(global_cells_.y);
  if (z >= static_cast<std::uint64_t>(global_cells_.z))
    return -1;
  return owner_coordinate(static_cast<int>(x), global_cells_.x,
                          process_grid_.x) +
         process_grid_.x *
             (owner_coordinate(static_cast<int>(y), global_cells_.y,
                               process_grid_.y) +
              process_grid_.y * owner_coordinate(static_cast<int>(z),
                                                 global_cells_.z,
                                                 process_grid_.z));
}
int ParcelMigrationPlan::auditor(ParcelId id) const noexcept {
  // Routing only; the full 128-bit key, never this hash, establishes
  // uniqueness.
  return static_cast<int>((id.high ^ id.low) %
                          static_cast<std::uint64_t>(ranks_));
}
Status ParcelMigrationPlan::consensus(Status local,
                                      int &failing_rank) const noexcept {
  return collective(communicator_, rank_, local, failing_rank);
}
Status ParcelMigrationPlan::counts(std::size_t capacity, int lanes,
                                   std::size_t &received,
                                   int &failing_rank) noexcept {
  if (MPI_Alltoall(send_counts_.data(), 1, MPI_INT, receive_counts_.data(), 1,
                   MPI_INT, communicator_) != MPI_SUCCESS)
    return mpi_failure();
  std::size_t sends = 0, receives = 0;
  Status local;
  for (int rank = 0; rank < ranks_; ++rank) {
    const int s = send_counts_[rank], r = receive_counts_[rank];
    if (s < 0 || r < 0 ||
        static_cast<std::size_t>(r) > capacity - std::min(capacity, receives)) {
      local = invalid(MigrationDetail::capacity_exceeded);
      break;
    }
    send_offsets_[rank] = static_cast<int>(sends * lanes);
    receive_offsets_[rank] = static_cast<int>(receives * lanes);
    sends += static_cast<std::size_t>(s);
    receives += static_cast<std::size_t>(r);
    if (sends > static_cast<std::size_t>(INT32_MAX / lanes) ||
        receives > capacity ||
        receives > static_cast<std::size_t>(INT32_MAX / lanes)) {
      local = invalid(MigrationDetail::capacity_exceeded);
      break;
    }
    send_counts_[rank] *= lanes;
    receive_counts_[rank] *= lanes;
    cursor_[rank] = send_offsets_[rank];
  }
  Status status = consensus(local, failing_rank);
  if (status)
    received = receives;
  return status;
}

Status ParcelMigrationPlan::prepare(Span<const ParcelMigrationValue> input,
                                    ParcelMigrationReport &report) noexcept {
  discard();
  report = {};
  if (communicator_ == MPI_COMM_NULL)
    return report.status;
  Status local;
  if ((input.size != 0 && input.data == nullptr) ||
      input.size > local_capacity_)
    local = invalid(MigrationDetail::capacity_exceeded);
  std::fill(send_counts_.begin(), send_counts_.end(), 0);
  std::uint64_t outgoing = 0;
  if (local)
    for (std::size_t i = 0; i < input.size; ++i) {
      if (!valid(input.data[i])) {
        local = invalid(MigrationDetail::invalid_value);
        break;
      }
      const int destination = owner(input.data[i].parcel.owner_global_cell);
      if (destination < 0) {
        local = invalid(MigrationDetail::invalid_owner);
        break;
      }
      destinations_[i] = destination;
      ++send_counts_[destination];
      outgoing += destination != rank_;
    }
  int failing_rank = -1;
  auto fail = [&](Status status) noexcept {
    discard();
    report = {};
    report.status = status;
    report.lowest_failing_rank = failing_rank;
    return status;
  };
  Status status = consensus(local, failing_rank);
  if (!status)
    return fail(status);
  std::size_t received = 0;
  status = counts(local_capacity_, kLanes, received, failing_rank);
  if (!status)
    return fail(status);
  const std::uint64_t incoming =
      received - static_cast<std::size_t>(receive_counts_[rank_] / kLanes);
  for (std::size_t i = 0; i < input.size; ++i) {
    int &offset = cursor_[destinations_[i]];
    encode(input.data[i], send_wire_.data() + offset);
    offset += kLanes;
  }
  if (MPI_Alltoallv(send_wire_.data(), send_counts_.data(),
                    send_offsets_.data(), MPI_UINT64_T, receive_wire_.data(),
                    receive_counts_.data(), receive_offsets_.data(),
                    MPI_UINT64_T, communicator_) != MPI_SUCCESS)
    return fail(mpi_failure());
  candidates_.resize(received);
  std::fill(send_counts_.begin(), send_counts_.end(), 0);
  for (std::size_t i = 0; i < received; ++i) {
    candidates_[i] = decode(receive_wire_.data() + i * kLanes);
    if (!valid(candidates_[i]) ||
        owner(candidates_[i].parcel.owner_global_cell) != rank_) {
      local = invalid(MigrationDetail::invalid_owner);
      break;
    }
    ++send_counts_[auditor(candidates_[i].parcel.id)];
  }
  status = consensus(local, failing_rank);
  if (!status)
    return fail(status);
  std::size_t audited = 0;
  status = counts(audit_capacity_, 2, audited, failing_rank);
  if (!status)
    return fail(status);
  for (const auto &candidate : candidates_) {
    const ParcelId id = candidate.parcel.id;
    int &offset = cursor_[auditor(id)];
    send_wire_[offset++] = id.high;
    send_wire_[offset++] = id.low;
  }
  if (MPI_Alltoallv(send_wire_.data(), send_counts_.data(),
                    send_offsets_.data(), MPI_UINT64_T, audit_wire_.data(),
                    receive_counts_.data(), receive_offsets_.data(),
                    MPI_UINT64_T, communicator_) != MPI_SUCCESS)
    return fail(mpi_failure());
  audit_ids_.resize(audited);
  for (std::size_t i = 0; i < audited; ++i)
    audit_ids_[i] = {audit_wire_[2 * i], audit_wire_[2 * i + 1]};
  std::sort(audit_ids_.begin(), audit_ids_.end());
  for (std::size_t i = 1; i < audited; ++i)
    if (audit_ids_[i - 1] == audit_ids_[i]) {
      local = invalid(MigrationDetail::duplicate_id);
      break;
    }
  status = consensus(local, failing_rank);
  if (!status)
    return fail(status);
  std::uint64_t local_totals[]{static_cast<std::uint64_t>(input.size),
                               static_cast<std::uint64_t>(received)};
  std::uint64_t global_totals[2]{};
  if (MPI_Allreduce(local_totals, global_totals, 2, MPI_UINT64_T, MPI_SUM,
                    communicator_) != MPI_SUCCESS)
    return fail(mpi_failure());
  if (global_totals[0] != global_totals[1])
    return fail(invalid(MigrationDetail::invalid_owner));
  std::sort(
      candidates_.begin(), candidates_.end(),
      [](const auto &a, const auto &b) { return a.parcel.id < b.parcel.id; });
  available_ = true;
  report.status = {};
  report.available = true;
  report.input_parcels = input.size;
  report.output_parcels = received;
  report.outgoing_parcels = outgoing;
  report.incoming_parcels = incoming;
  report.global_parcels = global_totals[0];
  return {};
}

Status
ParcelMigrationPlan::prepare_checked(Span<const ParcelMigrationValue> input,
                                     const ParcelLocationProvider &location,
                                     portable::Revision revision,
                                     ParcelMigrationReport &report) noexcept {
  discard();
  report = {};
  if (communicator_ == MPI_COMM_NULL)
    return report.status;
  Status local;
  if ((input.size && !input.data) || input.size > local_capacity_ ||
      revision.algorithm_version != 1U)
    local = invalid(MigrationDetail::invalid_plan);
  // All ranks must request the same accepted clock/geometry revision before
  // invoking geometry; no rank advances to the exchange on local success only.
  std::uint64_t clock[]{revision.accepted_step, revision.input_revision,
                        revision.algorithm_version},
      minimum[3]{}, maximum[3]{};
  if (MPI_Allreduce(clock, minimum, 3, MPI_UINT64_T, MPI_MIN, communicator_) !=
          MPI_SUCCESS ||
      MPI_Allreduce(clock, maximum, 3, MPI_UINT64_T, MPI_MAX, communicator_) !=
          MPI_SUCCESS)
    local = mpi_failure();
  if (!std::equal(minimum, minimum + 3, maximum))
    local = invalid(MigrationDetail::stale_revision);
  if (local)
    for (std::size_t i = 0; i < input.size; ++i) {
      if (!valid(input.data[i])) {
        local = invalid(MigrationDetail::invalid_value);
        break;
      }
      ParcelLocation located;
      if (!location.locate(input.data[i].parcel.position_m, revision,
                           located)) {
        local = invalid(MigrationDetail::location_unavailable);
        break;
      }
      if (located.revision != revision) {
        local = invalid(MigrationDetail::stale_revision);
        break;
      }
      if (located.global_cell != input.data[i].parcel.owner_global_cell ||
          located.owner_rank < 0 ||
          located.owner_rank != owner(located.global_cell)) {
        local = invalid(MigrationDetail::location_mismatch);
        break;
      }
    }
  int failing = -1;
  auto status = consensus(local, failing);
  if (!status) {
    report.status = status;
    report.lowest_failing_rank = failing;
    return status;
  }
  status = prepare(input, report);
  if (!status)
    return status;
  // A receiver's geometry must agree too; an incoming candidate is not yet
  // an accepted parcel and is withdrawn globally if any receiver rejects it.
  for (const auto &value : candidates_) {
    ParcelLocation located;
    if (!location.locate(value.parcel.position_m, revision, located)) {
      local = invalid(MigrationDetail::location_unavailable);
      break;
    }
    if (located.revision != revision) {
      local = invalid(MigrationDetail::stale_revision);
      break;
    }
    if (located.global_cell != value.parcel.owner_global_cell ||
        located.owner_rank != rank_) {
      local = invalid(MigrationDetail::location_mismatch);
      break;
    }
  }
  status = consensus(local, failing);
  if (!status) {
    discard();
    report = {};
    report.status = status;
    report.lowest_failing_rank = failing;
  }
  return status;
}

Status CartesianParcelLocationProvider::locate(
    const Vector3 &position, portable::Revision revision,
    ParcelLocation &candidate) const noexcept {
  candidate = {};
  if (revision != input_.revision)
    return invalid(MigrationDetail::stale_revision);
  MeshPatch patch;
  if (revision.algorithm_version != 1U ||
      !hundun::v04::detail::make_mesh_patch(0, input_.rank_count,
                                            input_.global_cells, patch))
    return invalid(MigrationDetail::invalid_plan);
  const int counts[]{input_.global_cells.x, input_.global_cells.y,
                     input_.global_cells.z};
  const int partitions[]{patch.process_grid.x, patch.process_grid.y,
                         patch.process_grid.z};
  int cell[3]{}, owner_cell[3]{};
  for (std::size_t d = 0; d < 3; ++d) {
    if (!std::isfinite(position[d]) || !std::isfinite(input_.origin_m[d]) ||
        !std::isfinite(input_.cell_width_m[d]) || input_.cell_width_m[d] <= 0.0)
      return invalid(MigrationDetail::invalid_value);
    const double coordinate =
        (position[d] - input_.origin_m[d]) / input_.cell_width_m[d];
    if (!std::isfinite(coordinate) || coordinate < 0 || coordinate >= counts[d])
      return invalid(MigrationDetail::location_unavailable);
    cell[d] = static_cast<int>(std::floor(coordinate));
    owner_cell[d] = owner_coordinate(cell[d], counts[d], partitions[d]);
  }
  ParcelLocation result;
  result.revision = revision;
  result.global_cell = static_cast<std::uint64_t>(cell[0]) +
                       static_cast<std::uint64_t>(counts[0]) *
                           (static_cast<std::uint64_t>(cell[1]) +
                            static_cast<std::uint64_t>(counts[1]) * cell[2]);
  result.owner_rank =
      owner_cell[0] +
      partitions[0] * (owner_cell[1] + partitions[1] * owner_cell[2]);
  candidate = result;
  return {};
}

} // namespace hundun::v04::spray::detail
