// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_ibm.hpp"

#include "field_view_interval_detail.hpp"

#include <mpi.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <new>
#include <utility>
#include <vector>

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kDonorInput = 14501U;
constexpr std::uint32_t kDonorCollective = 14502U;
constexpr std::uint32_t kDonorLayout = 14503U;
constexpr std::uint32_t kDonorExchange = 14504U;
constexpr std::uint64_t kFnvOffset = UINT64_C(1469598103934665603);
constexpr std::uint64_t kFnvPrime = UINT64_C(1099511628211);

std::uint64_t mix(std::uint64_t hash, std::uint64_t value) noexcept {
  hash ^= value;
  hash *= kFnvPrime;
  return hash;
}

bool same(Int3 left, Int3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

bool inside(Int3 cell, Int3 shape) noexcept {
  return cell.x >= 0 && cell.y >= 0 && cell.z >= 0 && cell.x < shape.x &&
         cell.y < shape.y && cell.z < shape.z;
}

bool in_padded(Int3 cell, Int3 shape, Int3 ghosts) noexcept {
  return cell.x >= -ghosts.x && cell.y >= -ghosts.y &&
         cell.z >= -ghosts.z && cell.x < shape.x + ghosts.x &&
         cell.y < shape.y + ghosts.y && cell.z < shape.z + ghosts.z;
}

bool periodic_extent_valid(Int3 global, unsigned reach,
                           const QuadraticStencilPlan& reconstruction) noexcept {
  const std::int32_t extents[3]{global.x, global.y, global.z};
  for (int axis = 0; axis < 3; ++axis) {
    if (reconstruction.periodic_axis(static_cast<CartesianAxis>(axis)) &&
        static_cast<std::int64_t>(extents[axis]) <=
            2 * static_cast<std::int64_t>(reach)) {
      return false;
    }
  }
  return true;
}

// A reconstruction donor carries its canonical global cell id separately
// from its raw logical index.  The latter may be a one-period image in the
// certified padded field; no arbitrary modulo or multi-period aliases are
// accepted here.
bool canonicalize_raw(Int3 raw, Int3 global, unsigned reach,
                      const QuadraticStencilPlan& reconstruction,
                      Int3& canonical) noexcept {
  const std::int32_t raw_values[3]{raw.x, raw.y, raw.z};
  const std::int32_t extents[3]{global.x, global.y, global.z};
  std::int32_t canonical_values[3]{};
  for (int axis = 0; axis < 3; ++axis) {
    const std::int64_t value = raw_values[axis];
    const std::int64_t extent = extents[axis];
    if (value >= 0 && value < extent) {
      canonical_values[axis] = raw_values[axis];
      continue;
    }
    if (!reconstruction.periodic_axis(static_cast<CartesianAxis>(axis)) ||
        extent <= 2 * static_cast<std::int64_t>(reach)) {
      return false;
    }
    if (value < 0) {
      if (value < -static_cast<std::int64_t>(reach)) return false;
      const std::int64_t wrapped = value + extent;
      if (wrapped < 0 || wrapped >= extent) return false;
      canonical_values[axis] = static_cast<std::int32_t>(wrapped);
    } else {
      if (value >= extent + static_cast<std::int64_t>(reach)) return false;
      const std::int64_t wrapped = value - extent;
      if (wrapped < 0 || wrapped >= extent) return false;
      canonical_values[axis] = static_cast<std::int32_t>(wrapped);
    }
  }
  canonical = {canonical_values[0], canonical_values[1], canonical_values[2]};
  return true;
}

bool raw_index(Int3 local, MeshPatch patch, Int3& raw) noexcept {
  const std::int64_t values[3]{
      static_cast<std::int64_t>(local.x) + patch.begin.x,
      static_cast<std::int64_t>(local.y) + patch.begin.y,
      static_cast<std::int64_t>(local.z) + patch.begin.z};
  for (const std::int64_t value : values) {
    if (value < std::numeric_limits<std::int32_t>::min() ||
        value > std::numeric_limits<std::int32_t>::max())
      return false;
  }
  raw = {static_cast<std::int32_t>(values[0]),
         static_cast<std::int32_t>(values[1]),
         static_cast<std::int32_t>(values[2])};
  return true;
}

bool valid_shape(Int3 shape) noexcept {
  return shape.x > 0 && shape.y > 0 && shape.z > 0;
}

bool valid_patch(Int3 global, MeshPatch patch, int rank, int size) noexcept {
  if (!valid_shape(global) || !valid_shape(patch.cells) ||
      !valid_shape(patch.process_grid) || patch.process_coord.x < 0 ||
      patch.process_coord.y < 0 || patch.process_coord.z < 0 ||
      patch.process_coord.x >= patch.process_grid.x ||
      patch.process_coord.y >= patch.process_grid.y ||
      patch.process_coord.z >= patch.process_grid.z)
    return false;
  const std::int64_t ranks =
      static_cast<std::int64_t>(patch.process_grid.x) *
      patch.process_grid.y * patch.process_grid.z;
  const int expected_rank =
      patch.process_coord.x +
      patch.process_grid.x *
          (patch.process_coord.y +
           patch.process_grid.y * patch.process_coord.z);
  return ranks == size && expected_rank == rank && patch.begin.x >= 0 &&
         patch.begin.y >= 0 && patch.begin.z >= 0 &&
         patch.begin.x + patch.cells.x <= global.x &&
         patch.begin.y + patch.cells.y <= global.y &&
         patch.begin.z + patch.cells.z <= global.z;
}

Int3 decode(GlobalCellId id, Int3 global, bool& valid) noexcept {
  const std::uint64_t gx = static_cast<std::uint64_t>(global.x);
  const std::uint64_t gy = static_cast<std::uint64_t>(global.y);
  const std::uint64_t gz = static_cast<std::uint64_t>(global.z);
  const std::uint64_t total = gx * gy * gz;
  valid = id < total;
  if (!valid) return {};
  const std::uint64_t plane = gx * gy;
  const std::uint64_t z = id / plane;
  const std::uint64_t remainder = id - z * plane;
  const std::uint64_t y = remainder / gx;
  const std::uint64_t x = remainder - y * gx;
  return {static_cast<std::int32_t>(x), static_cast<std::int32_t>(y),
          static_cast<std::int32_t>(z)};
}

std::int32_t owner_coordinate(std::int32_t index, std::int32_t global,
                              std::int32_t partitions) noexcept {
  const std::int32_t base = global / partitions;
  const std::int32_t remainder = global % partitions;
  const std::int32_t wide_end = (base + 1) * remainder;
  return index < wide_end ? index / (base + 1)
                          : remainder + (index - wide_end) / base;
}

int owner_rank(Int3 global_index, Int3 global, Int3 grid) noexcept {
  const Int3 owner{owner_coordinate(global_index.x, global.x, grid.x),
                   owner_coordinate(global_index.y, global.y, grid.y),
                   owner_coordinate(global_index.z, global.z, grid.z)};
  return owner.x + grid.x * (owner.y + grid.y * owner.z);
}

struct Need {
  int owner{};
  GlobalCellId global{};
  Int3 local{};
};

bool need_less(const Need& left, const Need& right) noexcept {
  if (left.owner != right.owner) return left.owner < right.owner;
  if (left.global != right.global) return left.global < right.global;
  if (left.local.x != right.local.x) return left.local.x < right.local.x;
  if (left.local.y != right.local.y) return left.local.y < right.local.y;
  return left.local.z < right.local.z;
}

bool multiply(std::size_t left, std::size_t right,
              std::size_t& out) noexcept {
  if (left != 0U && right > std::numeric_limits<std::size_t>::max() / left)
    return false;
  out = left * right;
  return true;
}

bool add(std::uint64_t left, std::uint64_t right,
         std::uint64_t& out) noexcept {
  if (right > std::numeric_limits<std::uint64_t>::max() - left) return false;
  out = left + right;
  return true;
}

bool mpi_live() noexcept {
  int initialized = 0;
  int finalized = 0;
  return MPI_Initialized(&initialized) == MPI_SUCCESS && initialized != 0 &&
         MPI_Finalized(&finalized) == MPI_SUCCESS && finalized == 0;
}

}  // namespace

struct RemoteDonorExchangePlan::Impl {
  MPI_Comm communicator{MPI_COMM_NULL};
  MeshPatch patch{};
  Int3 global{};
  std::vector<RemoteDonorFieldSpec> fields;
  std::vector<Int3> receive_targets;
  std::vector<Int3> send_sources;
  std::vector<double> receive_buffer;
  std::vector<double> send_buffer;
  std::vector<MPI_Request> requests;
  std::vector<int> receive_counts;
  std::vector<int> receive_displacements;
  std::vector<int> send_counts;
  std::vector<int> send_displacements;
  std::size_t values_per_cell{};
  std::uint8_t reach{};
  StageId stage{};
  RemoteDonorExchangeStats stats{};
  RemoteDonorExchangeCounters counters{};
  PlanFingerprint fingerprint{};
  bool ready{};

  ~Impl() noexcept {
    if (!mpi_live()) return;
    for (MPI_Request& request : requests) {
      if (request != MPI_REQUEST_NULL) (void)MPI_Request_free(&request);
    }
    if (communicator != MPI_COMM_NULL) (void)MPI_Comm_free(&communicator);
  }
};

RemoteDonorExchangePlan::~RemoteDonorExchangePlan() noexcept { release(); }

RemoteDonorExchangePlan::RemoteDonorExchangePlan(
    RemoteDonorExchangePlan&& other) noexcept
    : implementation_(std::exchange(other.implementation_, nullptr)) {}

RemoteDonorExchangePlan& RemoteDonorExchangePlan::operator=(
    RemoteDonorExchangePlan&& other) noexcept {
  if (this != &other) {
    release();
    implementation_ = std::exchange(other.implementation_, nullptr);
  }
  return *this;
}

void RemoteDonorExchangePlan::release() noexcept {
  Impl* implementation = std::exchange(implementation_, nullptr);
  if (implementation == nullptr) return;
  delete implementation;
}

Status RemoteDonorExchangePlan::analyze(
    MPI_Comm communicator, Int3 global_cells, MeshPatch patch,
    const QuadraticStencilPlan& reconstruction,
    Span<const RemoteDonorFieldSpec> fields, StageId stage,
    RemoteDonorExchangePlan& out) noexcept try {
  if (communicator == MPI_COMM_NULL || out.implementation_ != nullptr ||
      reconstruction.fingerprint() == 0U ||
      fields.data == nullptr || fields.size == 0U || stage == 0U)
    return {StatusCode::invalid_plan, kDonorInput};
  int rank = -1;
  int size = 0;
  if (MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(communicator, &size) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kDonorCollective};
  const unsigned local_reach = reconstruction.maximum_halo_reach();
  unsigned global_reach = 0U;
  if (MPI_Allreduce(&local_reach, &global_reach, 1, MPI_UNSIGNED, MPI_MAX,
                    communicator) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kDonorCollective};
  if (global_reach == 0U || global_reach > UINT8_MAX)
    return {StatusCode::invalid_plan, kDonorInput};
  bool local_valid = valid_patch(global_cells, patch, rank, size);
  std::size_t values_per_cell = 0U;
  for (std::size_t index = 0U; index < fields.size; ++index) {
    const RemoteDonorFieldSpec field = fields.data[index];
    local_valid = local_valid && field.field != 0U && field.components != 0U &&
                  values_per_cell <=
                      std::numeric_limits<std::size_t>::max() -
                          field.components;
    for (std::size_t prior = 0U; prior < index; ++prior)
      local_valid = local_valid && fields.data[prior].field != field.field;
    values_per_cell += field.components;
  }
  const auto globals = reconstruction.donor_global_cells();
  const auto locals = reconstruction.donor_local_indices();
  local_valid = local_valid && globals.size == locals.size;
  int valid_integer = local_valid ? 1 : 0;
  int globally_valid = 0;
  if (MPI_Allreduce(&valid_integer, &globally_valid, 1, MPI_INT, MPI_MIN,
                    communicator) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kDonorCollective};
  if (globally_valid == 0)
    return {StatusCode::invalid_plan, kDonorInput};
  local_valid = periodic_extent_valid(global_cells, global_reach,
                                      reconstruction);
  valid_integer = local_valid ? 1 : 0;
  if (MPI_Allreduce(&valid_integer, &globally_valid, 1, MPI_INT, MPI_MIN,
                    communicator) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kDonorCollective};
  if (globally_valid == 0)
    return {StatusCode::invalid_plan, kDonorInput};

  std::vector<Need> needs;
  needs.reserve(globals.size);
  for (std::size_t index = 0U; index < globals.size; ++index) {
    bool decoded = false;
    const Int3 canonical = decode(globals.data[index], global_cells, decoded);
    Int3 raw{};
    const Int3 local = locals.data[index];
    Int3 canonical_from_raw{};
    if (!decoded || !raw_index(local, patch, raw) ||
        !in_padded(local, patch.cells,
                   {static_cast<std::int32_t>(global_reach),
                    static_cast<std::int32_t>(global_reach),
                    static_cast<std::int32_t>(global_reach)}) ||
        !canonicalize_raw(raw, global_cells, global_reach, reconstruction,
                          canonical_from_raw) ||
        !same(canonical, canonical_from_raw)) {
      local_valid = false;
      break;
    }
    if (!inside(local, patch.cells))
      needs.push_back(
          {owner_rank(canonical, global_cells, patch.process_grid),
           globals.data[index], local});
  }
  std::sort(needs.begin(), needs.end(), need_less);
  // A canonical donor may serve distinct raw ghost targets.  The same
  // canonical/target pair can also be referenced by several reconstruction
  // groups, so retain one request for that storage slot deterministically.
  // Per-group canonical duplicates are rejected by collect_donors before this
  // exchange plan is constructed.
  const auto unique_end = std::unique(
      needs.begin(), needs.end(), [](const Need& left, const Need& right) {
        return left.global == right.global && same(left.local, right.local);
      });
  needs.erase(unique_end, needs.end());
  for (const Need& need : needs)
    local_valid = local_valid && need.owner >= 0 && need.owner < size;
  valid_integer = local_valid ? 1 : 0;
  if (MPI_Allreduce(&valid_integer, &globally_valid, 1, MPI_INT, MPI_MIN,
                    communicator) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kDonorCollective};
  if (globally_valid == 0)
    return {StatusCode::invalid_plan, kDonorLayout};

  std::vector<int> request_counts(static_cast<std::size_t>(size), 0);
  for (const Need& need : needs) {
    int& count = request_counts[static_cast<std::size_t>(need.owner)];
    if (count == std::numeric_limits<int>::max())
      return {StatusCode::invalid_plan, kDonorLayout};
    ++count;
  }
  std::vector<int> supply_counts(static_cast<std::size_t>(size), 0);
  if (MPI_Alltoall(request_counts.data(), 1, MPI_INT, supply_counts.data(), 1,
                   MPI_INT, communicator) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kDonorCollective};
  std::vector<int> request_displacements(static_cast<std::size_t>(size), 0);
  std::vector<int> supply_displacements(static_cast<std::size_t>(size), 0);
  int request_total = 0;
  int supply_total = 0;
  for (int peer = 0; peer < size; ++peer) {
    request_displacements[static_cast<std::size_t>(peer)] = request_total;
    supply_displacements[static_cast<std::size_t>(peer)] = supply_total;
    if (request_counts[static_cast<std::size_t>(peer)] >
            std::numeric_limits<int>::max() - request_total ||
        supply_counts[static_cast<std::size_t>(peer)] >
            std::numeric_limits<int>::max() - supply_total)
      return {StatusCode::invalid_plan, kDonorLayout};
    request_total += request_counts[static_cast<std::size_t>(peer)];
    supply_total += supply_counts[static_cast<std::size_t>(peer)];
  }
  std::vector<GlobalCellId> requested(static_cast<std::size_t>(request_total));
  for (std::size_t index = 0U; index < needs.size(); ++index)
    requested[index] = needs[index].global;
  std::vector<GlobalCellId> supplied(static_cast<std::size_t>(supply_total));
  if (MPI_Alltoallv(requested.data(), request_counts.data(),
                    request_displacements.data(), MPI_UINT64_T,
                    supplied.data(), supply_counts.data(),
                    supply_displacements.data(), MPI_UINT64_T,
                    communicator) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kDonorCollective};

  std::vector<Int3> send_sources;
  send_sources.reserve(supplied.size());
  for (GlobalCellId id : supplied) {
    bool decoded = false;
    const Int3 global = decode(id, global_cells, decoded);
    const Int3 local{global.x - patch.begin.x, global.y - patch.begin.y,
                     global.z - patch.begin.z};
    if (!decoded || !inside(local, patch.cells) ||
        owner_rank(global, global_cells, patch.process_grid) != rank) {
      local_valid = false;
      break;
    }
    send_sources.push_back(local);
  }
  valid_integer = local_valid ? 1 : 0;
  if (MPI_Allreduce(&valid_integer, &globally_valid, 1, MPI_INT, MPI_MIN,
                    communicator) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kDonorCollective};
  if (globally_valid == 0)
    return {StatusCode::invalid_plan, kDonorLayout};

  std::unique_ptr<Impl> candidate{new (std::nothrow) Impl};
  if (candidate == nullptr)
    return {StatusCode::allocation_failure, kDonorLayout};
  candidate->patch = patch;
  candidate->global = global_cells;
  candidate->fields.assign(fields.data, fields.data + fields.size);
  candidate->receive_targets.reserve(needs.size());
  for (const Need& need : needs)
    candidate->receive_targets.push_back(need.local);
  candidate->send_sources = std::move(send_sources);
  candidate->values_per_cell = values_per_cell;
  candidate->reach = static_cast<std::uint8_t>(global_reach);
  candidate->stage = stage;
  candidate->receive_counts = request_counts;
  candidate->receive_displacements = request_displacements;
  candidate->send_counts = supply_counts;
  candidate->send_displacements = supply_displacements;
  std::size_t receive_values = 0U;
  std::size_t send_values = 0U;
  if (!multiply(candidate->receive_targets.size(), values_per_cell,
                receive_values) ||
      !multiply(candidate->send_sources.size(), values_per_cell,
                send_values)) {
    return {StatusCode::invalid_plan, kDonorLayout};
  }
  std::uint32_t peer_messages = 0U;
  for (int peer = 0; peer < size; ++peer) {
    const int receive_cells = request_counts[static_cast<std::size_t>(peer)];
    const int send_cells = supply_counts[static_cast<std::size_t>(peer)];
    if ((receive_cells > 0 &&
         values_per_cell >
             static_cast<std::size_t>(std::numeric_limits<int>::max() /
                                      receive_cells)) ||
        (send_cells > 0 &&
         values_per_cell >
             static_cast<std::size_t>(std::numeric_limits<int>::max() /
                                      send_cells))) {
      return {StatusCode::invalid_plan, kDonorLayout};
    }
    peer_messages += receive_cells > 0 ? 1U : 0U;
    peer_messages += send_cells > 0 ? 1U : 0U;
  }
  candidate->stats.received_cells = candidate->receive_targets.size();
  candidate->stats.supplied_cells = candidate->send_sources.size();
  candidate->stats.bytes_per_exchange =
      (receive_values + send_values) * sizeof(double);
  candidate->stats.peer_messages = peer_messages;
  std::uint64_t hash = kFnvOffset;
  hash = mix(hash, reconstruction.fingerprint());
  hash = mix(hash, stage);
  hash = mix(hash, values_per_cell);
  hash = mix(hash, candidate->reach);
  hash = mix(hash, candidate->stats.received_cells);
  hash = mix(hash, candidate->stats.supplied_cells);
  for (const RemoteDonorFieldSpec field : candidate->fields) {
    hash = mix(hash, field.field);
    hash = mix(hash, field.components);
  }
  for (const Need& need : needs) hash = mix(hash, need.global);
  candidate->fingerprint = hash == 0U ? 1U : hash;
  out.implementation_ = candidate.release();
  return {};
} catch (const std::bad_alloc&) {
  return {StatusCode::allocation_failure, kDonorLayout};
} catch (...) {
  return {StatusCode::invalid_plan, kDonorInput};
}

Status RemoteDonorExchangePlan::bind(MPI_Comm communicator) noexcept {
  if (implementation_ == nullptr || communicator == MPI_COMM_NULL ||
      implementation_->fingerprint == 0U || implementation_->ready ||
      implementation_->communicator != MPI_COMM_NULL)
    return {StatusCode::invalid_plan, kDonorInput};
  Impl& plan = *implementation_;
  bool allocated = true;
  try {
    plan.receive_buffer.resize(plan.stats.received_cells *
                               plan.values_per_cell);
    plan.send_buffer.resize(plan.stats.supplied_cells * plan.values_per_cell);
    plan.requests.reserve(plan.stats.peer_messages);
  } catch (const std::bad_alloc&) {
    allocated = false;
  } catch (...) {
    allocated = false;
  }
  const int local_allocated = allocated ? 1 : 0;
  int all_allocated = 0;
  if (MPI_Allreduce(&local_allocated, &all_allocated, 1, MPI_INT, MPI_MIN,
                    communicator) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kDonorCollective};
  if (all_allocated == 0) {
    plan.receive_buffer.clear();
    plan.send_buffer.clear();
    plan.requests.clear();
    return {StatusCode::allocation_failure, kDonorLayout};
  }
  if (MPI_Comm_dup(communicator, &plan.communicator) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kDonorCollective};
  const int tag = 23000 + static_cast<int>(plan.stage % 1000U);
  bool local_ready = true;
  const int size = static_cast<int>(plan.receive_counts.size());
  for (int peer = 0; peer < size && local_ready; ++peer) {
    const int receive_cells = plan.receive_counts[peer];
    const int send_cells = plan.send_counts[peer];
    if (receive_cells > 0) {
      MPI_Request request = MPI_REQUEST_NULL;
      const std::size_t begin =
          static_cast<std::size_t>(plan.receive_displacements[peer]) *
          plan.values_per_cell;
      local_ready = MPI_Recv_init(
                        plan.receive_buffer.data() + begin,
                        static_cast<int>(receive_cells * plan.values_per_cell),
                        MPI_DOUBLE, peer, tag, plan.communicator, &request) ==
                    MPI_SUCCESS;
      if (request != MPI_REQUEST_NULL) plan.requests.push_back(request);
    }
    if (send_cells > 0 && local_ready) {
      MPI_Request request = MPI_REQUEST_NULL;
      const std::size_t begin =
          static_cast<std::size_t>(plan.send_displacements[peer]) *
          plan.values_per_cell;
      local_ready = MPI_Send_init(
                        plan.send_buffer.data() + begin,
                        static_cast<int>(send_cells * plan.values_per_cell),
                        MPI_DOUBLE, peer, tag, plan.communicator, &request) ==
                    MPI_SUCCESS;
      if (request != MPI_REQUEST_NULL) plan.requests.push_back(request);
    }
  }
  const int local_bound = local_ready ? 1 : 0;
  int all_bound = 0;
  if (MPI_Allreduce(&local_bound, &all_bound, 1, MPI_INT, MPI_MIN,
                    communicator) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kDonorCollective};
  if (all_bound == 0 || plan.requests.size() != plan.stats.peer_messages)
    return {StatusCode::mpi_failure, kDonorCollective};
  plan.ready = true;
  return {};
}

Status RemoteDonorExchangePlan::compile(
    MPI_Comm communicator, Int3 global_cells, MeshPatch patch,
    const QuadraticStencilPlan& reconstruction,
    Span<const RemoteDonorFieldSpec> fields, StageId stage,
    RemoteDonorExchangePlan& out) noexcept {
  if (out.implementation_ != nullptr)
    return {StatusCode::invalid_plan, kDonorInput};
  RemoteDonorExchangePlan candidate;
  Status status = analyze(communicator, global_cells, patch, reconstruction,
                          fields, stage, candidate);
  if (status) status = candidate.bind(communicator);
  if (!status) return status;
  out = std::move(candidate);
  return {};
}

Status RemoteDonorExchangePlan::preflight_exchange(
    StageId stage, Span<const FieldView> fields) const noexcept {
  if (implementation_ == nullptr || stage != implementation_->stage ||
      !implementation_->ready ||
      fields.data == nullptr || fields.size != implementation_->fields.size())
    return {StatusCode::invalid_plan, kDonorExchange};
  const Impl& plan = *implementation_;
  for (std::size_t index = 0U; index < fields.size; ++index) {
    const FieldView field = fields.data[index];
    detail::FieldStorageInterval interval{};
    if (field.field != plan.fields[index].field ||
        field.components != plan.fields[index].components ||
        !same(field.interior, plan.patch.cells) ||
        field.ghosts.x < plan.reach || field.ghosts.y < plan.reach ||
        field.ghosts.z < plan.reach || field.revision == 0U ||
        !detail::field_storage_interval(field, interval))
      return {StatusCode::invalid_plan, kDonorExchange};
  }
  RemoteDonorExchangeCounters next = plan.counters;
  if (!add(next.exchange_calls, 1U, next.exchange_calls) ||
      !add(next.peer_messages, plan.stats.peer_messages,
           next.peer_messages) ||
      !add(next.bytes, plan.stats.bytes_per_exchange, next.bytes))
    return {StatusCode::invalid_plan, kDonorExchange};
  return {};
}

Status RemoteDonorExchangePlan::exchange(StageId stage,
                                         Span<FieldView> fields) noexcept {
  if (implementation_ == nullptr ||
      implementation_->communicator == MPI_COMM_NULL ||
      !implementation_->ready)
    return {StatusCode::invalid_plan, kDonorExchange};
  Impl& plan = *implementation_;
  const Status preflight = preflight_exchange(
      stage, {fields.data, fields.size});
  const int local_preflight = preflight ? 1 : 0;
  int global_preflight = 0;
  if (MPI_Allreduce(&local_preflight, &global_preflight, 1, MPI_INT, MPI_MIN,
                    plan.communicator) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kDonorCollective};
  if (global_preflight == 0)
    return {StatusCode::invalid_plan, kDonorExchange};
  RemoteDonorExchangeCounters next = plan.counters;
  if (!add(next.exchange_calls, 1U, next.exchange_calls) ||
      !add(next.peer_messages, plan.stats.peer_messages,
           next.peer_messages) ||
      !add(next.bytes, plan.stats.bytes_per_exchange, next.bytes))
    return {StatusCode::invalid_plan, kDonorExchange};
  std::size_t packed = 0U;
  for (const Int3 source : plan.send_sources) {
    for (std::size_t field_index = 0U; field_index < fields.size;
         ++field_index) {
      const FieldView field = fields.data[field_index];
      for (std::uint8_t component = 0U; component < field.components;
           ++component)
        plan.send_buffer[packed++] = field.unchecked(source, component);
    }
  }
  if (!plan.requests.empty() &&
      (MPI_Startall(static_cast<int>(plan.requests.size()),
                    plan.requests.data()) != MPI_SUCCESS ||
       MPI_Waitall(static_cast<int>(plan.requests.size()),
                   plan.requests.data(), MPI_STATUSES_IGNORE) != MPI_SUCCESS))
    return {StatusCode::mpi_failure, kDonorExchange};
  std::size_t unpacked = 0U;
  for (const Int3 target : plan.receive_targets) {
    for (std::size_t field_index = 0U; field_index < fields.size;
         ++field_index) {
      FieldView& field = fields.data[field_index];
      for (std::uint8_t component = 0U; component < field.components;
           ++component)
        field.unchecked(target, component) = plan.receive_buffer[unpacked++];
    }
  }
  plan.counters = next;
  return {};
}

std::uint8_t RemoteDonorExchangePlan::reach() const noexcept {
  return implementation_ == nullptr ? 0U : implementation_->reach;
}

RemoteDonorExchangeStats RemoteDonorExchangePlan::stats() const noexcept {
  return implementation_ == nullptr ? RemoteDonorExchangeStats{}
                                    : implementation_->stats;
}

RemoteDonorExchangeCounters
RemoteDonorExchangePlan::runtime_counters() const noexcept {
  return implementation_ == nullptr ? RemoteDonorExchangeCounters{}
                                    : implementation_->counters;
}

PlanFingerprint RemoteDonorExchangePlan::fingerprint() const noexcept {
  return implementation_ == nullptr ? 0U : implementation_->fingerprint;
}

bool RemoteDonorExchangePlan::ready() const noexcept {
  return implementation_ != nullptr && implementation_->fingerprint != 0U &&
         implementation_->ready;
}

}  // namespace hundun::v04
