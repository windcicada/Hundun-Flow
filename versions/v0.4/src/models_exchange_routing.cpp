// SPDX-License-Identifier: Apache-2.0
#include "models_exchange_routing_detail.hpp"
#include <algorithm>
#include <climits>
#include <cstring>
#include <stdexcept>
namespace hundun::v04::portable {
namespace {
constexpr std::size_t lanes = 17;
std::uint64_t bits(double x) noexcept {
  std::uint64_t u;
  std::memcpy(&u, &x, 8);
  return u;
}
double value(std::uint64_t u) noexcept {
  double x;
  std::memcpy(&x, &u, 8);
  return x;
}
bool agree(MPI_Comm comm, bool local) noexcept {
  int a = local, b = 0;
  return MPI_Allreduce(&a, &b, 1, MPI_INT, MPI_MIN, comm) == MPI_SUCCESS && b;
}
} // namespace
ExchangeRoutingWorkspace::ExchangeRoutingWorkspace(std::size_t nc,
                                                   std::size_t ns,
                                                   std::size_t nr)
    : cell_capacity_(nc), segment_capacity_(ns), rank_capacity_(nr) {
  if (!nr || nr > INT_MAX || nc > std::size_t(INT_MAX) / nr ||
      ns > std::size_t(INT_MAX) / lanes)
    throw std::invalid_argument("invalid exchange routing capacity");
  counts_.resize(nr);
  offsets_.resize(nr);
  local_ids_.resize(nc);
  global_ids_.resize(nc * nr);
  send_.resize(ns * lanes);
  receive_.resize(ns * lanes);
  candidates_.resize(ns);
}
ExchangeRoutingReport ExchangeRoutingWorkspace::prepare(
    MPI_Comm comm, const ExchangeCell *cells, std::size_t nc,
    const ExchangeSegment *segments, std::size_t ns) noexcept {
  ExchangeRoutingReport out;
  if (comm == MPI_COMM_NULL)
    return out;
  int nr = 0;
  if (MPI_Comm_size(comm, &nr) != MPI_SUCCESS)
    return out;
  if (!agree(comm, nr > 0 && std::size_t(nr) <= rank_capacity_ &&
                       nc <= cell_capacity_ && ns <= segment_capacity_ &&
                       (!nc || cells) && (!ns || segments))) {
    out.status = Status::capacity_exceeded;
    return out;
  }
  auto gather = [&](int local, std::size_t capacity,
                    std::size_t &total) noexcept {
    if (MPI_Allgather(&local, 1, MPI_INT, counts_.data(), 1, MPI_INT, comm) !=
        MPI_SUCCESS)
      return false;
    total = 0;
    bool valid = true;
    for (int r = 0; r < nr; ++r) {
      if (counts_[r] < 0 || std::size_t(counts_[r]) > capacity - total) {
        valid = false;
        break;
      }
      offsets_[r] = static_cast<int>(total);
      total += counts_[r];
    }
    return agree(comm, valid);
  };
  std::size_t total_cells = 0, total_lanes = 0;
  if (!gather(static_cast<int>(nc), global_ids_.size(), total_cells)) {
    out.status = Status::capacity_exceeded;
    return out;
  }
  for (std::size_t i = 0; i < nc; ++i)
    local_ids_[i] = cells[i].global_cell;
  if (MPI_Allgatherv(local_ids_.data(), static_cast<int>(nc), MPI_UINT64_T,
                     global_ids_.data(), counts_.data(), offsets_.data(),
                     MPI_UINT64_T, comm) != MPI_SUCCESS)
    return out;
  std::sort(global_ids_.begin(), global_ids_.begin() + total_cells);
  if (std::adjacent_find(global_ids_.begin(),
                         global_ids_.begin() + total_cells) !=
      global_ids_.begin() + total_cells)
    return out;
  if (!gather(static_cast<int>(ns * lanes), receive_.size(), total_lanes)) {
    out.status = Status::capacity_exceeded;
    return out;
  }
  for (std::size_t i = 0; i < ns; ++i) {
    const auto &s = segments[i];
    auto *w = send_.data() + i * lanes;
    w[0] = s.revision.accepted_step;
    w[1] = s.revision.input_revision;
    w[2] = s.revision.algorithm_version;
    w[3] = s.global_cell;
    w[4] = s.parcel_id.high;
    w[5] = s.parcel_id.low;
    w[6] = s.segment_ordinal;
    w[7] = bits(s.deposition_weight);
    w[8] = static_cast<std::uint64_t>(s.channel);
    w[9] = bits(s.delta.mass_delta_kg);
    for (int d = 0; d < 3; ++d)
      w[10 + d] = bits(s.delta.momentum_delta_kg_m_per_s[d]);
    w[13] = bits(s.delta.thermochemical_enthalpy_delta_j);
    w[14] = bits(s.delta.kinetic_energy_delta_j);
    w[15] = bits(s.delta.thermal_exchange_to_gas_j);
    w[16] = s.vapor_species_index;
  }
  if (MPI_Allgatherv(send_.data(), static_cast<int>(ns * lanes), MPI_UINT64_T,
                     receive_.data(), counts_.data(), offsets_.data(),
                     MPI_UINT64_T, comm) != MPI_SUCCESS)
    return out;
  std::size_t count = 0;
  for (std::size_t i = 0; i < total_lanes / lanes; ++i) {
    const auto *w = receive_.data() + i * lanes;
    if (!std::binary_search(global_ids_.begin(),
                            global_ids_.begin() + total_cells, w[3]) ||
        w[2] > UINT32_MAX ||
        w[8] > static_cast<std::uint64_t>(ExchangeChannel::external))
      return out;
    if (std::find(local_ids_.begin(), local_ids_.begin() + nc, w[3]) ==
        local_ids_.begin() + nc)
      continue;
    ExchangeSegment s;
    s.revision = {w[0], w[1], static_cast<std::uint32_t>(w[2])};
    s.global_cell = w[3];
    s.parcel_id = {w[4], w[5]};
    s.segment_ordinal = w[6];
    s.deposition_weight = value(w[7]);
    s.channel = static_cast<ExchangeChannel>(w[8]);
    s.delta.mass_delta_kg = value(w[9]);
    for (int d = 0; d < 3; ++d)
      s.delta.momentum_delta_kg_m_per_s[d] = value(w[10 + d]);
    s.delta.thermochemical_enthalpy_delta_j = value(w[13]);
    s.delta.kinetic_energy_delta_j = value(w[14]);
    s.delta.thermal_exchange_to_gas_j = value(w[15]);
    s.vapor_species_index = w[16];
    candidates_[count++] = s;
  }
  out.status = Status::success;
  out.available = true;
  out.segments = candidates_.data();
  out.count = count;
  out.global_count = total_lanes / lanes;
  return out;
}
Status ExchangeRoutingWorkspace::audit_ids(MPI_Comm comm,
                                           const spray::ParcelId *ids,
                                           std::size_t count) noexcept {
  if (comm == MPI_COMM_NULL)
    return Status::invalid_input;
  int nr = 0;
  MPI_Comm_size(comm, &nr);
  if (!agree(comm, nr > 0 && std::size_t(nr) <= rank_capacity_ &&
                       count <= segment_capacity_ && (!count || ids)))
    return Status::capacity_exceeded;
  int local = static_cast<int>(count * 2);
  if (MPI_Allgather(&local, 1, MPI_INT, counts_.data(), 1, MPI_INT, comm) !=
      MPI_SUCCESS)
    return Status::provider_failure;
  std::size_t total = 0;
  bool valid = true;
  for (int r = 0; r < nr; ++r) {
    if (counts_[r] < 0 ||
        std::size_t(counts_[r]) > 2 * segment_capacity_ - total) {
      valid = false;
      break;
    }
    offsets_[r] = static_cast<int>(total);
    total += counts_[r];
  }
  if (!agree(comm, valid))
    return Status::capacity_exceeded;
  for (std::size_t i = 0; i < count; ++i) {
    send_[2 * i] = ids[i].high;
    send_[2 * i + 1] = ids[i].low;
  }
  if (MPI_Allgatherv(send_.data(), local, MPI_UINT64_T, receive_.data(),
                     counts_.data(), offsets_.data(), MPI_UINT64_T,
                     comm) != MPI_SUCCESS)
    return Status::provider_failure;
  for (std::size_t i = 0; i < total / 2; ++i) {
    candidates_[i].parcel_id = {receive_[2 * i], receive_[2 * i + 1]};
    if (!receive_[2 * i] && !receive_[2 * i + 1])
      return Status::invalid_input;
  }
  std::sort(
      candidates_.begin(), candidates_.begin() + total / 2,
      [](const auto &a, const auto &b) { return a.parcel_id < b.parcel_id; });
  for (std::size_t i = 1; i < total / 2; ++i)
    if (candidates_[i - 1].parcel_id == candidates_[i].parcel_id)
      return Status::invalid_input;
  return Status::success;
}
} // namespace hundun::v04::portable
