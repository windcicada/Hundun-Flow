// SPDX-License-Identifier: Apache-2.0
#include "mesh_focus_detail.hpp"
#include "models_exchange_owner_detail.hpp"
#include <algorithm>
#include <climits>
#include <cmath>
#include <cstring>
#include <limits>
#include <tuple>

namespace hundun::v04::portable {
namespace {
constexpr int lanes = 17;
using NativeStatus = hundun::v04::Status;
bool same(Int3 a, Int3 b) noexcept {
  return a.x == b.x && a.y == b.y && a.z == b.z;
}
std::uint64_t bits(double d) noexcept {
  std::uint64_t v;
  std::memcpy(&v, &d, 8);
  return v;
}
double real(std::uint64_t v) noexcept {
  double d;
  std::memcpy(&d, &v, 8);
  return d;
}
void encode(const ExchangeSegment &s, std::uint64_t *w) noexcept {
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
ExchangeSegment decode(const std::uint64_t *w) noexcept {
  ExchangeSegment s;
  s.revision = {w[0], w[1], static_cast<std::uint32_t>(w[2])};
  s.global_cell = w[3];
  s.parcel_id = {w[4], w[5]};
  s.segment_ordinal = w[6];
  s.deposition_weight = real(w[7]);
  s.channel = static_cast<ExchangeChannel>(w[8]);
  s.delta.mass_delta_kg = real(w[9]);
  for (int d = 0; d < 3; ++d)
    s.delta.momentum_delta_kg_m_per_s[d] = real(w[10 + d]);
  s.delta.thermochemical_enthalpy_delta_j = real(w[13]);
  s.delta.kinetic_energy_delta_j = real(w[14]);
  s.delta.thermal_exchange_to_gas_j = real(w[15]);
  s.vapor_species_index = w[16];
  return s;
}
auto key(const ExchangeSegment &s) noexcept {
  return std::tie(s.parcel_id.high, s.parcel_id.low, s.segment_ordinal);
}
bool same_budget(const ExchangeSegment &a, const ExchangeSegment &b) noexcept {
  return a.channel == b.channel &&
         a.vapor_species_index == b.vapor_species_index &&
         a.delta.mass_delta_kg == b.delta.mass_delta_kg &&
         a.delta.momentum_delta_kg_m_per_s ==
             b.delta.momentum_delta_kg_m_per_s &&
         a.delta.thermochemical_enthalpy_delta_j ==
             b.delta.thermochemical_enthalpy_delta_j &&
         a.delta.kinetic_energy_delta_j == b.delta.kinetic_energy_delta_j &&
         a.delta.thermal_exchange_to_gas_j == b.delta.thermal_exchange_to_gas_j;
}
int owner_axis(int i, int cells, int parts) noexcept {
  const int base = cells / parts, extra = cells % parts;
  const std::int64_t extended = std::int64_t(base + 1) * extra;
  return i < extended ? i / (base + 1) : extra + int((i - extended) / base);
}
} // namespace
NativeStatus OwnerExchangeRoutingPlan::configure(
    MPI_Comm comm, Int3 global, MeshPatch patch, std::size_t local_capacity,
    std::size_t audit_capacity, std::size_t species,
    std::uint64_t maximum_bytes) noexcept {
  if (comm == MPI_COMM_NULL)
    return {StatusCode::invalid_plan, 26200};
  int rank = 0, ranks = 0;
  if (MPI_Comm_rank(comm, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(comm, &ranks) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, 26201};
  OwnerExchangeRoutingPlan fresh;
  fresh.comm_ = comm;
  fresh.rank_ = rank;
  fresh.ranks_ = ranks;
  MeshPatch actual;
  Status local = Status::success;
  if (!hundun::v04::detail::make_mesh_patch(rank, ranks, global, actual) ||
      !same(actual.begin, patch.begin) || !same(actual.cells, patch.cells) ||
      !same(actual.process_grid, patch.process_grid) ||
      !same(actual.process_coord, patch.process_coord) ||
      local_capacity > std::size_t(INT_MAX / lanes) ||
      audit_capacity > std::size_t(INT_MAX / lanes) || !species ||
      species > 255)
    local = Status::invalid_input;
  if (fresh.agree(local) != Status::success)
    return {StatusCode::invalid_plan, 26200};
  std::uint64_t signature[]{std::uint64_t(global.x), std::uint64_t(global.y),
                            std::uint64_t(global.z), species},
      low[4]{}, high[4]{};
  if (MPI_Allreduce(signature, low, 4, MPI_UINT64_T, MPI_MIN, comm) !=
          MPI_SUCCESS ||
      MPI_Allreduce(signature, high, 4, MPI_UINT64_T, MPI_MAX, comm) !=
          MPI_SUCCESS)
    return {StatusCode::mpi_failure, 26201};
  if (!std::equal(low, low + 4, high))
    return {StatusCode::invalid_plan, 26200};
  const auto capacity = std::max(local_capacity, audit_capacity);
  const auto send_words = std::max<std::size_t>(1, local_capacity * lanes);
  const auto receive_words = std::max<std::size_t>(1, capacity * lanes);
  // Bounds above ensure these products fit uint64 even at maximal MPI int
  // sizes.
  const std::uint64_t bytes = sizeof(fresh) +
                              std::uint64_t(ranks) * 5 * sizeof(int) +
                              std::uint64_t(send_words + receive_words) * 8 +
                              std::uint64_t(capacity) * sizeof(ExchangeSegment);
  if (bytes > maximum_bytes)
    local = Status::capacity_exceeded;
  if (fresh.agree(local) != Status::success)
    return {StatusCode::allocation_failure, 26202};
  try {
    fresh.send_counts_.resize(ranks);
    fresh.receive_counts_.resize(ranks);
    fresh.send_offsets_.resize(ranks);
    fresh.receive_offsets_.resize(ranks);
    fresh.cursor_.resize(ranks);
    fresh.send_.resize(send_words);
    fresh.receive_.resize(receive_words);
    fresh.values_.resize(capacity);
  } catch (...) {
    local = Status::capacity_exceeded;
  }
  if (fresh.agree(local) != Status::success)
    return {StatusCode::allocation_failure, 26203};
  fresh.global_ = global;
  fresh.grid_ = patch.process_grid;
  fresh.local_capacity_ = local_capacity;
  fresh.audit_capacity_ = audit_capacity;
  fresh.species_ = species;
  fresh.owned_bytes_ = bytes;
  // Reconfiguration invalidates old borrowed certificates even when no payload
  // is published.
  fresh.generation_ = generation_ + 1;
  *this = std::move(fresh);
  return {};
}
Status OwnerExchangeRoutingPlan::agree(Status local) const noexcept {
  int first = local == Status::success ? INT_MAX : rank_, failing = 0;
  if (MPI_Allreduce(&first, &failing, 1, MPI_INT, MPI_MIN, comm_) !=
      MPI_SUCCESS)
    return Status::provider_failure;
  if (failing == INT_MAX)
    return Status::success;
  int code = static_cast<int>(local);
  if (MPI_Bcast(&code, 1, MPI_INT, failing, comm_) != MPI_SUCCESS)
    return Status::provider_failure;
  return static_cast<Status>(code);
}
int OwnerExchangeRoutingPlan::owner(std::uint64_t id) const noexcept {
  const int x = int(id % std::uint64_t(global_.x));
  id /= std::uint64_t(global_.x);
  const int y = int(id % std::uint64_t(global_.y));
  id /= std::uint64_t(global_.y);
  if (id >= std::uint64_t(global_.z))
    return -1;
  return owner_axis(x, global_.x, grid_.x) +
         grid_.x * (owner_axis(y, global_.y, grid_.y) +
                    grid_.y * owner_axis(int(id), global_.z, grid_.z));
}
int OwnerExchangeRoutingPlan::auditor(const ExchangeSegment &s) const noexcept {
  // Only bucket selection is hashed; the complete integer tuple is compared.
  std::uint64_t h = UINT64_C(14695981039346656037);
  for (auto v : {s.parcel_id.high, s.parcel_id.low, s.segment_ordinal})
    h = (h ^ v) * UINT64_C(1099511628211);
  return int(h % std::uint64_t(ranks_));
}
Status OwnerExchangeRoutingPlan::exchange(const ExchangeSegment *input,
                                          std::size_t n, bool audit,
                                          std::size_t &received) noexcept {
  std::fill(send_counts_.begin(), send_counts_.end(), 0);
  for (std::size_t i = 0; i < n; ++i)
    ++send_counts_[audit ? auditor(input[i]) : owner(input[i].global_cell)];
  if (MPI_Alltoall(send_counts_.data(), 1, MPI_INT, receive_counts_.data(), 1,
                   MPI_INT, comm_) != MPI_SUCCESS)
    return Status::provider_failure;
  std::size_t sent = 0;
  received = 0;
  const auto capacity = audit ? audit_capacity_ : local_capacity_;
  Status local = Status::success;
  for (int r = 0; r < ranks_; ++r) {
    if (receive_counts_[r] < 0 ||
        std::size_t(receive_counts_[r]) > capacity - received) {
      local = Status::capacity_exceeded;
      break;
    }
    send_offsets_[r] = int(sent * lanes);
    receive_offsets_[r] = int(received * lanes);
    sent += send_counts_[r];
    received += receive_counts_[r];
    send_counts_[r] *= lanes;
    receive_counts_[r] *= lanes;
    cursor_[r] = send_offsets_[r];
  }
  auto status = agree(local);
  if (status != Status::success)
    return status;
  for (std::size_t i = 0; i < n; ++i) {
    const auto target = audit ? auditor(input[i]) : owner(input[i].global_cell);
    encode(input[i], send_.data() + cursor_[target]);
    cursor_[target] += lanes;
  }
  if (MPI_Alltoallv(send_.data(), send_counts_.data(), send_offsets_.data(),
                    MPI_UINT64_T, receive_.data(), receive_counts_.data(),
                    receive_offsets_.data(), MPI_UINT64_T,
                    comm_) != MPI_SUCCESS)
    return Status::provider_failure;
  for (std::size_t i = 0; i < received; ++i)
    values_[i] = decode(receive_.data() + i * lanes);
  return Status::success;
}
Status OwnerExchangeRoutingPlan::prepare(Revision revision,
                                         const ExchangeSegment *input,
                                         std::size_t n) noexcept {
  discard();
  if (comm_ == MPI_COMM_NULL)
    return Status::invalid_input;
  Status local = Status::success;
  if (n > local_capacity_)
    local = Status::capacity_exceeded;
  else if ((n && !input) || revision.algorithm_version != 1 ||
           !revision.input_revision)
    local = Status::invalid_input;
  else
    for (std::size_t i = 0; i < n; ++i) {
      const auto &s = input[i];
      if (s.revision != revision) {
        local = Status::stale_revision;
        break;
      }
      std::uint64_t wire[lanes];
      encode(s, wire);
      bool finite = true;
      for (int j = 9; j <= 15; ++j)
        finite &= std::isfinite(real(wire[j]));
      if (owner(s.global_cell) < 0 || (!s.parcel_id.high && !s.parcel_id.low) ||
          !finite || !std::isfinite(s.deposition_weight) ||
          s.deposition_weight <= 0 || s.deposition_weight > 1 ||
          s.channel > ExchangeChannel::external ||
          s.vapor_species_index >= species_ ||
          (s.channel != ExchangeChannel::interphase &&
           s.deposition_weight != 1)) {
        local = Status::invalid_input;
        break;
      }
    }
  auto status = agree(local);
  if (status != Status::success)
    return status;
  std::uint64_t clock[]{revision.accepted_step, revision.input_revision},
      low[2]{}, high[2]{};
  if (MPI_Allreduce(clock, low, 2, MPI_UINT64_T, MPI_MIN, comm_) !=
          MPI_SUCCESS ||
      MPI_Allreduce(clock, high, 2, MPI_UINT64_T, MPI_MAX, comm_) !=
          MPI_SUCCESS)
    return Status::provider_failure;
  if (!std::equal(low, low + 2, high))
    return Status::stale_revision;
  std::size_t received = 0;
  status = exchange(input, n, true, received);
  if (status != Status::success)
    return status;
  std::sort(values_.begin(), values_.begin() + received,
            [](const auto &a, const auto &b) {
              return key(a) != key(b) ? key(a) < key(b)
                                      : a.global_cell < b.global_cell;
            });
  for (std::size_t i = 0; i < received;) {
    std::size_t end = i + 1;
    while (end < received && key(values_[i]) == key(values_[end]))
      ++end;
    long double weight = 0;
    for (std::size_t j = i; j < end; ++j) {
      if (!same_budget(values_[i], values_[j]) ||
          (j > i && (values_[j - 1].global_cell == values_[j].global_cell ||
                     values_[j].channel != ExchangeChannel::interphase))) {
        local = Status::invalid_input;
        break;
      }
      weight += values_[j].deposition_weight;
    }
    if (local != Status::success)
      break;
    if (std::abs(weight - 1.L) > 32 * std::numeric_limits<double>::epsilon()) {
      local = Status::conservation_failure;
      break;
    }
    i = end;
  }
  status = agree(local);
  if (status != Status::success)
    return status;
  status = exchange(input, n, false, count_);
  if (status != Status::success)
    return status;
  // Stable summation order is independent of sender rank and rank count.
  std::sort(values_.begin(), values_.begin() + count_,
            [](const auto &a, const auto &b) {
              return a.global_cell != b.global_cell
                         ? a.global_cell < b.global_cell
                         : key(a) < key(b);
            });
  for (std::size_t i = 0; i < count_; ++i)
    if (owner(values_[i].global_cell) != rank_)
      local = Status::invalid_input;
  status = agree(local);
  if (status != Status::success)
    return status;
  revision_ = revision;
  available_ = true;
  return Status::success;
}
RoutedExchangeSegments OwnerExchangeRoutingPlan::candidates() const noexcept {
  RoutedExchangeSegments result;
  if (available_) {
    result.data_ = values_.data();
    result.size_ = count_;
    result.revision_ = revision_;
    result.generation_ = &generation_;
    result.captured_generation_ = generation_;
  }
  return result;
}
} // namespace hundun::v04::portable
