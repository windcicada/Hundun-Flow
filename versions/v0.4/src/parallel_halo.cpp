// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_parallel.hpp"

#include "parallel_halo_detail.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <new>
#include <utility>
#include <vector>

namespace hundun::v04 {
namespace {

constexpr std::uint64_t kMaximumChunkDoubles =
    static_cast<std::uint64_t>(std::numeric_limits<int>::max());

enum class Face : std::uint8_t { xm, xp, ym, yp, zm, zp };

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
using FailurePoint = detail::HaloFailurePoint;
#else
enum class FailurePoint : std::uint8_t {
  none,
  pack,
  start,
  completion,
  unpack,
  reserve_before_contract,
  reserve_before_alltoall
};
#endif

struct FieldSlot {
  HaloFieldSpec spec{};
  FieldView begin_view{};
  FieldView finish_view{};
  RevisionToken pending_revision{};
  RevisionToken ghost_revision{};
  std::uint64_t seen_generation{};
};

struct Section {
  std::size_t field_slot{};
  Face face{};
  Face canonical_destination{};
  std::uint64_t offset{};
  std::uint64_t count{};
};

struct Peer {
  int rank{-1};
  bool local{};
  std::uint64_t send_offset{};
  std::uint64_t receive_offset{};
  std::uint64_t count{};
  std::vector<Section> sends;
  std::vector<Section> receives;
};

struct Chunk {
  std::size_t peer{};
  std::uint64_t offset{};
  int count{};
  int tag{};
};

struct HaloImplData {
  MPI_Comm control{MPI_COMM_NULL};
  MPI_Comm transport{MPI_COMM_NULL};
  MeshPatch patch{};
  HaloTopology topology{};
  std::vector<FieldSlot> fields;
  std::vector<std::int32_t> field_lookup;
  std::vector<FieldView> candidate_begin_views;
  std::vector<Peer> peers;
  std::vector<Chunk> chunks;
  std::vector<double> send_storage;
  std::vector<double> receive_storage;
  std::vector<MPI_Request> requests;
  std::vector<MPI_Status> statuses;
  HaloPlanStats stats{};
  HaloRuntimeCounters counters{};
  std::uint64_t identity{};
  std::uint64_t generation{};
  std::uint64_t field_visit_generation{};
  std::uint64_t maximum_chunk_doubles{kMaximumChunkDoubles};
  int maximum_tag{};
  int rank{-1};
  int size{};
  int lowest_failing_rank{-1};
  bool exchange_active{};
  bool requests_active{};
  bool poisoned{};
  bool prepared_epoch_active{};
  bool prepared_exchange_views_valid{};
};

std::atomic<std::uint64_t> g_next_halo_identity{1U};

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
std::atomic<detail::HaloFailurePoint> g_failure_point{
    detail::HaloFailurePoint::none};
std::atomic<int> g_failure_rank{-1};
std::atomic<std::uint64_t> g_test_maximum_chunk_doubles{
    kMaximumChunkDoubles};
std::atomic<int> g_test_tag_upper_bound{std::numeric_limits<int>::min()};
#endif

bool checked_add(std::uint64_t left, std::uint64_t right,
                 std::uint64_t& out) noexcept {
  if (left > std::numeric_limits<std::uint64_t>::max() - right) {
    return false;
  }
  out = left + right;
  return true;
}

bool checked_multiply(std::uint64_t left, std::uint64_t right,
                      std::uint64_t& out) noexcept {
  if (left != 0U && right > std::numeric_limits<std::uint64_t>::max() / left) {
    return false;
  }
  out = left * right;
  return true;
}

bool positive(Int3 value) noexcept {
  return value.x > 0 && value.y > 0 && value.z > 0;
}

int face_axis(Face face) noexcept {
  return (face == Face::xm || face == Face::xp)
             ? 0
             : ((face == Face::ym || face == Face::yp) ? 1 : 2);
}

Face opposite(Face face) noexcept {
  switch (face) {
    case Face::xm:
      return Face::xp;
    case Face::xp:
      return Face::xm;
    case Face::ym:
      return Face::yp;
    case Face::yp:
      return Face::ym;
    case Face::zm:
      return Face::zp;
    case Face::zp:
      return Face::zm;
  }
  return Face::xm;
}

bool periodic(HaloTopology topology, int axis) noexcept {
  return axis == 0 ? topology.periodic_x
                   : (axis == 1 ? topology.periodic_y : topology.periodic_z);
}

std::uint64_t hash_mix(std::uint64_t hash, std::uint64_t value) noexcept {
  hash ^= value;
  hash *= 1099511628211ULL;
  return hash;
}

bool valid_status_code(StatusCode code) noexcept {
  return static_cast<std::uint16_t>(code) <=
         static_cast<std::uint16_t>(StatusCode::io_failure);
}

Status raw_collective_status(MPI_Comm communicator, int rank, int size,
                             Status local,
                             int& lowest_failing_rank) noexcept {
  if (!valid_status_code(local.code)) {
    local = {StatusCode::invalid_plan, detail::halo_detail_state};
  }
  const std::uint64_t candidate =
      local.code == StatusCode::ok
          ? std::numeric_limits<std::uint64_t>::max()
          : static_cast<std::uint64_t>(rank);
  std::uint64_t selected = std::numeric_limits<std::uint64_t>::max();
  if (MPI_Allreduce(&candidate, &selected, 1, MPI_UINT64_T, MPI_MIN,
                    communicator) != MPI_SUCCESS) {
    lowest_failing_rank = -1;
    return {StatusCode::mpi_failure, detail::halo_detail_collective};
  }
  if (selected == std::numeric_limits<std::uint64_t>::max()) {
    lowest_failing_rank = -1;
    return {};
  }
  if (selected >= static_cast<std::uint64_t>(size)) {
    lowest_failing_rank = -1;
    return {StatusCode::mpi_failure, detail::halo_detail_collective};
  }
  const int selected_rank = static_cast<int>(selected);
  std::array<std::uint64_t, 2U> wire{};
  if (rank == selected_rank) {
    wire[0] = static_cast<std::uint64_t>(local.code);
    wire[1] = local.detail;
  }
  if (MPI_Bcast(wire.data(), static_cast<int>(wire.size()), MPI_UINT64_T,
                selected_rank, communicator) != MPI_SUCCESS ||
      wire[0] > static_cast<std::uint64_t>(StatusCode::io_failure) ||
      wire[1] > std::numeric_limits<std::uint32_t>::max()) {
    lowest_failing_rank = -1;
    return {StatusCode::mpi_failure, detail::halo_detail_collective};
  }
  lowest_failing_rank = selected_rank;
  return {static_cast<StatusCode>(wire[0]),
          static_cast<std::uint32_t>(wire[1])};
}

Status collective_status(HaloImplData& implementation, Status local) noexcept {
  ++implementation.counters.control_consensus_calls;
  return raw_collective_status(implementation.control, implementation.rank,
                               implementation.size, local,
                               implementation.lowest_failing_rank);
}

void retain_deferred(Status local, Status& deferred) noexcept {
  if (deferred && !local) {
    deferred = local;
  }
}

bool injected(FailurePoint point, int rank) noexcept {
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  return g_failure_point.load(std::memory_order_relaxed) == point &&
         g_failure_rank.load(std::memory_order_relaxed) == rank;
#else
  (void)point;
  (void)rank;
  return false;
#endif
}

std::uint64_t configured_maximum_chunk_doubles() noexcept {
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  return g_test_maximum_chunk_doubles.load(std::memory_order_relaxed);
#else
  return kMaximumChunkDoubles;
#endif
}

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
bool configured_tag_upper_bound(int& upper_bound) noexcept {
  const int configured =
      g_test_tag_upper_bound.load(std::memory_order_relaxed);
  if (configured == std::numeric_limits<int>::min()) {
    return false;
  }
  upper_bound = configured;
  return true;
}
#endif

void invalidate_certificates(HaloImplData& implementation) noexcept {
  for (FieldSlot& field : implementation.fields) {
    field.ghost_revision = RevisionToken{0U};
  }
}

Status validate_patch(const MeshPatch& patch, int rank, int size) noexcept {
  if (!positive(patch.cells) || !positive(patch.process_grid) ||
      patch.begin.x < 0 || patch.begin.y < 0 || patch.begin.z < 0 ||
      patch.process_coord.x < 0 || patch.process_coord.y < 0 ||
      patch.process_coord.z < 0 ||
      patch.process_coord.x >= patch.process_grid.x ||
      patch.process_coord.y >= patch.process_grid.y ||
      patch.process_coord.z >= patch.process_grid.z) {
    return {StatusCode::invalid_plan, detail::halo_detail_topology};
  }
  std::uint64_t xy = 0U;
  std::uint64_t ranks = 0U;
  if (!checked_multiply(static_cast<std::uint64_t>(patch.process_grid.x),
                        static_cast<std::uint64_t>(patch.process_grid.y), xy) ||
      !checked_multiply(xy,
                        static_cast<std::uint64_t>(patch.process_grid.z),
                        ranks) ||
      ranks != static_cast<std::uint64_t>(size)) {
    return {StatusCode::invalid_plan, detail::halo_detail_topology};
  }
  const std::uint64_t expected_rank =
      static_cast<std::uint64_t>(patch.process_coord.x) +
      static_cast<std::uint64_t>(patch.process_grid.x) *
          (static_cast<std::uint64_t>(patch.process_coord.y) +
           static_cast<std::uint64_t>(patch.process_grid.y) *
               static_cast<std::uint64_t>(patch.process_coord.z));
  if (expected_rank != static_cast<std::uint64_t>(rank)) {
    return {StatusCode::invalid_plan, detail::halo_detail_topology};
  }
  return {};
}

Status validate_specs(Span<const HaloFieldSpec> fields,
                      const MeshPatch& patch) noexcept {
  if (fields.data == nullptr || fields.size == 0U ||
      fields.size > static_cast<std::size_t>(
                        std::numeric_limits<FieldId>::max()) +
                        1U) {
    return {StatusCode::invalid_plan, detail::halo_detail_field};
  }
  for (std::size_t index = 0U; index < fields.size; ++index) {
    const HaloFieldSpec spec = fields.data[index];
    if (spec.width == 0U || spec.components == 0U ||
        static_cast<std::int32_t>(spec.width) > patch.cells.x ||
        static_cast<std::int32_t>(spec.width) > patch.cells.y ||
        static_cast<std::int32_t>(spec.width) > patch.cells.z) {
      return {StatusCode::invalid_plan, detail::halo_detail_field};
    }
    for (std::size_t prior = 0U; prior < index; ++prior) {
      if (fields.data[prior].field == spec.field) {
        return {StatusCode::invalid_plan, detail::halo_detail_field};
      }
    }
  }
  return {};
}

int neighbor_rank(const MeshPatch& patch, HaloTopology topology, Face face,
                  bool& exists) noexcept {
  const int axis = face_axis(face);
  const std::size_t axis_index = static_cast<std::size_t>(axis);
  const bool minus = face == Face::xm || face == Face::ym || face == Face::zm;
  std::array<std::int32_t, 3U> coordinate{patch.process_coord.x,
                                         patch.process_coord.y,
                                         patch.process_coord.z};
  const std::array<std::int32_t, 3U> grid{patch.process_grid.x,
                                         patch.process_grid.y,
                                         patch.process_grid.z};
  if (minus) {
    if (coordinate[axis_index] == 0) {
      if (!periodic(topology, axis)) {
        exists = false;
        return -1;
      }
      coordinate[axis_index] = grid[axis_index] - 1;
    } else {
      --coordinate[axis_index];
    }
  } else if (coordinate[axis_index] + 1 == grid[axis_index]) {
    if (!periodic(topology, axis)) {
      exists = false;
      return -1;
    }
    coordinate[axis_index] = 0;
  } else {
    ++coordinate[axis_index];
  }
  exists = true;
  const std::int64_t result = coordinate[0] +
                              static_cast<std::int64_t>(grid[0]) *
                                  (coordinate[1] +
                                   static_cast<std::int64_t>(grid[1]) *
                                       coordinate[2]);
  return static_cast<int>(result);
}

Status section_size(const MeshPatch& patch, HaloFieldSpec field, Face face,
                    std::uint64_t& out) noexcept {
  const std::uint64_t width = field.width;
  const std::uint64_t components = field.components;
  const std::uint64_t nx = static_cast<std::uint64_t>(patch.cells.x);
  const std::uint64_t ny = static_cast<std::uint64_t>(patch.cells.y);
  const std::uint64_t nz = static_cast<std::uint64_t>(patch.cells.z);
  std::uint64_t area = 0U;
  if (face_axis(face) == 0) {
    if (!checked_multiply(ny, nz, area)) {
      return {StatusCode::invalid_plan, detail::halo_detail_overflow};
    }
  } else if (face_axis(face) == 1) {
    if (!checked_multiply(nx, nz, area)) {
      return {StatusCode::invalid_plan, detail::halo_detail_overflow};
    }
  } else if (!checked_multiply(nx, ny, area)) {
    return {StatusCode::invalid_plan, detail::halo_detail_overflow};
  }
  std::uint64_t layer = 0U;
  if (!checked_multiply(area, width, layer) ||
      !checked_multiply(layer, components, out) || out == 0U) {
    return {StatusCode::invalid_plan, detail::halo_detail_overflow};
  }
  return {};
}

Status validate_collective_contract(HaloImplData& implementation) noexcept {
  std::uint64_t contract = 1469598103934665603ULL;
  contract = hash_mix(contract,
                      static_cast<std::uint32_t>(implementation.patch.process_grid.x));
  contract = hash_mix(contract,
                      static_cast<std::uint32_t>(implementation.patch.process_grid.y));
  contract = hash_mix(contract,
                      static_cast<std::uint32_t>(implementation.patch.process_grid.z));
  contract = hash_mix(contract, implementation.topology.periodic_x ? 1U : 0U);
  contract = hash_mix(contract, implementation.topology.periodic_y ? 1U : 0U);
  contract = hash_mix(contract, implementation.topology.periodic_z ? 1U : 0U);
  contract = hash_mix(contract, implementation.maximum_chunk_doubles);
  contract = hash_mix(contract,
                      static_cast<std::uint32_t>(implementation.maximum_tag));
  contract = hash_mix(contract, implementation.fields.size());
  for (const FieldSlot& field : implementation.fields) {
    contract = hash_mix(contract, field.spec.field);
    contract = hash_mix(contract, field.spec.width);
    contract = hash_mix(contract, field.spec.components);
  }
  std::uint64_t minimum = contract;
  std::uint64_t maximum = contract;
  const int minimum_result =
      MPI_Allreduce(MPI_IN_PLACE, &minimum, 1, MPI_UINT64_T, MPI_MIN,
                    implementation.control);
  const int maximum_result =
      MPI_Allreduce(MPI_IN_PLACE, &maximum, 1, MPI_UINT64_T, MPI_MAX,
                    implementation.control);
  if (minimum_result != MPI_SUCCESS || maximum_result != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, detail::halo_detail_collective};
  }
  return minimum == maximum
             ? Status{}
             : Status{StatusCode::invalid_plan,
                      detail::halo_detail_topology};
}

Status build_peers(HaloImplData& implementation) {
  constexpr std::array<Face, 6U> faces{Face::xm, Face::xp, Face::ym,
                                       Face::yp, Face::zm, Face::zp};
  for (const Face face : faces) {
    bool exists = false;
    const int neighbor = neighbor_rank(implementation.patch,
                                       implementation.topology, face, exists);
    if (!exists) {
      continue;
    }
    auto peer = std::find_if(implementation.peers.begin(),
                             implementation.peers.end(),
                             [neighbor](const Peer& value) {
                               return value.rank == neighbor;
                             });
    if (peer == implementation.peers.end()) {
      implementation.peers.push_back(Peer{});
      peer = implementation.peers.end() - 1;
      peer->rank = neighbor;
      peer->local = neighbor == implementation.rank;
    }
    for (std::size_t field = 0U; field < implementation.fields.size(); ++field) {
      std::uint64_t count = 0U;
      const Status sized = section_size(implementation.patch,
                                        implementation.fields[field].spec,
                                        face, count);
      if (!sized) {
        return sized;
      }
      peer->sends.push_back(
          Section{field, face, opposite(face), 0U, count});
      peer->receives.push_back(Section{field, face, face, 0U, count});
    }
  }
  std::sort(implementation.peers.begin(), implementation.peers.end(),
            [](const Peer& left, const Peer& right) {
              return left.rank < right.rank;
            });
  const auto section_less = [&implementation](const Section& left,
                                               const Section& right) {
    if (left.canonical_destination != right.canonical_destination) {
      return left.canonical_destination < right.canonical_destination;
    }
    return implementation.fields[left.field_slot].spec.field <
           implementation.fields[right.field_slot].spec.field;
  };
  std::uint64_t send_total = 0U;
  std::uint64_t receive_total = 0U;
  for (std::size_t peer_index = 0U;
       peer_index < implementation.peers.size(); ++peer_index) {
    Peer& peer = implementation.peers[peer_index];
    std::sort(peer.sends.begin(), peer.sends.end(), section_less);
    std::sort(peer.receives.begin(), peer.receives.end(), section_less);
    peer.send_offset = send_total;
    peer.receive_offset = receive_total;
    std::uint64_t peer_send = 0U;
    for (Section& section : peer.sends) {
      section.offset = send_total;
      if (!checked_add(send_total, section.count, send_total) ||
          !checked_add(peer_send, section.count, peer_send)) {
        return {StatusCode::invalid_plan, detail::halo_detail_overflow};
      }
    }
    std::uint64_t peer_receive = 0U;
    for (Section& section : peer.receives) {
      section.offset = receive_total;
      if (!checked_add(receive_total, section.count, receive_total) ||
          !checked_add(peer_receive, section.count, peer_receive)) {
        return {StatusCode::invalid_plan, detail::halo_detail_overflow};
      }
    }
    if (peer_send != peer_receive) {
      return {StatusCode::invalid_plan, detail::halo_detail_topology};
    }
    peer.count = peer_send;
    if (!peer.local) {
      // Chunks are appended in increasing merged-payload offset.  Every rank
      // builds the reciprocal layout in that same canonical destination-face
      // order, and a peer-local unique tag binds every receive to its exact
      // offset independent of MPI_Startall's request activation order.
      std::uint64_t chunk_ordinal = 0U;
      for (std::uint64_t offset = 0U; offset < peer.count;) {
        const std::uint64_t remaining = peer.count - offset;
        const std::uint64_t count =
            std::min(remaining, implementation.maximum_chunk_doubles);
        if (chunk_ordinal >
            static_cast<std::uint64_t>(implementation.maximum_tag)) {
          return {StatusCode::invalid_plan, detail::halo_detail_overflow};
        }
        implementation.chunks.push_back(
            Chunk{peer_index, offset, static_cast<int>(count),
                  static_cast<int>(chunk_ordinal)});
        offset += count;
        ++chunk_ordinal;
      }
    }
  }
  if (send_total > std::numeric_limits<std::size_t>::max() ||
      receive_total > std::numeric_limits<std::size_t>::max() ||
      implementation.chunks.size() >
          static_cast<std::size_t>(std::numeric_limits<int>::max()) / 2U) {
    return {StatusCode::invalid_plan, detail::halo_detail_overflow};
  }
  implementation.send_storage.resize(static_cast<std::size_t>(send_total));
  implementation.receive_storage.resize(
      static_cast<std::size_t>(receive_total));
  implementation.requests.assign(implementation.chunks.size() * 2U,
                                 MPI_REQUEST_NULL);
  implementation.statuses.resize(implementation.requests.size());
  implementation.stats.send_capacity_doubles = send_total;
  implementation.stats.receive_capacity_doubles = receive_total;
  implementation.stats.maximum_messages_per_exchange =
      implementation.chunks.size();
  if (!checked_multiply(send_total, sizeof(double),
                        implementation.stats.maximum_bytes_per_exchange)) {
    return {StatusCode::invalid_plan, detail::halo_detail_overflow};
  }
  implementation.stats.transport_peer_count =
      static_cast<std::size_t>(std::count_if(
          implementation.peers.begin(), implementation.peers.end(),
          [](const Peer& peer) { return !peer.local; }));
  implementation.stats.local_peer_count =
      implementation.peers.size() - implementation.stats.transport_peer_count;
  implementation.stats.persistent_request_count =
      implementation.requests.size();
  return {};
}

Status allocate_peer_size_wires(HaloImplData& implementation,
                                std::vector<std::uint64_t>& sends,
                                std::vector<std::uint64_t>& receives) {
  sends.assign(static_cast<std::size_t>(implementation.size), 0U);
  receives.assign(static_cast<std::size_t>(implementation.size), 0U);
  for (const Peer& peer : implementation.peers) {
    sends[static_cast<std::size_t>(peer.rank)] = peer.count;
  }
  return {};
}

Status exchange_and_validate_peer_sizes(
    HaloImplData& implementation, std::vector<std::uint64_t>& sends,
    std::vector<std::uint64_t>& receives) noexcept {
  if (MPI_Alltoall(sends.data(), 1, MPI_UINT64_T, receives.data(), 1,
                   MPI_UINT64_T, implementation.control) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, detail::halo_detail_collective};
  }
  for (int rank = 0; rank < implementation.size; ++rank) {
    if (receives[static_cast<std::size_t>(rank)] !=
        sends[static_cast<std::size_t>(rank)]) {
      return {StatusCode::invalid_plan, detail::halo_detail_topology};
    }
  }
  return {};
}

Status create_requests(HaloImplData& implementation) noexcept {
  const std::size_t send_base = implementation.chunks.size();
  for (std::size_t index = 0U; index < implementation.chunks.size(); ++index) {
    const Chunk& chunk = implementation.chunks[index];
    const Peer& peer = implementation.peers[chunk.peer];
    double* const receive = implementation.receive_storage.data() +
                            static_cast<std::size_t>(peer.receive_offset +
                                                     chunk.offset);
    double* const send = implementation.send_storage.data() +
                         static_cast<std::size_t>(peer.send_offset +
                                                  chunk.offset);
    if (MPI_Recv_init(receive, chunk.count, MPI_DOUBLE, peer.rank, chunk.tag,
                      implementation.transport,
                      &implementation.requests[index]) != MPI_SUCCESS ||
        MPI_Send_init(send, chunk.count, MPI_DOUBLE, peer.rank, chunk.tag,
                      implementation.transport,
                      &implementation.requests[send_base + index]) !=
            MPI_SUCCESS) {
      return {StatusCode::mpi_failure, detail::halo_detail_request};
    }
  }
  implementation.stats.request_storage_address =
      reinterpret_cast<std::uintptr_t>(implementation.requests.data());
  implementation.stats.send_storage_address =
      reinterpret_cast<std::uintptr_t>(implementation.send_storage.data());
  implementation.stats.receive_storage_address =
      reinterpret_cast<std::uintptr_t>(implementation.receive_storage.data());
  return {};
}

void free_requests(HaloImplData& implementation) noexcept {
  for (MPI_Request& request : implementation.requests) {
    if (request != MPI_REQUEST_NULL) {
      (void)MPI_Request_free(&request);
    }
  }
}

Status recreate_requests(HaloImplData& implementation) noexcept {
  free_requests(implementation);
  std::fill(implementation.requests.begin(), implementation.requests.end(),
            MPI_REQUEST_NULL);
  implementation.requests_active = false;
  const Status created = create_requests(implementation);
  implementation.poisoned = !created;
  return created;
}

bool valid_view_shape(const FieldView& view, const FieldSlot& field,
                      const MeshPatch& patch) noexcept {
  if (view.base == nullptr || view.field != field.spec.field ||
      view.interior.x != patch.cells.x || view.interior.y != patch.cells.y ||
      view.interior.z != patch.cells.z ||
      view.ghosts.x < static_cast<std::int32_t>(field.spec.width) ||
      view.ghosts.y < static_cast<std::int32_t>(field.spec.width) ||
      view.ghosts.z < static_cast<std::int32_t>(field.spec.width) ||
      view.components != field.spec.components || view.revision == 0U ||
      view.storage_identity == 0U || view.revision_domain == 0U) {
    return false;
  }
  const std::uint64_t padded_x =
      static_cast<std::uint64_t>(view.interior.x) +
      2U * static_cast<std::uint64_t>(view.ghosts.x);
  const std::uint64_t padded_y =
      static_cast<std::uint64_t>(view.interior.y) +
      2U * static_cast<std::uint64_t>(view.ghosts.y);
  const std::uint64_t padded_z =
      static_cast<std::uint64_t>(view.interior.z) +
      2U * static_cast<std::uint64_t>(view.ghosts.z);
  std::uint64_t z_minimum = 0U;
  std::uint64_t component_minimum = 0U;
  std::uint64_t component_extent = 0U;
  return view.stride_y >= padded_x &&
         checked_multiply(view.stride_y, padded_y, z_minimum) &&
         view.stride_z >= z_minimum &&
         checked_multiply(view.stride_z, padded_z, component_minimum) &&
         view.component_stride >= component_minimum &&
         checked_multiply(view.component_stride,
                          static_cast<std::uint64_t>(view.components),
                          component_extent) &&
         component_extent <= static_cast<std::uint64_t>(
                                 std::numeric_limits<std::ptrdiff_t>::max());
}

Status validate_begin_views(HaloImplData& implementation,
                            Span<const FieldView> views) noexcept {
  if (views.data == nullptr || views.size != implementation.fields.size()) {
    return {StatusCode::invalid_plan, detail::halo_detail_field};
  }
  ++implementation.field_visit_generation;
  if (implementation.field_visit_generation == 0U) {
    implementation.field_visit_generation = 1U;
    for (FieldSlot& field : implementation.fields) {
      field.seen_generation = 0U;
    }
  }
  const std::uint64_t visit = implementation.field_visit_generation;
  for (std::size_t input = 0U; input < views.size; ++input) {
    const FieldView view = views.data[input];
    const std::int32_t lookup =
        implementation.field_lookup[static_cast<std::size_t>(view.field)];
    if (lookup < 0) {
      return {StatusCode::invalid_plan, detail::halo_detail_field};
    }
    FieldSlot& field =
        implementation.fields[static_cast<std::size_t>(lookup)];
    if (field.seen_generation == visit ||
        !valid_view_shape(view, field, implementation.patch)) {
      return {StatusCode::invalid_plan, detail::halo_detail_field};
    }
    field.seen_generation = visit;
    implementation.candidate_begin_views[static_cast<std::size_t>(lookup)] =
        view;
  }
  return {};
}

void publish_begin_views(HaloImplData& implementation) noexcept {
  for (std::size_t field = 0U; field < implementation.fields.size(); ++field) {
    implementation.fields[field].begin_view =
        implementation.candidate_begin_views[field];
    implementation.fields[field].pending_revision =
        implementation.candidate_begin_views[field].revision;
  }
}

Status bind_finish_views(HaloImplData& implementation,
                         Span<FieldView> views) noexcept {
  if (views.data == nullptr || views.size != implementation.fields.size()) {
    return {StatusCode::invalid_plan, detail::halo_detail_field};
  }
  ++implementation.field_visit_generation;
  if (implementation.field_visit_generation == 0U) {
    implementation.field_visit_generation = 1U;
    for (FieldSlot& field : implementation.fields) {
      field.seen_generation = 0U;
    }
  }
  const std::uint64_t visit = implementation.field_visit_generation;
  for (std::size_t input = 0U; input < views.size; ++input) {
    const FieldView view = views.data[input];
    const std::int32_t lookup =
        implementation.field_lookup[static_cast<std::size_t>(view.field)];
    if (lookup < 0) {
      return {StatusCode::invalid_plan, detail::halo_detail_field};
    }
    FieldSlot& field =
        implementation.fields[static_cast<std::size_t>(lookup)];
    const FieldView source = field.begin_view;
    if (field.seen_generation == visit ||
        !valid_view_shape(view, field, implementation.patch) ||
        view.base != source.base || view.stride_y != source.stride_y ||
        view.stride_z != source.stride_z ||
        view.component_stride != source.component_stride ||
        view.replica != source.replica || view.revision != source.revision ||
        view.storage_identity != source.storage_identity ||
        view.revision_domain != source.revision_domain) {
      return {StatusCode::invalid_plan, detail::halo_detail_field};
    }
    field.seen_generation = visit;
    field.finish_view = view;
  }
  return {};
}

void pack_section(const FieldView& view, HaloFieldSpec spec, Face face,
                  double* output) noexcept {
  std::size_t cursor = 0U;
  const std::int32_t nx = view.interior.x;
  const std::int32_t ny = view.interior.y;
  const std::int32_t nz = view.interior.z;
  const std::int32_t width = spec.width;
  for (std::uint8_t component_index = 0U;
       component_index < spec.components; ++component_index) {
    if (face == Face::xm || face == Face::xp) {
      const std::int32_t start = face == Face::xm ? 0 : nx - width;
      for (std::int32_t z = 0; z < nz; ++z) {
        for (std::int32_t y = 0; y < ny; ++y) {
          for (std::int32_t offset = 0; offset < width; ++offset) {
            output[cursor++] =
                view.unchecked(Int3{start + offset, y, z}, component_index);
          }
        }
      }
    } else if (face == Face::ym || face == Face::yp) {
      const std::int32_t start = face == Face::ym ? 0 : ny - width;
      for (std::int32_t z = 0; z < nz; ++z) {
        for (std::int32_t offset = 0; offset < width; ++offset) {
          for (std::int32_t x = 0; x < nx; ++x) {
            output[cursor++] =
                view.unchecked(Int3{x, start + offset, z}, component_index);
          }
        }
      }
    } else {
      const std::int32_t start = face == Face::zm ? 0 : nz - width;
      for (std::int32_t offset = 0; offset < width; ++offset) {
        for (std::int32_t y = 0; y < ny; ++y) {
          for (std::int32_t x = 0; x < nx; ++x) {
            output[cursor++] =
                view.unchecked(Int3{x, y, start + offset}, component_index);
          }
        }
      }
    }
  }
}

void unpack_section(FieldView view, HaloFieldSpec spec, Face face,
                    const double* input) noexcept {
  std::size_t cursor = 0U;
  const std::int32_t nx = view.interior.x;
  const std::int32_t ny = view.interior.y;
  const std::int32_t nz = view.interior.z;
  const std::int32_t width = spec.width;
  for (std::uint8_t component_index = 0U;
       component_index < spec.components; ++component_index) {
    if (face == Face::xm || face == Face::xp) {
      const std::int32_t start = face == Face::xm ? -width : nx;
      for (std::int32_t z = 0; z < nz; ++z) {
        for (std::int32_t y = 0; y < ny; ++y) {
          for (std::int32_t offset = 0; offset < width; ++offset) {
            view.unchecked(Int3{start + offset, y, z}, component_index) =
                input[cursor++];
          }
        }
      }
    } else if (face == Face::ym || face == Face::yp) {
      const std::int32_t start = face == Face::ym ? -width : ny;
      for (std::int32_t z = 0; z < nz; ++z) {
        for (std::int32_t offset = 0; offset < width; ++offset) {
          for (std::int32_t x = 0; x < nx; ++x) {
            view.unchecked(Int3{x, start + offset, z}, component_index) =
                input[cursor++];
          }
        }
      }
    } else {
      const std::int32_t start = face == Face::zm ? -width : nz;
      for (std::int32_t offset = 0; offset < width; ++offset) {
        for (std::int32_t y = 0; y < ny; ++y) {
          for (std::int32_t x = 0; x < nx; ++x) {
            view.unchecked(Int3{x, y, start + offset}, component_index) =
                input[cursor++];
          }
        }
      }
    }
  }
}

void copy_local_peers(HaloImplData& implementation) noexcept {
  for (const Peer& peer : implementation.peers) {
    if (peer.local && peer.count != 0U) {
      std::memcpy(implementation.receive_storage.data() +
                      static_cast<std::size_t>(peer.receive_offset),
                  implementation.send_storage.data() +
                      static_cast<std::size_t>(peer.send_offset),
                  static_cast<std::size_t>(peer.count) * sizeof(double));
    }
  }
}

void finish_exchange_state(HaloImplData& implementation) noexcept {
  implementation.exchange_active = false;
  implementation.requests_active = false;
}

void cancel_and_drain(HaloImplData& implementation) noexcept {
  if (!implementation.requests_active || implementation.requests.empty()) {
    implementation.requests_active = false;
    return;
  }
  for (MPI_Request& request : implementation.requests) {
    (void)MPI_Cancel(&request);
  }
  (void)MPI_Waitall(static_cast<int>(implementation.requests.size()),
                    implementation.requests.data(),
                    implementation.statuses.data());
  implementation.requests_active = false;
}

void release_data(HaloImplData& implementation) noexcept {
  int initialized = 0;
  int finalized = 0;
  (void)MPI_Initialized(&initialized);
  if (initialized != 0) {
    (void)MPI_Finalized(&finalized);
  }
  if (initialized != 0 && finalized == 0) {
    cancel_and_drain(implementation);
    free_requests(implementation);
    if (implementation.transport != MPI_COMM_NULL) {
      (void)MPI_Comm_free(&implementation.transport);
    }
    if (implementation.control != MPI_COMM_NULL) {
      (void)MPI_Comm_free(&implementation.control);
    }
  }
  implementation.transport = MPI_COMM_NULL;
  implementation.control = MPI_COMM_NULL;
}

}  // namespace

struct HaloEngine::Impl : HaloImplData {};

namespace detail {

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
void set_halo_failure_for_test(HaloFailurePoint point,
                               int failing_rank) noexcept {
  g_failure_rank.store(failing_rank, std::memory_order_relaxed);
  g_failure_point.store(point, std::memory_order_relaxed);
}

void clear_halo_failure_for_test() noexcept {
  g_failure_point.store(HaloFailurePoint::none, std::memory_order_relaxed);
  g_failure_rank.store(-1, std::memory_order_relaxed);
}

void set_halo_maximum_chunk_doubles_for_test(
    std::uint64_t maximum) noexcept {
  g_test_maximum_chunk_doubles.store(maximum, std::memory_order_relaxed);
}

void clear_halo_maximum_chunk_doubles_for_test() noexcept {
  g_test_maximum_chunk_doubles.store(kMaximumChunkDoubles,
                                     std::memory_order_relaxed);
}

void set_halo_tag_upper_bound_for_test(int upper_bound) noexcept {
  g_test_tag_upper_bound.store(upper_bound, std::memory_order_relaxed);
}

void clear_halo_tag_upper_bound_for_test() noexcept {
  g_test_tag_upper_bound.store(std::numeric_limits<int>::min(),
                               std::memory_order_relaxed);
}
#endif

}  // namespace detail

HaloEngine::~HaloEngine() noexcept { release(); }

HaloEngine::HaloEngine(HaloEngine&& other) noexcept
    : implementation_(other.implementation_) {
  other.implementation_ = nullptr;
}

HaloEngine& HaloEngine::operator=(HaloEngine&& other) noexcept {
  if (this != &other) {
    release();
    implementation_ = other.implementation_;
    other.implementation_ = nullptr;
  }
  return *this;
}

void HaloEngine::release() noexcept {
  if (implementation_ != nullptr) {
    release_data(*implementation_);
    delete implementation_;
    implementation_ = nullptr;
  }
}

Status HaloEngine::reserve(MPI_Comm communicator, const MeshPatch& patch,
                           Span<const HaloFieldSpec> fields) noexcept {
  return reserve(communicator, patch, fields, HaloTopology{});
}

Status HaloEngine::reserve(MPI_Comm communicator, const MeshPatch& patch,
                           Span<const HaloFieldSpec> fields,
                           HaloTopology topology) noexcept {
  if (communicator == MPI_COMM_NULL) {
    return {StatusCode::invalid_plan, detail::halo_detail_input};
  }
  int rank = -1;
  int size = 0;
  if (MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(communicator, &size) != MPI_SUCCESS || size <= 0) {
    return {StatusCode::mpi_failure, detail::halo_detail_collective};
  }
  std::unique_ptr<Impl> candidate(new (std::nothrow) Impl);
  Status local = candidate == nullptr
                     ? Status{StatusCode::allocation_failure, 0U}
                     : validate_patch(patch, rank, size);
  if (local && implementation_ != nullptr &&
      implementation_->exchange_active) {
    local = {StatusCode::invalid_plan, detail::halo_detail_state};
  }
  if (local) {
    local = validate_specs(fields, patch);
  }
  int lowest = -1;
  Status consensus =
      raw_collective_status(communicator, rank, size, local, lowest);
  if (!consensus) {
    if (implementation_ != nullptr) {
      implementation_->lowest_failing_rank = lowest;
    }
    return consensus;
  }
  candidate->rank = rank;
  candidate->size = size;
  candidate->patch = patch;
  candidate->topology = topology;
  local = MPI_Comm_dup(communicator, &candidate->control) == MPI_SUCCESS
              ? Status{}
              : Status{StatusCode::mpi_failure,
                       detail::halo_detail_collective};
  consensus = raw_collective_status(communicator, rank, size, local, lowest);
  if (!consensus) {
    release_data(*candidate);
    if (implementation_ != nullptr) {
      implementation_->lowest_failing_rank = lowest;
    }
    return consensus;
  }
  const auto reject_candidate =
      [this, &candidate](Status failure) noexcept -> Status {
    const int selected_rank = candidate->lowest_failing_rank;
    release_data(*candidate);
    if (implementation_ != nullptr) {
      implementation_->lowest_failing_rank = selected_rank;
    }
    return failure;
  };
  local = MPI_Comm_set_errhandler(candidate->control, MPI_ERRORS_RETURN) ==
                  MPI_SUCCESS
              ? Status{}
              : Status{StatusCode::mpi_failure,
                       detail::halo_detail_collective};
  consensus = raw_collective_status(communicator, rank, size, local, lowest);
  if (!consensus) {
    release_data(*candidate);
    if (implementation_ != nullptr) {
      implementation_->lowest_failing_rank = lowest;
    }
    return consensus;
  }
  local = MPI_Comm_dup(communicator, &candidate->transport) == MPI_SUCCESS
              ? Status{}
              : Status{StatusCode::mpi_failure,
                       detail::halo_detail_collective};
  consensus = raw_collective_status(communicator, rank, size, local, lowest);
  if (!consensus) {
    release_data(*candidate);
    if (implementation_ != nullptr) {
      implementation_->lowest_failing_rank = lowest;
    }
    return consensus;
  }
  local = MPI_Comm_set_errhandler(candidate->transport, MPI_ERRORS_RETURN) ==
                  MPI_SUCCESS
              ? Status{}
              : Status{StatusCode::mpi_failure,
                       detail::halo_detail_collective};
  consensus = raw_collective_status(communicator, rank, size, local, lowest);
  if (!consensus) {
    release_data(*candidate);
    if (implementation_ != nullptr) {
      implementation_->lowest_failing_rank = lowest;
    }
    return consensus;
  }
  int* maximum_tag_pointer = nullptr;
  int maximum_tag_available = 0;
  const int tag_query = MPI_Comm_get_attr(candidate->transport, MPI_TAG_UB,
                                          &maximum_tag_pointer,
                                          &maximum_tag_available);
  local = tag_query == MPI_SUCCESS
              ? Status{}
              : Status{StatusCode::mpi_failure,
                       detail::halo_detail_collective};
  int injected_tag_upper_bound = 0;
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  if (configured_tag_upper_bound(injected_tag_upper_bound)) {
    if (injected_tag_upper_bound < 0) {
      maximum_tag_pointer = nullptr;
      maximum_tag_available = 0;
    } else {
      maximum_tag_pointer = &injected_tag_upper_bound;
      maximum_tag_available = 1;
    }
  }
#else
  (void)injected_tag_upper_bound;
#endif
  // MPI guarantees tags through 32767.  Some older MPI builds lose this
  // predefined attribute on a duplicated communicator, so a successful query
  // without a value falls back to that portable bound.  Valid reported bounds
  // are globally minimized before entering the plan contract or requests.
  constexpr int kMinimumTagUpperBound = 32767;
  const int local_maximum_tag =
      maximum_tag_available != 0 && maximum_tag_pointer != nullptr &&
              *maximum_tag_pointer >= 0
          ? *maximum_tag_pointer
          : kMinimumTagUpperBound;
  int global_maximum_tag = kMinimumTagUpperBound;
  const int tag_reduce = MPI_Allreduce(
      &local_maximum_tag, &global_maximum_tag, 1, MPI_INT, MPI_MIN,
      candidate->control);
  if (tag_reduce != MPI_SUCCESS && local) {
    local = {StatusCode::mpi_failure, detail::halo_detail_collective};
  }
  if (local) {
    candidate->maximum_tag = global_maximum_tag;
    candidate->stats.maximum_tag = global_maximum_tag;
    candidate->maximum_chunk_doubles = configured_maximum_chunk_doubles();
    if (candidate->maximum_chunk_doubles == 0U ||
        candidate->maximum_chunk_doubles > kMaximumChunkDoubles) {
      local = {StatusCode::invalid_plan, detail::halo_detail_overflow};
    }
  }
  consensus = collective_status(*candidate, local);
  if (!consensus) {
    return reject_candidate(consensus);
  }
  // From here onward every potentially allocating local phase reaches a
  // control consensus before any rank enters the next MPI collective.
  try {
    candidate->fields.reserve(fields.size);
    for (std::size_t index = 0U; index < fields.size; ++index) {
      candidate->fields.push_back(FieldSlot{fields.data[index]});
    }
    std::sort(candidate->fields.begin(), candidate->fields.end(),
              [](const FieldSlot& left, const FieldSlot& right) {
                return left.spec.field < right.spec.field;
              });
    candidate->field_lookup.assign(
        static_cast<std::size_t>(std::numeric_limits<FieldId>::max()) + 1U,
        -1);
    candidate->candidate_begin_views.resize(candidate->fields.size());
    for (std::size_t index = 0U; index < candidate->fields.size(); ++index) {
      candidate->field_lookup[candidate->fields[index].spec.field] =
          static_cast<std::int32_t>(index);
    }
    if (injected(FailurePoint::reserve_before_contract, rank)) {
      local = {StatusCode::allocation_failure,
               detail::halo_detail_reserve_allocation};
    }
  } catch (const std::bad_alloc&) {
    local = {StatusCode::allocation_failure, 0U};
  } catch (...) {
    local = {StatusCode::invalid_plan, detail::halo_detail_state};
  }
  consensus = collective_status(*candidate, local);
  if (!consensus) {
    return reject_candidate(consensus);
  }
  // All ranks have completed allocation before these two contract Allreduces.
  local = validate_collective_contract(*candidate);
  consensus = collective_status(*candidate, local);
  if (!consensus) {
    return reject_candidate(consensus);
  }
  try {
    local = build_peers(*candidate);
  } catch (const std::bad_alloc&) {
    local = {StatusCode::allocation_failure, 0U};
  } catch (...) {
    local = {StatusCode::invalid_plan, detail::halo_detail_state};
  }
  consensus = collective_status(*candidate, local);
  if (!consensus) {
    return reject_candidate(consensus);
  }
  std::vector<std::uint64_t> peer_size_sends;
  std::vector<std::uint64_t> peer_size_receives;
  try {
    local = allocate_peer_size_wires(*candidate, peer_size_sends,
                                     peer_size_receives);
    if (local &&
        injected(FailurePoint::reserve_before_alltoall, rank)) {
      local = {StatusCode::allocation_failure,
               detail::halo_detail_reserve_allocation};
    }
  } catch (const std::bad_alloc&) {
    local = {StatusCode::allocation_failure, 0U};
  } catch (...) {
    local = {StatusCode::invalid_plan, detail::halo_detail_state};
  }
  consensus = collective_status(*candidate, local);
  if (!consensus) {
    return reject_candidate(consensus);
  }
  // No rank enters Alltoall until all rank-local wire allocations succeeded.
  local = exchange_and_validate_peer_sizes(*candidate, peer_size_sends,
                                           peer_size_receives);
  consensus = collective_status(*candidate, local);
  if (!consensus) {
    return reject_candidate(consensus);
  }
  local = create_requests(*candidate);
  consensus = collective_status(*candidate, local);
  if (!consensus) {
    return reject_candidate(consensus);
  }
  candidate->identity =
      g_next_halo_identity.fetch_add(1U, std::memory_order_relaxed);
  if (candidate->identity == 0U) {
    candidate->identity =
        g_next_halo_identity.fetch_add(1U, std::memory_order_relaxed);
  }
  if (candidate->identity == 0U) {
    local = {StatusCode::invalid_plan, detail::halo_detail_state};
  }
  consensus = collective_status(*candidate, local);
  if (!consensus) {
    return reject_candidate(consensus);
  }
  candidate->lowest_failing_rank = -1;
  release();
  implementation_ = candidate.release();
  return {};
}

Status HaloEngine::begin(StageId stage, Span<const FieldView> fields,
                         HaloTicket& ticket) noexcept {
  return begin(stage, fields, {}, ticket);
}

Status HaloEngine::begin(StageId stage, Span<const FieldView> fields,
                         Status prerequisite, HaloTicket& ticket) noexcept {
  if (implementation_ == nullptr) {
    return {StatusCode::invalid_plan, detail::halo_detail_state};
  }
  HaloImplData& implementation = *implementation_;
  ++implementation.counters.begin_calls;
  // A new idle attempt invalidates an older certificate immediately, even if
  // its pack/preflight later fails.  A re-entrant begin must leave the active
  // exchange's pending metadata and certificate boundary untouched.
  if (!implementation.exchange_active) {
    invalidate_certificates(implementation);
  }
  Status local =
      implementation.exchange_active || ticket.active_ ||
              implementation.prepared_epoch_active
          ? Status{StatusCode::invalid_plan, detail::halo_detail_state}
          : prerequisite;
  if (local && implementation.poisoned) {
    local = recreate_requests(implementation);
  }
  if (local) {
    local = validate_begin_views(implementation, fields);
  }
  if (local && injected(FailurePoint::pack, implementation.rank)) {
    local = {StatusCode::invalid_plan, detail::halo_detail_pack_failure};
  }
  if (local) {
    for (const Peer& peer : implementation.peers) {
      for (const Section& section : peer.sends) {
        const FieldSlot& field = implementation.fields[section.field_slot];
        pack_section(implementation.candidate_begin_views[section.field_slot],
                     field.spec, section.face,
                     implementation.send_storage.data() +
                         static_cast<std::size_t>(section.offset));
      }
    }
  }
  if (local && injected(FailurePoint::start, implementation.rank)) {
    local = {StatusCode::mpi_failure, detail::halo_detail_start_failure};
  }
  Status consensus = collective_status(implementation, local);
  if (!consensus) {
    if (consensus.detail == detail::halo_detail_start_failure) {
      implementation.counters.bytes_packed +=
          implementation.stats.send_capacity_doubles * sizeof(double);
    }
    ticket.active_ = false;
    return consensus;
  }
  implementation.counters.bytes_packed +=
      implementation.stats.send_capacity_doubles * sizeof(double);
  local = {};
  if (!implementation.requests.empty()) {
    implementation.requests_active = true;
    if (MPI_Startall(static_cast<int>(implementation.requests.size()),
                     implementation.requests.data()) != MPI_SUCCESS) {
      local = {StatusCode::mpi_failure, detail::halo_detail_start_failure};
    }
  }
  consensus = collective_status(implementation, local);
  if (!consensus) {
    const int selected_failure = implementation.lowest_failing_rank;
    cancel_and_drain(implementation);
    const Status restored = recreate_requests(implementation);
    const Status restore_consensus = collective_status(implementation, restored);
    if (!restore_consensus) {
      implementation.poisoned = true;
    }
    implementation.lowest_failing_rank = selected_failure;
    implementation.exchange_active = false;
    ticket.active_ = false;
    return consensus;
  }
  implementation.counters.messages_started += implementation.chunks.size();
  publish_begin_views(implementation);
  ++implementation.generation;
  if (implementation.generation == 0U) {
    ++implementation.generation;
  }
  implementation.exchange_active = true;
  ticket.engine_identity_ = implementation.identity;
  ticket.generation_ = implementation.generation;
  ticket.stage_ = stage;
  ticket.active_ = true;
  return {};
}

Status HaloEngine::finish(HaloTicket& ticket,
                          Span<FieldView> fields) noexcept {
  if (implementation_ == nullptr) {
    return {StatusCode::invalid_plan, detail::halo_detail_state};
  }
  HaloImplData& implementation = *implementation_;
  ++implementation.counters.finish_calls;
  Status local{};
  if (!implementation.exchange_active || implementation.poisoned ||
      !ticket.active_ || ticket.engine_identity_ != implementation.identity ||
      ticket.generation_ != implementation.generation) {
    local = {StatusCode::invalid_plan, detail::halo_detail_state};
  } else {
    local = bind_finish_views(implementation, fields);
  }
  if (local && injected(FailurePoint::completion, implementation.rank)) {
    local = {StatusCode::mpi_failure,
             detail::halo_detail_completion_failure};
  }
  const Status preflight = collective_status(implementation, local);
  const int preflight_failure_rank = implementation.lowest_failing_rank;
  Status completion{};
  if (implementation.requests_active) {
    if (MPI_Waitall(static_cast<int>(implementation.requests.size()),
                    implementation.requests.data(),
                    implementation.statuses.data()) != MPI_SUCCESS) {
      completion = {StatusCode::mpi_failure,
                    detail::halo_detail_completion_failure};
    } else {
      implementation.requests_active = false;
    }
  }
  if (completion && preflight &&
      injected(FailurePoint::unpack, implementation.rank)) {
    completion = {StatusCode::invalid_plan,
                  detail::halo_detail_unpack_failure};
  }
  const Status completion_consensus =
      collective_status(implementation, completion);
  const int completion_failure_rank = implementation.lowest_failing_rank;
  if (!completion_consensus &&
      completion_consensus.detail == detail::halo_detail_completion_failure &&
      completion_consensus.code == StatusCode::mpi_failure) {
    cancel_and_drain(implementation);
    const Status restored = recreate_requests(implementation);
    const Status restore_consensus = collective_status(implementation, restored);
    if (!restore_consensus) {
      implementation.poisoned = true;
    }
    implementation.lowest_failing_rank = completion_failure_rank;
  }
  if (!preflight || !completion_consensus) {
    const Status result = !preflight ? preflight : completion_consensus;
    const int selected_failure =
        !preflight ? preflight_failure_rank : completion_failure_rank;
    finish_exchange_state(implementation);
    ticket = HaloTicket{};
    implementation.lowest_failing_rank = selected_failure;
    return result;
  }
  copy_local_peers(implementation);
  for (const Peer& peer : implementation.peers) {
    for (const Section& section : peer.receives) {
      FieldSlot& field = implementation.fields[section.field_slot];
      unpack_section(field.finish_view, field.spec, section.face,
                     implementation.receive_storage.data() +
                         static_cast<std::size_t>(section.offset));
    }
  }
  implementation.counters.bytes_unpacked +=
      implementation.stats.receive_capacity_doubles * sizeof(double);
  for (FieldSlot& field : implementation.fields) {
    field.ghost_revision = field.pending_revision;
  }
  implementation.lowest_failing_rank = -1;
  finish_exchange_state(implementation);
  ticket = HaloTicket{};
  return {};
}

Status HaloEngine::enter_prepared_epoch() noexcept {
  if (implementation_ == nullptr) {
    return {StatusCode::invalid_plan, detail::halo_detail_state};
  }
  HaloImplData& implementation = *implementation_;
  // An attempted epoch is itself a certificate boundary.  Invalidate before
  // checking readiness so a rank-local rejected entry cannot retain an older
  // ghost certificate while its peers reject the same outer authority.
  invalidate_certificates(implementation);
  for (FieldSlot& field : implementation.fields) {
    field.pending_revision = 0U;
  }
  if (implementation.prepared_epoch_active || implementation.exchange_active ||
      implementation.requests_active || implementation.poisoned) {
    return {StatusCode::invalid_plan, detail::halo_detail_state};
  }
  implementation.prepared_exchange_views_valid = false;
  implementation.prepared_epoch_active = true;
  implementation.lowest_failing_rank = -1;
  return {};
}

Status HaloEngine::begin_prepared(StageId stage,
                                  Span<const FieldView> fields,
                                  Status& deferred,
                                  HaloTicket& ticket) noexcept {
  if (implementation_ == nullptr) {
    return {StatusCode::invalid_plan, detail::halo_detail_state};
  }
  HaloImplData& implementation = *implementation_;
  ++implementation.counters.begin_calls;
  Status local =
      !implementation.prepared_epoch_active ||
              implementation.exchange_active || ticket.active_ ||
              implementation.poisoned
          ? Status{StatusCode::invalid_plan, detail::halo_detail_state}
          : Status{};
  if (local) {
    local = validate_begin_views(implementation, fields);
  }
  const bool views_valid = static_cast<bool>(local);
  retain_deferred(local, deferred);
  if (views_valid) {
    if (injected(FailurePoint::pack, implementation.rank)) {
      retain_deferred(
          {StatusCode::invalid_plan, detail::halo_detail_pack_failure},
          deferred);
    }
    for (const Peer& peer : implementation.peers) {
      for (const Section& section : peer.sends) {
        const FieldSlot& field = implementation.fields[section.field_slot];
        pack_section(implementation.candidate_begin_views[section.field_slot],
                     field.spec, section.face,
                     implementation.send_storage.data() +
                         static_cast<std::size_t>(section.offset));
      }
    }
    implementation.counters.bytes_packed +=
        implementation.stats.send_capacity_doubles * sizeof(double);
  }
  if (injected(FailurePoint::start, implementation.rank)) {
    retain_deferred(
        {StatusCode::mpi_failure, detail::halo_detail_start_failure},
        deferred);
  }
  if (!implementation.requests.empty()) {
    implementation.requests_active = true;
    if (MPI_Startall(static_cast<int>(implementation.requests.size()),
                     implementation.requests.data()) != MPI_SUCCESS) {
      cancel_and_drain(implementation);
      implementation.poisoned = true;
      implementation.exchange_active = false;
      implementation.prepared_exchange_views_valid = false;
      ticket = HaloTicket{};
      return {StatusCode::mpi_failure, detail::halo_detail_start_failure};
    }
  }
  implementation.counters.messages_started += implementation.chunks.size();
  if (views_valid) {
    publish_begin_views(implementation);
  }
  ++implementation.generation;
  if (implementation.generation == 0U) {
    ++implementation.generation;
  }
  implementation.prepared_exchange_views_valid = views_valid;
  implementation.exchange_active = true;
  ticket.engine_identity_ = implementation.identity;
  ticket.generation_ = implementation.generation;
  ticket.stage_ = stage;
  ticket.active_ = true;
  return {};
}

Status HaloEngine::finish_prepared(HaloTicket& ticket,
                                   Span<FieldView> fields,
                                   Status& deferred) noexcept {
  if (implementation_ == nullptr) {
    return {StatusCode::invalid_plan, detail::halo_detail_state};
  }
  HaloImplData& implementation = *implementation_;
  ++implementation.counters.finish_calls;
  Status local =
      !implementation.prepared_epoch_active ||
              !implementation.exchange_active || implementation.poisoned ||
              !ticket.active_ ||
              ticket.engine_identity_ != implementation.identity ||
              ticket.generation_ != implementation.generation
          ? Status{StatusCode::invalid_plan, detail::halo_detail_state}
          : Status{};
  if (local && implementation.prepared_exchange_views_valid) {
    local = bind_finish_views(implementation, fields);
  }
  const bool views_valid = implementation.prepared_exchange_views_valid &&
                           static_cast<bool>(local);
  retain_deferred(local, deferred);
  if (injected(FailurePoint::completion, implementation.rank)) {
    retain_deferred(
        {StatusCode::mpi_failure, detail::halo_detail_completion_failure},
        deferred);
  }
  if (implementation.requests_active) {
    if (MPI_Waitall(static_cast<int>(implementation.requests.size()),
                    implementation.requests.data(),
                    implementation.statuses.data()) != MPI_SUCCESS) {
      cancel_and_drain(implementation);
      implementation.poisoned = true;
      finish_exchange_state(implementation);
      implementation.prepared_exchange_views_valid = false;
      ticket = HaloTicket{};
      return {StatusCode::mpi_failure,
              detail::halo_detail_completion_failure};
    }
    implementation.requests_active = false;
  }
  if (injected(FailurePoint::unpack, implementation.rank)) {
    retain_deferred(
        {StatusCode::invalid_plan, detail::halo_detail_unpack_failure},
        deferred);
  }
  if (views_valid) {
    copy_local_peers(implementation);
    for (const Peer& peer : implementation.peers) {
      for (const Section& section : peer.receives) {
        FieldSlot& field = implementation.fields[section.field_slot];
        unpack_section(field.finish_view, field.spec, section.face,
                       implementation.receive_storage.data() +
                           static_cast<std::size_t>(section.offset));
      }
    }
    implementation.counters.bytes_unpacked +=
        implementation.stats.receive_capacity_doubles * sizeof(double);
  }
  finish_exchange_state(implementation);
  implementation.prepared_exchange_views_valid = false;
  ticket = HaloTicket{};
  return {};
}

void HaloEngine::close_prepared_epoch(bool publish,
                                      int lowest_failing_rank) noexcept {
  if (implementation_ == nullptr) {
    return;
  }
  HaloImplData& implementation = *implementation_;
  if (!implementation.prepared_epoch_active) {
    return;
  }
  const bool clean = !implementation.exchange_active &&
                     !implementation.requests_active &&
                     !implementation.poisoned;
  if (publish && clean) {
    for (FieldSlot& field : implementation.fields) {
      field.ghost_revision = field.pending_revision;
    }
    implementation.lowest_failing_rank = -1;
  } else {
    invalidate_certificates(implementation);
    implementation.lowest_failing_rank = lowest_failing_rank;
  }
  implementation.prepared_epoch_active = false;
  implementation.prepared_exchange_views_valid = false;
}

Status HaloEngine::validate_contract(
    MPI_Comm communicator, const MeshPatch& patch,
    Span<const HaloFieldSpec> fields, HaloTopology topology) const noexcept {
  if (implementation_ == nullptr || communicator == MPI_COMM_NULL ||
      fields.data == nullptr || fields.size != implementation_->fields.size()) {
    return {StatusCode::invalid_plan, detail::halo_detail_input};
  }
  int initialized = 0;
  int finalized = 0;
  if (MPI_Initialized(&initialized) != MPI_SUCCESS || initialized == 0 ||
      MPI_Finalized(&finalized) != MPI_SUCCESS || finalized != 0) {
    return {StatusCode::invalid_plan, detail::halo_detail_state};
  }
  int relation = MPI_UNEQUAL;
  if (MPI_Comm_compare(implementation_->control, communicator, &relation) !=
      MPI_SUCCESS) {
    return {StatusCode::mpi_failure, detail::halo_detail_collective};
  }
  if (relation != MPI_IDENT && relation != MPI_CONGRUENT) {
    return {StatusCode::invalid_plan, detail::halo_detail_topology};
  }
  const MeshPatch& expected = implementation_->patch;
  if (patch.begin.x != expected.begin.x || patch.begin.y != expected.begin.y ||
      patch.begin.z != expected.begin.z || patch.cells.x != expected.cells.x ||
      patch.cells.y != expected.cells.y || patch.cells.z != expected.cells.z ||
      patch.process_grid.x != expected.process_grid.x ||
      patch.process_grid.y != expected.process_grid.y ||
      patch.process_grid.z != expected.process_grid.z ||
      patch.process_coord.x != expected.process_coord.x ||
      patch.process_coord.y != expected.process_coord.y ||
      patch.process_coord.z != expected.process_coord.z) {
    return {StatusCode::invalid_plan, detail::halo_detail_topology};
  }
  const HaloTopology expected_topology = implementation_->topology;
  if (topology.periodic_x != expected_topology.periodic_x ||
      topology.periodic_y != expected_topology.periodic_y ||
      topology.periodic_z != expected_topology.periodic_z) {
    return {StatusCode::invalid_plan, detail::halo_detail_topology};
  }
  for (std::size_t index = 0U; index < fields.size; ++index) {
    const HaloFieldSpec actual = fields.data[index];
    const HaloFieldSpec expected_field = implementation_->fields[index].spec;
    if (actual.field != expected_field.field ||
        actual.width != expected_field.width ||
        actual.components != expected_field.components) {
      return {StatusCode::invalid_plan, detail::halo_detail_field};
    }
  }
  return {};
}

RevisionToken HaloEngine::ghost_revision(FieldId field) const noexcept {
  if (implementation_ == nullptr) {
    return RevisionToken{0U};
  }
  const std::int32_t slot =
      implementation_->field_lookup[static_cast<std::size_t>(field)];
  return slot < 0
             ? RevisionToken{0U}
             : implementation_->fields[static_cast<std::size_t>(slot)]
                   .ghost_revision;
}

int HaloEngine::lowest_failing_rank() const noexcept {
  return implementation_ == nullptr ? -1
                                    : implementation_->lowest_failing_rank;
}

bool HaloEngine::active() const noexcept {
  return implementation_ != nullptr && implementation_->exchange_active;
}

bool HaloEngine::ready() const noexcept {
  return implementation_ != nullptr && !implementation_->poisoned;
}

HaloPlanStats HaloEngine::plan_stats() const noexcept {
  return implementation_ == nullptr ? HaloPlanStats{} : implementation_->stats;
}

HaloRuntimeCounters HaloEngine::runtime_counters() const noexcept {
  return implementation_ == nullptr ? HaloRuntimeCounters{}
                                    : implementation_->counters;
}

}  // namespace hundun::v04
