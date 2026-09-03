// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/diag_performance_correctness.hpp"
#include "hundun/mesh_topology.hpp"
#include "hundun/rt_exchange_plan.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "hundun/rt_structured_decomposition.hpp"
#include "yyjson.h"

#include <mpi.h>

#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

namespace {

using hundun::diagnostics::parse_performance_correctness;
using hundun::mesh::LocalCellId;
using hundun::mesh::MeshTopology;
using hundun::runtime::DecompositionOptions;
using hundun::runtime::ExchangePlan;
using hundun::runtime::Int3;
using hundun::runtime::MpiContext;
using hundun::runtime::MpiEnvironment;
using hundun::runtime::StructuredDecomposition;

// Independent one-scalar fixture accounting. The assembly owns sixteen
// cell-sized FP64 work arrays plus the accepted per-cell metadata (132
// bytes/cell). The Stage 3-capable flow scratch path has a 340-byte-per-cell
// transient peak, one six-byte-per-cell partition workspace per split axis,
// and a 24-byte rank-local control block. Besides the persisted linear
// SolveReport reductions, each measured step performs 48 flow reductions
// carrying 2360 logical FP64 payload bytes. The extra scratch and payload
// remain present when the additive Stage 3 features are disabled, so
// schema-v1 accounting records them explicitly.
constexpr std::uint64_t kAssemblyBytesPerOwnedCell = 132U;
constexpr std::uint64_t kBasePeakBytesPerOwnedCell = 340U;
constexpr std::uint64_t kPartitionPeakBytesPerOwnedCell = 6U;
constexpr std::uint64_t kPeakFixedBytesPerRank = 24U;
constexpr std::uint64_t kNonlinearFp64CallsPerMeasuredStep = 48U;
constexpr std::uint64_t kFp64PayloadBytesPerMeasuredStep = 2360U;

std::uint64_t add(std::uint64_t left, std::uint64_t right) {
  if (right > std::numeric_limits<std::uint64_t>::max() - left)
    throw std::runtime_error("Task 25 oracle addition overflow");
  return left + right;
}

std::uint64_t multiply(std::uint64_t left, std::uint64_t right) {
  if (left != 0U &&
      right > std::numeric_limits<std::uint64_t>::max() / left)
    throw std::runtime_error("Task 25 oracle multiplication overflow");
  return left * right;
}

std::string read_file(const std::filesystem::path& path) {
  std::ifstream stream(path, std::ios::binary);
  if (!stream)
    throw std::runtime_error("unable to read Task 25 artifact");
  return {std::istreambuf_iterator<char>(stream),
          std::istreambuf_iterator<char>()};
}

yyjson_val* required(yyjson_val* object, const char* key) {
  auto* value = yyjson_obj_get(object, key);
  if (value == nullptr)
    throw std::runtime_error(std::string("missing artifact member: ") + key);
  return value;
}

template <std::size_t Size>
void require_exact_keys(yyjson_val* object,
                        const std::array<const char*, Size>& keys) {
  if (!yyjson_is_obj(object) || yyjson_obj_size(object) != Size)
    throw std::runtime_error("Task 25 exact-counter key set differs");
  for (const char* key : keys) {
    if (yyjson_obj_get(object, key) == nullptr)
      throw std::runtime_error("Task 25 exact-counter key set differs");
  }
}

std::uint64_t u64(yyjson_val* object, const char* key) {
  auto* value = required(object, key);
  if (!yyjson_is_uint(value))
    throw std::runtime_error(std::string("artifact member is not u64: ") +
                             key);
  return yyjson_get_uint(value);
}

std::array<std::uint64_t, 3> triplet(yyjson_val* object, const char* key) {
  auto* value = required(object, key);
  if (!yyjson_is_arr(value) || yyjson_arr_size(value) != 3U)
    throw std::runtime_error("artifact extent/grid shape differs");
  std::array<std::uint64_t, 3> result{};
  for (std::size_t index = 0U; index < result.size(); ++index) {
    auto* item = yyjson_arr_get(value, index);
    if (!yyjson_is_uint(item))
      throw std::runtime_error("artifact extent/grid value differs");
    result[index] = yyjson_get_uint(item);
  }
  return result;
}

std::uint64_t box_values(hundun::runtime::Box3 box) {
  const auto x =
      static_cast<std::uint64_t>(box.end.x - box.begin.x);
  const auto y =
      static_cast<std::uint64_t>(box.end.y - box.begin.y);
  const auto z =
      static_cast<std::uint64_t>(box.end.z - box.begin.z);
  return multiply(multiply(x, y), z);
}

int owner_coordinate(int cell, int cells, int partitions) {
  const int base = cells / partitions;
  const int remainder = cells % partitions;
  const int large_span = (base + 1) * remainder;
  return cell < large_span
             ? cell / (base + 1)
             : remainder + (cell - large_span) / base;
}

struct Partition final {
  std::uint64_t send_values{};
  std::uint64_t receive_values{};
  std::uint64_t send_messages{};
  std::uint64_t receive_messages{};
};

Partition pressure_partition(const MpiContext& mpi,
                             const StructuredDecomposition& decomposition,
                             const MeshTopology& topology) {
  std::vector<std::uint64_t> requests(
      static_cast<std::size_t>(mpi.size()));
  std::vector<std::uint64_t> incoming(
      static_cast<std::size_t>(mpi.size()));
  const auto global = decomposition.global_extent();
  const auto grid = decomposition.process_grid();
  for (LocalCellId local = topology.owned_cell_count();
       local < topology.local_cell_count(); ++local) {
    const auto cell = topology.global_cell(local);
    std::array<int, 3> coordinate{
        owner_coordinate(cell.x, global.x, grid.x),
        owner_coordinate(cell.y, global.y, grid.y),
        owner_coordinate(cell.z, global.z, grid.z)};
    int owner = -1;
    if (MPI_Cart_rank(decomposition.comm(), coordinate.data(), &owner) !=
            MPI_SUCCESS ||
        owner < 0 || owner >= mpi.size())
      throw std::runtime_error("unable to derive pressure partition owner");
    ++requests[static_cast<std::size_t>(owner)];
  }
  if (MPI_Alltoall(requests.data(), 1, MPI_UINT64_T, incoming.data(), 1,
                   MPI_UINT64_T, decomposition.comm()) != MPI_SUCCESS)
    throw std::runtime_error("unable to exchange pressure partition counts");
  Partition result;
  for (int rank = 0; rank < mpi.size(); ++rank) {
    const auto index = static_cast<std::size_t>(rank);
    result.receive_values = add(result.receive_values, requests[index]);
    result.send_values = add(result.send_values, incoming[index]);
    if (requests[index] != 0U)
      ++result.receive_messages;
    if (incoming[index] != 0U)
      ++result.send_messages;
  }
  return result;
}

std::uint64_t rank_sum(const MpiContext& mpi, std::uint64_t local) {
  std::uint64_t global = 0U;
  if (MPI_Allreduce(&local, &global, 1, MPI_UINT64_T, MPI_SUM,
                    mpi.comm()) != MPI_SUCCESS)
    throw std::runtime_error("unable to reduce Task 25 oracle value");
  return global;
}

struct Exact final {
  std::uint64_t allocated{};
  std::uint64_t peak{};
  std::uint64_t halo_send_bytes{};
  std::uint64_t halo_receive_bytes{};
  std::uint64_t halo_pack_bytes{};
  std::uint64_t halo_unpack_bytes{};
  std::uint64_t halo_send_messages{};
  std::uint64_t halo_receive_messages{};
  std::uint64_t fp64_calls{};
  std::uint64_t fp64_bytes{};
  std::uint64_t linear_reductions{};
  std::uint64_t momentum_matvec{};
  std::uint64_t pressure_matvec{};
  std::uint64_t momentum_preconditioner{};
  std::uint64_t pressure_preconditioner{};
};

bool exact_equal(const Exact& left, const Exact& right) noexcept {
  return left.allocated == right.allocated &&
         left.peak == right.peak &&
         left.halo_send_bytes == right.halo_send_bytes &&
         left.halo_receive_bytes == right.halo_receive_bytes &&
         left.halo_pack_bytes == right.halo_pack_bytes &&
         left.halo_unpack_bytes == right.halo_unpack_bytes &&
         left.halo_send_messages == right.halo_send_messages &&
         left.halo_receive_messages == right.halo_receive_messages &&
         left.fp64_calls == right.fp64_calls &&
         left.fp64_bytes == right.fp64_bytes &&
         left.linear_reductions == right.linear_reductions &&
         left.momentum_matvec == right.momentum_matvec &&
         left.pressure_matvec == right.pressure_matvec &&
         left.momentum_preconditioner == right.momentum_preconditioner &&
         left.pressure_preconditioner == right.pressure_preconditioner;
}

struct WorkTotals final {
  std::uint64_t linear_reductions{};
  std::uint64_t momentum_matvec{};
  std::uint64_t pressure_matvec{};
  std::uint64_t momentum_preconditioner{};
  std::uint64_t pressure_preconditioner{};
  std::uint64_t pressure_applies{};
};

struct ArtifactInput final {
  std::uint64_t measured_steps{};
  std::uint64_t repetitions{};
  Exact actual;
  WorkTotals work;
};
static_assert(std::is_trivially_copyable_v<ArtifactInput>);

std::array<std::uint64_t, 2> combined_message_counts(
    std::uint64_t runtime_send, std::uint64_t runtime_receive,
    std::uint64_t pressure_send, std::uint64_t pressure_receive) {
  return {add(runtime_send, pressure_send),
          add(runtime_receive, pressure_receive)};
}

Exact expected_exact(const MpiContext& mpi,
                     const StructuredDecomposition& decomposition,
                     const MeshTopology& topology,
                     const ExchangePlan& exchange_plan,
                     std::uint64_t scalar_count,
                     std::uint64_t measured_steps,
                     std::uint64_t repetitions,
                     const WorkTotals& work) {
  if (scalar_count != 1U)
    throw std::runtime_error(
        "Task 25 frozen exact-counter fixture requires one scalar");
  std::uint64_t send_region_values = 0U;
  std::uint64_t receive_region_values = 0U;
  std::uint64_t active_regions = 0U;
  for (const auto& region : exchange_plan.regions()) {
    if (region.neighbor_rank == MPI_PROC_NULL)
      continue;
    send_region_values =
        add(send_region_values, box_values(region.send_box));
    receive_region_values =
        add(receive_region_values, box_values(region.receive_box));
    ++active_regions;
  }

  // These are fixture formulas, intentionally isolated from the product
  // driver. The scalar parameter makes the field-exchange contribution
  // explicit: the accepted one-scalar fixture has 32 calls and 98 scalar
  // components per measured step.
  const auto runtime_calls_per_step = add(30U, multiply(2U, scalar_count));
  const auto runtime_components_per_step =
      add(96U, multiply(2U, scalar_count));
  const auto measured_repetitions =
      multiply(measured_steps, repetitions);
  const auto runtime_calls =
      multiply(measured_repetitions, runtime_calls_per_step);
  const auto runtime_components =
      multiply(measured_repetitions, runtime_components_per_step);
  const auto runtime_send_bytes = multiply(
      multiply(send_region_values, runtime_components), sizeof(double));
  const auto runtime_receive_bytes = multiply(
      multiply(receive_region_values, runtime_components), sizeof(double));
  const auto runtime_send_messages =
      multiply(active_regions, runtime_calls);
  const auto runtime_receive_messages =
      multiply(active_regions, runtime_calls);

  Exact result;
  result.linear_reductions = work.linear_reductions;
  result.momentum_matvec = work.momentum_matvec;
  result.pressure_matvec = work.pressure_matvec;
  result.momentum_preconditioner = work.momentum_preconditioner;
  result.pressure_preconditioner = work.pressure_preconditioner;
  const auto partition =
      pressure_partition(mpi, decomposition, topology);
  const auto pressure_send_bytes = multiply(
      multiply(partition.send_values, work.pressure_applies), sizeof(double));
  const auto pressure_receive_bytes = multiply(
      multiply(partition.receive_values, work.pressure_applies),
      sizeof(double));
  const auto pressure_send_messages =
      multiply(partition.send_messages, work.pressure_applies);
  const auto pressure_receive_messages =
      multiply(partition.receive_messages, work.pressure_applies);

  const auto local_send_bytes =
      add(runtime_send_bytes, pressure_send_bytes);
  const auto local_receive_bytes =
      add(runtime_receive_bytes, pressure_receive_bytes);
  const auto local_messages = combined_message_counts(
      runtime_send_messages, runtime_receive_messages,
      pressure_send_messages, pressure_receive_messages);
  const auto global_send_bytes = rank_sum(mpi, local_send_bytes);
  const auto global_receive_bytes = rank_sum(mpi, local_receive_bytes);
  const auto global_send_messages = rank_sum(mpi, local_messages[0U]);
  const auto global_receive_messages = rank_sum(mpi, local_messages[1U]);
  if (global_send_bytes != global_receive_bytes ||
      global_send_messages != global_receive_messages)
    throw std::runtime_error("fixture Halo oracle is not symmetric");
  result.halo_send_bytes = global_send_bytes;
  result.halo_receive_bytes = global_receive_bytes;
  result.halo_pack_bytes = global_send_bytes;
  result.halo_unpack_bytes = global_receive_bytes;
  result.halo_send_messages = global_send_messages;
  result.halo_receive_messages = global_receive_messages;

  const auto local = decomposition.local_extent();
  const auto local_cells = multiply(
      multiply(static_cast<std::uint64_t>(local.x),
               static_cast<std::uint64_t>(local.y)),
      static_cast<std::uint64_t>(local.z));
  result.allocated =
      multiply(
          rank_sum(mpi,
                   multiply(local_cells, kAssemblyBytesPerOwnedCell)),
          repetitions);
  const auto grid = decomposition.process_grid();
  const std::uint64_t partitioned_axes =
      static_cast<std::uint64_t>(grid.x > 1) +
      static_cast<std::uint64_t>(grid.y > 1) +
      static_cast<std::uint64_t>(grid.z > 1);
  const auto peak_bytes_per_cell =
      add(kBasePeakBytesPerOwnedCell,
          multiply(kPartitionPeakBytesPerOwnedCell, partitioned_axes));
  result.peak = rank_sum(
      mpi, add(multiply(local_cells, peak_bytes_per_cell),
               kPeakFixedBytesPerRank));

  result.fp64_calls = add(
      result.linear_reductions,
      multiply(multiply(kNonlinearFp64CallsPerMeasuredStep, measured_steps),
               repetitions));
  result.fp64_bytes =
      multiply(multiply(kFp64PayloadBytesPerMeasuredStep, measured_steps),
               repetitions);
  return result;
}

ArtifactInput parse_artifact(
    std::string_view artifact, int ranks,
    const std::array<std::uint64_t, 3>& global,
    const std::array<std::uint64_t, 3>& grid) {
  yyjson_doc* document = yyjson_read(artifact.data(), artifact.size(), 0);
  if (document == nullptr)
    throw std::runtime_error("invalid Task 25 artifact JSON");
  try {
    auto* root = yyjson_doc_get_root(document);
    if (u64(root, "schema_version") != 1U ||
        u64(root, "ranks") != static_cast<std::uint64_t>(ranks) ||
        triplet(root, "global_owned_cell_extents") != global ||
        triplet(root, "process_grid") != grid)
      throw std::runtime_error("Task 25 artifact fixture binding differs");
    auto* measurement = required(root, "measurement");
    const auto warmup_steps = u64(measurement, "warmup_steps");
    const auto measured_steps = u64(measurement, "measured_steps");
    const auto repetitions = u64(measurement, "repetitions");
    auto* correctness_json = required(root, "correctness");
    auto* passed = required(correctness_json, "passed");
    auto* summary = required(correctness_json, "summary");
    if (!yyjson_is_bool(passed) || !yyjson_get_bool(passed) ||
        !yyjson_is_str(summary))
      throw std::runtime_error("Task 25 artifact correctness differs");
    const auto correctness =
        parse_performance_correctness(yyjson_get_str(summary));
    hundun::diagnostics::validate_performance_correctness_coverage(
        correctness, warmup_steps, measured_steps, repetitions);
    ArtifactInput result;
    result.measured_steps = measured_steps;
    result.repetitions = repetitions;
    for (const auto& item : correctness.work) {
      if (item.phase != 'M')
        continue;
      result.work.linear_reductions =
          add(result.work.linear_reductions, item.reduction);
      if (item.slot < 3U) {
        result.work.momentum_matvec =
            add(result.work.momentum_matvec, item.matvec);
        result.work.momentum_preconditioner =
            add(result.work.momentum_preconditioner, item.preconditioner);
      } else {
        result.work.pressure_matvec =
            add(result.work.pressure_matvec, item.matvec);
        result.work.pressure_preconditioner =
            add(result.work.pressure_preconditioner, item.preconditioner);
        result.work.pressure_applies =
            add(result.work.pressure_applies, add(item.matvec, 1U));
      }
    }
    auto* exact = required(root, "exact_counters");
    require_exact_keys(
        exact,
        std::array<const char*, 8>{
            "allocated_bytes", "halo_payload_bytes", "halo_messages",
            "collectives", "collective_logical_payload_bytes", "matvec",
            "preconditioner_applications", "logical_io_bytes"});
    auto* allocated = required(exact, "allocated_bytes");
    auto* payload = required(exact, "halo_payload_bytes");
    auto* messages = required(exact, "halo_messages");
    auto* collectives = required(exact, "collectives");
    auto* collective_bytes =
        required(exact, "collective_logical_payload_bytes");
    auto* matvec = required(exact, "matvec");
    auto* preconditioner =
        required(exact, "preconditioner_applications");
    auto* io = required(exact, "logical_io_bytes");
    require_exact_keys(
        allocated,
        std::array<const char*, 2>{"execution.allocated",
                                   "execution.peak-live"});
    require_exact_keys(
        payload,
        std::array<const char*, 4>{"pack", "receive", "send", "unpack"});
    require_exact_keys(messages,
                       std::array<const char*, 2>{"receive", "send"});
    require_exact_keys(
        collectives,
        std::array<const char*, 3>{"checkpoint", "fp64-reduction",
                                   "linear-reduction"});
    require_exact_keys(
        collective_bytes,
        std::array<const char*, 1>{"fp64-reduction"});
    require_exact_keys(matvec,
                       std::array<const char*, 2>{"momentum", "pressure"});
    require_exact_keys(
        preconditioner,
        std::array<const char*, 2>{"momentum", "pressure"});
    require_exact_keys(
        io, std::array<const char*, 2>{"checkpoint", "diagnostics"});
    result.actual.allocated = u64(allocated, "execution.allocated");
    result.actual.peak = u64(allocated, "execution.peak-live");
    result.actual.halo_send_bytes = u64(payload, "send");
    result.actual.halo_receive_bytes = u64(payload, "receive");
    result.actual.halo_pack_bytes = u64(payload, "pack");
    result.actual.halo_unpack_bytes = u64(payload, "unpack");
    result.actual.halo_send_messages = u64(messages, "send");
    result.actual.halo_receive_messages = u64(messages, "receive");
    if (u64(collectives, "checkpoint") != 0U ||
        u64(io, "checkpoint") != 0U || u64(io, "diagnostics") != 0U)
      throw std::runtime_error("portable artifact I/O counters are nonzero");
    result.actual.fp64_calls = u64(collectives, "fp64-reduction");
    result.actual.linear_reductions =
        u64(collectives, "linear-reduction");
    result.actual.fp64_bytes =
        u64(collective_bytes, "fp64-reduction");
    result.actual.momentum_matvec = u64(matvec, "momentum");
    result.actual.pressure_matvec = u64(matvec, "pressure");
    result.actual.momentum_preconditioner =
        u64(preconditioner, "momentum");
    result.actual.pressure_preconditioner =
        u64(preconditioner, "pressure");
    yyjson_doc_free(document);
    return result;
  } catch (...) {
    yyjson_doc_free(document);
    throw;
  }
}

void validate_artifact(const MpiContext& mpi,
                       const StructuredDecomposition& decomposition,
                       const MeshTopology& topology,
                       const ExchangePlan& exchange_plan,
                       const ArtifactInput& input,
                       std::uint64_t scalar_count) {
  const auto expected =
      expected_exact(mpi, decomposition, topology, exchange_plan,
                     scalar_count, input.measured_steps, input.repetitions,
                     input.work);
  if (!exact_equal(input.actual, expected))
    throw std::runtime_error(
        "Task 25 artifact differs from independent exact-counter oracle");
}

int parse_positive(const char* text) {
  std::size_t consumed = 0U;
  const int value = std::stoi(text, &consumed);
  if (consumed != std::string_view(text).size() || value <= 0)
    throw std::runtime_error("Task 25 oracle argument is invalid");
  return value;
}

}  // namespace

int main(int argc, char** argv) {
  MpiEnvironment environment(argc, argv);
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  try {
    if (argc == 2 && std::string_view(argv[1]) == "--self-test") {
      const auto asymmetric =
          combined_message_counts(3U, 5U, 7U, 11U);
      if (asymmetric != std::array<std::uint64_t, 2>{10U, 16U})
        throw std::runtime_error(
            "Task 25 asymmetric message oracle self-test failed");
      return 0;
    }
    auto mpi = MpiContext::duplicate(MPI_COMM_WORLD);
    if (argc != 9)
      throw std::runtime_error(
          "usage: task25_performance_oracle artifact gx gy gz px py pz scalars");
    const Int3 global{parse_positive(argv[2]), parse_positive(argv[3]),
                      parse_positive(argv[4])};
    const Int3 grid_value{parse_positive(argv[5]), parse_positive(argv[6]),
                          parse_positive(argv[7])};
    const auto scalar_count =
        static_cast<std::uint64_t>(parse_positive(argv[8]));
    const auto grid_product =
        multiply(multiply(static_cast<std::uint64_t>(grid_value.x),
                          static_cast<std::uint64_t>(grid_value.y)),
                 static_cast<std::uint64_t>(grid_value.z));
    if (grid_product != static_cast<std::uint64_t>(mpi.size()))
      throw std::runtime_error("Task 25 oracle process grid differs");
    const std::array<std::uint64_t, 3> global_array{
        static_cast<std::uint64_t>(global.x),
        static_cast<std::uint64_t>(global.y),
        static_cast<std::uint64_t>(global.z)};
    const std::array<std::uint64_t, 3> grid_array{
        static_cast<std::uint64_t>(grid_value.x),
        static_cast<std::uint64_t>(grid_value.y),
        static_cast<std::uint64_t>(grid_value.z)};
    ArtifactInput input;
    int artifact_ok = 1;
    if (mpi.rank() == 0) {
      try {
        input = parse_artifact(read_file(argv[1]), mpi.size(),
                               global_array, grid_array);
      } catch (const std::exception&) {
        artifact_ok = 0;
      }
    }
    if (MPI_Bcast(&artifact_ok, 1, MPI_INT, 0, mpi.comm()) != MPI_SUCCESS)
      throw std::runtime_error(
          "unable to broadcast Task 25 artifact status");
    if (artifact_ok == 0)
      throw std::runtime_error(
          "unable to read or parse Task 25 artifact on rank zero");
    static_assert(sizeof(ArtifactInput) <=
                  static_cast<std::size_t>(
                      std::numeric_limits<int>::max()));
    if (MPI_Bcast(&input, static_cast<int>(sizeof(input)), MPI_BYTE, 0,
                  mpi.comm()) != MPI_SUCCESS)
      throw std::runtime_error(
          "unable to broadcast Task 25 artifact input");
    const std::array<bool, 3> periodic{true, true, true};
    auto decomposition = StructuredDecomposition::create(
        mpi, global, periodic, DecompositionOptions{grid_value});
    MeshTopology topology(decomposition);
    const auto exchange_plan = ExchangePlan::create(
        decomposition, decomposition.local_extent(), 2);
    validate_artifact(mpi, decomposition, topology, exchange_plan, input,
                      scalar_count);
    return 0;
  } catch (const std::exception& error) {
    if (rank == 0)
      std::cerr << error.what() << '\n';
    return 1;
  }
}
