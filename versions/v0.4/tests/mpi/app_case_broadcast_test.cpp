// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_case.hpp"

#include "app_case_detail.hpp"

#include <mpi.h>

#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <string_view>

namespace {

namespace fs = std::filesystem;
using hundun::v04::CaseCompiler;
using hundun::v04::Status;
using hundun::v04::StatusCode;
using hundun::v04::ValidatedModel;

int g_open_count = 0;
int g_observed_rank = -1;

void observe_open(int rank) {
  ++g_open_count;
  g_observed_rank = rank;
}

void reset_observer() {
  g_open_count = 0;
  g_observed_rank = -1;
  hundun::v04::detail::set_file_open_observer_for_test(observe_open);
}

bool expect(bool condition, int rank, std::string_view message) {
  if (!condition) {
    std::cerr << "rank " << rank << " FAIL: " << message << '\n';
  }
  return condition;
}

void write_file(const fs::path& path, std::string_view content) {
  std::ofstream output(path, std::ios::binary);
  output.write(content.data(), static_cast<std::streamsize>(content.size()));
}

bool same_u64(std::uint64_t value, MPI_Comm communicator) {
  std::uint64_t minimum = 0U;
  std::uint64_t maximum = 0U;
  return MPI_Allreduce(&value, &minimum, 1, MPI_UINT64_T, MPI_MIN,
                       communicator) == MPI_SUCCESS &&
         MPI_Allreduce(&value, &maximum, 1, MPI_UINT64_T, MPI_MAX,
                       communicator) == MPI_SUCCESS &&
         minimum == maximum;
}

bool same_status(Status status, MPI_Comm communicator) {
  const std::uint64_t packed =
      (static_cast<std::uint64_t>(status.code) << 32U) | status.detail;
  return same_u64(packed, communicator);
}

bool collective_result(bool local, MPI_Comm communicator) {
  const int value = local ? 1 : 0;
  int result = 0;
  return MPI_Allreduce(&value, &result, 1, MPI_INT, MPI_MIN, communicator) ==
             MPI_SUCCESS &&
         result != 0;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    return 2;
  }

  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  fs::path root;
  if (rank == 0) {
    const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
    root = fs::temp_directory_path() /
           ("hundun-v04-case-broadcast-" + std::to_string(stamp));
    std::error_code error;
    fs::remove_all(root, error);
    fs::create_directories(root);
    write_file(root / "profile.d", "0 1\n1 2\n");
    write_file(root / "shape.stl", "solid shape\nendsolid shape\n");
    write_file(root / "case.json", R"json({
      "schema_version": 1,
      "units": "SI",
      "mesh": {
        "kind": "tensor_stretched",
        "domain": {"lower": [0, 0, 0], "upper": [2, 1, 1]},
        "exact_cells": [16, 10, 8],
        "base_spacing": [0.25, 0.2, 0.2],
        "minimum_spacing": [0.05, 0.05, 0.05],
        "max_growth_ratio": 1.2,
        "focus_regions": [
          {"lower": [-1, 0.2, 0.2], "upper": [0.5, 0.8, 0.8],
           "target_spacing": [0.1, 0.1, 0.1]},
          {"lower": [-2, 0.2, 0.2], "upper": [0.5, 0.8, 0.8],
           "target_spacing": [0.1, 0.1, 0.1]}
        ],
        "limits": {
          "max_global_cells": 100000,
          "max_memory_bytes_per_rank": 134217728
        },
        "data_files": ["profile.d"],
        "stl_file": "shape.stl"
      },
      "flow": {
        "model": "single_phase_low_mach_compressible",
        "pressure_closure": "local_absolute_pressure_drho_dp",
        "reacting": false
      },
      "solver": {"coupling": "PISO", "pressure_correctors": 2},
      "turbulence": {"model": "vreman_wall_function"},
      "time": {"control": "adaptive_flow"}
    })json");
  }
  MPI_Barrier(MPI_COMM_WORLD);

  const fs::path nonroot_poison =
      "/nonroot-must-not-open-or-canonicalize-this-case";
  reset_observer();
  ValidatedModel first;
  const Status first_status = CaseCompiler::load_and_compile(
      MPI_COMM_WORLD, rank == 0 ? root : nonroot_poison, first);

  bool passed = true;
  passed &= expect(static_cast<bool>(first_status), rank,
                   "valid root case compiles collectively");
  passed &= expect(same_status(first_status, MPI_COMM_WORLD), rank,
                   "success status is identical");
  passed &= expect(same_u64(first.fingerprint, MPI_COMM_WORLD), rank,
                   "fingerprint is identical");
  const std::uint64_t model_signature =
      (static_cast<std::uint64_t>(first.mesh.kind) << 24U) |
      (static_cast<std::uint64_t>(first.turbulence) << 16U) |
      (static_cast<std::uint64_t>(first.time_control) << 8U) |
      static_cast<std::uint64_t>(first.data_files.size()) |
      (first.stl_file.has_value() ? (1ULL << 40U) : 0U);
  passed &= expect(same_u64(model_signature, MPI_COMM_WORLD), rank,
                   "typed model is identical");
  passed &= expect(first.mesh.kind == hundun::v04::GeometryKind::tensor_stretched &&
                       first.mesh.lower.x == 0.0 && first.mesh.upper.x == 2.0 &&
                       first.mesh.has_exact_cells &&
                       first.mesh.exact_cells.x == 16 &&
                       first.mesh.exact_cells.y == 10 &&
                       first.mesh.exact_cells.z == 8 &&
                       first.mesh.has_base_spacing &&
                       first.mesh.base_spacing.x == 0.25 &&
                       first.mesh.minimum_spacing.x == 0.05 &&
                       first.mesh.max_growth_ratio == 1.2 &&
                       first.mesh.limits.max_global_cells == 100000 &&
                       first.mesh.limits.max_memory_bytes_per_rank == 134217728,
                   rank, "wire publishes the complete typed mesh");
  passed &= expect(first.mesh.focus_regions.size() == 1U &&
                       first.mesh.focus_regions[0].lower.x == 0.0 &&
                       first.mesh.focus_regions[0].upper.x == 0.5 &&
                       first.mesh.focus_regions[0].target_spacing.x == 0.1,
                   rank, "wire publishes canonical clipped focus regions");
  passed &= expect(first.data_files.size() == 1U &&
                       first.data_files[0] == "profile.d" &&
                       first.stl_file == fs::path("shape.stl"),
                   rank, "relative references survive the wire");
  passed &= expect(rank == 0 ? g_open_count == 3 : g_open_count == 0, rank,
                   "only rank zero opens JSON/data/STL");
  passed &= expect(rank == 0 ? g_observed_rank == 0 : g_observed_rank == -1,
                   rank, "observer records only rank zero");

  const std::uint64_t first_fingerprint = first.fingerprint;
  if (rank == 0) {
    write_file(root / "profile.d", "0 1\n1 3\n");
  }
  MPI_Barrier(MPI_COMM_WORLD);
  reset_observer();
  ValidatedModel mutated;
  const Status mutated_status = CaseCompiler::load_and_compile(
      MPI_COMM_WORLD, rank == 0 ? root : nonroot_poison, mutated);
  passed &= expect(static_cast<bool>(mutated_status), rank,
                   "mutated referenced data recompiles");
  passed &= expect(same_u64(mutated.fingerprint, MPI_COMM_WORLD), rank,
                   "mutated fingerprint is identical");
  passed &= expect(mutated.fingerprint != first_fingerprint, rank,
                   "referenced bytes affect fingerprint");
  passed &= expect(rank == 0 ? g_open_count == 3 : g_open_count == 0, rank,
                   "recompile remains root-only I/O");

  if (rank == 0) {
    write_file(root / "case.json", "{ invalid json");
  }
  MPI_Barrier(MPI_COMM_WORLD);
  reset_observer();
  ValidatedModel rejected;
  rejected.fingerprint = 777U;
  const Status rejected_status = CaseCompiler::load_and_compile(
      MPI_COMM_WORLD, rank == 0 ? root : nonroot_poison, rejected);
  passed &= expect(rejected_status.code == StatusCode::invalid_case, rank,
                   "invalid root JSON fails collectively");
  passed &= expect(same_status(rejected_status, MPI_COMM_WORLD), rank,
                   "failure category and detail are identical");
  passed &= expect(rejected.fingerprint == 777U, rank,
                   "failed compile does not publish a partial model");
  passed &= expect(
      hundun::v04::detail::last_lowest_failing_rank_for_test() == 0, rank,
      "lowest failing rank is broadcast");
  passed &= expect(rank == 0 ? g_open_count == 1 : g_open_count == 0, rank,
                   "failure parsing remains root-only I/O");

  passed = collective_result(passed, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0) {
    std::error_code error;
    fs::remove_all(root, error);
  }
  hundun::v04::detail::set_file_open_observer_for_test(nullptr);

  if (rank == 0 && !passed) {
    std::cerr << "v04 app case broadcast failed for " << size << " ranks\n";
  }
  MPI_Finalize();
  return passed ? 0 : 1;
}
