// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "hundun/v04_status.hpp"

#include <mpi.h>

#include <cstdint>
#include <filesystem>
#include <optional>
#include <vector>

namespace hundun::v04 {

enum class GeometryKind : std::uint8_t { uniform, tensor_stretched };
enum class TurbulenceKind : std::uint8_t {
  none,
  wale,
  vreman_wall_function
};
enum class TimeControlKind : std::uint8_t {
  fixed,
  adaptive_flow,
  adaptive_acoustic
};

struct CaseSpec {
  std::filesystem::path root;
};

struct FocusRegionSpec {
  Real3 lower{};
  Real3 upper{};
  Real3 target_spacing{};
};

struct MeshLimits {
  std::uint64_t max_global_cells{};
  std::uint64_t max_memory_bytes_per_rank{};
};

struct CartesianMeshSpec {
  GeometryKind kind{GeometryKind::uniform};
  Real3 lower{};
  Real3 upper{};
  bool has_exact_cells{};
  Int3 exact_cells{};
  bool has_base_spacing{};
  Real3 base_spacing{};
  Real3 minimum_spacing{};
  double max_growth_ratio{1.0};
  std::vector<FocusRegionSpec> focus_regions;
  MeshLimits limits{};
};

struct ValidatedModel {
  CartesianMeshSpec mesh;
  TurbulenceKind turbulence{TurbulenceKind::vreman_wall_function};
  TimeControlKind time_control{TimeControlKind::fixed};
  std::vector<std::filesystem::path> data_files;
  std::optional<std::filesystem::path> stl_file;
  PlanFingerprint fingerprint{};
};

class CaseCompiler {
 public:
  static Status load_and_compile(MPI_Comm communicator,
                                 const std::filesystem::path& case_root,
                                 ValidatedModel& out);
};

}  // namespace hundun::v04
