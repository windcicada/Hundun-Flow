// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_app.hpp"
#include "hundun/v04_case.hpp"
#include "hundun/v04_physics.hpp"

#include "app_evidence_detail.hpp"
#include "app_identity_detail.hpp"

#include <mpi.h>

#include <fcntl.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <array>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

namespace fs = std::filesystem;
using namespace hundun::v04;

constexpr std::uint32_t kRunnerInput = 10901U;
constexpr std::string_view kSpecMagic =
    "HUNDUN_V04_LITERATURE_STATISTICS_V1";
constexpr std::string_view kRunMagic =
    "HUNDUN_V04_THIN_DOMAIN_RUN_V1";
constexpr std::string_view kStatisticsMagic =
    "HUNDUN_V04_THIN_DOMAIN_STATISTICS_V1";
constexpr std::string_view kAccumulatorMagic =
    "HUNDUN_V04_THIN_DOMAIN_ACCUMULATOR_V1";
constexpr std::string_view kCheckpointMagic =
    "HUNDUN_V04_THIN_DOMAIN_CHECKPOINT_V1";

static_assert(kRuntimePressureEnergyRefinementCapacity ==
              kPressureEnergyRefinementCapacity);

struct Options {
  fs::path spec;
  fs::path case_root;
  fs::path run_root;
  fs::path restart_root;
  std::uint64_t steps{};
  std::uint64_t visit_interval{};
  bool have_steps{};
  bool have_visit_interval{};
  bool dry_plan{};
  bool self_test{};
};

struct StatisticsSpec {
  std::uint64_t development_steps{};
  std::uint64_t collection_end_step{};
  std::uint64_t checkpoint_interval{};
  double rho_ref{};
  double u_ref{};
  double diameter{};
  double span{};
  double cylinder_center_x{};
  std::vector<double> station_x_over_d;
};

struct Bracket {
  std::int32_t lower{};
  std::int32_t upper{};
  double upper_weight{};
};

struct RuntimeGeometry {
  std::vector<Bracket> station_brackets;
  Bracket centerline_bracket{};
};

struct Accumulator {
  std::vector<double> profile;
  std::vector<double> centerline;
  std::uint64_t sample_steps{};
};

struct ThermoExtrema {
  double pressure_min{};
  double pressure_max{};
  double temperature_min{};
  double temperature_max{};
  double density_min{};
  double density_max{};
};

struct RestartBinding {
  fs::path statistics;
  fs::path accumulator;
  std::string generation;
};

bool parse_u64(std::string_view text, std::uint64_t& out) {
  if (text.empty() || text.front() == '-') return false;
  errno = 0;
  char* end = nullptr;
  const std::string owned(text);
  const unsigned long long value = std::strtoull(owned.c_str(), &end, 10);
  if (errno != 0 || end != owned.c_str() + owned.size()) return false;
  out = static_cast<std::uint64_t>(value);
  return true;
}

bool parse_real(std::string_view text, double& out) {
  if (text.empty()) return false;
  errno = 0;
  char* end = nullptr;
  const std::string owned(text);
  const double value = std::strtod(owned.c_str(), &end);
  if (errno != 0 || end != owned.c_str() + owned.size() ||
      !std::isfinite(value))
    return false;
  out = value;
  return true;
}

bool parse_options(int argc, char** argv, Options& out) {
  for (int index = 1; index < argc; ++index) {
    const std::string_view token(argv[index]);
    if (token == "--dry-plan") {
      if (out.dry_plan) return false;
      out.dry_plan = true;
    } else if (token == "--self-test") {
      if (out.self_test) return false;
      out.self_test = true;
    } else {
      if (index + 1 >= argc) return false;
      const std::string_view value(argv[++index]);
      if (token == "--spec" && out.spec.empty()) {
        out.spec = value;
      } else if (token == "--case-root" && out.case_root.empty()) {
        out.case_root = value;
      } else if (token == "--run-root" && out.run_root.empty()) {
        out.run_root = value;
      } else if (token == "--restart-root" && out.restart_root.empty()) {
        out.restart_root = value;
      } else if (token == "--steps" && !out.have_steps &&
                 parse_u64(value, out.steps) && out.steps != 0U) {
        out.have_steps = true;
      } else if (token == "--visit-interval" &&
                 !out.have_visit_interval &&
                 parse_u64(value, out.visit_interval)) {
        out.have_visit_interval = true;
      } else {
        return false;
      }
    }
  }
  if (out.self_test)
    return !out.dry_plan && out.spec.empty() && out.case_root.empty() &&
           out.run_root.empty() && out.restart_root.empty() &&
           !out.have_steps && !out.have_visit_interval;
  if (out.spec.empty() || out.case_root.empty()) return false;
  if (out.dry_plan)
    return out.run_root.empty() && out.restart_root.empty() &&
           !out.have_steps;
  return !out.run_root.empty() && out.have_steps;
}

bool parse_spec(const fs::path& path, StatisticsSpec& out,
                std::string& error) {
  std::ifstream input(path);
  if (!input) {
    error = "cannot_open_spec";
    return false;
  }
  std::string line;
  if (!std::getline(input, line) || line != kSpecMagic) {
    error = "invalid_spec_magic";
    return false;
  }
  StatisticsSpec candidate;
  bool development = false;
  bool collection_end = false;
  bool checkpoint = false;
  bool rho = false;
  bool velocity = false;
  bool diameter = false;
  bool span = false;
  bool center = false;
  bool ended = false;
  while (std::getline(input, line)) {
    if (line.empty()) continue;
    if (ended) {
      error = "data_after_end";
      return false;
    }
    std::istringstream row(line);
    std::string key;
    std::string value;
    std::string extra;
    if (!(row >> key)) continue;
    if (key == "end") {
      if (row >> extra) {
        error = "invalid_end";
        return false;
      }
      ended = true;
      continue;
    }
    if (!(row >> value) || row >> extra) {
      error = "invalid_spec_row";
      return false;
    }
    if (key == "development_steps" && !development) {
      development = parse_u64(value, candidate.development_steps);
      if (!development) error = "invalid_development_steps";
    } else if (key == "collection_end_step" && !collection_end) {
      collection_end = parse_u64(value, candidate.collection_end_step);
      if (!collection_end) error = "invalid_collection_end_step";
    } else if (key == "checkpoint_interval" && !checkpoint) {
      checkpoint = parse_u64(value, candidate.checkpoint_interval) &&
                   candidate.checkpoint_interval != 0U;
      if (!checkpoint) error = "invalid_checkpoint_interval";
    } else if (key == "rho_ref" && !rho) {
      rho = parse_real(value, candidate.rho_ref) && candidate.rho_ref > 0.0;
      if (!rho) error = "invalid_rho_ref";
    } else if (key == "u_ref" && !velocity) {
      velocity = parse_real(value, candidate.u_ref) && candidate.u_ref > 0.0;
      if (!velocity) error = "invalid_u_ref";
    } else if (key == "diameter" && !diameter) {
      diameter = parse_real(value, candidate.diameter) &&
                 candidate.diameter > 0.0;
      if (!diameter) error = "invalid_diameter";
    } else if (key == "span" && !span) {
      span = parse_real(value, candidate.span) && candidate.span > 0.0;
      if (!span) error = "invalid_span";
    } else if (key == "cylinder_center_x" && !center) {
      center = parse_real(value, candidate.cylinder_center_x);
      if (!center) error = "invalid_cylinder_center_x";
    } else if (key == "station_x_over_d") {
      double station = 0.0;
      if (!parse_real(value, station)) {
        error = "invalid_station";
        return false;
      }
      candidate.station_x_over_d.push_back(station);
    } else {
      error = "unknown_or_duplicate_key";
      return false;
    }
    if (!error.empty()) return false;
  }
  if (!ended || !development || !collection_end || !checkpoint || !rho ||
      !velocity || !diameter || !span || !center ||
      candidate.station_x_over_d.empty() ||
      candidate.collection_end_step <= candidate.development_steps) {
    error = "incomplete_or_inconsistent_spec";
    return false;
  }
  for (std::size_t index = 1U; index < candidate.station_x_over_d.size();
       ++index) {
    if (!(candidate.station_x_over_d[index] >
          candidate.station_x_over_d[index - 1U])) {
      error = "stations_not_strictly_increasing";
      return false;
    }
  }
  out = std::move(candidate);
  return true;
}

std::uint64_t mix(std::uint64_t hash, std::uint64_t value) {
  hash ^= value + UINT64_C(0x9e3779b97f4a7c15) + (hash << 6U) +
          (hash >> 2U);
  return hash;
}

std::uint64_t bits(double value) {
  std::uint64_t result = 0U;
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

std::uint64_t spec_fingerprint(const StatisticsSpec& spec) {
  std::uint64_t hash = UINT64_C(0x5448494e53544154);
  hash = mix(hash, spec.development_steps);
  hash = mix(hash, spec.collection_end_step);
  hash = mix(hash, spec.checkpoint_interval);
  hash = mix(hash, bits(spec.rho_ref));
  hash = mix(hash, bits(spec.u_ref));
  hash = mix(hash, bits(spec.diameter));
  hash = mix(hash, bits(spec.span));
  hash = mix(hash, bits(spec.cylinder_center_x));
  hash = mix(hash, spec.station_x_over_d.size());
  for (double station : spec.station_x_over_d)
    hash = mix(hash, bits(station));
  return hash == 0U ? 1U : hash;
}

bool all_true(MPI_Comm communicator, bool local) {
  int value = local ? 1 : 0;
  return MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_INT, MPI_MIN,
                       communicator) == MPI_SUCCESS &&
         value == 1;
}

bool consensus_u64(MPI_Comm communicator, std::uint64_t value) {
  std::uint64_t minimum = value;
  std::uint64_t maximum = value;
  return MPI_Allreduce(MPI_IN_PLACE, &minimum, 1, MPI_UINT64_T, MPI_MIN,
                       communicator) == MPI_SUCCESS &&
         MPI_Allreduce(MPI_IN_PLACE, &maximum, 1, MPI_UINT64_T, MPI_MAX,
                       communicator) == MPI_SUCCESS &&
         minimum == maximum;
}

bool write_exclusive(const fs::path& path, std::string_view text) {
  const int descriptor =
      ::open(path.c_str(), O_WRONLY | O_CREAT | O_EXCL | O_CLOEXEC, 0444);
  if (descriptor < 0) return false;
  std::size_t cursor = 0U;
  bool okay = true;
  while (cursor < text.size()) {
    const ssize_t count =
        ::write(descriptor, text.data() + cursor, text.size() - cursor);
    if (count <= 0) {
      okay = false;
      break;
    }
    cursor += static_cast<std::size_t>(count);
  }
  if (okay && ::fsync(descriptor) != 0) okay = false;
  if (::close(descriptor) != 0) okay = false;
  return okay;
}

bool sync_directory(const fs::path& path) {
  const int descriptor =
      ::open(path.c_str(), O_RDONLY | O_DIRECTORY | O_CLOEXEC);
  if (descriptor < 0) return false;
  const bool okay = ::fsync(descriptor) == 0;
  return ::close(descriptor) == 0 && okay;
}

std::string checkpoint_name(std::uint64_t step, std::string_view suffix) {
  std::ostringstream name;
  name << "step-" << std::setw(20) << std::setfill('0') << step << suffix;
  return name.str();
}

bool bracket(Span<const double> coordinates, double target, Bracket& out) {
  if (coordinates.data == nullptr || coordinates.size == 0U ||
      !std::isfinite(target) || target < coordinates.data[0] ||
      target > coordinates.data[coordinates.size - 1U])
    return false;
  const double* begin = coordinates.data;
  const double* end = begin + coordinates.size;
  const double* found = std::lower_bound(begin, end, target);
  if (found != end && *found == target) {
    const auto index = static_cast<std::int32_t>(found - begin);
    out = {index, index, 0.0};
    return true;
  }
  if (found == begin || found == end) return false;
  const auto upper = static_cast<std::int32_t>(found - begin);
  const auto lower = upper - 1;
  const double width = coordinates.data[upper] - coordinates.data[lower];
  if (!(width > 0.0)) return false;
  const double weight = (target - coordinates.data[lower]) / width;
  if (!std::isfinite(weight) || weight < 0.0 || weight > 1.0) return false;
  out = {lower, upper, weight};
  return true;
}

const SnapshotFieldView* find_field(const CommittedOutputSnapshot& snapshot,
                                    std::string_view name,
                                    std::uint8_t components) {
  const SnapshotFieldView* result = nullptr;
  for (std::size_t index = 0U; index < snapshot.fields.size; ++index) {
    const SnapshotFieldView& field = snapshot.fields.data[index];
    if (field.stable_name == name) {
      if (result != nullptr || field.values.components != components)
        return nullptr;
      result = &field;
    }
  }
  return result;
}

bool prepare_geometry(const StatisticsSpec& spec,
                      const CommittedOutputSnapshot& snapshot,
                      RuntimeGeometry& out) {
  if (!snapshot.committed || snapshot.geometry == nullptr) return false;
  RuntimeGeometry candidate;
  const Span<const double> x = snapshot.geometry->x().centres();
  const Span<const double> y = snapshot.geometry->y().centres();
  candidate.station_brackets.resize(spec.station_x_over_d.size());
  for (std::size_t index = 0U; index < spec.station_x_over_d.size(); ++index) {
    const double target = spec.cylinder_center_x +
                          spec.station_x_over_d[index] * spec.diameter;
    if (!bracket(x, target, candidate.station_brackets[index])) return false;
  }
  if (!bracket(y, 0.0, candidate.centerline_bracket)) return false;
  out = std::move(candidate);
  return true;
}

bool owns(std::int32_t global, std::int32_t begin, std::int32_t cells) {
  return global >= begin && global < begin + cells;
}

bool add_sample(const StatisticsSpec& spec, const RuntimeGeometry& runtime,
                const CommittedOutputSnapshot& snapshot,
                Accumulator& accumulator) {
  const SnapshotFieldView* velocity = find_field(snapshot, "U", 3U);
  if (velocity == nullptr || snapshot.geometry == nullptr) return false;
  const ConstFieldView u = velocity->values;
  const MeshPatch patch = snapshot.patch;
  const Int3 global = snapshot.geometry->global_cells();
  const std::size_t expected_profile =
      spec.station_x_over_d.size() * static_cast<std::size_t>(global.y) * 6U;
  const std::size_t expected_center = static_cast<std::size_t>(global.x) * 3U;
  if (u.interior.x != patch.cells.x || u.interior.y != patch.cells.y ||
      u.interior.z != patch.cells.z || global.x <= 0 || global.y <= 0 ||
      global.z <= 0 || accumulator.profile.size() != expected_profile ||
      accumulator.centerline.size() != expected_center)
    return false;

  for (std::size_t station = 0U; station < runtime.station_brackets.size();
       ++station) {
    const Bracket selected = runtime.station_brackets[station];
    const std::array<std::int32_t, 2U> planes{{selected.lower,
                                               selected.upper}};
    const std::array<double, 2U> weights{{1.0 - selected.upper_weight,
                                          selected.upper_weight}};
    const std::size_t plane_count = selected.lower == selected.upper ? 1U : 2U;
    for (std::size_t side = 0U; side < plane_count; ++side) {
      const std::int32_t gx = planes[side];
      const double weight = weights[side];
      if (weight == 0.0 || !owns(gx, patch.begin.x, patch.cells.x)) continue;
      const std::int32_t lx = gx - patch.begin.x;
      for (std::int32_t lz = 0; lz < patch.cells.z; ++lz) {
        for (std::int32_t ly = 0; ly < patch.cells.y; ++ly) {
          const double ux = u.unchecked({lx, ly, lz}, 0U);
          const double uy = u.unchecked({lx, ly, lz}, 1U);
          if (!std::isfinite(ux) || !std::isfinite(uy)) return false;
          const std::size_t gy =
              static_cast<std::size_t>(patch.begin.y + ly);
          const std::size_t base =
              (station * static_cast<std::size_t>(global.y) + gy) * 6U;
          accumulator.profile[base] += weight;
          accumulator.profile[base + 1U] += weight * ux;
          accumulator.profile[base + 2U] += weight * uy;
          accumulator.profile[base + 3U] += weight * ux * ux;
          accumulator.profile[base + 4U] += weight * uy * uy;
          accumulator.profile[base + 5U] += weight * ux * uy;
        }
      }
    }
  }

  const Bracket selected = runtime.centerline_bracket;
  const std::array<std::int32_t, 2U> planes{{selected.lower, selected.upper}};
  const std::array<double, 2U> weights{{1.0 - selected.upper_weight,
                                        selected.upper_weight}};
  const std::size_t plane_count = selected.lower == selected.upper ? 1U : 2U;
  for (std::size_t side = 0U; side < plane_count; ++side) {
    const std::int32_t gy = planes[side];
    const double weight = weights[side];
    if (weight == 0.0 || !owns(gy, snapshot.patch.begin.y,
                               snapshot.patch.cells.y))
      continue;
    const std::int32_t ly = gy - snapshot.patch.begin.y;
    for (std::int32_t lz = 0; lz < snapshot.patch.cells.z; ++lz) {
      for (std::int32_t lx = 0; lx < snapshot.patch.cells.x; ++lx) {
        const double ux = u.unchecked({lx, ly, lz}, 0U);
        if (!std::isfinite(ux)) return false;
        const std::size_t gx =
            static_cast<std::size_t>(snapshot.patch.begin.x + lx);
        const std::size_t base = gx * 3U;
        accumulator.centerline[base] += weight;
        accumulator.centerline[base + 1U] += weight * ux;
        accumulator.centerline[base + 2U] += weight * ux * ux;
      }
    }
  }
  ++accumulator.sample_steps;
  return true;
}

bool collect_probes(MPI_Comm communicator, const RuntimeGeometry& runtime,
                    const CommittedOutputSnapshot& snapshot,
                    std::vector<double>& out) {
  const SnapshotFieldView* velocity = find_field(snapshot, "U", 3U);
  if (!all_true(communicator, velocity != nullptr)) return false;
  const ConstFieldView u = velocity->values;
  std::vector<double> sums(runtime.station_brackets.size() * 4U, 0.0);
  const Bracket y_selected = runtime.centerline_bracket;
  const std::array<std::int32_t, 2U> y_planes{{y_selected.lower,
                                                y_selected.upper}};
  const std::array<double, 2U> y_weights{{1.0 - y_selected.upper_weight,
                                          y_selected.upper_weight}};
  const std::size_t y_count =
      y_selected.lower == y_selected.upper ? 1U : 2U;
  bool local_valid = true;
  for (std::size_t station = 0U; station < runtime.station_brackets.size();
       ++station) {
    const Bracket x_selected = runtime.station_brackets[station];
    const std::array<std::int32_t, 2U> x_planes{{x_selected.lower,
                                                  x_selected.upper}};
    const std::array<double, 2U> x_weights{{1.0 - x_selected.upper_weight,
                                            x_selected.upper_weight}};
    const std::size_t x_count =
        x_selected.lower == x_selected.upper ? 1U : 2U;
    for (std::size_t xs = 0U; xs < x_count; ++xs) {
      if (!owns(x_planes[xs], snapshot.patch.begin.x,
                snapshot.patch.cells.x))
        continue;
      const std::int32_t lx = x_planes[xs] - snapshot.patch.begin.x;
      for (std::size_t ys = 0U; ys < y_count; ++ys) {
        if (!owns(y_planes[ys], snapshot.patch.begin.y,
                  snapshot.patch.cells.y))
          continue;
        const std::int32_t ly = y_planes[ys] - snapshot.patch.begin.y;
        const double weight = x_weights[xs] * y_weights[ys];
        if (weight == 0.0) continue;
        for (std::int32_t lz = 0; lz < snapshot.patch.cells.z; ++lz) {
          const Int3 cell{lx, ly, lz};
          const double ux = u.unchecked(cell, 0U);
          const double uy = u.unchecked(cell, 1U);
          const double uz = u.unchecked(cell, 2U);
          if (!std::isfinite(ux) || !std::isfinite(uy) ||
              !std::isfinite(uz)) {
            local_valid = false;
            continue;
          }
          const std::size_t base = station * 4U;
          sums[base] += weight;
          sums[base + 1U] += weight * ux;
          sums[base + 2U] += weight * uy;
          sums[base + 3U] += weight * uz;
        }
      }
    }
  }
  if (!all_true(communicator, local_valid)) return false;
  if (sums.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
      MPI_Allreduce(MPI_IN_PLACE, sums.data(), static_cast<int>(sums.size()),
                    MPI_DOUBLE, MPI_SUM, communicator) != MPI_SUCCESS)
    return false;
  out.resize(runtime.station_brackets.size() * 3U);
  for (std::size_t station = 0U; station < runtime.station_brackets.size();
       ++station) {
    const std::size_t source = station * 4U;
    const std::size_t target = station * 3U;
    if (!std::isfinite(sums[source]) || sums[source] <= 0.0) return false;
    for (std::size_t component = 0U; component < 3U; ++component) {
      const double value = sums[source + component + 1U] / sums[source];
      if (!std::isfinite(value)) return false;
      out[target + component] = value;
    }
  }
  return true;
}

bool reduce_vector(MPI_Comm communicator, int rank,
                   const std::vector<double>& local,
                   std::vector<double>& global) {
  if (local.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    return false;
  if (rank == 0) global.resize(local.size());
  return MPI_Reduce(local.data(), rank == 0 ? global.data() : nullptr,
                    static_cast<int>(local.size()), MPI_DOUBLE, MPI_SUM, 0,
                    communicator) == MPI_SUCCESS;
}

std::string encode_accumulator(std::uint64_t fingerprint,
                               PlanFingerprint plan,
                               PlanFingerprint schema,
                               std::uint64_t step, double time,
                               const Accumulator& accumulator) {
  std::ostringstream text;
  text.imbue(std::locale::classic());
  text << std::setprecision(17) << kAccumulatorMagic << '\n'
       << "spec_fingerprint " << fingerprint << '\n'
       << "plan " << plan << '\n'
       << "schema " << schema << '\n'
       << "step " << step << '\n'
       << "time " << time << '\n'
       << "sample_steps " << accumulator.sample_steps << '\n'
       << "profile_size " << accumulator.profile.size() << '\n'
       << "centerline_size " << accumulator.centerline.size() << '\n';
  for (std::size_t index = 0U; index < accumulator.profile.size(); ++index)
    text << "p " << index << ' ' << accumulator.profile[index] << '\n';
  for (std::size_t index = 0U; index < accumulator.centerline.size(); ++index)
    text << "c " << index << ' ' << accumulator.centerline[index] << '\n';
  text << "end\n";
  return text.str();
}

bool decode_accumulator(const fs::path& path, std::uint64_t fingerprint,
                        PlanFingerprint plan, PlanFingerprint schema,
                        std::uint64_t step, double time,
                        Accumulator& accumulator) {
  std::ifstream input(path);
  std::string line;
  if (!input || !std::getline(input, line) || line != kAccumulatorMagic)
    return false;
  bool have_spec = false;
  bool have_plan = false;
  bool have_schema = false;
  bool have_step = false;
  bool have_time = false;
  bool have_samples = false;
  bool have_profile_size = false;
  bool have_centerline_size = false;
  bool ended = false;
  std::vector<bool> profile_seen(accumulator.profile.size(), false);
  std::vector<bool> centerline_seen(accumulator.centerline.size(), false);
  while (std::getline(input, line)) {
    if (line.empty()) continue;
    std::istringstream row(line);
    row.imbue(std::locale::classic());
    std::string key;
    if (!(row >> key) || ended) return false;
    if (key == "end") {
      std::string extra;
      if (row >> extra) return false;
      ended = true;
      continue;
    }
    if (key == "p" || key == "c") {
      std::size_t index = 0U;
      double value = 0.0;
      std::string extra;
      if (!(row >> index >> value) || row >> extra || !std::isfinite(value))
        return false;
      std::vector<double>& target =
          key == "p" ? accumulator.profile : accumulator.centerline;
      std::vector<bool>& seen = key == "p" ? profile_seen : centerline_seen;
      if (index >= target.size() || seen[index]) return false;
      target[index] = value;
      seen[index] = true;
      continue;
    }
    std::string value;
    std::string extra;
    if (!(row >> value) || row >> extra) return false;
    std::uint64_t number = 0U;
    if (key == "spec_fingerprint" && !have_spec) {
      have_spec = parse_u64(value, number) && number == fingerprint;
      if (!have_spec) return false;
    } else if (key == "plan" && !have_plan) {
      have_plan = parse_u64(value, number) && number == plan;
      if (!have_plan) return false;
    } else if (key == "schema" && !have_schema) {
      have_schema = parse_u64(value, number) && number == schema;
      if (!have_schema) return false;
    } else if (key == "step" && !have_step) {
      have_step = parse_u64(value, number) && number == step;
      if (!have_step) return false;
    } else if (key == "time" && !have_time) {
      double parsed = 0.0;
      have_time = parse_real(value, parsed) && parsed == time;
      if (!have_time) return false;
    } else if (key == "sample_steps" && !have_samples) {
      have_samples = parse_u64(value, accumulator.sample_steps);
      if (!have_samples) return false;
    } else if (key == "profile_size" && !have_profile_size) {
      have_profile_size = parse_u64(value, number) &&
                          number == accumulator.profile.size();
      if (!have_profile_size) return false;
    } else if (key == "centerline_size" && !have_centerline_size) {
      have_centerline_size = parse_u64(value, number) &&
                             number == accumulator.centerline.size();
      if (!have_centerline_size) return false;
    } else {
      return false;
    }
  }
  return ended && have_spec && have_plan && have_schema && have_step &&
         have_time && have_samples && have_profile_size &&
         have_centerline_size &&
         std::all_of(profile_seen.begin(), profile_seen.end(),
                     [](bool value) { return value; }) &&
         std::all_of(centerline_seen.begin(), centerline_seen.end(),
                     [](bool value) { return value; });
}

std::string encode_statistics(const StatisticsSpec& spec,
                              std::uint64_t fingerprint,
                              const RuntimeGeometry& runtime,
                              const CommittedOutputSnapshot& snapshot,
                              const Accumulator& accumulator) {
  std::ostringstream text;
  text.imbue(std::locale::classic());
  text << std::setprecision(17)
       << "{\"schema\":\"" << kStatisticsMagic << "\""
       << ",\"spec_fingerprint\":" << fingerprint
       << ",\"plan\":" << snapshot.plan
       << ",\"field_schema\":" << snapshot.schema
       << ",\"snapshot_step\":" << snapshot.step
       << ",\"snapshot_time\":" << snapshot.time
       << ",\"sample_steps\":" << accumulator.sample_steps
       << ",\"profiles\":[";
  const Span<const double> xs = snapshot.geometry->x().centres();
  const Span<const double> ys = snapshot.geometry->y().centres();
  bool first = true;
  for (std::size_t station = 0U; station < runtime.station_brackets.size();
       ++station) {
    const Bracket selected = runtime.station_brackets[station];
    const double target = spec.cylinder_center_x +
                          spec.station_x_over_d[station] * spec.diameter;
    for (std::size_t y = 0U; y < ys.size; ++y) {
      const std::size_t base = (station * ys.size + y) * 6U;
      if (!first) text << ',';
      first = false;
      text << "{\"station\":" << station
           << ",\"x_over_d\":" << spec.station_x_over_d[station]
           << ",\"target_x\":" << target
           << ",\"lower_x\":" << xs.data[selected.lower]
           << ",\"upper_x\":" << xs.data[selected.upper]
           << ",\"upper_weight\":" << selected.upper_weight
           << ",\"y\":" << ys.data[y]
           << ",\"y_over_d\":" << ys.data[y] / spec.diameter
           << ",\"sum_weight\":" << accumulator.profile[base]
           << ",\"sum_u\":" << accumulator.profile[base + 1U]
           << ",\"sum_v\":" << accumulator.profile[base + 2U]
           << ",\"sum_uu\":" << accumulator.profile[base + 3U]
           << ",\"sum_vv\":" << accumulator.profile[base + 4U]
           << ",\"sum_uv\":" << accumulator.profile[base + 5U]
           << '}';
    }
  }
  const Bracket center = runtime.centerline_bracket;
  text << "],\"centerline_bracket\":{\"lower_y\":"
       << ys.data[center.lower] << ",\"upper_y\":"
       << ys.data[center.upper] << ",\"upper_weight\":"
       << center.upper_weight << "},\"centerline\":[";
  for (std::size_t x = 0U; x < xs.size; ++x) {
    const std::size_t base = x * 3U;
    if (x != 0U) text << ',';
    text << "{\"x\":" << xs.data[x] << ",\"x_over_d\":"
         << (xs.data[x] - spec.cylinder_center_x) / spec.diameter
         << ",\"sum_weight\":" << accumulator.centerline[base]
         << ",\"sum_u\":" << accumulator.centerline[base + 1U]
         << ",\"sum_uu\":" << accumulator.centerline[base + 2U]
         << '}';
  }
  text << "]}\n";
  return text.str();
}

bool read_current_generation(const fs::path& restart_root,
                             std::string& generation) {
  std::ifstream input(restart_root / "current");
  return static_cast<bool>(std::getline(input, generation)) &&
         !generation.empty() && generation.find('/') == std::string::npos &&
         generation.find("..") == std::string::npos;
}

bool read_restart_binding(const fs::path& restart_root,
                          const RestartImage& image,
                          std::uint64_t fingerprint,
                          RestartBinding& out) {
  std::string generation;
  if (!read_current_generation(restart_root, generation)) return false;
  const fs::path run_root = restart_root.parent_path();
  const fs::path marker =
      run_root / checkpoint_name(image.step, ".complete");
  std::ifstream input(marker);
  std::string line;
  if (!input || !std::getline(input, line) || line != kCheckpointMagic)
    return false;
  bool have_spec = false;
  bool have_step = false;
  bool have_time = false;
  bool have_plan = false;
  bool have_schema = false;
  bool have_generation = false;
  bool have_statistics = false;
  bool have_accumulator = false;
  bool ended = false;
  fs::path statistics;
  fs::path accumulator;
  while (std::getline(input, line)) {
    if (line.empty()) continue;
    std::istringstream row(line);
    std::string key;
    std::string value;
    std::string extra;
    if (!(row >> key) || ended) return false;
    if (key == "end") {
      if (row >> extra) return false;
      ended = true;
      continue;
    }
    if (!(row >> value) || row >> extra) return false;
    std::uint64_t number = 0U;
    if (key == "spec_fingerprint" && !have_spec) {
      have_spec = parse_u64(value, number) && number == fingerprint;
      if (!have_spec) return false;
    } else if (key == "step" && !have_step) {
      have_step = parse_u64(value, number) && number == image.step;
      if (!have_step) return false;
    } else if (key == "time" && !have_time) {
      double parsed = 0.0;
      have_time = parse_real(value, parsed) && parsed == image.time;
      if (!have_time) return false;
    } else if (key == "plan" && !have_plan) {
      have_plan = parse_u64(value, number) && number == image.plan;
      if (!have_plan) return false;
    } else if (key == "schema" && !have_schema) {
      have_schema = parse_u64(value, number) && number == image.schema;
      if (!have_schema) return false;
    } else if (key == "restart_generation" && !have_generation) {
      have_generation = value == generation;
      if (!have_generation) return false;
    } else if (key == "statistics" && !have_statistics) {
      statistics = run_root / value;
      have_statistics = statistics.filename() == value &&
                        fs::is_regular_file(statistics);
      if (!have_statistics) return false;
    } else if (key == "accumulator" && !have_accumulator) {
      accumulator = run_root / value;
      have_accumulator = accumulator.filename() == value &&
                         fs::is_regular_file(accumulator);
      if (!have_accumulator) return false;
    } else {
      return false;
    }
  }
  if (!ended || !have_spec || !have_step || !have_time || !have_plan ||
      !have_schema || !have_generation || !have_statistics ||
      !have_accumulator)
    return false;
  out = {statistics, accumulator, generation};
  return true;
}

bool write_checkpoint_observables(
    MPI_Comm communicator, int rank, const fs::path& run_root,
    const StatisticsSpec& spec, std::uint64_t fingerprint,
    const RuntimeGeometry& runtime,
    const CommittedOutputSnapshot& snapshot,
    const Accumulator& local_accumulator) {
  Accumulator global;
  global.sample_steps = local_accumulator.sample_steps;
  bool okay = reduce_vector(communicator, rank, local_accumulator.profile,
                            global.profile) &&
              reduce_vector(communicator, rank, local_accumulator.centerline,
                            global.centerline);
  if (!all_true(communicator, okay)) return false;
  if (rank == 0) {
    const fs::path statistics =
        run_root / checkpoint_name(snapshot.step, ".statistics.json");
    const fs::path accumulator =
        run_root / checkpoint_name(snapshot.step, ".accumulator");
    okay = write_exclusive(
                statistics,
                encode_statistics(spec, fingerprint, runtime, snapshot,
                                  global)) &&
            write_exclusive(accumulator,
                            encode_accumulator(fingerprint, snapshot.plan,
                                               snapshot.schema, snapshot.step,
                                               snapshot.time, global));
  }
  return all_true(communicator, okay);
}

bool checkpoint(MPI_Comm communicator, int rank, ProductDriver& driver,
                const fs::path& run_root, const StatisticsSpec& spec,
                std::uint64_t fingerprint, const RuntimeGeometry& runtime,
                const CommittedOutputSnapshot& snapshot,
                const Accumulator& accumulator) {
  if (!write_checkpoint_observables(communicator, rank, run_root, spec,
                                    fingerprint, runtime, snapshot,
                                    accumulator))
    return false;
  RestartSnapshot restart;
  Status status = driver.committed_restart_snapshot(restart);
  if (status)
    status = RestartWriter::write(communicator, run_root / "Restart", restart,
                                  {2U});
  if (!status) return false;
  bool okay = true;
  if (rank == 0) {
    std::string generation;
    okay = read_current_generation(run_root / "Restart", generation);
    if (okay) {
      std::ostringstream marker;
      marker.imbue(std::locale::classic());
      marker << std::setprecision(17) << kCheckpointMagic << '\n'
             << "spec_fingerprint " << fingerprint << '\n'
             << "step " << snapshot.step << '\n'
             << "time " << snapshot.time << '\n'
             << "plan " << snapshot.plan << '\n'
             << "schema " << snapshot.schema << '\n'
             << "restart_generation " << generation << '\n'
             << "statistics "
             << checkpoint_name(snapshot.step, ".statistics.json") << '\n'
             << "accumulator "
             << checkpoint_name(snapshot.step, ".accumulator") << '\n'
             << "end\n";
      okay = write_exclusive(
                 run_root / checkpoint_name(snapshot.step, ".complete"),
                 marker.str()) &&
             sync_directory(run_root);
    }
  }
  return all_true(communicator, okay);
}

DriverInitialState initial_state(const ValidatedModel& model,
                                 std::vector<double>& scalars) {
  scalars.assign(model.transported_scalars.size(), 0.0);
  DriverInitialState initial;
  initial.transported_scalars = {scalars.data(), scalars.size()};
  for (const BoundaryFaceSpec& boundary : model.boundaries) {
    if (std::isfinite(boundary.temperature) && boundary.temperature > 0.0)
      initial.temperature = boundary.temperature;
    if (model.pressure_reference == PressureReferenceKind::boundary_absolute &&
        boundary.flow_kind == BoundaryKind::pressure_outlet &&
        std::isfinite(boundary.pressure) && boundary.pressure > 0.0)
      initial.pressure_reference = boundary.pressure;
    if (boundary.flow_kind == BoundaryKind::velocity_inlet ||
        boundary.flow_kind == BoundaryKind::static_state_inlet ||
        boundary.flow_kind == BoundaryKind::total_state_inlet)
      initial.velocity = boundary.velocity;
  }
  return initial;
}

bool create_run_root(MPI_Comm communicator, int rank, const fs::path& root) {
  bool okay = true;
  if (rank == 0) {
    std::error_code error;
    okay = !fs::exists(root, error) && !error &&
            fs::create_directories(root, error) && !error &&
            sync_directory(root.parent_path());
  }
  return all_true(communicator, okay);
}

bool finite_force(const SurfaceForce& force) {
  const std::array<double, 12U> values{{
      force.pressure.x, force.pressure.y, force.pressure.z,
      force.viscous.x,  force.viscous.y,  force.viscous.z,
      force.total.x,    force.total.y,    force.total.z,
      force.moment.x,   force.moment.y,   force.moment.z}};
  return std::all_of(values.begin(), values.end(),
                     [](double value) { return std::isfinite(value); });
}

bool collect_thermo_extrema(MPI_Comm communicator,
                            const ValidatedModel& model,
                            const ThermodynamicsPlan& thermodynamics,
                            double pressure_reference,
                            const CommittedOutputSnapshot& snapshot,
                            ThermoExtrema& out) {
  const SnapshotFieldView* velocity = find_field(snapshot, "U", 3U);
  const SnapshotFieldView* pressure = find_field(snapshot, "pi", 1U);
  const SnapshotFieldView* enthalpy = find_field(snapshot, "h", 1U);
  bool local_valid =
      velocity != nullptr && pressure != nullptr && enthalpy != nullptr;
  std::vector<const SnapshotFieldView*> independent;
  for (const TransportedScalarSpec& scalar : model.transported_scalars) {
    if (scalar.role == TransportedScalarRole::species) {
      const SnapshotFieldView* field =
          find_field(snapshot, scalar.stable_name, 1U);
      if (field == nullptr) local_valid = false;
      independent.push_back(field);
    }
  }
  local_valid = local_valid &&
                independent.size() ==
                    thermodynamics.independent_species_count();
  if (!all_true(communicator, local_valid)) return false;
  const ConstFieldView u = velocity->values;
  const ConstFieldView pi = pressure->values;
  const ConstFieldView h = enthalpy->values;
  local_valid = u.interior.x == snapshot.patch.cells.x &&
                u.interior.y == snapshot.patch.cells.y &&
                u.interior.z == snapshot.patch.cells.z &&
                pi.interior.x == u.interior.x &&
                pi.interior.y == u.interior.y &&
                pi.interior.z == u.interior.z &&
                h.interior.x == u.interior.x &&
                h.interior.y == u.interior.y &&
                h.interior.z == u.interior.z;
  if (!all_true(communicator, local_valid)) return false;
  std::array<double, 3U> minimum{{
      std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::infinity()}};
  std::array<double, 3U> maximum{{
      -std::numeric_limits<double>::infinity(),
      -std::numeric_limits<double>::infinity(),
      -std::numeric_limits<double>::infinity()}};
  std::vector<double> fractions(independent.size(), 0.0);
  for (std::int32_t z = 0; z < u.interior.z; ++z) {
    for (std::int32_t y = 0; y < u.interior.y; ++y) {
      for (std::int32_t x = 0; x < u.interior.x; ++x) {
        const Int3 cell{x, y, z};
        for (std::size_t scalar = 0U; scalar < independent.size(); ++scalar)
          fractions[scalar] = independent[scalar]->values.unchecked(cell, 0U);
        const double perturbation = pi.unchecked(cell, 0U);
        const double absolute = pressure_reference + perturbation;
        const double enthalpy_value = h.unchecked(cell, 0U);
        const Real3 velocity_value{u.unchecked(cell, 0U),
                                   u.unchecked(cell, 1U),
                                   u.unchecked(cell, 2U)};
        ThermoState state;
        const Status status = thermodynamics.evaluate_from_reference_pressure(
            pressure_reference, perturbation, enthalpy_value,
            {fractions.data(), fractions.size()}, velocity_value, state);
        if (!status || !std::isfinite(absolute) || absolute <= 0.0 ||
            !std::isfinite(state.temperature) || state.temperature <= 0.0 ||
            !std::isfinite(state.rho) || state.rho <= 0.0) {
          local_valid = false;
          continue;
        }
        minimum[0U] = std::min(minimum[0U], absolute);
        minimum[1U] = std::min(minimum[1U], state.temperature);
        minimum[2U] = std::min(minimum[2U], state.rho);
        maximum[0U] = std::max(maximum[0U], absolute);
        maximum[1U] = std::max(maximum[1U], state.temperature);
        maximum[2U] = std::max(maximum[2U], state.rho);
      }
    }
  }
  if (!all_true(communicator, local_valid)) return false;
  if (MPI_Allreduce(MPI_IN_PLACE, minimum.data(), 3, MPI_DOUBLE, MPI_MIN,
                    communicator) != MPI_SUCCESS ||
      MPI_Allreduce(MPI_IN_PLACE, maximum.data(), 3, MPI_DOUBLE, MPI_MAX,
                    communicator) != MPI_SUCCESS)
    return false;
  const bool valid =
      std::all_of(minimum.begin(), minimum.end(), [](double value) {
        return std::isfinite(value) && value > 0.0;
      }) &&
      std::all_of(maximum.begin(), maximum.end(), [](double value) {
        return std::isfinite(value) && value > 0.0;
      });
  if (!valid) return false;
  out = {minimum[0U], maximum[0U], minimum[1U], maximum[1U], minimum[2U],
         maximum[2U]};
  return true;
}

RuntimePressureEnergyRefinementTermination refinement_termination(
    PressureEnergyRefinementTermination termination) {
  switch (termination) {
    case PressureEnergyRefinementTermination::component_residuals_converged:
      return RuntimePressureEnergyRefinementTermination::
          component_residuals_converged;
    case PressureEnergyRefinementTermination::iteration_capacity_exhausted:
      return RuntimePressureEnergyRefinementTermination::
          iteration_capacity_exhausted;
    case PressureEnergyRefinementTermination::rejected_candidate:
      return RuntimePressureEnergyRefinementTermination::rejected_candidate;
    case PressureEnergyRefinementTermination::none:
      return RuntimePressureEnergyRefinementTermination::none;
  }
  return RuntimePressureEnergyRefinementTermination::none;
}

Status collect_evidence_resources(MPI_Comm communicator,
                                  const DriverResourceReport& local,
                                  DriverResourceReport& global) {
  std::array<std::uint64_t, 16U> maxima{{
      local.structured_messages,
      local.structured_bytes,
      local.ibm_messages,
      local.ibm_bytes,
      local.reduction_collectives,
      local.reduction_nanoseconds,
      local.linear_iterations,
      local.exact_numeric_refills,
      local.hierarchy_rebuilds,
      local.preconditioner_applications,
      local.structured_exchanges,
      local.ibm_exchanges,
      local.reduction_logical_bytes,
      local.mg_blocking_collectives,
      local.mg_collective_logical_bytes,
      local.predictor_blocking_collectives,
  }};
  if (MPI_Allreduce(MPI_IN_PLACE, maxima.data(),
                    static_cast<int>(maxima.size()), MPI_UINT64_T, MPI_MAX,
                    communicator) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kRunnerInput};
  std::array<std::uint64_t, 4U> sums{{
      local.structured_messages, local.structured_bytes, local.ibm_messages,
      local.ibm_bytes}};
  if (MPI_Allreduce(MPI_IN_PLACE, sums.data(), static_cast<int>(sums.size()),
                    MPI_UINT64_T, MPI_SUM, communicator) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kRunnerInput};
  global = local;
  global.structured_messages = sums[0U];
  global.structured_bytes = sums[1U];
  global.ibm_messages = sums[2U];
  global.ibm_bytes = sums[3U];
  global.reduction_collectives = maxima[4U];
  global.reduction_nanoseconds = maxima[5U];
  global.linear_iterations = maxima[6U];
  global.exact_numeric_refills = maxima[7U];
  global.hierarchy_rebuilds = maxima[8U];
  global.preconditioner_applications = maxima[9U];
  global.structured_exchanges = maxima[10U];
  global.ibm_exchanges = maxima[11U];
  global.reduction_logical_bytes = maxima[12U];
  global.mg_blocking_collectives = maxima[13U];
  global.mg_collective_logical_bytes = maxima[14U];
  global.predictor_blocking_collectives = maxima[15U];
  return {};
}

Status collect_peak_rss(MPI_Comm communicator, std::uint64_t& maximum_rank,
                        std::uint64_t& maximum_node) {
  rusage usage{};
  if (::getrusage(RUSAGE_SELF, &usage) != 0 || usage.ru_maxrss < 0)
    return {StatusCode::io_failure, kRunnerInput};
  const auto kib = static_cast<std::uint64_t>(usage.ru_maxrss);
  if (kib > UINT64_MAX / 1024U)
    return {StatusCode::invalid_plan, kRunnerInput};
  const std::uint64_t local = kib * 1024U;
  if (MPI_Allreduce(&local, &maximum_rank, 1, MPI_UINT64_T, MPI_MAX,
                    communicator) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kRunnerInput};
  MPI_Comm shared = MPI_COMM_NULL;
  if (MPI_Comm_split_type(communicator, MPI_COMM_TYPE_SHARED, 0,
                          MPI_INFO_NULL, &shared) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kRunnerInput};
  std::uint64_t node = 0U;
  const int reduce = MPI_Allreduce(&local, &node, 1, MPI_UINT64_T, MPI_SUM,
                                   shared);
  const int free = MPI_Comm_free(&shared);
  if (reduce != MPI_SUCCESS || free != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kRunnerInput};
  if (MPI_Allreduce(&node, &maximum_node, 1, MPI_UINT64_T, MPI_MAX,
                    communicator) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kRunnerInput};
  return {};
}

Status collect_stage_timings(
    MPI_Comm communicator, const DriverStepReport& step,
    std::array<StageTimingRecord, kDriverTimedStageCapacity>& stages) {
  constexpr std::array<StageId, kDriverTimedStageCapacity> kStages{
      10U, 12U, 15U, 20U, 30U, 40U, 45U, 50U, 60U, 70U};
  std::array<std::uint64_t, kDriverTimedStageCapacity> local{};
  for (std::size_t index = 0U; index < step.stage_timing_count; ++index) {
    const DriverStageTiming timing = step.stage_timings[index];
    const auto found = std::find(kStages.begin(), kStages.end(), timing.stage);
    if (timing.stage == 0U || found == kStages.end())
      return {StatusCode::invalid_plan, kRunnerInput};
    const std::size_t target =
        static_cast<std::size_t>(found - kStages.begin());
    if (local[target] != 0U)
      return {StatusCode::invalid_plan, kRunnerInput};
    local[target] = timing.nanoseconds;
  }
  auto minimum = local;
  auto maximum = local;
  auto sum = local;
  if (MPI_Allreduce(MPI_IN_PLACE, minimum.data(),
                    static_cast<int>(minimum.size()), MPI_UINT64_T, MPI_MIN,
                    communicator) != MPI_SUCCESS ||
      MPI_Allreduce(MPI_IN_PLACE, maximum.data(),
                    static_cast<int>(maximum.size()), MPI_UINT64_T, MPI_MAX,
                    communicator) != MPI_SUCCESS ||
      MPI_Allreduce(MPI_IN_PLACE, sum.data(), static_cast<int>(sum.size()),
                    MPI_UINT64_T, MPI_SUM, communicator) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kRunnerInput};
  int ranks = 0;
  if (MPI_Comm_size(communicator, &ranks) != MPI_SUCCESS || ranks <= 0)
    return {StatusCode::mpi_failure, kRunnerInput};
  for (std::size_t index = 0U; index < stages.size(); ++index)
    stages[index] = {kStages[index], minimum[index],
                     sum[index] / static_cast<std::uint64_t>(ranks),
                     maximum[index]};
  return {};
}

Status append_evidence(
    MPI_Comm communicator, const fs::path& path,
    const IoServicePlan& services, const ValidatedModel& model,
    PlanFingerprint product, PlanFingerprint cpu, PlanFingerprint stl,
    const RuntimeCandidateIdentity& candidate_identity,
    PlanFingerprint build_identity, PlanFingerprint binary_identity,
    const RuntimeRunStartAnchor& run_start,
    const DriverStepReport& step, std::uint64_t maximum_step_nanoseconds,
    bool restart_recovery, DriverResourceReport& global_resources) {
  Status status =
      collect_evidence_resources(communicator, step.resources,
                                 global_resources);
  if (!status) return status;
  std::uint64_t maximum_rank_rss = 0U;
  std::uint64_t maximum_node_rss = 0U;
  status = collect_peak_rss(communicator, maximum_rank_rss,
                            maximum_node_rss);
  if (!status) return status;
  std::array<StageTimingRecord, kDriverTimedStageCapacity> stages{};
  status = collect_stage_timings(communicator, step, stages);
  if (!status) return status;

  RuntimeEvidenceRecord evidence;
  evidence.build = build_identity;
  evidence.binary = binary_identity;
  evidence.candidate_identity = candidate_identity;
  evidence.case_model = model.fingerprint;
  evidence.stl = stl;
  evidence.product = product;
  evidence.cpu_plan = cpu;
  evidence.run_start = run_start;
  evidence.step = step.accepted_step;
  evidence.previous_committed_time = step.proposal.time;
  evidence.time = step.accepted_time;
  evidence.requested_bdf_order = step.proposal.bdf.order;
  evidence.bdf_order = step.effective_bdf.order;
  evidence.thermophysical_predictor_calls =
      step.thermophysical_predictor_calls;
  evidence.temporal_method_fallback = step.temporal_method_fallback;
  evidence.maximum_rank_step_nanoseconds = maximum_step_nanoseconds;
  evidence.maximum_rank_rss_bytes = maximum_rank_rss;
  evidence.maximum_node_rss_bytes = maximum_node_rss;
  evidence.structured_messages = global_resources.structured_messages;
  evidence.structured_bytes = global_resources.structured_bytes;
  evidence.ibm_messages = global_resources.ibm_messages;
  evidence.ibm_bytes = global_resources.ibm_bytes;
  if (global_resources.mg_blocking_collectives >
      UINT64_MAX - global_resources.reduction_collectives)
    return {StatusCode::invalid_plan, kRunnerInput};
  evidence.blocking_collectives = global_resources.reduction_collectives +
                                  global_resources.mg_blocking_collectives;
  if (global_resources.predictor_blocking_collectives >
      UINT64_MAX - evidence.blocking_collectives)
    return {StatusCode::invalid_plan, kRunnerInput};
  evidence.blocking_collectives +=
      global_resources.predictor_blocking_collectives;
  evidence.reduction_nanoseconds = global_resources.reduction_nanoseconds;
  evidence.linear_iterations = global_resources.linear_iterations;
  evidence.exact_numeric_refills = global_resources.exact_numeric_refills;
  evidence.coarse_numeric_refills = global_resources.exact_numeric_refills;
  evidence.preconditioner_setups = global_resources.hierarchy_rebuilds;
  evidence.preconditioner_reuses =
      global_resources.exact_numeric_refills >=
              global_resources.hierarchy_rebuilds
          ? global_resources.exact_numeric_refills -
                global_resources.hierarchy_rebuilds
          : 0U;
  evidence.pressure = step.piso.pressure;
  evidence.pressure_solve_calls = step.piso.pressure_solve_calls;
  evidence.pressure_solve_contract =
      RuntimePressureSolveContract::continuity_energy_coupled;
  evidence.pressure_energy_refinement_solve_calls =
      step.piso.pressure_energy_refinement_solve_calls;
  evidence.pressure_energy_refinement_termination =
      refinement_termination(step.piso.pressure_energy_refinement_termination);
  for (std::uint8_t index = 0U;
       index < step.piso.pressure_energy_refinement_solve_calls &&
       index < evidence.pressure_energy_refinement.size();
       ++index) {
    const PisoPressureEnergyRefinementSolveReport& source =
        step.piso.pressure_energy_refinement[index];
    RuntimePressureEnergyRefinementSolve& destination =
        evidence.pressure_energy_refinement[index];
    destination.solve = source.solve;
    destination.target_generation = source.target_generation;
    destination.collective_lineage = source.collective_lineage;
    destination.ordinal = source.ordinal;
  }
  evidence.terminal_physical_audit.present =
      step.piso.final_flux_revision != 0U;
  evidence.terminal_physical_audit.final_flux_revision =
      step.piso.final_flux_revision;
  evidence.terminal_physical_audit.eos_residual = step.piso.eos_residual;
  evidence.terminal_physical_audit.eos_tolerance = model.solver.terminal.eos;
  evidence.terminal_physical_audit.continuity_residual =
      step.piso.continuity_residual;
  evidence.terminal_physical_audit.continuity_tolerance =
      model.solver.terminal.continuity;
  evidence.terminal_physical_audit.energy_residual =
      step.piso.energy_residual;
  evidence.terminal_physical_audit.energy_tolerance =
      model.solver.terminal.continuity;
  evidence.terminal_physical_audit.closed_mass_residual =
      step.piso.closed_mass_residual;
  evidence.terminal_physical_audit.closed_mass_tolerance =
      model.solver.terminal.closed_mass;
  evidence.terminal_physical_audit.gauge_residual =
      step.piso.gauge_residual;
  evidence.terminal_physical_audit.gauge_tolerance =
      model.solver.terminal.gauge;
  status = detail::runtime_committed_cfl(
      communicator, step.piso.committed_convective_cfl,
      evidence.committed_convective_cfl);
  if (!status) return status;
  evidence.momentum_predictor = step.momentum_predictor_solve.components;
  evidence.momentum_predictor_solve_calls =
      step.momentum_predictor_solve.solve_calls;
  evidence.predictor_enthalpy_endpoint =
      step.thermophysical_predictor.enthalpy_endpoint;
  evidence.predictor_enthalpy_endpoint_alpha =
      step.thermophysical_predictor.enthalpy_endpoint_alpha;
  evidence.predictor_bdf_endpoint_alpha =
      step.thermophysical_predictor.bdf_endpoint_alpha;
  evidence.predictor_source_endpoint_alpha =
      step.thermophysical_predictor.source_endpoint_alpha;
  evidence.predictor_enthalpy_solve_calls =
      step.thermophysical_predictor.enthalpy_solve_calls;
  evidence.momentum_predictor_theta = step.momentum_predictor_limiter.theta;
  evidence.momentum_predictor_activations =
      step.momentum_predictor_limiter.activations;
  evidence.momentum_correction_metrics_applicable =
      step.momentum_predictor_limiter.correction_metrics_applicable;
  evidence.momentum_minimum_face_alpha =
      step.momentum_predictor_limiter.minimum_face_alpha;
  evidence.momentum_active_correction_faces =
      step.momentum_predictor_limiter.active_correction_faces;
  evidence.momentum_limited_face_fraction =
      step.momentum_predictor_limiter.limited_face_fraction;
  evidence.momentum_predictor_limited =
      step.momentum_predictor_limiter.limited;
  status = detail::runtime_advective_cfl(
      communicator, step.momentum_predictor_limiter.advective_cfl,
      evidence.momentum_advective_cfl);
  if (!status) return status;
  evidence.predictor_theta = step.thermophysical_predictor.theta;
  evidence.predictor_mass_flux_scale =
      step.thermophysical_predictor.mass_flux_scale;
  evidence.predictor_low_margin = step.thermophysical_predictor.low_margin;
  evidence.predictor_high_margin = step.thermophysical_predictor.high_margin;
  evidence.predictor_low_order_transport_passes =
      step.thermophysical_predictor.low_order_transport_passes;
  evidence.predictor_low_order_substeps =
      step.thermophysical_predictor.low_order_substeps;
  evidence.predictor_low_order_halo_exchanges =
      step.thermophysical_predictor.low_order_halo_exchanges;
  evidence.predictor_blocking_collectives =
      step.thermophysical_predictor.blocking_collectives;
  evidence.predictor_limiting_cell_x =
      step.thermophysical_predictor.limiting_cell.x;
  evidence.predictor_limiting_cell_y =
      step.thermophysical_predictor.limiting_cell.y;
  evidence.predictor_limiting_cell_z =
      step.thermophysical_predictor.limiting_cell.z;
  evidence.predictor_limiting_rank =
      step.thermophysical_predictor.limiting_rank;
  evidence.predictor_constraint =
      static_cast<std::uint8_t>(step.thermophysical_predictor.constraint);
  evidence.predictor_low_state =
      static_cast<std::uint8_t>(step.thermophysical_predictor.low_state);
  evidence.predictor_limited = step.thermophysical_predictor.limited;
  evidence.stages = {stages.data(), stages.size()};
  evidence.startup = step.accepted_step == 1U;
  evidence.retry = step.attempts != 1U;
  evidence.restart_recovery = restart_recovery;
  // This quick thin-domain comparison is intentionally not a formal release
  // statistics candidate.  The full V6 numerical evidence is still emitted.
  evidence.statistics_eligible = false;
  return EvidenceWriter::append(communicator, path, services, evidence);
}

bool terminal_audit_valid(const ValidatedModel& model,
                          const DriverStepReport& step) {
  const double energy_tolerance = model.solver.terminal.continuity;
  const std::array<double, 10U> values{{
      step.piso.eos_residual,
      model.solver.terminal.eos,
      step.piso.continuity_residual,
      model.solver.terminal.continuity,
      step.piso.energy_residual,
      energy_tolerance,
      step.piso.closed_mass_residual,
      model.solver.terminal.closed_mass,
      step.piso.gauge_residual,
      model.solver.terminal.gauge}};
  return step.piso.final_flux_revision != 0U &&
         std::all_of(values.begin(), values.end(), [](double value) {
           return std::isfinite(value) && value >= 0.0;
         }) &&
         step.piso.eos_residual <= model.solver.terminal.eos &&
         step.piso.continuity_residual <=
             model.solver.terminal.continuity &&
         step.piso.energy_residual <= energy_tolerance &&
         step.piso.closed_mass_residual <=
             model.solver.terminal.closed_mass &&
         step.piso.gauge_residual <= model.solver.terminal.gauge &&
         step.momentum_predictor_limiter.advective_cfl.valid() &&
         step.piso.committed_convective_cfl.valid();
}

int run(MPI_Comm communicator, int rank, const Options& options) {
  StatisticsSpec spec;
  std::string parse_error;
  bool okay = parse_spec(options.spec, spec, parse_error);
  if (!all_true(communicator, okay)) {
    if (rank == 0) std::cerr << "spec_error=" << parse_error << '\n';
    return 3;
  }
  const std::uint64_t fingerprint = spec_fingerprint(spec);
  if (!consensus_u64(communicator, fingerprint)) {
    if (rank == 0) std::cerr << "spec_consensus_failure\n";
    return 3;
  }

  ValidatedModel model;
  Status status =
      CaseCompiler::load_and_compile(communicator, options.case_root, model);
  ThermodynamicsPlan thermodynamics;
  if (status)
    status = ThermodynamicsPlan::compile(model.thermophysics,
                                         {model.transported_scalars.data(),
                                          model.transported_scalars.size()},
                                         thermodynamics);
  CompiledCasePlan plan;
  if (status)
    status = ProductCompiler::compile(communicator, model, options.case_root,
                                      plan);
  if (!status) {
    if (rank == 0)
      std::cerr << "compile_status=" << static_cast<unsigned>(status.code)
                << '/' << status.detail << '\n';
    return 4;
  }
  const PlanFingerprint product_fingerprint = plan.fingerprint();
  const PlanFingerprint cpu_fingerprint = plan.cpu_plan_fingerprint();
  const PlanFingerprint stl_fingerprint = plan.stl_fingerprint();
  const PlanSummary summary = plan.summary();
  const IoServicePlan* sealed_services = plan.io_services();
  if (sealed_services == nullptr || sealed_services->fingerprint() == 0U)
    return 4;
  const IoServicePlan services = *sealed_services;
  ProductDriver driver;
  status = ProductDriver::create(communicator, std::move(plan), driver);
  if (!status) return 4;

  std::uint64_t starting_step = 0U;
  bool restarted = false;
  bool restart_requires_recovery = false;
  RuntimeRunStartAnchor run_start;
  RestartBinding restart_binding;
  if (!options.restart_root.empty()) {
    RestartExpected expected;
    status = driver.restart_expected(expected);
    RestartImage image;
    if (status)
      status = RestartReader::load(communicator, options.restart_root,
                                   expected, image);
    if (status) {
      starting_step = image.step;
      run_start.kind = RuntimeRunStartKind::restart;
      run_start.previous_step = image.step;
      run_start.previous_time = image.time;
      run_start.restart_manifest_sha256 = image.source_manifest_sha256;
      restart_requires_recovery = image.backward_euler_recovery;
      okay = rank != 0 ||
              read_restart_binding(options.restart_root, image, fingerprint,
                                   restart_binding);
      if (all_true(communicator, okay)) {
        status = driver.initialize_restart(image);
        restarted = true;
      } else {
        status = {StatusCode::invalid_case, kRunnerInput};
      }
    }
  } else if (status) {
    std::vector<double> scalars;
    const DriverInitialState initial = initial_state(model, scalars);
    status = driver.initialize(initial);
  }
  if (!status) {
    if (rank == 0)
      std::cerr << "initialize_status=" << static_cast<unsigned>(status.code)
                << '/' << status.detail << '\n';
    return 5;
  }

  CommittedOutputSnapshot snapshot;
  status = driver.committed_output_snapshot(snapshot);
  RuntimeGeometry runtime;
  okay = static_cast<bool>(status) && prepare_geometry(spec, snapshot, runtime) &&
         find_field(snapshot, "U", 3U) != nullptr &&
         find_field(snapshot, "pi", 1U) != nullptr &&
         find_field(snapshot, "h", 1U) != nullptr;
  if (!all_true(communicator, okay)) {
    if (rank == 0) std::cerr << "snapshot_geometry_failure\n";
    return 5;
  }

  if (options.dry_plan) {
    const std::array<std::int32_t, 6U> local_patch{{
        snapshot.patch.begin.x, snapshot.patch.begin.y, snapshot.patch.begin.z,
        snapshot.patch.cells.x, snapshot.patch.cells.y,
        snapshot.patch.cells.z}};
    int ranks = 0;
    MPI_Comm_size(communicator, &ranks);
    std::vector<std::int32_t> patches;
    if (rank == 0) patches.resize(static_cast<std::size_t>(ranks) * 6U);
    okay = MPI_Gather(local_patch.data(), 6, MPI_INT32_T,
                      rank == 0 ? patches.data() : nullptr, 6, MPI_INT32_T, 0,
                      communicator) == MPI_SUCCESS;
    if (!all_true(communicator, okay)) return 5;
    if (rank == 0) {
      std::cout << "schema=HUNDUN_V04_THIN_DOMAIN_DRY_PLAN_V1"
                << " spec=" << fingerprint << " case=" << model.fingerprint
                << " product=" << product_fingerprint
                << " cpu=" << cpu_fingerprint << " stl=" << stl_fingerprint
                << " global=" << snapshot.geometry->global_cells().x << ','
                << snapshot.geometry->global_cells().y << ','
                << snapshot.geometry->global_cells().z
                << " lower=" << std::setprecision(17)
                << snapshot.geometry->lower().x << ','
                << snapshot.geometry->lower().y << ','
                << snapshot.geometry->lower().z
                << " upper=" << snapshot.geometry->upper().x << ','
                << snapshot.geometry->upper().y << ','
                << snapshot.geometry->upper().z
                << " arena_doubles=" << summary.arena_doubles
                << " ranks=" << ranks << " step=" << snapshot.step << '\n';
      for (int source = 0; source < ranks; ++source) {
        const std::size_t base = static_cast<std::size_t>(source) * 6U;
        std::cout << "patch rank=" << source << " begin=" << patches[base]
                  << ',' << patches[base + 1U] << ',' << patches[base + 2U]
                  << " cells=" << patches[base + 3U] << ','
                  << patches[base + 4U] << ',' << patches[base + 5U] << '\n';
      }
    }
    return 0;
  }

  if (starting_step >= spec.collection_end_step ||
      options.steps > spec.collection_end_step - starting_step) {
    if (rank == 0) std::cerr << "requested_steps_outside_window\n";
    return 3;
  }
  RuntimeCandidateIdentity candidate_identity;
  status = detail::runtime_candidate_identity(communicator,
                                              candidate_identity);
  if (!status) {
    if (rank == 0)
      std::cerr << "candidate_identity_status="
                << static_cast<unsigned>(status.code) << '/' << status.detail
                << '\n';
    return 5;
  }
  const PlanFingerprint build_identity =
      detail::runtime_sha256_fingerprint(candidate_identity.build_manifest);
  const PlanFingerprint binary_identity =
      detail::runtime_sha256_fingerprint(candidate_identity.executable);
  if (build_identity == 0U || binary_identity == 0U) return 5;
  const Int3 global_cells = snapshot.geometry->global_cells();
  Accumulator accumulator;
  accumulator.profile.assign(spec.station_x_over_d.size() *
                                 static_cast<std::size_t>(global_cells.y) * 6U,
                             0.0);
  accumulator.centerline.assign(
      static_cast<std::size_t>(global_cells.x) * 3U, 0.0);
  if (restarted) {
    okay = true;
    if (rank == 0)
      okay = decode_accumulator(
          restart_binding.accumulator, fingerprint, snapshot.plan,
          snapshot.schema, snapshot.step, snapshot.time, accumulator);
    if (!all_true(communicator, okay)) {
      if (rank == 0) std::cerr << "restart_accumulator_failure\n";
      return 5;
    }
    if (MPI_Bcast(&accumulator.sample_steps, 1, MPI_UINT64_T, 0,
                  communicator) != MPI_SUCCESS)
      return 5;
  }
  if (!consensus_u64(communicator, accumulator.sample_steps)) return 5;

  if (!create_run_root(communicator, rank, options.run_root)) {
    if (rank == 0) std::cerr << "run_root_not_exclusive\n";
    return 6;
  }
  if (rank == 0) {
    std::ostringstream metadata;
    metadata.imbue(std::locale::classic());
    metadata << std::setprecision(17) << kRunMagic << '\n'
             << "spec_fingerprint " << fingerprint << '\n'
             << "case " << model.fingerprint << '\n'
             << "product " << product_fingerprint << '\n'
             << "cpu " << cpu_fingerprint << '\n'
             << "stl " << stl_fingerprint << '\n'
             << "starting_step " << starting_step << '\n'
             << "starting_time " << snapshot.time << '\n'
             << "starting_sample_steps " << accumulator.sample_steps << '\n'
             << "requested_steps " << options.steps << '\n'
             << "visit_interval " << options.visit_interval << '\n'
             << "restarted " << (restarted ? 1 : 0) << '\n'
             << "restart_requires_recovery "
             << (restart_requires_recovery ? 1 : 0) << '\n'
             << "candidate_head " << candidate_identity.head.data() << '\n'
             << "candidate_tree " << candidate_identity.tree.data() << '\n'
             << "build_manifest_sha256 "
             << candidate_identity.build_manifest.data() << '\n'
             << "executable_sha256 " << candidate_identity.executable.data()
             << '\n'
             << "candidate_identity " << candidate_identity.identity.data()
             << '\n'
             << "statistics_eligible 0\n"
             << "end\n";
    okay = write_exclusive(options.run_root / "RUN.meta", metadata.str());
  }
  if (!all_true(communicator, okay)) return 6;

  std::ofstream force;
  std::ofstream health;
  std::ofstream probe;
  if (rank == 0) {
    force.open(options.run_root / "force.csv", std::ios::out);
    health.open(options.run_root / "health.csv", std::ios::out);
    probe.open(options.run_root / "probe.csv", std::ios::out);
    if (force)
      force << "step,time,requested_bdf_order,bdf_order,attempts,"
               "pressure_fx,pressure_fy,pressure_fz,viscous_fx,viscous_fy,"
               "viscous_fz,total_fx,total_fy,total_fz,moment_x,moment_y,"
               "moment_z,cd,cl,cz,certificate_plan,"
               "certificate_terminal_state,certificate_final_flux,"
               "certificate_force,certificate_state,included,exclusion\n";
    if (health)
      health << "step,time,requested_bdf_order,bdf_order,attempts,retry,"
                "restart_recovery,temporal_fallback,terminal_audit_present,"
                "final_flux_revision,eos,eos_tolerance,continuity,"
                "continuity_tolerance,energy,energy_tolerance,closed_mass,"
                "closed_mass_tolerance,gauge,gauge_tolerance,p_min,p_max,"
                "T_min,T_max,rho_min,rho_max,"
                "provisional_cfl_out,provisional_cfl_abs,provisional_cfl_limit,"
                "committed_cfl_out,committed_cfl_abs,committed_cfl_limit,"
                "afc_metrics_applicable,afc_retained_l1_ratio,"
                "afc_minimum_face_alpha,afc_active_correction_faces,"
                "afc_limited_faces,afc_limited_face_fraction,afc_limited,"
                "thermo_predictor_theta,thermo_predictor_limited,"
                "pressure_solve_calls,refinement_solve_calls,linear_iterations,"
                "max_rank_step_ns,included,exclusion\n";
    if (probe)
      probe << "step,time,station,x_over_d,u,v,w,included,exclusion\n";
    okay = static_cast<bool>(force) && static_cast<bool>(health) &&
           static_cast<bool>(probe);
  }
  if (!all_true(communicator, okay)) return 6;

  const double force_scale = 0.5 * spec.rho_ref * spec.u_ref * spec.u_ref *
                             spec.diameter * spec.span;
  if (!std::isfinite(force_scale) || !(force_scale > 0.0)) return 3;
  const std::uint64_t target_step = starting_step + options.steps;
  for (std::uint64_t index = 0U; index < options.steps; ++index) {
    const auto begin = std::chrono::steady_clock::now();
    DriverStepReport step;
    status = driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, step);
    const auto end = std::chrono::steady_clock::now();
    const auto local_nanoseconds = static_cast<std::uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin)
            .count());
    std::uint64_t maximum_nanoseconds = 0U;
    okay = MPI_Allreduce(&local_nanoseconds, &maximum_nanoseconds, 1,
                         MPI_UINT64_T, MPI_MAX, communicator) == MPI_SUCCESS;
    if (!status || !okay || !step.accepted ||
        step.accepted_step != starting_step + index + 1U ||
        !std::isfinite(step.accepted_time) ||
        !terminal_audit_valid(model, step)) {
      if (rank == 0)
        std::cerr << "advance_failure expected_step="
                  << starting_step + index + 1U
                  << " reported_step=" << step.accepted_step
                  << " attempts=" << step.attempts << " status="
                  << static_cast<unsigned>(status.code) << '/' << status.detail
                  << " stage=" << step.failed_stage << '\n';
      return 7;
    }
    const bool startup = step.accepted_step == 1U;
    // Legacy V1 images intentionally take one BE recovery step.  Exact V2
    // images restore t_n/t_{n-1}, rates and flux history and therefore keep
    // BDF2 without losing the first sample of every segment.
    const bool recovery =
        restarted && restart_requires_recovery && index == 0U;
    if ((startup || recovery) && step.effective_bdf.order != 1U) {
      if (rank == 0) std::cerr << "backward_euler_provenance_failure\n";
      return 7;
    }
    if (restarted && !restart_requires_recovery && index == 0U &&
        (step.proposal.bdf.order != 2U || step.effective_bdf.order != 2U)) {
      if (rank == 0) std::cerr << "exact_restart_bdf2_provenance_failure\n";
      return 7;
    }

    SurfaceForce surface;
    FinalForceCertificate certificate;
    status = driver.committed_surface_force(surface, certificate);
    if (status) status = driver.committed_output_snapshot(snapshot);
    ThermoExtrema extrema;
    std::vector<double> probes;
    okay = static_cast<bool>(status) && certificate.valid() &&
           finite_force(surface) && snapshot.committed &&
           snapshot.step == step.accepted_step &&
           snapshot.plan == product_fingerprint;
    if (!all_true(communicator, okay)) {
      if (rank == 0) std::cerr << "committed_observable_failure\n";
      return 7;
    }
    okay = collect_thermo_extrema(communicator, model, thermodynamics,
                                  driver.pressure_reference(), snapshot,
                                  extrema);
    if (okay) okay = collect_probes(communicator, runtime, snapshot, probes);
    if (!all_true(communicator, okay)) {
      if (rank == 0) std::cerr << "committed_observable_failure\n";
      return 7;
    }

    const bool in_window = step.accepted_step > spec.development_steps &&
                           step.accepted_step <= spec.collection_end_step;
    const bool retry = step.attempts != 1U;
    const bool fallback = step.temporal_method_fallback;
    const bool included =
        in_window && !startup && !recovery && !retry && !fallback;
    const char* exclusion =
        included      ? "none"
        : startup     ? "startup"
        : recovery    ? "restart_recovery"
        : retry       ? "retry"
        : fallback    ? "temporal_fallback"
                      : "development";
    if (included) {
      okay = add_sample(spec, runtime, snapshot, accumulator);
      if (!all_true(communicator, okay)) {
        if (rank == 0) std::cerr << "nonfinite_snapshot_observable\n";
        return 7;
      }
    }

    const bool visit =
        options.visit_interval != 0U &&
        (step.accepted_step % options.visit_interval == 0U ||
         step.accepted_step == target_step);
    if (visit)
      status = VisitWriter::write(communicator,
                                  options.run_root / "Visit", services,
                                  snapshot);
    DriverResourceReport global_resources;
    if (status)
      status = append_evidence(
          communicator, options.run_root / "evidence.jsonl", services, model,
          product_fingerprint, cpu_fingerprint, stl_fingerprint,
          candidate_identity, build_identity, binary_identity, run_start,
          step, maximum_nanoseconds, recovery, global_resources);
    if (!status) {
      if (rank == 0)
        std::cerr << "evidence_status=" << static_cast<unsigned>(status.code)
                  << '/' << status.detail << '\n';
      return 6;
    }

    if (rank == 0) {
      force << step.accepted_step << ',' << std::setprecision(17)
            << step.accepted_time << ','
            << static_cast<unsigned>(step.proposal.bdf.order) << ','
            << static_cast<unsigned>(step.effective_bdf.order) << ','
            << step.attempts << ',' << surface.pressure.x << ','
            << surface.pressure.y << ',' << surface.pressure.z << ','
            << surface.viscous.x << ',' << surface.viscous.y << ','
            << surface.viscous.z << ',' << surface.total.x << ','
            << surface.total.y << ',' << surface.total.z << ','
            << surface.moment.x << ',' << surface.moment.y << ','
            << surface.moment.z << ',' << surface.total.x / force_scale << ','
            << surface.total.y / force_scale << ','
            << surface.total.z / force_scale << ',' << certificate.plan << ','
            << certificate.terminal_state << ',' << certificate.final_flux
            << ',' << certificate.force << ',' << certificate.state << ','
            << (included ? 1 : 0) << ',' << exclusion << '\n';
      health
          << step.accepted_step << ',' << std::setprecision(17)
          << step.accepted_time << ','
          << static_cast<unsigned>(step.proposal.bdf.order) << ','
          << static_cast<unsigned>(step.effective_bdf.order) << ','
          << step.attempts << ',' << (retry ? 1 : 0) << ','
          << (recovery ? 1 : 0) << ',' << (fallback ? 1 : 0) << ','
          << 1 << ',' << step.piso.final_flux_revision << ','
          << step.piso.eos_residual << ',' << model.solver.terminal.eos << ','
          << step.piso.continuity_residual << ','
          << model.solver.terminal.continuity << ','
          << step.piso.energy_residual << ','
          << model.solver.terminal.continuity << ','
          << step.piso.closed_mass_residual << ','
          << model.solver.terminal.closed_mass << ','
          << step.piso.gauge_residual << ',' << model.solver.terminal.gauge
          << ',' << extrema.pressure_min << ',' << extrema.pressure_max << ','
          << extrema.temperature_min << ',' << extrema.temperature_max << ','
          << extrema.density_min << ',' << extrema.density_max << ','
          << step.momentum_predictor_limiter.advective_cfl.out_max << ','
          << step.momentum_predictor_limiter.advective_cfl.absolute_max << ','
          << step.momentum_predictor_limiter.advective_cfl.limit << ','
          << step.piso.committed_convective_cfl_out_max << ','
          << step.piso.committed_convective_cfl_abs_max << ','
          << step.piso.committed_convective_cfl_limit << ','
          << (step.momentum_predictor_limiter.correction_metrics_applicable
                  ? 1
                  : 0)
          << ','
          << step.momentum_predictor_limiter.theta << ','
          << step.momentum_predictor_limiter.minimum_face_alpha << ','
          << step.momentum_predictor_limiter.active_correction_faces << ','
          << step.momentum_predictor_limiter.activations << ','
          << step.momentum_predictor_limiter.limited_face_fraction << ','
          << (step.momentum_predictor_limiter.limited ? 1 : 0) << ','
          << step.thermophysical_predictor.theta << ','
          << (step.thermophysical_predictor.limited ? 1 : 0) << ','
          << static_cast<unsigned>(step.piso.pressure_solve_calls) << ','
          << static_cast<unsigned>(
                 step.piso.pressure_energy_refinement_solve_calls)
          << ',' << global_resources.linear_iterations << ','
          << maximum_nanoseconds << ',' << (included ? 1 : 0) << ','
          << exclusion << '\n';
      for (std::size_t station = 0U;
           station < spec.station_x_over_d.size(); ++station) {
        const std::size_t base = station * 3U;
        probe << step.accepted_step << ',' << std::setprecision(17)
              << step.accepted_time << ',' << station << ','
              << spec.station_x_over_d[station] << ',' << probes[base] << ','
              << probes[base + 1U] << ',' << probes[base + 2U] << ','
              << (included ? 1 : 0) << ',' << exclusion << '\n';
      }
      okay = static_cast<bool>(force) && static_cast<bool>(health) &&
             static_cast<bool>(probe);
    }
    if (!all_true(communicator, okay)) return 6;

    const bool save =
        step.accepted_step % spec.checkpoint_interval == 0U ||
        step.accepted_step == target_step;
    if (save) {
      if (rank == 0) {
        force.flush();
        health.flush();
        probe.flush();
        okay = static_cast<bool>(force) && static_cast<bool>(health) &&
               static_cast<bool>(probe);
      }
      if (!all_true(communicator, okay) ||
          !consensus_u64(communicator, accumulator.sample_steps) ||
          !checkpoint(communicator, rank, driver, options.run_root, spec,
                      fingerprint, runtime, snapshot, accumulator)) {
        if (rank == 0) std::cerr << "checkpoint_failure\n";
        return 6;
      }
      if (rank == 0)
        std::cout << "CHECKPOINT step=" << step.accepted_step
                  << " time=" << std::setprecision(17) << step.accepted_time
                  << " samples=" << accumulator.sample_steps << '\n';
    }
  }
  if (rank == 0) {
    force.close();
    health.close();
    probe.close();
    okay = ::chmod((options.run_root / "force.csv").c_str(), 0444) == 0 &&
           ::chmod((options.run_root / "health.csv").c_str(), 0444) == 0 &&
           ::chmod((options.run_root / "probe.csv").c_str(), 0444) == 0 &&
           sync_directory(options.run_root);
  }
  if (!all_true(communicator, okay)) return 6;
  if (rank == 0)
    std::cout << "COMPLETED steps=" << options.steps
              << " final_step=" << target_step
              << " samples=" << accumulator.sample_steps << '\n';
  return 0;
}

bool self_test() {
  const std::array<double, 4U> coordinates{{-1.5, -0.5, 0.5, 1.5}};
  Bracket exact;
  Bracket middle;
  Bracket invalid;
  bool okay = bracket({coordinates.data(), coordinates.size()}, 0.5, exact) &&
              exact.lower == 2 && exact.upper == 2 &&
              exact.upper_weight == 0.0 &&
              bracket({coordinates.data(), coordinates.size()}, 0.0,
                      middle) &&
              middle.lower == 1 && middle.upper == 2 &&
              std::abs(middle.upper_weight - 0.5) < 1.0e-15 &&
              !bracket({coordinates.data(), coordinates.size()}, 2.0,
                       invalid);
  StatisticsSpec spec;
  spec.development_steps = 10U;
  spec.collection_end_step = 20U;
  spec.checkpoint_interval = 5U;
  spec.rho_ref = 2.0;
  spec.u_ref = 3.0;
  spec.diameter = 4.0;
  spec.span = 5.0;
  spec.cylinder_center_x = -1.0;
  spec.station_x_over_d = {1.0, 2.0};
  const std::uint64_t fingerprint = spec_fingerprint(spec);
  okay = okay && fingerprint != 0U &&
         0.5 * spec.rho_ref * spec.u_ref * spec.u_ref * spec.diameter *
                 spec.span ==
             180.0;

  const fs::path root = fs::temp_directory_path() /
                        ("hundun-v04-thin-runner-" +
                         std::to_string(::getpid()));
  std::error_code error;
  fs::remove_all(root, error);
  error.clear();
  okay = okay && fs::create_directory(root, error) && !error;
  const fs::path spec_path = root / "statistics.d";
  const std::string spec_text =
      std::string(kSpecMagic) +
      "\ndevelopment_steps 10\ncollection_end_step 20\n"
      "checkpoint_interval 5\nrho_ref 2\nu_ref 3\ndiameter 4\n"
      "span 5\ncylinder_center_x -1\nstation_x_over_d 1\n"
      "station_x_over_d 2\nend\n";
  okay = okay && write_exclusive(spec_path, spec_text);
  StatisticsSpec parsed;
  std::string parse_error;
  okay = okay && parse_spec(spec_path, parsed, parse_error) &&
         spec_fingerprint(parsed) == fingerprint;

  Accumulator source;
  source.profile = {1.0, 2.0, 4.0};
  source.centerline = {3.0, 5.0};
  source.sample_steps = 7U;
  const fs::path accumulator_path = root / "state.accumulator";
  okay = okay && write_exclusive(
                       accumulator_path,
                       encode_accumulator(fingerprint, 11U, 12U, 9U, 0.25,
                                          source));
  Accumulator restored;
  restored.profile.assign(source.profile.size(), 0.0);
  restored.centerline.assign(source.centerline.size(), 0.0);
  okay = okay && decode_accumulator(accumulator_path, fingerprint, 11U, 12U,
                                    9U, 0.25, restored) &&
         restored.profile == source.profile &&
         restored.centerline == source.centerline &&
         restored.sample_steps == source.sample_steps;
  Accumulator wrong_size;
  wrong_size.profile.assign(2U, 0.0);
  wrong_size.centerline.assign(2U, 0.0);
  okay = okay && !decode_accumulator(accumulator_path, fingerprint, 11U, 12U,
                                     9U, 0.25, wrong_size);
  fs::remove_all(root, error);
  return okay;
}

void usage(int rank) {
  if (rank != 0) return;
  std::cerr
      << "usage:\n"
      << "  v04_thin_domain_runner --self-test\n"
      << "  v04_thin_domain_runner --spec PATH --case-root PATH --dry-plan\n"
      << "  v04_thin_domain_runner --spec PATH --case-root PATH "
         "--run-root PATH [--restart-root PATH] --steps N "
         "[--visit-interval N]\n";
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 2;
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Options options;
  int result = 0;
  if (!parse_options(argc, argv, options)) {
    usage(rank);
    result = 2;
  } else if (options.self_test) {
    const bool local = self_test();
    const bool passed = all_true(MPI_COMM_WORLD, local);
    if (rank == 0)
      std::cout << (passed ? "PASS" : "FAIL")
                << " v04_thin_domain_runner_self_test\n";
    result = passed ? 0 : 1;
  } else {
    result = run(MPI_COMM_WORLD, rank, options);
  }
  MPI_Finalize();
  return result;
}
