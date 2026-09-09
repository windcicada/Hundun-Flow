// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

// Offline bridge from the audited COAST/VTK mean-field transfer arrays to a
// native HUNDUN v0.4 legacy restart.  This intentionally reconstructs h and
// rho with HUNDUN thermodynamics; it is not a COAST PDF restart continuation.

#include "hundun/v04_app.hpp"
#include "hundun/v04_io.hpp"
#include "hundun/v04_mesh.hpp"
#include "hundun/v04_physics.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

namespace {

using namespace hundun::v04;

constexpr double kDefaultTime = 0.015465333126485348;
constexpr double kDefaultDt = 5.456638518808177e-7;
constexpr double kDefaultPressureReference = 100000.0;
constexpr double kSolidTemperature = 295.0;
constexpr double kCarrierOxygenMassFraction = 0.23291751145757963;
constexpr double kCoastGasConstant = 8314.3;

enum Detail : std::uint32_t {
  kArguments = 24001U,
  kInputContract = 24002U,
  kInputFile = 24003U,
  kInputValue = 24004U,
  kSpeciesContract = 24005U,
  kThermodynamics = 24006U,
  kGeometry = 24007U,
  kRestartLayout = 24008U,
  kTopologyMismatch = 24009U,
  kReadback = 24010U,
  kAuditWrite = 24011U,
  kCollective = 24012U,
  kAllocation = 24013U
};

struct Options {
  std::filesystem::path case_root;
  std::filesystem::path transfer_root;
  std::filesystem::path output_root;
  double time{kDefaultTime};
  double dt{kDefaultDt};
  double pressure_reference{kDefaultPressureReference};
  std::uint64_t step{20000U};
};

struct SourceBlock {
  Int3 begin{};
  Int3 cells{};
  std::vector<float> u;
  std::vector<float> v;
  std::vector<float> w;
  std::vector<float> pressure_absolute;
  std::vector<float> temperature;
  std::vector<float> methane;
  std::vector<float> coast_density;
  std::vector<float> coast_flow_enthalpy;
  std::vector<std::uint8_t> fluid;

  std::size_t offset(Int3 global) const noexcept {
    const auto x = static_cast<std::size_t>(global.x - begin.x);
    const auto y = static_cast<std::size_t>(global.y - begin.y);
    const auto z = static_cast<std::size_t>(global.z - begin.z);
    return x + static_cast<std::size_t>(cells.x) *
                   (y + static_cast<std::size_t>(cells.y) * z);
  }
};

struct ReconstructedBlock {
  std::vector<double> enthalpy;
  std::vector<double> density;
  std::vector<double> oxygen;
};

struct Audit {
  std::uint64_t fluid_cells{};
  std::uint64_t solid_cells{};
  std::uint64_t native_flux_changed_owned_faces{};
  std::uint64_t native_internal_source_faces{};
  std::uint64_t bridge_payload_bytes_per_rank{};
  std::uint64_t sealed_arena_bytes_per_rank{};
  std::uint64_t maximum_workspace_bytes_per_rank{};
  std::uint64_t service_staging_bytes_per_rank{};
  std::uint64_t estimated_single_product_peak_bytes_per_rank{};
  long double coast_mass{};
  long double hundun_mass{};
  long double coast_energy{};
  long double hundun_energy{};
  long double coast_kinetic_energy{};
  long double hundun_kinetic_energy{};
  long double coast_total_energy{};
  long double hundun_total_energy{};
  long double same_density_enthalpy_change{};
  long double absolute_energy_change{};
  long double density_error_sum{};
  long double density_error_square_sum{};
  long double native_flux_absolute_change{};
  long double native_internal_source_absolute_mass_flow{};
  double configured_immersed_source_mass_flow{};
  double density_error_max{};
  double field_readback_max{};
  double flux_readback_max{};
};

bool checked_product(Int3 cells, std::size_t& out) noexcept {
  if (cells.x <= 0 || cells.y <= 0 || cells.z <= 0) return false;
  const auto x = static_cast<std::size_t>(cells.x);
  const auto y = static_cast<std::size_t>(cells.y);
  const auto z = static_cast<std::size_t>(cells.z);
  if (x > std::numeric_limits<std::size_t>::max() / y) return false;
  const std::size_t xy = x * y;
  if (xy > std::numeric_limits<std::size_t>::max() / z) return false;
  out = xy * z;
  return true;
}

bool same_int3(Int3 left, Int3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

bool same_patch(MeshPatch left, MeshPatch right) noexcept {
  return same_int3(left.begin, right.begin) &&
         same_int3(left.cells, right.cells) &&
         same_int3(left.process_grid, right.process_grid) &&
         same_int3(left.process_coord, right.process_coord);
}

std::uint64_t fnv_mix(std::uint64_t hash, std::uint8_t byte) noexcept {
  return (hash ^ byte) * UINT64_C(1099511628211);
}

std::uint64_t options_fingerprint(const Options& options) {
  std::uint64_t hash = UINT64_C(1469598103934665603);
  const auto add = [&](std::string_view value) {
    for (unsigned char byte : value) hash = fnv_mix(hash, byte);
    hash = fnv_mix(hash, 0xffU);
  };
  add(options.case_root.generic_string());
  add(options.transfer_root.generic_string());
  add(options.output_root.generic_string());
  std::array<std::uint64_t, 4U> words{};
  std::memcpy(&words[0U], &options.time, sizeof(double));
  std::memcpy(&words[1U], &options.dt, sizeof(double));
  std::memcpy(&words[2U], &options.pressure_reference, sizeof(double));
  words[3U] = options.step;
  for (std::uint64_t word : words)
    for (unsigned shift = 0U; shift < 64U; shift += 8U)
      hash = fnv_mix(hash, static_cast<std::uint8_t>(word >> shift));
  return hash == 0U ? 1U : hash;
}

Status consensus(MPI_Comm communicator, Status local) noexcept {
  int rank = -1;
  int size = 0;
  if (MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(communicator, &size) != MPI_SUCCESS || size <= 0)
    return {StatusCode::mpi_failure, kCollective};
  int failing = local ? size : rank;
  int lowest = size;
  if (MPI_Allreduce(&failing, &lowest, 1, MPI_INT, MPI_MIN, communicator) !=
      MPI_SUCCESS)
    return {StatusCode::mpi_failure, kCollective};
  if (lowest == size) return {};
  std::array<std::uint32_t, 2U> payload{};
  if (rank == lowest) {
    payload[0U] = static_cast<std::uint32_t>(local.code);
    payload[1U] = local.detail;
  }
  if (MPI_Bcast(payload.data(), static_cast<int>(payload.size()), MPI_UINT32_T,
                lowest, communicator) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kCollective};
  return {static_cast<StatusCode>(payload[0U]), payload[1U]};
}

template <class Function>
Status local_stage(MPI_Comm communicator, Function&& function) noexcept {
  Status local;
  try {
    local = function();
  } catch (const std::bad_alloc&) {
    local = {StatusCode::allocation_failure, kAllocation};
  } catch (...) {
    local = {StatusCode::invalid_plan, kInputContract};
  }
  return consensus(communicator, local);
}

bool parse_double(std::string_view text, double& value) {
  std::string owned(text);
  char* end = nullptr;
  errno = 0;
  const double parsed = std::strtod(owned.c_str(), &end);
  if (errno != 0 || end == owned.c_str() || *end != '\0' ||
      !std::isfinite(parsed))
    return false;
  value = parsed;
  return true;
}

bool parse_u64(std::string_view text, std::uint64_t& value) {
  std::string owned(text);
  char* end = nullptr;
  errno = 0;
  const unsigned long long parsed = std::strtoull(owned.c_str(), &end, 10);
  if (errno != 0 || end == owned.c_str() || *end != '\0') return false;
  value = static_cast<std::uint64_t>(parsed);
  return true;
}

bool parse_options(int argc, char** argv, Options& options) {
  for (int index = 1; index < argc; ++index) {
    const std::string_view argument(argv[index]);
    const auto next = [&]() -> const char* {
      return ++index < argc ? argv[index] : nullptr;
    };
    if (argument == "--case") {
      const char* value = next();
      if (value == nullptr) return false;
      options.case_root = value;
    } else if (argument == "--transfer") {
      const char* value = next();
      if (value == nullptr) return false;
      options.transfer_root = value;
    } else if (argument == "--output") {
      const char* value = next();
      if (value == nullptr) return false;
      options.output_root = value;
    } else if (argument == "--time") {
      const char* value = next();
      if (value == nullptr || !parse_double(value, options.time)) return false;
    } else if (argument == "--dt") {
      const char* value = next();
      if (value == nullptr || !parse_double(value, options.dt)) return false;
    } else if (argument == "--pressure-reference") {
      const char* value = next();
      if (value == nullptr ||
          !parse_double(value, options.pressure_reference))
        return false;
    } else if (argument == "--step") {
      const char* value = next();
      if (value == nullptr || !parse_u64(value, options.step)) return false;
    } else {
      return false;
    }
  }
  return !options.case_root.empty() && !options.transfer_root.empty() &&
         !options.output_root.empty() && options.time >= 0.0 &&
         options.dt > 0.0 && options.pressure_reference > 0.0 &&
         options.step > 0U;
}

std::string usage() {
  return
      "usage: v04_coast_restart_import --case CASE_ROOT --transfer "
      "TRANSFER_ROOT --output RESTART_ROOT [--step N] [--time S] [--dt S] "
      "[--pressure-reference PA]\n";
}

float decode_little_float(const std::uint8_t* bytes) noexcept {
  const std::uint32_t bits = static_cast<std::uint32_t>(bytes[0U]) |
                             (static_cast<std::uint32_t>(bytes[1U]) << 8U) |
                             (static_cast<std::uint32_t>(bytes[2U]) << 16U) |
                             (static_cast<std::uint32_t>(bytes[3U]) << 24U);
  float value = 0.0F;
  std::memcpy(&value, &bits, sizeof(value));
  return value;
}

Status check_file_size(const std::filesystem::path& path,
                       std::uintmax_t expected) {
  std::error_code error;
  const auto status = std::filesystem::symlink_status(path, error);
  if (error || !std::filesystem::is_regular_file(status))
    return {StatusCode::io_failure, kInputFile};
  const auto bytes = std::filesystem::file_size(path, error);
  if (error || bytes != expected)
    return {StatusCode::invalid_case, kInputFile};
  return {};
}

Status read_f32_box(const std::filesystem::path& path, Int3 global_cells,
                    Int3 begin, Int3 cells, std::vector<float>& values) {
  std::size_t global_count = 0U;
  std::size_t local_count = 0U;
  if (!checked_product(global_cells, global_count) ||
      !checked_product(cells, local_count) ||
      global_count > std::numeric_limits<std::uintmax_t>::max() / 4U)
    return {StatusCode::invalid_plan, kInputContract};
  Status status = check_file_size(path, static_cast<std::uintmax_t>(global_count) * 4U);
  if (!status) return status;
  values.resize(local_count);
  std::ifstream stream(path, std::ios::binary);
  if (!stream) return {StatusCode::io_failure, kInputFile};
  std::vector<std::uint8_t> row(static_cast<std::size_t>(cells.x) * 4U);
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y) {
      const auto gx = static_cast<std::uint64_t>(begin.x);
      const auto gy = static_cast<std::uint64_t>(begin.y + y);
      const auto gz = static_cast<std::uint64_t>(begin.z + z);
      const auto nx = static_cast<std::uint64_t>(global_cells.x);
      const auto ny = static_cast<std::uint64_t>(global_cells.y);
      const std::uint64_t element = gx + nx * (gy + ny * gz);
      const std::uint64_t byte = element * 4U;
      if (byte > static_cast<std::uint64_t>(
                     std::numeric_limits<std::streamoff>::max()))
        return {StatusCode::invalid_plan, kInputContract};
      stream.seekg(static_cast<std::streamoff>(byte));
      stream.read(reinterpret_cast<char*>(row.data()),
                  static_cast<std::streamsize>(row.size()));
      if (!stream) return {StatusCode::io_failure, kInputFile};
      const std::size_t destination =
          static_cast<std::size_t>(cells.x) *
          (static_cast<std::size_t>(y) +
           static_cast<std::size_t>(cells.y) * z);
      for (std::int32_t x = 0; x < cells.x; ++x)
        values[destination + static_cast<std::size_t>(x)] =
            decode_little_float(row.data() + static_cast<std::size_t>(x) * 4U);
    }
  return {};
}

Status read_u8_box(const std::filesystem::path& path, Int3 global_cells,
                   Int3 begin, Int3 cells,
                   std::vector<std::uint8_t>& values) {
  std::size_t global_count = 0U;
  std::size_t local_count = 0U;
  if (!checked_product(global_cells, global_count) ||
      !checked_product(cells, local_count))
    return {StatusCode::invalid_plan, kInputContract};
  Status status = check_file_size(path, global_count);
  if (!status) return status;
  values.resize(local_count);
  std::ifstream stream(path, std::ios::binary);
  if (!stream) return {StatusCode::io_failure, kInputFile};
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y) {
      const auto gx = static_cast<std::uint64_t>(begin.x);
      const auto gy = static_cast<std::uint64_t>(begin.y + y);
      const auto gz = static_cast<std::uint64_t>(begin.z + z);
      const auto nx = static_cast<std::uint64_t>(global_cells.x);
      const auto ny = static_cast<std::uint64_t>(global_cells.y);
      const std::uint64_t byte = gx + nx * (gy + ny * gz);
      if (byte > static_cast<std::uint64_t>(
                     std::numeric_limits<std::streamoff>::max()))
        return {StatusCode::invalid_plan, kInputContract};
      stream.seekg(static_cast<std::streamoff>(byte));
      const std::size_t destination =
          static_cast<std::size_t>(cells.x) *
          (static_cast<std::size_t>(y) +
           static_cast<std::size_t>(cells.y) * z);
      stream.read(reinterpret_cast<char*>(values.data() + destination),
                  static_cast<std::streamsize>(cells.x));
      if (!stream) return {StatusCode::io_failure, kInputFile};
    }
  return {};
}

Status read_source(const std::filesystem::path& root, Int3 global_cells,
                   MeshPatch patch, SourceBlock& out) {
  const Int3 end{
      std::min(global_cells.x, patch.begin.x + patch.cells.x + 1),
      std::min(global_cells.y, patch.begin.y + patch.cells.y + 1),
      std::min(global_cells.z, patch.begin.z + patch.cells.z + 1)};
  SourceBlock candidate;
  candidate.begin = {std::max(0, patch.begin.x - 1),
                     std::max(0, patch.begin.y - 1),
                     std::max(0, patch.begin.z - 1)};
  candidate.cells = {end.x - candidate.begin.x, end.y - candidate.begin.y,
                     end.z - candidate.begin.z};
  struct Field {
    const char* name;
    std::vector<float>* values;
  };
  const std::array<Field, 8U> fields{{
      {"u.f32", &candidate.u},
      {"v.f32", &candidate.v},
      {"w.f32", &candidate.w},
      {"p_abs.f32", &candidate.pressure_absolute},
      {"temperature.f32", &candidate.temperature},
      {"Y_CH4.f32", &candidate.methane},
      {"rho.f32", &candidate.coast_density},
      {"flow_h.f32", &candidate.coast_flow_enthalpy}}};
  Status status;
  for (const Field& field : fields) {
    status = read_f32_box(root / field.name, global_cells, candidate.begin,
                          candidate.cells, *field.values);
    if (!status) return status;
  }
  status = read_u8_box(root / "fluid_mask.u8", global_cells,
                       candidate.begin, candidate.cells, candidate.fluid);
  if (!status) return status;
  for (std::uint8_t marker : candidate.fluid)
    if (marker > 1U) return {StatusCode::invalid_case, kInputValue};
  out = std::move(candidate);
  return {};
}

Status validate_species_contract(const ValidatedModel& model) {
  if (model.transported_scalars.size() != 2U ||
      model.transported_scalars[0U].stable_name != "CH4" ||
      model.transported_scalars[1U].stable_name != "O2" ||
      model.transported_scalars[0U].role != TransportedScalarRole::species ||
      model.transported_scalars[1U].role != TransportedScalarRole::species ||
      model.thermophysics.species.size() != 3U ||
      model.thermophysics.species[0U].stable_name != "CH4" ||
      model.thermophysics.species[1U].stable_name != "O2" ||
      model.thermophysics.species[2U].stable_name != "N2")
    return {StatusCode::invalid_case, kSpeciesContract};
  return {};
}

Status reconstruct(const SourceBlock& source,
                   const ThermodynamicsPlan& thermodynamics,
                   ReconstructedBlock& out) {
  const std::size_t count = source.fluid.size();
  ReconstructedBlock candidate;
  candidate.enthalpy.resize(count);
  candidate.density.resize(count);
  candidate.oxygen.resize(count);
  const std::array<double, 2U> solid_species{{
      0.0, kCarrierOxygenMassFraction}};
  double solid_h = 0.0;
  double solid_cp = 0.0;
  double solid_r = 0.0;
  Status status = thermodynamics.mixture_enthalpy(
      kSolidTemperature, {solid_species.data(), solid_species.size()},
      solid_h, solid_cp, solid_r);
  if (!status || !(solid_r > 0.0) || !std::isfinite(solid_h))
    return {StatusCode::invalid_case, kThermodynamics};
  for (std::size_t cell = 0U; cell < count; ++cell) {
    if (source.fluid[cell] == 0U) {
      candidate.enthalpy[cell] = solid_h;
      candidate.density[cell] =
          kDefaultPressureReference / (solid_r * kSolidTemperature);
      candidate.oxygen[cell] = kCarrierOxygenMassFraction;
      continue;
    }
    const double u = source.u[cell];
    const double v = source.v[cell];
    const double w = source.w[cell];
    const double pressure = source.pressure_absolute[cell];
    const double temperature = source.temperature[cell];
    const double methane = source.methane[cell];
    const double coast_density = source.coast_density[cell];
    const double coast_h = source.coast_flow_enthalpy[cell];
    if (!std::isfinite(u) || !std::isfinite(v) || !std::isfinite(w) ||
        !std::isfinite(pressure) || !(pressure > 0.0) ||
        !std::isfinite(temperature) ||
        temperature < thermodynamics.minimum_temperature() ||
        temperature > thermodynamics.maximum_temperature() ||
        !std::isfinite(methane) || methane < 0.0 || methane > 1.0 ||
        !std::isfinite(coast_density) || !(coast_density > 0.0) ||
        !std::isfinite(coast_h))
      return {StatusCode::numerical_failure, kInputValue};
    const double oxygen = (1.0 - methane) * kCarrierOxygenMassFraction;
    const std::array<double, 2U> species{{methane, oxygen}};
    double enthalpy = 0.0;
    double cp = 0.0;
    double gas_constant = 0.0;
    status = thermodynamics.mixture_enthalpy(
        temperature, {species.data(), species.size()}, enthalpy, cp,
        gas_constant);
    if (!status || !std::isfinite(enthalpy) || !(cp > 0.0) ||
        !(gas_constant > 0.0))
      return {StatusCode::numerical_failure, kThermodynamics};
    const double density = pressure / (gas_constant * temperature);
    if (!std::isfinite(density) || !(density > 0.0))
      return {StatusCode::numerical_failure, kThermodynamics};
    candidate.enthalpy[cell] = enthalpy;
    candidate.density[cell] = density;
    candidate.oxygen[cell] = oxygen;
  }
  out = std::move(candidate);
  return {};
}

double source_velocity(const SourceBlock& source, Int3 global,
                       std::size_t axis) noexcept {
  const std::size_t cell = source.offset(global);
  if (source.fluid[cell] == 0U) return 0.0;
  if (axis == 0U) return source.u[cell];
  if (axis == 1U) return source.v[cell];
  return source.w[cell];
}

double face_flux(const SourceBlock& source,
                 const ReconstructedBlock& reconstructed,
                 const CartesianGeometryPlan& geometry, Int3 global_face,
                 std::size_t axis) noexcept {
  const Int3 global_cells = geometry.global_cells();
  const std::int32_t face = axis == 0U ? global_face.x
                            : axis == 1U ? global_face.y
                                        : global_face.z;
  const std::int32_t limit = axis == 0U ? global_cells.x
                             : axis == 1U ? global_cells.y
                                         : global_cells.z;
  Int3 left = global_face;
  Int3 right = global_face;
  if (axis == 0U) --left.x;
  else if (axis == 1U) --left.y;
  else --left.z;
  const auto momentum = [&](Int3 cell) {
    const std::size_t offset = source.offset(cell);
    return reconstructed.density[offset] *
           source_velocity(source, cell, axis);
  };
  double face_momentum = 0.0;
  if (face == 0) {
    const std::size_t cell = source.offset(right);
    if (source.fluid[cell] == 0U) return 0.0;
    face_momentum = momentum(right);
  } else if (face == limit) {
    const std::size_t cell = source.offset(left);
    if (source.fluid[cell] == 0U) return 0.0;
    face_momentum = momentum(left);
  } else {
    const std::size_t left_cell = source.offset(left);
    const std::size_t right_cell = source.offset(right);
    if (source.fluid[left_cell] == 0U || source.fluid[right_cell] == 0U)
      return 0.0;
    const AxisMetrics& metrics =
        axis == 0U ? geometry.x() : axis == 1U ? geometry.y() : geometry.z();
    const Span<const double> centres = metrics.centres();
    const Span<const double> faces = metrics.faces();
    const auto lower = static_cast<std::size_t>(face - 1);
    const auto upper = static_cast<std::size_t>(face);
    const double left_distance = faces.data[upper] - centres.data[lower];
    const double right_distance = centres.data[upper] - faces.data[upper];
    face_momentum =
        (right_distance * momentum(left) + left_distance * momentum(right)) /
        (left_distance + right_distance);
  }
  const auto gx = static_cast<std::size_t>(global_face.x);
  const auto gy = static_cast<std::size_t>(global_face.y);
  const auto gz = static_cast<std::size_t>(global_face.z);
  double area = 0.0;
  if (axis == 0U)
    area = geometry.y().widths().data[gy] * geometry.z().widths().data[gz];
  else if (axis == 1U)
    area = geometry.x().widths().data[gx] * geometry.z().widths().data[gz];
  else
    area = geometry.x().widths().data[gx] * geometry.y().widths().data[gy];
  return face_momentum * area;
}

Status fill_image(const RestartExpected& expected,
                  const ValidatedModel& model,
                  const CartesianGeometryPlan& geometry,
                  const SourceBlock& source,
                  const ReconstructedBlock& reconstructed,
                  const Options& options, RestartImage& out) {
  if (!same_int3(expected.global_cells, geometry.global_cells()))
    return {StatusCode::invalid_plan, kGeometry};
  std::size_t local_cells = 0U;
  if (!checked_product(expected.target_patch.cells, local_cells))
    return {StatusCode::invalid_plan, kRestartLayout};
  RestartImage candidate;
  candidate.global_cells = expected.global_cells;
  candidate.patch = expected.target_patch;
  candidate.plan = expected.plan;
  candidate.schema = expected.schema;
  candidate.geometry = expected.geometry;
  candidate.time = options.time;
  candidate.dt = options.dt;
  candidate.pressure_reference = options.pressure_reference;
  candidate.step = options.step;
  candidate.controller_state = 1U;
  candidate.backward_euler_recovery = true;
  candidate.source_format_version = 1U;
  std::size_t scalar = 0U;
  for (std::size_t field_index = 0U; field_index < expected.fields.size;
       ++field_index) {
    const RestartExpectedField descriptor = expected.fields.data[field_index];
    RestartImageField field;
    field.role = descriptor.role;
    field.field = descriptor.field;
    field.components = descriptor.components;
    field.values.resize(local_cells * descriptor.components);
    for (std::int32_t z = 0; z < candidate.patch.cells.z; ++z)
      for (std::int32_t y = 0; y < candidate.patch.cells.y; ++y)
        for (std::int32_t x = 0; x < candidate.patch.cells.x; ++x) {
          const Int3 global{candidate.patch.begin.x + x,
                            candidate.patch.begin.y + y,
                            candidate.patch.begin.z + z};
          const std::size_t input = source.offset(global);
          const std::size_t cell = static_cast<std::size_t>(x) +
              static_cast<std::size_t>(candidate.patch.cells.x) *
                  (static_cast<std::size_t>(y) +
                   static_cast<std::size_t>(candidate.patch.cells.y) * z);
          const bool fluid = source.fluid[input] != 0U;
          switch (descriptor.role) {
            case RestartFieldRole::velocity:
              if (descriptor.components != 3U)
                return {StatusCode::invalid_plan, kRestartLayout};
              field.values[cell * 3U] = fluid ? source.u[input] : 0.0;
              field.values[cell * 3U + 1U] = fluid ? source.v[input] : 0.0;
              field.values[cell * 3U + 2U] = fluid ? source.w[input] : 0.0;
              break;
            case RestartFieldRole::pressure_perturbation:
              if (descriptor.components != 1U)
                return {StatusCode::invalid_plan, kRestartLayout};
              field.values[cell] = fluid
                  ? static_cast<double>(source.pressure_absolute[input]) -
                        options.pressure_reference
                  : kDefaultPressureReference - options.pressure_reference;
              break;
            case RestartFieldRole::pressure_absolute:
              if (descriptor.components != 1U)
                return {StatusCode::invalid_plan, kRestartLayout};
              field.values[cell] = fluid ? source.pressure_absolute[input]
                                         : kDefaultPressureReference;
              break;
            case RestartFieldRole::enthalpy:
              if (descriptor.components != 1U)
                return {StatusCode::invalid_plan, kRestartLayout};
              field.values[cell] = reconstructed.enthalpy[input];
              break;
            case RestartFieldRole::independent_species:
              if (descriptor.components != 1U ||
                  scalar >= model.transported_scalars.size())
                return {StatusCode::invalid_plan, kRestartLayout};
              if (model.transported_scalars[scalar].stable_name == "CH4")
                field.values[cell] = fluid ? source.methane[input] : 0.0;
              else if (model.transported_scalars[scalar].stable_name == "O2")
                field.values[cell] = reconstructed.oxygen[input];
              else
                return {StatusCode::invalid_case, kSpeciesContract};
              break;
            default:
              return {StatusCode::invalid_plan, kRestartLayout};
          }
        }
    if (descriptor.role == RestartFieldRole::independent_species) ++scalar;
    candidate.fields.push_back(std::move(field));
  }
  if (scalar != model.transported_scalars.size())
    return {StatusCode::invalid_plan, kRestartLayout};
  const Int3 cells = candidate.patch.cells;
  for (std::size_t axis = 0U; axis < 3U; ++axis) {
    Int3 extents = cells;
    if (axis == 0U) ++extents.x;
    else if (axis == 1U) ++extents.y;
    else ++extents.z;
    std::size_t face_count = 0U;
    if (!checked_product(extents, face_count))
      return {StatusCode::invalid_plan, kRestartLayout};
    candidate.final_mass_flux[axis].resize(face_count);
    for (std::int32_t z = 0; z < extents.z; ++z)
      for (std::int32_t y = 0; y < extents.y; ++y)
        for (std::int32_t x = 0; x < extents.x; ++x) {
          const Int3 local{x, y, z};
          const Int3 global{candidate.patch.begin.x + x,
                            candidate.patch.begin.y + y,
                            candidate.patch.begin.z + z};
          const std::size_t face = static_cast<std::size_t>(x) +
              static_cast<std::size_t>(extents.x) *
                  (static_cast<std::size_t>(y) +
                   static_cast<std::size_t>(extents.y) * z);
          static_cast<void>(local);
          candidate.final_mass_flux[axis][face] =
              face_flux(source, reconstructed, geometry, global, axis);
        }
  }
  out = std::move(candidate);
  return {};
}

bool output_activity_matches(const CommittedOutputSnapshot& snapshot,
                             const SourceBlock& source) noexcept {
  const Int3 cells = snapshot.patch.cells;
  std::size_t expected = 0U;
  if (!checked_product(cells, expected)) return false;
  if (snapshot.cell_activity.size != 0U &&
      snapshot.cell_activity.size != expected)
    return false;
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const std::size_t cell = static_cast<std::size_t>(x) +
            static_cast<std::size_t>(cells.x) *
                (static_cast<std::size_t>(y) +
                 static_cast<std::size_t>(cells.y) * z);
        const Int3 global{snapshot.patch.begin.x + x,
                          snapshot.patch.begin.y + y,
                          snapshot.patch.begin.z + z};
        const std::uint8_t actual = snapshot.cell_activity.size == 0U
                                        ? 1U
                                        : snapshot.cell_activity.data[cell];
        if (actual != source.fluid[source.offset(global)]) return false;
      }
  return true;
}

std::vector<std::uint8_t> owned_activity(const MeshPatch& patch,
                                         const SourceBlock& source) {
  std::size_t count = 0U;
  if (!checked_product(patch.cells, count)) return {};
  std::vector<std::uint8_t> result(count);
  for (std::int32_t z = 0; z < patch.cells.z; ++z)
    for (std::int32_t y = 0; y < patch.cells.y; ++y)
      for (std::int32_t x = 0; x < patch.cells.x; ++x) {
        const std::size_t cell = static_cast<std::size_t>(x) +
            static_cast<std::size_t>(patch.cells.x) *
                (static_cast<std::size_t>(y) +
                 static_cast<std::size_t>(patch.cells.y) * z);
        result[cell] = source.fluid[source.offset(
            {patch.begin.x + x, patch.begin.y + y, patch.begin.z + z})];
      }
  return result;
}

bool output_activity_matches(
    const CommittedOutputSnapshot& snapshot,
    const std::vector<std::uint8_t>& expected_activity) noexcept {
  std::size_t count = 0U;
  if (!checked_product(snapshot.patch.cells, count) ||
      expected_activity.size() != count ||
      (snapshot.cell_activity.size != 0U &&
       snapshot.cell_activity.size != count))
    return false;
  for (std::size_t cell = 0U; cell < count; ++cell) {
    const std::uint8_t actual = snapshot.cell_activity.size == 0U
                                    ? 1U
                                    : snapshot.cell_activity.data[cell];
    if (actual != expected_activity[cell]) return false;
  }
  return true;
}

Status adopt_native_flux(const CartesianGeometryPlan& geometry,
                         const RestartSnapshot& native,
                         RestartImage& image, Audit& audit) {
  if (!same_patch(native.patch, image.patch))
    return {StatusCode::invalid_plan, kRestartLayout};
  const std::array<ConstFaceFieldView, 3U> source{{
      native.final_mass_flux.x, native.final_mass_flux.y,
      native.final_mass_flux.z}};
  for (std::size_t axis = 0U; axis < 3U; ++axis) {
    const ConstFaceFieldView view = source[axis];
    const std::size_t count = static_cast<std::size_t>(view.extents.x) *
                              view.extents.y * view.extents.z;
    if (image.final_mass_flux[axis].size() != count)
      return {StatusCode::invalid_plan, kRestartLayout};
    Int3 owned = image.patch.cells;
    const std::int32_t patch_end = axis == 0U
        ? image.patch.begin.x + image.patch.cells.x
        : axis == 1U ? image.patch.begin.y + image.patch.cells.y
                     : image.patch.begin.z + image.patch.cells.z;
    const std::int32_t global_end = axis == 0U
        ? geometry.global_cells().x
        : axis == 1U ? geometry.global_cells().y
                     : geometry.global_cells().z;
    if (patch_end == global_end) {
      if (axis == 0U) ++owned.x;
      else if (axis == 1U) ++owned.y;
      else ++owned.z;
    }
    for (std::int32_t z = 0; z < view.extents.z; ++z)
      for (std::int32_t y = 0; y < view.extents.y; ++y)
        for (std::int32_t x = 0; x < view.extents.x; ++x) {
          const std::size_t face = static_cast<std::size_t>(x) +
              static_cast<std::size_t>(view.extents.x) *
                  (static_cast<std::size_t>(y) +
                   static_cast<std::size_t>(view.extents.y) * z);
          const double before = image.final_mass_flux[axis][face];
          const double after = view.unchecked({x, y, z});
          if (!std::isfinite(after))
            return {StatusCode::numerical_failure, kReadback};
          const bool locally_owned = x < owned.x && y < owned.y && z < owned.z;
          if (locally_owned && before != after) {
            ++audit.native_flux_changed_owned_faces;
            audit.native_flux_absolute_change += std::abs(after - before);
          }
          const std::int32_t global_face = axis == 0U
              ? image.patch.begin.x + x
              : axis == 1U ? image.patch.begin.y + y
                           : image.patch.begin.z + z;
          if (locally_owned && global_face > 0 && global_face < global_end &&
              before == 0.0 && after != 0.0) {
            ++audit.native_internal_source_faces;
            audit.native_internal_source_absolute_mass_flow +=
                std::abs(after);
          }
          image.final_mass_flux[axis][face] = after;
        }
  }
  return {};
}

double compare_fields(const RestartImage& expected,
                      const CommittedOutputSnapshot& actual) noexcept {
  if (!same_patch(expected.patch, actual.patch) ||
      expected.fields.size() != actual.fields.size)
    return std::numeric_limits<double>::infinity();
  double error = 0.0;
  const Int3 cells = expected.patch.cells;
  for (std::size_t field_index = 0U; field_index < expected.fields.size();
       ++field_index) {
    const RestartImageField& left = expected.fields[field_index];
    const ConstFieldView right = actual.fields.data[field_index].values;
    if (left.field != right.field || left.components != right.components)
      return std::numeric_limits<double>::infinity();
    for (std::int32_t z = 0; z < cells.z; ++z)
      for (std::int32_t y = 0; y < cells.y; ++y)
        for (std::int32_t x = 0; x < cells.x; ++x) {
          const std::size_t cell = static_cast<std::size_t>(x) +
              static_cast<std::size_t>(cells.x) *
                  (static_cast<std::size_t>(y) +
                   static_cast<std::size_t>(cells.y) * z);
          for (std::uint8_t component = 0U; component < left.components;
               ++component)
            error = std::max(
                error,
                std::abs(left.values[cell * left.components + component] -
                         right.unchecked({x, y, z}, component)));
        }
  }
  return error;
}

double compare_flux(const RestartImage& expected,
                    const RestartSnapshot& actual) noexcept {
  const std::array<ConstFaceFieldView, 3U> right{{
      actual.final_mass_flux.x, actual.final_mass_flux.y,
      actual.final_mass_flux.z}};
  double error = 0.0;
  for (std::size_t axis = 0U; axis < 3U; ++axis) {
    const ConstFaceFieldView view = right[axis];
    if (expected.final_mass_flux[axis].size() !=
        static_cast<std::size_t>(view.extents.x) * view.extents.y *
            view.extents.z)
      return std::numeric_limits<double>::infinity();
    for (std::int32_t z = 0; z < view.extents.z; ++z)
      for (std::int32_t y = 0; y < view.extents.y; ++y)
        for (std::int32_t x = 0; x < view.extents.x; ++x) {
          const std::size_t face = static_cast<std::size_t>(x) +
              static_cast<std::size_t>(view.extents.x) *
                  (static_cast<std::size_t>(y) +
                   static_cast<std::size_t>(view.extents.y) * z);
          error = std::max(
              error, std::abs(expected.final_mass_flux[axis][face] -
                              view.unchecked({x, y, z})));
        }
  }
  return error;
}

Status local_audit(const CartesianGeometryPlan& geometry,
                   const MeshPatch& patch, const SourceBlock& source,
                   const ReconstructedBlock& reconstructed, Audit& audit) {
  Audit candidate;
  const auto dx = geometry.x().widths();
  const auto dy = geometry.y().widths();
  const auto dz = geometry.z().widths();
  for (std::int32_t z = 0; z < patch.cells.z; ++z)
    for (std::int32_t y = 0; y < patch.cells.y; ++y)
      for (std::int32_t x = 0; x < patch.cells.x; ++x) {
        const Int3 global{patch.begin.x + x, patch.begin.y + y,
                          patch.begin.z + z};
        const std::size_t input = source.offset(global);
        if (source.fluid[input] == 0U) {
          ++candidate.solid_cells;
          continue;
        }
        ++candidate.fluid_cells;
        const long double volume =
            static_cast<long double>(dx.data[static_cast<std::size_t>(global.x)]) *
            dy.data[static_cast<std::size_t>(global.y)] *
            dz.data[static_cast<std::size_t>(global.z)];
        const long double coast_rho = source.coast_density[input];
        const long double coast_h = source.coast_flow_enthalpy[input];
        const long double hundun_rho = reconstructed.density[input];
        const long double hundun_h = reconstructed.enthalpy[input];
        const long double coast_energy = volume * coast_rho * coast_h;
        const long double hundun_energy = volume * hundun_rho * hundun_h;
        const long double speed_square =
            static_cast<long double>(source.u[input]) * source.u[input] +
            static_cast<long double>(source.v[input]) * source.v[input] +
            static_cast<long double>(source.w[input]) * source.w[input];
        const long double coast_kinetic =
            0.5L * volume * coast_rho * speed_square;
        const long double hundun_kinetic =
            0.5L * volume * hundun_rho * speed_square;
        const long double pressure_work_inventory =
            volume * source.pressure_absolute[input];
        const long double difference = hundun_energy - coast_energy;
        candidate.coast_mass += volume * coast_rho;
        candidate.hundun_mass += volume * hundun_rho;
        candidate.coast_energy += coast_energy;
        candidate.hundun_energy += hundun_energy;
        candidate.coast_kinetic_energy += coast_kinetic;
        candidate.hundun_kinetic_energy += hundun_kinetic;
        candidate.coast_total_energy +=
            coast_energy - pressure_work_inventory + coast_kinetic;
        candidate.hundun_total_energy +=
            hundun_energy - pressure_work_inventory + hundun_kinetic;
        candidate.same_density_enthalpy_change +=
            volume * coast_rho * (hundun_h - coast_h);
        candidate.absolute_energy_change += std::abs(difference);
        const double density_error =
            reconstructed.density[input] - source.coast_density[input];
        candidate.density_error_sum += density_error;
        candidate.density_error_square_sum +=
            static_cast<long double>(density_error) * density_error;
        candidate.density_error_max =
            std::max(candidate.density_error_max, std::abs(density_error));
      }
  audit = candidate;
  return {};
}

Status reduce_audit(MPI_Comm communicator, Audit& audit) noexcept {
  std::array<std::uint64_t, 4U> counts{{
      audit.fluid_cells, audit.solid_cells,
      audit.native_flux_changed_owned_faces,
      audit.native_internal_source_faces}};
  if (MPI_Allreduce(MPI_IN_PLACE, counts.data(), static_cast<int>(counts.size()),
                    MPI_UINT64_T, MPI_SUM, communicator) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kCollective};
  std::array<long double, 14U> sums{{
      audit.coast_mass,
      audit.hundun_mass,
      audit.coast_energy,
      audit.hundun_energy,
      audit.same_density_enthalpy_change,
      audit.absolute_energy_change,
      audit.density_error_sum,
      audit.density_error_square_sum,
      audit.native_flux_absolute_change,
      audit.native_internal_source_absolute_mass_flow,
      audit.coast_kinetic_energy,
      audit.hundun_kinetic_energy,
      audit.coast_total_energy,
      audit.hundun_total_energy}};
  if (MPI_Allreduce(MPI_IN_PLACE, sums.data(), static_cast<int>(sums.size()),
                    MPI_LONG_DOUBLE, MPI_SUM, communicator) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kCollective};
  std::array<double, 3U> maxima{{audit.density_error_max,
                                audit.field_readback_max,
                                audit.flux_readback_max}};
  if (MPI_Allreduce(MPI_IN_PLACE, maxima.data(), static_cast<int>(maxima.size()),
                    MPI_DOUBLE, MPI_MAX, communicator) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kCollective};
  std::array<std::uint64_t, 5U> memory{{
      audit.bridge_payload_bytes_per_rank,
      audit.sealed_arena_bytes_per_rank,
      audit.maximum_workspace_bytes_per_rank,
      audit.service_staging_bytes_per_rank,
      audit.estimated_single_product_peak_bytes_per_rank}};
  if (MPI_Allreduce(MPI_IN_PLACE, memory.data(), static_cast<int>(memory.size()),
                    MPI_UINT64_T, MPI_MAX, communicator) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kCollective};
  audit.fluid_cells = counts[0U];
  audit.solid_cells = counts[1U];
  audit.native_flux_changed_owned_faces = counts[2U];
  audit.native_internal_source_faces = counts[3U];
  audit.coast_mass = sums[0U];
  audit.hundun_mass = sums[1U];
  audit.coast_energy = sums[2U];
  audit.hundun_energy = sums[3U];
  audit.same_density_enthalpy_change = sums[4U];
  audit.absolute_energy_change = sums[5U];
  audit.density_error_sum = sums[6U];
  audit.density_error_square_sum = sums[7U];
  audit.native_flux_absolute_change = sums[8U];
  audit.native_internal_source_absolute_mass_flow = sums[9U];
  audit.coast_kinetic_energy = sums[10U];
  audit.hundun_kinetic_energy = sums[11U];
  audit.coast_total_energy = sums[12U];
  audit.hundun_total_energy = sums[13U];
  audit.density_error_max = maxima[0U];
  audit.field_readback_max = maxima[1U];
  audit.flux_readback_max = maxima[2U];
  audit.bridge_payload_bytes_per_rank = memory[0U];
  audit.sealed_arena_bytes_per_rank = memory[1U];
  audit.maximum_workspace_bytes_per_rank = memory[2U];
  audit.service_staging_bytes_per_rank = memory[3U];
  audit.estimated_single_product_peak_bytes_per_rank = memory[4U];
  return {};
}

Status reduce_readback_errors(MPI_Comm communicator, Audit& audit) noexcept {
  std::array<double, 2U> maxima{{audit.field_readback_max,
                                audit.flux_readback_max}};
  if (MPI_Allreduce(MPI_IN_PLACE, maxima.data(), static_cast<int>(maxima.size()),
                    MPI_DOUBLE, MPI_MAX, communicator) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kCollective};
  audit.field_readback_max = maxima[0U];
  audit.flux_readback_max = maxima[1U];
  return {};
}

std::string json_string(std::string_view value) {
  std::ostringstream stream;
  stream << '"';
  for (unsigned char byte : value) {
    switch (byte) {
      case '"': stream << "\\\""; break;
      case '\\': stream << "\\\\"; break;
      case '\b': stream << "\\b"; break;
      case '\f': stream << "\\f"; break;
      case '\n': stream << "\\n"; break;
      case '\r': stream << "\\r"; break;
      case '\t': stream << "\\t"; break;
      default:
        if (byte < 0x20U) {
          stream << "\\u" << std::hex << std::setw(4) << std::setfill('0')
                 << static_cast<unsigned>(byte) << std::dec;
        } else {
          stream << static_cast<char>(byte);
        }
    }
  }
  stream << '"';
  return stream.str();
}

Status write_audit(const Options& options, const CartesianGeometryPlan& geometry,
                   const Audit& audit) {
  const std::filesystem::path target =
      options.output_root / "mean-field-transfer.json";
  const std::filesystem::path temporary =
      options.output_root / "mean-field-transfer.json.tmp";
  std::ofstream stream(temporary, std::ios::binary | std::ios::trunc);
  if (!stream) return {StatusCode::io_failure, kAuditWrite};
  const long double energy_change = audit.hundun_energy - audit.coast_energy;
  const long double mass_change = audit.hundun_mass - audit.coast_mass;
  const long double density_mean = audit.fluid_cells == 0U
      ? 0.0L
      : audit.density_error_sum / audit.fluid_cells;
  const long double density_rmse = audit.fluid_cells == 0U
      ? 0.0L
      : std::sqrt(audit.density_error_square_sum / audit.fluid_cells);
  stream << std::setprecision(17)
         << "{\n"
         << "  \"schema\": \"gtmc-hundun-native-mean-field-transfer/v1\",\n"
         << "  \"status\": \"ok\",\n"
         << "  \"scope\": \"HUNDUN development-state transfer; not COAST/PDF restart equivalence\",\n"
         << "  \"case_root\": "
         << json_string(options.case_root.generic_string()) << ",\n"
         << "  \"transfer_root\": "
         << json_string(options.transfer_root.generic_string()) << ",\n"
         << "  \"global_cells_xyz\": [" << geometry.global_cells().x << ", "
         << geometry.global_cells().y << ", " << geometry.global_cells().z
         << "],\n"
         << "  \"time_s\": " << options.time << ",\n"
         << "  \"dt_s\": " << options.dt << ",\n"
         << "  \"step\": " << options.step << ",\n"
         << "  \"pressure_reference_Pa\": "
         << options.pressure_reference << ",\n"
         << "  \"restart_format_version\": 1,\n"
         << "  \"backward_euler_recovery_verified\": true,\n"
         << "  \"mapping\": {\n"
         << "    \"velocity\": \"u/v/w.f32 on fluid cells; zero solid placeholders\",\n"
         << "    \"pressure\": \"pi=double(p_abs.f32)-pressure_reference; source pressure_correction_dp is not used\",\n"
         << "    \"composition\": \"CH4=Y_CH4; O2=(1-Y_CH4)*0.23291751145757963; N2 dependent\",\n"
         << "    \"enthalpy_density\": \"reconstructed from p_abs/T/composition with the compiled HUNDUN NASA model\",\n"
         << "    \"imported_face_flux\": \"metric linear interpolation of HUNDUN rho*U; one-sided at global boundaries; initially zero on source solid-incident faces\",\n"
         << "    \"persisted_face_flux\": \"ProductDriver applies compiled immersed-inlet source authority before the native snapshot is written; outer-boundary faces retain the imported estimate until the recovery step resolves boundary conditions\"\n"
         << "  },\n"
         << "  \"excluded_source_fields\": [\"pressure_correction_dp\", \"mixture_fraction_nvf\", \"flow_h_as_restart_authority\"],\n"
         << "  \"gas_constants_J_per_kmol_K\": {\"COAST\": "
         << kCoastGasConstant << ", \"HUNDUN\": " << kUniversalGasConstant
         << ", \"difference\": "
         << kUniversalGasConstant - kCoastGasConstant << "},\n"
         << "  \"cell_counts\": {\"fluid\": " << audit.fluid_cells
         << ", \"solid_placeholder\": " << audit.solid_cells << "},\n"
         << "  \"mass_kg\": {\"COAST_vtk_rho\": " << audit.coast_mass
         << ", \"HUNDUN_eos\": " << audit.hundun_mass
         << ", \"change\": " << mass_change << "},\n"
         << "  \"enthalpy_inventory_J\": {\"COAST_raw_rho_h\": "
         << audit.coast_energy << ", \"HUNDUN_eos_rho_h\": "
         << audit.hundun_energy << ", \"change\": " << energy_change
         << ", \"sum_absolute_cell_change\": "
         << audit.absolute_energy_change
         << ", \"same_COAST_density_h_change\": "
         << audit.same_density_enthalpy_change << "},\n"
         << "  \"kinetic_energy_J\": {\"COAST_vtk_rho\": "
         << audit.coast_kinetic_energy << ", \"HUNDUN_eos_rho\": "
         << audit.hundun_kinetic_energy << ", \"change\": "
         << audit.hundun_kinetic_energy - audit.coast_kinetic_energy
         << "},\n"
         << "  \"total_energy_rho_e_plus_kinetic_J\": {\"COAST\": "
         << audit.coast_total_energy << ", \"HUNDUN\": "
         << audit.hundun_total_energy << ", \"change\": "
         << audit.hundun_total_energy - audit.coast_total_energy << "},\n"
         << "  \"HUNDUN_rho_minus_COAST_rho_kg_m3\": {\"mean\": "
         << density_mean << ", \"rmse\": " << density_rmse
         << ", \"max_absolute\": " << audit.density_error_max << "},\n"
         << "  \"native_flux_adjustment\": {\"changed_owned_faces\": "
         << audit.native_flux_changed_owned_faces
         << ", \"sum_absolute_change_kg_s\": "
         << audit.native_flux_absolute_change
         << ", \"internal_zero_to_nonzero_source_faces\": "
         << audit.native_internal_source_faces
         << ", \"internal_source_absolute_mass_flow_kg_s\": "
         << audit.native_internal_source_absolute_mass_flow
         << ", \"configured_immersed_source_mass_flow_kg_s\": "
         << audit.configured_immersed_source_mass_flow << "},\n"
         << "  \"memory_estimate_max_per_rank_bytes\": {\"bridge_payload\": "
         << audit.bridge_payload_bytes_per_rank
         << ", \"sealed_arena\": " << audit.sealed_arena_bytes_per_rank
         << ", \"maximum_workspace\": "
         << audit.maximum_workspace_bytes_per_rank
         << ", \"service_staging\": "
         << audit.service_staging_bytes_per_rank
         << ", \"single_product_sum_upper_estimate\": "
         << audit.estimated_single_product_peak_bytes_per_rank << "},\n"
         << "  \"native_readback_max_absolute\": {\"fields\": "
         << audit.field_readback_max << ", \"mass_flux\": "
         << audit.flux_readback_max << "}\n"
         << "}\n";
  stream.flush();
  if (!stream) return {StatusCode::io_failure, kAuditWrite};
  stream.close();
  if (!stream) return {StatusCode::io_failure, kAuditWrite};
  std::error_code error;
  std::filesystem::rename(temporary, target, error);
  if (error) return {StatusCode::io_failure, kAuditWrite};
  return {};
}

int run(MPI_Comm communicator, const Options& options, int rank) {
  const std::uint64_t fingerprint = options_fingerprint(options);
  std::array<std::uint64_t, 2U> fingerprints{{fingerprint, fingerprint}};
  if (MPI_Allreduce(MPI_IN_PLACE, &fingerprints[0U], 1, MPI_UINT64_T,
                    MPI_MIN, communicator) != MPI_SUCCESS ||
      MPI_Allreduce(MPI_IN_PLACE, &fingerprints[1U], 1, MPI_UINT64_T,
                    MPI_MAX, communicator) != MPI_SUCCESS ||
      fingerprints[0U] != fingerprints[1U]) {
    if (rank == 0) std::cerr << "option_consensus_failure\n";
    return 3;
  }

  ValidatedModel model;
  Status status =
      CaseCompiler::load_and_compile(communicator, options.case_root, model);
  if (status) status = consensus(communicator, validate_species_contract(model));
  if (!status) {
    if (rank == 0)
      std::cerr << "case_status=" << static_cast<unsigned>(status.code) << '/'
                << status.detail << '\n';
    return 4;
  }

  ThermodynamicsPlan thermodynamics;
  status = ThermodynamicsPlan::compile(
      model.thermophysics,
      {model.transported_scalars.data(), model.transported_scalars.size()},
      thermodynamics);
  status = consensus(communicator, status);
  CartesianGeometryPlan geometry;
  MeshPatch geometry_patch;
  if (status)
    status = CartesianGeometryCompiler::compile(
        communicator, model.mesh, {}, geometry, geometry_patch);
  if (!status) {
    if (rank == 0)
      std::cerr << "thermo_geometry_status="
                << static_cast<unsigned>(status.code) << '/' << status.detail
                << '\n';
    return 5;
  }

  SourceBlock source;
  status = local_stage(communicator, [&] {
    return read_source(options.transfer_root, geometry.global_cells(),
                       geometry_patch, source);
  });
  ReconstructedBlock reconstructed;
  if (status)
    status = local_stage(communicator, [&] {
      return reconstruct(source, thermodynamics, reconstructed);
    });
  if (!status) {
    if (rank == 0)
      std::cerr << "transfer_status=" << static_cast<unsigned>(status.code)
                << '/' << status.detail << '\n';
    return 6;
  }

  Audit audit;
  status = local_stage(communicator, [&] {
    return local_audit(geometry, geometry_patch, source, reconstructed, audit);
  });
  std::vector<std::uint8_t> expected_activity;
  if (status)
    status = local_stage(communicator, [&] {
      expected_activity = owned_activity(geometry_patch, source);
      return expected_activity.empty()
          ? Status{StatusCode::allocation_failure, kAllocation}
          : Status{};
    });
  if (status && model.patch_inlets)
    for (const PatchInletSpec& patch : model.patch_inlets->patches)
      if (patch.immersed)
        audit.configured_immersed_source_mass_flow +=
            std::abs(patch.boundary.mass_flow_rate);
  if (!status) {
    if (rank == 0)
      std::cerr << "audit_preflight_status="
                << static_cast<unsigned>(status.code) << '/' << status.detail
                << '\n';
    return 6;
  }

  CompiledCasePlan plan;
  status = ProductCompiler::compile(communicator, model, options.case_root,
                                    plan);
  PlanSummary plan_summary;
  if (status) plan_summary = plan.summary();
  ProductDriver driver;
  if (status)
    status = ProductDriver::create(communicator, std::move(plan), driver);
  RestartExpected expected;
  if (status) status = driver.restart_expected(expected);
  if (status && (!same_patch(expected.target_patch, geometry_patch) ||
                 expected.geometry != geometry.fingerprint()))
    status = {StatusCode::invalid_plan, kGeometry};
  status = consensus(communicator, status);
  if (!status) {
    if (rank == 0)
      std::cerr << "product_status=" << static_cast<unsigned>(status.code)
                << '/' << status.detail << '\n';
    return 7;
  }

  RestartImage image;
  status = local_stage(communicator, [&] {
    return fill_image(expected, model, geometry, source, reconstructed,
                      options, image);
  });
  if (status) {
    const auto vector_bytes = [](std::size_t count, std::size_t element) {
      return count > std::numeric_limits<std::uint64_t>::max() / element
          ? std::numeric_limits<std::uint64_t>::max()
          : static_cast<std::uint64_t>(count * element);
    };
    std::uint64_t payload = vector_bytes(source.fluid.size(), 1U);
    const auto add = [&](std::uint64_t value) {
      payload = value > std::numeric_limits<std::uint64_t>::max() - payload
          ? std::numeric_limits<std::uint64_t>::max()
          : payload + value;
    };
    for (const std::vector<float>* values : {
             &source.u, &source.v, &source.w, &source.pressure_absolute,
             &source.temperature, &source.methane, &source.coast_density,
             &source.coast_flow_enthalpy})
      add(vector_bytes(values->size(), sizeof(float)));
    for (const std::vector<double>* values : {
             &reconstructed.enthalpy, &reconstructed.density,
             &reconstructed.oxygen})
      add(vector_bytes(values->size(), sizeof(double)));
    for (const RestartImageField& field : image.fields)
      add(vector_bytes(field.values.size(), sizeof(double)));
    for (const std::vector<double>& flux : image.final_mass_flux)
      add(vector_bytes(flux.size(), sizeof(double)));
    audit.bridge_payload_bytes_per_rank = payload;
    audit.sealed_arena_bytes_per_rank =
        vector_bytes(plan_summary.arena_doubles, sizeof(double));
    audit.maximum_workspace_bytes_per_rank =
        plan_summary.maximum_workspace_bytes;
    audit.service_staging_bytes_per_rank =
        plan_summary.service_staging_bytes;
    std::uint64_t estimate = payload;
    for (std::uint64_t value : {
             audit.sealed_arena_bytes_per_rank,
             audit.maximum_workspace_bytes_per_rank,
             audit.service_staging_bytes_per_rank})
      estimate = value > std::numeric_limits<std::uint64_t>::max() - estimate
          ? std::numeric_limits<std::uint64_t>::max()
          : estimate + value;
    audit.estimated_single_product_peak_bytes_per_rank = estimate;
  }
  if (status) status = driver.initialize_restart(image);
  CommittedOutputSnapshot first_output;
  if (status) status = driver.committed_output_snapshot(first_output);
  if (status && !output_activity_matches(first_output, source))
    status = {StatusCode::invalid_case, kTopologyMismatch};
  if (status && compare_fields(image, first_output) != 0.0)
    status = {StatusCode::invalid_plan, kReadback};
  RestartSnapshot committed;
  if (status) status = driver.committed_restart_snapshot(committed);
  if (status)
    status = local_stage(communicator, [&] {
      return adopt_native_flux(geometry, committed, image, audit);
    });
  status = consensus(communicator, status);
  if (!status) {
    if (rank == 0)
      std::cerr << "initialize_status=" << static_cast<unsigned>(status.code)
                << '/' << status.detail << '\n';
    return 8;
  }
  status = reduce_audit(communicator, audit);
  if (status && audit.configured_immersed_source_mass_flow > 0.0) {
    const long double expected_mdot =
        audit.configured_immersed_source_mass_flow;
    const long double error =
        std::abs(audit.native_internal_source_absolute_mass_flow -
                 expected_mdot);
    const long double tolerance =
        256.0L * std::numeric_limits<double>::epsilon() *
        std::max(1.0L, expected_mdot);
    if (audit.native_internal_source_faces == 0U || error > tolerance)
      status = {StatusCode::invalid_plan, kReadback};
  }
  status = consensus(communicator, status);
  if (!status) {
    if (rank == 0)
      std::cerr << "native_flux_status="
                << static_cast<unsigned>(status.code) << '/' << status.detail
                << '\n';
    return 8;
  }

  // The in-memory source is a genuine current-state-only V1 image.  The
  // driver synthesizes t_(n-1) solely so the next step can recover with BE;
  // never persist that synthetic state as exact V2/V3 history.
  RestartSnapshot legacy = committed;
  legacy.previous_fields = {};
  legacy.accepted_rate_fields = {};
  legacy.previous_rate_fields = {};
  legacy.previous_mass_flux = {};
  legacy.previous_pressure_reference = 0.0;
  legacy.closed_mass_target = 0.0;
  legacy.method_history_signature = 0U;
  RestartWriteReport write_report;
  status = RestartWriter::write(
      communicator, options.output_root, legacy, {1U, &write_report});
  if (!status) {
    if (rank == 0)
      std::cerr << "write_status=" << static_cast<unsigned>(status.code) << '/'
                << status.detail << '\n';
    return 9;
  }

  // The first product is no longer needed after its borrowed snapshot has
  // been synchronously persisted and copied into the owning expected image.
  // Release it before compiling the readback product: on the full GTMC mesh
  // this avoids overlapping two complete solver arenas.
  driver = ProductDriver{};
  source = SourceBlock{};
  reconstructed = ReconstructedBlock{};

  CompiledCasePlan readback_plan;
  status = ProductCompiler::compile(communicator, model, options.case_root,
                                    readback_plan);
  ProductDriver readback_driver;
  if (status)
    status = ProductDriver::create(communicator, std::move(readback_plan),
                                   readback_driver);
  RestartExpected readback_expected;
  if (status) status = readback_driver.restart_expected(readback_expected);
  RestartImage readback;
  RestartReadReport read_report;
  if (status)
    status = RestartReader::load(communicator, options.output_root,
                                 readback_expected, readback, &read_report);
  if (status &&
      (readback.source_format_version != 1U ||
       !readback.backward_euler_recovery ||
       readback.step != options.step || readback.time != options.time ||
       readback.dt != options.dt ||
       readback.pressure_reference != options.pressure_reference))
    status = {StatusCode::invalid_plan, kReadback};
  if (status) status = readback_driver.initialize_restart(readback);
  CommittedOutputSnapshot output;
  if (status) status = readback_driver.committed_output_snapshot(output);
  if (status && !output_activity_matches(output, expected_activity))
    status = {StatusCode::invalid_case, kTopologyMismatch};
  RestartSnapshot readback_snapshot;
  if (status)
    status = readback_driver.committed_restart_snapshot(readback_snapshot);
  status = consensus(communicator, status);
  if (!status) {
    if (rank == 0)
      std::cerr << "readback_status=" << static_cast<unsigned>(status.code)
                << '/' << status.detail << '\n';
    return 10;
  }

  if (status) {
    audit.field_readback_max = compare_fields(image, output);
    audit.flux_readback_max = compare_flux(image, readback_snapshot);
    if (!std::isfinite(audit.field_readback_max) ||
        !std::isfinite(audit.flux_readback_max))
      status = {StatusCode::invalid_plan, kReadback};
  }
  status = consensus(communicator, status);
  if (status) status = reduce_readback_errors(communicator, audit);
  if (status && (audit.field_readback_max != 0.0 ||
                 audit.flux_readback_max != 0.0))
    status = {StatusCode::invalid_plan, kReadback};
  if (status && audit.fluid_cells + audit.solid_cells == 0U)
    status = {StatusCode::invalid_case, kInputValue};
  status = consensus(communicator, status);
  if (status && rank == 0) status = write_audit(options, geometry, audit);
  status = consensus(communicator, status);
  if (!status) {
    if (rank == 0)
      std::cerr << "audit_status=" << static_cast<unsigned>(status.code) << '/'
                << status.detail << '\n';
    return 11;
  }

  if (rank == 0) {
    std::cout << std::setprecision(17)
              << "schema=GTMC_HUNDUN_NATIVE_MEAN_FIELD_TRANSFER_V1"
              << " status=ok"
              << " source_format=1"
              << " backward_euler_recovery=true"
              << " step=" << options.step
              << " time=" << options.time
              << " fluid_cells=" << audit.fluid_cells
              << " solid_cells=" << audit.solid_cells
              << " native_internal_source_faces="
              << audit.native_internal_source_faces
              << " native_internal_source_mdot="
              << audit.native_internal_source_absolute_mass_flow
              << " estimated_peak_bytes_per_rank="
              << audit.estimated_single_product_peak_bytes_per_rank
              << " mass_change_kg=" << audit.hundun_mass - audit.coast_mass
              << " enthalpy_inventory_change_J="
              << audit.hundun_energy - audit.coast_energy
              << " total_energy_change_J="
              << audit.hundun_total_energy - audit.coast_total_energy
              << " field_readback_max=" << audit.field_readback_max
              << " flux_readback_max=" << audit.flux_readback_max << '\n';
  }
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 2;
  int rank = -1;
  int code = 2;
  if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) == MPI_SUCCESS) {
    Options options;
    const bool parsed = parse_options(argc, argv, options);
    const Status parse_status = consensus(
        MPI_COMM_WORLD,
        parsed ? Status{} : Status{StatusCode::invalid_case, kArguments});
    if (!parse_status) {
      if (rank == 0) std::cerr << usage();
      code = 2;
    } else {
      code = run(MPI_COMM_WORLD, options, rank);
    }
  }
  if (MPI_Finalize() != MPI_SUCCESS) return 2;
  return code;
}
