// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

// Production-API endpoint comparator for a full-step/half-step temporal
// convergence experiment.  The implementation intentionally never borrows
// ProductDriver internals: every endpoint is restored through exactly the
// same public compile/load/initialize path used by an application restart.

#include "hundun/v04_app.hpp"
#include "hundun/v04_ibm.hpp"
#include "hundun/v04_physics.hpp"

#include <mpi.h>
#include <unistd.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <new>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

namespace fs = std::filesystem;
using namespace hundun::v04;

constexpr std::uint32_t kComparatorInput = 9901U;
constexpr std::uint32_t kComparatorAuthority = 9902U;
constexpr std::uint32_t kComparatorRestart = 9903U;
constexpr std::uint32_t kComparatorNumeric = 9904U;
constexpr std::uint32_t kComparatorIo = 9905U;
constexpr std::uint64_t kHashOffset = UINT64_C(1469598103934665603);
constexpr std::uint64_t kHashPrime = UINT64_C(1099511628211);

enum class ComparisonMode : std::uint8_t {
  full_half,
  equivalent_endpoints,
};

constexpr std::string_view comparison_mode_name(ComparisonMode mode) noexcept {
  return mode == ComparisonMode::equivalent_endpoints
             ? std::string_view{"equivalent_endpoints"}
             : std::string_view{"full_half_2_to_1"};
}

struct EndpointTiming {
  double initial_dt{};
  double minimum_dt{};
  double maximum_dt{};
  double accepted_dt{};
  std::uint64_t step{};
};

struct TimingIdentity {
  bool model_dt{};
  bool accepted_dt{};
  bool step{};

  bool valid() const noexcept { return model_dt && accepted_dt && step; }
};

struct Options {
  fs::path full_case;
  fs::path full_restart;
  fs::path half_case;
  fs::path half_restart;
  std::optional<fs::path> json;
  std::optional<double> maximum_relative_rms;
  std::optional<double> maximum_relative_linf;
  ComparisonMode comparison_mode{ComparisonMode::full_half};
  bool self_test{};
};

struct RawField {
  RestartFieldRole role{RestartFieldRole::velocity};
  std::uint8_t components{};
  std::vector<double> values;
};

struct EndpointRaw {
  fs::path case_root;
  ValidatedModel model;
  PlanSummary summary{};
  PlanFingerprint case_fingerprint{};
  PlanFingerprint product_fingerprint{};
  PlanFingerprint cpu_fingerprint{};
  PlanFingerprint stl_fingerprint{};
  PlanFingerprint schema_fingerprint{};
  PlanFingerprint geometry_fingerprint{};
  Int3 global_cells{};
  MeshPatch patch{};
  double time{};
  double dt{};
  double pressure_reference{};
  std::uint64_t step{};
  std::vector<RawField> fields;
  std::array<std::vector<double>, 3U> final_mass_flux;
};

struct PhysicalAuthority {
  CartesianGeometryPlan geometry;
  MeshPatch patch{};
  ThermodynamicsPlan thermodynamics;
  TransportPlan transport;
  std::optional<StlScanPlan> scan;
  std::optional<ImmersedSurfacePlan> surface;
  std::optional<EBTopology> topology;
};

struct AuthoritySignature {
  PlanFingerprint geometry{};
  PlanFingerprint thermodynamics{};
  PlanFingerprint transport{};
  PlanFingerprint source_stl{};
  PlanFingerprint surface{};
  PlanFingerprint topology{};
  PlanFingerprint interface_metric{};
};

std::uint64_t bits(double value) noexcept {
  std::uint64_t result{};
  static_assert(sizeof(result) == sizeof(value));
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

bool same_bits(double left, double right) noexcept {
  return bits(left) == bits(right);
}

bool same(Int3 left, Int3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

bool same(MeshPatch left, MeshPatch right) noexcept {
  return same(left.begin, right.begin) && same(left.cells, right.cells) &&
         same(left.process_grid, right.process_grid) &&
         same(left.process_coord, right.process_coord);
}

bool positive(Int3 value) noexcept {
  return value.x > 0 && value.y > 0 && value.z > 0;
}

bool checked_cells(Int3 cells, std::size_t& out) noexcept {
  if (!positive(cells)) return false;
  const std::size_t x = static_cast<std::size_t>(cells.x);
  const std::size_t y = static_cast<std::size_t>(cells.y);
  const std::size_t z = static_cast<std::size_t>(cells.z);
  if (x > std::numeric_limits<std::size_t>::max() / y) return false;
  const std::size_t xy = x * y;
  if (xy > std::numeric_limits<std::size_t>::max() / z) return false;
  out = xy * z;
  return true;
}

std::size_t cell_offset(Int3 cells, Int3 index) noexcept {
  return static_cast<std::size_t>(index.x) +
         static_cast<std::size_t>(cells.x) *
             (static_cast<std::size_t>(index.y) +
              static_cast<std::size_t>(cells.y) *
                  static_cast<std::size_t>(index.z));
}

std::size_t face_offset(Int3 extents, Int3 index) noexcept {
  return static_cast<std::size_t>(index.x) +
         static_cast<std::size_t>(extents.x) *
             (static_cast<std::size_t>(index.y) +
              static_cast<std::size_t>(extents.y) *
                  static_cast<std::size_t>(index.z));
}

Status collective_status(MPI_Comm communicator, Status local) noexcept {
  int rank = 0;
  int size = 0;
  if (MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(communicator, &size) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kComparatorInput};
  int source = local ? size : rank;
  int first = size;
  if (MPI_Allreduce(&source, &first, 1, MPI_INT, MPI_MIN, communicator) !=
      MPI_SUCCESS)
    return {StatusCode::mpi_failure, kComparatorInput};
  if (first == size) return {};
  std::array<std::uint32_t, 2U> wire{};
  if (rank == first) {
    wire[0U] = static_cast<std::uint32_t>(local.code);
    wire[1U] = local.detail;
  }
  if (MPI_Bcast(wire.data(), static_cast<int>(wire.size()), MPI_UINT32_T,
                first, communicator) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kComparatorInput};
  return {static_cast<StatusCode>(wire[0U]), wire[1U]};
}

void hash_word(std::uint64_t& hash, std::uint64_t word) noexcept {
  for (unsigned byte = 0U; byte < 8U; ++byte) {
    hash ^= (word >> (byte * 8U)) & UINT64_C(0xff);
    hash *= kHashPrime;
  }
}

void hash_text(std::uint64_t& hash, std::string_view text) noexcept {
  hash_word(hash, text.size());
  for (const char character : text) {
    const auto value = static_cast<unsigned char>(character);
    hash ^= value;
    hash *= kHashPrime;
  }
}

void hash_path(std::uint64_t& hash, const fs::path& path) {
  hash_text(hash, path.lexically_normal().generic_string());
}

void hash_real3(std::uint64_t& hash, Real3 value) noexcept {
  hash_word(hash, bits(value.x));
  hash_word(hash, bits(value.y));
  hash_word(hash, bits(value.z));
}

void hash_int3(std::uint64_t& hash, Int3 value) noexcept {
  hash_word(hash, static_cast<std::uint32_t>(value.x));
  hash_word(hash, static_cast<std::uint32_t>(value.y));
  hash_word(hash, static_cast<std::uint32_t>(value.z));
}

// Typed physical-case identity.  The three time-step bounds are deliberately
// absent; every other validated model member participates.
std::uint64_t physical_model_fingerprint(const ValidatedModel& model) {
  std::uint64_t hash = kHashOffset;
  const CartesianMeshSpec& mesh = model.mesh;
  hash_word(hash, static_cast<std::uint8_t>(mesh.kind));
  hash_real3(hash, mesh.lower);
  hash_real3(hash, mesh.upper);
  hash_word(hash, mesh.has_exact_cells);
  hash_int3(hash, mesh.exact_cells);
  hash_word(hash, mesh.has_base_spacing);
  hash_real3(hash, mesh.base_spacing);
  hash_real3(hash, mesh.minimum_spacing);
  hash_word(hash, bits(mesh.max_growth_ratio));
  hash_word(hash, mesh.focus_regions.size());
  for (const FocusRegionSpec& focus : mesh.focus_regions) {
    hash_real3(hash, focus.lower);
    hash_real3(hash, focus.upper);
    hash_real3(hash, focus.target_spacing);
  }
  hash_word(hash, mesh.limits.max_global_cells);
  hash_word(hash, mesh.limits.max_memory_bytes_per_rank);
  hash_word(hash, static_cast<std::uint8_t>(model.turbulence));
  hash_word(hash, static_cast<std::uint8_t>(model.pressure_reference));
  for (const BoundaryFaceSpec& face : model.boundaries) {
    hash_word(hash, static_cast<std::uint8_t>(face.flow_kind));
    hash_word(hash, static_cast<std::uint8_t>(face.thermal_kind));
    hash_real3(hash, face.velocity);
    hash_real3(hash, face.direction);
    hash_real3(hash, face.backflow_velocity);
    hash_word(hash, bits(face.mass_flow_rate));
    hash_word(hash, bits(face.pressure));
    hash_word(hash, bits(face.temperature));
    hash_word(hash, bits(face.total_pressure));
    hash_word(hash, bits(face.total_temperature));
    hash_word(hash, bits(face.backflow_temperature));
    hash_word(hash, bits(face.heat_flux));
    hash_word(hash, bits(face.relaxation));
    hash_word(hash, bits(face.mach_limit));
    hash_word(hash, face.allow_backflow);
    hash_word(hash, face.scalars.size());
    for (const ScalarBoundarySpec& scalar : face.scalars) {
      hash_text(hash, scalar.stable_name);
      hash_word(hash, static_cast<std::uint8_t>(scalar.kind));
      hash_word(hash, bits(scalar.value));
      hash_word(hash, static_cast<std::uint8_t>(scalar.backflow_kind));
      hash_word(hash, bits(scalar.backflow_value));
    }
  }
  const PressureLinearSolverSpec& linear = model.solver.pressure;
  hash_word(hash, bits(linear.absolute_tolerance));
  hash_word(hash, bits(linear.relative_tolerance));
  hash_word(hash, linear.maximum_iterations);
  hash_word(hash, linear.true_residual_interval);
  hash_word(hash, linear.krylov_restart);
  hash_word(hash, static_cast<std::uint8_t>(linear.algorithm));
  hash_word(hash, static_cast<std::uint8_t>(linear.mg_correction_scaling));
  hash_word(hash, bits(model.solver.terminal.eos));
  hash_word(hash, bits(model.solver.terminal.continuity));
  hash_word(hash, bits(model.solver.terminal.closed_mass));
  hash_word(hash, bits(model.solver.terminal.gauge));
  hash_word(hash, static_cast<std::uint8_t>(model.schemes.momentum));
  hash_word(hash, static_cast<std::uint8_t>(model.schemes.enthalpy));
  hash_word(hash, static_cast<std::uint8_t>(model.schemes.species));
  hash_word(hash, static_cast<std::uint8_t>(model.schemes.passive_scalar));
  hash_word(hash, static_cast<std::uint8_t>(model.schemes.diffusion));
  hash_word(hash, bits(model.schemes.limiter));
  const TimeControlSpec& time = model.time;
  hash_word(hash, static_cast<std::uint8_t>(time.control));
  hash_word(hash, static_cast<std::uint8_t>(time.scheme));
  hash_word(hash, bits(time.convective_cfl));
  hash_word(hash, bits(time.viscous_cfl));
  hash_word(hash, bits(time.thermal_cfl));
  hash_word(hash, bits(time.species_cfl));
  hash_word(hash, bits(time.acoustic_cfl));
  hash_word(hash, bits(time.maximum_growth));
  hash_word(hash, bits(time.retry_factor));
  hash_word(hash, time.maximum_retries);
  hash_word(hash, bits(time.minimum_bdf_ratio));
  hash_word(hash, bits(time.maximum_bdf_ratio));
  const ThermophysicalSpec& thermo = model.thermophysics;
  hash_path(hash, thermo.data_file);
  hash_word(hash, bits(thermo.minimum_temperature));
  hash_word(hash, bits(thermo.maximum_temperature));
  hash_word(hash, bits(thermo.temperature_relative_tolerance));
  hash_word(hash, thermo.maximum_temperature_iterations);
  hash_word(hash, bits(thermo.closed_mass_relative_tolerance));
  hash_word(hash, thermo.maximum_closed_mass_iterations);
  hash_word(hash, bits(thermo.maximum_closed_mass_relative_step));
  hash_word(hash, thermo.species.size());
  for (const SpeciesThermophysicalSpec& species : thermo.species) {
    hash_text(hash, species.stable_name);
    hash_word(hash, bits(species.molecular_weight));
    hash_word(hash, bits(species.temperature_switch));
    for (double value : species.nasa7_low) hash_word(hash, bits(value));
    for (double value : species.nasa7_high) hash_word(hash, bits(value));
    hash_word(hash, static_cast<std::uint8_t>(species.transport_law));
    hash_word(hash, bits(species.viscosity_reference));
    hash_word(hash, bits(species.transport_reference_temperature));
    hash_word(hash, bits(species.sutherland_temperature));
    hash_word(hash, bits(species.prandtl));
    hash_word(hash, bits(species.conductivity));
  }
  hash_word(hash, model.transported_scalars.size());
  for (const TransportedScalarSpec& scalar : model.transported_scalars) {
    hash_text(hash, scalar.stable_name);
    hash_word(hash, static_cast<std::uint8_t>(scalar.role));
    hash_word(hash, bits(scalar.molecular_schmidt));
    hash_word(hash, bits(scalar.turbulent_schmidt));
  }
  hash_word(hash, model.data_files.size());
  for (const fs::path& path : model.data_files) hash_path(hash, path);
  hash_word(hash, model.immersed_boundary.has_value());
  if (model.immersed_boundary.has_value()) {
    hash_path(hash, model.immersed_boundary->stl_file);
    hash_word(hash,
              static_cast<std::uint8_t>(model.immersed_boundary->fluid_side));
    hash_word(hash, static_cast<std::uint8_t>(
                        model.immersed_boundary->reconstruction_policy));
  }
  if (hash == 0U) hash = 1U;
  return hash;
}

bool same_ibm_audit(const IbmReconstructionAudit& left,
                    const IbmReconstructionAudit& right) noexcept {
  return left.valid == right.valid && left.policy == right.policy &&
         left.standard_reach == right.standard_reach &&
         left.group_count == right.group_count &&
         left.quadratic_groups == right.quadratic_groups &&
         left.linear_groups == right.linear_groups &&
         left.expanded_search_groups == right.expanded_search_groups &&
         left.rank_fallback_groups == right.rank_fallback_groups &&
         left.condition_fallback_groups == right.condition_fallback_groups &&
         left.coverage_fallback_groups == right.coverage_fallback_groups &&
         left.donor_fallback_groups == right.donor_fallback_groups &&
         same_bits(left.maximum_condition_estimate,
                   right.maximum_condition_estimate) &&
         same_bits(left.maximum_functional_l1,
                   right.maximum_functional_l1);
}

bool same_summary(const PlanSummary& left, const PlanSummary& right) noexcept {
  return same(left.global_cells, right.global_cells) &&
         same(left.local_cells, right.local_cells) &&
         left.field_count == right.field_count &&
         left.arena_doubles == right.arena_doubles &&
         left.graph_stage_count == right.graph_stage_count &&
         left.graph_node_count == right.graph_node_count &&
         left.maximum_workspace_bytes == right.maximum_workspace_bytes &&
         left.service_staging_bytes == right.service_staging_bytes &&
         left.pressure_correctors == right.pressure_correctors &&
         same_bits(left.pressure_absolute_tolerance,
                   right.pressure_absolute_tolerance) &&
         same_bits(left.pressure_relative_tolerance,
                   right.pressure_relative_tolerance) &&
         left.pressure_maximum_iterations ==
             right.pressure_maximum_iterations &&
         left.pressure_true_residual_interval ==
             right.pressure_true_residual_interval &&
         left.pressure_krylov_restart == right.pressure_krylov_restart &&
         same_bits(left.terminal_eos_tolerance,
                   right.terminal_eos_tolerance) &&
         same_bits(left.terminal_continuity_tolerance,
                   right.terminal_continuity_tolerance) &&
         same_bits(left.terminal_closed_mass_tolerance,
                   right.terminal_closed_mass_tolerance) &&
         same_bits(left.terminal_gauge_tolerance,
                   right.terminal_gauge_tolerance) &&
         left.immersed == right.immersed &&
         same_ibm_audit(left.ibm_boundary_reconstruction,
                        right.ibm_boundary_reconstruction) &&
         same_ibm_audit(left.ibm_surface_reconstruction,
                        right.ibm_surface_reconstruction) &&
         left.exact_numeric_certified == right.exact_numeric_certified &&
         left.preconditioner_setup_certified ==
             right.preconditioner_setup_certified &&
         left.sealed == right.sealed;
}

Status parse_options(int argc, char** argv, Options& out, std::string& detail) {
  if (argc == 2 && std::string_view(argv[1]) == "--self-test") {
    out.self_test = true;
    return {};
  }
  if (argc < 5) {
    detail = "four positional paths are required";
    return {StatusCode::invalid_case, kComparatorInput};
  }
  Options candidate;
  candidate.full_case = argv[1];
  candidate.full_restart = argv[2];
  candidate.half_case = argv[3];
  candidate.half_restart = argv[4];
  for (int index = 5; index < argc; ++index) {
    const std::string_view argument = argv[index];
    const auto take_value = [&](std::string_view name,
                                std::string_view& value) -> bool {
      if (argument == name && index + 1 < argc) {
        value = argv[++index];
        return true;
      }
      return false;
    };
    std::string_view value;
    try {
      if (take_value("--json", value)) {
        candidate.json = fs::path(value);
      } else if (take_value("--max-relative-rms", value)) {
        const std::string text(value);
        std::size_t consumed = 0U;
        const double parsed = std::stod(text, &consumed);
        if (consumed != text.size()) throw std::invalid_argument("trailing");
        candidate.maximum_relative_rms = parsed;
      } else if (take_value("--max-relative-linf", value)) {
        const std::string text(value);
        std::size_t consumed = 0U;
        const double parsed = std::stod(text, &consumed);
        if (consumed != text.size()) throw std::invalid_argument("trailing");
        candidate.maximum_relative_linf = parsed;
      } else if (argument == "--equivalent-endpoints") {
        candidate.comparison_mode = ComparisonMode::equivalent_endpoints;
      } else {
        detail = "unknown or incomplete option: " + std::string(argument);
        return {StatusCode::invalid_case, kComparatorInput};
      }
    } catch (...) {
      detail = "invalid numeric option: " + std::string(argument);
      return {StatusCode::invalid_case, kComparatorInput};
    }
  }
  if (candidate.full_case.empty() || candidate.full_restart.empty() ||
      candidate.half_case.empty() || candidate.half_restart.empty() ||
      !candidate.maximum_relative_rms.has_value() ||
      !candidate.maximum_relative_linf.has_value() ||
      !std::isfinite(*candidate.maximum_relative_rms) ||
      !std::isfinite(*candidate.maximum_relative_linf) ||
      *candidate.maximum_relative_rms < 0.0 ||
      *candidate.maximum_relative_linf < 0.0) {
    detail = "paths and finite non-negative comparison thresholds are required";
    return {StatusCode::invalid_case, kComparatorInput};
  }
  out = std::move(candidate);
  return {};
}

const RawField* find_field(const EndpointRaw& endpoint,
                           RestartFieldRole role,
                           std::size_t occurrence = 0U) noexcept {
  for (const RawField& field : endpoint.fields) {
    if (field.role == role) {
      if (occurrence == 0U) return &field;
      --occurrence;
    }
  }
  return nullptr;
}

bool snapshot_field_matches(const RestartImageField& image,
                            ConstFieldView snapshot) noexcept {
  std::size_t cells = 0U;
  if (!checked_cells(snapshot.interior, cells) ||
      snapshot.components != image.components ||
      image.values.size() != cells * image.components)
    return false;
  for (std::int32_t z = 0; z < snapshot.interior.z; ++z)
    for (std::int32_t y = 0; y < snapshot.interior.y; ++y)
      for (std::int32_t x = 0; x < snapshot.interior.x; ++x) {
        const Int3 index{x, y, z};
        const std::size_t cell = cell_offset(snapshot.interior, index);
        for (std::uint8_t component = 0U; component < image.components;
             ++component)
          if (!same_bits(snapshot.unchecked(index, component),
                         image.values[cell * image.components + component]))
            return false;
      }
  return true;
}

bool snapshot_flux_matches(const std::vector<double>& image,
                           ConstFaceFieldView snapshot) noexcept {
  std::size_t faces = 0U;
  if (!checked_cells(snapshot.extents, faces) || image.size() != faces)
    return false;
  for (std::int32_t z = 0; z < snapshot.extents.z; ++z)
    for (std::int32_t y = 0; y < snapshot.extents.y; ++y)
      for (std::int32_t x = 0; x < snapshot.extents.x; ++x) {
        const Int3 index{x, y, z};
        if (!same_bits(snapshot.unchecked(index),
                       image[face_offset(snapshot.extents, index)]))
          return false;
      }
  return true;
}

Status load_endpoint(MPI_Comm communicator, const fs::path& case_root,
                     const fs::path& restart_directory,
                     EndpointRaw& out) noexcept try {
  ValidatedModel model;
  Status status =
      CaseCompiler::load_and_compile(communicator, case_root, model);
  if (!status) return status;
  CompiledCasePlan plan;
  status = ProductCompiler::compile(communicator, model, case_root, plan);
  if (!status) return status;
  const PlanFingerprint product_fingerprint = plan.fingerprint();
  const PlanFingerprint cpu_fingerprint = plan.cpu_plan_fingerprint();
  const PlanFingerprint stl_fingerprint = plan.stl_fingerprint();
  const PlanSummary summary = plan.summary();
  ProductDriver driver;
  status = ProductDriver::create(communicator, std::move(plan), driver);
  if (!status) return status;
  RestartExpected expected;
  status = driver.restart_expected(expected);
  if (!status) return status;
  RestartImage image;
  status = RestartReader::load(communicator, restart_directory, expected,
                               image);
  if (!status) return status;
  status = driver.initialize_restart(image);
  if (!status) return status;
  RestartSnapshot snapshot;
  status = driver.committed_restart_snapshot(snapshot);
  if (!status) return status;
  Status local{};
  if (!driver.initialized() || !same(image.global_cells, expected.global_cells) ||
      !same(image.patch, expected.target_patch) ||
      image.plan != expected.plan || image.schema != expected.schema ||
      image.geometry != expected.geometry ||
      !same(snapshot.global_cells, image.global_cells) ||
      !same(snapshot.patch, image.patch) || snapshot.plan != image.plan ||
      snapshot.schema != image.schema || snapshot.geometry != image.geometry ||
      !same_bits(snapshot.time, image.time) ||
      !same_bits(snapshot.dt, image.dt) ||
      !same_bits(snapshot.pressure_reference, image.pressure_reference) ||
      snapshot.step != image.step || snapshot.fields.size != image.fields.size() ||
      !snapshot.final_mass_flux.certificate.valid() ||
      !snapshot.final_mass_flux.certificate.matches(
          snapshot.final_mass_flux))
    local = {StatusCode::invalid_plan, kComparatorRestart};
  for (std::size_t index = 0U; local && index < image.fields.size(); ++index) {
    if (snapshot.fields.data[index].role != image.fields[index].role ||
        !snapshot_field_matches(image.fields[index],
                                snapshot.fields.data[index].values))
      local = {StatusCode::invalid_plan, kComparatorRestart};
  }
  if (local &&
      (!snapshot_flux_matches(image.final_mass_flux[0U],
                              snapshot.final_mass_flux.x) ||
       !snapshot_flux_matches(image.final_mass_flux[1U],
                              snapshot.final_mass_flux.y) ||
       !snapshot_flux_matches(image.final_mass_flux[2U],
                              snapshot.final_mass_flux.z)))
    local = {StatusCode::invalid_plan, kComparatorRestart};
  status = collective_status(communicator, local);
  if (!status) return status;
  EndpointRaw candidate;
  Status packing{};
  try {
    candidate.case_root = case_root;
    candidate.case_fingerprint = model.fingerprint;
    candidate.product_fingerprint = product_fingerprint;
    candidate.cpu_fingerprint = cpu_fingerprint;
    candidate.stl_fingerprint = stl_fingerprint;
    candidate.summary = summary;
    candidate.global_cells = image.global_cells;
    candidate.patch = image.patch;
    candidate.schema_fingerprint = image.schema;
    candidate.geometry_fingerprint = image.geometry;
    candidate.time = image.time;
    candidate.dt = image.dt;
    candidate.pressure_reference = image.pressure_reference;
    candidate.step = image.step;
    candidate.fields.reserve(image.fields.size());
    for (RestartImageField& field : image.fields)
      candidate.fields.push_back(
          {field.role, field.components, std::move(field.values)});
    candidate.final_mass_flux = std::move(image.final_mass_flux);
    candidate.model = std::move(model);
  } catch (...) {
    packing = {StatusCode::allocation_failure, kComparatorRestart};
  }
  status = collective_status(communicator, packing);
  if (!status) return status;
  out = std::move(candidate);
  return {};
} catch (const std::bad_alloc&) {
  return collective_status(
      communicator, {StatusCode::allocation_failure, kComparatorRestart});
} catch (...) {
  return collective_status(communicator,
                           {StatusCode::invalid_case, kComparatorRestart});
}

Status compile_authority(MPI_Comm communicator, const EndpointRaw& endpoint,
                         PhysicalAuthority& out,
                         AuthoritySignature& signature) noexcept try {
  PhysicalAuthority candidate;
  Status status = CartesianGeometryCompiler::compile(
      communicator, endpoint.model.mesh, {}, candidate.geometry,
      candidate.patch);
  if (status)
    status = ThermodynamicsPlan::compile(endpoint.model.thermophysics,
                                         {endpoint.model.transported_scalars.data(),
                                          endpoint.model.transported_scalars.size()},
                                         candidate.thermodynamics);
  if (status)
    status = TransportPlan::compile(endpoint.model.thermophysics,
                                    candidate.thermodynamics,
                                    candidate.transport);
  if (status && endpoint.model.immersed_boundary.has_value()) {
    candidate.scan.emplace();
    const StlScanBudget scan_budget{
        endpoint.model.mesh.limits.max_memory_bytes_per_rank / 2U,
        endpoint.model.mesh.limits.max_memory_bytes_per_rank,
        UINT64_C(16000000), UINT64_C(100000), 1U};
    status = StlScanCompiler::compile(
        communicator, endpoint.case_root,
        endpoint.model.immersed_boundary->stl_file, candidate.geometry,
        candidate.patch, CartesianAxis::y, scan_budget, *candidate.scan);
    if (status) {
      candidate.surface.emplace();
      status = ImmersedSurfaceCompiler::compile(*candidate.scan,
                                                *candidate.surface);
    }
    if (status) {
      candidate.topology.emplace();
      status = EBTopologyCompiler::compile(
          communicator, candidate.geometry, candidate.patch, *candidate.scan,
          *candidate.surface,
          endpoint.model.immersed_boundary->fluid_side, ImmersedPlanLimits{},
          *candidate.topology);
    }
  }
  if (!status) return status;
  Status local{};
  if (!same(candidate.patch, endpoint.patch) ||
      candidate.geometry.fingerprint() != endpoint.geometry_fingerprint)
    local = {StatusCode::invalid_plan, kComparatorAuthority};
  status = collective_status(communicator, local);
  if (!status) return status;
  AuthoritySignature authority_signature;
  authority_signature.geometry = candidate.geometry.fingerprint();
  authority_signature.thermodynamics = candidate.thermodynamics.fingerprint();
  authority_signature.transport = candidate.transport.fingerprint();
  if (candidate.surface.has_value()) {
    authority_signature.source_stl =
        candidate.surface->source_triangle_fingerprint();
    authority_signature.surface = candidate.surface->fingerprint();
  }
  if (candidate.topology.has_value()) {
    authority_signature.topology = candidate.topology->fingerprint();
    authority_signature.interface_metric =
        candidate.topology->interface_metric().physical_fingerprint();
  }
  signature = authority_signature;
  out = std::move(candidate);
  return {};
} catch (const std::bad_alloc&) {
  return collective_status(
      communicator, {StatusCode::allocation_failure, kComparatorAuthority});
} catch (...) {
  return collective_status(communicator,
                           {StatusCode::invalid_case, kComparatorAuthority});
}

struct ErrorMetric {
  double weight{};
  double difference_squared{};
  double full_squared{};
  double half_squared{};
  double difference_linf{};
  double full_linf{};
  double half_linf{};
  double rms{};
  double relative_rms{};
  double relative_linf{};
};

enum class MetricId : std::size_t {
  velocity,
  pressure_perturbation,
  pressure_absolute,
  enthalpy,
  density,
  temperature,
  sensible_energy,
  momentum,
  mass_flux,
  mass_flux_density,
  count
};

constexpr std::array<std::string_view,
                     static_cast<std::size_t>(MetricId::count)>
    kMetricNames{"U", "pi", "p_abs", "h", "rho", "T", "q", "rhoU",
                 "phi", "phi_over_A"};
constexpr std::array<std::string_view, 6U> kBoundaryNames{
    "x_min", "x_max", "y_min", "y_max", "z_min", "z_max"};

struct EndpointIntegral {
  double volume{};
  double mass{};
  double energy{};
  double kinetic{};
  double minimum_pressure{std::numeric_limits<double>::infinity()};
  double maximum_pressure{-std::numeric_limits<double>::infinity()};
  double minimum_density{std::numeric_limits<double>::infinity()};
  double maximum_density{-std::numeric_limits<double>::infinity()};
  double minimum_temperature{std::numeric_limits<double>::infinity()};
  double maximum_temperature{-std::numeric_limits<double>::infinity()};
  double minimum_enthalpy{std::numeric_limits<double>::infinity()};
  double maximum_enthalpy{-std::numeric_limits<double>::infinity()};
  double minimum_energy{std::numeric_limits<double>::infinity()};
  double maximum_energy{-std::numeric_limits<double>::infinity()};
  double maximum_mach{};
  std::uint64_t nonfinite_cells{};
  std::uint64_t nonpositive_pressure{};
  std::uint64_t nonpositive_density{};
  std::uint64_t nonpositive_temperature{};
  std::uint64_t mach_limit_violations{};
  std::uint64_t solid_velocity_nonpositive_zero{};
};

struct BoundaryFaceIntegral {
  double net_outward{};
  double inflow{};
  double outflow{};
  double maximum_absolute_flux{};
  std::uint64_t inlet_reversal_faces{};
  std::uint64_t velocity_inlet_flux_mismatch_faces{};
  std::uint64_t outlet_backflow_faces{};
  std::uint64_t prohibited_backflow_faces{};
  std::uint64_t symmetry_nonpositive_zero{};
  std::uint64_t symmetry_negative_zero{};
};

struct EndpointFaceIntegral {
  std::array<BoundaryFaceIntegral, 6U> boundary{};
  std::uint64_t nonfinite_faces{};
  std::uint64_t ibm_nonpositive_zero{};
  std::uint64_t ibm_negative_zero{};
};

struct ComparisonReport {
  std::uint64_t physical_model_fingerprint{};
  AuthoritySignature authority{};
  std::array<ErrorMetric, static_cast<std::size_t>(MetricId::count)> metrics{};
  EndpointIntegral full{};
  EndpointIntegral half{};
  EndpointFaceIntegral full_faces{};
  EndpointFaceIntegral half_faces{};
  double time_tolerance{};
  double maximum_allowed_mach{};
  std::vector<std::string> failures;
};

bool same(AuthoritySignature left, AuthoritySignature right) noexcept {
  return left.geometry == right.geometry &&
         left.thermodynamics == right.thermodynamics &&
         left.transport == right.transport &&
         left.source_stl == right.source_stl && left.surface == right.surface &&
         left.topology == right.topology &&
         left.interface_metric == right.interface_metric;
}

void add_scalar(ErrorMetric& metric, double weight, double full,
                double half) noexcept {
  const double difference = half - full;
  metric.weight += weight;
  metric.difference_squared += weight * difference * difference;
  metric.full_squared += weight * full * full;
  metric.half_squared += weight * half * half;
  metric.difference_linf =
      std::max(metric.difference_linf, std::abs(difference));
  metric.full_linf = std::max(metric.full_linf, std::abs(full));
  metric.half_linf = std::max(metric.half_linf, std::abs(half));
}

void add_vector(ErrorMetric& metric, double weight, Real3 full,
                Real3 half) noexcept {
  const Real3 difference{half.x - full.x, half.y - full.y, half.z - full.z};
  const auto norm = [](Real3 value) noexcept {
    return std::sqrt(value.x * value.x + value.y * value.y +
                     value.z * value.z);
  };
  metric.weight += weight;
  metric.difference_squared +=
      weight * (difference.x * difference.x + difference.y * difference.y +
                difference.z * difference.z);
  metric.full_squared +=
      weight * (full.x * full.x + full.y * full.y + full.z * full.z);
  metric.half_squared +=
      weight * (half.x * half.x + half.y * half.y + half.z * half.z);
  metric.difference_linf = std::max(metric.difference_linf, norm(difference));
  metric.full_linf = std::max(metric.full_linf, norm(full));
  metric.half_linf = std::max(metric.half_linf, norm(half));
}

double relative_error(double difference, double full, double half) noexcept {
  const double scale = std::max(full, half);
  if (scale > 0.0) return difference / scale;
  return difference == 0.0 ? 0.0 : std::numeric_limits<double>::infinity();
}

bool finite(Real3 value) noexcept {
  return std::isfinite(value.x) && std::isfinite(value.y) &&
         std::isfinite(value.z);
}

Real3 field_vector(const RawField& field, std::size_t cell) noexcept {
  return {field.values[cell * 3U], field.values[cell * 3U + 1U],
          field.values[cell * 3U + 2U]};
}

double field_scalar(const RawField& field, std::size_t cell) noexcept {
  return field.values[cell];
}

double maximum_allowed_mach(const ValidatedModel& model) noexcept {
  double limit = std::numeric_limits<double>::infinity();
  for (const BoundaryFaceSpec& boundary : model.boundaries)
    if (std::isfinite(boundary.mach_limit) && boundary.mach_limit > 0.0)
      limit = std::min(limit, boundary.mach_limit);
  return limit;
}

Status validate_endpoint_layout(MPI_Comm communicator,
                                const EndpointRaw& full,
                                const EndpointRaw& half,
                                const PhysicalAuthority& authority,
                                std::string_view& local_detail) noexcept {
  Status local{};
  std::size_t cells = 0U;
  const RawField* full_u = find_field(full, RestartFieldRole::velocity);
  const RawField* half_u = find_field(half, RestartFieldRole::velocity);
  const RawField* full_pi =
      find_field(full, RestartFieldRole::pressure_perturbation);
  const RawField* half_pi =
      find_field(half, RestartFieldRole::pressure_perturbation);
  const RawField* full_h = find_field(full, RestartFieldRole::enthalpy);
  const RawField* half_h = find_field(half, RestartFieldRole::enthalpy);
  if (!same(full.global_cells, half.global_cells) ||
      !same(full.patch, half.patch) || !same(full.patch, authority.patch) ||
      !checked_cells(full.patch.cells, cells)) {
    local_detail = "endpoint decomposition or cell shape differs";
    local = {StatusCode::invalid_plan, kComparatorNumeric};
  } else if (full_u == nullptr || half_u == nullptr || full_pi == nullptr ||
             half_pi == nullptr || full_h == nullptr || half_h == nullptr ||
             full_u->components != 3U || half_u->components != 3U ||
             full_pi->components != 1U || half_pi->components != 1U ||
             full_h->components != 1U || half_h->components != 1U) {
    local_detail = "required U/pi/h Restart roles are absent or malformed";
    local = {StatusCode::invalid_plan, kComparatorNumeric};
  } else if (full.fields.size() != half.fields.size()) {
    local_detail = "Restart field catalogs differ";
    local = {StatusCode::invalid_plan, kComparatorNumeric};
  }
  for (std::size_t index = 0U; local && index < full.fields.size(); ++index) {
    if (full.fields[index].role != half.fields[index].role ||
        full.fields[index].components != half.fields[index].components ||
        full.fields[index].values.size() !=
            cells * full.fields[index].components ||
        half.fields[index].values.size() !=
            cells * half.fields[index].components) {
      local_detail = "Restart field role/component/layout differs";
      local = {StatusCode::invalid_plan, kComparatorNumeric};
    }
  }
  std::size_t independent = 0U;
  for (const RawField& field : full.fields)
    independent +=
        field.role == RestartFieldRole::independent_species ? 1U : 0U;
  if (local &&
      independent != authority.thermodynamics.independent_species_count()) {
    local_detail = "independent-species Restart catalog differs from EOS";
    local = {StatusCode::invalid_plan, kComparatorNumeric};
  }
  if (local && authority.topology.has_value() &&
      authority.topology->region().size != cells) {
    local_detail = "production IBM region mask does not match the patch";
    local = {StatusCode::invalid_plan, kComparatorNumeric};
  }
  const std::array<Int3, 3U> extents{
      Int3{full.patch.cells.x + 1, full.patch.cells.y, full.patch.cells.z},
      Int3{full.patch.cells.x, full.patch.cells.y + 1, full.patch.cells.z},
      Int3{full.patch.cells.x, full.patch.cells.y, full.patch.cells.z + 1}};
  for (std::size_t axis = 0U; local && axis < 3U; ++axis) {
    std::size_t faces = 0U;
    if (!checked_cells(extents[axis], faces) ||
        full.final_mass_flux[axis].size() != faces ||
        half.final_mass_flux[axis].size() != faces) {
      local_detail = "final mass-flux layout differs";
      local = {StatusCode::invalid_plan, kComparatorNumeric};
    }
  }
  return collective_status(communicator, local);
}

void update_integral(EndpointIntegral& integral, double volume,
                     double p_abs, double h, double q, Real3 velocity,
                     const ThermoState& thermo, double mach_limit) noexcept {
  integral.volume += volume;
  integral.mass += volume * thermo.rho;
  integral.energy += volume * q;
  integral.kinetic +=
      volume * 0.5 * thermo.rho *
      (velocity.x * velocity.x + velocity.y * velocity.y +
       velocity.z * velocity.z);
  integral.minimum_pressure = std::min(integral.minimum_pressure, p_abs);
  integral.maximum_pressure = std::max(integral.maximum_pressure, p_abs);
  integral.minimum_density = std::min(integral.minimum_density, thermo.rho);
  integral.maximum_density = std::max(integral.maximum_density, thermo.rho);
  integral.minimum_temperature =
      std::min(integral.minimum_temperature, thermo.temperature);
  integral.maximum_temperature =
      std::max(integral.maximum_temperature, thermo.temperature);
  integral.minimum_enthalpy = std::min(integral.minimum_enthalpy, h);
  integral.maximum_enthalpy = std::max(integral.maximum_enthalpy, h);
  integral.minimum_energy = std::min(integral.minimum_energy, q);
  integral.maximum_energy = std::max(integral.maximum_energy, q);
  integral.maximum_mach = std::max(integral.maximum_mach, thermo.mach);
  const bool all_finite =
      std::isfinite(p_abs) && std::isfinite(h) && std::isfinite(q) &&
      finite(velocity) && std::isfinite(thermo.rho) &&
      std::isfinite(thermo.temperature) && std::isfinite(thermo.mach);
  integral.nonfinite_cells += all_finite ? 0U : 1U;
  integral.nonpositive_pressure += p_abs > 0.0 ? 0U : 1U;
  integral.nonpositive_density += thermo.rho > 0.0 ? 0U : 1U;
  integral.nonpositive_temperature += thermo.temperature > 0.0 ? 0U : 1U;
  integral.mach_limit_violations +=
      thermo.mach <= mach_limit ? 0U : 1U;
}

Status reduce_cell_results(
    MPI_Comm communicator,
    std::array<ErrorMetric, static_cast<std::size_t>(MetricId::count)>& metrics,
    EndpointIntegral& full, EndpointIntegral& half) noexcept {
  constexpr std::size_t count = static_cast<std::size_t>(MetricId::count);
  std::array<double, count * 4U> local_sums{};
  std::array<double, count * 4U> global_sums{};
  std::array<double, count * 3U> local_maxima{};
  std::array<double, count * 3U> global_maxima{};
  for (std::size_t index = 0U; index < count; ++index) {
    local_sums[index * 4U] = metrics[index].weight;
    local_sums[index * 4U + 1U] = metrics[index].difference_squared;
    local_sums[index * 4U + 2U] = metrics[index].full_squared;
    local_sums[index * 4U + 3U] = metrics[index].half_squared;
    local_maxima[index * 3U] = metrics[index].difference_linf;
    local_maxima[index * 3U + 1U] = metrics[index].full_linf;
    local_maxima[index * 3U + 2U] = metrics[index].half_linf;
  }
  std::array<double, 8U> local_integrals{
      full.volume, full.mass, full.energy, full.kinetic,
      half.volume, half.mass, half.energy, half.kinetic};
  std::array<double, 8U> global_integrals{};
  std::array<double, 10U> local_minima{
      full.minimum_pressure, full.minimum_density, full.minimum_temperature,
      full.minimum_enthalpy, full.minimum_energy, half.minimum_pressure,
      half.minimum_density, half.minimum_temperature, half.minimum_enthalpy,
      half.minimum_energy};
  std::array<double, 10U> global_minima{};
  std::array<double, 12U> local_endpoint_maxima{
      full.maximum_pressure, full.maximum_density, full.maximum_temperature,
      full.maximum_enthalpy, full.maximum_energy, full.maximum_mach,
      half.maximum_pressure, half.maximum_density, half.maximum_temperature,
      half.maximum_enthalpy, half.maximum_energy, half.maximum_mach};
  std::array<double, 12U> global_endpoint_maxima{};
  std::array<std::uint64_t, 12U> local_counts{
      full.nonfinite_cells, full.nonpositive_pressure,
      full.nonpositive_density, full.nonpositive_temperature,
      full.mach_limit_violations, full.solid_velocity_nonpositive_zero,
      half.nonfinite_cells, half.nonpositive_pressure,
      half.nonpositive_density, half.nonpositive_temperature,
      half.mach_limit_violations, half.solid_velocity_nonpositive_zero};
  std::array<std::uint64_t, 12U> global_counts{};
  const bool okay =
      MPI_Allreduce(local_sums.data(), global_sums.data(), local_sums.size(),
                    MPI_DOUBLE, MPI_SUM, communicator) == MPI_SUCCESS &&
      MPI_Allreduce(local_maxima.data(), global_maxima.data(),
                    local_maxima.size(), MPI_DOUBLE, MPI_MAX, communicator) ==
          MPI_SUCCESS &&
      MPI_Allreduce(local_integrals.data(), global_integrals.data(),
                    local_integrals.size(), MPI_DOUBLE, MPI_SUM,
                    communicator) == MPI_SUCCESS &&
      MPI_Allreduce(local_minima.data(), global_minima.data(),
                    local_minima.size(), MPI_DOUBLE, MPI_MIN, communicator) ==
          MPI_SUCCESS &&
      MPI_Allreduce(local_endpoint_maxima.data(),
                    global_endpoint_maxima.data(),
                    local_endpoint_maxima.size(), MPI_DOUBLE, MPI_MAX,
                    communicator) == MPI_SUCCESS &&
      MPI_Allreduce(local_counts.data(), global_counts.data(),
                    local_counts.size(), MPI_UINT64_T, MPI_SUM,
                    communicator) == MPI_SUCCESS;
  if (!okay) return {StatusCode::mpi_failure, kComparatorNumeric};
  for (std::size_t index = 0U; index < count; ++index) {
    ErrorMetric& metric = metrics[index];
    metric.weight = global_sums[index * 4U];
    metric.difference_squared = global_sums[index * 4U + 1U];
    metric.full_squared = global_sums[index * 4U + 2U];
    metric.half_squared = global_sums[index * 4U + 3U];
    metric.difference_linf = global_maxima[index * 3U];
    metric.full_linf = global_maxima[index * 3U + 1U];
    metric.half_linf = global_maxima[index * 3U + 2U];
    if (metric.weight > 0.0) {
      metric.rms = std::sqrt(metric.difference_squared / metric.weight);
      const double full_rms = std::sqrt(metric.full_squared / metric.weight);
      const double half_rms = std::sqrt(metric.half_squared / metric.weight);
      metric.relative_rms = relative_error(metric.rms, full_rms, half_rms);
    }
    metric.relative_linf =
        relative_error(metric.difference_linf, metric.full_linf,
                       metric.half_linf);
  }
  full.volume = global_integrals[0U];
  full.mass = global_integrals[1U];
  full.energy = global_integrals[2U];
  full.kinetic = global_integrals[3U];
  half.volume = global_integrals[4U];
  half.mass = global_integrals[5U];
  half.energy = global_integrals[6U];
  half.kinetic = global_integrals[7U];
  full.minimum_pressure = global_minima[0U];
  full.minimum_density = global_minima[1U];
  full.minimum_temperature = global_minima[2U];
  full.minimum_enthalpy = global_minima[3U];
  full.minimum_energy = global_minima[4U];
  half.minimum_pressure = global_minima[5U];
  half.minimum_density = global_minima[6U];
  half.minimum_temperature = global_minima[7U];
  half.minimum_enthalpy = global_minima[8U];
  half.minimum_energy = global_minima[9U];
  full.maximum_pressure = global_endpoint_maxima[0U];
  full.maximum_density = global_endpoint_maxima[1U];
  full.maximum_temperature = global_endpoint_maxima[2U];
  full.maximum_enthalpy = global_endpoint_maxima[3U];
  full.maximum_energy = global_endpoint_maxima[4U];
  full.maximum_mach = global_endpoint_maxima[5U];
  half.maximum_pressure = global_endpoint_maxima[6U];
  half.maximum_density = global_endpoint_maxima[7U];
  half.maximum_temperature = global_endpoint_maxima[8U];
  half.maximum_enthalpy = global_endpoint_maxima[9U];
  half.maximum_energy = global_endpoint_maxima[10U];
  half.maximum_mach = global_endpoint_maxima[11U];
  full.nonfinite_cells = global_counts[0U];
  full.nonpositive_pressure = global_counts[1U];
  full.nonpositive_density = global_counts[2U];
  full.nonpositive_temperature = global_counts[3U];
  full.mach_limit_violations = global_counts[4U];
  full.solid_velocity_nonpositive_zero = global_counts[5U];
  half.nonfinite_cells = global_counts[6U];
  half.nonpositive_pressure = global_counts[7U];
  half.nonpositive_density = global_counts[8U];
  half.nonpositive_temperature = global_counts[9U];
  half.mach_limit_violations = global_counts[10U];
  half.solid_velocity_nonpositive_zero = global_counts[11U];
  return {};
}

Status compare_cells(
    MPI_Comm communicator, const EndpointRaw& full, const EndpointRaw& half,
    const PhysicalAuthority& authority,
    std::array<ErrorMetric, static_cast<std::size_t>(MetricId::count)>& metrics,
    EndpointIntegral& full_integral, EndpointIntegral& half_integral,
    double mach_limit) noexcept try {
  const RawField& full_u = *find_field(full, RestartFieldRole::velocity);
  const RawField& half_u = *find_field(half, RestartFieldRole::velocity);
  const RawField& full_pi =
      *find_field(full, RestartFieldRole::pressure_perturbation);
  const RawField& half_pi =
      *find_field(half, RestartFieldRole::pressure_perturbation);
  const RawField& full_h = *find_field(full, RestartFieldRole::enthalpy);
  const RawField& half_h = *find_field(half, RestartFieldRole::enthalpy);
  std::vector<const RawField*> full_species;
  std::vector<const RawField*> half_species;
  for (std::size_t index = 0U; index < full.fields.size(); ++index)
    if (full.fields[index].role == RestartFieldRole::independent_species) {
      full_species.push_back(&full.fields[index]);
      half_species.push_back(&half.fields[index]);
    }
  std::vector<double> full_y(full_species.size());
  std::vector<double> half_y(half_species.size());
  const Span<const std::uint8_t> region = authority.topology.has_value()
                                               ? authority.topology->region()
                                               : Span<const std::uint8_t>{};
  std::size_t cells = 0U;
  checked_cells(full.patch.cells, cells);
  Status local{};
  for (std::int32_t z = 0; z < full.patch.cells.z && local; ++z)
    for (std::int32_t y = 0; y < full.patch.cells.y && local; ++y)
      for (std::int32_t x = 0; x < full.patch.cells.x && local; ++x) {
        const Int3 local_index{x, y, z};
        const std::size_t cell = cell_offset(full.patch.cells, local_index);
        const bool fluid = !authority.topology.has_value() ||
                           region.data[cell] ==
                               static_cast<std::uint8_t>(RegionFlag::fluid);
        const Real3 full_velocity = field_vector(full_u, cell);
        const Real3 half_velocity = field_vector(half_u, cell);
        if (!fluid) {
          const auto positive_zero = [](Real3 value) noexcept {
            return bits(value.x) == 0U && bits(value.y) == 0U &&
                   bits(value.z) == 0U;
          };
          full_integral.solid_velocity_nonpositive_zero +=
              positive_zero(full_velocity) ? 0U : 1U;
          half_integral.solid_velocity_nonpositive_zero +=
              positive_zero(half_velocity) ? 0U : 1U;
          continue;
        }
        for (std::size_t species = 0U; species < full_species.size(); ++species) {
          full_y[species] = field_scalar(*full_species[species], cell);
          half_y[species] = field_scalar(*half_species[species], cell);
        }
        const double full_pi_value = field_scalar(full_pi, cell);
        const double half_pi_value = field_scalar(half_pi, cell);
        const double full_p = full.pressure_reference + full_pi_value;
        const double half_p = half.pressure_reference + half_pi_value;
        const double full_h_value = field_scalar(full_h, cell);
        const double half_h_value = field_scalar(half_h, cell);
        ThermoState full_thermo;
        ThermoState half_thermo;
        local = authority.thermodynamics.evaluate_from_reference_pressure(
            full.pressure_reference, full_pi_value, full_h_value,
            {full_y.data(), full_y.size()}, full_velocity, full_thermo);
        if (local)
          local = authority.thermodynamics.evaluate_from_reference_pressure(
              half.pressure_reference, half_pi_value, half_h_value,
              {half_y.data(), half_y.size()}, half_velocity, half_thermo);
        if (!local) break;
        const double full_q = full_thermo.rho * full_h_value - full_p;
        const double half_q = half_thermo.rho * half_h_value - half_p;
        const Real3 full_momentum{full_thermo.rho * full_velocity.x,
                                  full_thermo.rho * full_velocity.y,
                                  full_thermo.rho * full_velocity.z};
        const Real3 half_momentum{half_thermo.rho * half_velocity.x,
                                  half_thermo.rho * half_velocity.y,
                                  half_thermo.rho * half_velocity.z};
        const Int3 global{full.patch.begin.x + x, full.patch.begin.y + y,
                          full.patch.begin.z + z};
        const double volume =
            authority.geometry.x().widths().data[global.x] *
            authority.geometry.y().widths().data[global.y] *
            authority.geometry.z().widths().data[global.z];
        add_vector(metrics[static_cast<std::size_t>(MetricId::velocity)],
                   volume, full_velocity, half_velocity);
        add_scalar(metrics[static_cast<std::size_t>(
                       MetricId::pressure_perturbation)],
                   volume, full_pi_value, half_pi_value);
        add_scalar(metrics[static_cast<std::size_t>(
                       MetricId::pressure_absolute)],
                   volume, full_p, half_p);
        add_scalar(metrics[static_cast<std::size_t>(MetricId::enthalpy)],
                   volume, full_h_value, half_h_value);
        add_scalar(metrics[static_cast<std::size_t>(MetricId::density)],
                   volume, full_thermo.rho, half_thermo.rho);
        add_scalar(metrics[static_cast<std::size_t>(MetricId::temperature)],
                   volume, full_thermo.temperature, half_thermo.temperature);
        add_scalar(metrics[static_cast<std::size_t>(
                       MetricId::sensible_energy)],
                   volume, full_q, half_q);
        add_vector(metrics[static_cast<std::size_t>(MetricId::momentum)],
                   volume, full_momentum, half_momentum);
        update_integral(full_integral, volume, full_p, full_h_value, full_q,
                        full_velocity, full_thermo, mach_limit);
        update_integral(half_integral, volume, half_p, half_h_value, half_q,
                        half_velocity, half_thermo, mach_limit);
      }
  const Status status = collective_status(communicator, local);
  if (!status) return status;
  return reduce_cell_results(communicator, metrics, full_integral,
                             half_integral);
} catch (const std::bad_alloc&) {
  return collective_status(
      communicator, {StatusCode::allocation_failure, kComparatorNumeric});
} catch (...) {
  return collective_status(communicator,
                           {StatusCode::numerical_failure, kComparatorNumeric});
}

bool inlet_kind(BoundaryKind kind) noexcept {
  return kind == BoundaryKind::velocity_inlet ||
         kind == BoundaryKind::mass_flow_inlet ||
         kind == BoundaryKind::static_state_inlet ||
         kind == BoundaryKind::total_state_inlet ||
         kind == BoundaryKind::nscbc_inlet;
}

bool outlet_kind(BoundaryKind kind) noexcept {
  return kind == BoundaryKind::pressure_outlet ||
         kind == BoundaryKind::nscbc_outlet;
}

std::size_t boundary_index(std::size_t axis, bool maximum) noexcept {
  return axis * 2U + (maximum ? 1U : 0U);
}

void update_boundary(BoundaryFaceIntegral& result,
                     const BoundaryFaceSpec& spec, double outward,
                     double raw_flux) noexcept {
  result.net_outward += outward;
  result.inflow += std::max(-outward, 0.0);
  result.outflow += std::max(outward, 0.0);
  result.maximum_absolute_flux =
      std::max(result.maximum_absolute_flux, std::abs(raw_flux));
  const double sign_tolerance =
      64.0 * std::numeric_limits<double>::epsilon() *
      std::max(1.0, std::abs(raw_flux));
  if (inlet_kind(spec.flow_kind) && outward > sign_tolerance)
    ++result.inlet_reversal_faces;
  if (outlet_kind(spec.flow_kind) && outward < -sign_tolerance) {
    ++result.outlet_backflow_faces;
    if (!spec.allow_backflow) ++result.prohibited_backflow_faces;
  }
  const bool impermeable =
      spec.flow_kind == BoundaryKind::no_slip_wall ||
      spec.flow_kind == BoundaryKind::moving_wall ||
      spec.flow_kind == BoundaryKind::slip ||
      spec.flow_kind == BoundaryKind::symmetry;
  if (impermeable && bits(raw_flux) != 0U) {
    ++result.symmetry_nonpositive_zero;
    if (bits(raw_flux) == UINT64_C(0x8000000000000000))
      ++result.symmetry_negative_zero;
  }
}

double face_area(const CartesianGeometryPlan& geometry, std::size_t axis,
                 Int3 global) noexcept {
  const Span<const double> dx = geometry.x().widths();
  const Span<const double> dy = geometry.y().widths();
  const Span<const double> dz = geometry.z().widths();
  if (axis == 0U) return dy.data[global.y] * dz.data[global.z];
  if (axis == 1U) return dx.data[global.x] * dz.data[global.z];
  return dx.data[global.x] * dy.data[global.y];
}

Status configured_velocity_inlet_composition(
    const EndpointRaw& endpoint, const BoundaryFaceSpec& boundary,
    const ThermodynamicsPlan& thermodynamics,
    std::vector<double>& composition) {
  composition.clear();
  composition.reserve(thermodynamics.independent_species_count());
  for (const TransportedScalarSpec& declared :
       endpoint.model.transported_scalars) {
    if (declared.role != TransportedScalarRole::species) continue;
    const auto configured = std::find_if(
        boundary.scalars.begin(), boundary.scalars.end(),
        [&](const ScalarBoundarySpec& scalar) {
          return scalar.stable_name == declared.stable_name;
        });
    if (configured == boundary.scalars.end() ||
        configured->kind != ScalarBoundaryKind::dirichlet)
      return {StatusCode::invalid_plan, kComparatorNumeric};
    composition.push_back(configured->value);
  }
  if (composition.size() != thermodynamics.independent_species_count())
    return {StatusCode::invalid_plan, kComparatorNumeric};
  return {};
}

Int3 boundary_owner(std::size_t axis, bool maximum, Int3 face,
                    Int3 cells) noexcept {
  Int3 owner = face;
  if (axis == 0U)
    owner.x = maximum ? cells.x - 1 : 0;
  else if (axis == 1U)
    owner.y = maximum ? cells.y - 1 : 0;
  else
    owner.z = maximum ? cells.z - 1 : 0;
  return owner;
}

double axis_component(Real3 value, std::size_t axis) noexcept {
  return axis == 0U ? value.x : (axis == 1U ? value.y : value.z);
}

Status expected_velocity_inlet_flux(
    const EndpointRaw& endpoint, const PhysicalAuthority& authority,
    std::size_t boundary, std::size_t axis, bool maximum, Int3 face,
    double area, Span<const double> composition, double& expected) noexcept {
  if (boundary >= endpoint.model.boundaries.size() || !(area > 0.0) ||
      !std::isfinite(area))
    return {StatusCode::invalid_plan, kComparatorNumeric};
  const BoundaryFaceSpec& spec = endpoint.model.boundaries[boundary];
  const RawField* pressure =
      find_field(endpoint, RestartFieldRole::pressure_perturbation);
  if (spec.flow_kind != BoundaryKind::velocity_inlet || pressure == nullptr ||
      pressure->components != 1U || !(spec.temperature > 0.0) ||
      !std::isfinite(spec.temperature))
    return {StatusCode::invalid_plan, kComparatorNumeric};
  const Int3 owner = boundary_owner(axis, maximum, face, endpoint.patch.cells);
  const std::size_t cell = cell_offset(endpoint.patch.cells, owner);
  double enthalpy = 0.0;
  double cp = 0.0;
  double gas_constant = 0.0;
  Status status = authority.thermodynamics.mixture_enthalpy(
      spec.temperature, composition, enthalpy, cp, gas_constant);
  ThermoState thermo;
  if (status)
    status = authority.thermodynamics.evaluate_from_reference_pressure(
        endpoint.pressure_reference, field_scalar(*pressure, cell), enthalpy,
        composition, spec.velocity, thermo, spec.temperature);
  if (!status) return status;
  expected = thermo.rho * axis_component(spec.velocity, axis) * area;
  if (!std::isfinite(expected))
    return {StatusCode::numerical_failure, kComparatorNumeric};
  return {};
}

bool velocity_inlet_flux_matches(const EndpointRaw& endpoint, double actual,
                                 double expected) noexcept {
  const double relative_tolerance = std::max(
      {256.0 * std::numeric_limits<double>::epsilon(),
       endpoint.model.solver.terminal.eos,
       endpoint.model.solver.terminal.continuity});
  const double scale =
      std::max({1.0, std::abs(actual), std::abs(expected)});
  return std::isfinite(actual) && std::isfinite(expected) &&
         std::isfinite(relative_tolerance) && relative_tolerance >= 0.0 &&
         std::abs(actual - expected) <= relative_tolerance * scale;
}

double link_flux(const EndpointRaw& endpoint,
                 const ImmersedLink& link) noexcept {
  Int3 face = link.fluid_local_index;
  std::size_t axis = 0U;
  switch (link.direction) {
    case ImmersedFaceDirection::x_negative:
      axis = 0U;
      break;
    case ImmersedFaceDirection::x_positive:
      axis = 0U;
      ++face.x;
      break;
    case ImmersedFaceDirection::y_negative:
      axis = 1U;
      break;
    case ImmersedFaceDirection::y_positive:
      axis = 1U;
      ++face.y;
      break;
    case ImmersedFaceDirection::z_negative:
      axis = 2U;
      break;
    case ImmersedFaceDirection::z_positive:
      axis = 2U;
      ++face.z;
      break;
  }
  const std::array<Int3, 3U> extents{
      Int3{endpoint.patch.cells.x + 1, endpoint.patch.cells.y,
           endpoint.patch.cells.z},
      Int3{endpoint.patch.cells.x, endpoint.patch.cells.y + 1,
           endpoint.patch.cells.z},
      Int3{endpoint.patch.cells.x, endpoint.patch.cells.y,
           endpoint.patch.cells.z + 1}};
  return endpoint.final_mass_flux[axis][face_offset(extents[axis], face)];
}

Status reduce_face_results(MPI_Comm communicator, ErrorMetric& flux,
                           ErrorMetric& flux_density,
                           EndpointFaceIntegral& full,
                           EndpointFaceIntegral& half) noexcept {
  std::array<double, 8U> local_metric_sums{
      flux.weight, flux.difference_squared, flux.full_squared,
      flux.half_squared, flux_density.weight,
      flux_density.difference_squared, flux_density.full_squared,
      flux_density.half_squared};
  std::array<double, 8U> global_metric_sums{};
  std::array<double, 6U> local_metric_maxima{
      flux.difference_linf, flux.full_linf, flux.half_linf,
      flux_density.difference_linf, flux_density.full_linf,
      flux_density.half_linf};
  std::array<double, 6U> global_metric_maxima{};
  std::array<double, 36U> local_boundary_sums{};
  std::array<double, 36U> global_boundary_sums{};
  std::array<double, 12U> local_boundary_maxima{};
  std::array<double, 12U> global_boundary_maxima{};
  std::array<std::uint64_t, 80U> local_counts{};
  std::array<std::uint64_t, 80U> global_counts{};
  for (std::size_t endpoint = 0U; endpoint < 2U; ++endpoint) {
    const EndpointFaceIntegral& source = endpoint == 0U ? full : half;
    for (std::size_t face = 0U; face < 6U; ++face) {
      const BoundaryFaceIntegral& boundary = source.boundary[face];
      const std::size_t sum = endpoint * 18U + face * 3U;
      local_boundary_sums[sum] = boundary.net_outward;
      local_boundary_sums[sum + 1U] = boundary.inflow;
      local_boundary_sums[sum + 2U] = boundary.outflow;
      local_boundary_maxima[endpoint * 6U + face] =
          boundary.maximum_absolute_flux;
      const std::size_t count = endpoint * 40U + face * 6U;
      local_counts[count] = boundary.inlet_reversal_faces;
      local_counts[count + 1U] =
          boundary.velocity_inlet_flux_mismatch_faces;
      local_counts[count + 2U] = boundary.outlet_backflow_faces;
      local_counts[count + 3U] = boundary.prohibited_backflow_faces;
      local_counts[count + 4U] = boundary.symmetry_nonpositive_zero;
      local_counts[count + 5U] = boundary.symmetry_negative_zero;
    }
    const std::size_t tail = endpoint * 40U + 36U;
    local_counts[tail] = source.nonfinite_faces;
    local_counts[tail + 1U] = source.ibm_nonpositive_zero;
  }
  // Keep the exact IBM zero-sign diagnostic separate from boundary-kind
  // counters so the latter remain one compact fixed-stride block.
  std::array<std::uint64_t, 2U> local_ibm_negative{
      full.ibm_negative_zero, half.ibm_negative_zero};
  std::array<std::uint64_t, 2U> global_ibm_negative{};
  const bool okay =
      MPI_Allreduce(local_metric_sums.data(), global_metric_sums.data(),
                    local_metric_sums.size(), MPI_DOUBLE, MPI_SUM,
                    communicator) == MPI_SUCCESS &&
      MPI_Allreduce(local_metric_maxima.data(), global_metric_maxima.data(),
                    local_metric_maxima.size(), MPI_DOUBLE, MPI_MAX,
                    communicator) == MPI_SUCCESS &&
      MPI_Allreduce(local_boundary_sums.data(), global_boundary_sums.data(),
                    local_boundary_sums.size(), MPI_DOUBLE, MPI_SUM,
                    communicator) == MPI_SUCCESS &&
      MPI_Allreduce(local_boundary_maxima.data(),
                    global_boundary_maxima.data(),
                    local_boundary_maxima.size(), MPI_DOUBLE, MPI_MAX,
                    communicator) == MPI_SUCCESS &&
      MPI_Allreduce(local_counts.data(), global_counts.data(),
                    local_counts.size(), MPI_UINT64_T, MPI_SUM,
                    communicator) == MPI_SUCCESS &&
      MPI_Allreduce(local_ibm_negative.data(), global_ibm_negative.data(),
                    local_ibm_negative.size(), MPI_UINT64_T, MPI_SUM,
                    communicator) == MPI_SUCCESS;
  if (!okay) return {StatusCode::mpi_failure, kComparatorNumeric};
  const auto assign_metric = [&](ErrorMetric& metric, std::size_t sum,
                                 std::size_t maximum) {
    metric.weight = global_metric_sums[sum];
    metric.difference_squared = global_metric_sums[sum + 1U];
    metric.full_squared = global_metric_sums[sum + 2U];
    metric.half_squared = global_metric_sums[sum + 3U];
    metric.difference_linf = global_metric_maxima[maximum];
    metric.full_linf = global_metric_maxima[maximum + 1U];
    metric.half_linf = global_metric_maxima[maximum + 2U];
    if (metric.weight > 0.0) {
      metric.rms = std::sqrt(metric.difference_squared / metric.weight);
      metric.relative_rms = relative_error(
          metric.rms, std::sqrt(metric.full_squared / metric.weight),
          std::sqrt(metric.half_squared / metric.weight));
    }
    metric.relative_linf = relative_error(
        metric.difference_linf, metric.full_linf, metric.half_linf);
  };
  assign_metric(flux, 0U, 0U);
  assign_metric(flux_density, 4U, 3U);
  for (std::size_t endpoint = 0U; endpoint < 2U; ++endpoint) {
    EndpointFaceIntegral& destination = endpoint == 0U ? full : half;
    for (std::size_t face = 0U; face < 6U; ++face) {
      BoundaryFaceIntegral& boundary = destination.boundary[face];
      const std::size_t sum = endpoint * 18U + face * 3U;
      boundary.net_outward = global_boundary_sums[sum];
      boundary.inflow = global_boundary_sums[sum + 1U];
      boundary.outflow = global_boundary_sums[sum + 2U];
      boundary.maximum_absolute_flux =
          global_boundary_maxima[endpoint * 6U + face];
      const std::size_t count = endpoint * 40U + face * 6U;
      boundary.inlet_reversal_faces = global_counts[count];
      boundary.velocity_inlet_flux_mismatch_faces =
          global_counts[count + 1U];
      boundary.outlet_backflow_faces = global_counts[count + 2U];
      boundary.prohibited_backflow_faces = global_counts[count + 3U];
      boundary.symmetry_nonpositive_zero = global_counts[count + 4U];
      boundary.symmetry_negative_zero = global_counts[count + 5U];
    }
    const std::size_t tail = endpoint * 40U + 36U;
    destination.nonfinite_faces = global_counts[tail];
    destination.ibm_nonpositive_zero = global_counts[tail + 1U];
    destination.ibm_negative_zero = global_ibm_negative[endpoint];
  }
  return {};
}

Status compare_faces(
    MPI_Comm communicator, const EndpointRaw& full, const EndpointRaw& half,
    const PhysicalAuthority& authority, ErrorMetric& flux,
    ErrorMetric& flux_density, EndpointFaceIntegral& full_integral,
    EndpointFaceIntegral& half_integral) noexcept try {
  std::array<std::vector<double>, 6U> full_compositions;
  std::array<std::vector<double>, 6U> half_compositions;
  Status local_status{};
  for (std::size_t face = 0U; face < 6U && local_status; ++face) {
    if (full.model.boundaries[face].flow_kind ==
        BoundaryKind::velocity_inlet)
      local_status = configured_velocity_inlet_composition(
          full, full.model.boundaries[face], authority.thermodynamics,
          full_compositions[face]);
    if (local_status && half.model.boundaries[face].flow_kind ==
                            BoundaryKind::velocity_inlet)
      local_status = configured_velocity_inlet_composition(
          half, half.model.boundaries[face], authority.thermodynamics,
          half_compositions[face]);
  }
  const std::array<Int3, 3U> extents{
      Int3{full.patch.cells.x + 1, full.patch.cells.y, full.patch.cells.z},
      Int3{full.patch.cells.x, full.patch.cells.y + 1, full.patch.cells.z},
      Int3{full.patch.cells.x, full.patch.cells.y, full.patch.cells.z + 1}};
  const std::array<std::int32_t, 3U> patch_begin{
      full.patch.begin.x, full.patch.begin.y, full.patch.begin.z};
  const std::array<std::int32_t, 3U> global_axis_cells{
      full.global_cells.x, full.global_cells.y, full.global_cells.z};
  for (std::size_t axis = 0U; axis < 3U && local_status; ++axis) {
    for (std::int32_t z = 0; z < extents[axis].z && local_status; ++z)
      for (std::int32_t y = 0; y < extents[axis].y && local_status; ++y)
        for (std::int32_t x = 0; x < extents[axis].x; ++x) {
          const Int3 local{x, y, z};
          const std::array<std::int32_t, 3U> coordinate{x, y, z};
          const std::int32_t face_coordinate = coordinate[axis];
          // The lower duplicate of an MPI-shared face belongs to the rank on
          // its negative side; the global minimum face is owned explicitly.
          if (face_coordinate == 0 && patch_begin[axis] != 0) continue;
          Int3 global{full.patch.begin.x + x, full.patch.begin.y + y,
                      full.patch.begin.z + z};
          // Along the selected axis this is a face index; transverse entries
          // remain cell indices and therefore address axis widths directly.
          const double area = face_area(authority.geometry, axis, global);
          const std::size_t offset = face_offset(extents[axis], local);
          const double full_value = full.final_mass_flux[axis][offset];
          const double half_value = half.final_mass_flux[axis][offset];
          full_integral.nonfinite_faces += std::isfinite(full_value) ? 0U : 1U;
          half_integral.nonfinite_faces += std::isfinite(half_value) ? 0U : 1U;
          add_scalar(flux, area, full_value, half_value);
          add_scalar(flux_density, area, full_value / area, half_value / area);
          const std::int32_t global_face = patch_begin[axis] + face_coordinate;
          const bool minimum = global_face == 0;
          const bool maximum = global_face == global_axis_cells[axis];
          if (minimum || maximum) {
            const std::size_t face = boundary_index(axis, maximum);
            const double sign = minimum ? -1.0 : 1.0;
            update_boundary(full_integral.boundary[face],
                            full.model.boundaries[face], sign * full_value,
                            full_value);
            update_boundary(half_integral.boundary[face],
                            half.model.boundaries[face], sign * half_value,
                            half_value);
            if (full.model.boundaries[face].flow_kind ==
                BoundaryKind::velocity_inlet) {
              double expected = 0.0;
              local_status = expected_velocity_inlet_flux(
                  full, authority, face, axis, maximum, local, area,
                  {full_compositions[face].data(),
                   full_compositions[face].size()},
                  expected);
              if (local_status &&
                  !velocity_inlet_flux_matches(full, full_value, expected))
                ++full_integral.boundary[face]
                      .velocity_inlet_flux_mismatch_faces;
            }
            if (local_status && half.model.boundaries[face].flow_kind ==
                                    BoundaryKind::velocity_inlet) {
              double expected = 0.0;
              local_status = expected_velocity_inlet_flux(
                  half, authority, face, axis, maximum, local, area,
                  {half_compositions[face].data(),
                   half_compositions[face].size()},
                  expected);
              if (local_status &&
                  !velocity_inlet_flux_matches(half, half_value, expected))
                ++half_integral.boundary[face]
                      .velocity_inlet_flux_mismatch_faces;
            }
            if (!local_status) break;
          }
        }
  }
  if (authority.topology.has_value()) {
    const Span<const ImmersedLink> links = authority.topology->links();
    for (std::size_t index = 0U; index < links.size; ++index) {
      const ImmersedLink& link = links.data[index];
      const double full_value = link_flux(full, link);
      const double half_value = link_flux(half, link);
      if (bits(full_value) != 0U) {
        ++full_integral.ibm_nonpositive_zero;
        if (bits(full_value) == UINT64_C(0x8000000000000000))
          ++full_integral.ibm_negative_zero;
      }
      if (bits(half_value) != 0U) {
        ++half_integral.ibm_nonpositive_zero;
        if (bits(half_value) == UINT64_C(0x8000000000000000))
          ++half_integral.ibm_negative_zero;
      }
    }
  }
  const Status agreement = collective_status(communicator, local_status);
  if (!agreement) return agreement;
  return reduce_face_results(communicator, flux, flux_density, full_integral,
                             half_integral);
} catch (const std::bad_alloc&) {
  return collective_status(
      communicator, {StatusCode::allocation_failure, kComparatorNumeric});
} catch (...) {
  return collective_status(communicator,
                           {StatusCode::numerical_failure, kComparatorNumeric});
}

bool close_roundoff(double left, double right, double multiplier = 256.0) {
  return std::abs(left - right) <=
         multiplier * std::numeric_limits<double>::epsilon() *
             std::max({1.0, std::abs(left), std::abs(right)});
}

TimingIdentity endpoint_timing_identity(ComparisonMode mode,
                                        EndpointTiming first,
                                        EndpointTiming second) noexcept {
  if (mode == ComparisonMode::equivalent_endpoints) {
    return {same_bits(first.initial_dt, second.initial_dt) &&
                same_bits(first.minimum_dt, second.minimum_dt) &&
                same_bits(first.maximum_dt, second.maximum_dt),
            same_bits(first.accepted_dt, second.accepted_dt),
            first.step == second.step};
  }
  const auto half_of = [](double coarse, double fine) noexcept {
    return close_roundoff(coarse, 2.0 * fine, 1024.0);
  };
  return {half_of(first.initial_dt, second.initial_dt) &&
              half_of(first.minimum_dt, second.minimum_dt) &&
              half_of(first.maximum_dt, second.maximum_dt),
          half_of(first.accepted_dt, second.accepted_dt),
          first.step <= std::numeric_limits<std::uint64_t>::max() / 2U &&
              second.step == first.step * 2U};
}

double scalar_relative_difference(double left, double right) noexcept {
  return relative_error(std::abs(right - left), std::abs(left),
                        std::abs(right));
}

void append_failure(ComparisonReport& report, bool failure,
                    std::string message) {
  if (failure) report.failures.push_back(std::move(message));
}

void certify_report(const EndpointRaw& full, const EndpointRaw& half,
                    const Options& options, ComparisonReport& report) {
  for (std::size_t index = 0U; index < report.metrics.size(); ++index) {
    const ErrorMetric& metric = report.metrics[index];
    append_failure(report,
                   !std::isfinite(metric.rms) ||
                       !std::isfinite(metric.difference_linf),
                   std::string(kMetricNames[index]) +
                       " comparison is non-finite");
    append_failure(report,
                   metric.relative_rms > *options.maximum_relative_rms,
                   std::string(kMetricNames[index]) +
                       " relative RMS exceeds --max-relative-rms");
    append_failure(report,
                   metric.relative_linf > *options.maximum_relative_linf,
                   std::string(kMetricNames[index]) +
                       " relative Linf exceeds --max-relative-linf");
  }
  const auto endpoint_physics = [&](const EndpointIntegral& endpoint,
                                    std::string_view name) {
    append_failure(report, endpoint.nonfinite_cells != 0U,
                   std::string(name) + " has non-finite fluid cells");
    append_failure(report, endpoint.nonpositive_pressure != 0U,
                   std::string(name) + " has non-positive absolute pressure");
    append_failure(report, endpoint.nonpositive_density != 0U,
                   std::string(name) + " has non-positive density");
    append_failure(report, endpoint.nonpositive_temperature != 0U,
                   std::string(name) + " has non-positive temperature");
    append_failure(report, endpoint.mach_limit_violations != 0U,
                   std::string(name) + " exceeds configured Mach limit");
    append_failure(report, endpoint.solid_velocity_nonpositive_zero != 0U,
                   std::string(name) +
                       " IBM solid velocity is not exact positive zero");
    append_failure(report,
                   !std::isfinite(endpoint.volume) ||
                       !std::isfinite(endpoint.mass) ||
                       !std::isfinite(endpoint.energy) ||
                       !std::isfinite(endpoint.kinetic) ||
                       !(endpoint.volume > 0.0) || !(endpoint.mass > 0.0),
                   std::string(name) + " global integral is inadmissible");
  };
  endpoint_physics(report.full, "full");
  endpoint_physics(report.half, "half");
  append_failure(report,
                 scalar_relative_difference(report.full.mass,
                                            report.half.mass) >
                     *options.maximum_relative_linf,
                 "global mass difference exceeds --max-relative-linf");
  append_failure(report,
                 scalar_relative_difference(report.full.energy,
                                            report.half.energy) >
                     *options.maximum_relative_linf,
                 "global energy difference exceeds --max-relative-linf");
  append_failure(report,
                 scalar_relative_difference(report.full.kinetic,
                                            report.half.kinetic) >
                     *options.maximum_relative_linf,
                 "global kinetic-energy difference exceeds "
                 "--max-relative-linf");
  append_failure(report,
                 scalar_relative_difference(
                     report.full.energy + report.full.kinetic,
                     report.half.energy + report.half.kinetic) >
                     *options.maximum_relative_linf,
                 "global total-energy difference exceeds "
                 "--max-relative-linf");
  const auto face_physics = [&](const EndpointRaw& endpoint,
                                const EndpointFaceIntegral& faces,
                                std::string_view name) {
    append_failure(report, faces.nonfinite_faces != 0U,
                   std::string(name) + " has non-finite final face flux");
    append_failure(report, faces.ibm_nonpositive_zero != 0U,
                   std::string(name) +
                       " IBM interface flux is not exact positive zero");
    for (std::size_t face = 0U; face < faces.boundary.size(); ++face) {
      const BoundaryFaceIntegral& value = faces.boundary[face];
      append_failure(report, value.inlet_reversal_faces != 0U,
                     std::string(name) + " inlet has reversed final flux");
      append_failure(
          report, value.velocity_inlet_flux_mismatch_faces != 0U,
          std::string(name) +
              " velocity inlet final flux misses its EOS closure");
      append_failure(report, value.prohibited_backflow_faces != 0U,
                     std::string(name) +
                         " outlet has prohibited final backflow");
      append_failure(report, value.symmetry_nonpositive_zero != 0U,
                     std::string(name) +
                         " impermeable boundary flux is not exact positive "
                         "zero");
      const BoundaryFaceSpec& spec = endpoint.model.boundaries[face];
      if (spec.flow_kind == BoundaryKind::mass_flow_inlet) {
        const double target = std::abs(spec.mass_flow_rate);
        const double tolerance =
            std::max(256.0 * std::numeric_limits<double>::epsilon(),
                     endpoint.model.solver.terminal.continuity) *
            std::max(1.0, target);
        append_failure(report, std::abs(value.inflow - target) > tolerance,
                       std::string(name) +
                           " mass-flow inlet misses its final target");
      }
    }
  };
  face_physics(full, report.full_faces, "full");
  face_physics(half, report.half_faces, "half");
}

std::string json_escape(std::string_view input) {
  std::string output;
  output.reserve(input.size() + 8U);
  for (const char character : input) {
    const auto value = static_cast<unsigned char>(character);
    switch (value) {
      case '"':
        output += "\\\"";
        break;
      case '\\':
        output += "\\\\";
        break;
      case '\b':
        output += "\\b";
        break;
      case '\f':
        output += "\\f";
        break;
      case '\n':
        output += "\\n";
        break;
      case '\r':
        output += "\\r";
        break;
      case '\t':
        output += "\\t";
        break;
      default:
        if (value < 0x20U) {
          std::ostringstream escaped;
          escaped << "\\u" << std::hex << std::setw(4) << std::setfill('0')
                  << static_cast<unsigned>(value);
          output += escaped.str();
        } else {
          output.push_back(static_cast<char>(value));
        }
    }
  }
  return output;
}

void json_number(std::ostream& stream, double value) {
  if (std::isfinite(value))
    stream << std::setprecision(17) << value;
  else
    stream << "null";
}

void json_fingerprint(std::ostream& stream, PlanFingerprint value) {
  stream << '"' << "0x" << std::hex << std::setw(16) << std::setfill('0')
         << value << std::dec << std::setfill(' ') << '"';
}

void write_endpoint_integral(std::ostream& stream,
                             const EndpointRaw& endpoint,
                             const EndpointIntegral& integral,
                             const EndpointFaceIntegral& faces) {
  stream << "{\"case_fingerprint\":";
  json_fingerprint(stream, endpoint.case_fingerprint);
  stream << ",\"product_fingerprint\":";
  json_fingerprint(stream, endpoint.product_fingerprint);
  stream << ",\"time\":";
  json_number(stream, endpoint.time);
  stream << ",\"dt\":";
  json_number(stream, endpoint.dt);
  stream << ",\"step\":" << endpoint.step
         << ",\"pressure_reference\":";
  json_number(stream, endpoint.pressure_reference);
  stream << ",\"integrals\":{\"volume\":";
  json_number(stream, integral.volume);
  stream << ",\"mass\":";
  json_number(stream, integral.mass);
  stream << ",\"energy_q\":";
  json_number(stream, integral.energy);
  stream << ",\"kinetic\":";
  json_number(stream, integral.kinetic);
  stream << ",\"total_energy\":";
  json_number(stream, integral.energy + integral.kinetic);
  stream << "},\"extrema\":{\"p_abs\":[";
  json_number(stream, integral.minimum_pressure);
  stream << ',';
  json_number(stream, integral.maximum_pressure);
  stream << "],\"rho\":[";
  json_number(stream, integral.minimum_density);
  stream << ',';
  json_number(stream, integral.maximum_density);
  stream << "],\"T\":[";
  json_number(stream, integral.minimum_temperature);
  stream << ',';
  json_number(stream, integral.maximum_temperature);
  stream << "],\"h\":[";
  json_number(stream, integral.minimum_enthalpy);
  stream << ',';
  json_number(stream, integral.maximum_enthalpy);
  stream << "],\"q\":[";
  json_number(stream, integral.minimum_energy);
  stream << ',';
  json_number(stream, integral.maximum_energy);
  stream << "],\"maximum_mach\":";
  json_number(stream, integral.maximum_mach);
  stream << "},\"violations\":{\"nonfinite_cells\":"
         << integral.nonfinite_cells << ",\"nonpositive_pressure\":"
         << integral.nonpositive_pressure << ",\"nonpositive_density\":"
         << integral.nonpositive_density
         << ",\"nonpositive_temperature\":"
         << integral.nonpositive_temperature
         << ",\"mach_limit\":" << integral.mach_limit_violations
         << ",\"solid_velocity_nonpositive_zero\":"
         << integral.solid_velocity_nonpositive_zero
         << ",\"nonfinite_faces\":" << faces.nonfinite_faces
         << ",\"ibm_nonpositive_zero\":" << faces.ibm_nonpositive_zero
         << ",\"ibm_negative_zero\":" << faces.ibm_negative_zero
         << "},\"boundaries\":[";
  for (std::size_t face = 0U; face < faces.boundary.size(); ++face) {
    if (face != 0U) stream << ',';
    const BoundaryFaceIntegral& boundary = faces.boundary[face];
    stream << "{\"face\":\"" << kBoundaryNames[face]
           << "\",\"net_outward\":";
    json_number(stream, boundary.net_outward);
    stream << ",\"inflow\":";
    json_number(stream, boundary.inflow);
    stream << ",\"outflow\":";
    json_number(stream, boundary.outflow);
    stream << ",\"maximum_absolute_flux\":";
    json_number(stream, boundary.maximum_absolute_flux);
    stream << ",\"inlet_reversal_faces\":"
           << boundary.inlet_reversal_faces
           << ",\"velocity_inlet_flux_mismatch_faces\":"
           << boundary.velocity_inlet_flux_mismatch_faces
           << ",\"outlet_backflow_faces\":"
           << boundary.outlet_backflow_faces
           << ",\"prohibited_backflow_faces\":"
           << boundary.prohibited_backflow_faces
           << ",\"symmetry_nonpositive_zero\":"
           << boundary.symmetry_nonpositive_zero
           << ",\"symmetry_negative_zero\":"
           << boundary.symmetry_negative_zero << '}';
  }
  stream << "]}";
}

void write_json(std::ostream& stream, const Options& options,
                const EndpointRaw& full, const EndpointRaw& half,
                const ComparisonReport& report) {
  stream << "{\"schema\":\"HUNDUN_V04_FULL_HALF_RESTART_COMPARE_V1\","
            "\"passed\":"
         << (report.failures.empty() ? "true" : "false")
         << ",\"comparison_mode\":\""
         << comparison_mode_name(options.comparison_mode)
         << "\",\"comparison_thresholds\":{\"maximum_relative_rms\":";
  json_number(stream, *options.maximum_relative_rms);
  stream << ",\"maximum_relative_linf\":";
  json_number(stream, *options.maximum_relative_linf);
  stream << "},\"physical_model_fingerprint\":";
  json_fingerprint(stream, report.physical_model_fingerprint);
  stream << ",\"authority\":{\"geometry\":";
  json_fingerprint(stream, report.authority.geometry);
  stream << ",\"thermodynamics\":";
  json_fingerprint(stream, report.authority.thermodynamics);
  stream << ",\"transport\":";
  json_fingerprint(stream, report.authority.transport);
  stream << ",\"source_stl\":";
  json_fingerprint(stream, report.authority.source_stl);
  stream << ",\"surface\":";
  json_fingerprint(stream, report.authority.surface);
  stream << ",\"topology\":";
  json_fingerprint(stream, report.authority.topology);
  stream << ",\"interface_metric\":";
  json_fingerprint(stream, report.authority.interface_metric);
  stream << "},\"time_tolerance\":";
  json_number(stream, report.time_tolerance);
  stream << ",\"maximum_allowed_mach\":";
  json_number(stream, report.maximum_allowed_mach);
  stream << ",\"metrics\":{";
  for (std::size_t index = 0U; index < report.metrics.size(); ++index) {
    if (index != 0U) stream << ',';
    const ErrorMetric& metric = report.metrics[index];
    stream << '"' << kMetricNames[index] << "\":{\"weight\":";
    json_number(stream, metric.weight);
    stream << ",\"rms\":";
    json_number(stream, metric.rms);
    stream << ",\"linf\":";
    json_number(stream, metric.difference_linf);
    stream << ",\"relative_rms\":";
    json_number(stream, metric.relative_rms);
    stream << ",\"relative_linf\":";
    json_number(stream, metric.relative_linf);
    stream << '}';
  }
  stream << "},\"global_relative_difference\":{\"mass\":";
  json_number(stream,
              scalar_relative_difference(report.full.mass, report.half.mass));
  stream << ",\"energy_q\":";
  json_number(stream, scalar_relative_difference(report.full.energy,
                                                 report.half.energy));
  stream << ",\"kinetic\":";
  json_number(stream, scalar_relative_difference(report.full.kinetic,
                                                 report.half.kinetic));
  stream << ",\"total_energy\":";
  json_number(stream,
              scalar_relative_difference(
                  report.full.energy + report.full.kinetic,
                  report.half.energy + report.half.kinetic));
  stream << "},\"full\":";
  write_endpoint_integral(stream, full, report.full, report.full_faces);
  stream << ",\"half\":";
  write_endpoint_integral(stream, half, report.half, report.half_faces);
  stream << ",\"failures\":[";
  for (std::size_t index = 0U; index < report.failures.size(); ++index) {
    if (index != 0U) stream << ',';
    stream << '"' << json_escape(report.failures[index]) << '"';
  }
  stream << "]}\n";
}

Status publish_report(MPI_Comm communicator, const Options& options,
                      const EndpointRaw& full, const EndpointRaw& half,
                      const ComparisonReport& report) noexcept try {
  int rank = 0;
  if (MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kComparatorIo};
  Status local{};
  if (rank == 0) {
    if (options.json.has_value() && *options.json == fs::path("-")) {
      write_json(std::cout, options, full, half, report);
      if (!std::cout) local = {StatusCode::io_failure, kComparatorIo};
    } else if (options.json.has_value()) {
      std::ofstream output(*options.json, std::ios::binary | std::ios::trunc);
      if (!output) {
        local = {StatusCode::io_failure, kComparatorIo};
      } else {
        write_json(output, options, full, half, report);
        output.flush();
        if (!output) local = {StatusCode::io_failure, kComparatorIo};
      }
    } else {
      write_json(std::cout, options, full, half, report);
      if (!std::cout) local = {StatusCode::io_failure, kComparatorIo};
    }
    std::cerr << (report.failures.empty() ? "PASS" : "FAIL")
              << " full/half Restart comparison: "
              << report.failures.size() << " certificate failure(s)\n";
    for (const std::string& failure : report.failures)
      std::cerr << "  - " << failure << '\n';
  }
  return collective_status(communicator, local);
} catch (...) {
  return collective_status(communicator,
                           {StatusCode::io_failure, kComparatorIo});
}

void print_status(int rank, std::string_view operation, Status status) {
  if (rank == 0)
    std::cerr << operation << " failed: " << status_message(status) << " ("
              << static_cast<unsigned>(status.code) << ':' << status.detail
              << ")\n";
}

void print_usage(std::ostream& stream) {
  stream << "usage: v04_full_half_restart_compare FULL_CASE FULL_RESTART "
            "HALF_CASE HALF_RESTART --max-relative-rms VALUE "
            "--max-relative-linf VALUE [--json FILE|-] "
            "[--equivalent-endpoints]\n";
}

int run_self_test(MPI_Comm communicator);
int run_compare(MPI_Comm communicator, const Options& options);

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 70;
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Options options;
  std::string detail;
  const Status parsed = parse_options(argc, argv, options, detail);
  int result = 0;
  if (!parsed) {
    if (rank == 0) {
      print_usage(std::cerr);
      if (!detail.empty()) std::cerr << detail << '\n';
    }
    result = 64;
  } else if (options.self_test) {
    result = run_self_test(MPI_COMM_WORLD);
  } else {
    result = run_compare(MPI_COMM_WORLD, options);
  }
  MPI_Finalize();
  return result;
}

namespace {

bool replace_once(std::string& text, std::string_view from,
                  std::string_view to) {
  const std::size_t position = text.find(from);
  if (position == std::string::npos ||
      text.find(from, position + from.size()) != std::string::npos)
    return false;
  text.replace(position, from.size(), to);
  return true;
}

bool set_fixed_dt(const fs::path& case_root, std::string_view value) {
  const fs::path case_json = case_root / "case.json";
  std::ifstream input(case_json, std::ios::binary);
  if (!input) return false;
  const std::string original{std::istreambuf_iterator<char>(input),
                             std::istreambuf_iterator<char>()};
  if (!input.good() && !input.eof()) return false;
  std::string text = original;
  if (!replace_once(text, "\"initial_dt\":0.001",
                    "\"initial_dt\":" + std::string(value)) ||
      !replace_once(text, "\"minimum_dt\":1e-8",
                    "\"minimum_dt\":" + std::string(value)) ||
      !replace_once(text, "\"maximum_dt\":0.1",
                    "\"maximum_dt\":" + std::string(value)))
    return false;
  std::ofstream output(case_json, std::ios::binary | std::ios::trunc);
  output.write(text.data(), static_cast<std::streamsize>(text.size()));
  output.flush();
  return static_cast<bool>(output);
}

bool set_uniform_open_boundaries(const fs::path& case_root) {
  const fs::path case_json = case_root / "case.json";
  std::ifstream input(case_json, std::ios::binary);
  if (!input) return false;
  const std::string original{std::istreambuf_iterator<char>(input),
                             std::istreambuf_iterator<char>()};
  if (!input.good() && !input.eof()) return false;
  std::string text = original;
  if (!replace_once(text, "\"pressure_reference\": \"closed_mass\"",
                    "\"pressure_reference\": \"boundary_absolute\"") ||
      !replace_once(
          text,
          "\"x_min\": {\"flow_kind\":\"periodic\",\"thermal_kind\":\"none\","
          "\"velocity\":[0,0,0]",
          "\"x_min\": {\"flow_kind\":\"velocity_inlet\","
          "\"thermal_kind\":\"none\",\"velocity\":[0.1,0,0]") ||
      !replace_once(
          text,
          "\"x_max\": {\"flow_kind\":\"periodic\",\"thermal_kind\":\"none\"",
          "\"x_max\": {\"flow_kind\":\"pressure_outlet\","
          "\"thermal_kind\":\"none\""))
    return false;
  std::ofstream output(case_json, std::ios::binary | std::ios::trunc);
  output.write(text.data(), static_cast<std::streamsize>(text.size()));
  output.flush();
  return static_cast<bool>(output);
}

int run_self_test(MPI_Comm communicator) {
  int rank = 0;
  MPI_Comm_rank(communicator, &rank);
  std::array<std::string, 11U> cli_storage{
      "compare", "full-case", "full-restart", "half-case", "half-restart",
      "--json", "report.json", "--max-relative-rms", "1e-7",
      "--max-relative-linf", "2e-6"};
  std::array<char*, 11U> cli{};
  for (std::size_t index = 0U; index < cli.size(); ++index)
    cli[index] = cli_storage[index].data();
  Options parsed_options;
  std::string parsed_detail;
  Status cli_status = parse_options(static_cast<int>(cli.size()), cli.data(),
                                    parsed_options, parsed_detail);
  if (!cli_status || parsed_options.full_case != fs::path("full-case") ||
      parsed_options.full_restart != fs::path("full-restart") ||
      parsed_options.half_case != fs::path("half-case") ||
      parsed_options.half_restart != fs::path("half-restart") ||
      parsed_options.json != std::optional<fs::path>{"report.json"} ||
      !parsed_options.maximum_relative_rms.has_value() ||
      !parsed_options.maximum_relative_linf.has_value() ||
      !same_bits(*parsed_options.maximum_relative_rms, 1.0e-7) ||
      !same_bits(*parsed_options.maximum_relative_linf, 2.0e-6) ||
      parsed_options.comparison_mode != ComparisonMode::full_half)
    cli_status = {StatusCode::invalid_case, kComparatorInput};
  std::array<std::string, 12U> equivalent_cli_storage{
      "compare",          "control-case",      "control-restart",
      "rank-change-case", "rank-change-restart", "--max-relative-rms",
      "1e-7",             "--max-relative-linf", "2e-6",
      "--json",           "equivalent.json",   "--equivalent-endpoints"};
  std::array<char*, 12U> equivalent_cli{};
  for (std::size_t index = 0U; index < equivalent_cli.size(); ++index)
    equivalent_cli[index] = equivalent_cli_storage[index].data();
  Options equivalent_options;
  std::string equivalent_detail;
  const Status equivalent_cli_status = parse_options(
      static_cast<int>(equivalent_cli.size()), equivalent_cli.data(),
      equivalent_options, equivalent_detail);
  if (cli_status &&
      (!equivalent_cli_status ||
       equivalent_options.comparison_mode !=
           ComparisonMode::equivalent_endpoints))
    cli_status = {StatusCode::invalid_case, kComparatorInput};
  const EndpointTiming full_half_coarse{1.0e-3, 1.0e-8, 1.0e-1, 1.0e-3,
                                        1U};
  const EndpointTiming full_half_fine{5.0e-4, 5.0e-9, 5.0e-2, 5.0e-4, 2U};
  const EndpointTiming equivalent_timing{1.0e-3, 1.0e-8, 1.0e-1, 1.0e-3,
                                         7U};
  EndpointTiming unequal_model_dt = equivalent_timing;
  unequal_model_dt.initial_dt =
      std::nextafter(unequal_model_dt.initial_dt, 1.0);
  EndpointTiming unequal_accepted_dt = equivalent_timing;
  unequal_accepted_dt.accepted_dt =
      std::nextafter(unequal_accepted_dt.accepted_dt, 1.0);
  EndpointTiming unequal_step = equivalent_timing;
  ++unequal_step.step;
  if (cli_status &&
      (!endpoint_timing_identity(ComparisonMode::full_half, full_half_coarse,
                                 full_half_fine)
            .valid() ||
       !endpoint_timing_identity(ComparisonMode::equivalent_endpoints,
                                 equivalent_timing, equivalent_timing)
            .valid() ||
       endpoint_timing_identity(ComparisonMode::equivalent_endpoints,
                                equivalent_timing, unequal_model_dt)
           .valid() ||
       endpoint_timing_identity(ComparisonMode::equivalent_endpoints,
                                equivalent_timing, unequal_accepted_dt)
           .valid() ||
       endpoint_timing_identity(ComparisonMode::equivalent_endpoints,
                                equivalent_timing, unequal_step)
           .valid()))
    cli_status = {StatusCode::invalid_case, kComparatorInput};
  const auto rejected_cli = [](std::vector<std::string> storage) {
    std::vector<char*> arguments;
    arguments.reserve(storage.size());
    for (std::string& argument : storage) arguments.push_back(argument.data());
    Options rejected;
    std::string rejected_detail;
    return !parse_options(static_cast<int>(arguments.size()), arguments.data(),
                          rejected, rejected_detail);
  };
  if (cli_status &&
      (!rejected_cli({"compare", "full-case", "full-restart", "half-case",
                      "half-restart"}) ||
       !rejected_cli({"compare", "full-case", "full-restart", "half-case",
                      "half-restart", "--max-relative-rms", "1e-7"}) ||
       !rejected_cli({"compare", "full-case", "full-restart", "half-case",
                      "half-restart", "--max-relative-rms", "inf",
                      "--max-relative-linf", "2e-6"}) ||
       !rejected_cli({"compare", "full-case", "full-restart", "half-case",
                      "half-restart", "--max-relative-rms", "nan",
                      "--max-relative-linf", "2e-6"}) ||
       !rejected_cli({"compare", "full-case", "full-restart", "half-case",
                      "half-restart", "--max-relative-rms", "1e-7junk",
                      "--max-relative-linf", "2e-6"}) ||
       !rejected_cli({"compare", "full-case", "full-restart", "half-case",
                      "half-restart", "--max-relative-rms", "1e-7",
                      "--max-relative-linf", "-1"})))
    cli_status = {StatusCode::invalid_case, kComparatorInput};
  cli_status = collective_status(communicator, cli_status);
  if (!cli_status) return 2;
  std::string root_text;
  Status root_status{};
  if (rank == 0) {
    try {
      root_text =
          (fs::temp_directory_path() /
           ("hundun-v04-full-half-restart-compare-" +
            std::to_string(::getpid())))
              .string();
    } catch (...) {
      root_status = {StatusCode::allocation_failure, kComparatorIo};
    }
  }
  root_status = collective_status(communicator, root_status);
  if (!root_status) return 2;
  std::uint64_t length = root_text.size();
  if (MPI_Bcast(&length, 1, MPI_UINT64_T, 0, communicator) != MPI_SUCCESS)
    return 2;
  Status resize_status{};
  try {
    root_text.resize(static_cast<std::size_t>(length));
  } catch (...) {
    resize_status = {StatusCode::allocation_failure, kComparatorIo};
  }
  resize_status = collective_status(communicator, resize_status);
  if (!resize_status) return 2;
  if (MPI_Bcast(root_text.data(), static_cast<int>(root_text.size()), MPI_CHAR,
                0, communicator) != MPI_SUCCESS)
    return 2;
  const fs::path root = root_text;
  const fs::path full_case = root / "full-case";
  const fs::path half_case = root / "half-case";
  const fs::path full_run = root / "full-run";
  const fs::path half_run = root / "half-run";
  const fs::path report_path = root / "comparison.json";
  const fs::path equivalent_report_path = root / "equivalent-comparison.json";
  Status local{};
  if (rank == 0) {
    std::error_code error;
    fs::remove_all(root, error);
    if (error || !ApplicationService::initialize_case_directory(full_case) ||
        !ApplicationService::initialize_case_directory(half_case) ||
        !set_fixed_dt(full_case, "0.001") ||
        !set_fixed_dt(half_case, "0.0005") ||
        !set_uniform_open_boundaries(full_case) ||
        !set_uniform_open_boundaries(half_case))
      local = {StatusCode::io_failure, kComparatorIo};
  }
  Status status = collective_status(communicator, local);
  if (!status) {
    print_status(rank, "self-test case construction", status);
    return 2;
  }
  if (MPI_Barrier(communicator) != MPI_SUCCESS) return 2;
  ApplicationRunOptions full_options;
  full_options.case_root = full_case;
  full_options.run_directory = full_run;
  full_options.source_root = HUNDUN_V04_SOURCE_ROOT;
  full_options.steps = 1U;
  full_options.output_interval = 0U;
  full_options.restart_interval = 1U;
  ApplicationRunReport full_report;
  status = ApplicationService::run(communicator, full_options, full_report);
  if (!status || full_report.accepted_steps != 1U) {
    if (status)
      status = {StatusCode::numerical_failure, kComparatorNumeric};
    print_status(rank, "self-test full run", status);
    return 2;
  }
  ApplicationRunOptions half_options;
  half_options.case_root = half_case;
  half_options.run_directory = half_run;
  half_options.source_root = HUNDUN_V04_SOURCE_ROOT;
  half_options.steps = 2U;
  half_options.output_interval = 0U;
  half_options.restart_interval = 1U;
  ApplicationRunReport half_report;
  status = ApplicationService::run(communicator, half_options, half_report);
  if (!status || half_report.accepted_steps != 2U) {
    if (status)
      status = {StatusCode::numerical_failure, kComparatorNumeric};
    print_status(rank, "self-test half run", status);
    return 2;
  }
  Options compare;
  compare.full_case = full_case;
  compare.full_restart = full_run / "Restart";
  compare.half_case = half_case;
  compare.half_restart = half_run / "Restart";
  compare.json = report_path;
  compare.maximum_relative_rms = 1.0e-12;
  compare.maximum_relative_linf = 1.0e-12;
  const int compared = run_compare(communicator, compare);
  Options equivalent_compare = compare;
  equivalent_compare.half_case = full_case;
  equivalent_compare.half_restart = full_run / "Restart";
  equivalent_compare.json = equivalent_report_path;
  equivalent_compare.comparison_mode = ComparisonMode::equivalent_endpoints;
  const int equivalent_compared =
      run_compare(communicator, equivalent_compare);
  bool report_valid = compared == 0 && equivalent_compared == 0;
  if (report_valid) {
    EndpointRaw endpoint;
    Status closure_status =
        load_endpoint(communicator, full_case, full_run / "Restart", endpoint);
    PhysicalAuthority closure_authority;
    AuthoritySignature closure_signature;
    if (closure_status)
      closure_status = compile_authority(communicator, endpoint,
                                         closure_authority, closure_signature);
    bool local_closure_valid = static_cast<bool>(closure_status);
    std::uint64_t local_checked_faces = 0U;
    if (closure_status && endpoint.patch.begin.x == 0) {
      std::vector<double> composition;
      try {
        closure_status = configured_velocity_inlet_composition(
            endpoint, endpoint.model.boundaries[0U],
            closure_authority.thermodynamics, composition);
      } catch (...) {
        closure_status = {StatusCode::allocation_failure, kComparatorNumeric};
      }
      const Int3 extents{endpoint.patch.cells.x + 1, endpoint.patch.cells.y,
                         endpoint.patch.cells.z};
      for (std::int32_t z = 0; z < endpoint.patch.cells.z && closure_status;
           ++z)
        for (std::int32_t y = 0; y < endpoint.patch.cells.y && closure_status;
             ++y) {
          const Int3 face{0, y, z};
          const Int3 global{endpoint.patch.begin.x,
                            endpoint.patch.begin.y + y,
                            endpoint.patch.begin.z + z};
          const double area = face_area(closure_authority.geometry, 0U, global);
          double expected = 0.0;
          closure_status = expected_velocity_inlet_flux(
              endpoint, closure_authority, 0U, 0U, false, face, area,
              {composition.data(), composition.size()}, expected);
          const double actual = endpoint.final_mass_flux[0U]
              [face_offset(extents, face)];
          if (closure_status &&
              (!velocity_inlet_flux_matches(endpoint, actual, expected) ||
               velocity_inlet_flux_matches(endpoint, 0.0, expected)))
            local_closure_valid = false;
          ++local_checked_faces;
        }
    }
    if (closure_status && endpoint.patch.begin.x == 0 &&
        endpoint.patch.begin.y == 0 && endpoint.patch.begin.z == 0) {
      const Int3 x_extents{endpoint.patch.cells.x + 1,
                           endpoint.patch.cells.y, endpoint.patch.cells.z};
      endpoint.final_mass_flux[0U][face_offset(x_extents, {0, 0, 0})] = 0.0;
    }
    ErrorMetric corrupted_flux;
    ErrorMetric corrupted_flux_density;
    EndpointFaceIntegral corrupted_full_faces;
    EndpointFaceIntegral corrupted_half_faces;
    if (closure_status)
      closure_status = compare_faces(
          communicator, endpoint, endpoint, closure_authority, corrupted_flux,
          corrupted_flux_density, corrupted_full_faces, corrupted_half_faces);
    if (closure_status &&
        (corrupted_full_faces.boundary[0U]
                 .velocity_inlet_flux_mismatch_faces != 1U ||
         corrupted_half_faces.boundary[0U]
                 .velocity_inlet_flux_mismatch_faces != 1U))
      local_closure_valid = false;
    local_closure_valid = local_closure_valid && closure_status;
    int local_closure = local_closure_valid ? 1 : 0;
    int global_closure = 0;
    std::uint64_t global_checked_faces = 0U;
    if (MPI_Allreduce(&local_closure, &global_closure, 1, MPI_INT, MPI_MIN,
                      communicator) != MPI_SUCCESS ||
        MPI_Allreduce(&local_checked_faces, &global_checked_faces, 1,
                      MPI_UINT64_T, MPI_SUM, communicator) != MPI_SUCCESS)
      report_valid = false;
    else
      report_valid = global_closure != 0 && global_checked_faces != 0U;
  }
  if (rank == 0 && report_valid) {
    std::ostringstream threshold_json;
    threshold_json << "\"comparison_thresholds\":{\"maximum_relative_rms\":";
    json_number(threshold_json, *compare.maximum_relative_rms);
    threshold_json << ",\"maximum_relative_linf\":";
    json_number(threshold_json, *compare.maximum_relative_linf);
    threshold_json << '}';
    const auto valid_report = [&](const fs::path& path,
                                  std::string_view mode) {
      std::ifstream input(path, std::ios::binary);
      const std::string json{std::istreambuf_iterator<char>(input),
                             std::istreambuf_iterator<char>()};
      return static_cast<bool>(input.good() || input.eof()) &&
             json.find("HUNDUN_V04_FULL_HALF_RESTART_COMPARE_V1") !=
                 std::string::npos &&
             json.find("\"passed\":true") != std::string::npos &&
             json.find("\"comparison_mode\":\"" + std::string(mode) +
                       "\"") != std::string::npos &&
             json.find(threshold_json.str()) != std::string::npos &&
             json.find("\"phi_over_A\"") != std::string::npos;
    };
    report_valid = valid_report(report_path, "full_half_2_to_1") &&
                   valid_report(equivalent_report_path,
                                "equivalent_endpoints");
  }
  int local_valid = report_valid ? 1 : 0;
  int global_valid = 0;
  if (MPI_Allreduce(&local_valid, &global_valid, 1, MPI_INT, MPI_MIN,
                    communicator) != MPI_SUCCESS)
    return 2;
  if (global_valid != 0 && rank == 0) {
    std::error_code error;
    fs::remove_all(root, error);
    if (error) global_valid = 0;
  }
  if (MPI_Bcast(&global_valid, 1, MPI_INT, 0, communicator) != MPI_SUCCESS)
    return 2;
  if (global_valid == 0 && rank == 0)
    std::cerr << "full/half comparator self-test retained at " << root << '\n';
  return global_valid != 0 ? 0 : 1;
}

int run_compare(MPI_Comm communicator, const Options& options) {
  int rank = 0;
  MPI_Comm_rank(communicator, &rank);
  EndpointRaw full;
  Status status = load_endpoint(communicator, options.full_case,
                                options.full_restart, full);
  if (!status) {
    print_status(rank, "full endpoint restore", status);
    return 2;
  }
  AuthoritySignature full_signature;
  {
    // Destroy the first full geometry/IBM authority before opening the half
    // endpoint.  Only the compact Restart fields survive this scope.
    PhysicalAuthority full_authority;
    status = compile_authority(communicator, full, full_authority,
                               full_signature);
    if (!status) {
      print_status(rank, "full physical authority", status);
      return 2;
    }
  }
  EndpointRaw half;
  status = load_endpoint(communicator, options.half_case,
                         options.half_restart, half);
  if (!status) {
    print_status(rank, "half endpoint restore", status);
    return 2;
  }
  PhysicalAuthority authority;
  AuthoritySignature half_signature;
  status = compile_authority(communicator, half, authority, half_signature);
  if (!status) {
    print_status(rank, "half physical authority", status);
    return 2;
  }

  ComparisonReport report;
  Status identity_status{};
  try {
    const std::uint64_t full_physical =
        physical_model_fingerprint(full.model);
    const std::uint64_t half_physical =
        physical_model_fingerprint(half.model);
    report.physical_model_fingerprint = half_physical;
    report.authority = half_signature;
    report.maximum_allowed_mach = maximum_allowed_mach(half.model);
    const double step_scale = static_cast<double>(
        std::max<std::uint64_t>(full.step, half.step));
    report.time_tolerance =
        128.0 * std::numeric_limits<double>::epsilon() *
        std::max(1.0, step_scale) *
        std::max({1.0, std::abs(full.time), std::abs(half.time)});
    append_failure(report, full_physical != half_physical,
                   "validated physical models differ outside the three dt "
                   "fields");
    append_failure(report, !same_summary(full.summary, half.summary),
                   "compiled Product plan summaries differ");
    append_failure(report, full.cpu_fingerprint != half.cpu_fingerprint,
                   "CPU execution fingerprints differ");
    append_failure(report,
                   full.schema_fingerprint != half.schema_fingerprint,
                   "Restart field-schema fingerprints differ");
    append_failure(report,
                   full.geometry_fingerprint != half.geometry_fingerprint,
                   "Restart geometry fingerprints differ");
    append_failure(report, full.stl_fingerprint != half.stl_fingerprint,
                   "Product STL fingerprints differ");
    append_failure(report, full.stl_fingerprint != full_signature.source_stl ||
                               half.stl_fingerprint !=
                                   half_signature.source_stl,
                   "Product STL fingerprint does not replay through the "
                   "public IBM compiler");
    append_failure(report, !same(full_signature, half_signature),
                   "geometry/thermo/transport/IBM authority fingerprints "
                   "differ");
    append_failure(report,
                   !std::isfinite(full.time) || !std::isfinite(half.time) ||
                       std::abs(full.time - half.time) > report.time_tolerance,
                   "Restart endpoints are not at the same physical time");
    const EndpointTiming full_timing{
        full.model.time.initial_dt, full.model.time.minimum_dt,
        full.model.time.maximum_dt, full.dt, full.step};
    const EndpointTiming half_timing{
        half.model.time.initial_dt, half.model.time.minimum_dt,
        half.model.time.maximum_dt, half.dt, half.step};
    const TimingIdentity timing = endpoint_timing_identity(
        options.comparison_mode, full_timing, half_timing);
    if (options.comparison_mode == ComparisonMode::equivalent_endpoints) {
      append_failure(report, !timing.model_dt,
                     "equivalent endpoint case initial/minimum/maximum dt "
                     "bits differ");
      append_failure(report, !timing.accepted_dt,
                     "equivalent endpoint Restart accepted dt bits differ");
      append_failure(report, !timing.step,
                     "equivalent endpoint Restart accepted steps differ");
    } else {
      append_failure(report, !timing.model_dt,
                     "half case initial/minimum/maximum dt are not one half "
                     "of the full case");
      append_failure(report, !timing.accepted_dt,
                     "Restart accepted dt does not have the expected 2:1 "
                     "full/half ratio");
      append_failure(report, !timing.step,
                     "Restart accepted step counts do not have the expected "
                     "1:2 full/half ratio");
    }
  } catch (...) {
    identity_status = {StatusCode::allocation_failure, kComparatorAuthority};
  }
  identity_status = collective_status(communicator, identity_status);
  if (!identity_status) {
    print_status(rank, "physical identity", identity_status);
    return 2;
  }

  std::string_view layout_detail;
  status = validate_endpoint_layout(communicator, full, half, authority,
                                    layout_detail);
  if (!status) {
    if (rank == 0 && !layout_detail.empty())
      std::cerr << "endpoint layout: " << layout_detail << '\n';
    print_status(rank, "endpoint layout", status);
    return 2;
  }
  status = compare_cells(communicator, full, half, authority, report.metrics,
                         report.full, report.half,
                         report.maximum_allowed_mach);
  if (!status) {
    print_status(rank, "cell/EOS comparison", status);
    return 2;
  }
  status = compare_faces(
      communicator, full, half, authority,
      report.metrics[static_cast<std::size_t>(MetricId::mass_flux)],
      report.metrics[static_cast<std::size_t>(MetricId::mass_flux_density)],
      report.full_faces, report.half_faces);
  if (!status) {
    print_status(rank, "face-flux comparison", status);
    return 2;
  }
  Status certificate_status{};
  try {
    certify_report(full, half, options, report);
  } catch (...) {
    certificate_status = {StatusCode::allocation_failure, kComparatorNumeric};
  }
  certificate_status = collective_status(communicator, certificate_status);
  if (!certificate_status) {
    print_status(rank, "certificate construction", certificate_status);
    return 2;
  }
  status = publish_report(communicator, options, full, half, report);
  if (!status) {
    print_status(rank, "report publication", status);
    return 2;
  }
  return report.failures.empty() ? 0 : 1;
}

}  // namespace
