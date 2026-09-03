// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/cfg_resolved_case_loader.hpp"

#include "cfg_case_config_loader_detail.hpp"

#include "hundun/rt_error.hpp"
#include "yyjson.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <iterator>
#include <limits>
#include <memory>
#include <set>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace hundun::config {
namespace {

using runtime::ConfigError;
using runtime::Error;

std::string escape_pointer_token(std::string_view token) {
  std::string escaped;
  escaped.reserve(token.size());
  for (const char character : token) {
    if (character == '~') {
      escaped += "~0";
    } else if (character == '/') {
      escaped += "~1";
    } else {
      escaped += character;
    }
  }
  return escaped;
}

std::string child_pointer(std::string_view parent, std::string_view key) {
  std::string result(parent);
  result += '/';
  result += escape_pointer_token(key);
  return result;
}

std::string index_pointer(std::string_view parent, std::size_t index) {
  return child_pointer(parent, std::to_string(index));
}

std::string actual_type(yyjson_val* value) {
  return value == nullptr ? "missing" : yyjson_get_type_desc(value);
}

[[noreturn]] void throw_type_error(std::string pointer,
                                   std::string_view expected,
                                   yyjson_val* actual) {
  throw ConfigError(std::move(pointer),
                    "expected " + std::string(expected) + ", got " +
                        actual_type(actual));
}

yyjson_val* require_member(yyjson_val* object, const char* key,
                           std::string_view pointer) {
  yyjson_val* value = yyjson_obj_get(object, key);
  if (value == nullptr) {
    throw ConfigError(child_pointer(pointer, key), "required member is missing");
  }
  return value;
}

yyjson_val* require_object_member(yyjson_val* object, const char* key,
                                  std::string_view pointer) {
  yyjson_val* value = require_member(object, key, pointer);
  if (!yyjson_is_obj(value)) {
    throw_type_error(child_pointer(pointer, key), "object", value);
  }
  return value;
}

yyjson_val* require_array_member(yyjson_val* object, const char* key,
                                 std::string_view pointer) {
  yyjson_val* value = require_member(object, key, pointer);
  if (!yyjson_is_arr(value)) {
    throw_type_error(child_pointer(pointer, key), "array", value);
  }
  return value;
}

void reject_unknown_keys(
    yyjson_val* object, std::initializer_list<std::string_view> allowed,
    std::string_view pointer) {
  std::vector<std::string> seen;
  seen.reserve(yyjson_obj_size(object));
  yyjson_obj_iter iterator = yyjson_obj_iter_with(object);
  while (yyjson_val* key_value = yyjson_obj_iter_next(&iterator)) {
    const std::string key(yyjson_get_str(key_value), yyjson_get_len(key_value));
    if (std::find(allowed.begin(), allowed.end(), std::string_view(key)) ==
        allowed.end()) {
      throw ConfigError(child_pointer(pointer, key), "unknown member");
    }
    if (std::find(seen.begin(), seen.end(), key) != seen.end()) {
      throw ConfigError(child_pointer(pointer, key), "duplicate member");
    }
    seen.push_back(key);
  }
}

double read_number(yyjson_val* value, std::string pointer) {
  if (!yyjson_is_num(value)) {
    throw_type_error(std::move(pointer), "number", value);
  }
  const double result = yyjson_get_num(value);
  if (!std::isfinite(result)) {
    throw ConfigError(std::move(pointer), "expected a finite number");
  }
  return result;
}

int read_integer(yyjson_val* value, std::string pointer) {
  if (!yyjson_is_int(value)) {
    throw_type_error(std::move(pointer), "integer", value);
  }
  if (yyjson_is_uint(value)) {
    const std::uint64_t result = yyjson_get_uint(value);
    if (result >
        static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
      throw ConfigError(std::move(pointer), "integer exceeds the C++ int range");
    }
    return static_cast<int>(result);
  }
  const std::int64_t result = yyjson_get_sint(value);
  if (result < static_cast<std::int64_t>(std::numeric_limits<int>::min()) ||
      result > static_cast<std::int64_t>(std::numeric_limits<int>::max())) {
    throw ConfigError(std::move(pointer), "integer exceeds the C++ int range");
  }
  return static_cast<int>(result);
}

std::string read_string(yyjson_val* value, std::string pointer) {
  if (!yyjson_is_str(value)) {
    throw_type_error(std::move(pointer), "string", value);
  }
  return std::string(yyjson_get_str(value), yyjson_get_len(value));
}

bool read_boolean(yyjson_val* value, std::string pointer) {
  if (!yyjson_is_bool(value)) {
    throw_type_error(std::move(pointer), "boolean", value);
  }
  return yyjson_get_bool(value);
}

double require_number(yyjson_val* object, const char* key,
                      std::string_view pointer) {
  return read_number(require_member(object, key, pointer),
                     child_pointer(pointer, key));
}

int require_integer(yyjson_val* object, const char* key,
                    std::string_view pointer) {
  return read_integer(require_member(object, key, pointer),
                      child_pointer(pointer, key));
}

std::string require_string(yyjson_val* object, const char* key,
                           std::string_view pointer) {
  return read_string(require_member(object, key, pointer),
                     child_pointer(pointer, key));
}

bool require_boolean(yyjson_val* object, const char* key,
                     std::string_view pointer) {
  return read_boolean(require_member(object, key, pointer),
                      child_pointer(pointer, key));
}

runtime::Int3 require_int3(yyjson_val* object, const char* key,
                          std::string_view pointer) {
  const std::string field_pointer = child_pointer(pointer, key);
  yyjson_val* array = require_member(object, key, pointer);
  if (!yyjson_is_arr(array)) {
    throw_type_error(field_pointer, "array of three integers", array);
  }
  if (yyjson_arr_size(array) != 3U) {
    throw ConfigError(field_pointer, "expected exactly three entries");
  }
  return {read_integer(yyjson_arr_get(array, 0U),
                       index_pointer(field_pointer, 0U)),
          read_integer(yyjson_arr_get(array, 1U),
                       index_pointer(field_pointer, 1U)),
          read_integer(yyjson_arr_get(array, 2U),
                       index_pointer(field_pointer, 2U))};
}

runtime::Real3 require_real3(yyjson_val* object, const char* key,
                            std::string_view pointer) {
  const std::string field_pointer = child_pointer(pointer, key);
  yyjson_val* array = require_member(object, key, pointer);
  if (!yyjson_is_arr(array)) {
    throw_type_error(field_pointer, "array of three numbers", array);
  }
  if (yyjson_arr_size(array) != 3U) {
    throw ConfigError(field_pointer, "expected exactly three entries");
  }
  return {read_number(yyjson_arr_get(array, 0U),
                      index_pointer(field_pointer, 0U)),
          read_number(yyjson_arr_get(array, 1U),
                      index_pointer(field_pointer, 1U)),
          read_number(yyjson_arr_get(array, 2U),
                      index_pointer(field_pointer, 2U))};
}

std::optional<double> optional_number(yyjson_val* object, const char* key,
                                     std::string_view pointer) {
  yyjson_val* value = yyjson_obj_get(object, key);
  if (value == nullptr) {
    return std::nullopt;
  }
  return read_number(value, child_pointer(pointer, key));
}

std::string read_complete_file(const std::filesystem::path& path) {
  std::ifstream stream(path, std::ios::binary);
  if (!stream) {
    throw ConfigError("", "could not open case JSON: " + path.string());
  }
  const std::string contents((std::istreambuf_iterator<char>(stream)),
                             std::istreambuf_iterator<char>());
  if (stream.bad()) {
    throw ConfigError("", "could not read complete case JSON: " + path.string());
  }
  return contents;
}

using Document = std::unique_ptr<yyjson_doc, decltype(&yyjson_doc_free)>;

Document parse_document(std::string contents) {
  yyjson_read_err read_error{};
  Document document(
      yyjson_read_opts(contents.data(), contents.size(), YYJSON_READ_NOFLAG,
                       nullptr, &read_error),
      &yyjson_doc_free);
  if (!document) {
    const std::string detail =
        read_error.msg == nullptr ? "unknown parse error" : read_error.msg;
    throw ConfigError("", "invalid JSON at byte " +
                              std::to_string(read_error.pos) + ": " + detail);
  }
  return document;
}

bool has_parent_component(const std::filesystem::path& path) {
  return std::any_of(path.begin(), path.end(), [](const auto& component) {
    return component == "..";
  });
}

bool has_path_prefix(const std::filesystem::path& path,
                     const std::filesystem::path& prefix) {
  auto path_iterator = path.begin();
  for (auto prefix_iterator = prefix.begin(); prefix_iterator != prefix.end();
       ++prefix_iterator, ++path_iterator) {
    if (path_iterator == path.end() || *path_iterator != *prefix_iterator) {
      return false;
    }
  }
  return true;
}

std::filesystem::path normalize_case_path(
    const std::string& raw, const std::filesystem::path& case_path,
    std::string_view pointer) {
  if (raw.empty()) {
    throw ConfigError(std::string(pointer), "expected a non-empty relative path");
  }
  if (raw.find('\0') != std::string::npos) {
    throw ConfigError(std::string(pointer), "path must not contain a null byte");
  }
  const std::filesystem::path path(raw);
  if (path.is_absolute() || path.has_root_name() || path.has_root_directory()) {
    throw ConfigError(std::string(pointer), "expected an unrooted relative path");
  }
  if (has_parent_component(path)) {
    throw ConfigError(std::string(pointer),
                      "path must not contain a '..' component");
  }
  const std::filesystem::path normalized = path.lexically_normal();
  if (normalized.empty() || normalized == ".") {
    throw ConfigError(std::string(pointer), "expected a non-empty relative path");
  }
  try {
    const std::filesystem::path case_root =
        std::filesystem::absolute(case_path).lexically_normal().parent_path();
    const std::filesystem::path resolved =
        (case_root / normalized).lexically_normal();
    if (!has_path_prefix(resolved, case_root)) {
      throw ConfigError(std::string(pointer), "path escapes the case directory");
    }
  } catch (const std::filesystem::filesystem_error& error) {
    throw ConfigError(std::string(pointer),
                      "could not resolve path against the case directory: " +
                          std::string(error.what()));
  }
  return normalized;
}

void validate_case_path(const std::filesystem::path& path,
                        std::string_view pointer) {
  const std::string text = path.generic_string();
  if (text.empty() || text == ".") {
    throw ConfigError(std::string(pointer), "expected a non-empty relative path");
  }
  if (text.find('\0') != std::string::npos) {
    throw ConfigError(std::string(pointer), "path must not contain a null byte");
  }
  if (path.is_absolute() || path.has_root_name() || path.has_root_directory()) {
    throw ConfigError(std::string(pointer), "expected an unrooted relative path");
  }
  if (has_parent_component(path)) {
    throw ConfigError(std::string(pointer),
                      "path must not contain a '..' component");
  }
  if (text != path.lexically_normal().generic_string()) {
    throw ConfigError(std::string(pointer), "path must be lexically normalized");
  }
}

template <class Enum>
bool valid_enum(Enum value, int last) {
  const int encoded = static_cast<int>(value);
  return encoded >= 0 && encoded <= last;
}

DensityModel parse_density_model(const std::string& value,
                                 const std::string& pointer) {
  if (value == "constant") {
    return DensityModel::constant;
  }
  if (value == "material") {
    return DensityModel::material;
  }
  if (value == "ideal_gas") {
    return DensityModel::ideal_gas;
  }
  throw ConfigError(pointer, "unknown density model '" + value + "'");
}

MeshMapping parse_mesh_mapping(const std::string& value,
                               const std::string& pointer) {
  if (value == "uniform_box") {
    return MeshMapping::uniform_box;
  }
  if (value == "analytic_warped_box") {
    return MeshMapping::analytic_warped_box;
  }
  throw ConfigError(pointer, "unknown mesh mapping '" + value + "'");
}

TimeMode parse_time_mode(const std::string& value,
                         const std::string& pointer) {
  if (value == "fixed") {
    return TimeMode::fixed;
  }
  if (value == "adaptive") {
    return TimeMode::adaptive;
  }
  throw ConfigError(pointer, "unknown time mode '" + value + "'");
}

PatchName parse_patch_name(const std::string& value,
                           const std::string& pointer) {
  static constexpr std::array<std::string_view, 6> names = {
      "x_min", "x_max", "y_min", "y_max", "z_min", "z_max"};
  const auto found = std::find(names.begin(), names.end(), value);
  if (found == names.end()) {
    throw ConfigError(pointer, "unknown boundary patch '" + value + "'");
  }
  return static_cast<PatchName>(std::distance(names.begin(), found));
}

BoundaryType parse_boundary_type(const std::string& value,
                                 const std::string& pointer) {
  if (value == "periodic") {
    return BoundaryType::periodic;
  }
  if (value == "no_slip_wall") {
    return BoundaryType::no_slip_wall;
  }
  if (value == "symmetry") {
    return BoundaryType::symmetry;
  }
  if (value == "velocity_inlet") {
    return BoundaryType::velocity_inlet;
  }
  if (value == "pressure_outlet") {
    return BoundaryType::pressure_outlet;
  }
  throw ConfigError(pointer, "unknown boundary type '" + value + "'");
}

InletThermalAuthority parse_thermal_authority(
    const std::string& value, const std::string& pointer) {
  if (value == "temperature") {
    return InletThermalAuthority::temperature;
  }
  if (value == "enthalpy") {
    return InletThermalAuthority::enthalpy;
  }
  throw ConfigError(pointer, "unknown inlet thermal authority '" + value +
                                 "'");
}

std::uint64_t checked_product(runtime::Int3 values,
                              std::string_view pointer) {
  const std::array<int, 3> components{values.x, values.y, values.z};
  std::uint64_t product = 1U;
  for (const int component : components) {
    if (component <= 0) {
      throw ConfigError(std::string(pointer),
                        "process-grid entries must be positive");
    }
    const auto factor = static_cast<std::uint64_t>(component);
    if (product > std::numeric_limits<std::uint64_t>::max() / factor) {
      throw ConfigError(std::string(pointer), "process-grid product overflows");
    }
    product *= factor;
  }
  return product;
}

void require_positive(double value, std::string pointer) {
  if (!std::isfinite(value) || value <= 0.0) {
    throw ConfigError(std::move(pointer),
                      "expected a finite value greater than zero");
  }
}

void require_nonnegative(double value, std::string pointer) {
  if (!std::isfinite(value) || value < 0.0) {
    throw ConfigError(std::move(pointer),
                      "expected a finite non-negative value");
  }
}

void require_finite(double value, std::string pointer) {
  if (!std::isfinite(value)) {
    throw ConfigError(std::move(pointer), "expected a finite number");
  }
}

void require_derived_positive(double value, std::string pointer,
                              std::string_view quantity) {
  if (!std::isfinite(value) || value <= 0.0) {
    throw ConfigError(
        std::move(pointer),
        "derived " + std::string(quantity) +
            " is not a positive finite binary64 value");
  }
}

double scaled_positive_ratio(double numerator, double denominator_a,
                             double denominator_b) {
  int numerator_exponent = 0;
  int denominator_a_exponent = 0;
  int denominator_b_exponent = 0;
  const double numerator_fraction =
      std::frexp(numerator, &numerator_exponent);
  const double denominator_a_fraction =
      std::frexp(denominator_a, &denominator_a_exponent);
  const double denominator_b_fraction =
      std::frexp(denominator_b, &denominator_b_exponent);

  const double fraction =
      (numerator_fraction / denominator_a_fraction) /
      denominator_b_fraction;
  const std::int64_t exponent =
      static_cast<std::int64_t>(numerator_exponent) -
      static_cast<std::int64_t>(denominator_a_exponent) -
      static_cast<std::int64_t>(denominator_b_exponent);
  if (exponent > std::numeric_limits<int>::max()) {
    return std::numeric_limits<double>::infinity();
  }
  if (exponent < std::numeric_limits<int>::min()) {
    return 0.0;
  }
  return std::scalbn(fraction, static_cast<int>(exponent));
}

bool valid_scalar_name(std::string_view name) {
  if (name.empty()) {
    return false;
  }
  const auto is_alpha_or_underscore = [](char character) {
    return (character >= 'A' && character <= 'Z') ||
           (character >= 'a' && character <= 'z') || character == '_';
  };
  const auto is_alnum_or_underscore = [&](char character) {
    return is_alpha_or_underscore(character) ||
           (character >= '0' && character <= '9');
  };
  return is_alpha_or_underscore(name.front()) &&
         std::all_of(name.begin() + 1, name.end(), is_alnum_or_underscore);
}

bool relatively_equal(double expected, double actual, double tolerance) {
  const double denominator =
      std::max({std::abs(expected), std::abs(actual),
                std::numeric_limits<double>::min()});
  return std::abs(expected - actual) / denominator <= tolerance;
}

void reject_boundary_member(bool present, std::size_t index,
                            std::string_view member) {
  if (present) {
    throw ConfigError(child_pointer(index_pointer("/boundaries", index), member),
                      "member is forbidden for this boundary type");
  }
}

FlowCaseConfig canonicalize_and_validate(FlowCaseConfig config) {
  if (config.schema_version != 2) {
    throw ConfigError("/schema_version", "expected schema version 2");
  }
  if (config.case_name.empty()) {
    throw ConfigError("/case/name", "case name must not be empty");
  }
  if (config.simulation_type != SimulationType::variable_density_flow) {
    throw ConfigError("/simulation/type",
                      "expected 'variable_density_flow'");
  }
  if (!valid_enum(config.density_model,
                  static_cast<int>(DensityModel::ideal_gas))) {
    throw ConfigError("/simulation/density_model", "invalid density model");
  }

  if (config.resources.expected_ranks.has_value() &&
      *config.resources.expected_ranks < 1) {
    throw ConfigError("/resources/expected_ranks",
                      "expected an integer of at least one");
  }
  const std::array<int, 3> cells{config.mesh.cells.x, config.mesh.cells.y,
                                 config.mesh.cells.z};
  for (std::size_t index = 0; index < cells.size(); ++index) {
    if (cells[index] < 1) {
      throw ConfigError(index_pointer("/mesh/cells", index),
                        "cell count must be positive");
    }
  }
  if (config.resources.process_grid.has_value()) {
    const auto grid = *config.resources.process_grid;
    const std::array<int, 3> components{grid.x, grid.y, grid.z};
    for (std::size_t index = 0; index < components.size(); ++index) {
      if (components[index] <= 0) {
        throw ConfigError(index_pointer("/resources/process_grid", index),
                          "process-grid entry must be positive");
      }
      if (components[index] > cells[index]) {
        throw ConfigError(index_pointer("/resources/process_grid", index),
                          "process-grid entry exceeds the cell count");
      }
    }
    const std::uint64_t product =
        checked_product(grid, "/resources/process_grid");
    if (product >
        static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
      throw ConfigError("/resources/process_grid",
                        "process-grid product exceeds the MPI int domain");
    }
    if (config.resources.expected_ranks.has_value() &&
        static_cast<int>(product) != *config.resources.expected_ranks) {
      throw ConfigError("/resources/process_grid",
                        "process-grid product does not equal expected_ranks");
    }
  }

  const std::array<double, 3> origin{config.mesh.origin_m.x,
                                     config.mesh.origin_m.y,
                                     config.mesh.origin_m.z};
  const std::array<double, 3> lengths{config.mesh.length_m.x,
                                      config.mesh.length_m.y,
                                      config.mesh.length_m.z};
  for (std::size_t index = 0; index < origin.size(); ++index) {
    require_finite(origin[index], index_pointer("/mesh/origin_m", index));
    require_positive(lengths[index], index_pointer("/mesh/length_m", index));
  }
  if (!valid_enum(config.mesh.mapping,
                  static_cast<int>(MeshMapping::analytic_warped_box))) {
    throw ConfigError("/mesh/mapping", "invalid mesh mapping");
  }
  if (config.mesh.mapping == MeshMapping::analytic_warped_box) {
    if (!config.mesh.warp_amplitude.has_value()) {
      throw ConfigError("/mesh/warp_amplitude",
                        "member is required for analytic_warped_box");
    }
    const auto amplitude = *config.mesh.warp_amplitude;
    const std::array<double, 3> values{amplitude.x, amplitude.y, amplitude.z};
    for (std::size_t index = 0; index < values.size(); ++index) {
      if (!std::isfinite(values[index]) || std::abs(values[index]) > 0.02) {
        throw ConfigError(index_pointer("/mesh/warp_amplitude", index),
                          "warp amplitude must be finite and at most 0.02");
      }
    }
  } else if (config.mesh.warp_amplitude.has_value()) {
    throw ConfigError("/mesh/warp_amplitude",
                      "member is forbidden for uniform_box");
  }

  if (!valid_enum(config.time.mode, static_cast<int>(TimeMode::adaptive))) {
    throw ConfigError("/time/mode", "invalid time mode");
  }
  if (config.time.steps < 1) {
    throw ConfigError("/time/steps", "steps must be at least one");
  }
  require_positive(config.time.initial_dt_s, "/time/initial_dt_s");
  require_positive(config.time.min_dt_s, "/time/min_dt_s");
  require_positive(config.time.max_dt_s, "/time/max_dt_s");
  if (config.time.min_dt_s > config.time.initial_dt_s) {
    throw ConfigError("/time/min_dt_s", "min_dt_s exceeds initial_dt_s");
  }
  if (config.time.initial_dt_s > config.time.max_dt_s) {
    throw ConfigError("/time/max_dt_s", "initial_dt_s exceeds max_dt_s");
  }
  if (config.time.cfl_target != 0.5) {
    throw ConfigError("/time/cfl_target", "flow requires cfl_target 0.5");
  }
  if (config.time.diffusion_number_target != 0.25) {
    throw ConfigError("/time/diffusion_number_target",
                      "flow requires diffusion_number_target 0.25");
  }
  if (config.time.growth_factor != 1.25) {
    throw ConfigError("/time/growth_factor",
                      "flow requires growth_factor 1.25");
  }
  if (config.time.retry_factor != 0.5) {
    throw ConfigError("/time/retry_factor",
                      "flow requires retry_factor 0.5");
  }
  if (config.time.max_retries != 8) {
    throw ConfigError("/time/max_retries",
                      "flow requires max_retries 8");
  }

  require_positive(config.physics.rho_ref_kg_per_m3,
                   "/physics/rho_ref_kg_per_m3");
  require_positive(config.physics.dynamic_viscosity_pa_s,
                   "/physics/dynamic_viscosity_pa_s");
  require_positive(config.physics.inlet_consistency_rtol,
                   "/physics/inlet_consistency_rtol");
  const bool ideal_gas = config.density_model == DensityModel::ideal_gas;
  const auto require_gas_value = [&](const std::optional<double>& value,
                                     std::string_view pointer) {
    if (!value.has_value()) {
      throw ConfigError(std::string(pointer),
                        "member is required for ideal_gas");
    }
    require_positive(*value, std::string(pointer));
  };
  if (ideal_gas) {
    require_gas_value(config.physics.cp_J_per_kg_K,
                      "/physics/cp_J_per_kg_K");
    require_gas_value(config.physics.gas_constant_J_per_kg_K,
                      "/physics/gas_constant_J_per_kg_K");
    require_gas_value(config.physics.thermodynamic_pressure_pa,
                      "/physics/thermodynamic_pressure_pa");
  } else {
    if (config.physics.cp_J_per_kg_K.has_value()) {
      throw ConfigError("/physics/cp_J_per_kg_K",
                        "member is forbidden for this density model");
    }
    if (config.physics.gas_constant_J_per_kg_K.has_value()) {
      throw ConfigError("/physics/gas_constant_J_per_kg_K",
                        "member is forbidden for this density model");
    }
    if (config.physics.thermodynamic_pressure_pa.has_value()) {
      throw ConfigError("/physics/thermodynamic_pressure_pa",
                        "member is forbidden for this density model");
    }
  }

  static const std::set<std::string> reserved = {
      "rho", "u", "v", "w", "pi", "h", "T", "face_mass_flux"};
  std::set<std::string> scalar_names;
  for (std::size_t index = 0; index < config.scalars.size(); ++index) {
    const auto& scalar = config.scalars[index];
    if (!valid_scalar_name(scalar.name) || reserved.count(scalar.name) != 0U) {
      throw ConfigError(
          child_pointer(index_pointer("/scalars", index), "name"),
          "scalar name is empty, invalid, or reserved");
    }
    if (!scalar_names.insert(scalar.name).second) {
      throw ConfigError(
          child_pointer(index_pointer("/scalars", index), "name"),
          "duplicate scalar name");
    }
    require_nonnegative(
        scalar.diffusivity_m2_per_s,
        child_pointer(index_pointer("/scalars", index),
                      "diffusivity_m2_per_s"));
  }

  std::array<bool, 6> patch_seen{};
  std::array<std::size_t, 6> patch_entry_index{};
  int inlet_count = 0;
  int outlet_count = 0;
  for (std::size_t index = 0; index < config.boundaries.size(); ++index) {
    auto& boundary = config.boundaries[index];
    if (!valid_enum(boundary.patch, static_cast<int>(PatchName::z_max))) {
      throw ConfigError(child_pointer(index_pointer("/boundaries", index),
                                      "patch"),
                        "invalid boundary patch");
    }
    if (!valid_enum(boundary.type,
                    static_cast<int>(BoundaryType::pressure_outlet))) {
      throw ConfigError(child_pointer(index_pointer("/boundaries", index),
                                      "type"),
                        "invalid boundary type");
    }
    const std::size_t patch = static_cast<std::size_t>(boundary.patch);
    if (patch_seen[patch]) {
      throw ConfigError(child_pointer(index_pointer("/boundaries", index),
                                      "patch"),
                        "duplicate boundary patch");
    }
    patch_seen[patch] = true;
    patch_entry_index[patch] = index;

    const bool is_inlet = boundary.type == BoundaryType::velocity_inlet;
    const bool is_outlet = boundary.type == BoundaryType::pressure_outlet;
    if (is_inlet) {
      ++inlet_count;
      if (!boundary.velocity_m_per_s.has_value()) {
        throw ConfigError(child_pointer(index_pointer("/boundaries", index),
                                        "velocity_m_per_s"),
                          "member is required for velocity_inlet");
      }
      const auto velocity = *boundary.velocity_m_per_s;
      const std::array<double, 3> values{velocity.x, velocity.y, velocity.z};
      for (std::size_t component = 0; component < values.size(); ++component) {
        require_finite(values[component],
                       index_pointer(child_pointer(
                                         index_pointer("/boundaries", index),
                                         "velocity_m_per_s"),
                                     component));
      }
      if (!boundary.thermal_authority.has_value() ||
          !valid_enum(*boundary.thermal_authority,
                      static_cast<int>(InletThermalAuthority::enthalpy))) {
        throw ConfigError(child_pointer(index_pointer("/boundaries", index),
                                        "thermal_authority"),
                          "valid thermal authority is required");
      }
      if (!boundary.scalar_values.has_value()) {
        throw ConfigError(child_pointer(index_pointer("/boundaries", index),
                                        "scalar_values"),
                          "member is required for velocity_inlet");
      }
      std::set<std::string> inlet_names;
      for (std::size_t value_index = 0;
           value_index < boundary.scalar_values->size(); ++value_index) {
        const auto& scalar_value = (*boundary.scalar_values)[value_index];
        const std::string value_pointer =
            index_pointer(child_pointer(index_pointer("/boundaries", index),
                                        "scalar_values"),
                          value_index);
        if (scalar_names.count(scalar_value.name) == 0U) {
          throw ConfigError(child_pointer(value_pointer, "name"),
                            "inlet scalar is not declared");
        }
        if (!inlet_names.insert(scalar_value.name).second) {
          throw ConfigError(child_pointer(value_pointer, "name"),
                            "duplicate inlet scalar value");
        }
        require_finite(scalar_value.value,
                       child_pointer(value_pointer, "value"));
      }
      if (inlet_names != scalar_names) {
        throw ConfigError(child_pointer(index_pointer("/boundaries", index),
                                        "scalar_values"),
                          "inlet must provide exactly every declared scalar");
      }
      if (boundary.enthalpy_J_per_kg.has_value()) {
        require_finite(*boundary.enthalpy_J_per_kg,
                       child_pointer(index_pointer("/boundaries", index),
                                     "enthalpy_J_per_kg"));
      }
      if (boundary.temperature_K.has_value()) {
        require_positive(*boundary.temperature_K,
                         child_pointer(index_pointer("/boundaries", index),
                                       "temperature_K"));
      }
      if (boundary.density_kg_per_m3.has_value()) {
        require_positive(*boundary.density_kg_per_m3,
                         child_pointer(index_pointer("/boundaries", index),
                                       "density_kg_per_m3"));
      }
      if (ideal_gas) {
        const std::string boundary_pointer =
            index_pointer("/boundaries", index);
        double temperature = 0.0;
        if (*boundary.thermal_authority ==
            InletThermalAuthority::temperature) {
          if (!boundary.temperature_K.has_value()) {
            throw ConfigError(child_pointer(index_pointer("/boundaries", index),
                                            "temperature_K"),
                              "authoritative temperature is required");
          }
          temperature = *boundary.temperature_K;
        } else {
          if (!boundary.enthalpy_J_per_kg.has_value()) {
            throw ConfigError(child_pointer(index_pointer("/boundaries", index),
                                            "enthalpy_J_per_kg"),
                              "authoritative enthalpy is required");
          }
          temperature = *boundary.enthalpy_J_per_kg /
                        *config.physics.cp_J_per_kg_K;
          require_derived_positive(
              temperature,
              child_pointer(boundary_pointer, "enthalpy_J_per_kg"),
              "temperature");
        }
        const double enthalpy = *config.physics.cp_J_per_kg_K * temperature;
        require_derived_positive(
            enthalpy, child_pointer(boundary_pointer, "enthalpy_J_per_kg"),
            "enthalpy");
        const double density = scaled_positive_ratio(
            *config.physics.thermodynamic_pressure_pa,
            *config.physics.gas_constant_J_per_kg_K, temperature);
        require_derived_positive(
            density, child_pointer(boundary_pointer, "density_kg_per_m3"),
            "density");
        const double tolerance = config.physics.inlet_consistency_rtol;
        if (boundary.temperature_K.has_value() &&
            !relatively_equal(temperature, *boundary.temperature_K,
                              tolerance)) {
          throw ConfigError(child_pointer(index_pointer("/boundaries", index),
                                          "temperature_K"),
                            "temperature is inconsistent with authority");
        }
        if (boundary.enthalpy_J_per_kg.has_value() &&
            !relatively_equal(enthalpy, *boundary.enthalpy_J_per_kg,
                              tolerance)) {
          throw ConfigError(child_pointer(index_pointer("/boundaries", index),
                                          "enthalpy_J_per_kg"),
                            "enthalpy is inconsistent with authority");
        }
        if (boundary.density_kg_per_m3.has_value() &&
            !relatively_equal(density, *boundary.density_kg_per_m3,
                              tolerance)) {
          throw ConfigError(child_pointer(index_pointer("/boundaries", index),
                                          "density_kg_per_m3"),
                            "density is inconsistent with ideal gas closure");
        }
      } else {
        if (*boundary.thermal_authority !=
            InletThermalAuthority::enthalpy) {
          throw ConfigError(child_pointer(index_pointer("/boundaries", index),
                                          "thermal_authority"),
                            "constant/material inlet authority must be enthalpy");
        }
        if (!boundary.enthalpy_J_per_kg.has_value()) {
          throw ConfigError(child_pointer(index_pointer("/boundaries", index),
                                          "enthalpy_J_per_kg"),
                            "inlet enthalpy is required");
        }
        if (boundary.temperature_K.has_value()) {
          throw ConfigError(child_pointer(index_pointer("/boundaries", index),
                                          "temperature_K"),
                            "temperature is forbidden for this density model");
        }
        if (config.density_model == DensityModel::material &&
            !boundary.density_kg_per_m3.has_value()) {
          throw ConfigError(child_pointer(index_pointer("/boundaries", index),
                                          "density_kg_per_m3"),
                            "material-density inlet density is required");
        }
        if (config.density_model == DensityModel::constant &&
            boundary.density_kg_per_m3.has_value() &&
            !relatively_equal(config.physics.rho_ref_kg_per_m3,
                              *boundary.density_kg_per_m3,
                              config.physics.inlet_consistency_rtol)) {
          throw ConfigError(child_pointer(index_pointer("/boundaries", index),
                                          "density_kg_per_m3"),
                            "inlet density is inconsistent with rho_ref");
        }
      }
      reject_boundary_member(boundary.pressure_perturbation_pa.has_value(),
                             index, "pressure_perturbation_pa");
    } else if (is_outlet) {
      ++outlet_count;
      if (!boundary.pressure_perturbation_pa.has_value()) {
        throw ConfigError(child_pointer(index_pointer("/boundaries", index),
                                        "pressure_perturbation_pa"),
                          "member is required for pressure_outlet");
      }
      require_finite(*boundary.pressure_perturbation_pa,
                     child_pointer(index_pointer("/boundaries", index),
                                   "pressure_perturbation_pa"));
      reject_boundary_member(boundary.velocity_m_per_s.has_value(), index,
                             "velocity_m_per_s");
      reject_boundary_member(boundary.thermal_authority.has_value(), index,
                             "thermal_authority");
      reject_boundary_member(boundary.temperature_K.has_value(), index,
                             "temperature_K");
      reject_boundary_member(boundary.enthalpy_J_per_kg.has_value(), index,
                             "enthalpy_J_per_kg");
      reject_boundary_member(boundary.density_kg_per_m3.has_value(), index,
                             "density_kg_per_m3");
      reject_boundary_member(boundary.scalar_values.has_value(), index,
                             "scalar_values");
    } else {
      reject_boundary_member(boundary.velocity_m_per_s.has_value(), index,
                             "velocity_m_per_s");
      reject_boundary_member(boundary.thermal_authority.has_value(), index,
                             "thermal_authority");
      reject_boundary_member(boundary.temperature_K.has_value(), index,
                             "temperature_K");
      reject_boundary_member(boundary.enthalpy_J_per_kg.has_value(), index,
                             "enthalpy_J_per_kg");
      reject_boundary_member(boundary.density_kg_per_m3.has_value(), index,
                             "density_kg_per_m3");
      reject_boundary_member(boundary.scalar_values.has_value(), index,
                             "scalar_values");
      reject_boundary_member(boundary.pressure_perturbation_pa.has_value(),
                             index, "pressure_perturbation_pa");
    }
  }
  if (!std::all_of(patch_seen.begin(), patch_seen.end(), [](bool value) {
        return value;
      })) {
    throw ConfigError("/boundaries", "all six patches are required");
  }
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    const std::size_t lower = axis * 2U;
    const std::size_t upper = lower + 1U;
    const bool lower_periodic =
        config.boundaries[patch_entry_index[lower]].type ==
        BoundaryType::periodic;
    const bool upper_periodic =
        config.boundaries[patch_entry_index[upper]].type ==
        BoundaryType::periodic;
    if (lower_periodic != upper_periodic) {
      const std::size_t offender = lower_periodic
                                       ? patch_entry_index[lower]
                                       : patch_entry_index[upper];
      throw ConfigError(child_pointer(index_pointer("/boundaries", offender),
                                      "type"),
                        "periodic boundaries must be paired by axis");
    }
  }
  if (!((inlet_count == 0 && outlet_count == 0) ||
        (inlet_count == 1 && outlet_count == 1))) {
    throw ConfigError("/boundaries",
                      "open cases require exactly one inlet and one outlet");
  }

  if (config.restart.read) {
    if (!config.restart.read_directory.has_value()) {
      throw ConfigError("/restart/read_directory",
                        "member is required when read is true");
    }
    validate_case_path(*config.restart.read_directory,
                       "/restart/read_directory");
  } else if (config.restart.read_directory.has_value()) {
    throw ConfigError("/restart/read_directory",
                      "member is forbidden when read is false");
  }
  validate_case_path(config.restart.write_directory,
                     "/restart/write_directory");
  if (config.restart.write_interval < 1) {
    throw ConfigError("/restart/write_interval",
                      "write interval must be positive");
  }
  validate_case_path(config.diagnostics.directory,
                     "/diagnostics/directory");
  if (config.diagnostics.write_interval < 1) {
    throw ConfigError("/diagnostics/write_interval",
                      "write interval must be positive");
  }
  validate_case_path(config.performance.directory,
                     "/performance/directory");
  if (config.performance.warmup_steps < 0) {
    throw ConfigError("/performance/warmup_steps",
                      "warmup steps must be non-negative");
  }
  if (config.performance.measured_steps < 1) {
    throw ConfigError("/performance/measured_steps",
                      "measured steps must be positive");
  }
  if (config.performance.repetitions < 1) {
    throw ConfigError("/performance/repetitions",
                      "repetitions must be positive");
  }

  std::sort(config.scalars.begin(), config.scalars.end(),
            [](const auto& left, const auto& right) {
              return left.name < right.name;
            });
  for (auto& boundary : config.boundaries) {
    if (boundary.scalar_values.has_value()) {
      std::sort(boundary.scalar_values->begin(), boundary.scalar_values->end(),
                [](const auto& left, const auto& right) {
                  return left.name < right.name;
                });
    }
  }
  std::sort(config.boundaries.begin(), config.boundaries.end(),
            [](const auto& left, const auto& right) {
              return static_cast<int>(left.patch) <
                     static_cast<int>(right.patch);
            });
  return config;
}

FlowCaseConfig parse_flow_case(yyjson_val* root,
                               const std::filesystem::path& case_path) {
  if (!yyjson_is_obj(root)) {
    throw_type_error("", "object", root);
  }
  reject_unknown_keys(
      root,
      {"schema_version", "case", "simulation", "resources", "mesh", "time",
       "physics", "scalars", "boundaries", "restart", "diagnostics",
       "performance"},
      "");

  FlowCaseConfig config{};
  config.schema_version = require_integer(root, "schema_version", "");

  yyjson_val* case_object = require_object_member(root, "case", "");
  reject_unknown_keys(case_object, {"name"}, "/case");
  config.case_name = require_string(case_object, "name", "/case");

  yyjson_val* simulation = require_object_member(root, "simulation", "");
  reject_unknown_keys(simulation, {"type", "density_model"}, "/simulation");
  const std::string simulation_type =
      require_string(simulation, "type", "/simulation");
  if (simulation_type != "variable_density_flow") {
    throw ConfigError("/simulation/type",
                      "expected 'variable_density_flow'");
  }
  config.simulation_type = SimulationType::variable_density_flow;
  config.density_model = parse_density_model(
      require_string(simulation, "density_model", "/simulation"),
      "/simulation/density_model");

  if (yyjson_val* resources = yyjson_obj_get(root, "resources")) {
    if (!yyjson_is_obj(resources)) {
      throw_type_error("/resources", "object", resources);
    }
    reject_unknown_keys(resources, {"expected_ranks", "process_grid"},
                        "/resources");
    if (yyjson_val* expected = yyjson_obj_get(resources, "expected_ranks")) {
      config.resources.expected_ranks =
          read_integer(expected, "/resources/expected_ranks");
    }
    if (yyjson_obj_get(resources, "process_grid") != nullptr) {
      config.resources.process_grid =
          require_int3(resources, "process_grid", "/resources");
    }
  }

  yyjson_val* mesh = require_object_member(root, "mesh", "");
  reject_unknown_keys(mesh,
                      {"cells", "origin_m", "length_m", "mapping",
                       "warp_amplitude"},
                      "/mesh");
  config.mesh.cells = require_int3(mesh, "cells", "/mesh");
  config.mesh.origin_m = require_real3(mesh, "origin_m", "/mesh");
  config.mesh.length_m = require_real3(mesh, "length_m", "/mesh");
  config.mesh.mapping = parse_mesh_mapping(
      require_string(mesh, "mapping", "/mesh"), "/mesh/mapping");
  if (yyjson_obj_get(mesh, "warp_amplitude") != nullptr) {
    config.mesh.warp_amplitude =
        require_real3(mesh, "warp_amplitude", "/mesh");
  }

  yyjson_val* time = require_object_member(root, "time", "");
  reject_unknown_keys(time,
                      {"mode", "steps", "initial_dt_s", "min_dt_s",
                       "max_dt_s", "cfl_target", "diffusion_number_target",
                       "growth_factor", "retry_factor", "max_retries"},
                      "/time");
  config.time.mode = parse_time_mode(require_string(time, "mode", "/time"),
                                     "/time/mode");
  config.time.steps = require_integer(time, "steps", "/time");
  config.time.initial_dt_s = require_number(time, "initial_dt_s", "/time");
  config.time.min_dt_s = require_number(time, "min_dt_s", "/time");
  config.time.max_dt_s = require_number(time, "max_dt_s", "/time");
  config.time.cfl_target = require_number(time, "cfl_target", "/time");
  config.time.diffusion_number_target =
      require_number(time, "diffusion_number_target", "/time");
  config.time.growth_factor =
      require_number(time, "growth_factor", "/time");
  config.time.retry_factor = require_number(time, "retry_factor", "/time");
  config.time.max_retries = require_integer(time, "max_retries", "/time");

  yyjson_val* physics = require_object_member(root, "physics", "");
  reject_unknown_keys(
      physics,
      {"rho_ref_kg_per_m3", "dynamic_viscosity_pa_s",
       "inlet_consistency_rtol", "cp_J_per_kg_K",
       "gas_constant_J_per_kg_K", "thermodynamic_pressure_pa"},
      "/physics");
  config.physics.rho_ref_kg_per_m3 =
      require_number(physics, "rho_ref_kg_per_m3", "/physics");
  config.physics.dynamic_viscosity_pa_s =
      require_number(physics, "dynamic_viscosity_pa_s", "/physics");
  config.physics.inlet_consistency_rtol =
      require_number(physics, "inlet_consistency_rtol", "/physics");
  config.physics.cp_J_per_kg_K =
      optional_number(physics, "cp_J_per_kg_K", "/physics");
  config.physics.gas_constant_J_per_kg_K =
      optional_number(physics, "gas_constant_J_per_kg_K", "/physics");
  config.physics.thermodynamic_pressure_pa =
      optional_number(physics, "thermodynamic_pressure_pa", "/physics");

  yyjson_val* scalars = require_array_member(root, "scalars", "");
  config.scalars.reserve(yyjson_arr_size(scalars));
  for (std::size_t index = 0; index < yyjson_arr_size(scalars); ++index) {
    yyjson_val* scalar = yyjson_arr_get(scalars, index);
    const std::string pointer = index_pointer("/scalars", index);
    if (!yyjson_is_obj(scalar)) {
      throw_type_error(pointer, "object", scalar);
    }
    reject_unknown_keys(scalar, {"name", "diffusivity_m2_per_s"}, pointer);
    config.scalars.push_back(
        {require_string(scalar, "name", pointer),
         require_number(scalar, "diffusivity_m2_per_s", pointer)});
  }

  yyjson_val* boundaries = require_array_member(root, "boundaries", "");
  if (yyjson_arr_size(boundaries) != config.boundaries.size()) {
    throw ConfigError("/boundaries", "expected exactly six boundary entries");
  }
  for (std::size_t index = 0; index < config.boundaries.size(); ++index) {
    yyjson_val* boundary_object = yyjson_arr_get(boundaries, index);
    const std::string pointer = index_pointer("/boundaries", index);
    if (!yyjson_is_obj(boundary_object)) {
      throw_type_error(pointer, "object", boundary_object);
    }
    reject_unknown_keys(
        boundary_object,
        {"patch", "type", "velocity_m_per_s", "thermal_authority",
         "temperature_K", "enthalpy_J_per_kg", "density_kg_per_m3",
         "scalar_values", "pressure_perturbation_pa"},
        pointer);
    auto& boundary = config.boundaries[index];
    boundary.patch = parse_patch_name(
        require_string(boundary_object, "patch", pointer),
        child_pointer(pointer, "patch"));
    boundary.type = parse_boundary_type(
        require_string(boundary_object, "type", pointer),
        child_pointer(pointer, "type"));
    if (yyjson_obj_get(boundary_object, "velocity_m_per_s") != nullptr) {
      boundary.velocity_m_per_s =
          require_real3(boundary_object, "velocity_m_per_s", pointer);
    }
    if (yyjson_obj_get(boundary_object, "thermal_authority") != nullptr) {
      boundary.thermal_authority = parse_thermal_authority(
          require_string(boundary_object, "thermal_authority", pointer),
          child_pointer(pointer, "thermal_authority"));
    }
    boundary.temperature_K =
        optional_number(boundary_object, "temperature_K", pointer);
    boundary.enthalpy_J_per_kg =
        optional_number(boundary_object, "enthalpy_J_per_kg", pointer);
    boundary.density_kg_per_m3 =
        optional_number(boundary_object, "density_kg_per_m3", pointer);
    boundary.pressure_perturbation_pa = optional_number(
        boundary_object, "pressure_perturbation_pa", pointer);
    if (yyjson_val* scalar_values =
            yyjson_obj_get(boundary_object, "scalar_values")) {
      if (!yyjson_is_arr(scalar_values)) {
        throw_type_error(child_pointer(pointer, "scalar_values"), "array",
                         scalar_values);
      }
      boundary.scalar_values.emplace();
      boundary.scalar_values->reserve(yyjson_arr_size(scalar_values));
      for (std::size_t value_index = 0;
           value_index < yyjson_arr_size(scalar_values); ++value_index) {
        yyjson_val* scalar_value = yyjson_arr_get(scalar_values, value_index);
        const std::string value_pointer =
            index_pointer(child_pointer(pointer, "scalar_values"), value_index);
        if (!yyjson_is_obj(scalar_value)) {
          throw_type_error(value_pointer, "object", scalar_value);
        }
        reject_unknown_keys(scalar_value, {"name", "value"}, value_pointer);
        boundary.scalar_values->push_back(
            {require_string(scalar_value, "name", value_pointer),
             require_number(scalar_value, "value", value_pointer)});
      }
    }
  }

  yyjson_val* restart = require_object_member(root, "restart", "");
  reject_unknown_keys(restart,
                      {"read", "read_directory", "write_directory",
                       "write_interval"},
                      "/restart");
  config.restart.read = require_boolean(restart, "read", "/restart");
  if (yyjson_val* read_directory = yyjson_obj_get(restart, "read_directory")) {
    config.restart.read_directory = normalize_case_path(
        read_string(read_directory, "/restart/read_directory"), case_path,
        "/restart/read_directory");
  }
  config.restart.write_directory = normalize_case_path(
      require_string(restart, "write_directory", "/restart"), case_path,
      "/restart/write_directory");
  config.restart.write_interval =
      require_integer(restart, "write_interval", "/restart");

  yyjson_val* diagnostics = require_object_member(root, "diagnostics", "");
  reject_unknown_keys(diagnostics,
                      {"directory", "write_interval", "write_mesh"},
                      "/diagnostics");
  config.diagnostics.directory = normalize_case_path(
      require_string(diagnostics, "directory", "/diagnostics"), case_path,
      "/diagnostics/directory");
  config.diagnostics.write_interval =
      require_integer(diagnostics, "write_interval", "/diagnostics");
  config.diagnostics.write_mesh =
      require_boolean(diagnostics, "write_mesh", "/diagnostics");

  yyjson_val* performance = require_object_member(root, "performance", "");
  reject_unknown_keys(performance,
                      {"enabled", "directory", "warmup_steps",
                       "measured_steps", "repetitions"},
                      "/performance");
  config.performance.enabled =
      require_boolean(performance, "enabled", "/performance");
  config.performance.directory = normalize_case_path(
      require_string(performance, "directory", "/performance"), case_path,
      "/performance/directory");
  config.performance.warmup_steps =
      require_integer(performance, "warmup_steps", "/performance");
  config.performance.measured_steps =
      require_integer(performance, "measured_steps", "/performance");
  config.performance.repetitions =
      require_integer(performance, "repetitions", "/performance");

  return canonicalize_and_validate(std::move(config));
}

const char* density_model_name(DensityModel value) {
  switch (value) {
    case DensityModel::constant:
      return "constant";
    case DensityModel::material:
      return "material";
    case DensityModel::ideal_gas:
      return "ideal_gas";
  }
  throw ConfigError("/simulation/density_model", "invalid density model");
}

const char* mesh_mapping_name(MeshMapping value) {
  switch (value) {
    case MeshMapping::uniform_box:
      return "uniform_box";
    case MeshMapping::analytic_warped_box:
      return "analytic_warped_box";
  }
  throw ConfigError("/mesh/mapping", "invalid mesh mapping");
}

const char* time_mode_name(TimeMode value) {
  switch (value) {
    case TimeMode::fixed:
      return "fixed";
    case TimeMode::adaptive:
      return "adaptive";
  }
  throw ConfigError("/time/mode", "invalid time mode");
}

const char* patch_name(PatchName value) {
  static constexpr std::array<const char*, 6> names = {
      "x_min", "x_max", "y_min", "y_max", "z_min", "z_max"};
  if (!valid_enum(value, static_cast<int>(PatchName::z_max))) {
    throw ConfigError("/boundaries", "invalid boundary patch");
  }
  return names[static_cast<std::size_t>(value)];
}

const char* boundary_type_name(BoundaryType value) {
  switch (value) {
    case BoundaryType::periodic:
      return "periodic";
    case BoundaryType::no_slip_wall:
      return "no_slip_wall";
    case BoundaryType::symmetry:
      return "symmetry";
    case BoundaryType::velocity_inlet:
      return "velocity_inlet";
    case BoundaryType::pressure_outlet:
      return "pressure_outlet";
  }
  throw ConfigError("/boundaries", "invalid boundary type");
}

const char* thermal_authority_name(InletThermalAuthority value) {
  switch (value) {
    case InletThermalAuthority::temperature:
      return "temperature";
    case InletThermalAuthority::enthalpy:
      return "enthalpy";
  }
  throw ConfigError("/boundaries", "invalid inlet thermal authority");
}

void require_json_write(bool successful) {
  if (!successful) {
    throw Error("failed to allocate deterministic resolved JSON");
  }
}

yyjson_mut_val* add_object(yyjson_mut_doc* document, yyjson_mut_val* parent,
                           const char* key) {
  yyjson_mut_val* object = yyjson_mut_obj_add_obj(document, parent, key);
  if (object == nullptr) {
    throw Error("failed to allocate deterministic resolved JSON");
  }
  return object;
}

yyjson_mut_val* add_array(yyjson_mut_doc* document, yyjson_mut_val* parent,
                          const char* key) {
  yyjson_mut_val* array = yyjson_mut_obj_add_arr(document, parent, key);
  if (array == nullptr) {
    throw Error("failed to allocate deterministic resolved JSON");
  }
  return array;
}

void add_string(yyjson_mut_doc* document, yyjson_mut_val* object,
                const char* key, const std::string& value) {
  require_json_write(yyjson_mut_obj_add_strncpy(
      document, object, key, value.data(), value.size()));
}

void add_int3(yyjson_mut_doc* document, yyjson_mut_val* object,
              const char* key, runtime::Int3 value) {
  yyjson_mut_val* array = add_array(document, object, key);
  require_json_write(yyjson_mut_arr_add_int(document, array, value.x));
  require_json_write(yyjson_mut_arr_add_int(document, array, value.y));
  require_json_write(yyjson_mut_arr_add_int(document, array, value.z));
}

void add_real3(yyjson_mut_doc* document, yyjson_mut_val* object,
               const char* key, runtime::Real3 value) {
  yyjson_mut_val* array = add_array(document, object, key);
  require_json_write(yyjson_mut_arr_add_real(document, array, value.x));
  require_json_write(yyjson_mut_arr_add_real(document, array, value.y));
  require_json_write(yyjson_mut_arr_add_real(document, array, value.z));
}

std::string serialize_flow_case(const FlowCaseConfig& input) {
  const FlowCaseConfig config = canonicalize_and_validate(input);
  using MutableDocument =
      std::unique_ptr<yyjson_mut_doc, decltype(&yyjson_mut_doc_free)>;
  MutableDocument document(yyjson_mut_doc_new(nullptr), &yyjson_mut_doc_free);
  if (!document) {
    throw Error("failed to allocate deterministic resolved JSON");
  }
  yyjson_mut_val* root = yyjson_mut_obj(document.get());
  if (root == nullptr) {
    throw Error("failed to allocate deterministic resolved JSON");
  }
  yyjson_mut_doc_set_root(document.get(), root);

  require_json_write(yyjson_mut_obj_add_int(
      document.get(), root, "schema_version", config.schema_version));
  yyjson_mut_val* case_object = add_object(document.get(), root, "case");
  add_string(document.get(), case_object, "name", config.case_name);
  yyjson_mut_val* simulation =
      add_object(document.get(), root, "simulation");
  require_json_write(yyjson_mut_obj_add_str(
      document.get(), simulation, "type", "variable_density_flow"));
  require_json_write(yyjson_mut_obj_add_str(
      document.get(), simulation, "density_model",
      density_model_name(config.density_model)));

  yyjson_mut_val* resources = add_object(document.get(), root, "resources");
  if (config.resources.expected_ranks.has_value()) {
    require_json_write(yyjson_mut_obj_add_int(
        document.get(), resources, "expected_ranks",
        *config.resources.expected_ranks));
  }
  if (config.resources.process_grid.has_value()) {
    add_int3(document.get(), resources, "process_grid",
             *config.resources.process_grid);
  }

  yyjson_mut_val* mesh = add_object(document.get(), root, "mesh");
  add_int3(document.get(), mesh, "cells", config.mesh.cells);
  add_real3(document.get(), mesh, "origin_m", config.mesh.origin_m);
  add_real3(document.get(), mesh, "length_m", config.mesh.length_m);
  require_json_write(yyjson_mut_obj_add_str(
      document.get(), mesh, "mapping", mesh_mapping_name(config.mesh.mapping)));
  if (config.mesh.warp_amplitude.has_value()) {
    add_real3(document.get(), mesh, "warp_amplitude",
              *config.mesh.warp_amplitude);
  }

  yyjson_mut_val* time = add_object(document.get(), root, "time");
  require_json_write(yyjson_mut_obj_add_str(
      document.get(), time, "mode", time_mode_name(config.time.mode)));
  require_json_write(yyjson_mut_obj_add_int(document.get(), time, "steps",
                                            config.time.steps));
  require_json_write(yyjson_mut_obj_add_real(
      document.get(), time, "initial_dt_s", config.time.initial_dt_s));
  require_json_write(yyjson_mut_obj_add_real(
      document.get(), time, "min_dt_s", config.time.min_dt_s));
  require_json_write(yyjson_mut_obj_add_real(
      document.get(), time, "max_dt_s", config.time.max_dt_s));
  require_json_write(yyjson_mut_obj_add_real(
      document.get(), time, "cfl_target", config.time.cfl_target));
  require_json_write(yyjson_mut_obj_add_real(
      document.get(), time, "diffusion_number_target",
      config.time.diffusion_number_target));
  require_json_write(yyjson_mut_obj_add_real(
      document.get(), time, "growth_factor", config.time.growth_factor));
  require_json_write(yyjson_mut_obj_add_real(
      document.get(), time, "retry_factor", config.time.retry_factor));
  require_json_write(yyjson_mut_obj_add_int(
      document.get(), time, "max_retries", config.time.max_retries));

  yyjson_mut_val* physics = add_object(document.get(), root, "physics");
  require_json_write(yyjson_mut_obj_add_real(
      document.get(), physics, "rho_ref_kg_per_m3",
      config.physics.rho_ref_kg_per_m3));
  require_json_write(yyjson_mut_obj_add_real(
      document.get(), physics, "dynamic_viscosity_pa_s",
      config.physics.dynamic_viscosity_pa_s));
  require_json_write(yyjson_mut_obj_add_real(
      document.get(), physics, "inlet_consistency_rtol",
      config.physics.inlet_consistency_rtol));
  if (config.physics.cp_J_per_kg_K.has_value()) {
    require_json_write(yyjson_mut_obj_add_real(
        document.get(), physics, "cp_J_per_kg_K",
        *config.physics.cp_J_per_kg_K));
    require_json_write(yyjson_mut_obj_add_real(
        document.get(), physics, "gas_constant_J_per_kg_K",
        *config.physics.gas_constant_J_per_kg_K));
    require_json_write(yyjson_mut_obj_add_real(
        document.get(), physics, "thermodynamic_pressure_pa",
        *config.physics.thermodynamic_pressure_pa));
  }

  yyjson_mut_val* scalars = add_array(document.get(), root, "scalars");
  for (const auto& scalar : config.scalars) {
    yyjson_mut_val* object = yyjson_mut_obj(document.get());
    if (object == nullptr ||
        !yyjson_mut_arr_append(scalars, object)) {
      throw Error("failed to allocate deterministic resolved JSON");
    }
    add_string(document.get(), object, "name", scalar.name);
    require_json_write(yyjson_mut_obj_add_real(
        document.get(), object, "diffusivity_m2_per_s",
        scalar.diffusivity_m2_per_s));
  }

  yyjson_mut_val* boundaries = add_array(document.get(), root, "boundaries");
  for (const auto& boundary : config.boundaries) {
    yyjson_mut_val* object = yyjson_mut_obj(document.get());
    if (object == nullptr || !yyjson_mut_arr_append(boundaries, object)) {
      throw Error("failed to allocate deterministic resolved JSON");
    }
    require_json_write(yyjson_mut_obj_add_str(
        document.get(), object, "patch", patch_name(boundary.patch)));
    require_json_write(yyjson_mut_obj_add_str(
        document.get(), object, "type", boundary_type_name(boundary.type)));
    if (boundary.velocity_m_per_s.has_value()) {
      add_real3(document.get(), object, "velocity_m_per_s",
                *boundary.velocity_m_per_s);
    }
    if (boundary.thermal_authority.has_value()) {
      require_json_write(yyjson_mut_obj_add_str(
          document.get(), object, "thermal_authority",
          thermal_authority_name(*boundary.thermal_authority)));
    }
    if (boundary.temperature_K.has_value()) {
      require_json_write(yyjson_mut_obj_add_real(
          document.get(), object, "temperature_K", *boundary.temperature_K));
    }
    if (boundary.enthalpy_J_per_kg.has_value()) {
      require_json_write(yyjson_mut_obj_add_real(
          document.get(), object, "enthalpy_J_per_kg",
          *boundary.enthalpy_J_per_kg));
    }
    if (boundary.density_kg_per_m3.has_value()) {
      require_json_write(yyjson_mut_obj_add_real(
          document.get(), object, "density_kg_per_m3",
          *boundary.density_kg_per_m3));
    }
    if (boundary.scalar_values.has_value()) {
      yyjson_mut_val* values = add_array(document.get(), object, "scalar_values");
      for (const auto& scalar_value : *boundary.scalar_values) {
        yyjson_mut_val* value_object = yyjson_mut_obj(document.get());
        if (value_object == nullptr ||
            !yyjson_mut_arr_append(values, value_object)) {
          throw Error("failed to allocate deterministic resolved JSON");
        }
        add_string(document.get(), value_object, "name", scalar_value.name);
        require_json_write(yyjson_mut_obj_add_real(
            document.get(), value_object, "value", scalar_value.value));
      }
    }
    if (boundary.pressure_perturbation_pa.has_value()) {
      require_json_write(yyjson_mut_obj_add_real(
          document.get(), object, "pressure_perturbation_pa",
          *boundary.pressure_perturbation_pa));
    }
  }

  yyjson_mut_val* restart = add_object(document.get(), root, "restart");
  require_json_write(yyjson_mut_obj_add_bool(document.get(), restart, "read",
                                             config.restart.read));
  if (config.restart.read_directory.has_value()) {
    add_string(document.get(), restart, "read_directory",
               config.restart.read_directory->generic_string());
  }
  add_string(document.get(), restart, "write_directory",
             config.restart.write_directory.generic_string());
  require_json_write(yyjson_mut_obj_add_int(
      document.get(), restart, "write_interval",
      config.restart.write_interval));

  yyjson_mut_val* diagnostics =
      add_object(document.get(), root, "diagnostics");
  add_string(document.get(), diagnostics, "directory",
             config.diagnostics.directory.generic_string());
  require_json_write(yyjson_mut_obj_add_int(
      document.get(), diagnostics, "write_interval",
      config.diagnostics.write_interval));
  require_json_write(yyjson_mut_obj_add_bool(
      document.get(), diagnostics, "write_mesh",
      config.diagnostics.write_mesh));

  yyjson_mut_val* performance =
      add_object(document.get(), root, "performance");
  require_json_write(yyjson_mut_obj_add_bool(
      document.get(), performance, "enabled", config.performance.enabled));
  add_string(document.get(), performance, "directory",
             config.performance.directory.generic_string());
  require_json_write(yyjson_mut_obj_add_int(
      document.get(), performance, "warmup_steps",
      config.performance.warmup_steps));
  require_json_write(yyjson_mut_obj_add_int(
      document.get(), performance, "measured_steps",
      config.performance.measured_steps));
  require_json_write(yyjson_mut_obj_add_int(
      document.get(), performance, "repetitions",
      config.performance.repetitions));

  std::size_t length = 0U;
  using Buffer = std::unique_ptr<char, decltype(&std::free)>;
  Buffer buffer(
      yyjson_mut_write(document.get(), YYJSON_WRITE_NOFLAG, &length),
      &std::free);
  if (!buffer) {
    throw Error("failed to write deterministic resolved JSON");
  }
  return std::string(buffer.get(), length);
}

int dispatch_schema_version(yyjson_val* root) {
  if (!yyjson_is_obj(root)) {
    throw_type_error("", "object", root);
  }
  int occurrences = 0;
  yyjson_obj_iter iterator = yyjson_obj_iter_with(root);
  while (yyjson_val* key_value = yyjson_obj_iter_next(&iterator)) {
    const std::string_view key(yyjson_get_str(key_value),
                               yyjson_get_len(key_value));
    if (key == "schema_version") {
      ++occurrences;
      if (occurrences > 1) {
        throw ConfigError("/schema_version", "duplicate member");
      }
    }
  }
  return read_integer(require_member(root, "schema_version", ""),
                      "/schema_version");
}

}  // namespace

ResolvedCase load_resolved_case(const std::filesystem::path& path) {
  const std::string contents = read_complete_file(path);
  Document document = parse_document(contents);
  yyjson_val* root = yyjson_doc_get_root(document.get());
  const int version = dispatch_schema_version(root);
  if (version == 1) {
    return ResolvedCase(detail::load_case_config_from_json(contents, path));
  }
  if (version == 2) {
    return ResolvedCase(parse_flow_case(root, path));
  }
  throw ConfigError("/schema_version",
                    "unsupported schema version " + std::to_string(version));
}

std::string to_resolved_json(const ResolvedCase& resolved_case) {
  return std::visit(
      [](const auto& config) -> std::string {
        using Config = std::decay_t<decltype(config)>;
        if constexpr (std::is_same_v<Config, CaseConfig>) {
          return to_resolved_json(config);
        } else {
          return serialize_flow_case(config);
        }
      },
      resolved_case);
}

}  // namespace hundun::config
