// SPDX-License-Identifier: Apache-2.0

#include "hundun/cfg_resolved_case_v3_loader.hpp"

#include "cfg_case_config_loader_detail.hpp"
#include "cfg_resolved_case_v3_loader_detail.hpp"

#include "hundun/cfg_resolved_case_loader.hpp"
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
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace hundun::config {
namespace {

using runtime::ConfigError;
using runtime::Error;

using Document = std::unique_ptr<yyjson_doc, decltype(&yyjson_doc_free)>;
using MutableDocument =
    std::unique_ptr<yyjson_mut_doc, decltype(&yyjson_mut_doc_free)>;

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

std::string actual_type(yyjson_val *value) {
  return value == nullptr ? "missing" : yyjson_get_type_desc(value);
}

[[noreturn]] void throw_type_error(std::string pointer,
                                   std::string_view expected,
                                   yyjson_val *actual) {
  throw ConfigError(std::move(pointer), "expected " + std::string(expected) +
                                            ", got " + actual_type(actual));
}

yyjson_val *require_member(yyjson_val *object, const char *key,
                           std::string_view pointer) {
  yyjson_val *value = yyjson_obj_get(object, key);
  if (value == nullptr) {
    throw ConfigError(child_pointer(pointer, key),
                      "required member is missing");
  }
  return value;
}

yyjson_val *require_object_member(yyjson_val *object, const char *key,
                                  std::string_view pointer) {
  yyjson_val *value = require_member(object, key, pointer);
  if (!yyjson_is_obj(value)) {
    throw_type_error(child_pointer(pointer, key), "object", value);
  }
  return value;
}

yyjson_val *require_array_member(yyjson_val *object, const char *key,
                                 std::string_view pointer) {
  yyjson_val *value = require_member(object, key, pointer);
  if (!yyjson_is_arr(value)) {
    throw_type_error(child_pointer(pointer, key), "array", value);
  }
  return value;
}

void reject_unknown_keys(yyjson_val *object,
                         std::initializer_list<std::string_view> allowed,
                         std::string_view pointer) {
  std::vector<std::string> seen;
  seen.reserve(yyjson_obj_size(object));
  yyjson_obj_iter iterator = yyjson_obj_iter_with(object);
  while (yyjson_val *key_value = yyjson_obj_iter_next(&iterator)) {
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

double read_number(yyjson_val *value, std::string pointer) {
  if (!yyjson_is_num(value)) {
    throw_type_error(std::move(pointer), "number", value);
  }
  const double result = yyjson_get_num(value);
  if (!std::isfinite(result)) {
    throw ConfigError(std::move(pointer), "expected a finite number");
  }
  return result;
}

int read_integer(yyjson_val *value, std::string pointer) {
  if (!yyjson_is_int(value)) {
    throw_type_error(std::move(pointer), "integer", value);
  }
  if (yyjson_is_uint(value)) {
    const std::uint64_t result = yyjson_get_uint(value);
    if (result > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
      throw ConfigError(std::move(pointer),
                        "integer exceeds the C++ int range");
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

std::string read_string(yyjson_val *value, std::string pointer) {
  if (!yyjson_is_str(value)) {
    throw_type_error(std::move(pointer), "string", value);
  }
  return std::string(yyjson_get_str(value), yyjson_get_len(value));
}

bool read_boolean(yyjson_val *value, std::string pointer) {
  if (!yyjson_is_bool(value)) {
    throw_type_error(std::move(pointer), "boolean", value);
  }
  return yyjson_get_bool(value);
}

double require_number(yyjson_val *object, const char *key,
                      std::string_view pointer) {
  return read_number(require_member(object, key, pointer),
                     child_pointer(pointer, key));
}

int require_integer(yyjson_val *object, const char *key,
                    std::string_view pointer) {
  return read_integer(require_member(object, key, pointer),
                      child_pointer(pointer, key));
}

std::string require_string(yyjson_val *object, const char *key,
                           std::string_view pointer) {
  return read_string(require_member(object, key, pointer),
                     child_pointer(pointer, key));
}

bool require_boolean(yyjson_val *object, const char *key,
                     std::string_view pointer) {
  return read_boolean(require_member(object, key, pointer),
                      child_pointer(pointer, key));
}

runtime::Int3 require_int3(yyjson_val *object, const char *key,
                           std::string_view pointer) {
  const std::string field_pointer = child_pointer(pointer, key);
  yyjson_val *array = require_member(object, key, pointer);
  if (!yyjson_is_arr(array)) {
    throw_type_error(field_pointer, "array of three integers", array);
  }
  if (yyjson_arr_size(array) != 3U) {
    throw ConfigError(field_pointer, "expected exactly three entries");
  }
  return {
      read_integer(yyjson_arr_get(array, 0U), index_pointer(field_pointer, 0U)),
      read_integer(yyjson_arr_get(array, 1U), index_pointer(field_pointer, 1U)),
      read_integer(yyjson_arr_get(array, 2U),
                   index_pointer(field_pointer, 2U))};
}

runtime::Real3 require_real3(yyjson_val *object, const char *key,
                             std::string_view pointer) {
  const std::string field_pointer = child_pointer(pointer, key);
  yyjson_val *array = require_member(object, key, pointer);
  if (!yyjson_is_arr(array)) {
    throw_type_error(field_pointer, "array of three numbers", array);
  }
  if (yyjson_arr_size(array) != 3U) {
    throw ConfigError(field_pointer, "expected exactly three entries");
  }
  return {
      read_number(yyjson_arr_get(array, 0U), index_pointer(field_pointer, 0U)),
      read_number(yyjson_arr_get(array, 1U), index_pointer(field_pointer, 1U)),
      read_number(yyjson_arr_get(array, 2U), index_pointer(field_pointer, 2U))};
}

std::optional<double> optional_number(yyjson_val *object, const char *key,
                                      std::string_view pointer) {
  yyjson_val *value = yyjson_obj_get(object, key);
  if (value == nullptr) {
    return std::nullopt;
  }
  return read_number(value, child_pointer(pointer, key));
}

Document parse_document(std::string &contents) {
  yyjson_read_err read_error{};
  Document document(yyjson_read_opts(contents.data(), contents.size(),
                                     YYJSON_READ_NOFLAG, nullptr, &read_error),
                    &yyjson_doc_free);
  if (!document) {
    const std::string detail =
        read_error.msg == nullptr ? "unknown parse error" : read_error.msg;
    throw ConfigError("", "invalid JSON at byte " +
                              std::to_string(read_error.pos) + ": " + detail);
  }
  return document;
}

std::string read_complete_file(const std::filesystem::path &path) {
  std::ifstream stream(path, std::ios::binary);
  if (!stream) {
    throw ConfigError("", "could not open case JSON: " + path.string());
  }
  const std::string contents((std::istreambuf_iterator<char>(stream)),
                             std::istreambuf_iterator<char>());
  if (stream.bad()) {
    throw ConfigError("",
                      "could not read complete case JSON: " + path.string());
  }
  return contents;
}

bool has_parent_component(const std::filesystem::path &path) {
  return std::any_of(path.begin(), path.end(),
                     [](const auto &component) { return component == ".."; });
}

bool has_path_prefix(const std::filesystem::path &path,
                     const std::filesystem::path &prefix) {
  auto path_iterator = path.begin();
  for (auto prefix_iterator = prefix.begin(); prefix_iterator != prefix.end();
       ++prefix_iterator, ++path_iterator) {
    if (path_iterator == path.end() || *path_iterator != *prefix_iterator) {
      return false;
    }
  }
  return true;
}

std::filesystem::path
normalize_case_path(const std::string &raw,
                    const std::filesystem::path &case_path,
                    std::string_view pointer) {
  if (raw.empty()) {
    throw ConfigError(std::string(pointer),
                      "expected a non-empty relative path");
  }
  if (raw.find('\0') != std::string::npos) {
    throw ConfigError(std::string(pointer),
                      "path must not contain a null byte");
  }
  const std::filesystem::path path(raw);
  if (path.is_absolute() || path.has_root_name() || path.has_root_directory()) {
    throw ConfigError(std::string(pointer),
                      "expected an unrooted relative path");
  }
  if (has_parent_component(path)) {
    throw ConfigError(std::string(pointer),
                      "path must not contain a '..' component");
  }
  const std::filesystem::path normalized = path.lexically_normal();
  if (normalized.empty() || normalized == ".") {
    throw ConfigError(std::string(pointer),
                      "expected a non-empty relative path");
  }
  try {
    const std::filesystem::path case_root =
        std::filesystem::absolute(case_path).lexically_normal().parent_path();
    const std::filesystem::path resolved =
        (case_root / normalized).lexically_normal();
    if (!has_path_prefix(resolved, case_root)) {
      throw ConfigError(std::string(pointer),
                        "path escapes the case directory");
    }
  } catch (const std::filesystem::filesystem_error &error) {
    throw ConfigError(std::string(pointer),
                      "could not resolve path against the case directory: " +
                          std::string(error.what()));
  }
  return normalized;
}

void validate_case_path(const std::filesystem::path &path,
                        std::string_view pointer) {
  const std::string text = path.generic_string();
  if (text.empty() || text == ".") {
    throw ConfigError(std::string(pointer),
                      "expected a non-empty relative path");
  }
  if (text.find('\0') != std::string::npos) {
    throw ConfigError(std::string(pointer),
                      "path must not contain a null byte");
  }
  if (path.is_absolute() || path.has_root_name() || path.has_root_directory()) {
    throw ConfigError(std::string(pointer),
                      "expected an unrooted relative path");
  }
  if (has_parent_component(path)) {
    throw ConfigError(std::string(pointer),
                      "path must not contain a '..' component");
  }
  if (text != path.lexically_normal().generic_string()) {
    throw ConfigError(std::string(pointer),
                      "path must be lexically normalized");
  }
}

DensityModel parse_density_model(const std::string &value,
                                 const std::string &pointer) {
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

MeshMapping parse_mesh_mapping(const std::string &value,
                               const std::string &pointer) {
  if (value == "uniform_box") {
    return MeshMapping::uniform_box;
  }
  if (value == "analytic_warped_box") {
    return MeshMapping::analytic_warped_box;
  }
  throw ConfigError(pointer, "unknown mesh mapping '" + value + "'");
}

TimeMode parse_time_mode(const std::string &value, const std::string &pointer) {
  if (value == "fixed") {
    return TimeMode::fixed;
  }
  if (value == "adaptive") {
    return TimeMode::adaptive;
  }
  throw ConfigError(pointer, "unknown time mode '" + value + "'");
}

PatchName parse_patch_name(const std::string &value,
                           const std::string &pointer) {
  static constexpr std::array<std::string_view, 6> names = {
      "x_min", "x_max", "y_min", "y_max", "z_min", "z_max"};
  const auto found = std::find(names.begin(), names.end(), value);
  if (found == names.end()) {
    throw ConfigError(pointer, "unknown boundary patch '" + value + "'");
  }
  return static_cast<PatchName>(std::distance(names.begin(), found));
}

BoundaryType parse_boundary_type(const std::string &value,
                                 const std::string &pointer) {
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

InletThermalAuthority parse_thermal_authority(const std::string &value,
                                              const std::string &pointer) {
  if (value == "temperature") {
    return InletThermalAuthority::temperature;
  }
  if (value == "enthalpy") {
    return InletThermalAuthority::enthalpy;
  }
  throw ConfigError(pointer, "unknown inlet thermal authority '" + value + "'");
}

void canonicalize_common_order(FlowCaseConfig &flow) {
  std::sort(flow.scalars.begin(), flow.scalars.end(),
            [](const auto &left, const auto &right) {
              return left.name < right.name;
            });
  for (auto &boundary : flow.boundaries) {
    if (boundary.scalar_values.has_value()) {
      std::sort(boundary.scalar_values->begin(), boundary.scalar_values->end(),
                [](const auto &left, const auto &right) {
                  return left.name < right.name;
                });
    }
  }
  std::sort(flow.boundaries.begin(), flow.boundaries.end(),
            [](const auto &left, const auto &right) {
              return static_cast<int>(left.patch) <
                     static_cast<int>(right.patch);
            });
}

FlowCaseConfig parse_common_flow(yyjson_val *root,
                                 const std::filesystem::path &case_path,
                                 bool immersed_flow_root) {
  if (!yyjson_is_obj(root)) {
    throw_type_error("", "object", root);
  }
  if (immersed_flow_root) {
    reject_unknown_keys(root,
                        {"schema_version", "case", "simulation", "resources",
                         "mesh", "time", "physics", "scalars", "boundaries",
                         "restart", "diagnostics", "performance",
                         "immersed_boundary", "les"},
                        "");
  } else {
    reject_unknown_keys(root,
                        {"schema_version", "case", "simulation", "resources",
                         "mesh", "time", "physics", "scalars", "boundaries",
                         "restart", "diagnostics", "performance"},
                        "");
  }

  FlowCaseConfig config{};
  config.schema_version = 2;

  yyjson_val *case_object = require_object_member(root, "case", "");
  reject_unknown_keys(case_object, {"name"}, "/case");
  config.case_name = require_string(case_object, "name", "/case");

  yyjson_val *simulation = require_object_member(root, "simulation", "");
  reject_unknown_keys(simulation, {"type", "density_model"}, "/simulation");
  const std::string simulation_type =
      require_string(simulation, "type", "/simulation");
  if (simulation_type != "variable_density_flow") {
    throw ConfigError("/simulation/type", "expected 'variable_density_flow'");
  }
  config.simulation_type = SimulationType::variable_density_flow;
  config.density_model = parse_density_model(
      require_string(simulation, "density_model", "/simulation"),
      "/simulation/density_model");

  if (yyjson_val *resources = yyjson_obj_get(root, "resources")) {
    if (!yyjson_is_obj(resources)) {
      throw_type_error("/resources", "object", resources);
    }
    reject_unknown_keys(resources, {"expected_ranks", "process_grid"},
                        "/resources");
    if (yyjson_val *expected = yyjson_obj_get(resources, "expected_ranks")) {
      config.resources.expected_ranks =
          read_integer(expected, "/resources/expected_ranks");
    }
    if (yyjson_obj_get(resources, "process_grid") != nullptr) {
      config.resources.process_grid =
          require_int3(resources, "process_grid", "/resources");
    }
  }

  yyjson_val *mesh = require_object_member(root, "mesh", "");
  reject_unknown_keys(
      mesh, {"cells", "origin_m", "length_m", "mapping", "warp_amplitude"},
      "/mesh");
  config.mesh.cells = require_int3(mesh, "cells", "/mesh");
  config.mesh.origin_m = require_real3(mesh, "origin_m", "/mesh");
  config.mesh.length_m = require_real3(mesh, "length_m", "/mesh");
  config.mesh.mapping = parse_mesh_mapping(
      require_string(mesh, "mapping", "/mesh"), "/mesh/mapping");
  if (yyjson_obj_get(mesh, "warp_amplitude") != nullptr) {
    config.mesh.warp_amplitude = require_real3(mesh, "warp_amplitude", "/mesh");
  }

  yyjson_val *time = require_object_member(root, "time", "");
  reject_unknown_keys(time,
                      {"mode", "steps", "initial_dt_s", "min_dt_s", "max_dt_s",
                       "cfl_target", "diffusion_number_target", "growth_factor",
                       "retry_factor", "max_retries"},
                      "/time");
  config.time.mode =
      parse_time_mode(require_string(time, "mode", "/time"), "/time/mode");
  config.time.steps = require_integer(time, "steps", "/time");
  config.time.initial_dt_s = require_number(time, "initial_dt_s", "/time");
  config.time.min_dt_s = require_number(time, "min_dt_s", "/time");
  config.time.max_dt_s = require_number(time, "max_dt_s", "/time");
  config.time.cfl_target = require_number(time, "cfl_target", "/time");
  config.time.diffusion_number_target =
      require_number(time, "diffusion_number_target", "/time");
  config.time.growth_factor = require_number(time, "growth_factor", "/time");
  config.time.retry_factor = require_number(time, "retry_factor", "/time");
  config.time.max_retries = require_integer(time, "max_retries", "/time");

  yyjson_val *physics = require_object_member(root, "physics", "");
  reject_unknown_keys(physics,
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

  yyjson_val *scalars = require_array_member(root, "scalars", "");
  config.scalars.reserve(yyjson_arr_size(scalars));
  for (std::size_t index = 0; index < yyjson_arr_size(scalars); ++index) {
    yyjson_val *scalar = yyjson_arr_get(scalars, index);
    const std::string pointer = index_pointer("/scalars", index);
    if (!yyjson_is_obj(scalar)) {
      throw_type_error(pointer, "object", scalar);
    }
    reject_unknown_keys(scalar, {"name", "diffusivity_m2_per_s"}, pointer);
    config.scalars.push_back(
        {require_string(scalar, "name", pointer),
         require_number(scalar, "diffusivity_m2_per_s", pointer)});
  }

  yyjson_val *boundaries = require_array_member(root, "boundaries", "");
  if (yyjson_arr_size(boundaries) != config.boundaries.size()) {
    throw ConfigError("/boundaries", "expected exactly six boundary entries");
  }
  for (std::size_t index = 0; index < config.boundaries.size(); ++index) {
    yyjson_val *boundary_object = yyjson_arr_get(boundaries, index);
    const std::string pointer = index_pointer("/boundaries", index);
    if (!yyjson_is_obj(boundary_object)) {
      throw_type_error(pointer, "object", boundary_object);
    }
    reject_unknown_keys(boundary_object,
                        {"patch", "type", "velocity_m_per_s",
                         "thermal_authority", "temperature_K",
                         "enthalpy_J_per_kg", "density_kg_per_m3",
                         "scalar_values", "pressure_perturbation_pa"},
                        pointer);
    auto &boundary = config.boundaries[index];
    boundary.patch =
        parse_patch_name(require_string(boundary_object, "patch", pointer),
                         child_pointer(pointer, "patch"));
    boundary.type =
        parse_boundary_type(require_string(boundary_object, "type", pointer),
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
    boundary.pressure_perturbation_pa =
        optional_number(boundary_object, "pressure_perturbation_pa", pointer);
    if (yyjson_val *scalar_values =
            yyjson_obj_get(boundary_object, "scalar_values")) {
      if (!yyjson_is_arr(scalar_values)) {
        throw_type_error(child_pointer(pointer, "scalar_values"), "array",
                         scalar_values);
      }
      boundary.scalar_values.emplace();
      boundary.scalar_values->reserve(yyjson_arr_size(scalar_values));
      for (std::size_t value_index = 0;
           value_index < yyjson_arr_size(scalar_values); ++value_index) {
        yyjson_val *scalar_value = yyjson_arr_get(scalar_values, value_index);
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

  yyjson_val *restart = require_object_member(root, "restart", "");
  reject_unknown_keys(
      restart, {"read", "read_directory", "write_directory", "write_interval"},
      "/restart");
  config.restart.read = require_boolean(restart, "read", "/restart");
  if (yyjson_val *read_directory = yyjson_obj_get(restart, "read_directory")) {
    config.restart.read_directory = normalize_case_path(
        read_string(read_directory, "/restart/read_directory"), case_path,
        "/restart/read_directory");
  }
  config.restart.write_directory = normalize_case_path(
      require_string(restart, "write_directory", "/restart"), case_path,
      "/restart/write_directory");
  config.restart.write_interval =
      require_integer(restart, "write_interval", "/restart");

  yyjson_val *diagnostics = require_object_member(root, "diagnostics", "");
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

  yyjson_val *performance = require_object_member(root, "performance", "");
  reject_unknown_keys(
      performance,
      {"enabled", "directory", "warmup_steps", "measured_steps", "repetitions"},
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

  static_cast<void>(to_resolved_json(ResolvedCase(config)));
  canonicalize_common_order(config);
  return config;
}

ImmersedBoundaryModel parse_ibm_model(const std::string &value) {
  if (value == "none") {
    return ImmersedBoundaryModel::none;
  }
  if (value == "local_flow_pattern_ghost_cell") {
    return ImmersedBoundaryModel::local_flow_pattern_ghost_cell;
  }
  throw ConfigError("/immersed_boundary/model",
                    "unknown immersed-boundary model '" + value + "'");
}

ImmersedFluidSide parse_fluid_side(const std::string &value) {
  if (value == "inside") {
    return ImmersedFluidSide::inside;
  }
  if (value == "outside") {
    return ImmersedFluidSide::outside;
  }
  throw ConfigError("/immersed_boundary/geometry/fluid_side",
                    "unknown fluid side '" + value + "'");
}

LesModel parse_les_model(const std::string &value) {
  if (value == "none") {
    return LesModel::none;
  }
  if (value == "wale") {
    return LesModel::wale;
  }
  throw ConfigError("/les/model", "unknown LES model '" + value + "'");
}

void require_range(double value, double lower, double upper,
                   std::string pointer) {
  if (!std::isfinite(value) || value < lower || value > upper) {
    throw ConfigError(std::move(pointer),
                      "value is outside the approved inclusive range");
  }
}

ImmersedFlowCaseConfig canonicalize_immersed_flow(ImmersedFlowCaseConfig config) {
  if (config.schema_version != 3) {
    throw ConfigError("/schema_version", "expected schema version 3");
  }
  if (config.common_flow.schema_version != 2) {
    throw ConfigError("/schema_version",
                      "common flow must use internal schema version 2");
  }
  static_cast<void>(to_resolved_json(ResolvedCase(config.common_flow)));
  canonicalize_common_order(config.common_flow);

  switch (config.immersed_boundary.model) {
  case ImmersedBoundaryModel::none:
    if (config.immersed_boundary.geometry.has_value()) {
      throw ConfigError("/immersed_boundary/geometry",
                        "member is forbidden when the model is none");
    }
    if (config.immersed_boundary.wall.has_value()) {
      throw ConfigError("/immersed_boundary/wall",
                        "member is forbidden when the model is none");
    }
    break;
  case ImmersedBoundaryModel::local_flow_pattern_ghost_cell: {
    if (!config.immersed_boundary.geometry.has_value()) {
      throw ConfigError("/immersed_boundary/geometry",
                        "required member is missing");
    }
    if (!config.immersed_boundary.wall.has_value()) {
      throw ConfigError("/immersed_boundary/wall",
                        "required member is missing");
    }
    auto &geometry = *config.immersed_boundary.geometry;
    validate_case_path(geometry.file, "/immersed_boundary/geometry/file");
    if (!std::isfinite(geometry.length_scale_to_m) ||
        geometry.length_scale_to_m <= 0.0) {
      throw ConfigError("/immersed_boundary/geometry/length_scale_to_m",
                        "expected a finite value greater than zero");
    }
    if (geometry.fluid_side != ImmersedFluidSide::inside &&
        geometry.fluid_side != ImmersedFluidSide::outside) {
      throw ConfigError("/immersed_boundary/geometry/fluid_side",
                        "invalid fluid side");
    }
    auto &velocity = config.immersed_boundary.wall->velocity_m_per_s;
    const std::array<double, 3> components{velocity.x, velocity.y, velocity.z};
    for (std::size_t index = 0; index < components.size(); ++index) {
      if (!std::isfinite(components[index]) || components[index] != 0.0) {
        throw ConfigError(
            index_pointer("/immersed_boundary/wall/velocity_m_per_s", index),
            "immersed-flow requires a mathematically zero wall velocity");
      }
    }
    velocity = {0.0, 0.0, 0.0};
    break;
  }
  default:
    throw ConfigError("/immersed_boundary/model",
                      "invalid immersed-boundary model");
  }

  switch (config.les.model) {
  case LesModel::none:
    if (config.les.wale.has_value()) {
      throw ConfigError("/les/wale",
                        "member is forbidden when the model is none");
    }
    break;
  case LesModel::wale:
    if (!config.les.wale.has_value()) {
      throw ConfigError("/les/wale", "required member is missing");
    }
    require_range(config.les.wale->coefficient, 1.0e-6, 1.0,
                  "/les/wale/coefficient");
    require_range(config.les.wale->turbulent_prandtl, 0.1, 10.0,
                  "/les/turbulent_prandtl");
    require_range(config.les.wale->turbulent_schmidt, 0.1, 10.0,
                  "/les/turbulent_schmidt");
    break;
  default:
    throw ConfigError("/les/model", "invalid LES model");
  }

  if (config.immersed_boundary.model == ImmersedBoundaryModel::none &&
      config.les.model == LesModel::none) {
    throw ConfigError("/les/model", "schema version 3 requires IBM or LES");
  }
  return config;
}

ImmersedBoundaryConfig parse_ibm(yyjson_val *root,
                                 const std::filesystem::path &case_path) {
  yyjson_val *object = require_object_member(root, "immersed_boundary", "");
  const ImmersedBoundaryModel model =
      parse_ibm_model(require_string(object, "model", "/immersed_boundary"));

  ImmersedBoundaryConfig result{};
  result.model = model;
  if (model == ImmersedBoundaryModel::none) {
    reject_unknown_keys(object, {"model"}, "/immersed_boundary");
    return result;
  }

  reject_unknown_keys(object, {"model", "geometry", "wall"},
                      "/immersed_boundary");
  yyjson_val *geometry =
      require_object_member(object, "geometry", "/immersed_boundary");
  reject_unknown_keys(geometry,
                      {"format", "file", "length_scale_to_m", "fluid_side"},
                      "/immersed_boundary/geometry");
  const std::string format =
      require_string(geometry, "format", "/immersed_boundary/geometry");
  if (format != "stl") {
    throw ConfigError("/immersed_boundary/geometry/format", "expected 'stl'");
  }
  StlGeometryConfig geometry_config{};
  geometry_config.file = normalize_case_path(
      require_string(geometry, "file", "/immersed_boundary/geometry"),
      case_path, "/immersed_boundary/geometry/file");
  geometry_config.length_scale_to_m = require_number(
      geometry, "length_scale_to_m", "/immersed_boundary/geometry");
  if (geometry_config.length_scale_to_m <= 0.0) {
    throw ConfigError("/immersed_boundary/geometry/length_scale_to_m",
                      "expected a finite value greater than zero");
  }
  geometry_config.fluid_side = parse_fluid_side(
      require_string(geometry, "fluid_side", "/immersed_boundary/geometry"));
  result.geometry = std::move(geometry_config);

  yyjson_val *wall =
      require_object_member(object, "wall", "/immersed_boundary");
  reject_unknown_keys(wall, {"velocity_m_per_s", "enthalpy", "scalars"},
                      "/immersed_boundary/wall");
  StaticImmersedWallConfig wall_config{};
  wall_config.velocity_m_per_s =
      require_real3(wall, "velocity_m_per_s", "/immersed_boundary/wall");
  const std::array<double, 3> velocity{wall_config.velocity_m_per_s.x,
                                       wall_config.velocity_m_per_s.y,
                                       wall_config.velocity_m_per_s.z};
  for (std::size_t index = 0; index < velocity.size(); ++index) {
    if (velocity[index] != 0.0) {
      throw ConfigError(
          index_pointer("/immersed_boundary/wall/velocity_m_per_s", index),
          "immersed-flow requires a mathematically zero wall velocity");
    }
  }
  const std::string enthalpy =
      require_string(wall, "enthalpy", "/immersed_boundary/wall");
  if (enthalpy != "zero_normal_diffusive_flux") {
    throw ConfigError("/immersed_boundary/wall/enthalpy",
                      "expected 'zero_normal_diffusive_flux'");
  }
  const std::string scalars =
      require_string(wall, "scalars", "/immersed_boundary/wall");
  if (scalars != "zero_normal_diffusive_flux") {
    throw ConfigError("/immersed_boundary/wall/scalars",
                      "expected 'zero_normal_diffusive_flux'");
  }
  wall_config.velocity_m_per_s = {0.0, 0.0, 0.0};
  result.wall = wall_config;
  return result;
}

LesConfig parse_les(yyjson_val *root) {
  yyjson_val *object = require_object_member(root, "les", "");
  const LesModel model =
      parse_les_model(require_string(object, "model", "/les"));
  LesConfig result{};
  result.model = model;
  if (model == LesModel::none) {
    reject_unknown_keys(object, {"model"}, "/les");
    return result;
  }

  reject_unknown_keys(
      object, {"model", "wale", "turbulent_prandtl", "turbulent_schmidt"},
      "/les");
  yyjson_val *wale = require_object_member(object, "wale", "/les");
  reject_unknown_keys(wale, {"coefficient"}, "/les/wale");
  WaleConfig config{};
  config.coefficient = require_number(wale, "coefficient", "/les/wale");
  config.turbulent_prandtl =
      require_number(object, "turbulent_prandtl", "/les");
  config.turbulent_schmidt =
      require_number(object, "turbulent_schmidt", "/les");
  require_range(config.coefficient, 1.0e-6, 1.0, "/les/wale/coefficient");
  require_range(config.turbulent_prandtl, 0.1, 10.0, "/les/turbulent_prandtl");
  require_range(config.turbulent_schmidt, 0.1, 10.0, "/les/turbulent_schmidt");
  result.wale = config;
  return result;
}

ImmersedFlowCaseConfig parse_immersed_flow(yyjson_val *root,
                                    const std::filesystem::path &case_path) {
  const int version = require_integer(root, "schema_version", "");
  if (version != 3) {
    throw ConfigError("/schema_version", "expected schema version 3");
  }
  ImmersedFlowCaseConfig config{};
  config.schema_version = 3;
  config.common_flow = parse_common_flow(root, case_path, true);
  config.immersed_boundary = parse_ibm(root, case_path);
  config.les = parse_les(root);
  return canonicalize_immersed_flow(std::move(config));
}

int dispatch_schema_version(yyjson_val *root) {
  if (!yyjson_is_obj(root)) {
    throw_type_error("", "object", root);
  }
  int occurrences = 0;
  yyjson_obj_iter iterator = yyjson_obj_iter_with(root);
  while (yyjson_val *key_value = yyjson_obj_iter_next(&iterator)) {
    const std::string_view key(yyjson_get_str(key_value),
                               yyjson_get_len(key_value));
    if (key == "schema_version") {
      ++occurrences;
      if (occurrences > 1) {
        throw ConfigError("/schema_version", "duplicate member");
      }
    }
  }
  return require_integer(root, "schema_version", "");
}

void require_json_write(bool successful) {
  if (!successful) {
    throw Error("failed to allocate deterministic immersed-flow resolved JSON");
  }
}

yyjson_mut_val *add_object(yyjson_mut_doc *document, yyjson_mut_val *parent,
                           const char *key) {
  yyjson_mut_val *object = yyjson_mut_obj_add_obj(document, parent, key);
  if (object == nullptr) {
    throw Error("failed to allocate deterministic immersed-flow resolved JSON");
  }
  return object;
}

void add_string(yyjson_mut_doc *document, yyjson_mut_val *object,
                const char *key, std::string_view value) {
  require_json_write(yyjson_mut_obj_add_strncpy(document, object, key,
                                                value.data(), value.size()));
}

const char *fluid_side_name(ImmersedFluidSide side) {
  switch (side) {
  case ImmersedFluidSide::inside:
    return "inside";
  case ImmersedFluidSide::outside:
    return "outside";
  }
  throw ConfigError("/immersed_boundary/geometry/fluid_side",
                    "invalid fluid side");
}

std::string serialize_immersed_flow_tail(const ImmersedFlowCaseConfig &input) {
  const ImmersedFlowCaseConfig config = canonicalize_immersed_flow(input);
  MutableDocument document(yyjson_mut_doc_new(nullptr), &yyjson_mut_doc_free);
  if (!document) {
    throw Error("failed to allocate deterministic immersed-flow resolved JSON");
  }
  yyjson_mut_val *root = yyjson_mut_obj(document.get());
  if (root == nullptr) {
    throw Error("failed to allocate deterministic immersed-flow resolved JSON");
  }
  yyjson_mut_doc_set_root(document.get(), root);

  yyjson_mut_val *ibm = add_object(document.get(), root, "immersed_boundary");
  if (config.immersed_boundary.model == ImmersedBoundaryModel::none) {
    require_json_write(
        yyjson_mut_obj_add_str(document.get(), ibm, "model", "none"));
  } else {
    require_json_write(yyjson_mut_obj_add_str(document.get(), ibm, "model",
                                              "local_flow_pattern_ghost_cell"));
    const auto &geometry = *config.immersed_boundary.geometry;
    yyjson_mut_val *geometry_object =
        add_object(document.get(), ibm, "geometry");
    require_json_write(yyjson_mut_obj_add_str(document.get(), geometry_object,
                                              "format", "stl"));
    add_string(document.get(), geometry_object, "file",
               geometry.file.generic_string());
    require_json_write(yyjson_mut_obj_add_real(document.get(), geometry_object,
                                               "length_scale_to_m",
                                               geometry.length_scale_to_m));
    require_json_write(
        yyjson_mut_obj_add_str(document.get(), geometry_object, "fluid_side",
                               fluid_side_name(geometry.fluid_side)));

    yyjson_mut_val *wall = add_object(document.get(), ibm, "wall");
    yyjson_mut_val *velocity =
        yyjson_mut_obj_add_arr(document.get(), wall, "velocity_m_per_s");
    if (velocity == nullptr) {
      throw Error("failed to allocate deterministic immersed-flow resolved JSON");
    }
    require_json_write(yyjson_mut_arr_add_real(document.get(), velocity, 0.0));
    require_json_write(yyjson_mut_arr_add_real(document.get(), velocity, 0.0));
    require_json_write(yyjson_mut_arr_add_real(document.get(), velocity, 0.0));
    require_json_write(yyjson_mut_obj_add_str(document.get(), wall, "enthalpy",
                                              "zero_normal_diffusive_flux"));
    require_json_write(yyjson_mut_obj_add_str(document.get(), wall, "scalars",
                                              "zero_normal_diffusive_flux"));
  }

  yyjson_mut_val *les = add_object(document.get(), root, "les");
  if (config.les.model == LesModel::none) {
    require_json_write(
        yyjson_mut_obj_add_str(document.get(), les, "model", "none"));
  } else {
    require_json_write(
        yyjson_mut_obj_add_str(document.get(), les, "model", "wale"));
    yyjson_mut_val *wale = add_object(document.get(), les, "wale");
    require_json_write(yyjson_mut_obj_add_real(
        document.get(), wale, "coefficient", config.les.wale->coefficient));
    require_json_write(
        yyjson_mut_obj_add_real(document.get(), les, "turbulent_prandtl",
                                config.les.wale->turbulent_prandtl));
    require_json_write(
        yyjson_mut_obj_add_real(document.get(), les, "turbulent_schmidt",
                                config.les.wale->turbulent_schmidt));
  }

  std::size_t length = 0U;
  using Buffer = std::unique_ptr<char, decltype(&std::free)>;
  Buffer buffer(yyjson_mut_write(document.get(), YYJSON_WRITE_NOFLAG, &length),
                &std::free);
  if (!buffer || length < 2U || buffer.get()[0] != '{' ||
      buffer.get()[length - 1U] != '}') {
    throw Error("failed to write deterministic immersed-flow resolved JSON");
  }
  return std::string(buffer.get() + 1, length - 2U);
}

std::string serialize_immersed_flow(const ImmersedFlowCaseConfig &config) {
  const ImmersedFlowCaseConfig canonical = canonicalize_immersed_flow(config);
  std::string common = to_resolved_json(ResolvedCase(canonical.common_flow));
  constexpr std::string_view version2 = "{\"schema_version\":2,";
  if (common.compare(0U, version2.size(), version2) != 0 || common.empty() ||
      common.back() != '}') {
    throw Error("flow canonical JSON does not have the frozen envelope");
  }
  common[version2.size() - 2U] = '3';
  common.pop_back();
  common += ',';
  common += serialize_immersed_flow_tail(canonical);
  common += '}';
  return common;
}

} // namespace

namespace detail {

ResolvedCaseV3
parse_resolved_case_v3_json(std::string contents,
                            const std::filesystem::path &authoritative_path) {
  Document document = parse_document(contents);
  yyjson_val *root = yyjson_doc_get_root(document.get());
  const int version = dispatch_schema_version(root);
  if (version == 1) {
    return ResolvedCaseV3(
        load_case_config_from_json(contents, authoritative_path));
  }
  if (version == 2) {
    return ResolvedCaseV3(parse_common_flow(root, authoritative_path, false));
  }
  if (version == 3) {
    return ResolvedCaseV3(parse_immersed_flow(root, authoritative_path));
  }
  throw ConfigError("/schema_version",
                    "unsupported schema version " + std::to_string(version));
}

} // namespace detail

ResolvedCaseV3 load_resolved_case_v3(const std::filesystem::path &path) {
  std::string contents = read_complete_file(path);
  Document document = parse_document(contents);
  const int version =
      dispatch_schema_version(yyjson_doc_get_root(document.get()));
  if (version == 1 || version == 2) {
    const ResolvedCase legacy = load_resolved_case(path);
    return std::visit(
        [](const auto &value) -> ResolvedCaseV3 {
          return ResolvedCaseV3(value);
        },
        legacy);
  }
  if (version == 3) {
    return detail::parse_resolved_case_v3_json(std::move(contents), path);
  }
  throw ConfigError("/schema_version",
                    "unsupported schema version " + std::to_string(version));
}

std::string to_resolved_json_v3(const ResolvedCaseV3 &resolved_case) {
  return std::visit(
      [](const auto &value) -> std::string {
        using Config = std::decay_t<decltype(value)>;
        if constexpr (std::is_same_v<Config, ImmersedFlowCaseConfig>) {
          return serialize_immersed_flow(value);
        } else {
          return to_resolved_json(ResolvedCase(value));
        }
      },
      resolved_case);
}

} // namespace hundun::config
