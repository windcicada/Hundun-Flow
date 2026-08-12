// SPDX-License-Identifier: Apache-2.0

#include "hundun/cfg_resolved_case_v4_loader.hpp"

#include "cfg_resolved_case_v3_loader_detail.hpp"
#include "cfg_resolved_case_v4_loader_detail.hpp"

#include "hundun/cfg_resolved_case_loader.hpp"
#include "hundun/cfg_resolved_case_v3_loader.hpp"
#include "hundun/rt_error.hpp"
#include "yyjson.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <limits>
#include <memory>
#include <string_view>
#include <utility>

namespace hundun::config {
namespace {

using runtime::ConfigError;
using runtime::Error;
using Document = std::unique_ptr<yyjson_doc, decltype(&yyjson_doc_free)>;
using MutableDocument =
    std::unique_ptr<yyjson_mut_doc, decltype(&yyjson_mut_doc_free)>;

std::string child(std::string_view parent, std::string_view key) {
  return std::string(parent) + "/" + std::string(key);
}

std::string index_pointer(std::string_view parent, std::size_t index) {
  return child(parent, std::to_string(index));
}

std::string type_name(yyjson_val *value) {
  return value == nullptr ? "missing" : yyjson_get_type_desc(value);
}

[[noreturn]] void type_error(std::string pointer, std::string_view expected,
                             yyjson_val *actual) {
  throw ConfigError(std::move(pointer), "expected " + std::string(expected) +
                                            ", got " + type_name(actual));
}

yyjson_val *member(yyjson_val *object, const char *key,
                   std::string_view pointer) {
  yyjson_val *value = yyjson_obj_get(object, key);
  if (value == nullptr) {
    throw ConfigError(child(pointer, key), "required member is missing");
  }
  return value;
}

yyjson_val *object_member(yyjson_val *object, const char *key,
                          std::string_view pointer) {
  yyjson_val *value = member(object, key, pointer);
  if (!yyjson_is_obj(value)) {
    type_error(child(pointer, key), "object", value);
  }
  return value;
}

void reject_keys(yyjson_val *object,
                 std::initializer_list<std::string_view> allowed,
                 std::string_view pointer) {
  std::vector<std::string> seen;
  yyjson_obj_iter iterator = yyjson_obj_iter_with(object);
  while (yyjson_val *key_value = yyjson_obj_iter_next(&iterator)) {
    const std::string key(yyjson_get_str(key_value), yyjson_get_len(key_value));
    if (std::find(allowed.begin(), allowed.end(), key) == allowed.end()) {
      throw ConfigError(child(pointer, key), "unknown member");
    }
    if (std::find(seen.begin(), seen.end(), key) != seen.end()) {
      throw ConfigError(child(pointer, key), "duplicate member");
    }
    seen.push_back(key);
  }
}

std::string string_value(yyjson_val *value, std::string pointer) {
  if (!yyjson_is_str(value)) {
    type_error(std::move(pointer), "string", value);
  }
  return std::string(yyjson_get_str(value), yyjson_get_len(value));
}

std::string required_string(yyjson_val *object, const char *key,
                            std::string_view pointer) {
  return string_value(member(object, key, pointer), child(pointer, key));
}

double number_value(yyjson_val *value, std::string pointer) {
  if (!yyjson_is_num(value)) {
    type_error(std::move(pointer), "number", value);
  }
  const double result = yyjson_get_num(value);
  if (!std::isfinite(result)) {
    throw ConfigError(std::move(pointer), "expected a finite number");
  }
  return result;
}

double positive_number(yyjson_val *object, const char *key,
                       std::string_view pointer) {
  const std::string field = child(pointer, key);
  const double value = number_value(member(object, key, pointer), field);
  if (value <= 0.0) {
    throw ConfigError(field, "expected a finite value greater than zero");
  }
  return value;
}

int positive_integer(yyjson_val *object, const char *key,
                     std::string_view pointer) {
  const std::string field = child(pointer, key);
  yyjson_val *value = member(object, key, pointer);
  if (!yyjson_is_int(value)) {
    type_error(field, "integer", value);
  }
  const std::uint64_t converted = yyjson_get_uint(value);
  if (!yyjson_is_uint(value) || converted == 0U ||
      converted > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
    throw ConfigError(field, "expected a positive C++ int");
  }
  return static_cast<int>(converted);
}

Document parse(std::string &contents) {
  yyjson_read_err read_error{};
  Document document(yyjson_read_opts(contents.data(), contents.size(),
                                     YYJSON_READ_NOFLAG, nullptr, &read_error),
                    &yyjson_doc_free);
  if (!document) {
    throw ConfigError("",
                      "invalid JSON at byte " + std::to_string(read_error.pos));
  }
  return document;
}

std::string read_file(const std::filesystem::path &path) {
  std::ifstream input(path, std::ios::binary);
  if (!input) {
    throw ConfigError("", "could not open case JSON: " + path.string());
  }
  const std::string contents((std::istreambuf_iterator<char>(input)),
                             std::istreambuf_iterator<char>());
  if (input.bad()) {
    throw ConfigError("",
                      "could not read complete case JSON: " + path.string());
  }
  return contents;
}

bool parent_component(const std::filesystem::path &path) {
  return std::any_of(path.begin(), path.end(),
                     [](const auto &part) { return part == ".."; });
}

std::filesystem::path mechanism_path(const std::string &raw) {
  const std::filesystem::path path(raw);
  if (raw.empty() || raw.find('\0') != std::string::npos ||
      path.is_absolute() || path.has_root_name() || path.has_root_directory() ||
      parent_component(path) || path.lexically_normal().empty() ||
      path.lexically_normal() == ".") {
    throw ConfigError("/reacting/chemistry/mechanism/file",
                      "expected a case-root-confined relative path");
  }
  return path.lexically_normal();
}

bool lower_sha256(std::string_view value) {
  return value.size() == 64U &&
         std::all_of(value.begin(), value.end(), [](char character) {
           return (character >= '0' && character <= '9') ||
                  (character >= 'a' && character <= 'f');
         });
}

std::uint64_t species_fingerprint(const std::vector<std::string> &names) {
  std::uint64_t hash = UINT64_C(14695981039346656037);
  for (const auto &name : names) {
    for (const char character : name) {
      hash ^= static_cast<unsigned char>(character);
      hash *= UINT64_C(1099511628211);
    }
    hash ^= UINT64_C(255);
    hash *= UINT64_C(1099511628211);
  }
  return hash == 0U ? 1U : hash;
}

PatchName patch_name(const std::string &value, const std::string &pointer) {
  static constexpr std::array<std::string_view, 6> names = {
      "x_min", "x_max", "y_min", "y_max", "z_min", "z_max"};
  const auto found = std::find(names.begin(), names.end(), value);
  if (found == names.end()) {
    throw ConfigError(pointer, "unknown boundary patch '" + value + "'");
  }
  return static_cast<PatchName>(std::distance(names.begin(), found));
}

std::string mutable_json(MutableDocument &document) {
  std::size_t length = 0U;
  using Buffer = std::unique_ptr<char, decltype(&std::free)>;
  Buffer buffer(yyjson_mut_write(document.get(), YYJSON_WRITE_NOFLAG, &length),
                &std::free);
  if (!buffer) {
    throw Error("failed to allocate schema-v4 intermediate JSON");
  }
  return std::string(buffer.get(), length);
}

FlowCaseConfig parse_common_authority(yyjson_doc *source,
                                      const std::filesystem::path &path) {
  MutableDocument copy(yyjson_doc_mut_copy(source, nullptr),
                       &yyjson_mut_doc_free);
  if (!copy) {
    throw Error("failed to allocate schema-v4 common-flow projection");
  }
  yyjson_mut_val *root = yyjson_mut_doc_get_root(copy.get());
  yyjson_mut_obj_remove_str(root, "reacting");
  yyjson_mut_obj_remove_str(root, "immersed_boundary");
  yyjson_mut_obj_remove_str(root, "les");
  yyjson_mut_obj_remove_str(root, "schema_version");
  yyjson_mut_obj_add_int(copy.get(), root, "schema_version", 2);
  yyjson_mut_val *simulation = yyjson_mut_obj_get(root, "simulation");
  yyjson_mut_obj_remove_str(simulation, "type");
  yyjson_mut_obj_add_str(copy.get(), simulation, "type",
                         "variable_density_flow");
  yyjson_mut_obj_remove_str(simulation, "density_model");
  yyjson_mut_obj_add_str(copy.get(), simulation, "density_model", "constant");
  yyjson_mut_val *boundaries = yyjson_mut_obj_get(root, "boundaries");
  for (std::size_t index = 0; index < yyjson_mut_arr_size(boundaries);
       ++index) {
    yyjson_mut_obj_remove_str(yyjson_mut_arr_get(boundaries, index),
                              "reacting");
  }
  std::string common_json = mutable_json(copy);
  const ResolvedCaseV3 parsed =
      detail::parse_resolved_case_v3_json(std::move(common_json), path);
  FlowCaseConfig result = std::get<FlowCaseConfig>(parsed);
  result.density_model = DensityModel::ideal_gas;
  return result;
}

void parse_stage3_options(yyjson_doc *source, yyjson_val *root,
                          const std::filesystem::path &path,
                          ResolvedReactingCaseV4 &result) {
  yyjson_val *ibm = object_member(root, "immersed_boundary", "");
  yyjson_val *les = object_member(root, "les", "");
  const std::string ibm_model =
      required_string(ibm, "model", "/immersed_boundary");
  const std::string les_model = required_string(les, "model", "/les");
  if (ibm_model == "none" && les_model == "none") {
    reject_keys(ibm, {"model"}, "/immersed_boundary");
    reject_keys(les, {"model"}, "/les");
    result.immersed_boundary.model = ImmersedBoundaryModel::none;
    result.les.model = LesModel::none;
    return;
  }
  MutableDocument copy(yyjson_doc_mut_copy(source, nullptr),
                       &yyjson_mut_doc_free);
  if (!copy) {
    throw Error("failed to allocate schema-v4 Stage-3 projection");
  }
  yyjson_mut_val *mutable_root = yyjson_mut_doc_get_root(copy.get());
  yyjson_mut_obj_remove_str(mutable_root, "reacting");
  yyjson_mut_obj_remove_str(mutable_root, "schema_version");
  yyjson_mut_obj_add_int(copy.get(), mutable_root, "schema_version", 3);
  yyjson_mut_val *simulation = yyjson_mut_obj_get(mutable_root, "simulation");
  yyjson_mut_obj_remove_str(simulation, "type");
  yyjson_mut_obj_add_str(copy.get(), simulation, "type",
                         "variable_density_flow");
  yyjson_mut_obj_remove_str(simulation, "density_model");
  yyjson_mut_obj_add_str(copy.get(), simulation, "density_model", "constant");
  yyjson_mut_val *boundaries = yyjson_mut_obj_get(mutable_root, "boundaries");
  for (std::size_t index = 0; index < yyjson_mut_arr_size(boundaries);
       ++index) {
    yyjson_mut_obj_remove_str(yyjson_mut_arr_get(boundaries, index),
                              "reacting");
  }
  std::string projection = mutable_json(copy);
  const ResolvedCaseV3 parsed =
      detail::parse_resolved_case_v3_json(std::move(projection), path);
  const auto &stage3 = std::get<ImmersedFlowCaseConfig>(parsed);
  result.immersed_boundary = stage3.immersed_boundary;
  result.les = stage3.les;
}

void parse_boundaries(yyjson_val *root, ResolvedReactingCaseV4 &result) {
  yyjson_val *boundaries = member(root, "boundaries", "");
  for (std::size_t input_index = 0; input_index < yyjson_arr_size(boundaries);
       ++input_index) {
    yyjson_val *boundary = yyjson_arr_get(boundaries, input_index);
    const std::string pointer = index_pointer("/boundaries", input_index);
    const PatchName patch = patch_name(
        required_string(boundary, "patch", pointer), child(pointer, "patch"));
    auto &output = result.boundary_reacting[static_cast<std::size_t>(patch)];
    yyjson_val *reacting = yyjson_obj_get(boundary, "reacting");
    if (reacting == nullptr) {
      continue;
    }
    if (!yyjson_is_obj(reacting)) {
      type_error(child(pointer, "reacting"), "object", reacting);
    }
    const std::string reacting_pointer = child(pointer, "reacting");
    reject_keys(reacting, {"species", "thermal"}, reacting_pointer);
    const std::string species =
        required_string(reacting, "species", reacting_pointer);
    if (species != "non_catalytic_impermeable") {
      throw ConfigError(child(reacting_pointer, "species"),
                        "expected 'non_catalytic_impermeable'");
    }
    output.non_catalytic_impermeable = true;
    yyjson_val *thermal = object_member(reacting, "thermal", reacting_pointer);
    const std::string thermal_pointer = child(reacting_pointer, "thermal");
    reject_keys(thermal, {"mode", "temperature_k"}, thermal_pointer);
    ReactingThermalBoundaryConfig thermal_config{};
    const std::string mode = required_string(thermal, "mode", thermal_pointer);
    yyjson_val *temperature = yyjson_obj_get(thermal, "temperature_k");
    if (mode == "adiabatic") {
      if (temperature != nullptr) {
        throw ConfigError(child(thermal_pointer, "temperature_k"),
                          "member is forbidden for an adiabatic wall");
      }
      thermal_config.mode = ReactingThermalMode::adiabatic;
    } else if (mode == "isothermal") {
      if (temperature == nullptr) {
        throw ConfigError(child(thermal_pointer, "temperature_k"),
                          "required member is missing");
      }
      thermal_config.mode = ReactingThermalMode::isothermal;
      thermal_config.temperature_k =
          number_value(temperature, child(thermal_pointer, "temperature_k"));
      if (*thermal_config.temperature_k <= 0.0) {
        throw ConfigError(child(thermal_pointer, "temperature_k"),
                          "expected a finite value greater than zero");
      }
    } else {
      throw ConfigError(child(thermal_pointer, "mode"),
                        "expected 'adiabatic' or 'isothermal'");
    }
    output.thermal = thermal_config;
  }
  for (std::size_t index = 0; index < result.common_flow.boundaries.size();
       ++index) {
    if (result.common_flow.boundaries[index].type ==
            BoundaryType::no_slip_wall &&
        !result.boundary_reacting[index].thermal.has_value()) {
      throw ConfigError(index_pointer("/boundaries", index) + "/reacting",
                        "reacting wall contract is required");
    }
  }
}

void parse_reacting(yyjson_val *root, ResolvedReactingCaseV4 &result) {
  yyjson_val *reacting = object_member(root, "reacting", "");
  reject_keys(
      reacting,
      {"chemistry", "thermodynamics", "transport", "pressure_constraint"},
      "/reacting");
  yyjson_val *chemistry = object_member(reacting, "chemistry", "/reacting");
  reject_keys(chemistry,
              {"backend", "mechanism", "relative_tolerance",
               "absolute_tolerance", "maximum_internal_steps"},
              "/reacting/chemistry");
  if (required_string(chemistry, "backend", "/reacting/chemistry") !=
      "cantera") {
    throw ConfigError("/reacting/chemistry/backend", "expected 'cantera'");
  }
  yyjson_val *mechanism =
      object_member(chemistry, "mechanism", "/reacting/chemistry");
  reject_keys(mechanism, {"file", "sha256", "phase"},
              "/reacting/chemistry/mechanism");
  result.mechanism.file = mechanism_path(
      required_string(mechanism, "file", "/reacting/chemistry/mechanism"));
  result.mechanism.sha256 =
      required_string(mechanism, "sha256", "/reacting/chemistry/mechanism");
  if (!lower_sha256(result.mechanism.sha256)) {
    throw ConfigError("/reacting/chemistry/mechanism/sha256",
                      "expected 64 lowercase hexadecimal characters");
  }
  result.mechanism.phase =
      required_string(mechanism, "phase", "/reacting/chemistry/mechanism");
  if (result.mechanism.phase.empty()) {
    throw ConfigError("/reacting/chemistry/mechanism/phase",
                      "expected a nonempty string");
  }
  result.chemistry.relative_tolerance =
      positive_number(chemistry, "relative_tolerance", "/reacting/chemistry");
  result.chemistry.absolute_tolerance =
      positive_number(chemistry, "absolute_tolerance", "/reacting/chemistry");
  result.chemistry.maximum_internal_steps = positive_integer(
      chemistry, "maximum_internal_steps", "/reacting/chemistry");

  yyjson_val *thermodynamics =
      object_member(reacting, "thermodynamics", "/reacting");
  reject_keys(
      thermodynamics,
      {"initial_p0_pa", "initial_temperature_k", "initial_mass_fractions"},
      "/reacting/thermodynamics");
  result.initial_p0_pa = positive_number(thermodynamics, "initial_p0_pa",
                                         "/reacting/thermodynamics");
  result.initial_temperature_k = positive_number(
      thermodynamics, "initial_temperature_k", "/reacting/thermodynamics");
  yyjson_val *fractions = member(thermodynamics, "initial_mass_fractions",
                                 "/reacting/thermodynamics");
  if (!yyjson_is_arr(fractions)) {
    type_error("/reacting/thermodynamics/initial_mass_fractions", "array",
               fractions);
  }
  if (yyjson_arr_size(fractions) == 0U) {
    throw ConfigError("/reacting/thermodynamics/initial_mass_fractions",
                      "expected at least one species");
  }
  double total = 0.0;
  for (std::size_t index = 0; index < yyjson_arr_size(fractions); ++index) {
    yyjson_val *entry = yyjson_arr_get(fractions, index);
    const std::string pointer =
        index_pointer("/reacting/thermodynamics/initial_mass_fractions", index);
    if (!yyjson_is_obj(entry)) {
      type_error(pointer, "object", entry);
    }
    reject_keys(entry, {"species", "value"}, pointer);
    const std::string species = required_string(entry, "species", pointer);
    if (species.empty() ||
        std::find(result.species_names.begin(), result.species_names.end(),
                  species) != result.species_names.end()) {
      throw ConfigError(child(pointer, "species"),
                        "species names must be nonempty and unique");
    }
    const double fraction =
        number_value(member(entry, "value", pointer), child(pointer, "value"));
    if (fraction < 0.0 || fraction > 1.0) {
      throw ConfigError(child(pointer, "value"),
                        "mass fraction must be in [0,1]");
    }
    result.species_names.push_back(species);
    result.initial_mass_fractions.push_back(fraction);
    total += fraction;
  }
  if (!std::isfinite(total) ||
      std::abs(total - 1.0) >
          1.0e-12 * static_cast<double>(result.species_names.size())) {
    throw ConfigError("/reacting/thermodynamics/initial_mass_fractions",
                      "mass fractions must sum to one");
  }
  result.composition_fingerprint = species_fingerprint(result.species_names);

  yyjson_val *transport = object_member(reacting, "transport", "/reacting");
  reject_keys(transport, {"model"}, "/reacting/transport");
  if (required_string(transport, "model", "/reacting/transport") !=
      "mixture_averaged") {
    throw ConfigError("/reacting/transport/model",
                      "expected 'mixture_averaged'");
  }
  yyjson_val *pressure =
      object_member(reacting, "pressure_constraint", "/reacting");
  reject_keys(pressure, {"mode"}, "/reacting/pressure_constraint");
  const std::string mode =
      required_string(pressure, "mode", "/reacting/pressure_constraint");
  if (mode == "open_fixed_p0") {
    result.pressure_mode = PressureConstraintMode::open_fixed_p0;
  } else if (mode == "closed") {
    result.pressure_mode = PressureConstraintMode::closed;
  } else if (mode == "partially_closed") {
    result.pressure_mode = PressureConstraintMode::partially_closed;
  } else {
    throw ConfigError("/reacting/pressure_constraint/mode",
                      "unknown pressure-constraint mode '" + mode + "'");
  }
}

void require_write(bool ok) {
  if (!ok) {
    throw Error("failed to allocate deterministic schema-v4 JSON");
  }
}

const char *pressure_name(PressureConstraintMode mode) {
  switch (mode) {
  case PressureConstraintMode::open_fixed_p0:
    return "open_fixed_p0";
  case PressureConstraintMode::closed:
    return "closed";
  case PressureConstraintMode::partially_closed:
    return "partially_closed";
  }
  throw ConfigError("/reacting/pressure_constraint/mode", "invalid mode");
}

std::string serialize(ResolvedReactingCaseV4 config) {
  if (config.schema_version != 4 || config.common_flow.schema_version != 2) {
    throw ConfigError("/schema_version", "expected schema version 4");
  }
  FlowCaseConfig common_projection = config.common_flow;
  common_projection.density_model = DensityModel::constant;
  const bool has_stage3 =
      config.immersed_boundary.model != ImmersedBoundaryModel::none ||
      config.les.model != LesModel::none;
  std::string base;
  if (has_stage3) {
    ImmersedFlowCaseConfig stage3{};
    stage3.schema_version = 3;
    stage3.common_flow = common_projection;
    stage3.immersed_boundary = config.immersed_boundary;
    stage3.les = config.les;
    base = to_resolved_json_v3(ResolvedCaseV3(std::move(stage3)));
  } else {
    base = to_resolved_json(ResolvedCase(common_projection));
  }
  Document immutable = parse(base);
  MutableDocument document(yyjson_doc_mut_copy(immutable.get(), nullptr),
                           &yyjson_mut_doc_free);
  yyjson_mut_val *root = yyjson_mut_doc_get_root(document.get());
  yyjson_mut_obj_remove_str(root, "schema_version");
  require_write(
      yyjson_mut_obj_add_int(document.get(), root, "schema_version", 4));
  yyjson_mut_val *simulation = yyjson_mut_obj_get(root, "simulation");
  yyjson_mut_obj_remove_str(simulation, "type");
  require_write(yyjson_mut_obj_add_str(document.get(), simulation, "type",
                                       "reacting_flow"));
  yyjson_mut_obj_remove_str(simulation, "density_model");
  require_write(yyjson_mut_obj_add_str(document.get(), simulation,
                                       "density_model", "ideal_gas"));

  if (!has_stage3) {
    yyjson_mut_val *ibm =
        yyjson_mut_obj_add_obj(document.get(), root, "immersed_boundary");
    yyjson_mut_val *les = yyjson_mut_obj_add_obj(document.get(), root, "les");
    require_write(yyjson_mut_obj_add_str(document.get(), ibm, "model", "none"));
    require_write(yyjson_mut_obj_add_str(document.get(), les, "model", "none"));
  }

  yyjson_mut_val *boundaries = yyjson_mut_obj_get(root, "boundaries");
  for (std::size_t index = 0; index < config.boundary_reacting.size();
       ++index) {
    const auto &boundary = config.boundary_reacting[index];
    if (!boundary.thermal.has_value())
      continue;
    if (!boundary.non_catalytic_impermeable ||
        config.common_flow.boundaries[index].type !=
            BoundaryType::no_slip_wall) {
      throw ConfigError(index_pointer("/boundaries", index) + "/reacting",
                        "invalid reacting wall contract");
    }
    yyjson_mut_val *entry = yyjson_mut_arr_get(boundaries, index);
    yyjson_mut_val *reacting =
        yyjson_mut_obj_add_obj(document.get(), entry, "reacting");
    require_write(yyjson_mut_obj_add_str(document.get(), reacting, "species",
                                         "non_catalytic_impermeable"));
    yyjson_mut_val *thermal =
        yyjson_mut_obj_add_obj(document.get(), reacting, "thermal");
    if (boundary.thermal->mode == ReactingThermalMode::adiabatic) {
      if (boundary.thermal->temperature_k.has_value()) {
        throw ConfigError(index_pointer("/boundaries", index) +
                              "/reacting/thermal/temperature_k",
                          "member is forbidden for an adiabatic wall");
      }
      require_write(
          yyjson_mut_obj_add_str(document.get(), thermal, "mode", "adiabatic"));
    } else {
      if (!boundary.thermal->temperature_k.has_value() ||
          !std::isfinite(*boundary.thermal->temperature_k) ||
          *boundary.thermal->temperature_k <= 0.0) {
        throw ConfigError(index_pointer("/boundaries", index) +
                              "/reacting/thermal/temperature_k",
                          "required finite positive value");
      }
      require_write(yyjson_mut_obj_add_str(document.get(), thermal, "mode",
                                           "isothermal"));
      require_write(yyjson_mut_obj_add_real(document.get(), thermal,
                                            "temperature_k",
                                            *boundary.thermal->temperature_k));
    }
  }

  yyjson_mut_val *reacting =
      yyjson_mut_obj_add_obj(document.get(), root, "reacting");
  yyjson_mut_val *chemistry =
      yyjson_mut_obj_add_obj(document.get(), reacting, "chemistry");
  require_write(
      yyjson_mut_obj_add_str(document.get(), chemistry, "backend", "cantera"));
  yyjson_mut_val *mechanism =
      yyjson_mut_obj_add_obj(document.get(), chemistry, "mechanism");
  require_write(yyjson_mut_obj_add_strcpy(
      document.get(), mechanism, "file",
      config.mechanism.file.generic_string().c_str()));
  require_write(yyjson_mut_obj_add_strcpy(document.get(), mechanism, "sha256",
                                          config.mechanism.sha256.c_str()));
  require_write(yyjson_mut_obj_add_strcpy(document.get(), mechanism, "phase",
                                          config.mechanism.phase.c_str()));
  require_write(yyjson_mut_obj_add_real(document.get(), chemistry,
                                        "relative_tolerance",
                                        config.chemistry.relative_tolerance));
  require_write(yyjson_mut_obj_add_real(document.get(), chemistry,
                                        "absolute_tolerance",
                                        config.chemistry.absolute_tolerance));
  require_write(yyjson_mut_obj_add_int(
      document.get(), chemistry, "maximum_internal_steps",
      config.chemistry.maximum_internal_steps));
  yyjson_mut_val *thermo =
      yyjson_mut_obj_add_obj(document.get(), reacting, "thermodynamics");
  require_write(yyjson_mut_obj_add_real(document.get(), thermo, "initial_p0_pa",
                                        config.initial_p0_pa));
  require_write(yyjson_mut_obj_add_real(document.get(), thermo,
                                        "initial_temperature_k",
                                        config.initial_temperature_k));
  if (config.species_names.empty() ||
      config.species_names.size() != config.initial_mass_fractions.size()) {
    throw ConfigError("/reacting/thermodynamics/initial_mass_fractions",
                      "one value is required for every species");
  }
  yyjson_mut_val *fractions =
      yyjson_mut_obj_add_arr(document.get(), thermo, "initial_mass_fractions");
  for (std::size_t index = 0; index < config.species_names.size(); ++index) {
    yyjson_mut_val *entry = yyjson_mut_obj(document.get());
    require_write(yyjson_mut_arr_add_val(fractions, entry));
    require_write(yyjson_mut_obj_add_strcpy(
        document.get(), entry, "species", config.species_names[index].c_str()));
    require_write(yyjson_mut_obj_add_real(
        document.get(), entry, "value", config.initial_mass_fractions[index]));
  }
  yyjson_mut_val *transport =
      yyjson_mut_obj_add_obj(document.get(), reacting, "transport");
  require_write(yyjson_mut_obj_add_str(document.get(), transport, "model",
                                       "mixture_averaged"));
  yyjson_mut_val *pressure =
      yyjson_mut_obj_add_obj(document.get(), reacting, "pressure_constraint");
  require_write(yyjson_mut_obj_add_str(document.get(), pressure, "mode",
                                       pressure_name(config.pressure_mode)));

  std::string output = mutable_json(document);
  static_cast<void>(detail::parse_resolved_reacting_case_v4_json(
      output, "schema-v4-in-memory.json"));
  return output;
}

} // namespace

namespace detail {

ResolvedReactingCaseV4 parse_resolved_reacting_case_v4_json(
    std::string contents, const std::filesystem::path &authoritative_path) {
  Document document = parse(contents);
  yyjson_val *root = yyjson_doc_get_root(document.get());
  if (!yyjson_is_obj(root))
    type_error("", "object", root);
  reject_keys(root,
              {"schema_version", "case", "simulation", "resources", "mesh",
               "time", "physics", "scalars", "boundaries", "restart",
               "diagnostics", "performance", "immersed_boundary", "les",
               "reacting"},
              "");
  yyjson_val *version = member(root, "schema_version", "");
  if (!yyjson_is_int(version) || yyjson_get_sint(version) != 4) {
    throw ConfigError("/schema_version", "expected schema version 4");
  }
  yyjson_val *simulation = object_member(root, "simulation", "");
  if (required_string(simulation, "type", "/simulation") != "reacting_flow") {
    throw ConfigError("/simulation/type", "expected 'reacting_flow'");
  }
  if (required_string(simulation, "density_model", "/simulation") !=
      "ideal_gas") {
    throw ConfigError("/simulation/density_model", "expected 'ideal_gas'");
  }

  ResolvedReactingCaseV4 result{};
  result.common_flow =
      parse_common_authority(document.get(), authoritative_path);
  parse_stage3_options(document.get(), root, authoritative_path, result);
  parse_reacting(root, result);
  parse_boundaries(root, result);
  return result;
}

} // namespace detail

ResolvedReactingCaseV4
load_resolved_reacting_case_v4(const std::filesystem::path &path) {
  return detail::parse_resolved_reacting_case_v4_json(read_file(path), path);
}

std::string
to_resolved_reacting_json_v4(const ResolvedReactingCaseV4 &resolved_case) {
  return serialize(resolved_case);
}

} // namespace hundun::config
