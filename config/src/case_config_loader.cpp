// SPDX-License-Identifier: Apache-2.0

#include "hundun/config/case_config_loader.hpp"

#include "hundun/runtime/error.hpp"
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

void reject_unknown_keys(
    yyjson_val* object, std::initializer_list<std::string_view> allowed,
    std::string_view pointer) {
  std::vector<std::string> seen;
  seen.reserve(yyjson_obj_size(object));
  yyjson_obj_iter iterator = yyjson_obj_iter_with(object);
  while (yyjson_val* key_value = yyjson_obj_iter_next(&iterator)) {
    const std::string key(yyjson_get_str(key_value), yyjson_get_len(key_value));
    const bool is_allowed =
        std::find(allowed.begin(), allowed.end(), std::string_view(key)) !=
        allowed.end();
    if (!is_allowed) {
      throw ConfigError(child_pointer(pointer, key), "unknown member");
    }
    if (std::find(seen.begin(), seen.end(), key) != seen.end()) {
      throw ConfigError(child_pointer(pointer, key), "duplicate member");
    }
    seen.push_back(key);
  }
}

double read_number_value(yyjson_val* value, std::string pointer) {
  if (!yyjson_is_num(value)) {
    throw_type_error(std::move(pointer), "number", value);
  }
  const double result = yyjson_get_num(value);
  if (!std::isfinite(result)) {
    throw ConfigError(std::move(pointer), "expected a finite number");
  }
  return result;
}

int read_integer_value(yyjson_val* value, std::string pointer) {
  if (!yyjson_is_int(value)) {
    throw_type_error(std::move(pointer), "integer", value);
  }

  if (yyjson_is_uint(value)) {
    const std::uint64_t number = yyjson_get_uint(value);
    if (number > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
      throw ConfigError(std::move(pointer),
                        "integer is outside the C++ int range: " +
                            std::to_string(number));
    }
    return static_cast<int>(number);
  }

  const std::int64_t number = yyjson_get_sint(value);
  if (number < static_cast<std::int64_t>(std::numeric_limits<int>::min()) ||
      number > static_cast<std::int64_t>(std::numeric_limits<int>::max())) {
    throw ConfigError(std::move(pointer),
                      "integer is outside the C++ int range: " +
                          std::to_string(number));
  }
  return static_cast<int>(number);
}

double require_number(yyjson_val* object, const char* key,
                      std::string_view pointer) {
  return read_number_value(require_member(object, key, pointer),
                           child_pointer(pointer, key));
}

int require_integer(yyjson_val* object, const char* key,
                    std::string_view pointer) {
  return read_integer_value(require_member(object, key, pointer),
                            child_pointer(pointer, key));
}

std::string require_string(yyjson_val* object, const char* key,
                           std::string_view pointer) {
  yyjson_val* value = require_member(object, key, pointer);
  if (!yyjson_is_str(value)) {
    throw_type_error(child_pointer(pointer, key), "string", value);
  }
  return std::string(yyjson_get_str(value), yyjson_get_len(value));
}

runtime::Int3 require_int3(yyjson_val* object, const char* key,
                          std::string_view pointer) {
  const std::string field_pointer = child_pointer(pointer, key);
  yyjson_val* array = require_member(object, key, pointer);
  if (!yyjson_is_arr(array)) {
    throw_type_error(field_pointer, "array of three integers", array);
  }
  if (yyjson_arr_size(array) != 3) {
    throw ConfigError(field_pointer, "expected exactly three entries, got " +
                                         std::to_string(yyjson_arr_size(array)));
  }
  return {read_integer_value(yyjson_arr_get(array, 0),
                             index_pointer(field_pointer, 0)),
          read_integer_value(yyjson_arr_get(array, 1),
                             index_pointer(field_pointer, 1)),
          read_integer_value(yyjson_arr_get(array, 2),
                             index_pointer(field_pointer, 2))};
}

runtime::Real3 require_real3(yyjson_val* object, const char* key,
                            std::string_view pointer) {
  const std::string field_pointer = child_pointer(pointer, key);
  yyjson_val* array = require_member(object, key, pointer);
  if (!yyjson_is_arr(array)) {
    throw_type_error(field_pointer, "array of three numbers", array);
  }
  if (yyjson_arr_size(array) != 3) {
    throw ConfigError(field_pointer, "expected exactly three entries, got " +
                                         std::to_string(yyjson_arr_size(array)));
  }
  return {read_number_value(yyjson_arr_get(array, 0),
                            index_pointer(field_pointer, 0)),
          read_number_value(yyjson_arr_get(array, 1),
                            index_pointer(field_pointer, 1)),
          read_number_value(yyjson_arr_get(array, 2),
                            index_pointer(field_pointer, 2))};
}

std::array<bool, 3> require_bool3(yyjson_val* object, const char* key,
                                  std::string_view pointer) {
  const std::string field_pointer = child_pointer(pointer, key);
  yyjson_val* array = require_member(object, key, pointer);
  if (!yyjson_is_arr(array)) {
    throw_type_error(field_pointer, "array of three booleans", array);
  }
  if (yyjson_arr_size(array) != 3) {
    throw ConfigError(field_pointer, "expected exactly three entries, got " +
                                         std::to_string(yyjson_arr_size(array)));
  }

  std::array<bool, 3> result{};
  for (std::size_t index = 0; index < result.size(); ++index) {
    yyjson_val* value = yyjson_arr_get(array, index);
    if (!yyjson_is_bool(value)) {
      throw_type_error(index_pointer(field_pointer, index), "boolean", value);
    }
    result[index] = yyjson_get_bool(value);
  }
  return result;
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

std::filesystem::path normalize_output_directory(
    const std::string& raw, const std::filesystem::path& case_path) {
  constexpr const char* pointer = "/output/directory";
  if (raw.empty()) {
    throw ConfigError(pointer, "expected a non-empty relative path");
  }
  if (raw.find('\0') != std::string::npos) {
    throw ConfigError(pointer, "path must not contain a null byte");
  }

  const std::filesystem::path path(raw);
  if (path.is_absolute() || path.has_root_name() || path.has_root_directory()) {
    throw ConfigError(pointer, "expected an unrooted relative path");
  }
  if (has_parent_component(path)) {
    throw ConfigError(pointer, "path must not contain a '..' component");
  }

  const std::filesystem::path normalized = path.lexically_normal();
  if (normalized.empty()) {
    throw ConfigError(pointer, "expected a non-empty relative path");
  }

  try {
    const std::filesystem::path case_root =
        std::filesystem::absolute(case_path).lexically_normal().parent_path();
    const std::filesystem::path resolved =
        (case_root / normalized).lexically_normal();
    if (!has_path_prefix(resolved, case_root)) {
      throw ConfigError(pointer, "path escapes the case directory");
    }
  } catch (const std::filesystem::filesystem_error& error) {
    throw ConfigError(pointer,
                      "could not resolve path against the case directory: " +
                          std::string(error.what()));
  }
  return normalized;
}

void require_positive_component(double value, std::string pointer) {
  if (!std::isfinite(value) || value <= 0.0) {
    throw ConfigError(std::move(pointer),
                      "expected a finite value greater than zero, got " +
                          std::to_string(value));
  }
}

void require_finite_component(double value, std::string pointer) {
  if (!std::isfinite(value)) {
    throw ConfigError(std::move(pointer), "expected a finite number");
  }
}

void validate_output_path(const std::filesystem::path& path) {
  constexpr const char* pointer = "/output/directory";
  if (path.empty()) {
    throw ConfigError(pointer, "expected a non-empty relative path");
  }
  if (path.is_absolute() || path.has_root_name() || path.has_root_directory()) {
    throw ConfigError(pointer, "expected an unrooted relative path");
  }
  if (has_parent_component(path)) {
    throw ConfigError(pointer, "path must not contain a '..' component");
  }
  if (path.generic_string().find('\0') != std::string::npos) {
    throw ConfigError(pointer, "path must not contain a null byte");
  }
  if (path.generic_string() != path.lexically_normal().generic_string()) {
    throw ConfigError(pointer, "path must be lexically normalized");
  }
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

CaseConfig parse_case_config(yyjson_val* root,
                             const std::filesystem::path& case_path) {
  if (!yyjson_is_obj(root)) {
    throw_type_error("", "object", root);
  }
  reject_unknown_keys(root,
                      {"schema_version", "case", "resources", "mesh", "time",
                       "transport", "initial_condition", "output"},
                      "");

  CaseConfig config{};
  config.schema_version = require_integer(root, "schema_version", "");

  yyjson_val* case_object = require_object_member(root, "case", "");
  reject_unknown_keys(case_object, {"name"}, "/case");
  config.case_name = require_string(case_object, "name", "/case");

  if (yyjson_val* resources = yyjson_obj_get(root, "resources")) {
    if (!yyjson_is_obj(resources)) {
      throw_type_error("/resources", "object", resources);
    }
    reject_unknown_keys(resources, {"expected_ranks"}, "/resources");
    if (yyjson_val* expected_ranks =
            yyjson_obj_get(resources, "expected_ranks")) {
      config.expected_ranks =
          read_integer_value(expected_ranks, "/resources/expected_ranks");
    }
  }

  yyjson_val* mesh = require_object_member(root, "mesh", "");
  reject_unknown_keys(mesh, {"cells", "origin_m", "length_m", "periodic"},
                      "/mesh");
  config.mesh.cells = require_int3(mesh, "cells", "/mesh");
  config.mesh.origin_m = require_real3(mesh, "origin_m", "/mesh");
  config.mesh.length_m = require_real3(mesh, "length_m", "/mesh");
  config.mesh.periodic = require_bool3(mesh, "periodic", "/mesh");

  yyjson_val* time = require_object_member(root, "time", "");
  reject_unknown_keys(time, {"dt_s", "steps"}, "/time");
  config.time.dt_s = require_number(time, "dt_s", "/time");
  config.time.steps = require_integer(time, "steps", "/time");

  yyjson_val* transport = require_object_member(root, "transport", "");
  reject_unknown_keys(transport,
                      {"velocity_m_per_s", "diffusivity_m2_per_s"},
                      "/transport");
  config.transport.velocity_m_per_s =
      require_real3(transport, "velocity_m_per_s", "/transport");
  config.transport.diffusivity_m2_per_s =
      require_number(transport, "diffusivity_m2_per_s", "/transport");

  yyjson_val* initial_condition =
      require_object_member(root, "initial_condition", "");
  reject_unknown_keys(initial_condition, {"type"}, "/initial_condition");
  config.initial_condition =
      require_string(initial_condition, "type", "/initial_condition");

  yyjson_val* output = require_object_member(root, "output", "");
  reject_unknown_keys(output,
                      {"directory", "write_interval", "restart_interval"},
                      "/output");
  config.output.directory = normalize_output_directory(
      require_string(output, "directory", "/output"), case_path);
  config.output.write_interval =
      require_integer(output, "write_interval", "/output");
  config.output.restart_interval =
      require_integer(output, "restart_interval", "/output");

  validate_case_config(config);
  return config;
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

void add_int3(yyjson_mut_doc* document, yyjson_mut_val* parent,
              const char* key, runtime::Int3 value) {
  yyjson_mut_val* array = add_array(document, parent, key);
  require_json_write(yyjson_mut_arr_add_int(document, array, value.x));
  require_json_write(yyjson_mut_arr_add_int(document, array, value.y));
  require_json_write(yyjson_mut_arr_add_int(document, array, value.z));
}

void add_real3(yyjson_mut_doc* document, yyjson_mut_val* parent,
               const char* key, runtime::Real3 value) {
  yyjson_mut_val* array = add_array(document, parent, key);
  require_json_write(yyjson_mut_arr_add_real(document, array, value.x));
  require_json_write(yyjson_mut_arr_add_real(document, array, value.y));
  require_json_write(yyjson_mut_arr_add_real(document, array, value.z));
}

}  // namespace

CaseConfig load_case_config(const std::filesystem::path& path) {
  std::string contents = read_complete_file(path);
  yyjson_read_err read_error{};
  using Document = std::unique_ptr<yyjson_doc, decltype(&yyjson_doc_free)>;
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
  return parse_case_config(yyjson_doc_get_root(document.get()), path);
}

void validate_case_config(const CaseConfig& config) {
  if (config.schema_version != 1) {
    throw ConfigError(
        "/schema_version",
        "expected schema version 1, got " +
            std::to_string(config.schema_version));
  }
  if (config.expected_ranks.has_value() && *config.expected_ranks < 1) {
    throw ConfigError("/resources/expected_ranks",
                      "expected an integer of at least 1, got " +
                          std::to_string(*config.expected_ranks));
  }

  const std::array<int, 3> cells = {config.mesh.cells.x, config.mesh.cells.y,
                                    config.mesh.cells.z};
  for (std::size_t index = 0; index < cells.size(); ++index) {
    if (cells[index] < 1) {
      throw ConfigError(index_pointer("/mesh/cells", index),
                        "expected an integer of at least 1, got " +
                            std::to_string(cells[index]));
    }
  }

  const std::array<double, 3> origin = {config.mesh.origin_m.x,
                                        config.mesh.origin_m.y,
                                        config.mesh.origin_m.z};
  const std::array<double, 3> length = {config.mesh.length_m.x,
                                        config.mesh.length_m.y,
                                        config.mesh.length_m.z};
  for (std::size_t index = 0; index < origin.size(); ++index) {
    require_finite_component(origin[index],
                             index_pointer("/mesh/origin_m", index));
    require_positive_component(length[index],
                               index_pointer("/mesh/length_m", index));
  }

  require_positive_component(config.time.dt_s, "/time/dt_s");
  if (config.time.steps < 0) {
    throw ConfigError("/time/steps", "expected an integer of at least 0, got " +
                                         std::to_string(config.time.steps));
  }

  const std::array<double, 3> velocity = {
      config.transport.velocity_m_per_s.x,
      config.transport.velocity_m_per_s.y,
      config.transport.velocity_m_per_s.z};
  for (std::size_t index = 0; index < velocity.size(); ++index) {
    require_finite_component(
        velocity[index], index_pointer("/transport/velocity_m_per_s", index));
  }
  if (!std::isfinite(config.transport.diffusivity_m2_per_s) ||
      config.transport.diffusivity_m2_per_s < 0.0) {
    throw ConfigError("/transport/diffusivity_m2_per_s",
                      "expected a finite non-negative number, got " +
                          std::to_string(config.transport.diffusivity_m2_per_s));
  }

  if (config.initial_condition != "sine_x") {
    throw ConfigError("/initial_condition/type",
                      "expected 'sine_x', got '" + config.initial_condition +
                          "'");
  }

  validate_output_path(config.output.directory);
  if (config.output.write_interval < 1) {
    throw ConfigError("/output/write_interval",
                      "expected an integer of at least 1, got " +
                          std::to_string(config.output.write_interval));
  }
  if (config.output.restart_interval < 1) {
    throw ConfigError("/output/restart_interval",
                      "expected an integer of at least 1, got " +
                          std::to_string(config.output.restart_interval));
  }
}

std::string to_resolved_json(const CaseConfig& config) {
  validate_case_config(config);

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
  require_json_write(yyjson_mut_obj_add_strncpy(
      document.get(), case_object, "name", config.case_name.data(),
      config.case_name.size()));

  if (config.expected_ranks.has_value()) {
    yyjson_mut_val* resources = add_object(document.get(), root, "resources");
    require_json_write(yyjson_mut_obj_add_int(
        document.get(), resources, "expected_ranks", *config.expected_ranks));
  }

  yyjson_mut_val* mesh = add_object(document.get(), root, "mesh");
  add_int3(document.get(), mesh, "cells", config.mesh.cells);
  add_real3(document.get(), mesh, "origin_m", config.mesh.origin_m);
  add_real3(document.get(), mesh, "length_m", config.mesh.length_m);
  yyjson_mut_val* periodic = add_array(document.get(), mesh, "periodic");
  for (const bool value : config.mesh.periodic) {
    require_json_write(
        yyjson_mut_arr_add_bool(document.get(), periodic, value));
  }

  yyjson_mut_val* time = add_object(document.get(), root, "time");
  require_json_write(yyjson_mut_obj_add_real(document.get(), time, "dt_s",
                                             config.time.dt_s));
  require_json_write(yyjson_mut_obj_add_int(document.get(), time, "steps",
                                            config.time.steps));

  yyjson_mut_val* transport = add_object(document.get(), root, "transport");
  add_real3(document.get(), transport, "velocity_m_per_s",
            config.transport.velocity_m_per_s);
  require_json_write(yyjson_mut_obj_add_real(
      document.get(), transport, "diffusivity_m2_per_s",
      config.transport.diffusivity_m2_per_s));

  yyjson_mut_val* initial_condition =
      add_object(document.get(), root, "initial_condition");
  require_json_write(yyjson_mut_obj_add_strncpy(
      document.get(), initial_condition, "type",
      config.initial_condition.data(), config.initial_condition.size()));

  yyjson_mut_val* output = add_object(document.get(), root, "output");
  const std::string directory = config.output.directory.generic_string();
  require_json_write(yyjson_mut_obj_add_strncpy(
      document.get(), output, "directory", directory.data(), directory.size()));
  require_json_write(yyjson_mut_obj_add_int(
      document.get(), output, "write_interval", config.output.write_interval));
  require_json_write(yyjson_mut_obj_add_int(
      document.get(), output, "restart_interval",
      config.output.restart_interval));

  std::size_t length = 0;
  using Buffer = std::unique_ptr<char, decltype(&std::free)>;
  Buffer buffer(
      yyjson_mut_write(document.get(), YYJSON_WRITE_NOFLAG, &length), &std::free);
  if (!buffer) {
    throw Error("failed to write deterministic resolved JSON");
  }
  return std::string(buffer.get(), length);
}

}  // namespace hundun::config
